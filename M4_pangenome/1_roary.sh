#!/bin/bash
#SBATCH -A SNIC_PROJECT_ID
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 0-4:00:00
#SBATCH -J roary
#SBATCH -o SNIC_PROJECT/logs/slurm-%A_%a.out
#SBATCH --mail-type=ALL

module load bioinfo-tools Roary prokka FastTree
SPECIES="Rothia"
REFS=$(find SNIC_PROJECT/ADRIAN/genomes/s__$SPECIES/*.fna -type f)
MAGS=$(find SNIC_PROJECT/ADRIAN/calculus/Gorilla/mag/metaWRAP/full-decontam-set/phylophlan/$SPECIES/*.fa -type f | grep -v -E "BL|BE|BS|ERR|EXB|LIB")
PROKKA=SNIC_PROJECT/ADRIAN/calculus/Gorilla/mag/prokka/$SPECIES

OUTDIR=SNIC_PROJECT/ADRIAN/calculus/Gorilla/mag/roary/$SPECIES
mkdir -p $OUTDIR

threads=$SLURM_CPU_ON_NODE

INDIR=$OUTDIR/input
mkdir -p $INDIR

# prokka on the references
# skip if already done
for i in $REFS
  do
    # extract name
    name=$(basename "$i" _genomic.fna)
    # if file exists, skip
    if [[ -f $PROKKA/$name/$name.gff ]]
      then
        echo "prokka already done for $name"
      else
        # run prokka
        prokka --proteins --cpus "$threads" --force --outdir $PROKKA/"$name" --prefix "$name" "$i"
    fi
  done    

# prokka on the MAGs
for i in $MAGS
  do
    # extract name
    name=$(basename "$i" _mappedcontigs.sorted.fa)
    if [[ -f $PROKKA/$name/$name.gff ]]
      then
        echo "prokka already done for $name"
      else
        prokka --cpus "$threads" --force --metagenome --outdir $PROKKA/"$name" --prefix "$name" "$i"
    fi
  done

cp $PROKKA/*/*.gff $INDIR

# ROARY
# Options: -p INT    number of threads [1]
#          -o STR    clusters output filename [clustered_proteins]
#          -f STR    output directory [.]
#          -e        create a multiFASTA alignment of core genes using PRANK
#          -n        fast core gene alignment with MAFFT, use with -e
#          -i        minimum percentage identity for blastp [95]
#          -cd FLOAT percentage of isolates a gene must be in to be core [99]
#          -qc       generate QC report with Kraken
#          -k STR    path to Kraken database for QC, use with -qc
#          -a        check dependancies and print versions
#          -g INT    maximum number of clusters [50000]
#          -r        create R plots, requires R and ggplot2
#          -s        dont split paralogs
#          -ap       allow paralogs in core alignment
#          -z        dont delete intermediate files
#          -v        verbose output to STDOUT
#          -w        print version and exit
#          -y        add gene inference information to spreadsheet, doesnt work with -e
#          -iv STR   Change the MCL inflation value [1.5]
#          -h        this help message

# should do additional QC steps to make sure that taxonomy matches what is expected
# roary -qc -k $krakenDB $PROKKA/*.gff
# assess coverage here to?

# pangenome with roary
roary -s -p "$threads" -i 95 -cd 75 -f $OUTDIR -e -n $PROKKA/*/*.gff

rm -fr $INDIR

# FastTree –nt –gtr $OUTDIR*_*/core_gene_alignment.aln > $OUTDIR*_*/core_gene_tree.newick