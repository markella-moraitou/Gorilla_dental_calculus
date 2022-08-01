#!/bin/bash -l
#SBATCH -A SNIC_PROJECT_ID
#SBATCH -p core
#SBATCH --ntasks 8 --mem-per-cpu=1G
#SBATCH -t 00:10:00
#SBATCH -J nitrate-tblastn
#SBATCH --mail-type=FAIL
#SBATCH -o SNIC_PROJECT/logs/slurm-%A_%a.out

# confirm the presence of genes in a genome assembly using blast
#
# Dependencies:
#   blast
#
# Input:
#   <genome_assembly>
#   <blast_db>
#
# Output:
#   <genome_assembly>.blast_confirm.txt
#
#   <genome_assembly>.blast_confirm.txt
#     - tab-delimited file with the following columns:
#       - qseqid
#       - sseqid
#       - stitle
#       - pident
#       - qcovs
#       - length
#       - mismatch
#       - gapopen
#       - qstart
#       - qend
#       - sstart
#       - send
#       - qframe
#       - sframe
#       - frames
#       - evalue
#       - bitscore
#       - qseq
#       - sseq

# load modules
module load bioinfo-tools blast R R_packages

# set variables
threads=$SLURM_NTASKS
species_list="Rothia
Neisseria
Veillonella"

for SPECIES in $species_list
  do

  OUTDIR=SNIC_PROJECT/ADRIAN/calculus/Gorilla/mag/nitrate_genes/blast/$SPECIES
  mkdir -p "$OUTDIR"

  # multifasta of nitrate/molybdenum genes
  GOI=SNIC_PROJECT/ADRIAN/genes/$SPECIES/nitrate.faa
  
  GENOME_DIR=SNIC_PROJECT/ADRIAN/genomes/s__$SPECIES

  # clean up genome directory first
  # for i in $(find $GENOME_DIR/*.fna -type f)
  #   do
  #     name=$(basename $i .fna)
  #     mkdir -p $GENOME_DIR/fasta/$name
  #     mv $i* $GENOME_DIR/fasta/$name
  #   done

  # clean up MAG directory first
  # for i in $(find $MAG_DIR/*.fa -type f)
  #   do
  #     name=$(basename $i .fa)
  #     mkdir -p $MAG_DIR/$name
  #     mv $i* $MAG_DIR/$name
  #   done

  INDIR=SNIC_PROJECT/ADRIAN/genomes/s__$SPECIES/fasta
  REFS=$(find "$INDIR"/*/*.fna -type f)

  MAG_DIR=SNIC_PROJECT/ADRIAN/calculus/Gorilla/mag/metaWRAP/full-decontam-set/phylophlan/$SPECIES
  MAGS=$(find "$MAG_DIR"/*/*.fa -type f | grep -v -E "BL|BE|BS|ERR|EXB|LIB")


  # make a blast database for each of the genomes of interest
  for i in $REFS
    do
      # extract directory name
      # use sed to cut at last /
      name=$(dirname "$i" | sed 's/.*\///')
      # if file exists, skip
      if [[ -f $INDIR/$name/blastdb/$name.nhr ]]
        then
          echo "blastdb already done for $name"
        else
          # make blastdb
          makeblastdb -in "$i" -dbtype nucl -out "$INDIR"/"$name"/blastdb/"$name"

      fi
      # run blast
      tblastn -query "$GOI" -db "$INDIR"/"$name"/blastdb/"$name" -outfmt "6 qseqid sseqid stitle pident qcovs length mismatch gapopen qstart qend sstart send qframe sframe frames evalue bitscore qseq sseq" -task tblastn -num_threads "$threads" -out "$OUTDIR"/"${name}"_blast_confirm.txt
    done

  # do the same for MAGs
  for i in $MAGS
    do
      # extract directory name
      # use sed to cut at last /
      name=$(dirname "$i" | sed 's/.*\///')
      # if file exists, skip
      if [[ -f $MAG_DIR/$name/blastdb/$name.nhr ]]
        then
          echo "blastdb already done for $name"
        else
          # make blastdb
          makeblastdb -in "$i" -dbtype nucl -out "$MAG_DIR"/"$name"/blastdb/"$name"

      fi
      # run blast
      tblastn -query "$GOI" -db "$MAG_DIR"/"$name"/blastdb/"$name" -outfmt "6 qseqid sseqid stitle pident qcovs length mismatch gapopen qstart qend sstart send qframe sframe frames evalue bitscore qseq sseq" -task tblastn -num_threads "$threads" -out "$OUTDIR"/"${name}"_blast_confirm.txt
    done
done

cd SNIC_PROJECT/ADRIAN/calculus/Gorilla/mag/nitrate_genes/blast || exit
Rscript SNIC_PROJECT/ADRIAN/scripts/MARKELLA/mag/summarise-nitrate-blast.R
