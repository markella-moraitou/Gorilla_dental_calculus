#!/bin/bash
#SBATCH -A SNIC_PROJECT
#SBATCH -p core
#SBATCH --ntasks 8 --mem-per-cpu=2G
#SBATCH -t 02:00:00
#SBATCH -M snowy
#SBATCH -J dietreads
#SBATCH --mail-type=FAIL
#SBATCH -o SNIC_PROJECT/ADRIAN/logs/slurm-%A.out

###### fastq decontamination using KrakenTools
# extract_kraken_reads.py usage/options
# python extract_kraken_reads.py ...
#     -k, --kraken MYFILE.KRAKEN.............Kraken output file
#     -s, -s1, -1, -U SEQUENCE.FILE..........FASTA/FASTQ sequence file (may be gzipped)
#     -s2, -2 SEQUENCE2.FILE.................FASTA/FASTQ sequence file (for paired reads, may be gzipped)
#     -o, --output2 OUTPUT.FASTA.............output FASTA/Q file with extracted seqs
#     -t, --taxid TID TID2 etc...............list of taxonomy IDs to extract (separated by spaces)
#
# Optional:
#     -o2, --output2 OUTPUT.FASTA.............second output FASTA/Q file with extracted seqs (for paired reads)
#     --fastq-output..........................Instead of producing FASTA files, print FASTQ files (requires FASTQ input)
#     --exclude...............................Instead of finding reads matching specified taxids, finds reads NOT matching specified taxids.
#     -r, --report MYFILE.KREPORT.............Kraken report file (required if specifying --include-children or --include-parents)
#     --include-children......................include reads classified at more specific levels than specified taxonomy ID levels.
#     --include-parents.......................include reads classified at all taxonomy levels between root and the specified taxonomy ID levels.
#     --max #.................................maximum number of reads to save.
#     --append................................if output file exists, appends reads
#     --noappend..............................[default] rewrites existing output file
######
# load modules
module load bioinfo-tools biopython/1.68-py3 bwa samtools mapDamage/2.0.9
echo $(module list)

ymd=$(date +%y%m%d)
echo $ymd
threads=8

# set variables
GENUS="Galium"
REF=SNIC_PROJECT/ADRIAN/genomes/$GENUS/*.fna
OUTDIR=SNIC_PROJECT/ADRIAN/D4_diet_damage/$GENUS
mkdir -p $OUTDIR

name="MTM009"
READS=SNIC_PROJECT/MARKELLA/D1_reads4Diet/${name}_m__bact_arch_vir_removed.fastq.gz
KRAKEN_OUTPUT=SNIC_PROJECT/MARKELLA/D2_kraken2_full_db/${name}_m_bact_arch_vir_removed_kraken2_output.txt
KRAKEN_REPORT=SNIC_PROJECT/MARKELLA/D2_kraken2_full_db/${name}_m_bact_arch_vir_removed_kraken2_report.txt
#taxa=$(grep $GENUS $KRAKEN_REPORT | awk '{print $5}' | sed 's/^/"/;s/$/"/' | tr '\n' ' ')
taxa=$(grep $GENUS $KRAKEN_REPORT | awk '{print $5}' | tr '\n' ' ')
echo $taxa

# subset fastq to genus reads
python SNIC_PROJECT/ADRIAN/bin/KrakenTools/extract_kraken_reads.py -k $KRAKEN_OUTPUT -s $READS -o ${OUTDIR}/${name}_$GENUS.fastq -t $taxa --fastq-output
pigz -f -p $threads ${OUTDIR}/${name}_$GENUS.fastq

# map to Galium reference genome
bwa index $REF
bwa mem -t $threads $REF ${OUTDIR}/${name}_$GENUS.fastq.gz | samtools view -hb -@ $threads -q 30 -F 4 -F 256 - | samtools sort -n -@ $threads -o $OUTDIR/${name}_$GENUS.bam -

# run mapdamage
cd $OUTDIR

# use forked version of mapDamage, for custom plots
module purge
module load R/4.1.1 R_packages/4.1.1
python3 SNIC_PROJECT/ADRIAN/bin/mapDamage/mapdamage/main.py -i $OUTDIR/${name}_$GENUS.bam -r $REF --merge-reference-sequences

