#!/bin/bash -l

#SBATCH -A SNIC_PROJECT_ID
#SBATCH -p core
#SBATCH -n 5
#SBATCH -t 20:00:00
#SBATCH -J phixRem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=USER_EMAIL

# analysis step 6: Map to PhiX and remove mapped reads
# Adapted from SNIC_PROJECT/JAELLE/DENTAL_CALCULUS_SECONDSCREEN_190219/SCRIPTS/rmphix_merged_190606.sh

#load modules
module load bioinfo-tools samtools bwa BEDTools

echo "$SLURM_JOB_NAME"
echo $(module list)

#define shortcuts

INDEXDIR=SNIC_PROJECT/JAELLE/REFERENCES/phix
DATADIR=6_deduplicate
OUTDIR=7_phixRem

cd $DATADIR || exit

#Align every fastq file against PhiX reference using BWA MEM
#Then output only the unmapped reads as BAM with suffix rmphix.bam
#Then convert this bam file to fastq
#Finally compress
find . -name "*fastq.gz" | while read i
do
    bwa mem -t "$SLURM_CPUS_ON_NODE" $INDEXDIR/phix "$i" | samtools view -Sb -f 4 - > $OUTDIR/"${i%dedup.fastq.gz}"rmphix.bam
    bedtools bamtofastq -i $OUTDIR/"${i%dedup.fastq.gz}"rmphix.bam -fq $OUTDIR/"${i%dedup.fastq.gz}"rmphix.fastq
    pigz -p "$SLURM_CPUS_ON_NODE" $OUTDIR/"${i%dedup.fastq.gz}"rmphix.fastq;
done
