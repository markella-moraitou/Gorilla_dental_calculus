#!/bin/bash -l

#SBATCH -A SNIC_PROJECT_ID
#SBATCH -p devcore
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -J read_count
#SBATCH --mail-type BEGIN
#SBATCH --mail-type END
#SBATCH --mail-type FAIL
#SBATCH --mail-user USER_EMAIL

#Counts the reads in the fastq files produced after decontamination with KrakenTools

echo "$SLURM_JOB_NAME"
echo $( date )

date=$(date +%y%m%d)

DIR=RD1_krakentools

cd $DIR || exit

echo "sample, decontam_krakentools" > readcount_decontam_kt_"${date}".txt

find . -name "*decontam1.fastq.gz" | while read i
do
    name=${i%_m_decontam1.fastq.gz}
    rc=$(( $(zcat "$i" | wc -l)/4 ))
    echo "$name, $rc">> readcount_decontam_kt_"${date}".txt
done
