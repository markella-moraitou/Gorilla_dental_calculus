#!/bin/bash -l

#SBATCH -A snic2020-5-528
#SBATCH -p devcore
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -J read_count
#SBATCH --mail-type BEGIN
#SBATCH --mail-type END
#SBATCH --mail-type FAIL
#SBATCH --mail-user Markella.Moraitou.0437@student.uu.se

#Counts the reads in the fastq files produced after decontamination with KrakenTools

echo $SLURM_JOB_NAME
echo $( date )

date=$(date +%y%m%d)

DIR=/proj/sllstore2017021/nobackup/MARKELLA/RD1_krakentools

cd $DIR

echo "sample, decontam_krakentools" > readcount_decontam_kt_${date}.txt

find . -name "*decontam1.fastq.gz" | while read i
do
    name=${i%_m_decontam1.fastq.gz}
    rc=$(( $(zcat $i | wc -l)/4 ))
    echo "$name, $rc">> readcount_decontam_kt_${date}.txt
done
