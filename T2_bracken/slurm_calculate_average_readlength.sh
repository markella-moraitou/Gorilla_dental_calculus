#!/bin/bash -l

#SBATCH -A SNIC_PROJECT_ID
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 20:00:00
#SBATCH -J read_length

#Script for calculating average read length for the sequences used as input for Kraken2
#For merged reads only
#Adapted from SNIC_PROJECT/JAELLE/DENTAL_CALCULUS_SECONDSCREEN_190219/SCRIPTS/calculate_average_readlength_190731.sh

ymd=$(date +%y%m%d)

DATADIR=8_humanHostFilt/unmapped
outfile=T2_bracken/log_readlength_hostfilt_${ymd}.txt

echo "#sample, ave read length" > $outfile

cd $DATADIR

find . -name "*m_host_unmapped.fastq.gz" | while read i
do
   c=$(zcat $i | awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' -)
   echo "${i%.gz}, $c" >> $outfile
done

#Calculate average read length
echo "Average read length:"

count=0; total=0; for i in $( awk '{ print $2; }' $outfile );\
do total=$(echo $total+$i | bc ); \
((count++)); done; echo "scale=2; $total / $count" | bc