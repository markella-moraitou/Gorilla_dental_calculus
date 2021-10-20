#!/bin/bash -l

#SBATCH -A snic2020-5-528
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5:00:00
#SBATCH -J read_count
#SBATCH --mail-type BEGIN
#SBATCH --mail-type END
#SBATCH --mail-type FAIL
#SBATCH --mail-user Markella.Moraitou.0437@student.uu.se

#Counts the reads in the fastq files produced after AdapterRemoval, after singleton removal, after barcode trimming (applicable only for Uppsala samples) and after merging.

echo $SLURM_JOB_NAME
echo $( date )

date=$(date +%y%m%d)

echo "sample, adapterRem_m, barcodeTrim_m, qualityFilter_m, deduplicate_m, phiXrem_m, hostHumFilt_m" > readcount_${date}.txt

GTRIMDIR=/proj/sllstore2017021/nobackup/MARKELLA/2_polyGrem
BCTRIMDIR=/proj/sllstore2017021/nobackup/MARKELLA/4_barcodeTrim
QUALFILDIR=/proj/sllstore2017021/nobackup/MARKELLA/5_qualityFilt
DEDUPDIR=/proj/sllstore2017021/nobackup/MARKELLA/6_deduplicate
PHIXDIR=/proj/sllstore2017021/nobackup/MARKELLA/7_phixRem
HOSTHUMDIR=/proj/sllstore2017021/nobackup/MARKELLA/8_humanHostFilt/unmapped

cd $BCTRIMDIR

find . -name "*_m.fastq.gz" | while read i
do
    name=${i%_m.fastq.gz}
    gtrim_um
    adapterRem_m=$(( $(zcat $i | wc -l)/4 ))    
    if [[ -f ${name}_m_bctrim.fastq.gz ]] #there are no *bctrim.fastq.gz files for the Jena samples, so print the same count as adapteRem
    then
        barcodeTrim_m=$(( $(zcat ${name}_m_bctrim.fastq.gz | wc -l)/4 ))
    else
        barcodeTrim_m=$adapterRem_m
    fi
    
    qualityFilter_m=$(( $(zcat $QUALFILDIR/${name}_m_passed.fastq.gz | wc -l)/4 ))

    deduplicate_m=$(( $(zcat $DEDUPDIR/${name}_m_dedup.fastq.gz | wc -l)/4 ))

    phiXrem_m=$(( $(zcat $PHIXDIR/${name}_m_rmphix.fastq.gz | wc -l)/4 ))

    hostHumFilt_m=$(( $(zcat $HOSTHUMDIR/${name}_m_host_unmapped.fastq.gz | wc -l)/4 ))
    echo "$name, $adapterRem_m, $barcodeTrim_m, $qualityFilter_m, $deduplicate_m, $phiXrem_m, $hostHumFilt_m" >> ../readcount_${date}.txt
done
