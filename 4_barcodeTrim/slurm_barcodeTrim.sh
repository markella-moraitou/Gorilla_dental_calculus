#!/bin/bash -l

#SBATCH -A SNIC_PROJECT_ID
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 15:00:00
#SBATCH -J barcodeTrim
#SBATCH --mail-type=ALL
#SBATCH --mail-user=USER_EMAIL

# analysis step 3: trimming barcodes and 3' adapters 
# adapted from: /proj/SNIC_PROJECT/nobackup/JAELLE/DENTAL_CALCULUS_SECONDSCREEN_190219/SCRIPTS/adrm_2_run1_190527.sh (merged)
# and: /proj/SNIC_PROJECT/nobackup/JAELLE/DENTAL_CALCULUS_SECONDSCREEN_190219/SCRIPTS/adrm_3_run1_190528.sh (unmerged)
# on merged and unmerged reads
# 3' adapters could be in the sequence if a read is very short and the sequencing continues past the insert

#This script shouldn't produce any output for the samples that have no barcodes (the ones left out of the if-statements).
#In the next processing step (5_qualityFilt), these samples should be taken straight from the previous step (3_adapterRem)

#load modules

module load bioinfo-tools AdapterRemoval bbmap pigz

echo "$SLURM_JOB_NAME"
echo $(module list)

#define shortcuts

OUTDIR=4_barcodeTrim
DATADIR=3_adapterRem

#Merging collapsed and collapsed.truncated output
cd $DATADIR || exit

find . -name '*collapsed.gz' -size +0 | while read i
do
    cat "$i" "${i%gz}"truncated.gz > $OUTDIR/"${i%.collapsed.gz}".fastq.gz;
done


#Removing barcodes from merged reads (that contain barcodes)
cd $OUTDIR || exit

find . -name '*_m.fastq.gz' -size +0 | while read i
do
    n=$(basename "$i")
    if [[ "${n:0:2}" == "G0" || "${n:0:2}" == "Gb" || "${n:0:2}" == "BL" || "${n:0:2}" == "BE" || "${n:0:2}" == "24" || "${n:0:2}" == "BS" || "${n:0:5}" == "ERR28" ]] #selecting the samples that contain barcodes, respectively DC2 (Uppsala), DC1 (Uppsala), library blanks, extraction blanks, ERR2868193
    then
        python2.7 trim_barcodes_henrique.py "$i" "${i%.fastq.gz}"
    fi
done

## I am not using the part of the script for unmerged reads, because I will only retain merged reads for downstream analysis

#Removing 5' barcodes from unmerged forward reads (that contain barcodes)
#cd $DATADIR
#find . -name '*pair1.truncated.gz' -size +0 | while read i
#    do
#    n=$(basename $i)
#    if [[ "${n:0:2}" == "G0" || "${n:0:2}" == "Gb" || "${n:0:2}" == "BL" || "${n:0:2}" == "BE" || "${n:0:2}" == "[0-9][0-9]" || "${n:0:2}" == "BS" || "${n:0:5}" == "ERR28" ]] #selecting the samples that contain barcodes, respectively DC2 and DC3(Uppsala), DC1 (Uppsala), library blanks, extraction blanks, ERR2868193
#    then
#        cd $OUTDIR
#        python2.7 trim_fwd_barcodes_henrique.py $DATADIR/$i ${i%.pair1.truncated.gz}_F;
#    fi
#done

cd $OUTDIR || exit

# output is not compressed, compress with gzip in parallel (pigz)
ls *fastq | while read i
do
    pigz -p "$SLURM_CPUS_ON_NODE" "$i"
done

# Remove 3' adapters+barcodes (from unmerged reads - if reads were long enough to be merged they don't have 3' barcodes/adapters)
# Use full P5 adapter sequence with barcodes
# Minlength of 30 bp (as now have removed barcodes so should be correct)
# Mismatch of 2 as Henrique previously tested was optimal of this step
# And rename files with fastq.gz extension for consistency
# Also do some read counts

#find . -name '*_F_bctrim.fastq.gz' | while read i #Out of the samples that got the 5' barcodes trimmed, I set the appropriate adapters (the adapters themselves have probably been removed at the adapter removal step, but the barcodes might be different so I am setting the adapters correctly to be sure)
#do
#    n=$(basename $i)
#    if [[ "${n:0:2}" == "G0" || "${n:0:2}" == "BL" || "${n:0:2}" == "BE" || "${n:0:2}" == "BS" ]] #DC2 samples
#    then
#        adapter_list=4_barcodeTrim/p7_p5_revcomp_dc2_dc3.txt
#    elif [[ "${n:0:2}" == "Gb" ]] #DC1 samples
#    then
#        adapter_list=4_barcodeTrim/p7_p5_revcomp_dc1.txt
#    elif [[ "${n:0:5}" == "ERR28" ]]
#    then
#        adapter_list=4_barcodeTrim/p7_p5_revcomp_ena.txt
#    fi
#    AdapterRemoval --file1 $i --basename ${i%.fastq.gz}_adrm --adapter-list $adapter_list --gzip --trimns --trimqualities --minquality 30 \
#    --minlength 30 --mm 2 --threads $SLURM_CPUS_ON_NODE
#    mv ${i%.fastq.gz}_adrm.truncated.gz ${i%.fastq.gz}_adrm.fastq.gz
#done
