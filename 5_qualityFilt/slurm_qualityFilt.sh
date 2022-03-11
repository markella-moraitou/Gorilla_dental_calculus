#!/bin/bash -l

#SBATCH -A SNIC_PROJECT_ID
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -J qualityFilt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=USER_EMAIL

# analysis step 4: Quality filtering with PrinSeq-Lite
# for collapsed reads as well as the forward reads that couldn't be merged
# Filtering outreads with mean base quality <30 and length <30bp
# Adapted from /proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_SECONDSCREEN_190219/SCRIPTS/prinseq_merged_190604.sh
# and /proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_SECONDSCREEN_190219/SCRIPTS/prinseq_unmerged_190604.sh
# Not all of the samples went through barcode trimming, so this script takes input from both 3_adapterRem and 4_barcodeTrim to process the total number of samples

echo "$SLURM_JOB_NAME"
echo $(module list)

#Defining shortcuts

DATADIR1=4_barcodeTrim
DATADIR2=3_adapterRem
OUTDIR=5_qualityFilt
PRINSEQ=software/prinseq-lite-0.20.4

cd $DATADIR1 || exit #To get the files that went through barcode trimming

#Merged reads
find . -name "*_m_bctrim.fastq.gz" | while read i
do
    echo "$i"
    zcat "$i" | perl $PRINSEQ/prinseq-lite.pl -fastq stdin -min_qual_mean 30 -out_good $OUTDIR/"${i%_bctrim.fastq.gz}"_passed -out_bad $OUTDIR/"${i%_bctrim.fastq.gz}"_failed -log $OUTDIR/log_prinseq_merged.log;
done

#Forward reads - not processing them for the purpose of this analysis
#find . -name "*_F_bctrim_adrm.fastq.gz" | while read i
#do
#    echo "$i"
#    zcat $i | perl $PRINSEQ/prinseq-lite.pl -fastq stdin -min_qual_mean 30 -out_good $OUTDIR/${i%_bctrim_adrm.fastq.gz}_passed -out_bad $OUTDIR/${i%_bctrim_adrm.fastq.gz}_failed -log $OUTDIR/log_prinseq_forward.log;
#done

#Staying in DATADIR1 even for the merged files that didn't go through barcode trimming, because the concatenated .collaped.gz and .collaped.truncated.gz outputs are located there

find . -name "*_m.fastq.gz" | while read i
do
    if [[ ! -f $OUTDIR/${i%.fastq.gz}_passed.fastq ]] #to make sure I don't overwrite the prinSeq output of barcode trimmed samples with the filtered but not trimmed output
    then
        echo "$i"
        zcat "$i" | perl $PRINSEQ/prinseq-lite.pl -fastq stdin -min_qual_mean 30 -out_good $OUTDIR/"${i%.fastq.gz}"_passed -out_bad $OUTDIR/"${i%.fastq.gz}"_failed -log $OUTDIR/log_prinseq_merged.log;
    fi
done

#cd $DATADIR2 #The forward reads are in the 3_adapterRem directory - not processing them for the purpose of this analysis

#find . -name "*_m.pair1.truncated.gz" | while read i
#do
#    if [[ ! -f $OUTDIR/${i%.pair1.truncated.gz}_F_passed.fastq ]] #to make sure I don't overwrite the prinSeq output of barcode trimmed samples with the filtered but not trimmed output
#    then
#        echo "$i"
#        zcat $i | perl $PRINSEQ/prinseq-lite.pl -fastq stdin -min_qual_mean 30 -out_good $OUTDIR/${i%.pair1.truncated.gz}_F_passed -out_bad $OUTDIR/${i%.pair1.truncated.gz}_F_failed -log $OUTDIR/log_prinseq_forward.log;
#    fi
#done

# output is not compressed, so compress with pigz in parallel

cd $OUTDIR || exit

find . -name "*.fastq" | while read i
do
    gzip "$i"
done
