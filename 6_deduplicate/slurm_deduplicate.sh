#!/bin/bash -l

#SBATCH -A snic2020-5-528
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 15:00:00
#SBATCH -J deduplicate
#SBATCH -C mem256GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Markella.Moraitou.0437@student.uu.se

# analysis step 5: Removing duplicates 
# with Tom's in house script
# Adapted from /proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_SECONDSCREEN_190219/SCRIPTS/rmdup_merged_190606.sh

# may need to increase memory depending on input file size

echo $SLURM_JOB_NAME
echo $(module list)

# define shortcuts

DATADIR=/proj/sllstore2017021/nobackup/MARKELLA/5_qualityFilt
OUTDIR=/proj/sllstore2017021/nobackup/MARKELLA/6_deduplicate
rmdup_script=/proj/sllstore2017021/nobackup/EXAMPLES_for_new_users/remove_duplicates_single_end.py

#run deduplication
cd $DATADIR

find . -name "*passed.fastq.gz" | while read i
do
    cd $OUTDIR
    python2.7 $rmdup_script $DATADIR/$i ${i%_passed.fastq.gz};
done

# output is not compressed, compress with gzip
cd $OUTDIR

for i in *.fastq;
do
    gzip $i;
done

