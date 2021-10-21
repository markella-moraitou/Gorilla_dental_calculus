#!/bin/bash -l

#SBATCH -A snic2020-5-528
#SBATCH -p node
#SBATCH -n 1
#SBATCH -C mem256GB
#SBATCH -t 10:00:00
#SBATCH -J kraken2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Markella.Moraitou.0437@student.uu.se

# Microbial assignments using kraken for file to be used as sources in Feast
# Standard kraken database (RefSeq bacteria+archaea+viral+human)
# build by UPPMAX $KRAKEN_DB (updated started of each month)
# Adapted from /proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_SECONDSCREEN_190219/SCRIPTS/kraken2_190618.sh

# Merged data reads that couldn't be merged

module load bioinfo-tools Kraken2
echo $(module list)
echo $SLURM_JOB_NAME
echo $(readlink -f $KRAKEN2_DEFAULT_DB)

DATADIR=/proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_FIRSTSCREEN_181220/SOURCES_P6_18119
OUTDIR=/proj/sllstore2017021/nobackup/MARKELLA/T1_kraken2/Sources

# Set up database
 MY_DB_DIR=$SNIC_TMP/Kraken2
 MY_DB=$MY_DB_DIR/${KRAKEN2_DEFAULT_DB##*/}
 mkdir -p $MY_DB
 cp -av $KRAKEN2_DEFAULT_DB/* $MY_DB/

# Database version
echo "Kraken default database version:"
readlink $KRAKEN2_DEFAULT_DB

cd $DATADIR

#Only for merged reads

find . -name "*.gz" | while read i
do
    echo "Sample: ${i%.*gz}"
    kraken2 --db $MY_DB $i --threads 20 --report $OUTDIR/${i%.*gz}_kraken2_report.txt --report-zero-counts --output $OUTDIR/${i%.gz}_kraken2_output.txt;
done
