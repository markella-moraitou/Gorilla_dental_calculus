#!/bin/bash -l

#SBATCH -A SNIC_PROJECT_ID
#SBATCH -p node
#SBATCH -n 1
#SBATCH -C mem256GB
#SBATCH -t 7:00:00
#SBATCH -J kraken2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=USER_EMAIL

# Microbial assignments using kraken for merged data
# Standard kraken database (RefSeq bacteria+archaea+viral+human)
# build by UPPMAX $KRAKEN_DB (updated started of each month)
# Adapted from /proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_SECONDSCREEN_190219/SCRIPTS/kraken2_190618.sh

# Merged data reads that couldn't be merged

module load bioinfo-tools Kraken2
echo $(module list)
echo "$SLURM_JOB_NAME"
echo $(readlink -f "$KRAKEN2_DEFAULT_DB")

DATADIR=8_humanHostFilt/unmapped
OUTDIR=T1_kraken2/new_kraken_version

# Set up database
 MY_DB_DIR=$SNIC_TMP/Kraken2
 MY_DB=$MY_DB_DIR/${KRAKEN2_DEFAULT_DB##*/}
 mkdir -p "$MY_DB"
 cp -av "$KRAKEN2_DEFAULT_DB"/* "$MY_DB"/

cd $DATADIR || exit

#Only for merged reads

find . -name "*_m_host_unmapped.fastq.gz" | while read i
do
  if [ -f $OUTDIR/"${i%host_unmapped.fastq.gz}"kraken2_output.txt ] 
  then
      continue
  else
    echo "Sample: ${i%_host_unmapped.fastq.gz}"
    kraken2 --db "$MY_DB" "$i" --threads "$SLURM_CPUS_ON_NODE" --report $OUTDIR/"${i%host_unmapped.fastq.gz}"kraken2_report.txt --report-zero-counts --output $OUTDIR/"${i%host_unmapped.fastq.gz}"kraken2_output.txt;
  fi
done
