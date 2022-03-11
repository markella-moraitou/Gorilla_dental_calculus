#!/bin/bash -l

#SBATCH -A SNIC_PROJECT_ID
#SBATCH -p node
#SBATCH -n 1
#SBATCH -C mem512GB
#SBATCH -t 40:00:00
#SBATCH -J kraken2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=USER_EMAIL

# Microbial assignments using kraken for merged data
# Using the NCBI nt database
# build by UPPMAX $KRAKEN_DB (updated monthly)
# Adapted from /proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_SECONDSCREEN_190219/SCRIPTS/kraken2_190618.sh

# Merged data reads that couldn't be merged

module load bioinfo-tools Kraken2
echo $(module list)
echo "$SLURM_JOB_NAME"

KRAKEN2_NT_DB=/sw/data/uppnex/Kraken2_data/latest_nt/
echo $(readlink -f $KRAKEN2_NT_DB)

DATADIR=D1_reads4Diet
OUTDIR=D2_kraken2_full_db

# Set up database - I am using the NCBI nt database because I need to detect eukaryotic taxa 
 MY_DB_DIR=$SNIC_TMP/Kraken2
 MY_DB=$MY_DB_DIR/${KRAKEN2_NT_DB##*/}
 mkdir -p "$MY_DB"
 cp -av $KRAKEN2_NT_DB/* "$MY_DB"/

cd $DATADIR || exit

#Only for merged reads

find . -name "*__bact_arch_vir_removed.fastq.gz" | while read i
do
    if [[ -f $OUTDIR/${i%_bact_arch_vir_removed.fastq.gz}bact_arch_vir_removed_kraken2_report.txt ]]
    then
      continue
    else
      echo "Sample: ${i%__bact_arch_vir_removed.fastq.gz}"
      kraken2 --db "$MY_DB" "$i" --threads 20 --report $OUTDIR/"${i%_bact_arch_vir_removed.fastq.gz}"bact_arch_vir_removed_kraken2_report.txt --report-zero-counts --output $OUTDIR/"${i%_bact_arch_vir_removed.fastq.gz}"bact_arch_vir_removed_kraken2_output.txt;
    fi
done
