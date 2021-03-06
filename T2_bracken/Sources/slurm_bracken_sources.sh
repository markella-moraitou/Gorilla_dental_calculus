#!/bin/bash -l

#SBATCH -A SNIC_PROJECT_ID
#SBATCH -p node
#SBATCH -N 1-1
#SBATCH -t 04:00:00
#SBATCH -J bracken
#SBATCH --mail-type=ALL
#SBATCH --mail-user=USER_EMAIL

# For sources to be used in FEAST
# Building Bracken db file from Kraken2 default database
# For average read length of 55 bp -- the average of the mean read length of every sample in the dataset (see: log_readlength_hostfilt_05022021.txt)
# Using kmer length of 35, as recommended by Bracken for Kraken2
# Adapted from /proj/SNIC_PROJECT/nobackup/JAELLE/DENTAL_CALCULUS_SECONDSCREEN_190219/SCRIPTS/bracken_unmerged_190801.sh

# Bracken on merged samples
# species level (-l S)
# no threshold i.e. 0 reads required before reestimation (-t 0)
# will apply abundance threshold independently

module load bioinfo-tools Kraken2
echo $(module list)
echo "$SLURM_JOB_NAME"

echo "Kraken default database version:"
echo $(readlink -f "$KRAKEN2_DEFAULT_DB")

DATADIR=T1_kraken2/Sources
OUTDIR=T2_bracken/Sources

# Set up database
 MY_DB_DIR=$SNIC_TMP/Kraken2
 MY_DB=$MY_DB_DIR/${KRAKEN2_DEFAULT_DB##*/}
 mkdir -p "$MY_DB"
 cp -av "$KRAKEN2_DEFAULT_DB"/* "$MY_DB"/

cd $DATADIR || exit
echo "Bracken-build"
#Bracken-build step
~/bin/Bracken/bracken-build -d "$MY_DB" -t "$SLURM_CPUS_ON_NODE" -k 35 -l 55

#Run Bracken for merged reads, species level
find . -name "*report.txt" -maxdepth 1 | while read i
do
   ~/bin/Bracken/bracken -d "$MY_DB" -i "$i" -o $OUTDIR/"${i%report.txt}"s_bracken.txt -r 55 -l S -t 0;
done

mv *report_bracken_species.txt $OUTDIR
