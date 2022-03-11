#!/bin/bash
#SBATCH -A SNIC_PROJECT_ID
#SBATCH -p node
#SBATCH -n 1
#SBATCH -C mem256GB
#SBATCH -t 06:00:00
#SBATCH -J gtdbtk
#SBATCH -o /proj/SNIC_PROJECT/nobackup/ADRIAN/logs/slurm-%A_%a.out
#SBATCH --mail-type=FAIL

module load bioinfo-tools prodigal hmmer pplacer FastTree gsl GTDB-Tk/0.3.2
echo $(module list)

ymd=$(date +%y%m%d)
echo "$ymd"

threads=$SLURM_CPUS_ON_NODE
BASEDIR=/home/adrianf/project-folder/nobackup/ADRIAN/calculus/Gorilla/mag/metaWRAP/full-decontam-set

BIN_REASSEMBLY=${BASEDIR}/MAG_binreassembly_metaWRAP
DATADIR=$BIN_REASSEMBLY/reassembled_bins
OUTDIR=${BASEDIR}/MAG_binclassification_gtdbtk
rm -fr $OUTDIR
mkdir -p $OUTDIR

# the batch file HAS TO CONTAIN EXACTLY 2 COLUMNS
# the output directory CANNOT HAVE ANYTHING IN IT
find $DATADIR/*.fa | awk '{print $1}' |  sed 's!.*/!!' | sed 's/\.[^\\.]*$//' > $OUTDIR/sample-names.txt
find $DATADIR/*.fa > $OUTDIR/file-paths.txt
paste $OUTDIR/file-paths.txt $OUTDIR/sample-names.txt > $OUTDIR/bin_list.txt


# Infer taxonomy of top quality bins with GTDB-TK
gtdbtk classify_wf --force --cpus "$threads" --batchfile $OUTDIR/bin_list.txt  --out_dir $OUTDIR


