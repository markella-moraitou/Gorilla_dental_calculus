#!/bin/bash -l

#SBATCH -A snic2020-5-528
#SBATCH -p node
#SBATCH -C mem256GB
#SBATCH -N 1-1
#SBATCH -t 10:00:00
#SBATCH -J gtdbtk
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Markella.Moraitou.0437@student.uu.se

#Infer taxonomy of top quality bins with GTDB-TK
#Adapted from: /proj/sllstore2017021/nobackup/ADRIAN/scripts/mag/mag_bin_taxonomy_gtdbtk_210111.sh
module load bioinfo-tools prodigal hmmer pplacer FastTree gsl GTDB-Tk

echo $(module list)

OUTDIR=/proj/sllstore2017021/nobackup/MARKELLA/M4_MAGtaxonomy
bin_list=/proj/sllstore2017021/nobackup/MARKELLA/M3_binning/medium_quality_MAGs.txt

gtdbtk classify_wf --cpus $SLURM_CPUS_ON_NODE --batchfile $bin_list --out_dir $OUTDIR --prefix mq


bin_list=/proj/sllstore2017021/nobackup/MARKELLA/M3_binning/high_quality_MAGs.txt

gtdbtk classify_wf --cpus $SLURM_CPUS_ON_NODE --batchfile $bin_list --out_dir $OUTDIR --prefix hq
