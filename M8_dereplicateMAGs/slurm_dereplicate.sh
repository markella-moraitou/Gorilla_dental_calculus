#!/bin/bash -l

#SBATCH -A snic2020-5-528
#SBATCH -p devcore
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -J dRep
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Markella.Moraitou.0437@student.uu.se

#Dereplicate MQ and HQ MAGs with ANI larger than 95%

module load python bioinfo-tools mash prodigal CheckM centrifuge FastANI/1.33

echo $SLURM_JOB_NAME
echo $(module list)

outdir=/proj/sllstore2017021/nobackup/MARKELLA/M8_dereplicateMAGs
datadir=/proj/sllstore2017021/nobackup/MARKELLA/M3_binning
dRep=/home/markmora/.local/bin/dRep 

conda activate /home/markmora/.conda/envs/dRep

cd $datadir

#Get CheckM contamination and completeness info
echo "genome,completeness,contamination" > $outdir/genome_info.txt
awk -v FS="\t" 'BEGIN{OFS=","}; NR>1 { print $1 ".fa", $12, $13}' $datadir/CheckM_results/all_checkm_results.tab >> $outdir/genome_info.txt

#List of medium quality and high quality MAGs
$dRep dereplicate $outdir -g $(cat high_quality_MAGs.txt medium_quality_MAGs.txt | awk '{print $1}') --genomeInfo $outdir/genome_info.txt -sa 0.95

