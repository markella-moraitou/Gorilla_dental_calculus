#!/bin/bash -l

#SBATCH -A snic2020-5-528
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 10:00:00
#SBATCH -J metabat
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Markella.Moraitou.0437@student.uu.se

#MetaBat binning with minContig 1500
#CheckM to evaluate
#Adapted from /proj/sllstore2017021/nobackup/ADRIAN/scripts/mag/mag_binning_metabat2_210111.sh

module load bioinfo-tools MetaBat CheckM
echo $SLURM_JOB_NAME
echo $(module list)

ymd=$(date +%y%m%d)

ASDIR=/proj/sllstore2017021/nobackup/MARKELLA/M1_MAGassemblies
COVDIR=/proj/sllstore2017021/nobackup/MARKELLA/M2_readmapping
OUTDIR=/proj/sllstore2017021/nobackup/MARKELLA/M3_binning
CHKMDIR=$OUTDIR/CheckM_results

cd $COVDIR

ls *sorted.bam | while read i;
do
   echo $i
   jgi_summarize_bam_contig_depths --outputDepth $OUTDIR/${i%.sorted.bam}_depth.txt $i
   metabat2 -i $ASDIR/${i%_mappedcontigs.sorted.bam}/${i%_mappedcontigs.sorted.bam}.contigs.fa -a $OUTDIR/${i%.sorted.bam}_depth.txt -o $OUTDIR/${i%_mappedcontigs.sorted.bam} --minContig 1500;
done

checkm lineage_wf -t 10 -x fa --tab_table -f $CHKMDIR/all_checkm_results.tab $OUTDIR $CHKMDIR 

#Get list of genomes with high and medium quality
cd $OUTDIR
bash get_qual_genomes.sh