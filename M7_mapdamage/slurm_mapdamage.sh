#!/bin/bash -l

#SBATCH -A snic2020-5-528
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -J mapDamage
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Markella.Moraitou.0437@student.uu.se

#Calculate damage patterns of HQ and MQ MAGs using mapDamage

module load bioinfo-tools bbmap mapDamage samtools

echo $SLURM_JOB_NAME
echo $(module list)

#BAM files to be used
MAGdir=/proj/sllstore2017021/nobackup/MARKELLA/M3_binning
readdir=/proj/sllstore2017021/nobackup/MARKELLA/RD2_mapping/decontam_files
outdir=/proj/sllstore2017021/nobackup/MARKELLA/M7_mapdamage


awk '{print $2}' $MAGdir/high_quality_MAGs.txt | while read i
do
    #Map using bbmap
    bbwrap.sh ref=$MAGdir/${i}.fa in=$readdir/${i%.*}_m_decontam.fastq.gz mapper=bbmap out=$outdir/${i}.sam  maxindel=80 pigz=t unpigz=t nodisk
    #Save mapped reads into BAM file
    samtools view -Sb -F 4 -@ $SLURM_CPUS_ON_NODE $outdir/${i}.sam -o $outdir/${i}.bam
    cd $outdir
    #Run mapDamage
    mapDamage -i ${i}.bam -r $MAGdir/${i}.fa
done

awk '{print $2}' $MAGdir/medium_quality_MAGs.txt | while read i
do
    #Map using bbmap
    bbwrap.sh ref=$MAGdir/${i}.fa in=$readdir/${i%.*}_m_decontam.fastq.gz mapper=bbmap out=$outdir/${i}.sam  ma$
    #Save mapped reads into BAM file
    samtools view -Sb -F 4 -@ $SLURM_CPUS_ON_NODE $outdir/${i}.sam -o $outdir/${i}.bam
    cd $outdir
    #Run mapDamage
    mapDamage -i ${i}.bam -r $MAGdir/${i}.fa
done

