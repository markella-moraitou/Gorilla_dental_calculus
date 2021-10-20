#!/bin/bash -l

#SBATCH -A snic2020-5-528
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 5:00:00
#SBATCH -J mapDamage_mtDNA
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Markella.Moraitou.0437@student.uu.se

#Calculate damage patterns of host mtDNA using mapDamage(to compare with MAG damage patterns)

module load bioinfo-tools bbmap mapDamage samtools BEDTools

echo $SLURM_JOB_NAME
echo $(module list)

#BAM files to be used

readdir=/proj/sllstore2017021/nobackup/MARKELLA/8_humanHostFilt/mapped
outdir=/proj/sllstore2017021/nobackup/MARKELLA/M7_mapdamage/mtDNA
hq_mags=/proj/sllstore2017021/nobackup/MARKELLA/M3_binning/high_quality_MAGs.txt
reference=$outdir/gorilla_gorilla_chrMT.fa

cd $readdir

ls *_m_host_mapped_sorted.bam | while read i
do 
    #check if sample has produced any MAGs
    if grep -q $(basename ${i%_m_host_mapped_sorted.bam}) $hq_mags
    then 
        cd $outdir
        #Map using bbmap - after piping from bamtofastq
        bedtools bamtofastq -i $readdir/$i -fq /dev/stdout | bbwrap.sh ref=$reference int=t in=stdin.fq mapper=bbmap out=${i%m_host_mapped_sorted.bam}mtDNA.sam  maxindel=80 pigz=t unpigz=t nodisk
        #Save mapped reads into BAM file
        samtools view -Sb -F 4 -@ $SLURM_CPUS_ON_NODE ${i%m_host_mapped_sorted.bam}mtDNA.sam -o ${i%m_host_mapped_sorted.bam}mtDNA.bam
        #Run mapDamage
        mapDamage -i ${i%m_host_mapped_sorted.bam}mtDNA.bam -r $reference
    fi
done
