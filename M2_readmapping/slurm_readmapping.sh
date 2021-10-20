#!/bin/bash -l

#SBATCH -A snic2020-5-528
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 12:00:00
#SBATCH -J mag_readmapping
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Markella.Moraitou.0437@student.uu.se

#Map reads back to Megahit assemblies and get some mapping stats
#Using bbmap instead of bwa since need to build a ref for bwa
#Adapted from /proj/sllstore2017021/nobackup/ADRIAN/scripts/mag/mag_assembly_readmap_210111.sh


module load bioinfo-tools samtools bbmap
echo $SLURM_JOB_NAME
echo $(module list)

ymd=$(date +%y%m%d)

#Read dir
RDIR=/proj/sllstore2017021/nobackup/MARKELLA/RD2_mapping/decontam_files
#Contig dir
CDIR=/proj/sllstore2017021/nobackup/MARKELLA/M1_MAGassemblies
OUTDIR=/proj/sllstore2017021/nobackup/MARKELLA/M2_readmapping

cd $RDIR

find . -name "*fastq.gz" | while read i;
do
    echo ${i%_m_decontam.fastq.gz}
    #align reads to contigs with bbmap, generate coverage report
    bbwrap.sh ref=$CDIR/${i%_m_decontam.fastq.gz}/${i%_m_decontam.fastq.gz}.contigs.fa in=$i mapper=bbmap out=$OUTDIR/${i%_m_decontam.fastq.gz}_mappedcontigs.sam maxindel=80 pigz=t unpigz=t nodisk
    #view all mapped reads | sort ready for binning
    samtools view -Sb -F 4 -@ $SLURM_CPUS_ON_NODE $OUTDIR/${i%_m_decontam.fastq.gz}_mappedcontigs.sam | samtools sort -@ $SLURM_CPUS_ON_NODE -o $OUTDIR/${i%_m_decontam.fastq.gz}_mappedcontigs.sorted.bam
    ct=$(samtools view -c -@ $SLURM_CPUS_ON_NODE $OUTDIR/${i%_m_decontam.fastq.gz}_mappedcontigs.sorted.bam)
    echo "${i%.fastq.gz}, $ct" >> $OUTDIR/log_map_readct_${ymd}.txt;
done
