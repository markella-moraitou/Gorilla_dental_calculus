#!/bin/bash -l

#SBATCH -A snic2020-5-528
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 24:00:00
#SBATCH -J decontam_reads
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Markella.Moraitou.0437@student.uu.se

#Removing contamination from reads - Step 1
#Copy of script for processing blanks
#Preprocessing for MAGs and functional profiling
#This script uses Kraken2 output to extract the reads that have been identified as exogenous during community-level analysis

module load bioinfo-tools pigz
echo $SLURM_JOB_NAME
echo $(module list)

DATADIR=/proj/sllstore2017021/nobackup/MARKELLA/T1_kraken2
SEQDIR=/proj/sllstore2017021/nobackup/MARKELLA/8_humanHostFilt/unmapped
OUTDIR=/proj/sllstore2017021/nobackup/MARKELLA/RD1_krakentools
LISTDIR=/proj/sllstore2017021/nobackup/MARKELLA/RD2_mapping
SOFTWARE=/proj/sllstore2017021/nobackup/MARKELLA/software/KrakenTools

#Activate python environment
source $OUTDIR/krakentools/bin/activate

#Remove reads assigned to species in exogenous_id_list.txt according to Kraken2 outputs
cd $DATADIR
ls [EB][ELRS]*output.txt | while read i
do
    if [[ -f $OUTDIR/${i%kraken2_output.txt}decontam1.fastq.gz ]]; then #Check if output already exists
        continue
    else
        $SOFTWARE/extract_kraken_reads.py -k $i -s $SEQDIR/${i%kraken2_output.txt}host_unmapped.fastq.gz \
        -o $OUTDIR/${i%kraken2_output.txt}decontam1.fastq --fastq-output --exclude -t $(cat $LISTDIR/exogenous_id_list.txt)
    fi
done

ls [EL][XI]*output.txt | while read i
do
    if [[ -f $OUTDIR/${i%kraken2_output.txt}decontam1.fastq.gz ]]; then #Check if output already exists
        continue
    else
        $SOFTWARE/extract_kraken_reads.py -k $i -s $SEQDIR/${i%kraken2_output.txt}host_unmapped.fastq.gz \
        -o $OUTDIR/${i%kraken2_output.txt}decontam1.fastq --fastq-output --exclude -t  $(cat $LISTDIR/exogenous_id_list.txt)
    fi
done


#Deactivate python environment
deactivate

#compress output files
cd $OUTDIR
pigz -p $SLURM_CPUS_ON_NODE *decontam1.fastq
