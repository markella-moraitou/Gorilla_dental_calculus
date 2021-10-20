#!/bin/bash -l

#SBATCH -A snic2020-5-528
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 3:00:00
#SBATCH -J reads_for_diet_profiling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Markella.Moraitou.0437@student.uu.se

# Preprocessing for diet profiling 
# This script uses Kraken2 output to extract the reads that have been classified as bacterial, viral or archeal
# Only processes blanks and environmental controls

module load bioinfo-tools pigz
echo $SLURM_JOB_NAME
echo $(module list)

DATADIR=/proj/sllstore2017021/nobackup/MARKELLA/T1_kraken2
SEQDIR=/proj/sllstore2017021/nobackup/MARKELLA/8_humanHostFilt/unmapped
OUTDIR=/proj/sllstore2017021/nobackup/MARKELLA/D1_reads4Diet
SOFTWARE=/proj/sllstore2017021/nobackup/MARKELLA/software/KrakenTools

#Activate python environment
source /proj/sllstore2017021/nobackup/MARKELLA/RD1_krakentools/krakentools/bin/activate

#Remove reads assigned to species in exogenous_id_list.txt according to Kraken2 outputs
cd $DATADIR
ls [BE][ELRS]*output.txt | while read i
do
    if [[ -f $OUTDIR/${i%kraken2_output.txt}_bact_arch_vir_removed.fastq.gz ]]; then #Check if output already exists
        continue
    else
    #TaxIDs to be excluded: 2 (Bacteria), 2157 (Archaea), 10239 (Viruses)
        $SOFTWARE/extract_kraken_reads.py -k $i -s $SEQDIR/${i%kraken2_output.txt}host_unmapped.fastq.gz  --report ${i%_output.txt}_report.txt \
        -o $OUTDIR/${i%kraken2_output.txt}_bact_arch_vir_removed.fastq --fastq-output --exclude -t 2 2157 10239 --include-children
    fi
done

ls [EL][XI]*output.txt | while read i
do
    if [[ -f $OUTDIR/${i%kraken2_output.txt}_bact_arch_vir_removed.fastq.gz ]]; then #Check if output already exists
        continue
    else
    #TaxIDs to be excluded: 2 (Bacteria), 2157 (Archaea), 10239 (Viruses)
        $SOFTWARE/extract_kraken_reads.py -k $i -s $SEQDIR/${i%kraken2_output.txt}host_unmapped.fastq.gz  --report ${i%_output.txt}_report.txt \
        -o $OUTDIR/${i%kraken2_output.txt}_bact_arch_vir_removed.fastq --fastq-output --exclude -t 2 2157 10239 --include-children
    fi
done

#deactivate python env
deactivate

#compress output files
cd $OUTDIR
pigz -p $SLURM_CPUS_ON_NODE *_bact_arch_vir_removed.fastq
