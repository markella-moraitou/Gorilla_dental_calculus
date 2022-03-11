#!/bin/bash -l

#SBATCH -A SNIC_PROJECT_ID
#SBATCH -p node
#SBATCH -N 1-1
#SBATCH -t 3-00:00:00
#SBATCH -J HUMAnN2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=USER_EMAIL

# Humann2 pipeline with Metaphlan2 for tax profiling and mapping/DIAMOND for functional assignments
# Using contam-filtered reads based on Kraken2 assignments (../M1_decontam/)
# Adapted from /proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_FIRSTSCREEN_181220/SCRIPTS/func_humann2_190415.sh

#Activate conda environment

conda activate /home/markmora/.conda/envs/humann2

#error loading humann2 straight, had to load metaphlan2 first
module load bioinfo-tools metaphlan2 biopython humann2

echo "$SLURM_JOB_NAME"
echo $(module list)

DATADIR=RD2_mapping/decontam_files
OUTDIR=F1_HUMAnN2

#Install database
#cd /home/markmora/.conda/envs/humann2/bin/
#/home/markmora/.conda/envs/humann2/bin/metaphlan2.py --install

cd $DATADIR || exit
ls *m_decontam.fastq.gz | while read i 
do
    if [[ -f $OUTDIR/${i%_m_decontam.fastq.gz}_genefamilies.tsv ]]
    then
        continue
    else
        humann2 --nucleotide-database /sw/apps/bioinfo/humann2/data/chocophlan/ \
        --protein-database /sw/apps/bioinfo/humann2/data/uniref \
        --metaphlan-options "--bowtie2db /home/markmora/.conda/envs/humann2/bin/metaphlan_databases" \
        --input "$i" \
        --output $OUTDIR \
        --output-basename "${i%_m_decontam.fastq.gz}" \
        --threads "$SLURM_CPUS_ON_NODE"
    fi
done

conda deactivate
