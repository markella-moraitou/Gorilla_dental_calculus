#!/bin/bash -l

#SBATCH -A snic2020-5-528
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 3:00:00
#SBATCH -J index_combined_refs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Markella.Moraitou.0437@student.uu.se

#This script concatenates all the contaminant genomes together (contam_references)
#and all the non contaminant genomes together (non_contam_references)
#and then indexes the combined references
#BWA index requires 5.37N memory, where N is the size of the database

module load bioinfo-tools bwa
echo $SLURM_JOB_NAME
echo $(module list)

REFDIR=/proj/sllstore2017021/nobackup/MARKELLA/RD2_mapping/refs

###First, create the contaminant combined reference
cd $REFDIR/contam_references

cat **/**[!m]_genomic.fna.gz > contam_references.fna.gz

#Add identifier in sequence headers
gunzip contam_references.fna.gz
sed -i 's/>/>contam_/g' contam_references.fna
gzip contam_references.fna

###Second, create the non-contaminant combined reference
cd $REFDIR/noncontam_references

cat **/**[!m]_genomic.fna.gz > noncontam_references.fna.gz

#Add identifier in sequence headers
gunzip noncontam_references.fna.gz
sed -i 's/>/>endogenous_/g' noncontam_references.fna
gzip noncontam_references.fna

###Merge the two combined references into one
cd $REFDIR

cat contam_references/contam_references.fna.gz noncontam_references/noncontam_references.fna.gz > gorilla_dencalc_reference.fna.gz

#Index the final reference
bwa index -p gorilla_dencalc_reference -a bwtsw gorilla_dencalc_reference.fna.gz
