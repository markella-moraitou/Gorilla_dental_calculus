#!/bin/bash -l

#SBATCH -A snic2020-5-528
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 60:00:00
#SBATCH -J get_accession_nums
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Markella.Moraitou.0437@student.uu.se

#This script parses the assembly summary files for bacteria and archaea(ftp://ftp.ncbi.nih.gov/genomes/genbank/bacteria/assembly_summary.txt)
#and saves the species name for each assembly accession number

echo $SLURM_JOB_NAME
echo $(module list)

summary_file=/proj/sllstore2017021/nobackup/MARKELLA/RD2_mapping/assembly_summary_bacteria.txt

touch species_names.txt

#Keep count so that the number of jobs running in parallel does not exceed the threshold

count=0
threshold=50

cat mag_tree_accnum.txt | while read j
do
    i=$(echo $j | tr -d '\r')
    #Check if the accession number already exists
    grep -q $i species_names.txt && continue
    #Check if either the accession number ($1) or the corresponding RefSeq accession number ($18) match, while incrementing the count variable
    count=$(expr $count + 1)
    awk -vLOOKUPVAL=$i -v OFS='\t' '
    BEGIN {FS="\t"}; $1 == LOOKUPVAL || $18 == LOOKUPVAL {print LOOKUPVAL, $8}' < $summary_file >> species_names.txt &
    #Wait once more than 50 jobs are running
    echo "jobs $(jobs -p | grep "^[0-9]" | wc -l)"
    if [[ $(jobs -p | grep "^[0-9]" | wc -l) -gt $threshold ]]
    then
        echo "waiting"
        wait $(jobs -p | head -1)
    fi
done
