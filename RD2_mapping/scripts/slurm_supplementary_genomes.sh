#!/bin/bash -l

#SBATCH -A snic2020-5-528
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J get_suppl_genomes
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Markella.Moraitou.0437@student.uu.se

#Second round of selecting reference genomes - for species that had no available genomes
#A list with sister species is provided, all the available genomes are listed out of which one genome per original species is selected.
#Modified from slurm_get_contam_genomes.sh and slurm_get_noncontam_genomes.sh

module load R_packages

echo $SLURM_JOB_NAME
echo $(module list)

REFDIR=/proj/sllstore2017021/nobackup/MARKELLA/RD2_mapping/refs
MAINDIR=/proj/sllstore2017021/nobackup/MARKELLA/RD2_mapping

cd $MAINDIR

#Get lists of taxa that are not covered by the current contaminant and noncontaminant lists
Rscript $MAINDIR/scripts/reference_lists.R

#The first part of the script is very similar to the first round: looking for available genomes with gradually less stringent criteria
#But the output includes one more column regarding the original species that a genome is a substitute candidate for.


####For contaminants####
touch contaminant_genomes_list_suppl.txt

awk '{print $2}' exogenous_id_suppl.txt | while read j
do
    i=$(echo $j | tr -d '\r')
    original=$( awk -vLOOKUPVAL=$i 'BEGIN {FS="\t"}; $2 == LOOKUPVAL {print $1}' exogenous_id_suppl.txt ) 
    #save the TaxID of the original species in a variable
    #1 - RefSeq category: reference, version status: latest, assembly level: complete genome, genome representation: full
    grep -qP "${i}\t" contaminant_genomes_list_suppl.txt || ( awk -vLOOKUPVAL=$i -v original="$original" -v OFS='\t' '
    BEGIN {FS="\t"}; $7 == LOOKUPVAL && $5 == "reference" && $11 == "latest"  && $12 == "Complete genome" && $14 == "Full" {
    print original, $7, $5, $12, $14, $15, $20}' < assembly_summary_combined.txt >> contaminant_genomes_list_suppl.txt && 
    
    #2 - allow assembly level: chromosome
    grep -qP "${i}\t" contaminant_genomes_list_suppl.txt ) || ( awk -vLOOKUPVAL=$i -v original="$original" -v OFS='\t' '
    BEGIN {FS="\t"}; $7 == LOOKUPVAL && $5 == "reference" && $11 == "latest"  && $12 == "Chromosome" && $14 == "Full" {
    print original, $7, $5, $12, $14, $15, $20}' < assembly_summary_combined.txt >> contaminant_genomes_list_suppl.txt &&
    
    #3 - allow assembly level: scaffold
    grep -qP "${i}\t" contaminant_genomes_list_suppl.txt ) || ( awk -vLOOKUPVAL=$i -v original="$original" -v OFS='\t' '
    BEGIN {FS="\t"}; $7 == LOOKUPVAL && $5 == "reference" && $11 == "latest"  && $12 == "Scaffold" && $14 == "Full" {
    print original, $7, $5, $12, $14, $15, $20}' < assembly_summary_combined.txt >> contaminant_genomes_list_suppl.txt &&
    
    #4 - allow assembly level: contig
    grep -qP "${i}\t" contaminant_genomes_list_suppl.txt ) || ( awk -vLOOKUPVAL=$i -v original="$original" -v OFS='\t' '
    BEGIN {FS="\t"}; $7 == LOOKUPVAL && $5 == "reference" && $11 == "latest"  && $12 == "Contig" && $14 == "Full" {
    print original, $7, $5, $12, $14, $15, $20}' < assembly_summary_combined.txt >> contaminant_genomes_list_suppl.txt &&

    #5 - allow RefSeq category: NA, but go back to complete genomes
    grep -qP "${i}\t" contaminant_genomes_list_suppl.txt ) || ( awk -vLOOKUPVAL=$i -v original="$original" -v OFS='\t' '
    BEGIN {FS="\t"}; $7 == LOOKUPVAL && $5 == "na" && $11 == "latest"  && $12 == "Complete genome" && $14 == "Full" {
    print original, $7, $5, $12, $14, $15, $20}' < assembly_summary_combined.txt >> contaminant_genomes_list_suppl.txt &&
    
    #6 - allow RefSeq category: NA and assembly level: chromosome
    grep -qP "${i}\t" contaminant_genomes_list_suppl.txt ) || ( awk -vLOOKUPVAL=$i -v original="$original" -v OFS='\t' '
    BEGIN {FS="\t"}; $7 == LOOKUPVAL && $5 == "na" && $11 == "latest"  && $12 == "Chromosome" && $14 == "Full" {
    print original, $7, $5, $12, $14, $15, $20}' < assembly_summary_combined.txt >> contaminant_genomes_list_suppl.txt &&
    
    #7 - allow RefSeq category: NA and assembly level: scaffold
    grep -qP "${i}\t" contaminant_genomes_list_suppl.txt ) || ( awk -vLOOKUPVAL=$i -v original="$original" -v OFS='\t' '
    BEGIN {FS="\t"}; $7 == LOOKUPVAL && $5 == "na" && $11 == "latest"  && $12 == "Scaffold" && $14 == "Full" {
    print original, $7, $5, $12, $14, $15, $20}' < assembly_summary_combined.txt >> contaminant_genomes_list_suppl.txt &&
    
    #8 - allow RefSeq category: NA and assembly level: contig
    grep -qP "${i}\t" contaminant_genomes_list_suppl.txt ) || ( awk -vLOOKUPVAL=$i -v original="$original" -v OFS='\t' '
    BEGIN {FS="\t"}; $7 == LOOKUPVAL && $5 == "na" && $11 == "latest"  && $12 == "Contig" && $14 == "Full" {
    print original, $7, $5, $12, $14, $15, $20}' < assembly_summary_combined.txt >> contaminant_genomes_list_suppl.txt &&
    
    #9 - allow partial, chromosome-level assemblies
    grep -qP "${i}\t" contaminant_genomes_list_suppl.txt ) || ( awk -vLOOKUPVAL=$i -v original="$original" -v OFS='\t' '
    BEGIN {FS="\t"};  $7 == LOOKUPVAL && $5 == "na" && $11 == "latest"  && $12 == "Chromosome" && $14 == "Partial" {
    print original, $7, $5, $12, $14, $15, $20}' < assembly_summary_combined.txt >> contaminant_genomes_list_suppl.txt &&
    
    #10 - allow partial, scaffold-level assemblies
    grep -qP "${i}\t" contaminant_genomes_list_suppl.txt ) || ( awk -vLOOKUPVAL=$i -v original="$original" -v OFS='\t' '
    BEGIN {FS="\t"};  $7 == LOOKUPVAL && $5 == "na" && $11 == "latest"  && $12 == "Scaffold" && $14 == "Partial" {
    print original, $7, $5, $12, $14, $15, $20}' < assembly_summary_combined.txt >> contaminant_genomes_list_suppl.txt &&
    
    #11 - allow partial, contig-level assemblies
    grep -qP "${i}\t" contaminant_genomes_list_suppl.txt ) || ( awk -vLOOKUPVAL=$i -v original="$original" -v OFS='\t' '
    BEGIN {FS="\t"};  $7 == LOOKUPVAL && $5 == "na" && $11 == "latest"  && $12 == "Contig" && $14 == "Partial" {
    print original, $7, $5, $12, $14, $15, $20}' < assembly_summary_combined.txt >> contaminant_genomes_list_suppl.txt
    
    #12 - finally, any match
    grep -qP "${i}\t" contaminant_genomes_list_suppl.txt ) ||  awk -vLOOKUPVAL=$i -v original="$original" -v OFS='\t' '
    BEGIN {FS="\t"};  $7 == LOOKUPVAL {
    print original, $7, $5, $12, $14, $15, $20}' < assembly_summary_combined.txt >> contaminant_genomes_list_suppl.txt
done

#Produce a new file with only the most recently published genomes per taxon

touch contaminant_genomes_list_recent_suppl.txt

#I will select the most recent entry per original species (not substitute species)

for i in $(awk '{print $1}' contaminant_genomes_list_suppl.txt | sort | uniq )
do
    grep -P "${i}\t" contaminant_genomes_list_suppl.txt | sort -t$'\t' -k6 | tail -1 >> contaminant_genomes_list_recent_suppl.txt
done

#Download genomes and place in directory named after the TaxID (inside the contam_references folder
awk -v OFS='\t' 'BEGIN {FS="\t"}; {print $7}' contaminant_genomes_list_recent_suppl.txt | while read i
do
    cd $REFDIR/contam_references
    dirname=$( grep "$i" $MAINDIR/contaminant_genomes_list_recent_suppl.txt | awk '{print $2}' )
    echo "Downloading $dirname ($i) as a substitute for $( grep "$i" $MAINDIR/contaminant_genomes_list_recent_suppl.txt | awk '{print $1}' )"
    mkdir $dirname
    cd $dirname
    wget -np "$i"/*[!m]_genomic.fna.gz
    cd $MAINDIR
done

echo "Total number of genomes downloaded (including first round of selections):"
ls $REFDIR/contam_references/**/**_genomic.fna.gz | wc -l


####For non-contaminants ####
touch noncontaminant_genomes_list_suppl.txt

awk '{print $2}' abundant_id_suppl.txt | while read j
do
    i=$(echo $j | tr -d '\r')
    original=$( awk -vLOOKUPVAL=$i 'BEGIN {FS="\t"}; $2 == LOOKUPVAL {print $1}' exogenous_id_suppl.txt ) 
    #save the TaxID of the original species in a variable
    #1 - RefSeq category: reference, version status: latest, assembly level: complete genome, genome representation: full
    grep -qP "${i}\t" contaminant_genomes_list_suppl.txt || ( awk -vLOOKUPVAL=$i -v original="$original" -v OFS='\t' '
    BEGIN {FS="\t"}; $7 == LOOKUPVAL && $5 == "reference" && $11 == "latest"  && $12 == "Complete genome" && $14 == "Full" {
    print original, $7, $5, $12, $14, $15, $20}' < assembly_summary_combined.txt >> contaminant_genomes_list_suppl.txt && 
    
    #2 - allow assembly level: chromosome
    grep -qP "${i}\t" contaminant_genomes_list_suppl.txt ) || ( awk -vLOOKUPVAL=$i -v original="$original" -v OFS='\t' '
    BEGIN {FS="\t"}; $7 == LOOKUPVAL && $5 == "reference" && $11 == "latest"  && $12 == "Chromosome" && $14 == "Full" {
    print original, $7, $5, $12, $14, $15, $20}' < assembly_summary_combined.txt >> contaminant_genomes_list_suppl.txt &&
    
    #3 - allow assembly level: scaffold
    grep -qP "${i}\t" contaminant_genomes_list_suppl.txt ) || ( awk -vLOOKUPVAL=$i -v original="$original" -v OFS='\t' '
    BEGIN {FS="\t"}; $7 == LOOKUPVAL && $5 == "reference" && $11 == "latest"  && $12 == "Scaffold" && $14 == "Full" {
    print original, $7, $5, $12, $14, $15, $20}' < assembly_summary_combined.txt >> contaminant_genomes_list_suppl.txt &&
    
    #4 - allow assembly level: contig
    grep -qP "${i}\t" contaminant_genomes_list_suppl.txt ) || ( awk -vLOOKUPVAL=$i -v original="$original" -v OFS='\t' '
    BEGIN {FS="\t"}; $7 == LOOKUPVAL && $5 == "reference" && $11 == "latest"  && $12 == "Contig" && $14 == "Full" {
    print original, $7, $5, $12, $14, $15, $20}' < assembly_summary_combined.txt >> contaminant_genomes_list_suppl.txt &&

    #5 - allow RefSeq category: NA, but go back to complete genomes
    grep -qP "${i}\t" contaminant_genomes_list_suppl.txt ) || ( awk -vLOOKUPVAL=$i -v original="$original" -v OFS='\t' '
    BEGIN {FS="\t"}; $7 == LOOKUPVAL && $5 == "na" && $11 == "latest"  && $12 == "Complete genome" && $14 == "Full" {
    print original, $7, $5, $12, $14, $15, $20}' < assembly_summary_combined.txt >> contaminant_genomes_list_suppl.txt &&
    
    #6 - allow RefSeq category: NA and assembly level: chromosome
    grep -qP "${i}\t" contaminant_genomes_list_suppl.txt ) || ( awk -vLOOKUPVAL=$i -v original="$original" -v OFS='\t' '
    BEGIN {FS="\t"}; $7 == LOOKUPVAL && $5 == "na" && $11 == "latest"  && $12 == "Chromosome" && $14 == "Full" {
    print original, $7, $5, $12, $14, $15, $20}' < assembly_summary_combined.txt >> contaminant_genomes_list_suppl.txt &&
    
    #7 - allow RefSeq category: NA and assembly level: scaffold
    grep -qP "${i}\t" contaminant_genomes_list_suppl.txt ) || ( awk -vLOOKUPVAL=$i -v original="$original" -v OFS='\t' '
    BEGIN {FS="\t"}; $7 == LOOKUPVAL && $5 == "na" && $11 == "latest"  && $12 == "Scaffold" && $14 == "Full" {
    print original, $7, $5, $12, $14, $15, $20}' < assembly_summary_combined.txt >> contaminant_genomes_list_suppl.txt &&
    
    #8 - allow RefSeq category: NA and assembly level: contig
    grep -qP "${i}\t" contaminant_genomes_list_suppl.txt ) || ( awk -vLOOKUPVAL=$i -v original="$original" -v OFS='\t' '
    BEGIN {FS="\t"}; $7 == LOOKUPVAL && $5 == "na" && $11 == "latest"  && $12 == "Contig" && $14 == "Full" {
    print original, $7, $5, $12, $14, $15, $20}' < assembly_summary_combined.txt >> contaminant_genomes_list_suppl.txt &&
    
    #9 - allow partial, chromosome-level assemblies
    grep -qP "${i}\t" contaminant_genomes_list_suppl.txt ) || ( awk -vLOOKUPVAL=$i -v original="$original" -v OFS='\t' '
    BEGIN {FS="\t"};  $7 == LOOKUPVAL && $5 == "na" && $11 == "latest"  && $12 == "Chromosome" && $14 == "Partial" {
    print original, $7, $5, $12, $14, $15, $20}' < assembly_summary_combined.txt >> contaminant_genomes_list_suppl.txt &&
    
    #10 - allow partial, scaffold-level assemblies
    grep -qP "${i}\t" contaminant_genomes_list_suppl.txt ) || ( awk -vLOOKUPVAL=$i -v original="$original" -v OFS='\t' '
    BEGIN {FS="\t"};  $7 == LOOKUPVAL && $5 == "na" && $11 == "latest"  && $12 == "Scaffold" && $14 == "Partial" {
    print original, $7, $5, $12, $14, $15, $20}' < assembly_summary_combined.txt >> contaminant_genomes_list_suppl.txt &&
    
    #11 - allow partial, contig-level assemblies
    grep -qP "${i}\t" contaminant_genomes_list_suppl.txt ) || ( awk -vLOOKUPVAL=$i -v original="$original" -v OFS='\t' '
    BEGIN {FS="\t"};  $7 == LOOKUPVAL && $5 == "na" && $11 == "latest"  && $12 == "Contig" && $14 == "Partial" {
    print original, $7, $5, $12, $14, $15, $20}' < assembly_summary_combined.txt >> contaminant_genomes_list_suppl.txt
    
    #12 - finally, any match
    grep -qP "${i}\t" contaminant_genomes_list_suppl.txt ) ||  awk -vLOOKUPVAL=$i -v original="$original" -v OFS='\t' '
    BEGIN {FS="\t"};  $7 == LOOKUPVAL {
    print original, $7, $5, $12, $14, $15, $20}' < assembly_summary_combined.txt >> contaminant_genomes_list_suppl.txt
done

#Produce a new file with only the most recently published genomes per taxon

touch noncontaminant_genomes_list_recent_suppl.txt

#I will select the most recent entry per original species (not substitute species)

for i in $(awk '{print $1}' noncontaminant_genomes_list_suppl.txt | sort | uniq )
do
    grep -P "${i}\t" noncontaminant_genomes_list_suppl.txt | sort -t$'\t' -k6 | tail -1 >> noncontaminant_genomes_list_recent_suppl.txt
done

#Download genomes and place in directory named after the TaxID (inside the contam_references folder
awk -v OFS='\t' 'BEGIN {FS="\t"}; {print $7}' noncontaminant_genomes_list_recent_suppl.txt | while read i
do
    cd $REFDIR/noncontam_references
    dirname=$( grep "$i" $MAINDIR/noncontaminant_genomes_list_recent_suppl.txt | awk '{print $2}' )
    echo "Downloading $dirname ($i) as a substitute for $( grep "$i" $MAINDIR/noncontaminant_genomes_list_recent_suppl.txt | awk '{print $1}' )"
    mkdir $dirname
    cd $dirname
    wget -np "$i"/*[!m]_genomic.fna.gz
    cd $MAINDIR
done

echo "Total number of genomes downloaded (including first round of selections):"
ls $REFDIR/noncontam_references/**/**_genomic.fna.gz | wc -l


