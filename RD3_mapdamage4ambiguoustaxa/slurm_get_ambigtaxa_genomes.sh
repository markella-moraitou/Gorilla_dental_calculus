#!/bin/bash -l

#SBATCH -A snic2020-5-528
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 8:00:00
#SBATCH -J get_ambigtaxa_genomes
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Markella.Moraitou.0437@student.uu.se

#This script parses the assembly summary files for bacteria and archaea(ftp://ftp.ncbi.nih.gov/genomes/genbank/bacteria/assembly_summary.txt) and viruses (ftp://ftp.ncbi.nih.gov/genomes/genbank/bacteria/assembly_summary_viral.txt)
#and saves ftp links to genome assemblies for each of the provided TaxIDs
#Looks genomes up with increasily relaxed criteria: Starts looking for complete reference genomes 
#and finishes with non-reference genomes at the scaffold-level
#At every step, it checks if a genome has already been listed for that taxon, so that it doesn't continue searching with less stringent criteria

echo $SLURM_JOB_NAME
echo $(module list)

OUTDIR=/proj/sllstore2017021/nobackup/MARKELLA/RD3_mapdamage4ambiguoustaxa
REFDIR=$OUTDIR/genomes
assembly_summary=/proj/sllstore2017021/nobackup/MARKELLA/RD2_mapping/assembly_summary_combined.txt

cd $OUTDIR

touch ambiguous_genomes_list.txt

cat top_ambiguous_taxa.txt | awk '{print $2}'| while read j
do
    i=$(echo $j | tr -d '\r')
    #1 - RefSeq category: reference, version status: latest, assembly level: complete genome, genome representation: full
    grep -qP "${i}\t" ambiguous_genomes_list.txt || ( awk -vLOOKUPVAL=$i -v OFS='\t' '
    BEGIN {FS="\t"}; $7 == LOOKUPVAL && $5 == "reference" && $11 == "latest"  && $12 == "Complete genome" && $14 == "Full" {
    print $7, $5, $12, $14, $15, $20}' < $assembly_summary >> ambiguous_genomes_list.txt && 
    
    #2 - allow assembly level: chromosome
    grep -qP "${i}\t" ambiguous_genomes_list.txt ) || ( awk -vLOOKUPVAL=$i -v OFS='\t' '
    BEGIN {FS="\t"}; $7 == LOOKUPVAL && $5 == "reference" && $11 == "latest"  && $12 == "Chromosome" && $14 == "Full" {
    print $7, $5, $12, $14, $15, $20}' < $assembly_summary >> ambiguous_genomes_list.txt &&
    
    #3 - allow assembly level: scaffold
    grep -qP "${i}\t" ambiguous_genomes_list.txt ) || ( awk -vLOOKUPVAL=$i -v OFS='\t' '
    BEGIN {FS="\t"}; $7 == LOOKUPVAL && $5 == "reference" && $11 == "latest"  && $12 == "Scaffold" && $14 == "Full" {
    print $7, $5, $12, $14, $15, $20}' < $assembly_summary >> ambiguous_genomes_list.txt &&
    
    #4 - allow assembly level: contig
    grep -qP "${i}\t" ambiguous_genomes_list.txt ) || ( awk -vLOOKUPVAL=$i -v OFS='\t' '
    BEGIN {FS="\t"}; $7 == LOOKUPVAL && $5 == "reference" && $11 == "latest"  && $12 == "Contig" && $14 == "Full" {
    print $7, $5, $12, $14, $15, $20}' < $assembly_summary >> ambiguous_genomes_list.txt &&

    #5 - allow RefSeq category: NA, but go back to complete genomes
    grep -qP "${i}\t" ambiguous_genomes_list.txt ) || ( awk -vLOOKUPVAL=$i -v OFS='\t' '
    BEGIN {FS="\t"}; $7 == LOOKUPVAL && $5 == "na" && $11 == "latest"  && $12 == "Complete genome" && $14 == "Full" {
    print $7, $5, $12, $14, $15, $20}' < $assembly_summary >> ambiguous_genomes_list.txt &&
    
    #6 - allow RefSeq category: NA and assembly level: chromosome
    grep -qP "${i}\t" ambiguous_genomes_list.txt ) || ( awk -vLOOKUPVAL=$i -v OFS='\t' '
    BEGIN {FS="\t"}; $7 == LOOKUPVAL && $5 == "na" && $11 == "latest"  && $12 == "Chromosome" && $14 == "Full" {
    print $7, $5, $12, $14, $15, $20}' < $assembly_summary >> ambiguous_genomes_list.txt &&
    
    #7 - allow RefSeq category: NA and assembly level: scaffold
    grep -qP "${i}\t" ambiguous_genomes_list.txt ) || ( awk -vLOOKUPVAL=$i -v OFS='\t' '
    BEGIN {FS="\t"}; $7 == LOOKUPVAL && $5 == "na" && $11 == "latest"  && $12 == "Scaffold" && $14 == "Full" {
    print $7, $5, $12, $14, $15, $20}' < $assembly_summary >> ambiguous_genomes_list.txt &&
    
    #8 - allow RefSeq category: NA and assembly level: contig
    grep -qP "${i}\t" ambiguous_genomes_list.txt ) || ( awk -vLOOKUPVAL=$i -v OFS='\t' '
    BEGIN {FS="\t"}; $7 == LOOKUPVAL && $5 == "na" && $11 == "latest"  && $12 == "Contig" && $14 == "Full" {
    print $7, $5, $12, $14, $15, $20}' < $assembly_summary >> ambiguous_genomes_list.txt &&
    
    #9 - allow partial, chromosome-level assemblies
    grep -qP "${i}\t" ambiguous_genomes_list.txt ) || ( awk -vLOOKUPVAL=$i -v OFS='\t' '
    BEGIN {FS="\t"};  $7 == LOOKUPVAL && $5 == "na" && $11 == "latest"  && $12 == "Chromosome" && $14 == "Partial" {
    print $7, $5, $12, $14, $15, $20}' < $assembly_summary >> ambiguous_genomes_list.txt &&
    
    #10 - allow partial, scaffold-level assemblies
    grep -qP "${i}\t" ambiguous_genomes_list.txt ) || ( awk -vLOOKUPVAL=$i -v OFS='\t' '
    BEGIN {FS="\t"};  $7 == LOOKUPVAL && $5 == "na" && $11 == "latest"  && $12 == "Scaffold" && $14 == "Partial" {
    print $7, $5, $12, $14, $15, $20}' < $assembly_summary >> ambiguous_genomes_list.txt &&
    
    #11 - allow partial, contig-level assemblies
    grep -qP "${i}\t" ambiguous_genomes_list.txt ) || ( awk -vLOOKUPVAL=$i -v OFS='\t' '
    BEGIN {FS="\t"};  $7 == LOOKUPVAL && $5 == "na" && $11 == "latest"  && $12 == "Contig" && $14 == "Partial" {
    print $7, $5, $12, $14, $15, $20}' < $assembly_summary >> ambiguous_genomes_list.txt
    
    #12 - finally, any match
    grep -qP "${i}\t" ambiguous_genomes_list.txt ) ||  awk -vLOOKUPVAL=$i -v OFS='\t' '
    BEGIN {FS="\t"};  $7 == LOOKUPVAL {
    print $7, $5, $12, $14, $15, $20}' < $assembly_summary >> ambiguous_genomes_list.txt
done

#Produce a new file with only the most recently published genomes per taxon

touch ambiguous_genomes_list_recent.txt

for i in $(awk '{print $1}' ambiguous_genomes_list.txt | sort | uniq )
do
    grep -P "^${i}\t" ambiguous_genomes_list.txt | sort -t$'\t' -k5 | tail -1 >> ambiguous_genomes_list_recent.txt
done

#Download genomes and place in directory named after the TaxID
awk -v OFS='\t' 'BEGIN {FS="\t"}; {print $6}' ambiguous_genomes_list_recent.txt | while read i
do
    cd $REFDIR
    #create directory name (unless already there)
    dirname=$( grep "$i" $OUTDIR/ambiguous_genomes_list_recent.txt | awk '{print $1}' )
    [[ -d $REFDIR/$dirname ]] || mkdir $dirname
    cd $dirname
    #download file unless already there
    ls | grep -q *[!m]_genomic.fna.gz || ( echo "Downloading $dirname ($i)" && wget -np "$i"/*[!m]_genomic.fna.gz )
    cd $OUTDIR
done

echo "Number of genomes downloaded:"
ls $REFDIR/**/**_genomic.fna.gz | wc -l
