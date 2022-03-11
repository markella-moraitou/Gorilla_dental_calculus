#!/bin/bash -l

#SBATCH -A SNIC_PROJECT_ID
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 30:00:00
#SBATCH -J mapDamage
#SBATCH --mail-type=ALL
#SBATCH --mail-user=USER_EMAIL

#Calculate damage patterns of abundant ambiguous taxa(found in both contaminant and oral lists) 
#The sample used for the mapping is the one where the taxon is the most abundant (excluding bad samples according to FEAST)

module load bioinfo-tools bbmap mapDamage samtools pigz

echo "$SLURM_JOB_NAME"
echo $(module list)

#BAM files to be used
readdir=8_humanHostFilt/unmapped
outdir=RD3_mapdamage4ambiguoustaxa
refdir=$outdir/genomes

cd $outdir || exit

#Create file to save base change stats from each taxon
echo -e "TaxID\tChr\tEnd\tStd\tPos\tA\tC\tG\tT\tTotal\tG>A\tC>T\tA>G\tT>C\tA>C\tA>T\tC>G\tC>A\tT>G\tT>A\tG>C\tG>T\tA>-\tT>-\tC>-\tG>-\t->A\t->T\t->C\t->G\tS" > $outdir/edge_misincorporation.txt

#For every taxon-sample pair in the list
cat top_ambiguous_taxa.txt | while read j
do
    i=$(echo "$j" | tr -d '\r')
    #Get sample name
    sample=$(echo "$i" | awk '{print $1}')
    #Get taxon name
    taxon=$(echo "$i" | awk '{print $2}')
    #Get the reference genome path
    refgen=$(ls $refdir/"${taxon}"/*.fna.gz)
    #Check if references are there (some might be missing)
    ls genomes | grep -qw "$taxon" && (
    #First remove all spaces from header
    pigz -d -p "$SLURM_CPUS_ON_NODE" "$refgen"
    sed -i 's/ /_/g' "${refgen%.gz}"
    #Map using bbmap - change minimum identity to 0.95 (ANI threshold for the same species is 0.95)
    bbwrap.sh ref="${refgen%.gz}" in=$readdir/"${sample}"_m_host_unmapped.fastq.gz mapper=bbmap out=$outdir/"${taxon}".sam  maxindel=80 minid=0.95 pigz=t unpigz=t nodisk
    #Save mapped reads into BAM file
    samtools view -Sb -F 4 -@ "$SLURM_CPUS_ON_NODE" $outdir/"${taxon}".sam -o $outdir/"${taxon}".bam
  
    #Run mapDamage
    cd mapdamage || exit
    mapDamage -i $outdir/"${taxon}".bam -r "${refgen%.gz}"
    #Re-zip reference
    pigz -p "$SLURM_CPUS_ON_NODE" "${refgen%.gz}"
    
    #Get stats on damage patterns
    cd results_${taxon} || exit
    #From the misincorporation plot, get info on the first three positions from the 3' and 5' edge
    awk -v taxon="$taxon" -v OFS="\t" 'BEGIN {FS="\t"}; $3 == "+" && $4 < 4 {print taxon, $0}' misincorporation.txt >> $outdir/edge_misincorporation.txt
    cd $outdir || exit
    
    rm $outdir/"${taxon}".sam )
done

#Combine all PDFs in one
pdfunite $(ls -rt mapdamage/results_**/Fragmisincorporation_plot.pdf) Fragmisincorporation_plot_combined.pdf

rm mapdamage/results_**/Fragmisincorporation_plot.pdf
