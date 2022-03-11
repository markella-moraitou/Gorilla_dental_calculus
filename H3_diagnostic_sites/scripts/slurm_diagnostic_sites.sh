#!/bin/bash -l

#SBATCH -A SNIC_PROJECT_ID
#SBATCH -p devcore
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -J haplotype_net
#SBATCH --mail-user=USER_EMAIL
#SBATCH --mail-type=ALL

#Get diagnostic sites to distinguish between subspecies

module load bioinfo-tools python3 ANGSD R_packages samtools
module list

DATADIR=H1_mtgenomes
OUTDIR=H3_diagnostic_sites
REFLIST=H2_haplotype_network/fasta_list.txt #List of full mitochondrial references

## First, create concatenated fasta files
#For the references
grep "1G[bg][bg]" $REFLIST | xargs zcat > $OUTDIR/concatenated_mt_refs.fa
#For the samples to be investigated
zcat $DATADIR/*.gz > $OUTDIR/concatenated_mt_samples.fa

cd $OUTDIR || exit

mkdir $OUTDIR/output
## Find diagnostic sites using the references - The scripts compares population pairs, so first I have to find the differences between western and eastern gorillas and then between Grauer's and mountain gorillas

#Create popfile to use with the python script (should have the sample/header name in the first column and the population assignment in the second)
conda activate diagnostic_sites

#Diagnostic sites between eastern and western gorillas
fastafile=H3_diagnostic_sites/concatenated_mt_refs.fa
popfile=H3_diagnostic_sites/popfile_east_west.txt
outfile=H3_diagnostic_sites/output/diagnostic_sites_east_west.txt

python3 $OUTDIR/scripts/FindDiagnosticSites_edited.py --input $fastafile --popfile $popfile --output $outfile

#Diagnostic sites between Grauer's and Mountain gorillas
popfile=H3_diagnostic_sites/popfile_mountain_grauers.txt
outfile=H3_diagnostic_sites/output/diagnostic_sites_mountain_grauers.txt

python3 $OUTDIR/scripts/FindDiagnosticSites_edited.py --input $fastafile --popfile $popfile --output $outfile

## Measure distance to references: count how many diagnostic site each sample has for each lineage

#For eastern-western diagnostic sites
fastafile=H3_diagnostic_sites/concatenated_mt_samples.fa
reffile=H3_diagnostic_sites/output/diagnostic_sites_east_west.txt
outfile=H3_diagnostic_sites/output/distance_to_refs_east_west.txt

python3 $OUTDIR/scripts/DistanceToReferences.py --input $fastafile --reference $reffile --output $outfile

#For Grauer's-Mountain diagnostic sites
reffile=H3_diagnostic_sites/output/diagnostic_sites_mountain_grauers.txt
outfile=H3_diagnostic_sites/output/distance_to_refs_mountain_grauers.txt

python3 $OUTDIR/scripts/DistanceToReferences.py --input $fastafile --reference $reffile --output $outfile

gzip *fa
rm *fai

#Obtain allele counts at all bases at all sites where coverage > 0
#Code below adapted from /proj/sllstore2017021/nobackup/SAM/Bear_Lineage_Project/scripts/1run_angsd_on_mt_copy.sh

cd $DATADIR || exit
ls *MT_RG.bam | while read i
do
    output=${i%.bam}
    cd $OUTDIR/output || exit
    angsd -out "$output" -doFasta 2 -doCounts 1 -dumpCounts 3 -i $DATADIR/"$i"
    zcat "${output}".pos.gz | tail -n +2 | cut -f 2 > tmp.positions
    #gunzip countsfile for pasting
    gunzip "${output}".counts.gz
    
    #creating a temporary file per base with position in first column, allele depth in second and allele in third (this is messy but should work)
    for n in {1..4}
    do
        allele=$(head -n 1 "${output}".counts | cut -f "$n")
        cat "${output}".counts | cut -f "$n" | tail -n +2 | paste tmp.positions - > tmp.counts_"$allele"
        awk -v allele="${allele/tot/}" ' {print $0, allele} ' tmp.counts_"$allele" > tmp.counts_"$allele"_2
        mv tmp.counts_"$allele"_2 tmp.counts_"$allele"
        cat tmp.counts_"$allele"
    done > "${output}".allele_counts

    #remove tempfiles
    rm -f tmp.*
    
    #Plot read coverage
    Rscript $OUTDIR/scripts/plot_readCoverage.r "${output}".allele_counts
    gzip "${output}".counts
    
    #Extract reads spanning diagnostic sites
    samtools index $DATADIR/"$i"
    #First for east-west diagnostic sites
    sites_e_w=$(tail -n +2 diagnostic_sites_east_west.txt | awk '{print $1}' | while read i
            do
                echo "gorilla_NC_011120.1:${i}-${i}"
            done | sed -z 's/\n/ /g')
    diagnostic_reads_e_w=$(samtools view $DATADIR/"$i" "$sites_e_w" | sort | uniq | wc -l)
    #Write number of diagnostic reads to file
    awk -v sample="${i%%_*}" -v nreads="$diagnostic_reads_e_w" -v OFS="\t" '$1==sample {print $0, nreads}' distance_to_refs_east_west.txt >> distance_to_refs_east_west_copy.txt
    #Second for mountain-Grauer's diagnostic sites
    sites_m_g=$(tail -n +2 diagnostic_sites_mountain_grauers.txt | awk '{print $1}' | while read i
            do
                echo "gorilla_NC_011120.1:${i}-${i}"
            done | sed -z 's/\n/ /g')
    diagnostic_reads_m_g=$(samtools view $DATADIR/"$i" "$sites_m_g" | sort | uniq | wc -l)
    #Write number of diagnostic reads to file
    awk -v sample="${i%%_*}" -v nreads="$diagnostic_reads_m_g" -v OFS="\t" '$1==sample {print $0, nreads}' distance_to_refs_mountain_grauers.txt >> distance_to_refs_mountain_grauers_copy.txt
done

cd $OUTDIR/output || exit
#Update files and add headers
echo -e "$(head -1 distance_to_refs_east_west.txt)\treads.spanning.sites" > distance_to_refs_east_west.txt
cat distance_to_refs_east_west_copy.txt >> distance_to_refs_east_west.txt
rm distance_to_refs_east_west_copy.txt


echo -e "$(head -1 distance_to_refs_mountain_grauers.txt)\treads.spanning.sites" > distance_to_refs_mountain_grauers.txt
cat distance_to_refs_mountain_grauers_copy.txt >> distance_to_refs_mountain_grauers.txt
rm distance_to_refs_mountain_grauers_copy.txt
