#!/bin/bash -l

#SBATCH -A snic2020-5-528
#SBATCH -p devel
#SBATCH -N 1-1
#SBATCH -t 1:00:00
#SBATCH -J haplotype_net
#SBATCH --mail-user=markella.moraitou.0437@student.uu.se
#SBATCH --mail-type=ALL

#This script concatenates mitochondrial genomes assembled from DC (with completeness above a certain threshold) and downloaded from NCBI,
#and prepares the input needed for PopArt to plot the network
DATADIR=/proj/sllstore2017021/nobackup/MARKELLA/H1_mtgenomes
OUTDIR=/proj/sllstore2017021/nobackup/MARKELLA/H2_haplotype_network
REFDIR=/proj/sllstore2017021/nobackup/MARKELLA/downloaded_mtgenomes
DBPREFIX=/proj/sllstore2017021/nobackup/JAELLE/REFERENCES/gorilla_gorilla_human/calculushost
coveragefile=$DATADIR/mt_coverage.csv

module load bioinfo-tools samtools python R_packages FastANI

### References + DC genomes above the threshold will be used to construct a network with Popart ###
threshold=80

#Bring all genomes to the same length
mkdir $OUTDIR/genomes_for_popart

#First, get the downloaded genomes selected to be included in the haplotype network
cd $REFDIR

#Trim reference genomes to remove hypervariable region 
cat gorilla_mt_genomes_to_download.txt | awk '{print $1}' | while read j
do
    ls ${j}*fasta.gz
done | uniq | while read i
do
   gunzip $i
   #Get chromosome name
   chr=$(head -1 ${i%.gz} | awk '{print $1}' | sed "s/>//g" )
   samtools faidx ${i%.gz} "${chr}:1-15446" -o $OUTDIR/genomes_for_popart/${i%fasta.gz}trimmed.fasta
   sed -i 's/:1-15446//g' $OUTDIR/genomes_for_popart/${i%fasta.gz}trimmed.fasta #Remove suffix from header
   gzip ${i%.gz}
   gzip $OUTDIR/genomes_for_popart/${i%fasta.gz}trimmed.fasta
done

cd $OUTDIR

#Only complete genomes with be used as references for FastANI
ls genomes_for_popart/*fasta.gz > $OUTDIR/fastani_ref_list.txt

#Also select the DC-assembled mt genomes above the threshold
cd $DATADIR
awk -v FS=", " -v threshold=$threshold 'NR > 1 && $2 > threshold {print $1}' $coveragefile | while read j
do
    ls ${j}*fa.gz
done | uniq | while read i
do
    gunzip $i
    #Get chromosome name
    chr=$(head -1 ${i%.gz} | awk '{print $1}' | sed "s/>//g" )
    samtools faidx ${i%.gz} "${chr}:1-15446" -o $OUTDIR/genomes_for_popart/${i%fa.gz}trimmed.fasta
    sed -i 's/:1-15446//g' $OUTDIR/genomes_for_popart/${i%fa.gz}trimmed.fasta #Remove suffix from header
    gzip ${i%.gz}
    gzip $OUTDIR/genomes_for_popart/${i%fa.gz}trimmed.fasta
done

cd $OUTDIR

#List all files that are as long as the reference
echo "" > $OUTDIR/fasta_list.txt

cd $OUTDIR/genomes_for_popart/
ls *.fasta.gz | while read i
do
    #Get sequence length (number of characters excluding title, minus newlines)
    length=$(zcat $i | tail -n +2 | wc | awk '{print $3-$1}')
    #If the sequence is smaller
    if [[ 15445 -gt $length ]]
    then
        continue
    else
        echo $OUTDIR/genomes_for_popart/$i
    fi
done  >> $OUTDIR/fasta_list.txt

cd $OUTDIR

#Concatenate all fasta files that were listed, in preparation for multiple sequence alignment
rm  mt_genomes.fa*
cat $OUTDIR/fasta_list.txt | xargs cat > mt_genomes.fa.gz

#decompress
pigz -d -p $SLURM_CPUS_ON_NODE mt_genomes.fa.gz

#Get NEXUS format
conda activate /home/markmora/.conda/envs/seqmagick
/home/markmora/.local/bin/seqmagick convert --output-format nexus --alphabet dna mt_genomes.fa mt_genomes.nex

#Get haplotype network with pegas
#Rscript plot_hapnets.r

#Append traits block on nexus file
cp mt_genomes.nex mt_genomes_with_traits.nex

#First, append the header part of the traits block
echo "BEGIN TRAITS;
Dimensions NTRAITS=5;
Format labels=yes missing=? separator=Comma;
TraitLabels gorilla graueri beringei DC reference;
Matrix" >> mt_genomes_with_traits.nex

#Then select the rows of the traits matrix that match the sequences
grep ">" mt_genomes.fa | sed 's/>//g' > mt_genomes_samples.txt
grep -f mt_genomes_samples.txt traits_matrix.txt >> mt_genomes_with_traits.nex

#Finally, append the closing part
echo ";
end;" >> mt_genomes_with_traits.nex
#Use this output on PopArt (desktop)

#compress
pigz -p $SLURM_CPUS_ON_NODE mt_genomes.fa

### DC genomes below the threshold will be placed near the closest haplotype (genome with the smallest ANI) ###

mkdir $OUTDIR/genomes_for_fastani

cd $DATADIR
#Also use a lower threshold
awk -v FS=", " -v threshold=$threshold '$2 < threshold && $2 > 8 {print $1}' $coveragefile | while read j
do
    ls ${j}*fa.gz
done | uniq | while read i
do
    gunzip $i
    #Get chromosome name
    chr=$(head -1 ${i%.gz} | awk '{print $1}' | sed "s/>//g" )
    samtools faidx ${i%.gz} "${chr}:1-15446" -o $OUTDIR/genomes_for_fastani/${i%fa.gz}trimmed.fasta
    sed -i 's/:1-15446//g' $OUTDIR/genomes_for_fastani/${i%fa.gz}trimmed.fasta #Remove suffix from header
    gzip $OUTDIR/genomes_for_fastani/${i%fa.gz}trimmed.fasta
    gzip ${i%.gz}
done

#List all incomplete genomes to be used as queries with fastani
cd $OUTDIR

ls genomes_for_fastani/*fasta.gz > $OUTDIR/fastani_query_list.txt

#Run FastANI
fastANI --ql fastani_query_list.txt --rl fastani_ref_list.txt -o fastani_output.txt
sed -i 's/genomes_for_fastani//g' fastani_output.txt
sed -i 's/genomes_for_popart//g' fastani_output.txt

#Get the refernces with the highest ANI per query
Rscript find_closest_haplotype.R

rm -r $OUTDIR/genomes_for_popart
rm -r $OUTDIR/genomes_for_fastani
