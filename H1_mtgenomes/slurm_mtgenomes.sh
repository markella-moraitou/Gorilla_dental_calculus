#!/bin/bash -l

#SBATCH -A SNIC_PROJECT_ID
#SBATCH -p devel
#SBATCH -N 1-1
#SBATCH -t 1:00:00
#SBATCH -J mtgenomes
#SBATCH --mail-user=USER_EMAIL
#SBATCH --mail-type=ALL

# This script takes fasta reads and a reference genome, mappes the reads to the reference (indexes the reference
# if needed), creates a bam file that is then subsetted to only include mitochondrial genome, adds read group data,
# and uses angsd to create a consensus fasta file.

#Adapted from: /proj/sllstore2017021/nobackup/SAM/Bear_Lineage_Project/scripts/bwa.RG.dup.sh

module load bioinfo-tools BEDTools bwa samtools picard ANGSD

echo "$SLURM_JOB_NAME"
echo $(module list)

OUTDIR=H1_mtgenomes
DATADIR=8_humanHostFilt/mapped
DBPREFIX=/proj/sllstore2017021/nobackup/JAELLE/REFERENCES/gorilla_gorilla_human/calculushost
#Samples that were retained in the analysis
retainedsamples=RD2_mapping/retained_samples.txt

cd $DATADIR/ || exit

#Save coverage info
echo "sample, coverage, Nreads" > $OUTDIR/mt_coverage.csv

ls *m_host_mapped_sorted.bam | grep -f $retainedsamples | while read i
do
    #Get fastq file
    bedtools bamtofastq -i "$i" -fq $OUTDIR/"${i%.bam}".fastq
    #The file that will be mapped is the fastq
    file=$OUTDIR/${i%.bam}.fastq

    #Get sample name
    name=$(basename "$file")
    SM=$(echo "$name" | cut -d _ -f 1)

    echo "Starting BWA"
    bwa mem -t "$SLURM_CPUS_ON_NODE" -M $DBPREFIX "$file" | samtools view -S -b -h -q 30 - > $OUTDIR/"$SM".bam

    cd $OUTDIR || exit
    #Map reads, filter q > 30, sort and index
    echo "Starting processing with samtools"
    samtools sort $OUTDIR/"$SM".bam > $OUTDIR/"$SM"_sort.bam
    samtools index $OUTDIR/"$SM"_sort.bam

    #Get only reads mapping to the mitochondrial genome
    #exlude d-loop (15447-16364) because of the hypervariable region
    samtools view -h -b $OUTDIR/"$SM"_sort.bam "gorilla_NC_011120.1:1-15446"  > $OUTDIR/"$SM"_sort_MT.bam
    samtools index $OUTDIR/"$SM"_sort_MT.bam

    #How many reads map to the mt genome?

    echo 'Reads mapping to the mt genome:'
    echo $(samtools view -c "$SM"_sort_MT.bam)
    
    # set up Read Group metadata for marking duplicates
    ID=$(samtools view $DATADIR/"$i" | awk 'NR==1 {print $1}' | cut -d ":" -f 1,2-4 | sed "s|":"|"."|g")
    PU=$(samtools view $DATADIR/"$i" |  awk 'NR==1 {print $1}' | cut -d ":" -f 2-4 | cut -d "@" -f 2 | sed "s|":"|"."|g")
    SM=$(echo "$name" | cut -d _ -f 1)
    PL=ILLUMINA
    LB="$SM"
    echo "Starting Picard"
    java -jar "$PICARD_ROOT"/picard.jar AddOrReplaceReadGroups \
            -I $OUTDIR/"$SM"_sort_MT.bam \
            -OUTPUT $OUTDIR/"$SM"_sort_MT_RG.bam \
            -ID "$ID" \
            -LB "$LB" \
            -PL ILLUMINA \
            -PU "$PU" \
            -SM "$SM"

    #I could use gatk MarkDuplicatesSpark here, but it removes way too many reads

    echo "Starting angsd"
    angsd -i $OUTDIR/"$SM"_sort_MT_RG.bam -out $OUTDIR/"$SM"_sort_MT_RG_angsd -doCounts 1 -doFasta 2
    
    #Rename headers
    pigz -p "$SLURM_CPUS_ON_NODE" -d "$SM"_sort_MT_RG_angsd.fa.gz
    sed -i "s/>gorilla_NC_011120.1/>$SM/g" "$SM"_sort_MT_RG_angsd.fa
    pigz -p "$SLURM_CPUS_ON_NODE" "$SM"_sort_MT_RG_angsd.fa
      
    cd $DATADIR || exit
    
    #Calculate mitochondrial genome coverage
    cov=$(samtools coverage $OUTDIR/"$SM"_sort_MT.bam -r gorilla_NC_011120.1| awk 'NR > 1{print $6}')
    nreads=$(samtools view -c $OUTDIR/"$SM"_sort_MT.bam)
    echo "$SM, $cov, $nreads" >> $OUTDIR/mt_coverage.csv
    
    #Remove intermediate files
    rm $OUTDIR/"${i%.bam}".fastq
    rm $OUTDIR/"$SM"_sort.bam
    rm $OUTDIR/"$SM"_sort.bam.bai
    rm $OUTDIR/"$SM"_sort_MT.bam
    rm $OUTDIR/"$SM"_sort_MT.bam.bai
    rm $OUTDIR/"$SM"_sort_MT_RG_angsd.arg
done
