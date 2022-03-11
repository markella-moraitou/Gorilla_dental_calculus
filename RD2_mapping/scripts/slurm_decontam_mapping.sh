#!/bin/bash -l

#SBATCH -A SNIC_PROJECT_ID
#SBATCH -p core
#SBATCH -n 5
#SBATCH -t 48:00:00
#SBATCH -J decontam_mapping
#SBATCH --mail-type=ALL
#SBATCH --mail-user=USER_EMAIL

# Map to contaminant genomes
# Save save unmapped reads for downstream analysis
# Based on /proj/SNIC_PROJECT/nobackup/JAELLE/DENTAL_CALCULUS_SECONDSCREEN_190219/SCRIPTS/rmhost_merged_190606.sh,
# but adapted for filter decontamination

module load bioinfo-tools samtools bwa BEDTools

# save versions of tools
echo $(module list)

REFDIR=RD2_mapping/refs
DATADIR=RD1_krakentools
OUTDIR=RD2_mapping/decontam_files

echo "$SLURM_JOB_NAME"

mkdir $OUTDIR

cd $DATADIR || exit

ymd=$(date +%y%m%d)

echo "sample, unmapped, mapped, contam_mapped, noncontam_mapped, decontam_retained" > $OUTDIR/readcount_decontam_"${ymd}".txt

find . -maxdepth 1 -name "*decontam1.fastq.gz" | while read i
do
    if [[ -f $OUTDIR/${i%1.fastq.gz}2.fastq.gz ]]
    then
        continue
    #process only if the file name matches with the list of retained samples or if it is a control or blank
    else
        echo "${i%_m_decontam1.fastq.gz}"
        bwa mem -t "$SLURM_CPUS_ON_NODE" $REFDIR/gorilla_dencalc_reference "$i" > $OUTDIR/"${i%decontam1.fastq.gz}"temp.sam   
        cd $OUTDIR || exit
        #save unmapped reads and get fastq file
        samtools view -Sb -f 4 -@ "$SLURM_CPUS_ON_NODE" "${i%decontam1.fastq.gz}"temp.sam > "${i%decontam1.fastq.gz}"unmapped.bam
        bedtools bamtofastq -i "${i%decontam1.fastq.gz}"unmapped.bam -fq "${i%decontam1.fastq.gz}"unmapped.fastq
        #save mapped reads      
        samtools view -Sb -F 4 -@ "$SLURM_CPUS_ON_NODE" "${i%decontam1.fastq.gz}"temp.sam > "${i%decontam1.fastq.gz}"mapped.bam
        #extract names of reads that primarily map to non-contaminants
        samtools view -F 4 -F 256 -@ "$SLURM_CPUS_ON_NODE" "${i%decontam1.fastq.gz}"mapped.bam | awk '$3 ~ /^endogenous/ {print $1}' > "${i%decontam1.fastq.gz}"endogenous_mapped.txt
        #extract those reads into a separate bam file
        zcat $DATADIR/"$i" | grep -A 3 --no-group-separator -Ff "${i%decontam1.fastq.gz}"endogenous_mapped.txt > "${i%decontam1.fastq.gz}"endogenous_mapped.fastq
        #merge unmapped reads and reads that mapped to noncontaminants to produce the decontaminated file
        cat "${i%decontam1.fastq.gz}"endogenous_mapped.fastq "${i%decontam1.fastq.gz}"unmapped.fastq > "${i%decontam1.fastq.gz}"decontam.fastq
        
        #get some stats
        unmapct=$(samtools view -c -@ "$SLURM_CPUS_ON_NODE" $OUTDIR/"${i%decontam1.fastq.gz}"unmapped.bam)
        mapct=$(samtools view -c -@ "$SLURM_CPUS_ON_NODE" $OUTDIR/"${i%decontam1.fastq.gz}"mapped.bam)
        contamct=$(samtools view -F 4 -F 256 -@ "$SLURM_CPUS_ON_NODE" $OUTDIR/"${i%decontam1.fastq.gz}"mapped.bam | awk '$3 ~ /^contam/'  | wc -l)
        noncontamct=$(samtools view -F 4 -F 256 -@ "$SLURM_CPUS_ON_NODE" $OUTDIR/"${i%decontam1.fastq.gz}"mapped.bam | awk '$3 ~ /^endogenous/'  | wc -l)
        decontam_retained=$(grep "@" "${i%decontam1.fastq.gz}"decontam.fastq | wc -l)
        echo "${i%_m_decontam1.fastq.gz}, $unmapct, $mapct, $contamct, $noncontamct, $decontam_retained" >> $OUTDIR/readcount_decontam_"${ymd}".txt

        #compress output
        pigz -p "$SLURM_CPUS_ON_NODE" "${i%decontam1.fastq.gz}"decontam.fastq
        #remove unnecessary file
        rm $OUTDIR/"${i%decontam1.fastq.gz}"unmapped.bam
        rm $OUTDIR/"${i%decontam1.fastq.gz}"mapped.bam
        rm $OUTDIR/"${i%decontam1.fastq.gz}"unmapped.fastq
        rm $OUTDIR/"${i%decontam1.fastq.gz}"endogenous_mapped.fastq
        rm $OUTDIR/"${i%decontam1.fastq.gz}"temp.sam
        cd $DATADIR || exit
    fi
done
