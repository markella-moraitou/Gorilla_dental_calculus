#!/bin/bash -l

#SBATCH -A SNIC_PROJECT_ID
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 10:00:00
#SBATCH -J HumanHostRem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=USER_EMAIL

# Filtering step 4: Map to host and human genomes
# Evaluate host content and human contamination
# Save mapped reads for future analysis
# Save unmapped reads for microbial analysis
# Adapted from /proj/SNIC_PROJECT/nobackup/JAELLE/DENTAL_CALCULUS_SECONDSCREEN_190219/SCRIPTS/rmhost_merged_190606.sh

module load bioinfo-tools samtools bwa BEDTools

# save versions of tools
echo $(module list)

REFDIR=/proj/SNIC_PROJECT/nobackup/JAELLE/REFERENCES
DATADIR=7_phixRem
UNMAPDIR=8_humanHostFilt/unmapped
MAPDIR=8_humanHostFilt/mapped

mkdir $UNMAPDIR
mkdir $MAPDIR

echo "$SLURM_JOB_NAME"
cd $DATADIR || exit

echo "sample, unmapped read #, mapped read #, host-sp read #, human-sp read #" >> $UNMAPDIR/readcount_hostmapping.txt

find . -maxdepth 1 -name "*fastq.gz" | while read i
do
    # Save sam for whole alignment
    host="${i:2:2}"
        if [[ "$host" = "G0" || "$host" = "Gb" || "$host" = "ER" ]];
        then
            INDEX=$REFDIR/gorilla_gorilla_human/calculushost
        elif [[ "$host" = "BE" || "$host" = "BL" || "$host" = "BS" || "$host"="LI" || "$host"="EX" ]] #assume no reference, only human (e.g. blanks)
        then
            INDEX=$REFDIR/human_genome_new/Hsap
        else
            INDEX=$REFDIR/gorilla_gorilla_human/calculushost #Jena gorillas (harder to select the file names)
        fi
        if [[ -f $UNMAPDIR/${i%rmphix.fastq.gz}host_unmapped.fastq.gz ]]
        then
            continue
        else
            bwa mem -t 10 $INDEX "$i" > $MAPDIR/"${i%rmphix.fastq.gz}"temp.sam
    
    # save mapped bam of host reads after sorting and index
            samtools view -Sb -F 4 -@ 10 $MAPDIR/"${i%rmphix.fastq.gz}"temp.sam | samtools sort -o $MAPDIR/"${i%rmphix.fastq.gz}"host_mapped_sorted.bam -@ 10 -
            samtools index $MAPDIR/"${i%rmphix.fastq.gz}"host_mapped_sorted.bam

    # save unmapped reads and get fastq.gz reads
            samtools view -Sb -f 4 -@ 10 $MAPDIR/"${i%rmphix.fastq.gz}"temp.sam > $UNMAPDIR/"${i%rmphix.fastq.gz}"host_unmapped.bam
            bedtools bamtofastq -i $UNMAPDIR/"${i%rmphix.fastq.gz}"host_unmapped.bam -fq $UNMAPDIR/"${i%rmphix.fastq.gz}"host_unmapped.fastq
            pigz -p 10 $UNMAPDIR/"${i%rmphix.fastq.gz}"host_unmapped.fastq

    # get some stats
            mapct=$(samtools view -c $MAPDIR/"${i%rmphix.fastq.gz}"host_mapped_sorted.bam)
            humanct=$(samtools view -F 4 -F 256 -q 30 $MAPDIR/"${i%rmphix.fastq.gz}"host_mapped_sorted.bam | awk '$3 ~ /^human/' | wc -l)
            unmapct=$(samtools view -c $UNMAPDIR/"${i%rmphix.fastq.gz}"host_unmapped.bam)
            hostct=$(samtools view -F 4 -F 256 -q 30 $MAPDIR/"${i%rmphix.fastq.gz}"host_mapped_sorted.bam | awk '$3 ~ /^gorilla/' | wc -l)

            echo "${i%_rmphix.fastq.gz}, $unmapct, $mapct, $hostct, $humanct" >> $UNMAPDIR/readcount_hostmapping.txt

    # remove temp sam file
            rm $MAPDIR/"${i%rmphix.fastq.gz}"temp.sam;
        fi
done
