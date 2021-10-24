#!/bin/bash -l

#SBATCH -A snic2020-5-528
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 30:00:00
#SBATCH -J samples_to_MAGs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Markella.Moraitou.0437@student.uu.se

#Map reads to medium and high quality MAG assemblies
#Count number of reads per sample mapping primarily to each MAG
#Using bbmap instead of bwa since need to build a ref for bwa
#Adapted from /proj/sllstore2017021/nobackup/ADRIAN/scripts/mag/mag_assembly_readmap_210111.sh

module load bioinfo-tools samtools bbmap QualiMap
echo $SLURM_JOB_NAME
echo $(module list)

ymd=$(date +%y%m%d)

#Read dir
RDIR=/proj/sllstore2017021/nobackup/MARKELLA/RD2_mapping/decontam_files
#Contig dir
CDIR=/proj/sllstore2017021/nobackup/MARKELLA/M3_binning

OUTDIR=/proj/sllstore2017021/nobackup/MARKELLA/M5_samplestoMAGs

noncontamref=/proj/sllstore2017021/nobackup/MARKELLA/RD2_mapping/refs/noncontam_references/noncontam_references.fna.gz

cd $OUTDIR

## Concatenate medium and high quality drafts
cat $CDIR/high_and_medium_quality_MAGs.txt | while read i
do
    #get label to use as identifier
    label=$(echo $i | awk '{print ">" $2 "_"}')
    #Add identifier and append to the concatenated reference fasta
    echo $i | awk '{print $1}' | xargs cat | sed "s/>/${label}/g" >> MAGs_concatenated_ref.fa
done

pigz -p $SLURM_CPUS_ON_NODE MAGs_concatenated_ref.fa

#Also add non contaminant references, to provide more options for mapping
cat MAGs_concatenated_ref.fa.gz $noncontamref > MAG_noncontam_ref.fa.gz

cd $RDIR

#Start mapping samples to the concatenated MAG fasta
find . -name "*fastq.gz" | while read i;
do
    echo ${i%_m_decontam.fastq.gz}
    #align reads to contigs with bbmap, generate coverage report
    bbwrap.sh ref=$OUTDIR/MAG_noncontam_ref.fa.gz in=$i mapper=bbmap out=$OUTDIR/${i%_m_decontam.fastq.gz}_mappedMAGs.sam  maxindel=80 pigz=t unpigz=t nodisk
    #view all mapped reads | sort ready for binning
    samtools view -Sb -F 4 -@ $SLURM_CPUS_ON_NODE $OUTDIR/${i%_m_decontam.fastq.gz}_mappedMAGs.sam | samtools sort -@ $SLURM_CPUS_ON_NODE -o $OUTDIR/${i%_m_decontam.fastq.gz}_mappedMAGs.sorted.bam
    cd $OUTDIR
done

#Start file where info will be logged
echo "sample, MAG, read_count, coverage_breadth" > $OUTDIR/log_map_readct_${ymd}.txt

#Get 20 groups of samples to be run in parallel
samples_per_thread=$( echo $(find . -name "*fastq.gz" | wc -l) / 20 | bc)

#For each group
for j in {1..20}
do
  #Get the list of samples in that group(e.g 1-6, 7-12 etc)
  ls *fastq.gz | awk -v samples=$samples_per_thread -v multiply=$j 'NR > samples*(multiply-1) && NR <= samples*multiply {print $0}' | while read i;
  do
    #Create temporary depth file including depth per position
    #By default, reads that have any of the flags UNMAP, SECONDARY, QCFAIL, or DUP set are skipped
    samtools depth -a $OUTDIR/${i%m_decontam.fastq.gz}mappedMAGs.sorted.bam > $OUTDIR/${i%m_decontam.fastq.gz}depth_temp.txt
    #Get read depth and coverage for each MAG
    cat $CDIR/high_and_medium_quality_MAGs.txt | awk '{print $2}' | while read k
    do
        cd $OUTDIR
        #Keep coverage data only for this specific MAG
        awk -v mag="${k}" '$1 ~ mag' ${i%m_decontam.fastq.gz}depth_temp.txt > ${i%m_decontam.fastq.gz}${k}_depth.txt
        #Prints some statistics:
        # Number of reads
        # Coverage breadth: adapted from https://sarahpenir.github.io/bioinformatics/awk/calculating-mapping-stats-from-a-bam-file-using-samtools-and-awk/
        echo ${i%.fastq.gz}, $k, $(samtools view -F 4 -F 256 ${i%m_decontam.fastq.gz}mappedMAGs.sorted.bam | awk -v mag="$k" '$3 ~ mag' | wc -l), $(awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' ${i%m_decontam.fastq.gz}${k}_depth.txt) >> log_map_readct_${ymd}.txt; 
        #rm $OUTDIR/${i%m_decontam.fastq.gz}${k}_depth.txt
        cd $RDIR
    done
  done &
done
wait

rm $OUTDIR/**depth_temp.txt
rm $OUTDIR/**_mappedMAGs.sam
