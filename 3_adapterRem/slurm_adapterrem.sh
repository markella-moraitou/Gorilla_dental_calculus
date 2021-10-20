#!/bin/bash -l

#SBATCH -A snic2020-5-528
#SBATCH -p node
#SBATCH -N 1-1
#SBATCH -t 15:00:00
#SBATCH -J AdapterRemoval
#SBATCH --mail-type BEGIN
#SBATCH --mail-type END
#SBATCH --mail-type FAIL
#SBATCH --mail-user Markella.Moraitou.0437@student.uu.se

# analysis step 2: Adapter Removal
# adapted from: /proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_SECONDSCREEN_190219/SCRIPTS/adrm_run1_190501.sh
# collapsing reads
# trims @ < 30 and Ns (default mismatch rate mm 3)
# removes sequences < 30 bp in length

#load modules
module load bioinfo-tools
module load AdapterRemoval/2.2.2
module load bbmap
module load pigz

echo $SLURM_JOB_NAME
echo $(module list)
echo "Number of cores ${SLURM_CPUS_ON_NODE}"

#define shortcuts

OUTDIR=/proj/sllstore2017021/nobackup/MARKELLA/3_adapterRem
DATADIR=/proj/sllstore2017021/nobackup/MARKELLA/2_polyGrem

# adapters - should match the END of fwd and rev reads, respectively
#adapter1 is rev comp of my recorded P7 full-length adapter
#adapter2 is rev comp of my recorded P5 full-length adapter
adapter1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNATCTCGTATGCCGTCTTCTGCTTG"
adapter2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT"

#Starting with DC2, DC3 and Jena samples for which polyG removal was performed.
#All contain the same set of adapters, except for the samples from Jena, which contain no adapters so nothing will be trimmed.

cd $DATADIR

find . -name '*F_concat_Gtrimmed.fastq.gz' -size +0 | while read file 
do
    cd $OUTDIR
    #First, re-pair reads
    repair.sh in=$DATADIR/$file in2=$DATADIR/${file%F_concat_Gtrimmed.fastq.gz}R_concat_Gtrimmed.fastq.gz out=$OUTDIR/${file%concat_Gtrimmed.fastq.gz}nosingletons.fastq out2=$OUTDIR/${file%F_concat_Gtrimmed.fastq.gz}R_nosingletons.fastq outsingle=$OUTDIR/${file%F_concat_Gtrimmed.fastq.gz}singletons.fastq overwrite=t
    #Next, remove adapters
    AdapterRemoval --file1 ${file%concat_Gtrimmed.fastq.gz}nosingletons.fastq --file2 ${file%F_concat_Gtrimmed.fastq.gz}R_nosingletons.fastq \
    --basename ${file%F_concat_Gtrimmed.fastq.gz}m --adapter1 $adapter1 --adapter2 $adapter2 --trimns --trimqualities --minquality 30 \
    --minlength 30 --gzip --mm 3 --minalignmentlength 11 --collapse --threads $SLURM_CPUS_ON_NODE
done

cd $OUTDIR 

# repair output is not compressed, compress with gzip in parallel (pigz)
ls *fastq | while read i
do
    pigz -p $SLURM_CPUS_ON_NODE $i
done

#DC1 samples also have the same set of adapters but are located in a different directory (polyG removal was not needed since they were sequenced on Illumina HiSeq)

SEQDIR=/proj/sllstore2017021/nobackup/JAELLE

find $SEQDIR/DENTAL_CALCULUS_FIRSTSCREEN_**/P1_** -name "G*wash*.fastq.gz" | while read i
do
    name=${i##/**/} #remove path
    echo $name
done | sort | uniq | while read n #iterate through list of unique names
do #concatenate them when they have the same name (save in polyGrem directory for consistency)
    ls -d -1 $SEQDIR/DENTAL_CALCULUS_FIRSTSCREEN_**/P1_**/$n | xargs cat > $DATADIR/${n%.fastq.gz}_concat.fastq.gz
done

cd $DATADIR

find . -name "G*wash*F_concat.fastq.gz" | while read i
do
    cd $OUTDIR
    AdapterRemoval --file1 $DATADIR/$i --file2 $DATADIR/${i%F_concat.fastq.gz}R_concat.fastq.gz --basename ${i%F_concat.fastq.gz}m --adapter1 $adapter1 \
    --adapter2 $adapter2 --trimns --trimqualities --minquality 30 --minlength 30 --gzip --mm 3 --minalignmentlength 11 --collapse --threads $SLURM_CPUS_ON_NODE
done

#The museum controls have the standard Illumina adapter and are also located in a different directory
#identified by running --identify-adapter on adapterRemoval (/proj/sllstore2017021/nobackup/MARKELLA/contamination_controls/less slurm-18000861.out)

SEQDIR=/proj/sllstore2017021/nobackup/MARKELLA/contamination_controls

cd $SEQDIR

find . -name 'ERR2503700_F.fastq.gz' -size +0 | while read i 
do
    cd $OUTDIR
    adapter1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNATCTCGTATGCCGTCTTCTGCTTG" 
    adapter2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
    AdapterRemoval --file1 $SEQDIR/$i --file2 $SEQDIR/${i%F.fastq.gz}R.fastq.gz --basename ${i%F.fastq.gz}m --adapter1 $adapter1 \
    --adapter2 $adapter2 --trimns --trimqualities --minquality 30 \
    --minlength 30 --gzip --mm 3 --minalignmentlength 11 --collapse --threads $SLURM_CPUS_ON_NODE
done

find . -name 'ERR2868193_F.fastq.gz' -size +0 | while read i 
do
    cd $OUTDIR
    adapter1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG"
    adapter2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
    AdapterRemoval --file1 $SEQDIR/$i --file2 $SEQDIR/${i%F.fastq.gz}R.fastq.gz --basename ${i%F.fastq.gz}m --adapter1 $adapter1 \
    --adapter2 $adapter2 --trimns --trimqualities --minquality 30 \
    --minlength 30 --gzip --mm 3 --minalignmentlength 11 --collapse --threads $SLURM_CPUS_ON_NODE
done

#Compress output
ls *fastq | while read i
do
    pigz -p $SLURM_CPUS_ON_NODE $i
done