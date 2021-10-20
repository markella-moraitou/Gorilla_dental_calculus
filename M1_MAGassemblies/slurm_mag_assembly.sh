#!/bin/bash
#SBATCH -A snic2020-5-528
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 2-00:00:00
#SBATCH -J mag_assembly

# MAG assembly of metagenomes with MEGAHIT2
# Going with defaults but note tips online say
# increase --k-min to e.g. 27 for ultra complex data like soil
# decrease --k-step to e.g. 10 low-coverage data
# --kmin-1pass mode is more memory efficient for ultra low-depth data like soil
# Adapted from /proj/sllstore2017021/nobackup/ADRIAN/scripts/mag/megahit_stats_210111.sh

# also generate some basic stats
module load bioinfo-tools megahit
echo $SLURM_JOB_NAME
echo $(module list)

ymd=$(date +%y%m%d)
threads=$SLURM_CPUS_ON_NODE
DATADIR=/proj/sllstore2017021/nobackup/MARKELLA/RD2_mapping/decontam_files
OUTDIR=/proj/sllstore2017021/nobackup/MARKELLA/M1_MAGassemblies

echo "sample, contig_n, min_size, max_size, average_size, contig1.5kb_n, contig20kb_n" > $OUTDIR/log_mag_assembly_$ymd.txt

# run megahit
cd $DATADIR
find . -maxdepth 1 -name "*fastq.gz" | while read i
do
    cd $OUTDIR
    name=${i%_m_decontam.fastq.gz}
    megahit -r $DATADIR/$i -t $threads -o $name --out-prefix $name
    fin=${name}/${name}.contigs.fa
    ct=$(grep "^>" -c $fin)
    g1kb=$(awk '-F[=]' '$4 > 1500' $fin | wc -l)
    g20kb=$(awk '-F[=]' '$4 > 20000' $fin | wc -l)
    cmax=$(grep "^>" $fin | awk '-F[=]' '{print $4}' | sort -nr | head -1)
    cmin=$(grep "^>" $fin | awk '-F[=]' '{print $4}' | sort -n | head -1)
    cmean=$(grep "^>" $fin | awk '-F[=]' '{sum+=$4} END { print sum/NR}')
    echo "${name}, $ct, $cmin, $cmax, $cmean, $g1kb, $g20kb" >> $OUTDIR/log_mag_assembly_$ymd.txt
done
