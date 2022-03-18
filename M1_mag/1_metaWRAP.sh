#!/bin/bash
#SBATCH -A SNIC_PROJECT_ID
#SBATCH -p node
#SBATCH -N 3
#SBATCH -t 0-24:00:00
#SBATCH -J magWRAP
#SBATCH -o slurm-%A_%a.out
#SBATCH --mail-type=FAIL

#### ASSEMBLY, BINNING, AND BINNING REFINEMENT
module load bioinfo-tools metaWRAP bbmap seqtk gnuparallel

# source /sw/apps/conda/latest/rackham/etc/profile.d/conda.sh
# conda activate /home/adrianf/project-folder/nobackup/ADRIAN/bin/bmtagger-env

ymd=$(date +%y%m%d)
echo "$ymd"

# set root for checkm
# echo $(readlink -f $CHECKM_DATA)
checkm data setRoot /home/adrianf/project-folder/nobackup/SHARED_REFERENCES/CHECKM/2015_01_16

# on a single node
threads=60
mem=384

### DIR with decontaminated READS
INTERLEAVED=/home/adrianf/project-folder/nobackup/MARKELLA/RD2_mapping/decontam_files
BASEDIR=/home/adrianf/project-folder/nobackup/ADRIAN/calculus/Gorilla/mag/metaWRAP/full-decontam-set
mkdir -p $BASEDIR

# generate a list of MAGs for the species of interest
# MAGLIST=$(grep $SSPECIES /home/adrianf/project-folder/nobackup/MARKELLA/M4_MAGtaxonomy/classify/hq.bac120.summary.tsv | awk '{print $1}' | cut -d . -f 1 )

### bin/ DIR for metaWRAP
SOFT=/home/adrianf/project-folder/nobackup/ADRIAN/bin/metaWRAP/bin/metawrap-scripts

ASSEMBLY_OUT=${BASEDIR}/MAG_assembly_metaWRAP
INITIAL_BINNING=${BASEDIR}/MAG_initialbinning_metaWRAP
BIN_REFINEMENT=${BASEDIR}/MAG_binsrefined_metaWRAP
BIN_REASSEMBLY=${BASEDIR}/MAG_binreassembly_metaWRAP
BIN_CLASSIFICATION=${BASEDIR}/MAG_binclassification_metaWRAP
QUANT_BINS=${BASEDIR}/MAG_quantbins_metaWRAP

DATADIR=/scratch/adrianf
mkdir -p $DATADIR/all

rm -f $DATADIR/*.fastq

echo "########################
#### START PIPELINE ####
########################"

# copy
find $INTERLEAVED/*_m_decontam.fastq.gz -size +4000 | parallel -j 47 cp {} $DATADIR/

# unzip
find $DATADIR/*_m_decontam.fastq.gz -size +4000 | parallel -j 47 "gunzip {}"

find $DATADIR/*_m_decontam.fastq -size +4000 | parallel -j 47 reformat.sh in={} out1={.}_1.fastq out2={.}_2.fastq tossbrokenreads addslash overwrite

cat $DATADIR/*_1.fastq >> $DATADIR/all/ALL_READS_1.fastq
cat $DATADIR/*_2.fastq >> $DATADIR/all/ALL_READS_2.fastq
echo "done cat-ing"

#### ASSEMBLY ####
echo "ASSEMBLY START"

rm -fr $ASSEMBLY_OUT
mkdir -p $ASSEMBLY_OUT

metawrap assembly -m $mem -t $threads --megahit -1 $DATADIR/all/ALL_READS_1.fastq -2 $DATADIR/all/ALL_READS_2.fastq -o $ASSEMBLY_OUT

echo "ASSEMBLY DONE"

#### BINNING ####
echo "INITIAL BINNING START"

rm -fr $INITIAL_BINNING
mkdir -p $INITIAL_BINNING/work_files

cp $ASSEMBLY_OUT/final_assembly.fasta $INITIAL_BINNING/work_files/assembly.fa

metawrap binning -o $INITIAL_BINNING -t "$threads" -a $ASSEMBLY_OUT --metabat2 --maxbin2 --concoct --interleaved $DATADIR/*_m_decontam.fastq

echo "INITIAL BINNING END"

#### BIN REFINEMENT ####
# -c = minimum completion
# -x = maximum contamination
# --keep-ambiguous	for contigs that end up in more than one bin, keep them in all bins default: keeps them only in the best bin
# --remove-ambiguous	for contigs that end up in more than one bin, remove them in all bins default: keeps them only in the best bin

echo "BIN REFINEMENT START"

rm -fr $BIN_REFINEMENT
mkdir $BIN_REFINEMENT

metawrap bin_refinement -o $BIN_REFINEMENT -t "$threads" -m $mem -A $INITIAL_BINNING/metabat2_bins -B $INITIAL_BINNING/maxbin2_bins -C $INITIAL_BINNING/concoct_bins -c 50 -x 10

echo "BIN REFINEMENT END"

#### QUANT BINS
# We would like to know how the extracted genomes are distributed across the samples, and in what abundances each bin is present in each sample. The Quant_bin module can give us this information. It used Salmon - a tool conventionally used for transcript quantitation - to estimate the abundance of each scaffold in each sample, and then computes the average bin abundances.

# NOTE: In order to run this module, it is highly recomended to use the non-reassembled bins (the bins produced by the Bin_Refinment module, not the Reassemble_bins module) and provide the entire non-binned assembly with the -a option. This will give more accurate bin abundances that are in context of the entire community.

echo "BIN QUANT START"
rm -fr $QUANT_BINS
mkdir $QUANT_BINS
metawrap quant_bins -b $BIN_REFINEMENT/metawrap_50_10_bins -o "$QUANT_BINS" -a "$ASSEMBLY_OUT"/final_assembly.fasta $DATADIR/*_[1,2].fastq

echo "BIN QUANT END"

#### RE-ASSEMBLY
# -c INT		minimum desired bin completion % (default=70)"
# -x INT		maximum desired bin contamination % (default=10)"
# -l INT		minimum contig length to be included in reassembly (default=500)"

echo "BIN RE-ASSEMBLY START"

rm -fr $BIN_REASSEMBLY
mkdir $BIN_REASSEMBLY

metawrap reassemble_bins -o $BIN_REASSEMBLY -1 $DATADIR/*_1.fastq -2 $DATADIR/*_2.fastq -t "$threads" -m $mem -c 50 -x 10 -l 500 -b $BIN_REFINEMENT/metawrap_50_10_bins 

echo "BIN RE-ASSEMBLY END"

#### Taxonomic Classification
# use gtdbtk here instead
# metawrap classify_bins -b $BIN_REASSEMBLY/reassembled_bins -o $BIN_CLASSIFICATION -t "$threads"

echo "
###							 ###
### END PIPELINE ###
###							 ###
"
