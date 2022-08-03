#!/bin/bash -l

#SBATCH -A SNIC_PROJECT_ID
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 0-1:00:00
#SBATCH -J panaroo
#SBATCH -o SNIC_PROJECT/logs/slurm-%A_%a.out
#SBATCH --array=1-3%3
#SBATCH --mail-type=FAIL

# install panaroo
source /sw/apps/conda/latest/rackham/etc/profile.d/conda.sh
# mamba create -p SNIC_PROJECT/ADRIAN/bin/panaroo-env -c conda-forge -c bioconda -c defaults panaroo

# activate environment
conda activate SNIC_PROJECT/ADRIAN/bin/panaroo-env

# load modules
module load bioinfo-tools prokka

species_list="Neisseria
Rothia
Veillonella"

# iterate over species_list using slurm array
SPECIES=$(echo "$species_list" | cut -f "$SLURM_ARRAY_TASK_ID" -d ' ')

threads=20

BASEDIR=SNIC_PROJECT/ADRIAN/calculus/Gorilla/mag/panaroo/$SPECIES
INDIR=$BASEDIR/input
REANNOTATED=$BASEDIR/reannotation
OUTDIR=$BASEDIR/results
PROKKA=SNIC_PROJECT/ADRIAN/calculus/Gorilla/mag/prokka/$SPECIES
ANNOTATED=$BASEDIR/annotated

# run panaroo
panaroo -t $threads -i "$ANNOTATED"/*/*.gff -o "$OUTDIR" --clean-mode strict --remove-invalid-genes --merge_paralogs
