#!/bin/bash -l
#SBATCH -A SNIC_PROJECT_ID
#SBATCH -p core -n 2
#SBATCH -t 1:30:00
#SBATCH -J magtree
#SBATCH --mail-type=FAIL
#SBATCH -o SNIC_PROJECT/logs/slurm-%A_%a.out

module load bioinfo-tools PhyloPhlAn diamond MAFFT trimAl FastTree

SPECIES="Veillonella"
SSPECIES="dispar"

threads=$SLURM_NTASKS
# threads=2

BASEDIR=SNIC_PROJECT/ADRIAN/calculus/Gorilla/mag/metaWRAP/decontam-set/phylophlan/$SPECIES

INPUT=${BASEDIR}/input
OUTDIR=$BASEDIR/output
TREEDIR=${OUTDIR}/tree

rm -fr $OUTDIR/input*/*.tre

mkdir -p "$INPUT"
mkdir -p "$OUTDIR/phylophlan_configs"
mkdir -p "$TREEDIR"


# Automatically retrieve the set of core UniRef90 proteins for the species:
# had to use a well annotated complete reference genome here
# phylophlan_setup_database \
#     -g s__${SPECIES}_$SSPECIES \
#     -o $BASEDIR \
#     --verbose 2>&1 | tee $OUTDIR/phylophlan_setup_database.log

# files have to have ".fna" suffix!
rename ".fa" ".fna" $INPUT/*.fa

cd $BASEDIR || exit

# configuration file generated using the following command:
# rm -f ${OUTDIR}/references_config.cfg
# phylophlan_write_config_file \
#     -o $OUTDIR/references_config.cfg \
#     -d a \
#     --db_aa diamond \
#     --map_aa diamond \
#     --map_dna diamond \
#     --msa mafft \
#     --overwrite \
#     --trim trimal \
#     --tree1 fasttree
    # --tree2 raxml

# bug with uppmax installation, raxml not the right version?
# printf '%s\n%s\n' '[tree2]' 'input = -s' 'database = -t' 'threads = -T' 'output = -n' "params = -p 1989 -m PROTCATLG" 'version = -v' 'command_line = #program_name# #params# #threads# #database# #output_path# #input# #output#' 'program_name = raxmlHPC-PTHREADS-AVX' 'output_path = -w' >> ${OUTDIR}/references_config.cfg

# and we can now build a phylogenetic tree of all our genomes:
# make sure to specify where the database file ends up!

phylophlan \
    -i ${INPUT} \
    --output_folder ${OUTDIR} \
    -d s__${SPECIES}_$SSPECIES \
    -t a \
    -f ${OUTDIR}/references_config.cfg \
    --nproc "$threads" \
    --diversity low \
    --fast \
    --verbose 2>&1 | tee $OUTDIR/phylophlan_$SPECIES.log

# mv ${OUTDIR}/tree/bins_${SPECIES}/${SPECIES}_RAxML_bestTree.bins_refined.tre ${OUTDIR}/tree/bins_${SPECIES}/${SPECIES}_RAxML_bestTree.bins_refined.tre
