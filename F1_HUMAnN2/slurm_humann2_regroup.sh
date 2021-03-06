#!/bin/bash -l

#SBATCH -A SNIC_PROJECT_ID
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH -J humann2_regroup
#SBATCH --mail-type=ALL
#SBATCH --mail-user=USER_EMAIL

#Script for doing some additional manipulations to HUMANn2 output
conda activate /home/markmora/.conda/envs/humann2

#error loading humann2 straight, had to load metaphlan2 first
module load bioinfo-tools metaphlan2 biopython humann2/0.11.2

echo "$SLURM_JOB_NAME"
echo $(module list)

OUTDIR=F1_HUMAnN2
MAPPING_DIR=software/utility_mapping/ #tables downloaded with humann2_databases --download utility_mapping full $DIR

#postprocessing:
#regroup gene families into different functional categories (KEGG orthogroups and GO terms)
#renormalization
#joining tables

cd $OUTDIR || exit
for i in *genefamilies.tsv;
do
    #Get KO groups
    humann2_regroup_table --input "$i" --custom  $MAPPING_DIR/map_ko_uniref90.txt.gz --output "${i%genefamilies.tsv}"KOgroups.tsv
    humann2_renorm_table --input "${i%genefamilies.tsv}"KOgroups.tsv --output "${i%genefamilies.tsv}"KOgroups_cpm.tsv --units cpm
    #Normalize path abundances
    humann2_renorm_table --input "${i%genefamilies.tsv}"pathabundance.tsv --output "${i%genefamilies.tsv}"pathabundance_cpm.tsv --units cpm;
    #Get GO terms
    humann2_regroup_table --input "$i" --custom $MAPPING_DIR/map_go_uniref90.txt.gz --output "${i%genefamilies.tsv}"go.tsv
    humann2_renorm_table --input "${i%genefamilies.tsv}"go.tsv --output "${i%genefamilies.tsv}"go_cpm.tsv --units cpm
done

# join all normalized output files into one table
humann2_join_tables --input $OUTDIR --output all_KOgroups_cpm.tsv --file_name KOgroups_cpm
humann2_join_tables --input $OUTDIR --output all_GOterms_cpm.tsv --file_name go_cpm
humann2_join_tables --input $OUTDIR --output all_pathabund_cpm.tsv --file_name pathabundance_cpm
humann2_join_tables --input $OUTDIR --output all_pathcov.tsv --file_name pathcoverage

# get unstratified info
grep "|" -v $OUTDIR/all_KOgroups_cpm.tsv > $OUTDIR/all_KOgroups_cpm_unstratified.tsv
grep "|" -v $OUTDIR/all_GOterms_cpm.tsv > $OUTDIR/all_GOterms_cpm_unstratified.tsv
grep "|" -v $OUTDIR/all_pathabund_cpm.tsv > $OUTDIR/all_pathabund_cpm_unstratified.tsv
grep "|" -v $OUTDIR/all_pathcov.tsv > $OUTDIR/all_pathcov_unstratified.tsv


#get human_readable names
humann2_rename_table --input all_GOterms_cpm_unstratified.tsv -n go --output all_GOterms_cpm_unstratified_renamed.tsv
humann2_rename_table --input all_GOterms_cpm.tsv -n go --output all_GOterms_cpm_renamed.tsv

conda deactivate
