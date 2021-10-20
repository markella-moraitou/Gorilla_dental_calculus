#!/bin/bash -l

#SBATCH -A snic2020-5-528
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 6:00:00
#SBATCH -J tax_analysis
#SBATCH --mail-type BEGIN
#SBATCH --mail-type END
#SBATCH --mail-type FAIL
#SBATCH --mail-user Markella.Moraitou.0437@student.uu.se

module load R_packages

Rscript 1_phyloseq.R > log1_phyloseq.txt
Rscript 2_FEAST.R > log2_FEAST.txt
Rscript 3_decontam.R > log3_decontam.txt
Rscript 4_abundance_filtering.R > log4_abundance_filtering.txt
Rscript 5_abundance_based_decontam.R > log5_abundance_based_decontam.txt
Rscript 6a_reference_based_decontam.R > log6a_reference_based_decontam.txt
Rscript 6b_reference_based_decontam.R > log6b_reference_based_decontam.txt
Rscript 7_permanova.R > log7_permanova.txt
Rscript 8_alpha.R > log8_alpha.txt
Rscript 10_ancom.R > log10_ancom.txt
