# Studying the gorilla oral microbiome using dental calculus from natural history museum specimens

Studying the oral microbiome of three gorilla subspecies using a shotgun metagenomic approach on dental calculus. The workflow is adapted from Brealey et al. 2020.

## Pre-processing of reads

Script for demultiplexing is not included here.

- 2_polyGrem: Removal of PolyG tails
- 3_adapterRem: Adapter removal and merging of paired-end reads
- 4_barcodeTrim: Barcode trimming
- 5_qualityFilt: Quality filtering
- 6_deduplicate: Removal of duplicates
- 7_phixRem: Removal of sequences mapping to the PhiX genome
- 8_humanHostFilt: Removal of sequences mapping to the human and gorilla genome

## Host mtDNA analysis

- H1_mtgenomes: Obtaining mt genomes from dental calculus by mapping to the gorilla mitochondrial genome reference. Script for obtaining mitochondrial sequences using a reference genome were adapted from Samantha Lopez Clinton's work.
- H2_haplotype_network: Plot haplotype network using downloaded complete mt genomes and dental calculus genomes > 80% complete. For dental calculus genomes < 80% and > 8% complete, find the most closely related genome from the haplotype network based on ANI values. Haplotype network script plot_hapnets.r was adapted from Zach Nolen's work. (not used in final report)
- H3_diagnostic_sites: Identifying and using mitochondrial diagnostic sites to verify species and subspecies assignments. Python and R scripts (DistanceToReferences.py, FindDiagnosticSites_edited.py, plot_readCoverage.r) were adapted from Axel Jensen's work.

## Taxonomic classification and analysis

- T1_kraken2: Taxonomic classification using Kraken2.
- T2_bracken: Re-estimation of taxa abundances using Bracken
- T3_community-level: Rscripts for the analysis of the taxonomic classification output (decontamination, abundancefiltering, PERMANOVA, identification of differentially abundant taxa, etc.)

## Read-level decontamination and MapDamage

- RD3: The scripts in this directory are run as part of the taxonomic-level analysis scripts in T3_community-level, in order to identify contaminants based on damage patterns. Takes list of taxa that are found in both oral and contaminant lists (produced by T3_community-level/6a_reference_based_decontam.R) and are represented by a sufficient number of reads, downloads one genome per taxon and then runs mapDamage. The output is then processed with T3_community-level/6b_reference_based_decontam.R
- RD1: First step of read-level decontamination of the FastQ files. Removes reads that were assigned to identified contaminants by Kraken2.
- RD2: Second step of read-level decontamination of the FastQ files. Removed contaminant reads by mapping to a concatenated file of all contaminant genomes.

## Functional profiling and analysis

- F1_HUMAnN2: Functional profiling using HUMAnN2
- F2_functional_stats: Rscripts for the analysis of the functional profiling output (identify subspecies-associated biological processes, quantify contributions of different taxa to biological processes, etc.)

## Metagenomes-assembled genomes (MAGs)

- M1_metawrap: assembly, binning, refinement, taxonomic classification, and phylogenies of MAGs
- M2_isolation-sources: Isolation source information from multiple databases, categorization and visualization of isolation sources for each MAG
- M3_trees: scripts for producing core-gene trees for select oral taxa, using Phylophlan.
- M4_pangenome: scripts for running pangenomic analyses

## Analyses of dietary components

- D1_reads4Diet: Preprocessing for diet profiling, uses Kraken2 output to extract the reads that have been classified as bacterial, viral or archeal
- D2_kraken2_full_db: Run of Kraken2 on the NCBI full-nt database, with genus and species level assignments retained.
- D3_diet_stats: Statistical analyses with kraken output of dietary components. Constructing of reference database of dietary components. Preprocessing and filtering of abundance table for dietary taxa. ANCOM and PERMANOVAs on the filtered dataset of dietary taxa.
- D4_diet_damage: extracts reads used to assign a given dietary taxon, maps these reads to a reference genome, and runs mapDamage to get the misincorporation rates.

## Misc. Analyses, plots, tables, other stats

- frag_length_vs_samp_age: script for plotting the relationship between the average fragment length for each sample and the sample age
- microbial_vs_diet_MRM: script for plotting the relationship between microbial and dietary distances
- plots: scripts for producing most figures presented in the report and the supplement
- tables_and_stats: script producing most information provided in supplementary tables
