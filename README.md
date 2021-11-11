# Studying the gorilla oral microbiome using dental calculus from natural history museum specimens
Studying the oral microbiome of three gorilla subspecies using a shotgun metagenomic approach on dental calculus. The workflow is adapted from Brealey et al. 2020.
## Pre-processing of reads
* 2_polyGrem: Removal of PolyG tails 
* 3_adapterRem: Adapter removal and merging of paired-end reads 
* 4_barcodeTrim: Barcode trimming 
* 5_qualityFilt: Quality filtering 
* 6_deduplicate: Removal of duplicates 
* 7_phixRem: Removal of sequences mapping to the PhiX genome 
* 8_humanHostFilt: Removal of sequences mapping to the human and gorilla genome
## Host mtDNA analysis
* H1_mtgenomes: Obtaining mt genomes from dental calculus by mapping to a the gorilla mitochondrial genome reference 
* H2_haplotype_network: Plot haplotype network using downloaded complete mt genomes and dental calculus genomes > 80% complete. For dental calculus genomes < 80% and > 8% complete, find the most closely related genome from the haplotype network based on ANI values
## Taxonomic classification and analysis
* T1_kraken2: Taxonomic classification using Kraken2 
* T2_bracken: Re-estimation of taxa abundances using Bracken 
* T3_community-level: Rscripts for the analysis of the taxonomic classification output (decontamination, abundancefiltering, PERMANOVA, identification of differentially abundant taxa, etc.)
## Read-level decontamination and MapDamage
* RD3: The scripts in this directory are run as part of the taxonomic-level analysis scripts in T3_community-level, in order to identify contaminants based on damage patterns. Takes list of taxa that are found in both oral and contaminant lists (produced by T3_community-level/6a_reference_based_decontam.R) and are represented by a sufficient number of reads, downloads one genome per taxon and then runs mapDamage. The output is then processed with T3_community-level/6b_reference_based_decontam.R 
* RD1: First step of read-level decontamination of the FastQ files. Removes reads that were assigned to identified contaminants by Kraken2. 
* RD2: Second step of read-level decontamination of the FastQ files. Removed contaminant reads by mapping to a concatenated file of all contaminant genomes.
## Functional profiling and analysis
* F1_HUMAnN2: Functional profiling using HUMAnN2 
* F2_functional_stats: Rscripts for the analysis of the functional profiling output (identify subspecies associated biological processes, quantify contributions of different taxa to biological processes, etc.)
## Metagenomes-assembled genomes (MAGs)
### This part will probably be replaced by a new pipeline
* M1_MAGassemblies: Contig assembly using MEGAHIT 
* M2_readmapping: Mapping reads back to assemblies in preparation for binning 
* M3_binning: MAG binning script using MetaBAT 
* M4_MAGtaxonomy: Get MAG taxonomy 
* M5_samplestoMAGs: Map samples to MAGs to quantify coverage and abundance per subspecies 
* M6_MAGstats: Rscripts for the analysis of MAGs (phylogenetic trees, coverage and abundance per subspecies) 
* M7_mapdamage: Authentication using mapDamage 
* M8_dereplicateMAGs: Dereplication of MAGs
