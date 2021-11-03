# Studying the gorilla oral microbiome using dental calculus from natural history museum specimens
Studying the oral microbiome of three gorilla subspecies using a shotgun metagenomic approach on dental calculus.
The workflow is adapted from Brealey et al. 2020.

## Pre-processing of reads
* 2_polyGrem: Removal of PolyG tails
* 3_adapterRem: Adapter removal and merging of paired-end reads
* 4_barcodeTrim: Barcode trimming
* 5_qualityFilt: Quality filtering
* 6_deduplicate: Removal of duplicates
* 7_phixRem: Removal of sequences mapping to the PhiX genome
* 8_humanHostFilt: Removal of sequences mapping to the human and gorilla genome

## Taxonomic classification and analysis
* T1_kraken2: Taxonomic classification using Kraken2
* T2_bracken: Re-estimation of taxa abundances using Bracken
* T3_community-level: Rscript for the analysis of the taxonomic classification output (decontamination, abundancefiltering, PERMANOVA, identification of differentially abundant taxa etc.)

