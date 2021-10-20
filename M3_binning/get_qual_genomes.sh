#!/bin/bash -l

#Select genomes of medium or high quality based on CheckM output
#Medium-quality draft: >50% completion, <10% contamination
#High-quality draft: >90% completion, <5% contamination 
#According to MIMAG high quality genes must also fulfill the following criteria: Gaps span repetitive regions. Presence of the 23S, 16S, and 5S rRNA genes and at least 18 tRNAs
#Output formatted so it can be used with GTDB-Tk with the --batchfile option

#Create a list with all MAGs meeting the medium-quality draft criteria (but not high-quality draft criteria)
awk -v OFS='\t' -v path="${PWD}/" 'BEGIN {FS="\t"};(( $12>=50 && $12<90 ) && $13<10 ) || ( $12>=50 && ( $13<10 && $13>=5 )) {print path $1 ".fa", $1}' CheckM_results/all_checkm_results.tab > medium_quality_MAGs.txt

#Create a list with all MAGs meeting the high-quality draft criteria (only regarding coverage and contamination)
awk -v OFS='\t' -v path="${PWD}/" 'BEGIN {FS="\t"}; $12>=90 && $13<5 {print path $1 ".fa", $1}' CheckM_results/all_checkm_results.tab > high_quality_MAGs.txt
