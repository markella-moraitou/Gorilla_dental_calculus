#!/bin/bash

#Adapted from SNIC_PROJECT/JAELLE/DENTAL_CALCULUS_SECONDSCREEN_190219/SCRIPTS/kraken2_otu_tables.sh
#Retaining assignments at the genus and species level
kraken-biom --max G --min S --fmt tsv -o otu_table_kraken_fulldb.txt *report.txt
