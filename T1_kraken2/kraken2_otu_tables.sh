#!/bin/bash

#Adapted from SNIC_PROJECT/JAELLE/DENTAL_CALCULUS_SECONDSCREEN_190219/SCRIPTS/kraken2_otu_tables.sh

kraken-biom --max G --min S --fmt tsv -o otu_table_kraken.txt *report.txt
