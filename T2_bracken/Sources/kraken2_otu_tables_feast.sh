#!/bin/bash

#Adapted from /proj/SNIC_PROJECT/nobackup/JAELLE/DENTAL_CALCULUS_SECONDSCREEN_190219/SCRIPTS/kraken2_otu_tables.sh
#Run it after activating conda environment kraken-biom

#Include samples (sinks) to the table
cp ../*_m_* .

kraken-biom --max G --min S --fmt tsv -o species_table_feast.txt *_bracken_species.txt

#Remove the samples (sinks)
rm -I *_m_*
