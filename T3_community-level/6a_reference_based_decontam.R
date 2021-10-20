#Start logging
sink(file = "log6a_reference_based_decontam.txt")

#packages
library(dplyr)
library(tidyverse)
library(phyloseq)
library(ggpubr)
library(reshape2)
library(readr)
library(phyloseq)
library(tibble)
library(tidytree)
library(ggtree)
library(FEAST)
library(MutationalPatterns)

load(".RData")

#Community-level taxonomic analysis - Script 6a (before running map damage)
#Removal of taxa based on references (using mapDamage to sort out ambiguities)

#### Compare dataset with common contaminants ####

#
print("How many of the species in our data belong to known contaminants")
length(which(tax_table(spe_data_envrem)[,7] %in% contaminant_list))

print("How many genera do they belong to?")
length(intersect(tax_table(spe_data_envrem)[,7], contaminant_list))

#Species that are found in the contaminant list but not the oral lists
unambiguous_contam <- taxa_names(spe_data_envrem)[which(tax_table(spe_data_envrem)[,7] %in% contaminant_list &
                                                          !(tax_table(spe_data_envrem)[,7] %in% append(str_remove(core_micr$Taxon, " .*"), homd$Genus)))]
#Note: core_micr contains some family level taxonomic classifications - these are being ignored

#
print("How many of the species are unambigous contaminants (they only appear in contaminant list)")
length(unambiguous_contam)

print("How many different genera do they belong to?")
tax_table(spe_data_envrem)[which(taxa_names(spe_data_envrem) %in% unambiguous_contam),7] %>%
  unique %>% length


#Species that are found in both the contaminant lists and the oral list
ambiguous_contam <- taxa_names(spe_data_envrem)[which(tax_table(spe_data_envrem)[,7] %in% contaminant_list &
                                                        tax_table(spe_data_envrem)[,7] %in% append(str_remove(core_micr$Taxon, " .*"), homd$Genus))]
#
print("How many of the species are found in both contaminant and oral lists?")
length(ambiguous_contam)

print("How many different genera do they belong to?")
tax_table(spe_data_envrem)[which(taxa_names(spe_data_envrem) %in% ambiguous_contam),7] %>%
  unique %>% length

#### Check number of reads per taxon ####
#Use number of reads instead of taxon abundance - therefore using Kraken2 output
spe_kraken <- read.table(file="/proj/sllstore2017021/nobackup/MARKELLA/T1_kraken2/otu_table_kraken.txt",
                         header = TRUE, comment.char = "", skip=1, sep="\t")

spe_kraken$X.OTU.ID <- as.character(spe_kraken$X.OTU.ID)
colnames(spe_kraken) <- str_remove(colnames(spe_kraken), "_m_kraken2_report")

#Get read count per sample for all ambiguous taxa
ambiguous_abund <-
  spe_kraken[which(spe_kraken$X.OTU.ID %in% ambiguous_contam),] %>%
  melt(id.vars="X.OTU.ID")

colnames(ambiguous_abund) <- c("taxon", "sample", "read_count")

rownames(ambiguous_abund) <- NULL

#Rank by abundance and then keep the sample for which the abundance is highest
#After removing the "bad" samples
ambiguous_abund <-
  ambiguous_abund %>% group_by(taxon) %>% filter(!(sample %in% bad_samples))%>% mutate(max_read_count = max(read_count)) %>% 
  ungroup %>% filter(read_count==max_read_count, read_count > 5000)

#Order by abundance
ambiguous_abund <- ambiguous_abund[order(ambiguous_abund$read_count, decreasing=TRUE),]

#Get taxon name
ambiguous_abund$taxon_name <- taxonomy_species[match(ambiguous_abund$taxon, 
                                                     rownames(taxonomy_species)),8]

#Save the list - these taxa will be checked with mapDamage on UPPMAX
#include sample name also (this will be the sample that will be used for mapping)
write.table(ambiguous_abund %>% select(sample, taxon, taxon_name), file="/proj/sllstore2017021/nobackup/MARKELLA/RD3_mapdamage4ambiguoustaxa/top_ambiguous_taxa.txt",
            col.names = FALSE, quote = FALSE, row.names = FALSE)

sessionInfo()
#Stop logging
sink(file=NULL)
        
save.image()
#### Run map damage for top ambiguous taxa ####
#system("cd /proj/sllstore2017021/nobackup/MARKELLA/RD3_mapdamage4ambiguoustaxa; sbatch slurm_get_ambigtaxa_genomes.sh)
