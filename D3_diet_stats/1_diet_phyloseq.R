#Start logging
sink(file = "log1_diet_phyloseq.txt")

#load packages
library(readxl)
library(phyloseq)
library(taxize)
library(curl)
library(compositions)
library(stringr)
library(ggplot2)
library(microbiome)
library(tidyverse)

load("T3_community-level/.RData")

#Diet analysis - Script 1
#Set up phyloseq object for analysis

#Import kraken-biom otu table for eukaryotic taxa into R
#Based on Kraken2 output instead of Bracken
euk_table <- read.table("D2_kraken2_full_db/otu_table_kraken_fulldb.txt", skip=1, sep="\t", comment.char="", header=T)

#Modifications on otu table

print("There are a few samples/blanks with no taxa identified. Those will be removed")
colnames(euk_table)[which(colSums(euk_table)==0)]

euk_table <- euk_table[,which(colSums(euk_table)!=0)] #remove columns with sum 0
euk_table <- t(euk_table) #transpose table. Vegdist documentation refers to columns as species and rows as sites.
colnames(euk_table)<-as.character(euk_table[1,]) #set first row as column (taxa) names
euk_table <- euk_table[-1,] #remove original first row

#Shorten sample names
rownames(euk_table) <- str_remove(rownames(euk_table), "_m_bact_arch_vir_removed_kraken2_report")

#transpose again to use phyloseq
euk_table=t(euk_table) 

#### Create Taxonomic Matrix ####

#Download taxonomy for every taxon
#taxonomy_euk_original <- list()
#for (i in 1:nrow(euk_table)){
#  taxon <- rownames(euk_table)[i]
#  taxonomy_euk_original[taxon] <- classification(sci_id=taxon, db="ncbi")
#}

#saveRDS(taxonomy_euk_original, file = "D3_diet_stats/taxonomy_euk_original_kraken.rds")
taxonomy_euk_original <- readRDS(file = "D3_diet_stats/taxonomy_euk_original_kraken.rds")

taxonomy_euk <- matrix(nrow = nrow(euk_table), ncol = 9)
colnames(taxonomy_euk) <- c("taxonID","superkingdom", "clade", "phylum", "class", "order", "family", "genus", "species")
#for loop to get a table format
for (i in c(1:nrow(euk_table)) ){
  taxon=rownames(euk_table)[i]
  classification <- taxonomy_euk_original[[taxon]]
   for (j in 2:length(colnames(taxonomy_euk))){ #for every taxonomic ranking in the matrix's header
    for (k in 1:length(classification$rank)){ #find the field in the classification object that matches
      if (as.character(classification$rank[k]) == as.character(colnames(taxonomy_euk)[j])){
        taxonomy_euk[i,j] = classification$name[k] #use name of ranking to fill in the matrix
        taxonomy_euk[i,1] = taxon
      }
    }
  }
}

taxonomy_euk <- as.data.frame(taxonomy_euk)

#Check if there are differences between the taxa in species_table and taxonomy_species
unwantedeuk <- taxonomy_euk %>% filter(!(taxonID %in% row.names(euk_table))) %>% pull(taxonID) %>% unique

#remove the unwanted taxa
taxonomy_euk=taxonomy_euk[! (taxonomy_euk$taxonID %in% unwantedeuk),]

missingeuk=unique(row.names(euk_table)[! row.names(euk_table) %in% taxonomy_euk$taxonID])

rownames(taxonomy_euk) <- taxonomy_euk[,1]
taxonomy_euk <- taxonomy_euk[,-1]

#Turn to matrix
taxonomy_euk <- as.matrix(taxonomy_euk)

#### Create phyloseq objects ####

#create phyloseq objects from: 
#sample_data object based on metadata
#otu_table based on species_table and species_table_norm respectively
euk_data <- phyloseq(otu_table(euk_table, taxa_are_rows = TRUE), sample_data(metadata), tax_table(taxonomy_euk))

#Get taxa richness data
sample_data(euk_data)$initial_richness <- sapply(row.names(sample_data(euk_data)), function(x) { #Taxa per sample (all ranks)
  estimate_richness(euk_data, measures = "Observed")[x,]})

#### Aggregate to genus level ####
#Aggregate taxa to genus level, after removing class taxonomic info (discrepancies lead to introduction of duplicates genera)
#Remove 'clade' column from taxonomy table
tax_table(euk_data)[,2] <- NA

#euk_genus <- tax_glom(euk_data, taxrank = "genus", NArm = FALSE, bad_empty = NA)
#saveRDS(euk_genus, file="euk_genus.RDS")
euk_genus <- readRDS("euk_genus.RDS")

print("How many genera?")
ntaxa(euk_genus)
save.image()

#Just keep non-NA classifications
euk_genus <- subset_taxa(euk_genus, !(is.na(tax_table(euk_genus)[,7])))
print("How many genera after removing NAs?")
ntaxa(euk_genus)

sample_data(euk_genus)$genus_richness <- sapply(row.names(sample_data(euk_genus)), function(x) { #Genera per sample
  estimate_richness(euk_genus, measures = "Observed")[x,]})
  
print("How many non eukaryotes?")
length(which(tax_table(euk_genus)[,1]!="Eukaryota"))
table(tax_table(euk_genus)[,1])

#Rename with genus names
taxa_names(euk_genus) <- make.names(tax_table(euk_genus)[,7], unique=TRUE)

print("How many of the non-eukaryotes were there in my prokaryote dataset?")
length(intersect(tax_table(euk_genus)[which(tax_table(euk_genus)[,1]!="Eukaryota"),7], tax_table(spe_data)[,7]))

#What about the remaining?
print("How abundant are the missing non-eukaryotes?")
extra_prok <- setdiff(tax_table(euk_genus)[which(tax_table(euk_genus)[,1]!="Eukaryota"),7], tax_table(spe_data)[,7])
summary(taxa_sums(subset_taxa(euk_genus, tax_table(euk_genus)[,7] %in% extra_prok)))

summary(taxa_sums(euk_genus))

sessionInfo()
#Stop logging
sink(file = NULL)

save.image()
