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

load("/proj/sllstore2017021/nobackup/MARKELLA/T3_community-level/.RData")

#Diet analysis - Script 1
#Set up phyloseq object for analysis

#Import kraken-biom otu table for eukaryotic taxa into R
#Based on Kraken2 output instead of Bracken
euk_table <- read.table("/proj/sllstore2017021/nobackup/MARKELLA/D2_kraken2_full_db/otu_table_kraken_fulldb.txt", skip=1, sep="\t", comment.char="", header=T)

#Modifications on otu table
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
}
#"1940210" missing -- it was deleted on March 1, 2021
#"242820" missing too

#saveRDS(taxonomy_euk_original, file = "/proj/sllstore2017021/nobackup/MARKELLA/D2_kraken2_full_db/taxonomy_euk_original_kraken.rds")
taxonomy_euk_original <- readRDS(file = "/proj/sllstore2017021/nobackup/MARKELLA/D3_diet_stats/taxonomy_euk_original_kraken.rds")

taxonomy_euk <- matrix(nrow = nrow(euk_table), ncol = 9)
colnames(taxonomy_euk) <- c("taxonID","superkingdom", "clade", "phylum", "class", "order", "family", "genus", "species")
#for loop to get a table format
for (i in c(1:59309, 59311:74448, 74450:nrow(euk_table)) ){
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


#The two missing taxa had to be added manually
taxonomy_euk[which(rownames(euk_table)=="1940210"),] <- c("1940210", NA, NA, NA, NA, NA, NA, "Anthracinomyces", NA)
taxonomy_euk[which(rownames(euk_table)=="242820"),] <- c("242820", "Eukaryota", NA, "Chordata", "Actinopterygii", "Cichliformes", "Cichlidae", "Andinoacara", "Andinoacara rivulatus")

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

#Remove mislabelled blanks
euk_data <- subset_samples(euk_data, !(sample_names(euk_data) %in% c("BE113", "BL109")))
euk_data <- subset_taxa(euk_data, taxa_sums(euk_data)>0)

#Get taxa richness data
sample_data(euk_data)$initial_richness <- sapply(row.names(sample_data(euk_data)), function(x) { #Taxa per sample (all ranks)
  estimate_richness(euk_data, measures = "Observed")[x,]})

#### Aggregate to genus level ####
#Aggregate taxa to genus level, after removing class taxonomic info (discrepancies lead to introduction of duplicates genera)
#Remove 'clade' column from taxonomy table
tax_table(euk_genus)[,2] <- NA

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
table(tax_table(euk_genus)[,1])

taxa_names(euk_genus) <- make.names(tax_table(euk_genus)[,7], unique=TRUE)

print("Which of the non-eukaryotes were the in my prokaryote dataset?")
length(intersect(tax_table(euk_genus)[which(tax_table(euk_genus)[,1]!="Eukaryota"),7], tax_table(spe_data)[,7]))
length(unique(tax_table(euk_genus)[which(tax_table(euk_genus)[,1]!="Eukaryota"),7]))
#2439 out of the 3330 non-eukaryotic genera were in my prokaryote dataset

#What about the remaining 891?
print("How abundant are the missing non-eukaryotes?")
extra_prok <- setdiff(tax_table(euk_genus)[which(tax_table(euk_genus)[,1]!="Eukaryota"),7], tax_table(spe_data)[,7])
summary(taxa_sums(subset_taxa(euk_genus, tax_table(euk_genus)[,7] %in% extra_prok)))

summary(taxa_sums(euk_genus))

print("Which are the most abundant missing non-eukaryotes?")
sort(taxa_sums(subset_taxa(euk_genus, tax_table(euk_genus)[,7] %in% extra_prok)))
taxonomy_euk[c("651822", "1461583"),]
#The two most abundant missing taxa are in the genus Fretibacterium

#Stop logging
sink(file = NULL)

save.image()