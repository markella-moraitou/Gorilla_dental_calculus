#Start logging
sink(file = "log1_phyloseq.txt")

#load packages
library(readxl)
library(phyloseq)
library(taxize)
library(curl)
library(dplyr)
library(compositions)
library(stringr)
library(ggplot2)
library(microbiome)
library(ggpubr)

#Community-level taxonomic analysis - Script 1
#Set up phyloseq object for analysis

### Importing Data ###

#Import kraken-biom species tables into R
species_table <- read.table("T2_bracken/species_table.txt", skip=1, sep="\t", comment.char="", header=T)

#Modifications on otu tables
species_table <- species_table[,which(colSums(species_table)!=0)] #remove columns with sum 0
species_table <- t(species_table) #transpose table. Vegdist documentation refers to columns as species and rows as sites.
colnames(species_table)<-as.character(species_table[1,]) #set first row as column (taxa) names
species_table <- species_table[-1,] #remove original first row


#for-loop to keep only the base name for sample names
for(i in 1:length(row.names(species_table))){
  j=str_remove(row.names(species_table)[i], "_m_kraken2_report_bracken_species")
  row.names(species_table)[i]<-j
}
species_table=t(species_table) #transpose again to use phyloseq

species_table_norm=microbiome::transform(species_table, transform = "clr") #CLR normalization

#Import metadata
metadata <- read_excel("metadata_gorilla_analysis.xlsx", sheet = "Sheet1")

#Modifications on metadata table
metadata <- as.data.frame(metadata) #turn from tibble to data frame
row.names(metadata) <- metadata[,1] #Set row names
metadata <- metadata[,-1]
metadata <- Filter(function(x)!all(is.na(x)), metadata) #remove columns with all values = NA
metadata <- subset(metadata, row.names(metadata) %in% colnames(species_table)) #subset metadata to have only the samples that exist in the OTU table
metadata <- metadata[order(row.names(metadata)),] #order by row name

#Make subspecies an ordered factor
metadata$Spec.subspecies <- factor(metadata$Spec.subspecies, levels = c("gorilla", "graueri", "beringei"))

#Add read count before Kraken
##Read counts after human-host mapping
rc_human_host <- read.table("8_humanHostFilt/unmapped/readcount_hostmapping.txt", sep=",", comment.char="", header=T)

#Fix row names
rc_human_host$sample <- str_remove(rc_human_host$sample, "_m") %>% str_remove("./")

#Add columns
metadata$readcount.m.before.Kraken <- rc_human_host$unmapped.read..[match(rownames(metadata), rc_human_host$sample)]

rm(rc_human_host)

#Save metadata object
saveRDS(metadata, file="metadata.RDS")

#Create taxonomic table for species-level classifications
#Already created and it takes time, so I am just getting it from the save file
taxonomy_species_original <- read_excel("T3_community-level/taxonomy_species_original.xlsx", sheet = "taxonomy_species_original")

#taxonomy_species_original <- matrix( nrow = length(row.names(species_table_norm)), ncol = 9) #create empty taxonomic matrix
#colnames(taxonomy_species_original)=c("ID", "superkingdom", "clade", "phylum", "class", "order", "family", "genus", "species") #set names of taxonomic rankings

#for loop to search the classification of every taxonomic ID in the species table
#for (i in 1:length(row.names(species_table))){
#  taxon=row.names(species_table)[i] 
#  classification = classification(sci_id=taxon, db="ncbi") #search for classification info of every taxon ID
#  for (j in 1:length(colnames(taxonomy_species_original))){ #for every taxonomic ranking in the matrix's header
#    for (k in 1:length(classification[[taxon]]$rank)){ #find the field in the classification object that matches
#      if (as.character(classification[[taxon]]$rank[k]) == as.character(colnames(taxonomy_species_original)[j])){
#        taxonomy_species_original[i,j] = classification[[taxon]]$name[k] #use name of ranking to fill in the matrix
#        taxonomy_species_original[i,1] = classification[[taxon]]$id[k] #use taxon ID to fill in the matrix
#      }
#    }
#  }
#}

taxonomy_species = as.data.frame(taxonomy_species_original) #create a copy of the original matrix as a data frame

#Turn string "NA"s into real NAs
taxonomy_species[taxonomy_species=="NA"] <- NA

#Check if there are differences between the taxa in species_table and taxonomy_species
unwantedspecies=unique(taxonomy_species$ID[! taxonomy_species$ID %in% row.names(species_table)])

#remove the unwanted taxa
taxonomy_species=taxonomy_species[! (taxonomy_species$ID %in%  unwantedspecies),]

missingspecies=unique(row.names(species_table)[! row.names(species_table) %in% taxonomy_species$ID])

#add missing species (this might not be necessary if there aren't missing species)
taxonomy_species_extension <- matrix(nrow = length(missingspecies), ncol = 9)
colnames(taxonomy_species_extension)=c("ID", "superkingdom", "clade", "phylum", "class", "order", "family", "genus", "species") #set names of taxonomic rankings

for (i in 1:length(missingspecies)){
  taxon=missingspecies[i]
  classification = classification(sci_id=taxon, db="ncbi")
  for (j in 1:length(colnames(taxonomy_species_extension))){
    for (k in 1:length(classification[[taxon]]$rank)){
      if (as.character(classification[[taxon]]$rank[k]) == as.character(colnames(taxonomy_species_extension)[j])){
        taxonomy_species_extension[i,j] = classification[[taxon]]$name[k]
        taxonomy_species_extension[i,1] = taxon #using taxonomy ID from my OTU table, because there are sometimes different IDs that cause discrepancies between the two tables
      }
    }
  }
}

taxonomy_species_extension = as.data.frame(taxonomy_species_extension)

#merge the two taxonomy tables and do a few modifications
taxonomy_species <- rbind(taxonomy_species, taxonomy_species_extension)

#Save it as text for back-up, since it takes long to generate
write.table(taxonomy_species, "T3_community-level/taxonomy_species_original.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

taxonomy_species <- taxonomy_species[order(taxonomy_species$ID),] #order by ID
taxonomy_species <- unique(taxonomy_species) #keep only unique entries
row.names(taxonomy_species)<- taxonomy_species$ID #set row names
taxonomy_species=taxonomy_species[,-1]
taxonomy_species=as.matrix(taxonomy_species) #Apparently phyloseq likes taxonomy tables to be matrices
#Again, make sure "NA"s are NAs
taxonomy_species[taxonomy_species=="NA"] <- NA


### Create phyloseq objects ###

#create phyloseq objects from: 
#sample_data object based on metadata
#otu_table based on species_table and species_table_norm respectively
spe_data <- phyloseq(otu_table(species_table, taxa_are_rows = TRUE), sample_data(metadata), tax_table(taxonomy_species))

#Log-transform read counts
spe_data@sam_data$readcount_log <- log(spe_data@sam_data$readcount.m.before.Kraken)
#Check distribution of read
ggsave(file="read_count_raw_hist.png",
            gghistogram(spe_data@sam_data, "readcount_log", fill="Sample.type", position="stack") +
            geom_vline(xintercept= log(300000), size=1, colour="red"),
            device="png")
spe_data@sam_data$readcount_log <- NULL

#Remove samples with very few reads, according to the distribition
low_content_samples <- subset_samples(spe_data, readcount.m.before.Kraken < 300000 & Sample.type=="sample") %>% sample_names
spe_data <- subset_samples(spe_data, !(sample_names(spe_data) %in% low_content_samples))

print("Sample excluded because they were low in content")
low_content_samples

#Add metadata about species richness
sample_data(spe_data)$richness_initial <- sapply(sample_names(spe_data)[match(sample_names(spe_data), sample_names(spe_data))], function(x) { #Species richness
  estimate_richness(spe_data, measures = "Observed")[x,]})

spe_data_norm <- phyloseq(otu_table(species_table_norm, taxa_are_rows = TRUE), sample_data(metadata), tax_table(taxonomy_species))


### Plot ###

### PCoA using Jaccard distances ###

#PCoA using jaccard for species table
jaccard_1 <- ordinate(spe_data, method="PCoA", distance="jaccard")


#Plot and save ordinations
pdf(file = "T3_community-level/jaccard_1.pdf")
par(mfrow=c(1,3))
plot_ordination(spe_data, jaccard_1, color=c("Seq.centre"), shape="Spec.subspecies", title="PCoA on jaccard distances based on species level assignments") #to highlight seq. centre
plot_ordination(spe_data, jaccard_1, color=c("plot.label"), title="PCoA on jaccard distances based on species level assignments") #to highlight subspecies
plot_ordination(spe_data, jaccard_1, color=c("Seq.centre"), shape="Duplic.pair", title="PCoA on jaccard distances based on species level assignments") + #to highlight samples sequences both in Jena and Uppsala
  scale_shape_manual(values=c(0, 1, 5, 6, 16))
dev.off()

### PCoA using Aitchison distances ###

#PCoA using euclidean distances of CLR-normalized abundances on the species table 
clr_1 <- ordinate(spe_data_norm, method="PCoA", distance="euclidean")

#Plot and save ordinations
pdf(file = "T3_community-level/clr_1.pdf")
par(mfrow=c(1,3))
plot_ordination(spe_data_norm, clr_1, color=c("Seq.centre"), shape="Spec.subspecies", title="PCoA on Aitchison distances based on species level assignments") #to highlight seq. centre
plot_ordination(spe_data_norm, clr_1, color=c("plot.label"), title="PCoA on Aitchison distances based on species level assignments") #to highlight subspecies
plot_ordination(spe_data_norm, clr_1, color=c("Seq.centre"), shape="Duplic.pair", title="PCoA on Aitchison distances based on species level assignments") + #to highlight samples sequences both in Jena and Uppsala
  scale_shape_manual(values=c(0, 1, 5, 6, 16))
dev.off()

sessionInfo()

#Stop logging
sink(file = NULL)

save.image() 
