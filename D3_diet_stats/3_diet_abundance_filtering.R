#Start logging
sink(file = "log3_diet_abundance_filtering.txt")

#load packages
library(phyloseq)
library(vegan)
library(data.table)
library(compositions)
library(ggplot2)
library(ggpubr)
library(stringr)
library(funrar)
library(decontam)
library(dplyr)

load(".RData")

#Diet analysis - Script 2
#Abundance filtering and decontam

#### Decontam ####

sample_data(euk_genus)$is.neg=FALSE
sample_data(euk_genus)[which(sample_data(euk_genus)$plot.label == "blank/control"),]$is.neg=TRUE

#Create separate phyloseq objects for Jena and for Uppsala (including extraction and library blanks in Uppsala dataset)
uppsala_euk <- subset_samples(euk_genus, euk_genus@sam_data$Seq.centre=="Uppsala" & euk_genus@sam_data$Sample.type!="control") #only Uppsala samples, except for environmental controls
jena_euk <- subset_samples(euk_genus, euk_genus@sam_data$Seq.centre=="Jena") #only Jena samples

#Run decontam
#In this analysis, the null hypothesis is that taxa are NOT contaminants

#For Uppsala samples - taxa identified by either prevalence or frequency
isContaminant(t(uppsala_euk), method=c("either"), conc="pre.seq.conc.copies.ul", neg="is.neg", normalize = FALSE, threshold = uppsala_threshold) -> upp_contam_euk

#For Jena samples
isContaminant(t(jena_euk), method=c("either"), conc="pre.seq.conc.copies.ul",  neg="is.neg", normalize = FALSE, threshold = jena_threshold) -> jena_contam_euk

#How many contaminants identified?
print("Contaminants in the Uppsala dataset")
table(upp_contam_euk$contaminant)
print("Contaminants in the Jena dataset")
table(jena_contam_euk$contaminant)

print("How much do the contaminants from the different datasets overlap?")
row.names(upp_contam_euk[which(upp_contam_euk$contaminant==TRUE),]) %in% row.names(jena_contam_euk[which(jena_contam_euk$contaminant==TRUE),]) %>% table
row.names(jena_contam_euk[which(jena_contam_euk$contaminant==TRUE),]) %in% row.names(upp_contam_euk[which(upp_contam_euk$contaminant==TRUE),]) %>% table

#Do these species even exist in the other dataset?
print("Which Uppsala contaminants are found in the Jena dataset?")
table(rowSums(otu_table(jena_euk)[row.names(upp_contam_euk[which(upp_contam_euk$contaminant==TRUE),]),])>0)
print("Which Jena contaminants are found in the Uppsala dataset?")
table(rowSums(otu_table(uppsala_euk)[row.names(jena_contam_euk[which(jena_contam_euk$contaminant==TRUE),]),])>0)

#Remove identified contaminants
#combine the two contaminant sets
contam_euk <- rbind(jena_contam_euk[which(jena_contam_euk$contaminant==TRUE),], upp_contam_euk[which(upp_contam_euk$contaminant==TRUE),])

#Remove contaminants from full dataset
euk_genus_decontam <- subset_taxa(euk_genus, !(taxa_names(euk_genus) %in% row.names(contam_euk)))

#Remove blanks from dataset
euk_genus_decontam <- subset_samples(euk_genus_decontam, !(Sample.type %in% c("ExtBlank", "LibBlank")))

#remove empty rows from decontam object
euk_genus_decontam <- subset_taxa(euk_genus_decontam, rowSums(otu_table(euk_genus_decontam))>0)

print("How many of the contaminants were in the diet references list?")
intersect(rownames(contam_euk), rownames(diet_ref_gen))

#Include metadata column with species richness after decontam
sample_data(euk_genus_decontam)$richness_after_decontam <- sapply(row.names(sample_data(euk_genus_decontam)), function(x) { #Species per sample after decontam
  estimate_richness(euk_genus_decontam, measures = "Observed")[x,]})

#### Filter for abundance ####

#Remove all taxa that have under 10 reads in a given species (turn to zero)
euk_genus_filt <- euk_genus_decontam
otu_table(euk_genus_filt)[otu_table(euk_genus_filt)<=10] <- 0

#Remove taxa with empty rows
euk_genus_filt <- prune_taxa(taxa_sums(euk_genus_filt)>0, euk_genus_filt)

#Include metadata column with species richness after abundance filtering
sample_data(euk_genus_filt)$richness_after_filt <- sapply(row.names(sample_data(euk_genus_filt)), function(x) { #Genera per sample after filtering
  estimate_richness(euk_genus_filt, measures = "Observed")[x,]})

#### Remove taxa on environmental controls #### 

#taxa present in environmental controls
env_taxa_euk <- euk_genus_filt %>% subset_samples(Sample.type=="control") %>% otu_table %>% as.data.frame %>% mutate(taxa_sums=rowSums(.)) %>% filter(taxa_sums>0) %>% rownames

euk_genus_decontam2 <- subset_taxa(euk_genus_filt, !(taxa_names(euk_genus_filt) %in% env_taxa_euk))

#remove environmental controls
euk_genus_decontam2 <- subset_samples(euk_genus_decontam2, Sample.type!="control")

#Include metadata column with species richness after decontamination using museum taxa
sample_data(euk_genus_decontam2)$richness_after_envrem <- sapply(row.names(sample_data(euk_genus_decontam2)), function(x) { #Species per sample after abundance-based decontam
  estimate_richness(euk_genus_decontam2, measures = "Observed")[x,]})


print("Have a look at the dataset")
table(tax_table(euk_genus_decontam2)[,1])
table(tax_table(subset_taxa(euk_genus_decontam2, tax_table(euk_genus_decontam2)[,1]=="Eukaryota"))[,3])

#Stop logging
sink(file = NULL)

save.image() 

