#Start logging
sink(file = "log2_diet_abundance_filtering.txt")

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

load("./RData")

#Diet analysis - Script 2
#Abundance filtering and decontam

#I am doing this before decontamination with decontam. I am not going to be looking at community-level diversity
#Instead I will try to extract dietary elements.
#Filtering is needed to remove spurious assignments

#### Decontam ####

#Create separate phyloseq objects for Jena and for Uppsala (including extraction and library blanks in Uppsala dataset)
uppsala_euk <- subset_samples(euk_genus, grepl("^[GB][0bEL]", sample_names(euk_genus))) #only Uppsala samples
jena_euk <- subset_samples(euk_genus, grepl("[A-Z][A-Q,S-Z][A-Z]", sample_names(euk_genus))) #only Jena samples

#Add information about negative controls in Uppsala metadata
sample_data(uppsala_euk)$is.neg=FALSE
sample_data(uppsala_euk)[which(sample_data(uppsala_euk)$plot.label=="Blank"),]$is.neg=TRUE

#Run decontam
#In this analysis, the null hypothesis is that taxa are NOT contaminants

#For Uppsala samples - taxa identified by either prevalence or frequency
isContaminant(t(uppsala_euk), method=c("either"), conc="pre.seq.conc.copies.ul", neg="is.neg", normalize = FALSE, threshold = c(freq_threshold, prev_threshold)) -> upp_contam_euk

#For Jena samples
isContaminant(t(jena_euk), method=c("frequency"), conc="pre.seq.conc.copies.ul", normalize = FALSE, threshold = freq_threshold) -> jena_contam_euk

#How many contaminants identified?
print("Contaminants in the Uppsala dataset")
table(upp_contam_euk$contaminant) #5366 in the Uppsala dataset
print("Contaminants in the Jena dataset")
table(jena_contam_euk$contaminant) #413

print("How much do the contaminants from the different datasets overlap?")
row.names(upp_contam_euk[which(upp_contam_euk$contaminant==TRUE),]) %in% row.names(jena_contam_euk[which(jena_contam_euk$contaminant==TRUE),]) %>% table
row.names(jena_contam_euk[which(jena_contam_euk$contaminant==TRUE),]) %in% row.names(upp_contam_euk[which(upp_contam_euk$contaminant==TRUE),]) %>% table
#130 are the same

#Do these species even exist in the other dataset?
print("Which Uppsala contaminants are found in the Jena dataset?")
table(rowSums(otu_table(jena_euk)[row.names(upp_contam_euk[which(upp_contam_euk$contaminant==TRUE),]),])>0)
#3423 out of 5366 Uppsala contaminants are found in the Jena dataset
print("Which Jena contaminants are found in the Uppsala dataset?")
table(rowSums(otu_table(uppsala_euk)[row.names(jena_contam_euk[which(jena_contam_euk$contaminant==TRUE),]),])>0)
#365 out of 413 Jena contaminants are not found in the Uppsala dataset

#Remove identified contaminants
#combine the two contaminant sets
contam_euk <- rbind(jena_contam_euk[which(jena_contam_euk$contaminant==TRUE),], upp_contam_euk[which(upp_contam_euk$contaminant==TRUE),])

#Remove contaminants from full dataset
euk_genus_decontam <- subset_taxa(euk_genus, !(taxa_names(euk_genus) %in% row.names(contam_euk)))

#Remove blanks from dataset
euk_genus_decontam <- subset_samples(euk_genus_decontam, !grepl("^[B][EL]", sample_names(euk_genus_decontam)))

#remove empty rows from decontam object
euk_genus_decontam <- subset_taxa(euk_genus_decontam, rowSums(otu_table(euk_genus_decontam))>0)

print("How many of the contaminants were in the diet references list?")
intersect(rownames(contam_euk), rownames(diet_ref_gen))
#1

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
env_taxa_euk <- euk_genus_filt %>% otu_table %>% as.matrix %>% as.data.frame %>% filter(ERR2503700>0 | ERR2868193 > 0) %>% rownames
euk_genus_decontam2 <- subset_taxa(euk_genus_filt, !(taxa_names(euk_genus_filt) %in% env_taxa_euk))

#remove environmental controls
euk_genus_decontam2 <- subset_samples(euk_genus_filt, !(grepl("ERR", sample_names(euk_genus_filt))))

#Remove environmental taxa
euk_genus_decontam2 <- subset_taxa(euk_genus_decontam2, !(taxa_names(euk_genus_decontam2) %in% env_taxa_euk))

#Include metadata column with species richness after decontamination using museum taxa
sample_data(euk_genus_decontam2)$richness_after_envrem <- sapply(row.names(sample_data(euk_genus_decontam2)), function(x) { #Species per sample after abundance-based decontam
  estimate_richness(euk_genus_decontam2, measures = "Observed")[x,]})


print("Have a look at the dataset")
table(tax_table(euk_genus_decontam2)[,1])
table(tax_table(subset_taxa(euk_genus_decontam2, tax_table(euk_genus_decontam2)[,1]=="Eukaryota"))[,3])

#Stop logging
sink(file = NULL)

save.image() 

