#Start logging
sink(file = "log3_decontam.txt")

#load packages
library(decontam)
library(phyloseq)
library(compositions)
library(ggplot2)
library(dplyr)
library(stringr)
library(tibble)

load(".RData")

#Community-level taxonomic analysis - Script 3
#Decontamination using extraction\library blanks

sample_data(spe_data)$is.neg=FALSE
sample_data(spe_data)[which(sample_data(spe_data)$plot.label == "blank/control"),]$is.neg=TRUE

#Create separate phyloseq objects for Jena and for Uppsala (including extraction and library blanks in Uppsala dataset)
uppsala <- subset_samples(spe_data, spe_data@sam_data$Seq.centre=="Uppsala" & spe_data@sam_data$Sample.type!="control") #only Uppsala samples, except for environmental controls
jena <- subset_samples(spe_data, spe_data@sam_data$Seq.centre=="Jena") #only Jena samples

#### Test Decontam with different threshold ####

#Load contaminant lists

salter_contam <- read.table("T3_community-level/contaminant_list_Salter.txt") %>%
  pull(V1)

#From Weyrich et al I am retaining only genus level
weyrich_contam <- read.table("T3_community-level/contaminant_list_Weyrich.txt",
                             sep="\t", header = TRUE) %>%
  pull(Genera.taxonomy) %>% str_remove(pattern=".*g__") %>% unique

#Merge into one
contaminant_list <- append(salter_contam, weyrich_contam) %>% unique

#Also append Streptomyces and Mycobacterium (they are never part of the oral microbiome)
contaminant_list <- append(contaminant_list, c("Mycobacterium", "Streptomyces"))

#Load data downloaded from https://github.com/jfy133/Hominid_Calculus_Microbiome_Evolution/blob/master/06-additional_data_files/Data_R22_S22A_CoreMicrobiome/coremicrobiome_presenceabsence_upsettable_allsoftware_maltdbnt_0.04_fracinds0.5_fracpops0.66_singleindpopsdroppedF_hostgenus_species_20190902.tsv
core_micr <- read.table("T3_community-level/core_microbiome.txt", sep = "\t", header = TRUE)
#Substitute underscores with spaces in Taxon columns
core_micr$Taxon <- sub("_", " ", core_micr$Taxon)
core_micr$Taxon <- sub("_", " ", core_micr$Taxon)
core_micr$Taxon <- sub("_", " ", core_micr$Taxon) 
core_micr$Taxon <- sub("_", " ", core_micr$Taxon) #Couldn't find better way to remove underscores

#Load data from HOMD
homd <- read.table("T3_community-level/homd_taxonomy_table.txt", sep = "\t", header = TRUE, fill=T)

#Keep only explicitly oral taxa
homd <- homd %>% filter(Body_site=="Oral")

#Choose thresholds for decontam by taking into account the proportion of the dataset that consists of oral and known contaminants

#Table is already produced, just read from file
decontam_test <- read.table(file="/crexT3_community-level/decontam_test.txt", sep=",",
                            header=TRUE)

decontam_test <- data.frame()

#Run decontam with a random threshold
id.cont.uppsala <- isContaminant(t(uppsala), method=c("combined"),
                       conc="pre.seq.conc.copies.ul", neg="is.neg", normalize = TRUE) %>% rownames_to_column
id.cont.jena <- isContaminant(t(jena), method=c("combined"),
                       conc="pre.seq.conc.copies.ul", neg="is.neg", normalize = TRUE) %>% rownames_to_column      
#Test for each threshold combination, to identify the best threshold per dataset             
for (u_threshold in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)) {
  for (j_threshold in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)) {
    #Identify contaminants in each threshold combination
    id.cont <- rbind(filter(id.cont.uppsala, p < u_threshold),
               filter(id.cont.jena, p < j_threshold)) %>%
      #Get taxonomy info of unique entries
      select(rowname) %>%
      mutate(tax_genus=taxonomy_species[match(rowname, rownames(taxonomy_species)),7],
             tax_name=taxonomy_species[match(rowname, rownames(taxonomy_species)),8]) %>% unique %>% 
             remove_rownames %>% column_to_rownames(var = "rowname")    
      #Proportion of the dataset explained by oral taxa/known contaminants etc in this subset.
      stats <- data.frame(n_contam_total=nrow(id.cont),
                n_contam_list=which(id.cont$tax_genus %in% contaminant_list) %>% length,
                n_homd=which(row.names(id.cont) %in% homd$NCBI_taxon_id) %>% length,
                n_core_micr=which(id.cont$tax_name %in% core_micr$Taxon) %>% length,
                n_streptomyces=which(id.cont$tax_genus=="Streptomyces") %>% length,
                n_mycobacterium=which(id.cont$tax_genus=="Mycobacterium") %>% length) %>%
      #add proportion
      mutate(prop_core_micr=n_core_micr/n_contam_total,
             prop_homd = n_homd/n_contam_total) %>%
      #add threshold info
      mutate(uppsala_threshold=u_threshold, jena_threshold=j_threshold)
    decontam_test <- rbind(decontam_test, stats)
  }  
}

#Save table
write.table(decontam_test, file="/crexT3_community-level/decontam_test.txt", sep=",",
            quote=FALSE, row.names=FALSE)

#Plot some heatmaps
library(microbiome)
library(ggpubr)

#Oral-contaminant ratio for different settings
#Proportion of oral taxa
decontam_prop_core_micr <- heat(decontam_test, Xvar="uppsala_threshold", Yvar="jena_threshold", fill="prop_core_micr",
     order.rows=FALSE, order.cols=FALSE) +
  scale_fill_gradient2(guide = guide_colorbar(barheight = 10, title="core microbiome\ntaxa proportion")) +
  ylab("jena_threshold") + xlab("uppsala_threshold")
  
decontam_prop_homd <- heat(decontam_test, Xvar="uppsala_threshold", Yvar="jena_threshold", fill="prop_homd",
     order.rows=FALSE, order.cols=FALSE) +
  scale_fill_gradient2(guide = guide_colorbar(barheight = 10, title="HOMD\ntaxa proportion")) +
  ylab("jena_threshold") + xlab("uppsala_threshold")

#Number of contaminants for different settings
decontam_test_contam <- heat(decontam_test, Xvar="uppsala_threshold", Yvar="jena_threshold", fill="n_contam_total",
     order.rows=FALSE, order.cols=FALSE) +
  scale_fill_gradient2(guide = guide_colorbar(barheight = 10, title="Number of\nidentified contaminants")) +
  ylab("jena_threshold") + xlab("uppsala_threshold")

#Number of taxa from contaminant list for different settings
decontam_test_contam_list <- heat(decontam_test, Xvar="uppsala_threshold", Yvar="jena_threshold", fill="n_contam_list",
     order.rows=FALSE, order.cols=FALSE) +
  scale_fill_gradient2(guide = guide_colorbar(barheight = 10, title="Number of\nknown contaminants")) +
  ylab("jena_threshold") + xlab("uppsala_threshold")

#Number of HOMD taxa for different settings
decontam_test_homd <- heat(decontam_test, Xvar="uppsala_threshold", Yvar="jena_threshold", fill="n_homd",
     order.rows=FALSE, order.cols=FALSE) +
  scale_fill_gradient2(guide = guide_colorbar(barheight = 10, title="Number of\nHOMD taxa")) +
  ylab("jena_threshold") + xlab("uppsala_threshold")

#Number of core microbiome taxa for different settings
decontam_test_core <- heat(decontam_test, Xvar="uppsala_threshold", Yvar="jena_threshold", fill="n_core_micr",
     order.rows=FALSE, order.cols=FALSE) +
  scale_fill_gradient2(guide = guide_colorbar(barheight = 10, title="Number of\ncore micr taxa")) +
  ylab("jena_threshold") + xlab("uppsala_threshold")

#Number of Streptomyces and Mycobacterium for different setting
decontam_test_strept <- heat(decontam_test, Xvar="uppsala_threshold", Yvar="jena_threshold", fill="n_streptomyces",
     order.rows=FALSE, order.cols=FALSE) +
  scale_fill_gradient2(guide = guide_colorbar(barheight = 10, title="Number of\nStreptomyces")) +
  ylab("jena_threshold") + xlab("uppsala_threshold")

decontam_test_mycob <- heat(decontam_test, Xvar="uppsala_threshold", Yvar="jena_threshold", fill="n_mycobacterium",
     order.rows=FALSE, order.cols=FALSE) +
  scale_fill_gradient2(guide = guide_colorbar(barheight = 10, title="Number of\nMycobacterium taxa")) +
  ylab("jena_threshold") + xlab("uppsala_threshold")

pdf(file="T3_community-level/decontam_test_heatmaps.pdf")
decontam_prop_core_micr
decontam_prop_homd
decontam_test_contam
decontam_test_contam_list
decontam_test_homd
decontam_test_core
decontam_test_strept
decontam_test_mycob
dev.off()

#I will use the max threshold for which we don't start removing oral taxa (the smallest oral to contam ratio)
uppsala_threshold <- decontam_test %>% filter(prop_core_micr==min(prop_core_micr)) %>% pull(uppsala_threshold) %>% unique %>% max
print("Uppsala threshold:")
uppsala_threshold

jena_threshold <- decontam_test %>% filter(prop_core_micr==min(prop_core_micr)) %>% pull(jena_threshold) %>% unique %>% max
print("Jena threshold:")
jena_threshold 

#### Run decontam ####
#In this analysis, the null hypothesis is that taxa are NOT contaminants

#For Uppsala samples - taxa identified by either prevalence or frequency
isContaminant(t(uppsala), method=c("combined"), conc="pre.seq.conc.copies.ul",
              neg="is.neg", normalize = TRUE, threshold = uppsala_threshold) -> upp_contam

#For Jena samples
isContaminant(t(jena), method=c("combined"),  conc="pre.seq.conc.copies.ul",
              neg="is.neg", normalize = TRUE, threshold = jena_threshold) -> jena_contam

#Look at taxonomic details of the identified contaminants
#For Uppsala
print("Uppsala contaminants:")
print("How many:")
taxonomy_species[which(row.names(taxonomy_species) %in% row.names(upp_contam[which(upp_contam$contaminant==TRUE),])),] %>% nrow
print("Details:")
taxonomy_species[which(row.names(taxonomy_species) %in% row.names(upp_contam[which(upp_contam$contaminant==TRUE),])),]

#For Jena
print("Jena contaminants:")
print("How many:")
taxonomy_species[which(row.names(taxonomy_species) %in% row.names(jena_contam[which(jena_contam$contaminant==TRUE),])),] %>% nrow
print("Details:")
taxonomy_species[which(row.names(taxonomy_species) %in% row.names(jena_contam[which(jena_contam$contaminant==TRUE),])),]

#How much do the contaminants from the different datasets overlap?
print("How many of the uppsala contaminants are also found in the jena contaminants:")
row.names(upp_contam[which(upp_contam$contaminant==TRUE),]) %in% row.names(jena_contam[which(jena_contam$contaminant==TRUE),]) %>%
  table
print("How many of the jena contaminants are also found in the uppsala contaminants:")
row.names(jena_contam[which(jena_contam$contaminant==TRUE),]) %in% row.names(upp_contam[which(upp_contam$contaminant==TRUE),]) %>%
  table
  
#Do these species even exist in the other dataset?
print("Which uppsala contaminants exist in the jena dataset?")
(rowSums(otu_table(jena)[row.names(upp_contam[which(upp_contam$contaminant==TRUE),]),])>0) %>% table
print("Which jena contaminants exist in the uppsala dataset?")
(rowSums(otu_table(uppsala)[row.names(jena_contam[which(jena_contam$contaminant==TRUE),]),])>0) %>% table

#Remove identified contaminants
#combine the two contaminant sets
contam <- rbind(jena_contam[which(jena_contam$contaminant==TRUE),], upp_contam[which(upp_contam$contaminant==TRUE),])

#Remove contaminants from full dataset
spe_data_decontam <- subset_taxa(spe_data, !(taxa_names(spe_data) %in% row.names(contam)))

#Remove blanks from dataset
spe_data_decontam <- subset_samples(spe_data_decontam, !(Sample.type %in% c("ExtBlank", "LibBlank")))

#Remove bad samples from dataset
spe_data_decontam <- subset_samples(spe_data_decontam, !(sample_names(spe_data_decontam) %in% bad_samples))

#remove empty rows from decontam object
spe_data_decontam <- subset_taxa(spe_data_decontam, rowSums(otu_table(spe_data_decontam))>0)

#Normalize decontaminated data
spe_data_decontam_norm <- phyloseq(otu_table(microbiome::transform(otu_table(spe_data_decontam), transform = "clr"), taxa_are_rows = TRUE), tax_table(spe_data_decontam), sample_data(spe_data_decontam))

 
#### Plot ####

### PCoA using Jaccard distances ###

#PCoA using jaccard for species table
jaccard_2 <- ordinate(spe_data_decontam, method="PCoA", distance="jaccard")

#Plot ordination
pdf(file = "T3_community-level/jaccard_2.pdf")
par(mfrow=c(1,3))
plot_ordination(spe_data_decontam, jaccard_2, color=c("Seq.centre"), shape="Spec.subspecies", title="PCoA on jaccard distances based on species level assignments \nafter decontam using lab blanks") #to highlight seq. centre
plot_ordination(spe_data_decontam, jaccard_2, color=c("Spec.subspecies"), title="PCoA on jaccard distances based on species level assignments \nafter decontam using lab blanks") #to highlight subspecies
plot_ordination(spe_data_decontam, jaccard_2, color=c("Seq.centre"), shape="Duplic.pair", title="PCoA on jaccard distances based on species level assignments \nafter decontam using lab blanks") + #to highlight samples sequences both in Jena and Uppsala
  scale_shape_manual(values=c(0, 1, 5, 6, 16))
dev.off()

### PCoA using Aitchison distances ###

#PCoA using euclidean distances of CLR-normalized abundances on the species table 
clr_2 <- ordinate(spe_data_decontam_norm, method="PCoA", distance="euclidean")

#Plot ordination
pdf(file = "T3_community-level/clr_2.pdf")
par(mfrow=c(1,3))
plot_ordination(spe_data_decontam_norm, clr_2, color=c("Seq.centre"), shape="Spec.subspecies", title="PCoA on Aitchison distances based on species level assignments \nafter decontam") #to highlight seq. centre
plot_ordination(spe_data_decontam_norm, clr_2, color=c("Spec.subspecies"), title="PCoA on Aitchison distances based on species level assignments \nafter decontam") #to highlight subspecies
plot_ordination(spe_data_decontam_norm, clr_2, color=c("Seq.centre"), shape="Duplic.pair", title="PCoA on Aitchison distances based on species level assignments \nafter decontam") + #to highlight samples sequences both in Jena and Uppsala
  scale_shape_manual(values=c(0, 1, 5, 6, 16))
dev.off()

sessionInfo()
#Stop logging
sink(file = NULL)

save.image()
