#Start logging
sink(file = "log4_abundance_filtering.txt")

#load packages
library(phyloseq)
library(vegan)
library(data.table)
library(compositions)
library(ggplot2)
library(ggpubr)
library(stringr)

load(".RData")

#Community-level taxonomic analysis - Script 4
#Abundance filtering

#### Abundance filtering ####

#Add metadata about species richness
sample_data(spe_data_decontam)$richness_after_decontam <- sapply(row.names(sample_data(spe_data_decontam)), function(x) { #Species per sample after decontam
  estimate_richness(spe_data_decontam, measures = "Observed")[x,]})

# Abundance filtering function based on Jaelle's Rscript_DC2_prelim_results_200911.R script
abundance_filter_f <- function(count.mat, cutoff) {
  #count.mat must be in format of samples = rows and taxa = columns
  #row.names must be sample IDs
  #cutoff must be proportion (e.g. 0.001 = 0.1% relative abundance)
  count.filt <- count.mat
  prop.mat <- as.data.frame(prop.table(as.matrix(count.filt), margin = 1))
  prop.mat$Sample <- row.names(prop.mat)
  prop.mat.m <- melt(prop.mat, by="Sample", value.name = "Proportion", variable.name = "Taxon")
  samples <- as.character(prop.mat$Sample)
  for (s in samples) {
    exclude.taxa <- as.vector(subset(prop.mat.m, Sample==s & Proportion < cutoff)$Taxon)
    count.filt[s,exclude.taxa] = 0
  }
  return(count.filt)
}


#Create abundance-filtered phyloseq objects
species_table_filt <- t(abundance_filter_f(as.data.frame(t(otu_table(spe_data_decontam))), 5e-04)) #filter OTU table of decontaminated phyloseq (the latest processing step)
species_table_filt <- species_table_filt[rowSums(species_table_filt)>0,] #removal of empty rows                                                           

#Create phyloseq object
spe_data_filt <- phyloseq(otu_table(species_table_filt, taxa_are_rows = TRUE), sample_data(spe_data_decontam), tax_table(spe_data_decontam))

#Include metadata column with species richness after abundance filtering
sample_data(spe_data_filt)$richness_after_filt <- sapply(row.names(sample_data(spe_data_filt)), function(x) { #Species per sample before filtering
  estimate_richness(spe_data_filt, measures = "Observed")[x,]})

#CLR-normalization and creation of normalized phyloseq objects
spe_data_filt_norm <- phyloseq(otu_table(microbiome::transform(species_table_filt, transform = "clr"), taxa_are_rows = TRUE), sample_data(spe_data_decontam_norm), tax_table(spe_data_decontam_norm))

#How many are oral?
#core hominid microbiome
print("How many of the retained taxa are from the core hominid microbiome?")
length(taxa_names(spe_data_filt)[which(tax_table(spe_data_filt)[,8] %in% core_micr$Taxon)])


#HOMD
print("How many of the retained taxa are from HOMD?")
length(taxa_names(spe_data_filt)[which(taxa_names(spe_data_filt) %in% homd$NCBI_taxon_id)])

# contaminants
print("How many of the retained taxa are contaminants?")
length(taxa_names(spe_data_filt)[which(taxa_names(spe_data_filt) %in% unambiguous_contam)])

#### Plot ####

### PCoA using Jaccard distances ###

#PCoA using jaccard for species table
jaccard_3 <- ordinate(spe_data_filt, method="PCoA", distance="jaccard")

#Plot ordination
pdf(file = "T3_community-level/jaccard_3.pdf")
par(mfrow=c(1,3))
plot_ordination(spe_data_filt, jaccard_3, color=c("Seq.centre"), shape="Spec.subspecies", title="PCoA on jaccard distances based on species level assignments \nafter abudance filtering") #to highlight seq. centre
plot_ordination(spe_data_filt, jaccard_3, color=c("Spec.subspecies"), title="PCoA on jaccard distances based on species level assignments \nafter abudance filtering") #to highlight subspecies
plot_ordination(spe_data_filt, jaccard_3, color=c("Seq.centre"), shape="Duplic.pair", title="PCoA on jaccard distances based on species level assignments \nafter abudance filtering") + #to highlight samples sequences both in Jena and Uppsala
  scale_shape_manual(values=c(0, 1, 5, 6, 16))
dev.off()

### PCoA using Aitchison distances ###

#PCoA using euclidean distances of CLR-normalized abundances on the species table 
clr_3 <- ordinate(spe_data_filt_norm, method="PCoA", distance="euclidean")

#Plot ordination
pdf(file = "T3_community-level/clr_3.pdf")
par(mfrow=c(1,3))
plot_ordination(spe_data_filt_norm, clr_3, color=c("Seq.centre"), shape="Spec.subspecies", title="PCoA on Aitchison distances based on species level assignments \nafter abudance filtering") #to highlight seq. centre
plot_ordination(spe_data_filt_norm, clr_3, color=c("Spec.subspecies"), title="PCoA on Aitchison distances based on species level assignments \nafter abudance filtering") #to highlight subspecies

plot_ordination(spe_data_filt_norm, clr_3, color=c("Seq.centre"), shape="Duplic.pair", title="PCoA on Aitchison distances based on species level assignments \nafter abudance filtering") + #to highlight samples sequences both in Jena and Uppsala
  scale_shape_manual(values=c(0, 1, 5, 6, 16))
dev.off()

sessionInfo()
#Stop logging
sink(file = NULL)

save.image()
