#Start logging
sink(file = "log5_env_taxa_removal.txt")

#load packages
library(phyloseq)
library(funrar)
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(reshape2)
library(dplyr)
library(stringr)

load(".RData")

#Community-level taxonomic analysis - Script 5
#Decontamination method based on the relative abundance 
#of taxa insamples and environmental controls

#Store previously discovered env taxa in a separate object, to compare later
#env_taxa_old <- env_taxa

#### Calculate relative abundance ratios ####

#Define function for calculating relative abundance ratios
rel.abund.ratios <- function(taxa_table_ra, sample_type_vector) {
  #Calculate rel. abundance ratios between each sample and one environmental control
  
  samplenames <- colnames(taxa_table_ra)[which(sample_type_vector=="sample")]
  controls <- colnames(taxa_table_ra)[which(sample_type_vector=="control")]

  #Get long table with columns containing: taxon, sample, rel. abund in sample, rel. abund in control1, rel. abund in control 2 etc...
  taxa_table_ra_m <- reshape(taxa_table_ra, idvar = "taxon", ids=row.names(taxa_table_ra),
                                  times=names(taxa_table_ra)[-c(which(sample_type_vector=="control"))],
                                  timevar = "sample", varying = list(names(taxa_table_ra)[-c(which(sample_type_vector=="control"))]),
                                  direction = "long")
  k <- ncol(taxa_table_ra_m)
  colnames(taxa_table_ra_m)[k-1] <- "rel.abund"
  #Reorder of columns
  taxa_table_ra_m <- taxa_table_ra_m[,c(k,k-2,k-1,1:(k-3))]

  #Make long again
  taxa_table_ra_m <- reshape(taxa_table_ra_m, idvar="taxon.sample", ids=rownames(taxa_table_ra_m),
                                  times=names(taxa_table_ra_m)[4:ncol(taxa_table_ra_m)], timevar = "control",
                                  varying = list(names(taxa_table_ra_m)[4:ncol(taxa_table_ra_m)]),
                                  direction = "long")
  
  #Tidy up table a bit
  taxa_table_ra_m[,5] <- as.numeric(taxa_table_ra_m[,5])
  row.names(taxa_table_ra_m) <- NULL
  colnames(taxa_table_ra_m) <- c("taxon", "sample", "rel.abund", "control", "rel.abund.in.control")
  taxa_table_ra_m <- taxa_table_ra_m[,-6]


  #Remove rows where the control abundances is 0 (will produce infinite values when calculating the ratios)
  taxa_table_ra_m <- taxa_table_ra_m[which(taxa_table_ra_m$rel.abund.in.control>0),]

  #Create a matrix with the relative abundance ratio for each sample and control combination
  ra_ratio <- taxa_table_ra_m[,c("taxon","sample","control")]
  ra_ratio$ratio <- taxa_table_ra_m$rel.abund/taxa_table_ra_m$rel.abund.in.control
  return(ra_ratio)
}

#Transform abundances into relative abundances (using otu table after 0.01% filtering)
species_table_filt_ra <- t(make_relative(t(otu_table(spe_data_filt))))
species_table_filt_ra <- as.data.frame(species_table_filt_ra)

#Create vector that indicates in a sample is control or not
control.or.not <- ifelse(grepl("ERR", colnames(species_table_filt_ra)) | grepl("BS", colnames(species_table_filt_ra)),
                          "control", "sample")

#Calculate relative abundance ratios
ra_ratio <- rel.abund.ratios(species_table_filt_ra, control.or.not)

#### Scatterplots ####
#Plot separately for each sample
ra_ratio_plots <- list()

samplenames <- colnames(species_table_filt_ra)[which(control.or.not=="sample")]

for (i in 1:length(samplenames)) {
  #Get a subset per sample
  subset <- ra_ratio[ra_ratio$sample==samplenames[i],]
  #Get the abundance rank of every taxon in the sample
  order <- species_table_filt_ra[order(species_table_filt_ra[,samplenames[i]], decreasing=TRUE),]
  order <- cbind(rownames(order), 1:nrow(order))
  #Add rank info to subset
  subset$rank <- as.numeric(order[match(subset$taxon, order[,1]),2])
  #Remove taxa with zero abundance in samples
  subset <- subset[which(subset$ratio > 0),]
  #Log transform
  subset$ratio <- log(subset$ratio)
  #Plot
  plot <- ggscatter(subset, "rank", "ratio", color = "control") +
    annotate("segment", size = 0.5, y = 0, yend = 0, x=0, xend=max(subset$rank), colour = "black") +
    xlab("Abundance rank in sample") + ylab("Rel.abund ration")
  ra_ratio_plots[[i]] <- plot
}

#Plot grid
n <- length(ra_ratio_plots)
nCol <- floor(sqrt(n))

#Save plot
ggsave("relabund_ratio_per_sample.png",
  plot = grid.arrange(grobs = ra_ratio_plots, ncol = nCol),
  device = "png",
  path = "T3_community-level",
  width = 100, height = 100, units = "cm")


#### Get "environmental" taxa ####

#Find taxa that are consistently more relatively abundant in controls

#Create a table that shows if a taxon is less abundant in the sample (for every sample-control comparison)
less_in_samples <- ra_ratio[,c(1,2,3)]
less_in_samples$less.abundant.in.sample <- sapply(ra_ratio$ratio, function(x) {x<1})

#Create a table that shows if a taxon passes the criterion 
#of being less abundant in all samples than at least one control
#If at least one sample per group is FALSE (meaning: more abundant in the sample), then the whole group is FALSE (meaning: not a contaminant)
always_less <- less_in_samples %>% group_by(taxon, control) %>%
  dplyr::summarize(always.more.in.env = all(less.abundant.in.sample == "TRUE")) %>%
  #if at least one of the comparisons yield TRUE the taxon gets marked as TRUE
  group_by(taxon) %>% dplyr::summarise(environmental = any(always.more.in.env == "TRUE"))
print("How many taxa among those investigated are environmental?")
table(always_less$environmental)


#Get species names of the taxa identified as environmental
env_taxa <- always_less$taxon[which(always_less$environmental=="TRUE")]
env_taxa <- as.data.frame(env_taxa)
env_taxa$species <- taxonomy_species[match(env_taxa$env_taxa, rownames(taxonomy_species)),8]
colnames(env_taxa) <- c("TaxID", "Species")

#Are any of these oral?

#core hominid microbiome
print("How many of these are core hominid microbiome taxa?")
env_taxa[which(env_taxa$Species %in% core_micr$Taxon),]


#HOMD
print("How many of these are HOMD taxa?")
env_taxa[which(env_taxa$TaxID %in% homd$NCBI_taxon_id),]


### Boxplot per taxon ###
#Boxplot of rel. abund. ratios per taxon

#Log transform data with pseudocount
ra_ratio_log <- ra_ratio
ra_ratio_log$ratio <- log(ra_ratio$ratio + 0.0001) 
colnames(ra_ratio_log)[4] <- "log.ratio"

#Look at the relative abundance of taxa flagged as environmental
env_taxa_relabund <- ggscatter(ra_ratio_log[which(ra_ratio_log$taxon %in% env_taxa$TaxID),],
                               x = "taxon", y = "log.ratio", color = "control",
                               bxp.errorbar = FALSE) +
  annotate("segment", size = 1, y = 0, yend = 0, x=0, xend=nrow(env_taxa), colour = "black")

env_taxa_relabund

#### Dataset after removal of environmental taxa ####

#remove environmental controls
spe_data_envrem <- subset_samples(spe_data_filt, !(grepl("ERR", sample_names(spe_data_filt)) | grepl("BS", sample_names(spe_data_filt))))

#keep only present taxa
spe_data_envrem <- subset_taxa(spe_data_envrem, rowSums(otu_table(spe_data_envrem))>0)

#Remove environmental taxa
spe_data_envrem <- subset_taxa(spe_data_envrem, !(taxa_names(spe_data_envrem) %in% env_taxa$TaxID))


#Include metadata column with species richness after abundance ratio-based decontam
sample_data(spe_data_envrem)$richness_after_envrem <- sapply(row.names(sample_data(spe_data_envrem)), function(x) { #Species per sample before filtering
  estimate_richness(spe_data_envrem, measures = "Observed")[x,]})

spe_data_envrem_norm <- phyloseq(otu_table(microbiome::transform(otu_table(spe_data_envrem), transform = "clr"), taxa_are_rows = TRUE), tax_table(spe_data_envrem), sample_data(spe_data_envrem))

#How many are oral?
#core hominid microbiome
print("How many of the retained taxa are from the core hominid microbiome?")
length(taxa_names(spe_data_envrem)[which(tax_table(spe_data_envrem)[,8] %in% core_micr$Taxon)])


#HOMD
print("How many of the retained taxa are from HOMD?")
length(taxa_names(spe_data_envrem)[which(taxa_names(spe_data_envrem) %in% homd$NCBI_taxon_id)])



#### Plot ####

### PCoA using Jaccard distances ###

#PCoA using jaccard for species table
jaccard_4 <- ordinate(spe_data_envrem, method="PCoA", distance="jaccard")

#Plot ordination
pdf(file = "T3_community-level/jaccard_4.pdf")
par(mfrow=c(1,3))
plot_ordination(spe_data_envrem, jaccard_4, color=c("Seq.centre"), shape="Spec.subspecies", title="PCoA on jaccard distances based on species level assignments \nafter abudance filtering and removing \"bad\" samples and environmental taxa") #to highlight seq. centre
plot_ordination(spe_data_envrem, jaccard_4, color=c("Spec.subspecies"), title="PCoA on jaccard distances based on species level assignments \nafter abudance filtering and removing \"bad\" samples and environmental taxa") #to highlight subspecies
plot_ordination(spe_data_envrem, jaccard_4, color=c("Seq.centre"), shape="Duplic.pair", title="PCoA on jaccard distances based on species level assignments \nafter abudance filtering and removing \"bad\" samples and environmental taxa") + #to highlight samples sequences both in Jena and Uppsala
  scale_shape_manual(values=c(0, 1, 5, 6, 16))
dev.off()

### PCoA using Aitchison distances ###

#PCoA using euclidean distances of CLR-normalized abundances on the species table 
clr_4 <- ordinate(spe_data_envrem_norm, method="PCoA", distance="euclidean")

#Plot ordination
pdf(file = "T3_community-level/clr_4.pdf")
par(mfrow=c(1,3))
plot_ordination(spe_data_envrem_norm, clr_4, color=c("Seq.centre"), shape="Spec.subspecies", title="PCoA on Aitchison distances based on species level assignments \nafter abudance filtering and removing \"bad\" samples and environmental taxa") #to highlight seq. centre
plot_ordination(spe_data_envrem_norm, clr_4, color=c("Spec.subspecies"), title="PCoA on Aitchison distances based on species level assignments \nafter abudance filtering and removing \"bad\" samples and environmental taxa") #to highlight subspecies
plot_ordination(spe_data_envrem_norm, clr_4, color=c("Seq.centre"), shape="Duplic.pair", title="PCoA on Aitchison distances based on species level assignments \nafter abudance filtering and removing \"bad\" samples and environmental taxa") + #to highlight samples sequences both in Jena and Uppsala
  scale_shape_manual(values=c(0, 1, 5, 6, 16))
dev.off()

sessionInfo()
#Stop logging
sink(file = NULL)
save.image()
