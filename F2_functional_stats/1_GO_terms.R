sink(file="log1_GO_terms.txt")

#load packages
library(readr)
library(stringr)
library(ape)
library(vegan)
library(phyloseq)
library(ggpubr)
library(tidyverse)
library(microbiome)
library(reshape2)
library(microbiome)
library(EcolUtils)

load("/proj/sllstore2017021/nobackup/MARKELLA/T3_community-level/.RData")

#Functional analysis - Script 1
#Gene ontology terms

#Import unstratified GO table
GO_unstr <- read_tsv("/proj/sllstore2017021/nobackup/MARKELLA/F1_HUMAnN2/all_GOterms_cpm_unstratified_renamed.tsv")
GO_unstr <- as.data.frame(GO_unstr)

#Rename columns
for(i in 1:length(colnames(GO_unstr))){
  j=str_remove(colnames(GO_unstr)[i], "_Abundance-RPKs")
  colnames(GO_unstr)[i]<-j
}

#Turn first column to row names
row.names(GO_unstr) <- GO_unstr$`# Gene Family`
GO_unstr <- GO_unstr[,-1]

#Remove unmapped and ungrouped
GO_unstr <- GO_unstr[-c(1,2),]

#Keep only samples that were retained as 'good' samples under community-level taxonomic analysis
GO_unstr <- GO_unstr[,which(colnames(GO_unstr) %in% sample_names(spe_data_final))]

#Remove empty rows
GO_unstr <- GO_unstr[which(rowSums(GO_unstr)>0),]

### Create phyloseq object ###
GO_phyloseq <- phyloseq(otu_table(GO_unstr, taxa_are_rows = TRUE), sample_data(spe_data_final))

#Create a subset of the phyloseq object that contains only biological processes
#Not normalizing on purpose
GO_BP_phyloseq <- phyloseq(otu_table(GO_unstr[which(grepl("BP", rownames(GO_unstr))),], taxa_are_rows = TRUE), sample_data(spe_data_final))

#save object
saveRDS(GO_BP_phyloseq, "GO_BP_phyloseq")

### PCoA ### 

#Calculate euclidean distance matrix
GO_eucl <- ordinate(GO_BP_phyloseq, method="PCoA", distance="euclidean")

#Plot and save ordinations
pdf(file = "/proj/sllstore2017021/nobackup/MARKELLA/F2_functional_stats/GO_abund_ordination.pdf")
plot_ordination(GO_BP_phyloseq, GO_eucl, color="Spec.subspecies", shape = "Seq.centre", title="PCoA plot based on normalized relative abundance\nof biological processes") + 
  theme_bw() + geom_point(size=4)  + theme(legend.text = element_text(size=15), legend.title = element_text(size=15)) +
  scale_color_manual(values=c("darkgoldenrod2", "dark grey", "steelblue1"))
dev.off()

### PERMANOVA ###

GO_BP_model <- adonis(t(otu_table(GO_BP_phyloseq)) ~ 
                     sample_data(GO_BP_phyloseq)$readcount.m.before.Kraken + 
                     as.factor(sample_data(GO_BP_phyloseq)$Seq.centre) +
                     as.factor(sample_data(GO_BP_phyloseq)$Spec.subspecies), permutations = 10000, method = "euclidean")
GO_BP_model

#### Pairwise PERMANOVA ####

#using a Jaccard dissimilarity matrix
vegdist(t(otu_table(GO_BP_phyloseq)), method="euclidean") %>%
  adonis.pair(Factor=sample_data(GO_BP_phyloseq)$Spec.subspecies)

sink(file=NULL)
save.image()
