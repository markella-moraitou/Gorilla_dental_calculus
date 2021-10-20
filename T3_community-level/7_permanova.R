#Start logging
sink(file = "log7_permanova.txt")

#load packages
library(vegan)
library(readxl)
library(phyloseq)
library(ggplot2)
library(EcolUtils)
library(dplyr)

load(".RData")

#Community-level taxonomic analysis - Script 7
#PERMANOVA after decontaminating and filtering for abundance

read_count <- sample_data(spe_data_final)$readcount.m.before.Kraken
spec.subspecies <- as.factor(sample_data(spe_data_final)$Spec.subspecies)
seq_centre <- as.factor(sample_data(spe_data_final)$Seq.centre)


####Model 1 - 1. Read count, 2. Spec.subspecies, 3.Seq. centre ####
print("1. Read count, 2. Spec.subspecies, 3. Seq.centre")
#using a Jaccard dissimilarity matrix
model1_jaccard <- adonis(t(otu_table(spe_data_final)) ~ read_count + spec.subspecies + seq_centre, permutations = 10000, method = "jaccard")

model1_jaccard

#using an Aitchison disimilarity matrix
model1_clr <- adonis(t(otu_table(spe_data_final_norm)) ~ read_count + spec.subspecies + seq_centre, permutations = 10000, method = "euclidean")

model1_clr

####Model 2 - 1. Read count, 2. Seq.centre, 3. Spec.subspecies ####
print("1. Read count, 2. Seq.centre, 3. Spec.subspecies")
#using a Jaccard dissimilarity matrix
model2_jaccard <- adonis(t(otu_table(spe_data_final)) ~ read_count + seq_centre + spec.subspecies, permutations = 10000, method = "jaccard")

model2_jaccard

#using a aitchison disimilarity matrix
model2_clr <- adonis(t(otu_table(spe_data_final_norm)) ~ read_count + seq_centre  + spec.subspecies, permutations = 10000, method = "euclidean")

model2_clr

#### Pairwise PERMANOVA ####
print("Post hoc analysis: Pairwise PERMANOVA")

print("using a Jaccard dissimilarity matrix")
vegdist(t(otu_table(spe_data_final)), method="jaccard") %>%
  adonis.pair(Factor=spec.subspecies)

print("using an Aitchison dissimilarity matrix")
vegdist(t(otu_table(spe_data_final_norm)), method="euclidean") %>%
  adonis.pair(Factor=spec.subspecies)

sessionInfo()
#Stop logging
sink(file = NULL)

save.image()
