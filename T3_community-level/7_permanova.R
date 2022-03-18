#Start logging
sink(file = "log7_permanova.txt")

#load packages
library(vegan)
library(readxl)
library(phyloseq)
library(ggplot2)
library(EcolUtils)
library(dplyr)
library(tidyr)
library(stringr)

load(".RData")

#Community-level taxonomic analysis - Script 7
#PERMANOVA after decontaminating and filtering for abundance

#Add sexing output
gorilla_sexing <- read_xlsx("/crexGorilla_sexing_full_summary.xlsx", skip=1)

#Keep only columns with sex assignments then get a summary of all the sex assignments per sample (hopefully on test does not contradict the other)
gorilla_sexing <- gorilla_sexing[,1:7] %>%
  #merge all sex assigment columns 
  unite("sex", allreads...3:s5000...7, na.rm = TRUE, sep="") %>%
  #ignore unassigned and get the unique assigments per row
  mutate(sex=str_remove(sex, "U+") %>% str_remove("U"))

uniqchars <- function(x) unique(strsplit(x, "")[[1]]) 
gorilla_sexing$sex <- sapply(gorilla_sexing$sex, function(x) {uniqchars(x)[1]})

#Turn empty fields into NAs
gorilla_sexing$sex[gorilla_sexing$sex==""] <- NA

#Check if the result is the same in both reference
gorilla_sexing <- gorilla_sexing %>% group_by(sample) %>% summarise(consistent_result=all(sex=="M" | is.na(sex)) | all(sex=="F" | is.na(sex)), sex=sex) %>% unique %>% filter(!(is.na(sex)))


#Add sexing info to phyloseq metadata
spe_data_final@sam_data$Sex <- gorilla_sexing$sex[match(sample_names(spe_data_final), gorilla_sexing$sample)]
spe_data_final_norm@sam_data$Sex <- spe_data_final@sam_data$Sex

print("Sex-subspecies table")
table(spe_data_final@sam_data$Spec.subspecies, spe_data_final@sam_data$Sex

read_count <- as.numeric(sample_data(spe_data_final)$readcount.m.before.Kraken)
spec.subspecies <- as.factor(sample_data(spe_data_final)$Spec.subspecies)
seq_centre <- as.factor(sample_data(spe_data_final)$Seq.centre)
sex <- as.factor(sample_data(spe_data_final)$Sex)

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

####Model 3 - 1. Read count, 2. Seq.centre, 3. Sex, 4. Spec.subspecies ####
read_count <- subset_samples(spe_data_final, !is.na(Sex))@sam_data$readcount.m.before.Kraken %>% as.numeric
seq_centre <- subset_samples(spe_data_final, !is.na(Sex))@sam_data$Seq.centre %>% as.factor
spec.subspecies <- subset_samples(spe_data_final, !is.na(Sex))@sam_data$Spec.subspecies %>% as.factor
sex <- subset_samples(spe_data_final, !is.na(Sex))@sam_data$Sex %>% as.factor

print("1. Read count, 2. Seq.centre, 3. Sex, 4. Spec.subspecies -- with a subset containing sex assignments")
model3_jaccard <- adonis(t(otu_table(subset_samples(spe_data_final, !is.na(Sex)))) ~ read_count + seq_centre + sex + spec.subspecies, permutations = 10000, method = "jaccard")
model3_jaccard

model3_clr <- adonis(t(otu_table(subset_samples(spe_data_final_norm, !is.na(Sex)))) ~ read_count + seq_centre + sex + spec.subspecies, permutations = 10000, method = "euclidean")
model3_clr

sessionInfo()
#Stop logging
sink(file = NULL)

save.image()
