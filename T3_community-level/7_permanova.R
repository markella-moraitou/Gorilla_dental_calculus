#Start logging
sink(file = "log7_permanova.txt")

#load packages
library(tidyverse)
library(vegan)
library(readxl)
library(phyloseq)
library(EcolUtils)
library(microbiome)

load(".RData")

#Community-level taxonomic analysis - Script 7
#PERMANOVA after decontaminating and filtering for abundance

#Add sexing information
meta <- read_tsv("tables_and_stats/S1_Sample_metadata.tsv", skip=1)
gorilla_sexing <- meta %>% select(Sample_ID,contains("Sex")) %>% 
  mutate(`Sex based on museum records`=str_split(`Sex based on museum records`," ",simplify = T)[,1],
         `Sex based on museum records`=ifelse(`Sex based on museum records`=="Unknown",NA,`Sex based on museum records`),
         `Sex assigned molecularly (using sexassign)`=ifelse(`Sex assigned molecularly (using sexassign)`=="Unknown",NA,`Sex assigned molecularly (using sexassign)`))

gorilla_age <- meta %>% select(Sample_ID,Age) %>% mutate(Age=ifelse(Age=="Unknown",NA,Age),
                                                         Age=fct_recode(Age,"Juvenile"="Subadult","Adult"="Adult (Old)"))

#Add sexing info to phyloseq metadata
spe_data_final@sam_data$museum_sex <- gorilla_sexing$`Sex based on museum records`[match(sample_names(spe_data_final), gorilla_sexing$Sample_ID)]
spe_data_final@sam_data$molecular_sex <- gorilla_sexing$`Sex assigned molecularly (using sexassign)`[match(sample_names(spe_data_final), gorilla_sexing$Sample_ID)]
spe_data_final_norm@sam_data$museum_sex <- gorilla_sexing$`Sex based on museum records`[match(sample_names(spe_data_final_norm), gorilla_sexing$Sample_ID)]
spe_data_final_norm@sam_data$molecular_sex <- gorilla_sexing$`Sex assigned molecularly (using sexassign)`[match(sample_names(spe_data_final_norm), gorilla_sexing$Sample_ID)]

#Add age info to phyloseq metadata
spe_data_final@sam_data$age <- gorilla_age$Age[match(sample_names(spe_data_final), gorilla_age$Sample_ID)]
spe_data_final_norm@sam_data$age <- gorilla_age$Age[match(sample_names(spe_data_final_norm), gorilla_age$Sample_ID)]

print("Sex-subspecies table")
table(spe_data_final@sam_data$Spec.subspecies, spe_data_final@sam_data$museum_sex)
table(spe_data_final@sam_data$Spec.subspecies, spe_data_final@sam_data$molecular_sex)

print("age-subspecies table")
table(spe_data_final@sam_data$Spec.subspecies, spe_data_final@sam_data$age)

read_count <- as.numeric(sample_data(spe_data_final)$readcount.m.before.Kraken)
spec.subspecies <- as.factor(sample_data(spe_data_final)$Spec.subspecies)
seq_centre <- as.factor(sample_data(spe_data_final)$Seq.centre)

####Model 1 - 1. Read count, 2. Spec.subspecies, 3.Seq. centre ####
print("1. Read count, 2. Spec.subspecies, 3. Seq.centre")
#using a Jaccard dissimilarity matrix
model1_jaccard <- adonis2(t(abundances(spe_data_final)) ~ read_count + spec.subspecies + seq_centre, permutations = 10000, method = "jaccard")
model1_jaccard

#using an Aitchison disimilarity matrix
model1_clr <- adonis2(t(abundances(spe_data_final_norm)) ~ read_count + spec.subspecies + seq_centre, permutations = 10000, method = "euclidean")
model1_clr

####Model 2 - 1. Read count, 2. Seq.centre, 3. Spec.subspecies ####
print("1. Read count, 2. Seq.centre, 3. Spec.subspecies")
#using a Jaccard dissimilarity matrix
model2_jaccard <- adonis2(t(abundances(spe_data_final)) ~ read_count + seq_centre + spec.subspecies, permutations = 10000, method = "jaccard")
model2_jaccard

#using a aitchison disimilarity matrix
model2_clr <- adonis2(t(abundances(spe_data_final_norm)) ~ read_count + seq_centre  + spec.subspecies, permutations = 10000, method = "euclidean")
model2_clr

#### Pairwise PERMANOVA ####
print("Post hoc analysis: Pairwise PERMANOVA")

print("using a Jaccard dissimilarity matrix")
vegdist(t(abundances(spe_data_final)), method="jaccard") %>%
  adonis.pair(Factor=spec.subspecies)

print("using an Aitchison dissimilarity matrix")
vegdist(t(abundances(spe_data_final_norm)), method="euclidean") %>%
  adonis.pair(Factor=spec.subspecies)

####Model 3 - 1. Read count, 2. Seq.centre, 3. sex, 4. Spec.subspecies ####
# museum sex assignments
read_count <- subset_samples(spe_data_final, !is.na(museum_sex))@sam_data$readcount.m.before.Kraken %>% as.numeric
seq_centre <- subset_samples(spe_data_final, !is.na(museum_sex))@sam_data$Seq.centre %>% as.factor
spec.subspecies <- subset_samples(spe_data_final, !is.na(museum_sex))@sam_data$Spec.subspecies %>% as.factor
sex <- subset_samples(spe_data_final, !is.na(museum_sex))@sam_data$museum_sex %>% as.factor

print("1. Read count, 2. Seq.centre, 3. museum_sex, 4. Spec.subspecies -- with a subset containing sex assignments")

model3_jaccard <- adonis2(t(abundances(subset_samples(spe_data_final, !is.na(museum_sex)))) ~ read_count + seq_centre + sex + spec.subspecies, permutations = 10000, method = "jaccard")
model3_jaccard

model3_clr <- adonis2(t(abundances(subset_samples(spe_data_final_norm, !is.na(museum_sex)))) ~ read_count + seq_centre + sex + spec.subspecies, permutations = 10000, method = "euclidean")
model3_clr

# molecular sex assignments
read_count <- subset_samples(spe_data_final, !is.na(molecular_sex))@sam_data$readcount.m.before.Kraken %>% as.numeric
seq_centre <- subset_samples(spe_data_final, !is.na(molecular_sex))@sam_data$Seq.centre %>% as.factor
spec.subspecies <- subset_samples(spe_data_final, !is.na(molecular_sex))@sam_data$Spec.subspecies %>% as.factor
sex <- subset_samples(spe_data_final, !is.na(molecular_sex))@sam_data$molecular_sex %>% as.factor

model3_jaccard <- adonis2(t(abundances(subset_samples(spe_data_final, !is.na(molecular_sex)))) ~ read_count + seq_centre + sex + spec.subspecies, permutations = 10000, method = "jaccard")
model3_jaccard

model3_clr <- adonis2(t(abundances(subset_samples(spe_data_final_norm, !is.na(molecular_sex)))) ~ read_count + seq_centre + sex + spec.subspecies, permutations = 10000, method = "euclidean")
model3_clr


####Model 4 - 1. Read count, 2. Seq.centre, 4. age, 4. Spec.subspecies ####
read_count <- subset_samples(spe_data_final, !is.na(age))@sam_data$readcount.m.before.Kraken %>% as.numeric
seq_centre <- subset_samples(spe_data_final, !is.na(age))@sam_data$Seq.centre %>% as.factor
spec.subspecies <- subset_samples(spe_data_final, !is.na(age))@sam_data$Spec.subspecies %>% as.factor
age <- subset_samples(spe_data_final, !is.na(age))@sam_data$age %>% as.factor

print("1. Read count, 2. Seq.centre, 4. age, 4. Spec.subspecies -- with a subset containing age assignments")

model4_jaccard <- adonis2(t(abundances(subset_samples(spe_data_final, !is.na(age)))) ~ read_count + seq_centre + age + spec.subspecies, permutations = 10000, method = "jaccard")
model4_jaccard

cat(round(model4_jaccard$R2,3),sep = "\n")
cat(round(model4_jaccard$`Pr(>F)`,3),sep = "\n")

model4_clr <- adonis2(t(abundances(subset_samples(spe_data_final_norm, !is.na(age)))) ~ read_count + seq_centre + age + spec.subspecies, permutations = 10000, method = "euclidean")
model4_clr

cat(round(model4_clr$R2,3),sep = "\n")
cat(round(model4_clr$`Pr(>F)`,3),sep = "\n")

####Model 5 - 1. Read count, 2. Seq.centre, 3. sex, 4. age, 4. Spec.subspecies ####
# museum sex assignments
read_count <- subset_samples(spe_data_final, !is.na(museum_sex) & !is.na(age))@sam_data$readcount.m.before.Kraken %>% as.numeric
seq_centre <- subset_samples(spe_data_final, !is.na(museum_sex) & !is.na(age))@sam_data$Seq.centre %>% as.factor
spec.subspecies <- subset_samples(spe_data_final, !is.na(museum_sex) & !is.na(age))@sam_data$Spec.subspecies %>% as.factor
sex <- subset_samples(spe_data_final, !is.na(museum_sex) & !is.na(age))@sam_data$museum_sex %>% as.factor
age <- subset_samples(spe_data_final, !is.na(museum_sex) & !is.na(age))@sam_data$age %>% as.factor

print("1. Read count, 2. Seq.centre, 5. museum_sex, 4. Spec.subspecies -- with a subset containing sex assignments")

model5_jaccard <- adonis2(t(abundances(subset_samples(spe_data_final, !is.na(museum_sex) & !is.na(age)))) ~ read_count + seq_centre + sex + age + spec.subspecies, permutations = 10000, method = "jaccard")
model5_jaccard

cat(round(model5_jaccard$R2,3),sep = "\n")
cat(round(model5_jaccard$`Pr(>F)`,3),sep = "\n")

model5_clr <- adonis2(t(abundances(subset_samples(spe_data_final_norm, !is.na(museum_sex) & !is.na(age)))) ~ read_count + seq_centre + sex + age + spec.subspecies, permutations = 10000, method = "euclidean")
model5_clr

cat(round(model5_clr$R2,3),sep = "\n")
cat(round(model5_clr$`Pr(>F)`,3),sep = "\n")

# molecular sex assignments
read_count <- subset_samples(spe_data_final, !is.na(molecular_sex) & !is.na(age))@sam_data$readcount.m.before.Kraken %>% as.numeric
seq_centre <- subset_samples(spe_data_final, !is.na(molecular_sex) & !is.na(age))@sam_data$Seq.centre %>% as.factor
spec.subspecies <- subset_samples(spe_data_final, !is.na(molecular_sex) & !is.na(age))@sam_data$Spec.subspecies %>% as.factor
sex <- subset_samples(spe_data_final, !is.na(molecular_sex) & !is.na(age))@sam_data$molecular_sex %>% as.factor
age <- subset_samples(spe_data_final, !is.na(molecular_sex) & !is.na(age))@sam_data$age %>% as.factor

model5_jaccard <- adonis2(t(abundances(subset_samples(spe_data_final, !is.na(molecular_sex) & !is.na(age)))) ~ read_count + sex + age + spec.subspecies, permutations = 10000, method = "jaccard")
model5_jaccard

cat(round(model5_jaccard$R2,3),sep = "\n")
cat(round(model5_jaccard$`Pr(>F)`,3),sep = "\n")

model5_clr <- adonis2(t(abundances(subset_samples(spe_data_final_norm, !is.na(molecular_sex) & !is.na(age)))) ~ read_count + sex + age + spec.subspecies, permutations = 10000, method = "euclidean")
model5_clr

cat(round(model5_clr$R2,3),sep = "\n")
cat(round(model5_clr$`Pr(>F)`,3),sep = "\n")

sessionInfo()
#Stop logging
sink(file = NULL)

save.image()
