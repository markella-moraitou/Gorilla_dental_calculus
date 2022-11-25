#load packages
library(vegan)
library(readxl)
library(phyloseq)
library(ggplot2)
library(EcolUtils)
library(dplyr)
library(tidyr)
library(stringr)
library(microbiome)

euk_genus_diet<-readRDS("D3_diet_stats/euk_genus_diet.rds")
euk_genus_diet_norm<-readRDS("D3_diet_stats/euk_genus_diet_norm.rds")

#Add sexing information
meta <- read_tsv("tables_and_stats/S1_Sample_metadata.tsv", skip=1)
gorilla_sexing <- meta %>% select(Sample_ID,contains("Sex")) %>% 
  mutate(`Sex based on museum records`=str_split(`Sex based on museum records`," ",simplify = T)[,1],
         `Sex based on museum records`=ifelse(`Sex based on museum records`=="Unknown",NA,`Sex based on museum records`),
         `Sex assigned molecularly (using sexassign)`=ifelse(`Sex assigned molecularly (using sexassign)`=="Unknown",NA,`Sex assigned molecularly (using sexassign)`))

gorilla_age <- meta %>% select(Sample_ID,Age) %>% mutate(Age=ifelse(Age=="Unknown",NA,Age),
                                                         Age=fct_recode(Age,"Juvenile"="Subadult","Adult"="Adult (Old)"))

#Add sexing info to phyloseq metadata
euk_genus_diet@sam_data$museum_sex <- gorilla_sexing$`Sex based on museum records`[match(sample_names(euk_genus_diet), gorilla_sexing$Sample_ID)]
euk_genus_diet@sam_data$molecular_sex <- gorilla_sexing$`Sex assigned molecularly (using sexassign)`[match(sample_names(euk_genus_diet), gorilla_sexing$Sample_ID)]
euk_genus_diet_norm@sam_data$museum_sex <- gorilla_sexing$`Sex based on museum records`[match(sample_names(euk_genus_diet_norm), gorilla_sexing$Sample_ID)]
euk_genus_diet_norm@sam_data$molecular_sex <- gorilla_sexing$`Sex assigned molecularly (using sexassign)`[match(sample_names(euk_genus_diet_norm), gorilla_sexing$Sample_ID)]

#Add age info to phyloseq metadata
euk_genus_diet@sam_data$age <- gorilla_age$Age[match(sample_names(euk_genus_diet), gorilla_age$Sample_ID)]
euk_genus_diet_norm@sam_data$age <- gorilla_age$Age[match(sample_names(euk_genus_diet_norm), gorilla_age$Sample_ID)]

####Model 3 - 1. Read count, 2. Seq.centre, 3. museum_sex, 4. Spec.subspecies ####
print("1. Read count, 2. Seq.centre, 3. museum_sex, 4. Spec.subspecies -- with a subset containing sex assignments")

euk_genus_diet_clean<-euk_genus_diet@sam_data %>% data.frame() %>% filter(!is.na(museum_sex)) %>% rownames(.) %>% 
  prune_samples(.,euk_genus_diet) %>% 
  filter_taxa(., function(x) sum(x) > 0, TRUE) %>% 
  prune_samples(sample_names(.)[!sample_names(.) %in% c("ABM006","DJA002","EBO001","GDC005")],.)

read_count <- prune_samples(euk_genus_diet_clean@sam_data %>% data.frame() %>% filter(!is.na(museum_sex)) %>% rownames(.),euk_genus_diet_clean)@sam_data$readcount.m.before.Kraken %>% as.numeric
seq_centre <- prune_samples(euk_genus_diet_clean@sam_data %>% data.frame() %>% filter(!is.na(museum_sex)) %>% rownames(.),euk_genus_diet_clean)@sam_data$Seq.centre %>% as.factor
spec.subspecies <- prune_samples(euk_genus_diet_clean@sam_data %>% data.frame() %>% filter(!is.na(museum_sex)) %>% rownames(.),euk_genus_diet_clean)@sam_data$Spec.subspecies %>% as.factor
sex <- prune_samples(euk_genus_diet_clean@sam_data %>% data.frame() %>% filter(!is.na(museum_sex)) %>% rownames(.),euk_genus_diet_clean)@sam_data$museum_sex %>% as.factor

model3_jaccard <- adonis2(t(abundances(subset_samples(euk_genus_diet_clean, !is.na(museum_sex)))) ~ read_count + seq_centre + sex + spec.subspecies, permutations = 10000, method = "jaccard")
model3_jaccard

cat(round(model3_jaccard$R2,3),sep = "\n")
cat(round(model3_jaccard$`Pr(>F)`,3),sep = "\n")

euk_genus_diet_norm_clean<-euk_genus_diet_norm@sam_data %>% data.frame() %>% filter(!is.na(museum_sex)) %>% rownames(.) %>% 
  prune_samples(.,euk_genus_diet) %>% 
  filter_taxa(., function(x) sum(x) > 0, TRUE) %>% 
  prune_samples(sample_names(.)[!sample_names(.) %in% c("ABM006","DJA002","EBO001","GDC005")],.)

read_count <- prune_samples(euk_genus_diet_norm_clean@sam_data %>% data.frame() %>% filter(!is.na(museum_sex)) %>% rownames(.),euk_genus_diet_norm_clean)@sam_data$readcount.m.before.Kraken %>% as.numeric
seq_centre <- prune_samples(euk_genus_diet_norm_clean@sam_data %>% data.frame() %>% filter(!is.na(museum_sex)) %>% rownames(.),euk_genus_diet_norm_clean)@sam_data$Seq.centre %>% as.factor
spec.subspecies <- prune_samples(euk_genus_diet_norm_clean@sam_data %>% data.frame() %>% filter(!is.na(museum_sex)) %>% rownames(.),euk_genus_diet_norm_clean)@sam_data$Spec.subspecies %>% as.factor
sex <- prune_samples(euk_genus_diet_norm_clean@sam_data %>% data.frame() %>% filter(!is.na(museum_sex)) %>% rownames(.),euk_genus_diet_norm_clean)@sam_data$museum_sex %>% as.factor

model3_clr <- adonis2(t(abundances(subset_samples(euk_genus_diet_clean, !is.na(museum_sex)))) ~ read_count + seq_centre + sex + spec.subspecies, permutations = 10000, method = "euclidean")
model3_clr

cat(round(model3_clr$R2,3),sep = "\n")
cat(round(model3_clr$`Pr(>F)`,3),sep = "\n")

print("1. Read count, 2. Seq.centre, 3. molecular_sex, 4. Spec.subspecies -- with a subset containing sex assignments")

euk_genus_diet_clean<-euk_genus_diet@sam_data %>% data.frame() %>% filter(!is.na(molecular_sex)) %>% rownames(.) %>% 
  prune_samples(.,euk_genus_diet) %>% 
  filter_taxa(., function(x) sum(x) > 0, TRUE) %>% 
  prune_samples(sample_names(.)[!sample_names(.) %in% c("ABM006","DJA002","EBO001","GDC005")],.)

read_count <- prune_samples(euk_genus_diet_clean@sam_data %>% data.frame() %>% filter(!is.na(molecular_sex)) %>% rownames(.),euk_genus_diet_clean)@sam_data$readcount.m.before.Kraken %>% as.numeric
seq_centre <- prune_samples(euk_genus_diet_clean@sam_data %>% data.frame() %>% filter(!is.na(molecular_sex)) %>% rownames(.),euk_genus_diet_clean)@sam_data$Seq.centre %>% as.factor
spec.subspecies <- prune_samples(euk_genus_diet_clean@sam_data %>% data.frame() %>% filter(!is.na(molecular_sex)) %>% rownames(.),euk_genus_diet_clean)@sam_data$Spec.subspecies %>% as.factor
sex <- prune_samples(euk_genus_diet_clean@sam_data %>% data.frame() %>% filter(!is.na(molecular_sex)) %>% rownames(.),euk_genus_diet_clean)@sam_data$molecular_sex %>% as.factor

model3_jaccard <- adonis2(t(abundances(subset_samples(euk_genus_diet_clean, !is.na(molecular_sex)))) ~ read_count + seq_centre + sex + spec.subspecies, permutations = 10000, method = "jaccard")
model3_jaccard

cat(round(model3_jaccard$R2,3),sep = "\n")
cat(round(model3_jaccard$`Pr(>F)`,3),sep = "\n")

euk_genus_diet_norm_clean<-euk_genus_diet_norm@sam_data %>% data.frame() %>% filter(!is.na(molecular_sex)) %>% rownames(.) %>% 
  prune_samples(.,euk_genus_diet) %>% 
  filter_taxa(., function(x) sum(x) > 0, TRUE) %>% 
  prune_samples(sample_names(.)[!sample_names(.) %in% c("ABM006","DJA002","EBO001","GDC005")],.)

read_count <- prune_samples(euk_genus_diet_norm_clean@sam_data %>% data.frame() %>% filter(!is.na(molecular_sex)) %>% rownames(.),euk_genus_diet_norm_clean)@sam_data$readcount.m.before.Kraken %>% as.numeric
seq_centre <- prune_samples(euk_genus_diet_norm_clean@sam_data %>% data.frame() %>% filter(!is.na(molecular_sex)) %>% rownames(.),euk_genus_diet_norm_clean)@sam_data$Seq.centre %>% as.factor
spec.subspecies <- prune_samples(euk_genus_diet_norm_clean@sam_data %>% data.frame() %>% filter(!is.na(molecular_sex)) %>% rownames(.),euk_genus_diet_norm_clean)@sam_data$Spec.subspecies %>% as.factor
sex <- prune_samples(euk_genus_diet_norm_clean@sam_data %>% data.frame() %>% filter(!is.na(molecular_sex)) %>% rownames(.),euk_genus_diet_norm_clean)@sam_data$molecular_sex %>% as.factor

model3_clr <- adonis2(t(abundances(subset_samples(euk_genus_diet_clean, !is.na(molecular_sex)))) ~ read_count + seq_centre + sex + spec.subspecies, permutations = 10000, method = "euclidean")
model3_clr

cat(round(model3_clr$R2,3),sep = "\n")
cat(round(model3_clr$`Pr(>F)`,3),sep = "\n")

print("1. Read count, 2. Seq.centre, 3. age, 4. Spec.subspecies -- with a subset containing age assignments")

euk_genus_diet_clean<-euk_genus_diet@sam_data %>% data.frame() %>% filter(!is.na(age)) %>% rownames(.) %>% 
  prune_samples(.,euk_genus_diet) %>% 
  filter_taxa(., function(x) sum(x) > 0, TRUE) %>% 
  prune_samples(sample_names(.)[!sample_names(.) %in% c("ABM006","DJA002","EBO001","GDC005")],.)

read_count <- prune_samples(euk_genus_diet_clean@sam_data %>% data.frame() %>% filter(!is.na(age)) %>% rownames(.),euk_genus_diet_clean)@sam_data$readcount.m.before.Kraken %>% as.numeric
seq_centre <- prune_samples(euk_genus_diet_clean@sam_data %>% data.frame() %>% filter(!is.na(age)) %>% rownames(.),euk_genus_diet_clean)@sam_data$Seq.centre %>% as.factor
spec.subspecies <- prune_samples(euk_genus_diet_clean@sam_data %>% data.frame() %>% filter(!is.na(age)) %>% rownames(.),euk_genus_diet_clean)@sam_data$Spec.subspecies %>% as.factor
age <- prune_samples(euk_genus_diet_clean@sam_data %>% data.frame() %>% filter(!is.na(age)) %>% rownames(.),euk_genus_diet_clean)@sam_data$age %>% as.factor

model4_jaccard <- adonis2(t(abundances(subset_samples(euk_genus_diet_clean,!is.na(age)))) ~ read_count + seq_centre + age + spec.subspecies, permutations = 10000, method = "jaccard")
model4_jaccard

cat(round(model4_jaccard$R2,3),sep = "\n")
cat(round(model4_jaccard$`Pr(>F)`,3),sep = "\n")

euk_genus_diet_norm_clean<-euk_genus_diet_norm@sam_data %>% data.frame() %>% filter(!is.na(age)) %>% rownames(.) %>% 
  prune_samples(.,euk_genus_diet) %>% 
  filter_taxa(., function(x) sum(x) > 0, TRUE) %>% 
  prune_samples(sample_names(.)[!sample_names(.) %in% c("ABM006","DJA002","EBO001","GDC005")],.)

read_count <- prune_samples(euk_genus_diet_norm_clean@sam_data %>% data.frame() %>% filter(!is.na(age)) %>% rownames(.),euk_genus_diet_norm_clean)@sam_data$readcount.m.before.Kraken %>% as.numeric
seq_centre <- prune_samples(euk_genus_diet_norm_clean@sam_data %>% data.frame() %>% filter(!is.na(age)) %>% rownames(.),euk_genus_diet_norm_clean)@sam_data$Seq.centre %>% as.factor
spec.subspecies <- prune_samples(euk_genus_diet_norm_clean@sam_data %>% data.frame() %>% filter(!is.na(age)) %>% rownames(.),euk_genus_diet_norm_clean)@sam_data$Spec.subspecies %>% as.factor
sex <- prune_samples(euk_genus_diet_norm_clean@sam_data %>% data.frame() %>% filter(!is.na(age)) %>% rownames(.),euk_genus_diet_norm_clean)@sam_data$age %>% as.factor

model4_clr <- adonis2(t(abundances(subset_samples(euk_genus_diet_clean,!is.na(age)))) ~ read_count + seq_centre + age + spec.subspecies, permutations = 10000, method = "euclidean")
model4_clr

cat(round(model4_clr$R2,3),sep = "\n")
cat(round(model4_clr$`Pr(>F)`,3),sep = "\n")

print("1. Read count, 2. Seq.centre, 3. museum_sex, 4. age, 5. Spec.subspecies -- with a subset for museum_sex and age assignments")

euk_genus_diet_clean<-euk_genus_diet@sam_data %>% data.frame() %>% filter(!is.na(museum_sex) & !is.na(age)) %>% rownames(.) %>% 
  prune_samples(.,euk_genus_diet) %>% 
  filter_taxa(., function(x) sum(x) > 0, TRUE) %>% 
  prune_samples(sample_names(.)[!sample_names(.) %in% c("ABM006","DJA002","EBO001","GDC005")],.)

read_count <- prune_samples(euk_genus_diet_clean@sam_data %>% data.frame() %>% filter(!is.na(museum_sex) & !is.na(age)) %>% rownames(.),euk_genus_diet_clean)@sam_data$readcount.m.before.Kraken %>% as.numeric
seq_centre <- prune_samples(euk_genus_diet_clean@sam_data %>% data.frame() %>% filter(!is.na(museum_sex) & !is.na(age)) %>% rownames(.),euk_genus_diet_clean)@sam_data$Seq.centre %>% as.factor
spec.subspecies <- prune_samples(euk_genus_diet_clean@sam_data %>% data.frame() %>% filter(!is.na(museum_sex) & !is.na(age)) %>% rownames(.),euk_genus_diet_clean)@sam_data$Spec.subspecies %>% as.factor
age <- prune_samples(euk_genus_diet_clean@sam_data %>% data.frame() %>% filter(!is.na(museum_sex) & !is.na(age)) %>% rownames(.),euk_genus_diet_clean)@sam_data$age %>% as.factor
sex <- prune_samples(euk_genus_diet_clean@sam_data %>% data.frame() %>% filter(!is.na(museum_sex) & !is.na(age)) %>% rownames(.),euk_genus_diet_clean)@sam_data$museum_sex %>% as.factor

model5_jaccard <- adonis2(t(abundances(subset_samples(euk_genus_diet_clean, !is.na(museum_sex) & !is.na(age)))) ~ read_count + seq_centre + sex + age + spec.subspecies, permutations = 10000, method = "jaccard")
model5_jaccard

cat(round(model5_jaccard$R2,3),sep = "\n")
cat(round(model5_jaccard$`Pr(>F)`,3),sep = "\n")

euk_genus_diet_norm_clean<-euk_genus_diet_norm@sam_data %>% data.frame() %>% filter(!is.na(museum_sex) & !is.na(age)) %>% rownames(.) %>% 
  prune_samples(.,euk_genus_diet) %>% 
  filter_taxa(., function(x) sum(x) > 0, TRUE) %>% 
  prune_samples(sample_names(.)[!sample_names(.) %in% c("ABM006","DJA002","EBO001","GDC005")],.)

read_count <- prune_samples(euk_genus_diet_norm_clean@sam_data %>% data.frame() %>% filter(!is.na(museum_sex) & !is.na(age)) %>% rownames(.),euk_genus_diet_norm_clean)@sam_data$readcount.m.before.Kraken %>% as.numeric
seq_centre <- prune_samples(euk_genus_diet_norm_clean@sam_data %>% data.frame() %>% filter(!is.na(museum_sex) & !is.na(age)) %>% rownames(.),euk_genus_diet_norm_clean)@sam_data$Seq.centre %>% as.factor
spec.subspecies <- prune_samples(euk_genus_diet_norm_clean@sam_data %>% data.frame() %>% filter(!is.na(museum_sex) & !is.na(age)) %>% rownames(.),euk_genus_diet_norm_clean)@sam_data$Spec.subspecies %>% as.factor
age <- prune_samples(euk_genus_diet_norm_clean@sam_data %>% data.frame() %>% filter(!is.na(museum_sex) & !is.na(age)) %>% rownames(.),euk_genus_diet_norm_clean)@sam_data$age %>% as.factor
sex <- prune_samples(euk_genus_diet_norm_clean@sam_data %>% data.frame() %>% filter(!is.na(museum_sex) & !is.na(age)) %>% rownames(.),euk_genus_diet_norm_clean)@sam_data$museum_sex %>% as.factor

model5_clr <- adonis2(t(abundances(subset_samples(euk_genus_diet_norm_clean, !is.na(museum_sex) & !is.na(age)))) ~ read_count + seq_centre + sex + age + spec.subspecies, permutations = 10000, method = "euclidean")
model5_clr

cat(round(model5_clr$R2,3),sep = "\n")
cat(round(model5_clr$`Pr(>F)`,3),sep = "\n")

print("1. Read count, 2. Seq.centre, 3. molecular_sex, 4. age, 5. Spec.subspecies -- with a subset for molecular_sex and age assignments")

euk_genus_diet_clean<-euk_genus_diet@sam_data %>% data.frame() %>% filter(!is.na(molecular_sex) & !is.na(age)) %>% rownames(.) %>% 
  prune_samples(.,euk_genus_diet) %>% 
  filter_taxa(., function(x) sum(x) > 0, TRUE) %>% 
  prune_samples(sample_names(.)[!sample_names(.) %in% c("ABM006","DJA002","EBO001","GDC005")],.)

read_count <- prune_samples(euk_genus_diet_clean@sam_data %>% data.frame() %>% filter(!is.na(molecular_sex) & !is.na(age)) %>% rownames(.),euk_genus_diet_clean)@sam_data$readcount.m.before.Kraken %>% as.numeric
seq_centre <- prune_samples(euk_genus_diet_clean@sam_data %>% data.frame() %>% filter(!is.na(molecular_sex) & !is.na(age)) %>% rownames(.),euk_genus_diet_clean)@sam_data$Seq.centre %>% as.factor
spec.subspecies <- prune_samples(euk_genus_diet_clean@sam_data %>% data.frame() %>% filter(!is.na(molecular_sex) & !is.na(age)) %>% rownames(.),euk_genus_diet_clean)@sam_data$Spec.subspecies %>% as.factor
age <- prune_samples(euk_genus_diet_clean@sam_data %>% data.frame() %>% filter(!is.na(molecular_sex) & !is.na(age)) %>% rownames(.),euk_genus_diet_clean)@sam_data$age %>% as.factor
sex <- prune_samples(euk_genus_diet_clean@sam_data %>% data.frame() %>% filter(!is.na(molecular_sex) & !is.na(age)) %>% rownames(.),euk_genus_diet_clean)@sam_data$molecular_sex %>% as.factor

model5_jaccard <- adonis2(t(abundances(subset_samples(euk_genus_diet_clean, !is.na(molecular_sex) & !is.na(age)))) ~ read_count + sex + age + spec.subspecies, permutations = 10000, method = "jaccard")
model5_jaccard

cat(round(model5_jaccard$R2,3),sep = "\n")
cat(round(model5_jaccard$`Pr(>F)`,3),sep = "\n")

euk_genus_diet_norm_clean<-euk_genus_diet_norm@sam_data %>% data.frame() %>% filter(!is.na(molecular_sex) & !is.na(age)) %>% rownames(.) %>% 
  prune_samples(.,euk_genus_diet) %>% 
  filter_taxa(., function(x) sum(x) > 0, TRUE) %>% 
  prune_samples(sample_names(.)[!sample_names(.) %in% c("ABM006","DJA002","EBO001","GDC005")],.)

read_count <- prune_samples(euk_genus_diet_norm_clean@sam_data %>% data.frame() %>% filter(!is.na(molecular_sex) & !is.na(age)) %>% rownames(.),euk_genus_diet_norm_clean)@sam_data$readcount.m.before.Kraken %>% as.numeric
seq_centre <- prune_samples(euk_genus_diet_norm_clean@sam_data %>% data.frame() %>% filter(!is.na(molecular_sex) & !is.na(age)) %>% rownames(.),euk_genus_diet_norm_clean)@sam_data$Seq.centre %>% as.factor
spec.subspecies <- prune_samples(euk_genus_diet_norm_clean@sam_data %>% data.frame() %>% filter(!is.na(molecular_sex) & !is.na(age)) %>% rownames(.),euk_genus_diet_norm_clean)@sam_data$Spec.subspecies %>% as.factor
age <- prune_samples(euk_genus_diet_norm_clean@sam_data %>% data.frame() %>% filter(!is.na(molecular_sex) & !is.na(age)) %>% rownames(.),euk_genus_diet_norm_clean)@sam_data$age %>% as.factor
sex <- prune_samples(euk_genus_diet_norm_clean@sam_data %>% data.frame() %>% filter(!is.na(molecular_sex) & !is.na(age)) %>% rownames(.),euk_genus_diet_norm_clean)@sam_data$molecular_sex %>% as.factor

model5_clr <- adonis2(t(abundances(subset_samples(euk_genus_diet_norm_clean, !is.na(molecular_sex) & !is.na(age)))) ~ read_count + sex + age + spec.subspecies, permutations = 10000, method = "euclidean")
model5_clr

cat(round(model5_clr$R2,3),sep = "\n")
cat(round(model5_clr$`Pr(>F)`,3),sep = "\n")
