#Start logging
sink(file = "log10_ancom.txt")

library(readr)
library(nlme)
library(tidyverse)
library(ggplot2)
library(compositions)
library(microbiome)
library(ggpubr)
library(readr)
library(tidyverse)
library(reshape2)
library(gridExtra)
library(dplyr)
library(phyloseq)
source("T3_community-level/ancom-functions.R")

load(".RData")

#Community-level taxonomic analysis - Script 10
#Detecting differentially abundant taxa with ANCOM

#Based on Adrian's ancom.R script
# feature_table_pre_process
# feature_table_pre_process(feature_table, meta_data, sample_var, group_var = NULL, out_cut = 0.05, zero_cut = 0.90, lib_cut, neg_lb)

  # Arguments
      # feature_table: Data frame or matrix representing observed OTU/SV table with taxa in rows (rownames) and samples in columns (colnames). Note that this is the absolute abundance table, do not transform it to relative abundance table (where the column totals are equal to 1).
      # meta_data: Data frame or matrix of all variables and covariates of interest.
      # sample_var: Character. The name of column storing sample IDs.
      # group_var: Character. The name of the group indicator. group_var is required for detecting structural zeros and outliers. For the definitions of different zeros (structural zero, outlier zero, and sampling zero), please refer to ANCOM-II.
      # out_cut Numerical fraction between 0 and 1. For each taxon, observations with proportion of mixture distribution less than out_cut will be detected as outlier zeros; while observations with proportion of mixture distribution greater than 1 - out_cut will be detected as outlier values.
      # zero_cut: Numerical fraction between 0 and 1. Taxa with proportion of zeroes greater than zero_cut are not included in the analysis.
      # lib_cut: Numeric. Samples with library size less than lib_cut are not included in the analysis.
      # neg_lb: Logical. TRUE indicates a taxon would be classified as a structural zero in the corresponding experimental group using its asymptotic lower bound.
        # More specifically, neg_lb = TRUE indicates you are using both criteria stated in section 3.2 of ANCOM-II to detect structural zeros;
        # Otherwise, neg_lb = FALSE will only use the equation 1 in section 3.2 of ANCOM-II for declaring structural zeros.

  # Values
    # feature_table: A data frame of pre-processed OTU table.
    # meta_data: A data frame of pre-processed metadata.
    # structure_zeros: A matrix consists of 0 and 1s with 1 indicating the taxon is identified as a structural zero in the corresponding group.

# ANCOM main function
# ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, alpha, adj_formula, rand_formula, ...)
  # Arguments
      # feature_table: Data frame representing OTU/SV table with taxa in rows (rownames) and samples in columns (colnames). It can be the output value from feature_table_pre_process. Note that this is the absolute abundance table, do not transform it to relative abundance table (where the column totals are equal to 1).
      # meta_data: Data frame of variables. Can be the output value from feature_table_pre_process.
      # struc_zero: A matrix consists of 0 and 1s with 1 indicating the taxon is identified as a structural zero in the corresponding group. Can be the output value from feature_table_pre_process.
      # main_var: Character. The name of the main variable of interest. ANCOM v2.1 currently supports categorical main_var.
      # p_adjust_method: Character. Specifying the method to adjust p-values for multiple comparisons. Default is bBHb (Benjamini-Hochberg procedure).
# alpha: Level of significance. Default is 0.05.
# adj_formula: Character string representing the formula for adjustment (see example).
# rand_formula: Character string representing the formula for random effects in lme (see example).

#### Step 1: Data preprocessing ####

# group_var 
# used to detect structural zeros

# out_cut
# default is 0.05
# if zero counts make up of less than 5% of the proposed Gaussian mixture model, 
# it will be detected as outlier zeros and replaced with NA.
# If you believe your dataset is exempt from erroneous data entries, you can also specify out_cut = 0 to disable this detection.

# zero_cut 
# used for filtering non-informative taxa
# with 100 samples and zero_cut = 0.90 (default value), taxa with more than 90 zero entries out of 100 taxa will be discarded

# lib_cut 
# used for filtering samples
# with lib_cut = 1000, any samples with library size (total observed counts) less than 1000 will be discarded.
# library size varies a lot across different studies, some may have a lot of samples with library size less than 1000.
# In such cases, sticking with the default value will lose a lot of power.
# If you do not want to filter out any sample based on library sizes, you can simply set lib_cut = 0 to disable this function.

feature_table <- data.frame(otu_table(spe_data_final))
meta_data <- data.frame(sample_data(spe_data_final))
#Include sample IDs as a column
meta_data$SampleID <- rownames(meta_data)
sample_var = "SampleID"
group_var = "Spec.subspecies"
out_cut = 0.05
zero_cut = 0.90
lib_cut = 0
neg_lb = FALSE

#prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var,out_cut, zero_cut, lib_cut, neg_lb)
#feature_table = prepro$feature_table # Preprocessed feature table
#meta_data = prepro$meta_data # Preprocessed metadata
#struc_zero_spe = prepro$structure_zeros # Structural zero info

#### Step 2: ANCOM ####
main_var = "Spec.subspecies" # taxa that are uniquely found in one disease state will be labelled as structural zeros
p_adj_method = "BH"
alpha = 0.05
adj_formula = "readcount.m.before.Kraken" #the read depth before the analysis will be used as a covariate
rand_formula = NULL

#species_ancom_res = ANCOM(feature_table, meta_data, struc_zero_spe, main_var, p_adj_method,alpha, adj_formula, rand_formula)

# save results
write_csv(cbind(species_ancom_res$out,struc_zero_spe), "T3_community-level/ancom_results.csv")

save.image()

#### Step 3: Volcano Plot ####
# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero_spe), nrow(feature_table), sum(apply(struc_zero_spe, 1, sum) == 0))

# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1), 0.5 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6", "detected_0.5")

# Annotation data
dat_ann = data.frame(x = min(species_ancom_res$fig$data$x), y = cut_off["detected_0.9"], label = "W[0.9]")

ancom_fig = species_ancom_res$fig +
  geom_hline(yintercept = cut_off["detected_0.9"], linetype = "dashed") +
  geom_text(data = dat_ann, aes(x = x, y = y, label = label),
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)+
  theme(legend.position = "top")

ggsave(plot=ancom_fig,filename = "T3_community-level/ancom_fig.png")

#### Rerun ANCOM for sequencing facility ####
feature_table <- data.frame(otu_table(spe_data_final))
meta_data <- data.frame(sample_data(spe_data_final))
#Include sample IDs as a column
meta_data$SampleID <- rownames(meta_data)
sample_var = "SampleID"
group_var = "Seq.centre" #differentially abundant taxa between sequencing centres will likely indicate contamination 
out_cut = 0.05
zero_cut = 0.90
lib_cut = 0
neg_lb = FALSE

#prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var,out_cut, zero_cut, lib_cut, neg_lb)
#feature_table = prepro$feature_table # Preprocessed feature table
#meta_data = prepro$meta_data # Preprocessed metadata
#struc_zero = prepro$structure_zeros # Structural zero info

#### Step 2: ANCOM ####
main_var = "Seq.centre" # taxa that are uniquely found in one seq facility will be marked as zero
p_adj_method = "BH"
alpha = 0.05
adj_formula = "Spec.subspecies" #Host subspecies will be used as a covariate
rand_formula = NULL

#species_ancom_res_seqc = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method,alpha, adj_formula, rand_formula)

save.image()

# save results
write_csv(cbind(species_ancom_res_seqc$out,struc_zero), "T3_community-level/ancom_results_seqc.csv")

#### Step 3: Volcano Plot ####
# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))

# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1), 0.5 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6", "detected_0.5")

# Annotation data
dat_ann = data.frame(x = min(species_ancom_res_seqc$fig$data$x), y = cut_off["detected_0.7"], label = "W[0.7]")

ancom_fig_seqc = species_ancom_res_seqc$fig +
  geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") +
  geom_text(data = dat_ann, aes(x = x, y = y, label = label),
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)+
  theme(legend.position = "top")

ggsave(plot=ancom_fig,filename = "T3_community-level/ancom_fig_seqc.png")


#### Look into differentially abundant species ####
#For host subspecies
ancom_species <- read.table("T3_community-level/ancom_results.csv", header = TRUE, sep = ",", dec = ".") %>% filter(detected_0.9==TRUE) %>% select(taxa_id)

#Add species names
ancom_species$taxa_names <- taxonomy_species[match(ancom_species$taxa_id, rownames(taxonomy_species)),8]

#
print("How many differentially abundant taxa between host subspecies?")
nrow(ancom_species)

#For sequencing facility
ancom_species_seqc <- read.table("T3_community-level/ancom_results_seqc.csv", header = TRUE, sep = ",", dec = ".") %>% filter(detected_0.9==TRUE) %>% select(taxa_id)

#Add species names
ancom_species_seqc$taxa_names <- taxonomy_species[match(ancom_species_seqc$taxa_id, rownames(taxonomy_species)),8]

#
print("How many differentially abundant taxa between sequencing facilities?")
nrow(ancom_species_seqc)


##Remove these that are differentially abundant between sequencing facilities at the genus level
ancom_species <-
  ancom_species %>%
  filter(!(str_remove(taxa_names, " .*") %in% str_remove(ancom_species_seqc$taxa_names, " .*")))

#
print("How many differentially abundant taxa remain between host subspecies after removing those that are due to the sequencing facility?")
nrow(ancom_species)

#Mark those that are structural zeros
ancom_species$is.zero <- ifelse(ancom_species$taxa_id %in% rownames(struc_zero_spe[which(rowSums(struc_zero_spe) > 0),]),
                                TRUE, FALSE)

print("How many of there are stuctural zeroes?")
table(ancom_species$is.zero)

#hominid
print("Which core hominid taxa are structural zeroes")
table(ancom_species$is.zero[which(ancom_species$taxa_names %in% core_micr$Taxon)])


#HOMD
print("Which HOMD hominid taxa are structural zeroes")
table(ancom_species$is.zero[which(ancom_species$taxa_id %in% homd$NCBI_taxon_id)])

#### Heatmap on all taxa ####
ancom_species_table <- otu_table(subset_taxa(spe_data_final_norm, taxa_names(spe_data_final_norm) %in% ancom_species$taxa_id))

#Set 0 values (the most negative value of CLR-normalized abundances) to min(clr.pseudocounts)*2 instead of NA (which is what Jaelle did). This is because the heatmap function that I use doesn't recognize NAs
clr.pseudozeros <- sapply(colnames(ancom_species_table), function(x){min(ancom_species_table[,x])})
for (s in colnames(ancom_species_table)) {
  ancom_species_table[which(ancom_species_table[,s] == clr.pseudozeros[s]),s] <- min(clr.pseudozeros)*2
}

#Melt
ancom_species_table <- reshape2::melt(ancom_species_table)
colnames(ancom_species_table) <- c("taxon", "sample", "clr-abundance")

#Add subspecies column
ancom_species_table$subspecies <- metadata$Spec.subspecies[match(ancom_species_table$sample, rownames(metadata))]

#Add taxon names column
ancom_species_table$taxon_name <- taxonomy_species[,8][match(ancom_species_table$taxon, rownames(taxonomy_species))]

#Add taxon order column: 
ancom_species_table <- 
  ancom_species_table %>% mutate(taxon_order = taxonomy_species[,5][match(ancom_species_table$taxon, rownames(taxonomy_species))])

#The following can be used if there are too many different orders to be displayed properly
#For infrequent orders (less that <5 species in the ancom table) change order to Other
ancom_species_table <- 
  ancom_species_table %>% left_join(as.data.frame((table(ancom_species_table$taxon_order)/nsamples(spe_data_final))), by=c("taxon_order" = "Var1")) %>% 
  mutate(taxon_order=ifelse(Freq<3, "Other", taxon_order))

#Add taxon phylum column
ancom_species_table$taxon_phylum <- taxonomy_species[,3][match(ancom_species_table$taxon, rownames(taxonomy_species))]

#Also domain - to split bacteria from archaea
#Add taxon phylum column
ancom_species_table$domain <- taxonomy_species[,1][match(ancom_species_table$taxon, rownames(taxonomy_species))]

#Order table based on subspecies and microbial order (also domain)
ancom_species_table$subspecies <- factor(ancom_species_table$subspecies, levels = c("gorilla", "graueri", "beringei"))

#levels <- ancom_species_table[order(ancom_species_table$Freq),] %>% pull(taxon_order) %>% unique %>% as.character
#levels <- levels[-which(levels=="Other")] %>% append("Other")
#levels <- levels[c(length(levels), 1:length(levels)-1)]

#ancom_species_table$taxon_order <- factor(ancom_species_table$taxon_order, levels = levels)
ancom_species_table <- ancom_species_table[order(ancom_species_table$subspecies, ancom_species_table$taxon_order),]

#Labeller for the facet labels
facet.labels <- list("gorilla"="western lowland", "graueri"="Grauer's", "beringei"="mountain")
facet.labeller <- function(variable,value){
  return(facet.labels[value])
}

#Plot
#This version be used if there are too many taxa for all labels to be displayed properly
#diff_abund_heat <-
#  heat(ancom_species_table, Yvar="taxon_name", fill="clr-abundance", Xvar = "sample", order.cols = FALSE, order.rows = FALSE)+
#  scale_fill_gradient2(low="blue", mid="white", high="red",
#    guide = guide_colorbar(barheight = 10, title = "clr-normalized\nabundance", draw.ulim=FALSE),
#    limits=c(min(clr.pseudozeros), NA), na.value="grey") +
#  facet_grid(~subspecies, scales="free", labeller=facet.labeller) +
#  theme(plot.title = element_text(size = 20, hjust = 0.8), axis.text.x = element_text(angle = 90, vjust = -0.2), 
#        axis.title = element_text(size=15), axis.text.y = element_blank())

#This version can be used if all taxa labels can be displayed
diff_abund_heat <-
  heat(ancom_species_table, Yvar="taxon_name", fill="clr-abundance", Xvar = "sample", order.cols = FALSE, order.rows = TRUE)+
  scale_fill_gradient2(low="blue", mid="white", high="red",
                       guide = guide_colorbar(barheight = 10, title = "clr-normalized\nabundance", draw.ulim=FALSE),
                       limits=c(min(clr.pseudozeros), NA), na.value="grey") +
  facet_grid(~subspecies, scales="free", labeller=facet.labeller) +
  theme(plot.title = element_text(size = 20, hjust = 0.8), axis.text.x = element_text(angle = 90, vjust = -0.2), 
        axis.title = element_text(size=15), axis.text.y = element_text())

levels <- ancom_species_table %>% pull(taxon_order) %>% unique %>% as.character
#Bring 'Other' to be the first in order
levels <- append("Other", levels[-which(levels=="Other")])

#### Colour the labels by order ####
heat_ancom_ylabs <- levels(diff_abund_heat$data$YYYY)
heat_ancom_ylabs <- cbind(heat_ancom_ylabs, as.character(ancom_species_table$taxon_order[match(heat_ancom_ylabs, ancom_species_table$taxon_name)]))
heat_ancom_ylabs <- as.data.frame(heat_ancom_ylabs)
colnames(heat_ancom_ylabs) <- c("species", "order")
heat_ancom_ylabs$order <- factor(heat_ancom_ylabs$order, levels=levels)

#Create a palette
ancom_spe_palette <- c("#999999", "cadetblue1", "navy", "chartreuse3", "darkolivegreen", "gold","lightsalmon", "firebrick1",
                       "darkorchid3", "hotpink", "plum2", "dodgerblue", "khaki", "darkorange", "violetred")
names(ancom_spe_palette) <- levels(heat_ancom_ylabs$order)

#Get species-level palette
ancom_spe_palette.s <- heat_ancom_ylabs
ancom_spe_palette.s$colour <- ancom_spe_palette[match(ancom_spe_palette.s$order, names(ancom_spe_palette))]
ancom_spe_palette.s <- t(ancom_spe_palette.s[, -which(colnames(ancom_spe_palette.s)=="order")])

#Set taxon names as palette names
colnames(ancom_spe_palette.s) <- ancom_spe_palette.s[1,]
ancom_spe_palette.s <- ancom_spe_palette.s[-1,]

#Update plot
diff_abund_heat <-
  diff_abund_heat +  theme(axis.ticks.y = element_line(colour = ancom_spe_palette.s, size=2.5),
                           axis.ticks.length.y = unit(0.5, "cm"),
                           axis.text.y = element_text()) +
  ylab("Differentially abundant microbial species")
diff_abund_heat

ggsave(diff_abund_heat, file="diff_abund_heat.png", device="png", height=9, width=7)

#### Heatmap on differentially abundant oral taxa ####
#Save the unique taxa IDs of differentially abundant oral taxa
ancom_oral <- rbind(ancom_species[which(ancom_species$taxa_names %in% core_micr$Taxon),],ancom_species[which(ancom_species$taxa_id %in% homd$NCBI_taxon_id),])

rownames(ancom_oral) <- NULL
ancom_oral <- unique(ancom_oral)
rownames(ancom_oral) <- NULL

write.table(ancom_oral, "T3_community-level/ancom_species_oral.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

#Plot abundances
ancom_oral_table <- otu_table(subset_taxa(spe_data_final_norm, taxa_names(spe_data_final_norm) %in% ancom_oral$taxa_id))

#Set 0 values (the most negative value of CLR-normalized abundances) to min(clr.pseudocounts)*2 instead of NA (which is what Jaelle did). This is because the heatmap function that I use doesn't recognize NAs
clr.pseudozeros.o <- sapply(colnames(ancom_oral_table), function(x){min(ancom_oral_table[,x])})
for (s in colnames(ancom_oral_table)) {
  ancom_oral_table[which(ancom_oral_table[,s] == clr.pseudozeros.o[s]),s] <- min(clr.pseudozeros.o)*2
}


#Melt
ancom_oral_table <- reshape2::melt(ancom_oral_table)
colnames(ancom_oral_table) <- c("taxon", "sample", "clr-abundance")

#Add subspecies column
ancom_oral_table$subspecies <- metadata$Spec.subspecies[match(ancom_oral_table$sample, rownames(metadata))]

#Add taxon names column
ancom_oral_table$taxon_name <- taxonomy_species[,8][match(ancom_oral_table$taxon, rownames(taxonomy_species))]

#Order table based on subspecies
ancom_oral_table$subspecies <- factor(ancom_oral_table$subspecies, levels = c("gorilla", "graueri", "beringei"))
ancom_oral_table <- ancom_oral_table[order(ancom_oral_table$subspecies),]

#Plot
heat_oral <-
  heat(ancom_oral_table, Yvar="taxon_name", fill="clr-abundance", Xvar = "sample", order.cols = FALSE) +
  scale_fill_gradient2(low="blue",mid="white" ,high="red" ,
                       guide = guide_colorbar(barheight = 10, title = "clr-normalized\nabundance", draw.ulim=TRUE),
                       limits=c(min(clr.pseudozeros), NA), na.value="grey") +
  facet_grid(~subspecies, scales="free", labeller=facet.labeller) +
  theme(plot.title = element_text(size = 20, hjust = 0.8), axis.text.x = element_text(angle = 90, vjust = -0.2), 
        axis.title = element_text(size=15), axis.text.y = element_text(size=10))

heat_oral
#### Structural zeros: how strong is the signal? ####
#I want to check the abundances of taxa that are flagged as structural zeros in one or more subspecies
#Could a single or a few samples where the taxon in present be driving the abundance of the whole subspecies

#I will use the raw abundances for this
#ancom_zeros_table <- otu_table(subset_taxa(spe_data_final, taxa_names(spe_data_final) %in% ancom_species$taxa_id[ancom_species$is.zero==TRUE]))

#Melt the table for plotting
#ancom_zeros_table <- melt(ancom_zeros_table)
#colnames(ancom_zeros_table) <- c("taxon", "sample", "abundance")

#log transform abundance with pseudocount
#ancom_zeros_table$abundance <- log(ancom_zeros_table$abundance + 1)

#heat_struc_zeros <- heat(ancom_zeros_table, Xvar="sample", Yvar="taxon", fill="abundance") + 
#  theme(axis.title = element_text(size=15), 
#        axis.text.y = element_blank(), axis.ticks.y = element_blank())  +
#  scale_fill_gradient2(high = "white", mid="navyblue", low = "black", guide = guide_colorbar(barheight = 10, title = "log-transformed\nraw abundance"))

#ggsave(heat_struc_zeros, filename = "T3_community-level/heat_struc_zeros.png")

##Separate struc zeros per subspecies
print("How many structural zeroes are there in each subspecies?")
mountain_zeros_table <- otu_table(subset_taxa(spe_data_final, taxa_names(spe_data_final) %in% rownames(struc_zero_spe)[which(struc_zero_spe[,1]==1)]))
western_zeros_table <- otu_table(subset_taxa(spe_data_final, taxa_names(spe_data_final) %in% rownames(struc_zero_spe)[which(struc_zero_spe[,2]==1)]))
grauer_zeros_table <- otu_table(subset_taxa(spe_data_final, taxa_names(spe_data_final) %in% rownames(struc_zero_spe)[which(struc_zero_spe[,3]==1)]))

#For the struct zero table of subspecies A, I will calculate the mean abundance of each taxon within A
#(it should be close to zero, but might not be exactly zero)
#The I will calculate the difference of every sample in subspecies B and C, and this mean(A) for every taxon
comp_heat4struczero <- function(subspecies) {
  result <- list()
  #Take a subset of taxa that are flagged as structural zeros in the reference subspecies
  subspecies_zeros_table <- otu_table(subset_taxa(spe_data_final, taxa_names(spe_data_final) %in% rownames(struc_zero_spe)[which(struc_zero_spe[,paste("structural_zero (", subspecies, ")", sep = "")]==1)]))
  subspecies_zeros_table <- as.data.frame(as.matrix(subspecies_zeros_table))
  
  #Save number of taxa
  numtaxa <- nrow(subspecies_zeros_table)
  
  #Calculate average abundance of every taxon in the reference subspecies - most likely all zeros, but just in case
  subspecies_zeros_table$reference_means <- rowMeans(subspecies_zeros_table[which(sample_data(spe_data_final)$Spec.subspecies==subspecies)])
  
  nonref <- which(sample_data(spe_data_final)$Spec.subspecies!=subspecies)
  
  #Save number of the first subspecies in order (for plotting a line to separate the subspecies on the heatmap)
  line1 <- as.data.frame(table(sample_data(spe_data_final)$Spec.subspecies[nonref]))[1,2]
  line2 <- sum(as.data.frame(table(sample_data(spe_data_final)$Spec.subspecies[nonref]))[,2])
  #Only keep non reference samples (and the last column which contains the means)
  subspecies_zeros_table <- subspecies_zeros_table[,c(nonref, 44)]
  
  subspecies_zeros_table$taxon <- rownames(subspecies_zeros_table)
  
  #Melt for plotting
  subspecies_zeros_table <- melt(subspecies_zeros_table)
  
  colnames(subspecies_zeros_table) <- c("taxon", "sample", "log.abundance")
  subspecies_zeros_table$log.abundance <- log(subspecies_zeros_table$log.abundance + 1)
  
  #Add subspecies column
  subspecies_zeros_table$subspecies <- metadata$Spec.subspecies[match(subspecies_zeros_table$sample, rownames(metadata))]
  
  #Order table based on subspecies
  subspecies_zeros_table <- subspecies_zeros_table[order(subspecies_zeros_table$subspecies),]
  
  #Plot
  heat <- heat(subspecies_zeros_table, Yvar = "taxon", Xvar = "sample", fill = "log.abundance",
               order.rows = FALSE, order.cols = FALSE) + 
    theme(axis.title = element_text(size=15), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    scale_fill_gradient2(guide = guide_colorbar(barheight = 10, title = "log diff in\nraw abundance")) +
    ggtitle(paste("reference subspecies:", subspecies)) +
    annotate("segment", size = 1, x = line1+0.5, xend = line1+0.5, y = 0.5, yend = numtaxa + 0.5, colour = "black") +
    annotate("segment", size = 1, x = line2+.5, xend = line2+0.5, y = 0.5, yend = numtaxa + 0.5, colour = "black")
  
  
  hist <- gghistogram(subspecies_zeros_table[which(subspecies_zeros_table$sample!="reference_means"),], x="log.abundance", bins=50) +
    xlab("log diff in raw abundance") + ggtitle(paste("number of taxa:", numtaxa))
  
  result[["hist"]] <- hist
  result[["heat"]] <- heat
  return(result)
}

#Run for each subspecies
comp_plot_struc_zeros <- list()
for (subspecies in c("beringei", "gorilla", "graueri")) {
  res <- comp_heat4struczero(subspecies)
  comp_plot_struc_zeros[[paste("heatmap_", subspecies, sep="")]] <- res[["heat"]]
  comp_plot_struc_zeros[[paste("histogram_", subspecies, sep="")]] <- res[["hist"]]
}

ggsave(grid.arrange(grobs=comp_plot_struc_zeros, ncol=2), file="T3_community-level/comp_plot_struc_zeros.png", width=15, height=15)

#### Compare abundances of differentially abundant and total taxa ####

#Get taxa sums for all taxa in the final dataset
comp_ancomVStotal <- taxa_sums(spe_data_final)
comp_ancomVStotal <- as.data.frame(cbind(comp_ancomVStotal, ifelse(names(comp_ancomVStotal) %in% ancom_species$taxa_id, "Diff. abundant", "Not diff. abundant")))
colnames(comp_ancomVStotal) <- c("log_abundance", "type")
comp_ancomVStotal$log_abundance <- log(as.numeric(comp_ancomVStotal$log_abundance) + 1)

#Plot
ggsave(
  ggboxplot(comp_ancomVStotal, x="type", y="log_abundance", fill="grey") +
    annotate("text", x=0.8, y=10, label=paste0("n = ", sum(comp_ancomVStotal$type=="Not diff. abundant"))) +
    annotate("text", x=1.8, y=10.2, label=paste0("n = ", sum(comp_ancomVStotal$type=="Diff. abundant"))) +
    ylab("log-transformed abundance"),
  file="comp_ancomVStotal_boxplot.png",
  device="png")

ggsave(
  gghistogram(comp_ancomVStotal, x="log_abundance", fill="type"),
  file="comp_ancomVStotal_histogram.png",
  device="png")

save.image()
q()

#### Maaslin ####
library(Maaslin2)

print("Run Maaslin - to compare")
#Also Run Maaslin 2
fit_species <- Maaslin2(input_data = t(otu_table(spe_data_final)), input_metadata = as.data.frame(as.matrix(sample_data(spe_data_final))),
                        fixed_effects = "Spec.subspecies", plot_heatmap = TRUE, n = "clr",
                        output = "T3_community-level/species_maaslin_output")

#Filter the features with the biggest effect
significant_species <- read_tsv("T3_community-level/species_maaslin_output/significant_results.tsv")

#Only keep species that are significant according to adjusted p-value
significant_species <- significant_species[significant_species$qval<0.05,]

#Fix taxa names
significant_species$feature <- str_remove(significant_species$feature, "X")

#Compare with the taxa identified by ANCOM
print("Which taxa identified by ANCOM are also identified by Maaslin2")
species_ancom_res$out[which(species_ancom_res$out$taxa_id %in% significant_species$feature),]