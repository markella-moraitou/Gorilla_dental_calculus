#Start logging
sink(file = "log5_ancom.txt")


library(readr)
library(nlme)
library(tidyverse)
library(ggplot2)
library(compositions)
library(phyloseq)
library(microbiome)
library(ggpubr)
library(readr)
library(tidyverse)
library(reshape2)
library(scales)
library(cowplot)
source("/proj/sllstore2017021/nobackup/MARKELLA/T3_community-level/ancom-functions.R")

#Diet analysis - Script 5
#Look for differentially abundant taxa

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
# p_adjust_method: Character. Specifying the method to adjust p-values for multiple comparisons. Default is “BH” (Benjamini-Hochberg procedure).
# alpha: Level of significance. Default is 0.05.
# adj_formula: Character string representing the formula for adjustment (see example).
# rand_formula: Character string representing the formula for random effects in lme (see example).

#### GENUS-LEVEL ASSIGNMENTS ####
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

#Remove environmental controls from dataset
euk_genus_diet <- subset_samples(euk_genus_diet, !is.na(Spec.subspecies))
euk_genus_diet_norm <- euk_genus_diet
otu_table(euk_genus_diet_norm) <- otu_table(microbiome::transform(otu_table(euk_genus_diet_norm), transform = "clr"),
                                            taxa_are_rows = TRUE)


feature_table <- data.frame(otu_table(euk_genus_diet))
meta_data <- data.frame(sample_data(euk_genus_diet))
#Include sample IDs as a column
meta_data$SampleID <- rownames(meta_data)
sample_var = "SampleID"
group_var = "Spec.subspecies"
out_cut = 0.05
zero_cut = 0.90
lib_cut = 0
neg_lb = FALSE

prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var,out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

#### Step 2: ANCOM ####
main_var = "Spec.subspecies" # taxa that are uniquely found in one subspecies state will be labelled as structural zeros
p_adj_method = "BH"
alpha = 0.05
adj_formula = "readcount.m.before.Kraken"#the read depth before the analysis will be used as a covariate
rand_formula = NULL

euk_ancom_res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method,alpha, adj_formula, rand_formula)
save.image()

# save results
write_csv(cbind(euk_ancom_res$out,struc_zero), "/proj/sllstore2017021/nobackup/MARKELLA/D3_diet_stats/euk_ancom_results.csv")

#### Step 3: Volcano Plot ####
# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))

# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1), 0.5 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6", "detected_0.5")

# Annotation data
dat_ann = data.frame(x = min(euk_ancom_res$fig$data$x), y = cut_off["detected_0.9"], label = "W[0.9]")

euk_ancom_fig = euk_ancom_res$fig +
  geom_hline(yintercept = cut_off["detected_0.9"], linetype = "dashed") +
  geom_text(data = dat_ann, aes(x = x, y = y, label = label),
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)+
  theme(legend.position = "top")

euk_ancom_fig

ggsave(plot=ancom_fig,filename = "/proj/sllstore2017021/nobackup/MARKELLA/D3_diet_stats/euk_ancom_fig.png")


#### Look into differentially abundant species ####
ancom_euk <- as.data.frame(euk_ancom_res$out$taxa_id[which(euk_ancom_res$out$detected_0.9==TRUE)])

colnames(ancom_euk) <- "genus"

#Add family names
ancom_euk$family <- tax_table(euk_genus_diet)[match(ancom_euk$genus, taxa_names(euk_genus_diet)),6]
#Add phylum names
ancom_euk$phylum <- tax_table(euk_genus_diet)[match(ancom_euk$genus, taxa_names(euk_genus_diet)),3]

dim(ancom_euk)
#37 differentially abundant plant taxa

table(ancom_euk$phylum)
#36 streptophyta, 1 fungus

#Mark those that are structural zeros
ancom_euk$is.zero <- ifelse(ancom_euk$genus %in% rownames(struc_zero[which(rowSums(struc_zero) > 0),]),
                                TRUE, FALSE)

table(ancom_euk$is.zero)
#35 are structural zeros in at least one group, 2 are differentially abundant, but with a high W statist

#### Plot heatmap ####
#Plot abundances
ancom_euk_abund <- otu_table(subset_taxa(euk_genus_diet_norm, taxa_names(euk_genus_diet_norm) %in% ancom_euk$genus))

#Set 0 values (the most negative value of CLR-normalized abundances) to min(clr.pseudocounts)*2 instead of NA (which is what Jaelle did). This is because the heatmap function that I use doesn't recognize NAs
clr.pseudozeros.dd <- sapply(colnames(ancom_euk_abund), function(x){min(ancom_euk_abund[,x])})
for (s in colnames(ancom_euk_abund)) {
  ancom_euk_abund[which(ancom_euk_abund[,s] == clr.pseudozeros.dd[s]),s] <- min(clr.pseudozeros.dd)*2
}

#Melt
ancom_euk_abund <- reshape2::melt(ancom_euk_abund)
colnames(ancom_euk_abund) <- c("genus", "sample", "clr-abundance")

#Add host subspecies column
ancom_euk_abund$host_subspecies <- metadata$Spec.subspecies[match(ancom_euk_abund$sample, rownames(metadata))]

#Add phylum column
ancom_euk_abund$phylum_name <- taxonomy_euk[,3][match(ancom_euk_abund$genus, taxonomy_euk[,7])]

#Add family columns
ancom_euk_abund$family_name <- taxonomy_euk[,6][match(ancom_euk_abund$genus, taxonomy_euk[,7])]

#Order table based on host subspecies, phylum and family
ancom_euk_abund$host_subspecies <- factor(ancom_euk_abund$host_subspecies, levels = c("gorilla", "graueri", "beringei"))
ancom_euk_abund <- ancom_euk_abund[order(ancom_euk_abund$host_subspecies, ancom_euk_abund$phylum_name, ancom_euk_abund$family_name),]


#Plot

#Labeller for the facet labels
facet.labels <- list("gorilla"="western lowland", "graueri"="Grauer's", "beringei"="mountain")
facet.labeller <- function(variable,value){
  return(facet.labels[value])
}

heat_diet <- heat(ancom_euk_abund, Yvar="genus", fill="clr-abundance", Xvar = "sample", order.cols = FALSE, order.rows = FALSE) +
  theme(legend.position = "left", plot.title=element_text(hjust = 0.5)) +
  scale_fill_gradient2(low="blue", mid="white", high="red",
    guide = guide_colorbar(barheight = 2, title = "clr-normalized\nabundance", draw.ulim=FALSE),
    limits=c(min(clr.pseudozeros.dd), NA), na.value="grey") +
    facet_grid(~host_subspecies, scales="free", labeller=facet.labeller) +
  ggtitle("Abundances of differentially abundant gorilla dietary genera across samples")


#### Colour the labels by phylum ####
heat_diet_ylabs <- levels(heat_diet$data$YYYY)
heat_diet_ylabs <- cbind(heat_diet_ylabs, tax_table(euk_genus_diet)[match(heat_diet_ylabs, taxa_names(euk_genus_diet)),3])
heat_diet_ylabs <- as.data.frame(heat_diet_ylabs)
colnames(heat_diet_ylabs) <- c("genus", "phylum")
heat_diet_ylabs$phylum <- factor(heat_diet_ylabs$phylum)

#Also add family infomation (will be needed later)
heat_diet_ylabs$family <- tax_table(euk_genus_diet)[match(heat_diet_ylabs$genus, taxa_names(euk_genus_diet)),6]

#Save numbers of genera per family
heat_diet_fam_freqs <- table(heat_diet_ylabs$family)
heat_diet_fam_freqs <- as.data.frame(heat_diet_fam_freqs)
heat_diet_fam_freqs <- heat_diet_fam_freqs[match(unique(heat_diet_ylabs$family), heat_diet_fam_freqs$Var1),]

#Subset diet_palette.g to get colours for the taxa in this heatmap
diet_palette.g <- diet_palette_full.g[match(heat_diet_ylabs$genus, colnames(diet_palette_full.g))]
names(diet_palette.g) <- heat_diet_ylabs$genus
diet_palette.g <- t(diet_palette.g)

#Also draw lines to show genera that belong to the same family
#From the sorted table with the number of genera per family, get the cumulative sum
#to determine where every line will be along the y axis
heat_diet_fam_freqs$Freq <- cumsum(heat_diet_fam_freqs$Freq)


#Update plot
heat_diet <-
  heat_diet + theme(axis.text.y = element_text(colour = diet_palette.g),
                    legend.position = "top") +
  geom_hline(yintercept=heat_diet_fam_freqs$Freq + 0.5, size=0.5, colour="white")
 

#### Create a sidebar ####
#based on the families in the studies
ref_sidebar <- heat_diet_ylabs[,c("genus", "family")]
ref_sidebar <- cbind(ref_sidebar, rownames_to_column(diet_ref_fams)[match(ref_sidebar$family, rownames_to_column(diet_ref_fams)$rowname),])
rownames(ref_sidebar) <- NULL
#Add info about genus mentions in literature
ref_sidebar <- cbind(ref_sidebar, rownames_to_column(diet_ref_gen)[match(ref_sidebar$genus, rownames_to_column(diet_ref_gen)$rowname),])
#A bunch of NAs because of the genera that are not mentioned in the references
#turn them into FALSE (because the species is absent)
ref_sidebar[is.na(ref_sidebar)] <- FALSE

#Remove unnecessary colums
ref_sidebar <- ref_sidebar[,-which(colnames(ref_sidebar) %in% c("family", "rowname")),]

#Replace FALSE -> 0, TRUE -> 1
ref_sidebar[ref_sidebar==FALSE] <- 0

#Add columns referring to the same diet. 
#This way if the exact genus is mentioned in the literature it takes a value of 2
#But if it is only the family, it takes a value of 1

ref_sidebar$western <- ref_sidebar$western + ref_sidebar$western.1
ref_sidebar$grauers <- ref_sidebar$grauers + ref_sidebar$grauers.1
ref_sidebar$mountain <- ref_sidebar$mountain + ref_sidebar$mountain.1

#Remove extra columns
ref_sidebar <- ref_sidebar[, -which(grepl(".1", colnames(ref_sidebar)))]
rownames(ref_sidebar) <- NULL

#Melt
ref_sidebar <- reshape::melt(ref_sidebar)

colnames(ref_sidebar) <- c("diet_gen", "host_subspecies", "presence/absence")

#Set presence/absence as numeric
ref_sidebar$`presence/absence` <- as.numeric(ref_sidebar$`presence/absence`)

#Plot heatmap
heat_sidebar <- heat(ref_sidebar, Yvar="diet_gen", fill="presence/absence", Xvar = "host_subspecies", order.cols = FALSE, order.rows = FALSE) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.title=element_text(size=10, hjust = 0.5)) +
  guides(fill=FALSE) + scale_fill_gradient2(low = "grey", mid="palegreen2", high = "palegreen4", midpoint = 1) +
  geom_hline(yintercept=heat_diet_fam_freqs$Freq + 0.5, size=0.5, colour="white") +
  ggtitle("Mentioned\nin literature")

#Combine the two plots
heat_diet_complete <- plot_grid(heat_diet, heat_sidebar, align = "h", ncol = 2, rel_widths = c(45, 5),
                                axis="tb")
heat_diet_complete

ggsave(heat_diet_complete, file="/proj/sllstore2017021/nobackup/MARKELLA/D3_diet_stats/heat_diet_complete.png", device="png", height=7, width=10)

#Get the total number of reads for differentially abundant taxa
ancom_taxa_rc <- taxa_sums(subset_taxa(euk_genus_diet, taxa_names(euk_genus_diet) %in% ancom_euk$genus))
summary(ancom_taxa_rc)
summary(taxa_sums(euk_genus_diet))


#Create lookup table to compare genera in the dataset with dietary genera
lookup_diet_genera <- ancom_euk[,1:2]

#For each of the differentially abundant dietary taxa
for (i in 1:nrow(lookup_diet_genera)) {
  #Collect the names of the sister taxa found in the diet of each host subspecies
  wes <- unique(remis_western$genus_name[which(remis_western$family_name == as.character(lookup_diet_genera$family[i]))])
  wes <- unique(append(wes,
                rogers_western$genus_name[which(rogers_western$family_name == as.character(lookup_diet_genera$family[i]))]))
  gr <- unique(yamagiwa_grauers$genus_name[which(yamagiwa_grauers$family_name == as.character(lookup_diet_genera$family[i]))])
  mo <- unique(rothman_mountain$genus_name[which(rothman_mountain$family_name == as.character(lookup_diet_genera$family[i]))])
  #Keep unique entries and format as a comma-separated list
  lookup_diet_genera$western_diet[i] <- ifelse(is_empty(wes), NA,
                                               ifelse(length(wes)==0, wes, paste(wes, spe=",")))
  lookup_diet_genera$grauers_diet[i] <- ifelse(is_empty(gr), NA, gr)
  lookup_diet_genera$mountain_diet[i] <- ifelse(is_empty(mo), NA, mo)
}

