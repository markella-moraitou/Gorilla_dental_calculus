#load packages
library(stringr)
library(reshape2)
library(ggpubr)
library(phyloseq)
library(dplyr)
library(tibble)

load("/proj/sllstore2017021/nobackup/MARKELLA/T3_community-level/.RData")

#Technical stats that will be needed for the dataset + tables for the Supplementary.

#### Table S2 - Read counts ####

##Read counts for read preprocessing
rc_prepro <- read.table("/proj/sllstore2017021/nobackup/MARKELLA/readcount_211016.txt", sep=",", comment.char="", header=T)

#Remove _m suffix
colnames(rc_prepro) <- str_remove(colnames(rc_prepro), "_m")

#Removed human host mapping columns (it is included in the next file)
rc_prepro <- rc_prepro[,-which(colnames(rc_prepro) == "hostHumFilt")]

##Read counts after human-host mapping
rc_human_host <- read.table("/proj/sllstore2017021/nobackup/MARKELLA/8_humanHostFilt/unmapped/readcount_hostmapping.txt", sep=",", comment.char="", header=T)

#Fix row names
rc_human_host$sample <- str_remove(rc_human_host$sample, "_m")

#Fix column names
colnames(rc_human_host) <- c("sample", "unmapped", "mapped", "mapped_host", "mapped_human")

##Read counts after decontaminating using Kraken-Tools
rc_decontam_kt <- read.table("/proj/sllstore2017021/nobackup/MARKELLA/RD1_krakentools/readcount_decontam_kt_211024.txt", sep=",", comment.char="", header=T)

colnames(rc_decontam_kt) <- c("sample", "after_kraken_tools")

##Read counts after decontaminating with mapping
rc_decontam_map <- read.table("/proj/sllstore2017021/nobackup/MARKELLA/RD2_mapping/decontam_files/readcount_decontam_211021.txt", sep=",", comment.char="", header=T)

##Merge tables together

readcounts <- merge(rc_prepro, rc_human_host)
readcounts$decontam_kt <- rc_decontam_kt[match(readcounts$sample, rc_decontam_kt$sample),2]
readcounts$decontam_map <- rc_decontam_map[match(readcounts$sample, rc_decontam_map$sample),6]

#Fix sample names
readcounts$sample <- str_remove(readcounts$sample, "./")

#Keep samples that were retained after Kraken2/Bracken
readcounts <- readcounts %>% filter(sample %in% rownames(metadata))

#Melt and plot (exluding details about human host mapping)
readcounts_m <- melt(readcounts[ , -which(colnames(readcounts) %in% 
                                            c("mapped", "mapped_host", "mapped_human"))])
colnames(readcounts_m) <- c("sample", "step", "read_count")

ggsave(
  ggboxplot(readcounts_m, y="read_count", x="step", fill="grey") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)),
  file="readcounts_boxplot.png",
  device="png")

#Add sample type info
readcounts$sample.type <- metadata$Sample.type[match(readcounts$sample, rownames(metadata))]

#Add 'retained' status
readcounts$retained <- 
  #A sample is considered retained if it is there in the final dataset
  ifelse(readcounts$sample.type=="sample", ifelse(readcounts$sample %in% sample_names(spe_data_final), "YES", "NO"),
  #Blanks/controls are not retained in the final dataset, so it is enough that they go through Kraken2/Bracken
  ifelse(readcounts$sample %in% sample_names(spe_data), "YES", "NO"))
#Add suspecies info
readcounts$subspecies <- metadata$Spec.subspecies[match(readcounts$sample, rownames(metadata))]

#Save table
write.table(readcounts, file="readcounts.csv", sep=",", quote=FALSE)

#### Read counts summary ####
#Get summary stats for the read count of samples and blanks separately
readcounts_summary <- list()
for (step in c("adapterRem", "unmapped", "decontam_map")) {
  table <- data.frame()
  for (group in c("sample", "blank/control")) {
    for (statistic in c("mean", "median", "SD", "maximum", "minimum")) {
      #Get a subset of the readcounts table referring to either samples or control/blanks
      expression <- (group=="sample")
      subset <- readcounts %>% filter((sample.type=="sample")==expression) %>% pull(step)
      #Calculate statistics
      table[group, statistic] <- 
        ifelse(statistic=="mean", mean(subset, na.rm=TRUE),
          ifelse(statistic=="median", median(subset, na.rm=TRUE),
            ifelse(statistic=="SD", sd(subset, na.rm=TRUE),
              ifelse(statistic=="maximum", max(subset, na.rm=TRUE),
                ifelse(statistic=="minimum", min(subset, na.rm=TRUE), NULL)))))    
    }
  }
  readcounts_summary[[step]] <- table
}

#### Table S2 - Taxa retained and oral proportion ####
taxa_n_retained <- data.frame(step=c("raw", "decontam", "abundance_filt", 
                                     "env_removal", "reference_based"),
                              species_n=rep(NA, 5), #Total number of species
                              genera_n=rep(NA, 5), #Number of genera
                              hominid_core_n=rep(NA, 5), #Number of core hominid microbiome taxa
                              homd_n=rep(NA, 5), #Number of HOMD taxa
                              contaminant_n=rep(NA, 5), #Number of known contaminants
                              
                              hominid_core_prop=rep(NA, 5), #Proportion of core hominid microbiome taxa
                              homd_prop=rep(NA, 5), #Proportion of HOMD taxa
                              contaminant_prop=rep(NA, 5), #Proportion of known contaminants
                              
                              hominid_core_rel_abund=rep(NA, 5), #Rel. abundance of core hominid microbiome taxa
                              homd_rel_abund=rep(NA, 5)) #Rel. abundance of HOMD taxa
i=1
for (n in c(spe_data, spe_data_decontam, spe_data_filt, spe_data_envrem, spe_data_final)) {
  speciesN <- ntaxa(n)
  generaN <- n@tax_table[,7] %>% unique %>% length 
  core_n <- length(intersect(tax_table(n)[,8], core_micr$Taxon))
  homd_n <- length(intersect(taxa_names(n), homd$NCBI_taxon_id))
  contaminant_n <- length(which(tax_table(n)[,7] %in% contaminant_list))
  total_abundance <- taxa_sums(n) %>% sum
  
  taxa_n_retained$species_n[i] <- speciesN
  taxa_n_retained$genera_n[i] <- generaN
  taxa_n_retained$hominid_core_n[i] <- core_n
  taxa_n_retained$homd_n[i] <- homd_n
  taxa_n_retained$contaminant_n[i] <- contaminant_n
  
  taxa_n_retained$hominid_core_prop[i] <- core_n/speciesN
  taxa_n_retained$homd_prop[i] <- homd_n/speciesN
  taxa_n_retained$contaminant_prop[i] <- contaminant_n / speciesN
  
  taxa_n_retained$hominid_core_rel_abund[i] <- (subset_taxa(n, species %in% core_micr$Taxon) %>% taxa_sums %>% sum ) / total_abundance
  taxa_n_retained$homd_rel_abund[i] <- (subset_taxa(n, taxa_names(n) %in% homd$NCBI_taxon_id) %>% taxa_sums %>% sum ) / total_abundance
  taxa_n_retained$contaminant_rel_abund[i] <- (subset_taxa(n, genus %in% contaminant_list) %>% taxa_sums %>% sum ) / total_abundance
    
  i<-i+1
}

#Taxa that are genus-level oral
#spe_data_final %>% subset_taxa(genus %in% homd$Genus | genus %in% str_remove(core_micr$Taxon, " .*")) %>% taxa_n

#Save table
write.table(taxa_n_retained, file="taxa_n_retained.txt", sep=",", quote=FALSE, row.names=FALSE, dec=".")

### Not in Supp - Taxa retained per sample in raw and final dataset
taxa_per_sample <- data.frame(row.names=c("readcount_before_Kraken", "taxa_total_raw", "taxa_total_final", "core_micr_raw", "core_micr_final", "homd_raw", "homd_final", "core_micr_prop_raw", "core_micr_prop_final", "homd_prop_raw", "homd_prop_final"))

for (i in 1:length(sample_names(spe_data_final))) {
  sample=sample_names(spe_data_final)[i]
  #Number of reads
  readcount_before_Kraken <- metadata %>% filter(rownames(metadata) == sample) %>% pull(readcount.m.before.Kraken)
  #Get list of taxa present in sample in the raw dataset
  taxa_list_raw <- prune_samples(sample_names(spe_data)==sample, spe_data)
  taxa_list_raw <- prune_taxa(taxa_sums(taxa_list_raw)>0, taxa_list_raw) %>% tax_table %>% as.data.frame %>% select(species)
  #Get list of taxa present in sample in the final dataset
  taxa_list_final <- prune_samples(sample_names(spe_data_final)==sample, spe_data_final)
  taxa_list_final <- prune_taxa(taxa_sums(taxa_list_final)>0, taxa_list_final) %>% tax_table %>% as.data.frame %>% select(species)
  #Total number of taxa
  taxa_total_raw <- nrow(taxa_list_raw)
  taxa_total_final <- nrow(taxa_list_final)
  #Number of core micr taxa
  core_micr_raw <- which(taxa_list_raw$species %in% core_micr$Taxon) %>% length
  core_micr_final <- which(taxa_list_final$species %in% core_micr$Taxon) %>% length
  #Number of HOMD taxa
  homd_raw <- which(rownames(taxa_list_raw) %in% homd$NCBI_taxon_id) %>% length
  homd_final <- which(rownames(taxa_list_final) %in% homd$NCBI_taxon_id) %>% length  
  #Proportion of core micr taxa
  core_micr_prop_raw <- core_micr_raw / taxa_total_raw
  core_micr_prop_final <- core_micr_final / taxa_total_final
  #Proportion of HOMD taxa
  homd_prop_raw <- homd_raw / taxa_total_raw
  homd_prop_final <- homd_final / taxa_total_final
  #Append row
  add <- c(readcount_before_Kraken, taxa_total_raw, taxa_total_final, core_micr_raw, core_micr_final, homd_raw, homd_final, core_micr_prop_raw, core_micr_prop_final, homd_prop_raw, homd_prop_final)
  taxa_per_sample <- cbind(taxa_per_sample, add)
  colnames(taxa_per_sample)[i] <- sample
 }
  
#transpose
taxa_per_sample <- t(taxa_per_sample)

#Save table
write.table(taxa_per_sample, file="taxa_per_sample.txt", sep=",", quote=FALSE, dec=".")


#### Table S5 - Subspecies verification ####
#Gather all data on mitochondrial genomes

#Get mitochondrial genome coverage per sample
mt_coverage <- read.table("/proj/sllstore2017021/nobackup/MARKELLA/H1_mtgenomes/mt_coverage.csv", header=TRUE, sep=",", dec=".")

#Add diagnostic sites per sample
diagnostic_sites <- read.table("/proj/sllstore2017021/nobackup/MARKELLA/H3_diagnostic_sites/output/distance_to_refs_east_west.txt", header=TRUE, sep="\t", row.names="Sample") %>% cbind(read.table("/proj/sllstore2017021/nobackup/MARKELLA/H3_diagnostic_sites/output/distance_to_refs_mountain_grauers.txt", header=TRUE, sep="\t", row.names="Sample"))

#Fix colnames
colnames(diagnostic_sites) <- make.unique(colnames(diagnostic_sites))

#Make sample IDs a column
diagnostic_sites <- rownames_to_column(diagnostic_sites, var="sample")

#Combine the two tables
mt_info <- left_join(mt_coverage, diagnostic_sites)

#Get subspecies classification based on presence of diagnostic sites

#Function to get the assigned species and associated p-value using chi-squared test


mt_info <- mt_info %>% mutate(inferred_subspecies = ifelse(Western > Eastern, "gorilla", 
                                                            ifelse(Grauers > Mountain, "graueri",
                                                                    ifelse(Mountain > Grauers, "beringei", "undetermined"))),
                         #Get subspecies according to metadata
                         metadata_subspecies = metadata$Spec.subspecies[match(mt_info$sample, rownames(metadata))])

write.table(mt_info, "mt_genomes_info.txt", sep=",", quote=FALSE, dec=".")


#### Table S6 - ANCOM taxa ####

#Show all differentially abundant species by host and facility

ancom_subspecies <-
 read.table("/proj/sllstore2017021/nobackup/MARKELLA/T3_community-level/ancom_results.csv", header = TRUE, sep = ",", dec = ".") %>% filter(detected_0.9==TRUE) %>% select(taxa_id, W, structural_zero..gorilla., structural_zero..graueri., structural_zero..beringei.) %>% mutate(taxa_names=taxonomy_species[match(taxa_id, rownames(taxonomy_species)),8]) %>%
 mutate(diff_abund_by_subspecies="YES")
 
ancom_facility <-
 read.table("/proj/sllstore2017021/nobackup/MARKELLA/T3_community-level/ancom_results_seqc.csv", header = TRUE, sep = ",", dec = ".") %>% filter(detected_0.9==TRUE) %>% select(taxa_id, W, structural_zero..Jena., structural_zero..Uppsala.) %>%
 mutate(diff_abund_by_facility="YES")
 
#Merge tables
ancom_full <- full_join(ancom_subspecies, ancom_facility, by="taxa_id")

#Some modifications
ancom_full[,which(grepl("structural.zero", colnames(ancom_full)))] <- sapply(ancom_full[,which(grepl("structural.zero", colnames(ancom_full)))], as.logical)
ancom_full$diff_abund_by_subspecies[which(is.na(ancom_full$diff_abund_by_subspecies))] <- "NO"
ancom_full$diff_abund_by_facility[which(is.na(ancom_full$diff_abund_by_facility))] <- "NO"

#Reorder columns
ancom_full <- ancom_full[, c(1, 6, 2, 3, 4, 5, 7, 8, 9, 10, 11)]

write.table(ancom_full, file="ancom_full.csv", sep=",", quote=FALSE)

#### Table S7 - Subspecies associated BP processes ####

load("/proj/sllstore2017021/nobackup/MARKELLA/F2_functional_stats/.RData")
write.table(most_signif_BP, file="most_signif_BPs.csv", sep=",", quote=FALSE)

#### Table S8 - MAG taxonomy ####
load("/proj/sllstore2017021/nobackup/MARKELLA/M6_MAGstats/.RData")

hq_mq_table <- hq_mq_mags %>% select(user_genome, phylum, family, genus, species, host_subspecies, draft_quality, fastani_reference)

write.table(hq_mq_table, file="hq_mq_table.csv", sep=",", quote=FALSE)

#### Table S9 - Dietary ANCOM taxa ####

#Show all differentially abundant species by host and facility

ancom_diet <-
 read.table("/proj/sllstore2017021/nobackup/MARKELLA/D3_diet_stats/euk_ancom_results.csv", header = TRUE, sep = ",", dec = ".") %>% filter(detected_0.9==TRUE) %>% select(taxa_id, W, structural_zero..gorilla., structural_zero..graueri., structural_zero..beringei.)

ancom_diet[,which(grepl("structural.zero", colnames(ancom_diet)))] <- sapply(ancom_diet[,which(grepl("structural.zero", colnames(ancom_diet)))], as.logical)
 
write.table(ancom_diet, file="ancom_diet.csv", sep=",", quote=FALSE)


#### MAG assembly ####

## contig assembly
mag_assemblies <- read.table("/proj/sllstore2017021/nobackup/MARKELLA/M1_MAGassemblies/log_mag_assembly_210719.txt", sep=",", comment.char="", header=T)
dim(mag_assemblies)

table(mag_assemblies$contig_n<2)

#remove samples with 0 or 1 contig
mag_assemblies <- mag_assemblies[which(mag_assemblies$contig_n>2),]

#Get some stats
summary(mag_assemblies$contig_n)
summary(mag_assemblies$min_size)
summary(mag_assemblies$max_size)
summary(mag_assemblies$contig1.5kb_n)
summary(mag_assemblies$contig20kb_n)

## binning
bins <- read.table("/proj/sllstore2017021/nobackup/MARKELLA/M3_binning/CheckM_results/all_checkm_results.tab", sep="\t", comment.char="", header=T)

dim(bins)
table(substring(bins$Bin.Id, 1, 6))
length(unique(substring(bins$Bin.Id, 1, 6)))

#get some stats
summary(bins$Completeness)
summary(bins$Contamination)

save.image()
