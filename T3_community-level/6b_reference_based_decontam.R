#Start logging
sink(file = "log6b_reference_based_decontam.txt")

#packages
library(dplyr)
library(tidyverse)
library(phyloseq)
library(ggpubr)
library(reshape2)
library(readr)
library(phyloseq)
library(tibble)
library(tidytree)
library(ggtree)
library(FEAST)
library(MutationalPatterns)

load(".RData")

#Community-level taxonomic analysis - Script 6b (after running map damage)

#### Investigate misincorporation in ambiguous taxa ####

#Load table
edge_misinc <- read_tsv(file="/proj/sllstore2017021/nobackup/MARKELLA/RD3_mapdamage4ambiguoustaxa/edge_misincorporation.txt")

#TaxID and position shouldn't be treated as numericals
edge_misinc$TaxID <- as.character(edge_misinc$TaxID)
edge_misinc$Pos <- as.factor(edge_misinc$Pos)

#Get average of all base changes per taxon and end
edge_misinc <- 
  edge_misinc %>% group_by(TaxID, End, Pos) %>% summarise_if(is.numeric, mean) %>% ungroup


#Calculate relative frequencies of each base change (change freq / reference base freq)
damage_assessment <- 
  edge_misinc %>% mutate_if(grepl("G>", names(.)), ~ ./G) %>% 
  mutate_if(grepl("C>", names(.)), ~ ./C) %>%
  mutate_if(grepl("T>", names(.)), ~ ./T) %>%
  mutate_if(grepl("A>", names(.)), ~ ./A) %>%
  #Get maximum frequency across all changes in a position
  group_by(TaxID, Pos, End) %>%
  mutate(max_freq=max(`G>A`, `C>T`, `A>G`, `T>C`, `A>C`, `A>T`, `C>G`, `C>A`, `T>G`, `T>A`, `G>C`, `G>T`)) %>% 
  #In 5p: are G>A changes the most frequent and is this frequency above 0.02? 
  #In 3p: are C>T changes the most frequent and this frequency above 0.02?
  mutate(target_change_highest=ifelse(End=="3p", as.logical(`G>A` >= max_freq), as.logical(`C>T` >= max_freq)),
         frequency_above_threshold=ifelse(End=="3p", as.logical(`G>A` >= 0.025), as.logical(`C>T` >= 0.025))) %>%
  #For how many out of 3 positions are the above criteria true?
  ungroup %>% group_by(TaxID, End) %>% 
  summarise(n_damaged_positions = sum(target_change_highest), n_highfreq_positions = sum(frequency_above_threshold)) %>% 
  #If the target base change is the highest in at least two positions and shows frequency above 0.02 in at least one
  #then that end shows damage
  ungroup %>% mutate(end_shows_damage = as.logical(n_damaged_positions >= 2 & n_highfreq_positions >= 1)) %>%
  #If both ends show damage, then flag the taxon as showing damage patterns
  group_by(TaxID) %>% summarise(damage = all(end_shows_damage==TRUE))

print("How many have sufficient damage?")
damage_assessment$damage %>% table

#Are the modern looking or the ancient looking taxa more abundant?
#Get abundances per sample for all ambiguous taxa
damage_assessment %>% left_join(otu_table(spe_data_envrem) %>% as.data.frame %>%
                                  mutate(TaxID=taxa_names(spe_data_envrem)), copy=TRUE, by="TaxID") %>%
  #Get mean abundances per taxon
  melt %>% group_by(TaxID, damage) %>% summarise(mean_abund=mean(value)) %>% ggboxplot(x="damage", y="mean_abund") %>%
  ggsave(device="png", file="/proj/sllstore2017021/nobackup/MARKELLA/T3_community-level/ambiguous_taxa_abundance.png")

#Ancient looking taxa tend to be on average more abundant. 
#However there are several highly abundant modern taxa

#### Final dataset ####

#remove unambiguous contaminants
spe_data_final <- subset_taxa(spe_data_envrem, 
                              !(taxa_names(spe_data_envrem) %in% unambiguous_contam))

#remove ambiguous contaminants without damage patterns
spe_data_final <- subset_taxa(spe_data_final, 
                              !(taxa_names(spe_data_final) %in% 
                                  (damage_assessment %>% filter(damage==FALSE) %>% pull(TaxID))))

#Make host subspecies a factor
spe_data_final@sam_data$Spec.subspecies <- factor(spe_data_final@sam_data$Spec.subspecies,
                                                  levels=c("gorilla", "graueri", "beringei"))
                                                  
#PCoA before removing the duplicates
jaccard_5_with_duplicates <- ordinate(spe_data_final, method="PCoA", distance="jaccard")

spe_data_final_norm <- microbiome::transform(spe_data_final, "clr")
clr_5_with_duplicates <- ordinate(spe_data_final_norm, method="PCoA", distance="euclidean")
                                                  
#Plot
jaccard_with_duplicates <- plot_ordination(spe_data_final, jaccard_5_with_duplicates, color=c("Seq.centre"), shape="Duplic.pair", title="PCoA on jaccard distances based on species level assignments,\nhighlighting samples sequenced in both facilities") + #to highlight samples sequences both in Jena and Uppsala
  scale_shape_manual(values=c(0, 1, 5, 16))
  
clr_with_duplicates <- plot_ordination(spe_data_final_norm, clr_5_with_duplicates, color=c("Seq.centre"), shape="Duplic.pair", title="PCoA on jaccard distances based on species level assignments,\nhighlighting samples sequenced in both facilities") + #to highlight samples sequences both in Jena and Uppsala
  scale_shape_manual(values=c(0, 1, 5, 16))

#Remove duplicates
spe_data_final <- subset_samples(spe_data_final, 
                                 !(sample_names(spe_data_final) %in% c("IBA001", "WAL001", "BIT001")))

#### Rerunning FEAST after all decontamination and filtering steps ####
#This isn't used to remove bad samples, but to show how the oral proportion has increased

#Create a new table including the filtered OTU table plus the sources (NOTE: this table includes only taxa that were retained after filtering)
feast_table_filt <- abundance_filter_f(as.data.frame(feast_table), 0.0001)

#Only keep sources - samples will be substituted
feast_table_filt <- feast_table_filt[which(grepl("soil|human", rownames(feast_table_filt))),]

#Add samples from the final OTU table
feast_table_filt <-
  feast_table_filt %>% t %>% as.data.frame %>%
  merge(as.data.frame(otu_table(spe_data_final)), by=0, all=TRUE) %>% 
  column_to_rownames(var="Row.names") %>% t

#Turn NA to 0
feast_table_filt[is.na(feast_table_filt)] <- 0

#Keep the metadata rows that are needed
feast_metadata_filt <- feast_metadata[which(rownames(feast_metadata) %in% rownames(feast_table_filt)),]
id = 1
for (i in 1:nrow(feast_metadata_filt)){ #change the id's so that they are consistent
  if (!is.na(feast_metadata_filt$id[i])){
    feast_metadata_filt$id[i] <- id
    id <- id + 1
  }
}  

#Rerun FEAST for filtered dataset
#feast_output_filt <- FEAST(feast_table_filt, feast_metadata_filt, different_sources_flag = FALSE,
#                           outfile = "species_dc_filt", dir_path = "/proj/sllstore2017021/nobackup/MARKELLA/T3_community-level")
#Already produced so I am reading off the saved file
feast_output_filt <- read.table("/proj/sllstore2017021/nobackup/MARKELLA/T3_community-level/species_dc_filt_source_contributions_matrix.txt", header = TRUE, sep = "\t", dec = ".")
                           
#Construct an edited matrix that shows the contributions per environment not per source
feast_output_env_filt <- matrix(nrow = nrow(feast_output_filt), ncol = 7)
rownames(feast_output_env_filt) <-  gsub("_[a-z]+", "", rownames(feast_output_filt))
colnames(feast_output_env_filt) <- c("human_calculus", "human_gut", "human_plaque", "human_skin", "labcontam_rmhuman", "soil_tundra", "Unknown")

for (i in 1:nrow(feast_output_env_filt)){
  for (j in 1:ncol(feast_output_env_filt)){
    #create a subset of the data that only contain the columns of the environment indicated by colnames(feast_output_env)[j]
    pattern <- row.names(feast_metadata_filt)[which(feast_metadata_filt$Env==colnames(feast_output_env_filt)[j])] #the sources from the environment of interest
    x = colnames(feast_output_filt) #all the sources
    pattern <- gsub("_[a-z]+", "", pattern) #Keep only IDs
    x <- gsub("_[a-z]+", "", x) #Keep only IDs
    subset <- feast_output_filt[,which(x %in% pattern)] #pick the data about the sources of the environment we are interested in
    feast_output_env_filt[i,j]=rowSums(subset)[i] #sum the contributions from the subset to take the value for the contribution of the environment
    if (colnames(feast_output_env_filt)[j]=="Unknown"){
      feast_output_env_filt[i,j]=feast_output_filt$Unknown[i]
    }
  }
}

#Row sums should be equal to 1
rowSums(feast_output_env_filt)

#Reorder rows
feast_output_env_filt <- feast_output_env_filt[,c(1,3,2,4,5,6,7)]

#Also plot contributions before decontam and filtering, but only for retained samples (for comparison purposes)
oral_proportion_filt <- feast_output_env_filt %>% as.data.frame %>% select(human_calculus, human_plaque) %>%
 mutate(oral_proportion_log = log(human_calculus + human_plaque + 0.01))
 
ggsave(gghistogram(oral_proportion_filt, "oral_proportion_log", fill="#C77CFF", position="stack") +
        #Include vertical line with the cutoff of 0.03 (3%)
        geom_vline(xintercept=log(0.03)),
 file="/proj/sllstore2017021/nobackup/MARKELLA/T3_community-level/oral_proporton_hist_filt.png")

#### Taxa lists for read-level decontam ####
#Extract a list of "exogenous" taxa -- based on this list, the taxa that have been removed in this analysis
#will also be removed in the read level using Kraken-tools
#This includes taxa marked as contaminants by decontam, abundance-based decontamination and reference-based decontamination

exogenous_id_list <- append(rownames(contam), env_taxa$TaxID) %>% 
  append(damage_assessment %>% filter(damage==FALSE) %>% pull(TaxID))

write.table(exogenous_id_list, file="/proj/sllstore2017021/nobackup/MARKELLA/RD2_mapping/exogenous_id_list.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

print("How many taxa are there in the exogenous taxa ID list")
length(exogenous_id_list)            

#do the same with retained samples
retained_samples <- sample_names(spe_data_final)
print("How many samples are retained?")
length(retained_samples)


write.table(retained_samples, file="/proj/sllstore2017021/nobackup/MARKELLA/RD2_mapping/retained_samples.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

#I will also extract a list of the taxa with the highest read count in the envrem dataset, 
#to use for comparative mapping

#Keep only samples from the retained samples list
abundant_id_list <-
  spe_kraken %>% melt %>% filter(variable %in% retained_samples) %>%
  #Get all taxa that have number of reads above 1.000
  filter(value>1000) %>% pull(X.OTU.ID) %>% unique

#Get the taxa that haven't been listed as exogenous
abundant_id_list <- setdiff(abundant_id_list, exogenous_id_list)

#Keep only the taxa that exist in the final dataset
abundant_id_list <- intersect(abundant_id_list, taxa_names(spe_data_final))

print("How many taxa are there in the abundant taxa ID list")
length(abundant_id_list)


#Plot contributions of retained samples after filtering
png(file = "/proj/sllstore2017021/nobackup/MARKELLA/T3_community-level/source_contribution_filt.png", width = 1000, height = 480)
feast_output_env_filt[which(rownames(feast_output_env_filt) %in% retained_samples),] %>% t %>%
  plot_contribution() + 
  theme(axis.text.x=element_text(angle = +90, hjust = 0)) +
  scale_fill_brewer(palette = "Spectral")
dev.off()

#Plot contributions before filtering and decontam(only of retained samples for comparison purposes)
png(file = "/proj/sllstore2017021/nobackup/MARKELLA/T3_community-level/source_contribution_retained_samples.png", width = 1000, height = 480)
#Filter out bad samples, blanks and controls
feast_output_env[which(rownames(feast_output_env) %in% retained_samples),] %>% t %>%
  plot_contribution() + 
  theme(axis.text.x=element_text(angle = +90, hjust = 0)) +
  scale_fill_brewer(palette = "Spectral")
dev.off()

#Normalize
spe_data_final_norm <- microbiome::transform(spe_data_final, "clr")

#
print("How does the abundance compare?")
print("Abundant taxa sums")
spe_data_envrem %>% subset_taxa(taxa_names(spe_data_envrem) %in% abundant_id_list) %>%
  taxa_sums() %>% summary
print("Exogenous taxa sums")
spe_data_envrem %>% subset_taxa(taxa_names(spe_data_envrem) %in% exogenous_id_list) %>%
  taxa_sums() %>% summary
#The exogenous taxa are less abundant


#
print("There should be no intersection with exogenous")
length(intersect(abundant_id_list, exogenous_id_list))

#Save file
write.table(abundant_id_list, file="/proj/sllstore2017021/nobackup/MARKELLA/RD2_mapping/abundant_id_list.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

#### Plot ####
#PCoA using jaccard for species table
jaccard_5 <- ordinate(spe_data_final, method="PCoA", distance="jaccard")

#Plot ordination
pdf(file = "/proj/sllstore2017021/nobackup/MARKELLA/T3_community-level/jaccard_5.pdf")
par(mfrow=c(1,3))
plot_ordination(spe_data_final, jaccard_5, color=c("Seq.centre"), shape="Spec.subspecies", title="PCoA on jaccard distances based on species level assignments \nafter reference-based decontamination") #to highlight seq. centre
plot_ordination(spe_data_final, jaccard_5, color=c("Spec.subspecies"), title="PCoA on jaccard distances based on species level assignments \nafter reference-based decontamination") #to highlight subspecies
dev.off()

### PCoA using Aitchison distances ###

#PCoA using euclidean distances of CLR-normalized abundances on the species table 
clr_5 <- ordinate(spe_data_final_norm, method="PCoA", distance="euclidean")

#Plot ordination
pdf(file = "/proj/sllstore2017021/nobackup/MARKELLA/T3_community-level/clr_5.pdf")
par(mfrow=c(1,3))
plot_ordination(spe_data_final_norm, clr_5, color=c("Seq.centre"), shape="Spec.subspecies", title="PCoA on Aitchison distances based on species level assignments \nafter reference-based decontamination") #to highlight seq. centre
plot_ordination(spe_data_final_norm, clr_5, color=c("Spec.subspecies"), title="PCoA on Aitchison distances based on species level assignments \nafter reference-based decontamination") #to highlight subspecies
dev.off()

sessionInfo()
#Stop logging
sink(file = NULL)

save.image()
