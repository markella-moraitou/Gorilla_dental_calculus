#Start logging
sink(file = "log2_FEAST.txt")

#load packages
Packages <- c("Rcpp", "RcppArmadillo", "vegan", "dplyr", "reshape2", "gridExtra", "ggplot2", "ggthemes")
lapply(Packages, library, character.only = TRUE)
library(FEAST)
library(stringr)
library(MutationalPatterns)
library(plyr)
library(ggpubr)
library(stringr)

load(".RData")

#Community-level taxonomic analysis - Script 2
#Identifying sources contributing to the dataset using FEAST

#### Prepare files for running FEAST ####

#Import kraken-biom of samples (sinks) and sources
feast_table <- read.table("T2_bracken/Sources/species_table_feast.txt", skip=1, sep="\t", comment.char="", header=T)

#Modifications on otu table
feast_table <- feast_table[,which(colSums(feast_table)!=0)] #remove columns with sum 0
feast_table <- t(feast_table) #transpose table. Vegdist documentation refers to columns as species and rows as sites.
colnames(feast_table)<-as.character(feast_table[1,]) #set first row as column (taxa) names
feast_table <- feast_table[-1,] #remove original first row
feast_table <- as.matrix(feast_table) #turn into matrix

#for-loop to keep only the base name for sample names
for(i in 1:length(row.names(feast_table))){
  j=str_remove(row.names(feast_table)[i], "_kraken2_report_bracken_species")
  row.names(feast_table)[i]<-j
}

#also to remove _m from sample names
for(i in 1:length(row.names(feast_table))){
  j=str_remove(row.names(feast_table)[i], "_m")
  row.names(feast_table)[i]<-j
}

#Remove samples that were removed from the phyloseq objects because they were low in content
feast_table <- feast_table[which(!(rownames(feast_table) %in% low_content_samples)),]

#create metadata matrix for FEAST
feast_metadata <- matrix( nrow = length(row.names(feast_table)), ncol = 3) #create empty matrix
colnames(feast_metadata) <- c("Env", "SourceSink", "id")
row.names(feast_metadata) <- row.names(feast_table)
feast_metadata <- as.data.frame(feast_metadata) #turn into dataframe

#populate feast_metadata dataframe
id=1 #will increment this by one for every sink
for (i in 1:length(row.names(feast_table))){
  if (row.names(feast_table)[i] %in% colnames(species_table)){
    feast_metadata$SourceSink[i] <- "Sink"
    feast_metadata$id[i] <- id #give a unique number from 1 to [number of sinks]
    id <- id + 1
    if (grepl("^B[LE]", row.names(feast_table)[i])){
      feast_metadata$Env[i] = "lab" 
    } else if (grepl("^ERR", row.names(feast_table)[i]) | grepl("^BS", row.names(feast_table)[i])){
      feast_metadata$Env[i] = "museum"
    } else {
      feast_metadata$Env[i] = "gorilla_calculus"
    }
  } else {
    feast_metadata$SourceSink[i] <- "Source"
    feast_metadata$Env[i] <- gsub("[SE]R[SR]\\d+_", "", row.names(feast_table)[i]) #use environment from name
  }
}

#### Run feast and investigate output ####
#feast_output <- FEAST(feast_table, feast_metadata, different_sources_flag = FALSE, outfile = "species_dc", dir_path = "T3_community-level")
#Already produced so I am reading off the saved file
feast_output <- read.table("T3_community-level/species_dc_source_contributions_matrix.txt", header = TRUE, sep = "\t", dec = ".")

#Construct an edited matrix that shows the contributions per environment not per source
feast_output_env <- matrix(nrow = nrow(feast_output), ncol = 7)
rownames(feast_output_env) <-  gsub("_[a-z]+", "", rownames(feast_output))
colnames(feast_output_env) <- c("human_calculus", "human_gut", "human_plaque", "human_skin", "labcontam_rmhuman", "soil_tundra", "Unknown")

for (i in 1:nrow(feast_output_env)){
  for (j in 1:ncol(feast_output_env)){
    #create a subset of the data that only contain the columns of the environment indicated by colnames(feast_output_env)[j]
    pattern <- row.names(feast_metadata)[which(feast_metadata$Env==colnames(feast_output_env)[j])] #the sources from the environment of interest
    x = colnames(feast_output) #all the sources
    pattern <- gsub("_[a-z]+", "", pattern) #Keep only IDs
    x <- gsub("_[a-z]+", "", x) #Keep only IDs
    subset <- feast_output[,which(x %in% pattern)] #pick the data about the sources of the environment we are interested in
    feast_output_env[i,j]=rowSums(subset)[i] #sum the contributions from the subset to take the value for the contribution of the environment
    if (colnames(feast_output_env)[j]=="Unknown"){
      feast_output_env[i,j]=feast_output$Unknown[i]
    }
  }
}

#Row sums should be equal to 1
rowSums(feast_output_env)

#Reorder rows
feast_output_env <- feast_output_env[,c(1,3,2,4,5,6,7)]

#Plot contributions
png(file = "T3_community-level/source_contribution.png", width = 1000, height = 480)
plot_contribution(t(feast_output_env)) +
  theme(axis.text.x=element_text(angle = +90, hjust = 0)) +
  scale_fill_brewer(palette = "Spectral")
dev.off()

#Deciding which samples to exclude due to low oral content

#Plot distribution of oral source contributions (summing human_calculus and human plaque)
oral_proportion <- feast_output_env %>% as.data.frame %>% select(human_calculus, human_plaque) %>%
 mutate(oral_proportion_log = log(human_calculus + human_plaque + 0.01)) %>%
 mutate(Sample.type=metadata$Sample.type[match(rownames(.), rownames(metadata))])
 
ggsave(gghistogram(oral_proportion, "oral_proportion_log", fill="Sample.type", position="stack") +
        #Include vertical line with the cutoff of 0.03 (3%)
        geom_vline(xintercept=log(0.03), size=1, colour="red"),
 file="T3_community-level/oral_proporton_hist.png")
 
#### Remove bad samples ####

#bad samples = samples below the 0.03 oral source threshold plus a dead in captivity sample
bad_samples <- oral_proportion %>% filter(oral_proportion_log < log(0.03), Sample.type=="sample") %>% rownames
bad_samples <- append(bad_samples, "G0006") #Add G0006 (Dead in captivity) to samples to be removed

print("Samples to be removed:")
bad_samples

#Bad samples will be removed after running decontam (3_decontam.R)

#Plot contributions only for retained samples (for comparison purposes)
png(file = "T3_community-level/source_contribution_retained_samples.png", width = 1000, height = 480)
#Filter out bad samples, blanks and controls
feast_output_env[which(!(rownames(feast_output_env) %in% bad_samples | grepl("BE|BL|BS|ER|EX|LI", rownames(feast_output_env)))),] %>% t %>%
  plot_contribution() + 
  theme(axis.text.x=element_text(angle = +90, hjust = 0)) +
  scale_fill_brewer(palette = "Spectral")
dev.off()

sessionInfo()
#Stop logging
sink(file = NULL)

save.image()
