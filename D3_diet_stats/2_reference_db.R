#Start logging
sink(file = "log2_reference_db.txt")


#load packages
library(phyloseq)
library(taxize)
library(reshape2)
library(microbiome)
library(splitstackshape)
library(stringr)

load(".RData")

#Diet analysis - Script 3
#Create a reference database based on the literature

####Load tables from previous publications ####

#Script for processing the databases

#Yamagiwa-2005
yamagiwa_grauers <- read.table("/proj/sllstore2017021/nobackup/MARKELLA/D3_diet_stats/Yamagiwa-2005.txt", sep="\t", fill = TRUE)
yamagiwa_grauers <- as.data.frame(unique(yamagiwa_grauers$V1))
colnames(yamagiwa_grauers) <- "taxon_name"

#Rothman-2014
rothman_mountain <- read.table("/proj/sllstore2017021/nobackup/MARKELLA/D3_diet_stats/Rothman-2014.txt", sep="\t", fill = TRUE)
colnames(rothman_mountain) <- "taxon_name"
rothman_mountain <- unique(rothman_mountain)

#Remis-2001
remis_western <- read.table("/proj/sllstore2017021/nobackup/MARKELLA/D3_diet_stats/Remis-2001.txt", sep="\t", fill = TRUE)
colnames(remis_western) <- "taxon_name"
remis_western <- unique(remis_western)

#Rogers-2004
rogers_western <- read.table("/proj/sllstore2017021/nobackup/MARKELLA/D3_diet_stats/Rogers-2004.txt", sep="\t", fill = TRUE)
colnames(rogers_western) <- "taxon_name"
rogers_western <- unique(rogers_western)

#Michel-in-prep
michel_grauers <- read.table("/proj/sllstore2017021/nobackup/MARKELLA/D3_diet_stats/Michel-in-prep.txt", sep="\t", fill = TRUE)
michel_grauers <- as.data.frame(michel_grauers[!is.na(michel_grauers)])
colnames(michel_grauers) <- "taxon_name"
michel_grauers <- unique(michel_grauers)



#Function for processing the references - 1st step: looking for taxIDs
diet_db_1 <- function(taxon_list) {
  #Get taxonomic ID were possible
  for(i in 1:nrow(taxon_list)) {
    taxon_list$TaxID[i] <- get_uid(taxon_list$taxon_name[i], db="ncbi")[[1]]
  }
  #Add column for updated names
  taxon_list$updated_name <- NA
  #Add column for GBIF tax IDs
  taxon_list$GBIF_taxon_ID <- NA
  return(taxon_list)
}

#Define a list with a references to be used
diet_ref_list <- c("rogers_western", "remis_western", "yamagiwa_grauers", "michel_grauers", "rothman_mountain")

#Run the first processing step for each reference and write out file
#for (i in diet_ref_list){
#  assign(i, diet_db_1(get(i)))
#  write.table(get(i)$taxon_name[which(is.na(get(i)$TaxID))],
#              file=paste("/proj/sllstore2017021/nobackup/MARKELLA/D3_diet_stats/",
#                         str_to_title(str_remove(i, "_.*")), "-ood_names.txt", sep = ""), 
#              sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
#  #Print table of missing tax IDs
#  print(table(is.na(get(i)$TaxID)))
#}
#This script prompts interaction when a genus name appears in multiple clades. The option "eudicots" or "monocots" was always selected

#Edits to be done: 
#Replace XXX sp. to XXX (reason: unspecified_species)
#manually check NCBI Taxonomy database to see if the name is updated but misspelled (reason: misspelled)
#Check for copying errors (reason: miscopied)
#Look for updated names at https://www.gbif.org/species/search?q=
#If taxon not found - reason: absent_gbif
#If taxon present in GBIF but absent from NCBI - reason: absent_ncbi and provide name of specific entry in GBIF and GBIF taxon ID
#If taxon found in both - reason: synonym

#### Update database ####

#Function for processing the references - 2nd step: loading updated names and GBIF TaxIDs when a taxon was absent from NCBI
diet_db_2 <- function(taxon_list, taxon_list_name){
  taxon_list_updates <- read.table(paste("/proj/sllstore2017021/nobackup/MARKELLA/D3_diet_stats/",
                                         str_to_title(str_remove(taxon_list_name, "_.*")), "-ood_names_updates.txt", sep = ""), 
                                   sep="\t", header = TRUE, fill = TRUE)
  ##Incorporate updated names and tax IDs in the main files
  for (j in 1:nrow(taxon_list_updates)){
    #For unspecified_species, misspelled, miscopied and synonym -> add updated name
    orig_name <- taxon_list_updates$original_name[j]
    if (taxon_list_updates$update_reason[j] %in% c("misspelled", "synonym", "unspecified_species")) {
      taxon_list$updated_name[which(taxon_list$taxon_name==orig_name)] <- taxon_list_updates$update_names[j] 
      #For absent_ncbi -> add GBIF taxon ID  
    } else if (taxon_list_updates$update_reason[j] == "absent_ncbi") {
      taxon_list$GBIF_taxon_ID[which(taxon_list$taxon_name==orig_name)] <- taxon_list_updates$GBIF_ID[j]}}
  #Get missing NCBI TaxIDs
  for(j in 1:nrow(taxon_list)) {
    #If the taxID is missing and the name has been updated, look again for a taxID using the new name
    if (!is.na(taxon_list$updated_name[j]) & is.na(taxon_list$TaxID[j])) {
      taxon_list$TaxID[j] <- get_uid(taxon_list$updated_name[j], db="ncbi")[[1]]}}
  #Get taxonomic information
  taxon_list$genus_id <- NA
  taxon_list$genus_name <- NA
  taxon_list$family_id <- NA
  taxon_list$family_name <- NA
  for (j in 1:nrow(taxon_list)){
    if (!is.na(taxon_list$TaxID[j])){
      sciid <- as.character(taxon_list$TaxID[j])
      classification <- classification(sci_id=sciid, db="ncbi")[[sciid]]
    } else if (!is.na(taxon_list$GBIF_taxon_ID[j])) {
      sciid <- as.character(taxon_list$GBIF_taxon_ID[j])
      classification <- classification(sci_id=sciid, db="gbif")[[sciid]]
    } else {
      next
    }
    if (taxon_list_name=="michel_grauers") { #For Michel's data only get info on family ranks
      taxon_list$family_id[j] <- classification$id[classification$rank=="family"]
      taxon_list$family_name[j] <- classification$name[classification$rank=="family"]
    } else {
    taxon_list$genus_id[j] <- classification$id[classification$rank=="genus"]
    taxon_list$genus_name[j] <- classification$name[classification$rank=="genus"]
    taxon_list$family_id[j] <- classification$id[classification$rank=="family"]
    taxon_list$family_name[j] <- classification$name[classification$rank=="family"]
    }
  }
  return(taxon_list)
}

#Update all references with correct names, NCBI taxIDs or GBIF taxIDs (where the NCBI one are not available)
#Finally, get taxonomic information for all taxa having at list one type of ID
#for (i in diet_ref_list){
#  taxon_list <- get(i)
#  taxon_list_name <- i
#  assign(i, diet_db_2(taxon_list, taxon_list_name))
#}

#### Look at the database and do some formatting ####
#How many taxa for each subspecies
length(unique(append(remis_western$genus_id, rogers_western$genus_name))) 
length(unique(append(remis_western$family_id, rogers_western$family_name))) 

length(unique(append(michel_grauers$family_id, yamagiwa_grauers$family_name)))

length(unique(rothman_mountain$genus_name)) 
length(unique(rothman_mountain$family_name))

diet_database <- data.frame()
#### Combine and save database ####
#for (i in diet_ref_list) {
# pub=str_remove(i, "_.*")
# subspecies=str_remove(i, ".*_")
# table <- get(i)
# table$publication <- pub
# table$subspecies <- subspecies
# diet_database <- rbind(diet_database, table)
#}

#write.table(diet_database, "diet_database.txt", quote=FALSE, sep=",", row.names=FALSE)
diet_database <- read.table("diet_database.txt", quote=FALSE, sep=",")

## Create joint table with all taxa mentioned in the reference
#first, collect the family name column from each table is a list
diet_ref_fams <- list()
for (i in diet_ref_list) {
  diet_ref_fams[[paste0(i)]] = unique(get(i)$family_name)
}

diet_ref_fams <- t(splitstackshape:::charMat(listOfValues = diet_ref_fams, fill = 0L))
colnames(diet_ref_fams) <- diet_ref_list
diet_ref_fams <- as.data.frame(diet_ref_fams)

#Merge columns refering to the same gorilla subspecies
diet_ref_fams$western <- as.logical(diet_ref_fams$rogers_western) | as.logical(diet_ref_fams$remis_western) 
diet_ref_fams$grauers <- as.logical(diet_ref_fams$yamagiwa_grauers) | as.logical(diet_ref_fams$michel_grauers) 
diet_ref_fams$mountain <- as.logical(diet_ref_fams$rothman_mountain)

diet_ref_fams <- diet_ref_fams[,6:8]

#Do the same at the genus level
diet_ref_gen <- list()
for (i in diet_ref_list) {
  diet_ref_gen[[paste0(i)]] = unique(get(i)$genus_name)
}

diet_ref_gen <- t(splitstackshape:::charMat(listOfValues = diet_ref_gen, fill = 0L))
colnames(diet_ref_gen) <- diet_ref_list
diet_ref_gen <- as.data.frame(diet_ref_gen)

#Merge columns refering to the same gorilla subspecies
diet_ref_gen$western <- as.logical(diet_ref_gen$rogers_western) | as.logical(diet_ref_gen$remis_western) 
diet_ref_gen$grauers <- as.logical(diet_ref_gen$yamagiwa_grauers) | as.logical(diet_ref_gen$michel_grauers) 
diet_ref_gen$mountain <- as.logical(diet_ref_gen$rothman_mountain)

diet_ref_gen <- diet_ref_gen[,6:8]


#Stop logging
sink(file = NULL)

save.image()
