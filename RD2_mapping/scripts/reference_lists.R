#load packages
library(stringr)
library(phyloseq)
library(taxize)
library(dplyr)

load("/proj/sllstore2017021/nobackup/MARKELLA/T3_community-level/.RData")

#Checking the reference genomes table produced on Rackham
#and choosing substitute species for the missing taxa

#### For the contaminants ####

contam_ref_recent <- read.table("/proj/sllstore2017021/nobackup/MARKELLA/RD2_mapping/contaminant_genomes_list_recent.txt", sep="\t", fill=TRUE, row.names = NULL, header = FALSE)
contam_ref_recent$V1 <- as.character(contam_ref_recent$V1)

#Keep only unique entries
contam_ref_recent <- unique(contam_ref_recent)

#Compare with contaminant list
print("How many contaminant taxa are represented in the genomes list?")
length(intersect(contam_ref_recent$V1, exogenous_id_list))
print("How many are missing?")
length(setdiff(exogenous_id_list, contam_ref_recent$V1))

#Make a list of taxa that require a genus level reference
#First get missing taxa IDs
contam_refs_to_add <- exogenous_id_list[which(exogenous_id_list %in%
                                                   setdiff(exogenous_id_list, contam_ref_recent$V1))]
#Then get species and genus info
contam_refs_to_add <- as.data.frame(contam_refs_to_add)
colnames(contam_refs_to_add) <- "TaxID"
contam_refs_to_add <- cbind(contam_refs_to_add, taxonomy_species[match(contam_refs_to_add$TaxID, rownames(taxonomy_species)),c(7,8)])

#Only keep taxa whose genus isn't already represented in contaminant ref genomes
contam_refs_to_add <- contam_refs_to_add[which(!(contam_refs_to_add$genus %in% taxonomy_species[match(contam_ref_recent$V1, rownames(taxonomy_species)),7])),]

print("Amongst the contaminant refs that are missing from the genome list, are there any genera that are part of the final dataset?")
#Check if one of these genera is also found in the final dataset
contam_refs_to_add[which(contam_refs_to_add$genus %in% tax_table(spe_data_final)[,7]),]

#Get all the species in that genus
contam_refs_to_add$all_species_in_genus <- ncbi_children(name=contam_refs_to_add$genus, out_type = "uid", ambiguous = FALSE)

#Extract elements from list into a separate object (which will be used to download the genomes)
exogenous_id_suppl <- c()

for (i in 1:length(contam_refs_to_add$all_species_in_genus)) {
  exogenous_id_suppl <- append(exogenous_id_suppl, contam_refs_to_add$all_species_in_genus[[i]])
}

#Remove the "original" species
exogenous_id_suppl <- exogenous_id_suppl[(which(!(exogenous_id_suppl %in% contam_refs_to_add$TaxID)))]

#Add a column with the species that is being "replaced"
exogenous_id_suppl <- as.data.frame(exogenous_id_suppl)
colnames(exogenous_id_suppl) <- "sister_species"
for (i in 1:nrow(exogenous_id_suppl)){
  exogenous_id_suppl$original_species[i] <- 
    contam_refs_to_add$TaxID[which(grepl(exogenous_id_suppl$sister_species[i], contam_refs_to_add$all_species_in_genus))]
}

exogenous_id_suppl <- exogenous_id_suppl[,c(2,1)]
print("How many additional contaminant taxa we will need substitutes for?")
nrow(exogenous_id_suppl)

#Save file
write.table(exogenous_id_suppl, file="/proj/sllstore2017021/nobackup/MARKELLA/RD2_mapping/exogenous_id_suppl.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep="\t")


#### For the noncontaminants ####

noncontam_ref_recent <- read.table("/proj/sllstore2017021/nobackup/MARKELLA/RD2_mapping/noncontaminant_genomes_list_recent.txt", sep="\t", row.names = NULL, header = FALSE)
noncontam_ref_recent$V1 <- as.character(noncontam_ref_recent$V1)

#Keep only unique entries
noncontam_ref_recent <- unique(noncontam_ref_recent)

#Compare with abundant taxa list
print("How many abundant noncontaminant taxa are represented in the genomes list?")
length(intersect(noncontam_ref_recent$V1, abundant_id_list))
print("How many are missing?")
length(setdiff(abundant_id_list, noncontam_ref_recent$V1))

#Make a list of taxa that require a genus level reference
#First get missing taxa IDs
noncontam_refs_to_add <- abundant_id_list[which(abundant_id_list %in%
                                                setdiff(abundant_id_list, noncontam_ref_recent$V1))]

#Then get species and genus info
noncontam_refs_to_add <- as.data.frame(noncontam_refs_to_add)
colnames(noncontam_refs_to_add) <- "TaxID"
noncontam_refs_to_add <- cbind(noncontam_refs_to_add, taxonomy_species[match(noncontam_refs_to_add$TaxID, rownames(taxonomy_species)),c(7,8)])

#Only keep taxa whose genus isn't already represented in contaminant ref genomes
noncontam_refs_to_add <- noncontam_refs_to_add[which(!(noncontam_refs_to_add$genus %in% taxonomy_species[match(noncontam_ref_recent$V1, rownames(taxonomy_species)),7])),]

#Remove Homo sapiens
noncontam_refs_to_add <- noncontam_refs_to_add %>% filter(genus!="Homo")

#Get all the species in that genus
noncontam_refs_to_add$all_species_in_genus <- ncbi_children(name=noncontam_refs_to_add$genus, out_type = "uid", ambiguous = FALSE)

#Extract elements from list into a separate object (which will be used to download the genomes)
abundant_id_suppl <- c()

for (i in 1:length(noncontam_refs_to_add$all_species_in_genus)) {
  abundant_id_suppl <- append(abundant_id_suppl, noncontam_refs_to_add$all_species_in_genus[[i]])
}

#Remove the "original" species
abundant_id_suppl <- abundant_id_suppl[(which(!(abundant_id_suppl %in% noncontam_refs_to_add$TaxID)))]

#Add a column with the species that is being "replaced"
abundant_id_suppl <- as.data.frame(abundant_id_suppl)
colnames(abundant_id_suppl) <- "sister_species"
for (i in 1:nrow(abundant_id_suppl)){
  abundant_id_suppl$original_species[i] <- 
    noncontam_refs_to_add$TaxID[which(grepl(abundant_id_suppl$sister_species[i], noncontam_refs_to_add$all_species_in_genus))]
}

abundant_id_suppl <- abundant_id_suppl[,c(2,1)]
print("Amongst the noncontaminant refs that are missing from the genome list, are there any genera that are part of the final dataset?")
nrow(abundant_id_suppl)

#There should be no intersect
print("There should be no taxa in common amongst the taxa substituting contaminants and noncontaminants") 
intersect(abundant_id_suppl$sister_species, exogenous_id_suppl$sister_species)

#Save file
write.table(abundant_id_suppl, file="/proj/sllstore2017021/nobackup/MARKELLA/RD2_mapping/abundant_id_suppl.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep="\t")

