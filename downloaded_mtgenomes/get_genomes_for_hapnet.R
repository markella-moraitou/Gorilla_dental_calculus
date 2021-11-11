#Get genomes from van der Valk 2018 that correspond to the samples of my analysis
#Then add one genome per haplotype
library(dplyr)
library(readxl)
library(tidyverse)

#Get list of DC assembled genomes (this analysis)
DC_samples <- read.table("mt_genomes_samples.txt", header=FALSE) %>% pull(V1)

#Get list of accession numbers of downloaded gorilla genomes (these include vdV2018)
downloaded_genomes <- read.table("gorilla_mt_genomes_to_download.txt", sep='\t', header=TRUE)

#Get full metadata list
metadata_full <-  read_excel("/crex/proj/sllstore2017021/nobackup/MARKELLA/Sample_list_dental_calculus v.22.06.2021.xlsx") %>%
  #Merge Jena and Uppsala names in one column
  mutate(Seq.label=ifelse(is.na(Seq.label), Jena.lab, Seq.label)) %>% mutate(Jena.lab=NULL)

#### Get matching genomes ####

#Get list of deep seq assembled genomes (van der Valk 2018)
deepseq_samples <-  read_excel("/proj/sllstore2017021/nobackup/MARKELLA/published_work/van der Valk (2018) Significant loss of mitochondrial diversity - Supplementary Tables.xlsx", sheet="S1", skip=1) %>% filter(!is.na(`Individual ID`))

deepseq_samples <- deepseq_samples %>%
  #Match seq labels from my dataset using museum accessions
  mutate(DC.label=metadata_full$Seq.label[match(str_remove(deepseq_samples$`Museum catalog number`, "[A-Z]+_"),
                                                str_remove(metadata_full$AccessionNo, "-M.{0,1}$"))])
                                                                                
#Match my labeling to Tom's labeling to accession number
duplicate_genomes <- deepseq_samples  %>%
  #Filter to keep only those that represent DC_samples
  filter(DC.label %in% DC_samples) %>% 
  #Get only the relevant columns
  select(`Individual ID`, `mtDNA haplotype`, DC.label)
  
#Add accession number
duplicate_genomes$accession <- downloaded_genomes$accession[match(duplicate_genomes$`Individual ID`, downloaded_genomes$individual)]

duplicate_genomes <- duplicate_genomes %>% filter(!is.na(accession))

write.table(duplicate_genomes, file="duplicate_genomes.txt", sep='\t', quote=FALSE, row.names=FALSE)

#### Update gorilla_mt_genomes_to_download.txt list with more detailed species classification from the vdV2018 Supp
downloaded_genomes$organism <- ifelse(downloaded_genomes$individual %in% deepseq_samples$`Individual ID`,
                                      deepseq_samples$Species[match(downloaded_genomes$individual, deepseq_samples$`Individual ID`)],
                                      downloaded_genomes$organism)
                                     
#Save updated version
write.table(downloaded_genomes, file="gorilla_mt_genomes_to_download.txt", sep='\t', row.names=FALSE, quote=FALSE)

#### Get one sample per haplotype ####

#Get one sample per haplotype, the one with the highest coverage
sample_per_haplotype <-
  deepseq_samples %>% mutate(accession=downloaded_genomes$accession[match(deepseq_samples$`Individual ID`, downloaded_genomes$individual)]) %>%
  #Remove rows with missing accessions
  filter(!is.na(accession)) %>%
  group_by(`mtDNA haplotype`) %>% mutate(max_coverage=max(`mtDNA coverage`)) %>% ungroup %>%
  filter(`mtDNA coverage`==max_coverage) %>% select(`Individual ID`, `mtDNA haplotype`, DC.label, accession)
  
#Add accession number
sample_per_haplotype$accession <- downloaded_genomes$accession[match(sample_per_haplotype$`Individual ID`, downloaded_genomes$individual)]

write.table(sample_per_haplotype, file="sample_per_haplotype.txt", sep='\t', quote=FALSE, row.names=FALSE)


  
  