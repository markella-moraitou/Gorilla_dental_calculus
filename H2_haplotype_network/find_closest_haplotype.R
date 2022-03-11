library(dplyr)
library(stringr)

#This script processes the output of FastANI to assign each incomplete genome to the most closely related haplotype

#Load FastANI output
ANItable <- read.table("H2_haplotype_network/fastani_output.txt", sep="\t", dec=".")
colnames(ANItable) <- c("query", "reference", "ANI", "bidirect_frag_mappings", "total_query_fragments")

#Filter for the highest ANI value per query
maxANItable <- ANItable %>% group_by(query) %>% mutate(maxANI=max(ANI)) %>% filter(ANI==maxANI)

#Remove unnecessary text
maxANItable$query <- str_remove(maxANItable$query, "_sort_MT_RG_angsd.trimmed.fasta.gz")
maxANItable$reference <- str_remove(maxANItable$reference, "_sort_MT_RG_angsd.trimmed.fasta.gz") %>% str_remove(".trimmed.fasta.gz")

write.table(maxANItable, file="H2_haplotype_network/fastani_output_maxANI.txt", 
            quote=FALSE, sep="\t", row.names=FALSE)
            