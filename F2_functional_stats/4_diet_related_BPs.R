sink(file="log4_diet_related_BPs.txt")

#Load packages
library(GOxploreR)
library(dplyr)
library(phyloseq)
library(stringr)
library(tibble)
library(reshape2)
library(microbiome)

load(".RData")

#Get all children node for the following general GO terms
#GO:0006109 - regulation of carbohydrate metabolic process
#GO:0051246 - regulation of protein metabolic process
#GO:0019216 - regulation of lipid metabolic process
#GO:0043610 - regulation of carbohydrate utilization
#GO:0043609 - regulation of carbon utilization
#GO:2000942 - regulation of amylopectin metabolic process
#GO:0051246 - regulation of protein metabolic process
#GO:0019740 - nitrogen utilization
#GO:0006791 - sulfur utilization
#GO:0009758 - carbohydrate utilization

BP_categories <- c("GO:0006109", "GO:0051246", "GO:0019216",
                   "GO:0043610", "GO:0043609", "GO:2000942",
                   "GO:0051246", "GO:0019740", "GO:0006791", 
                   "GO:0009758")

BP_diet_related <- c()
for (n in BP_categories) {
  BP_diet_related <- append(BP_diet_related,
                            GO2DecBP(goterm = n)) %>% unique
}

#Maaslin output
significant_BP$feature %>% 
  str_remove("...BP.*") %>%
  gsub(pattern=".", replacement=":", fixed = TRUE) %>% 
  intersect(BP_diet_related) -> significant_BP_diet

print("Which BP terms are related to dietary molecules?")

GO_BP_phyloseq %>% otu_table %>%
  as.data.frame %>%
  rownames_to_column %>%
  #Get separate columns with the GO ID and the GO human-readable name
  mutate(GO_ID=str_remove(rowname, ": .BP..*")) %>% 
  mutate(GO_name=str_remove(rowname, ".*BP.. *")) %>%
  #Keep only dietary-related GOs
  filter(GO_ID %in% significant_BP_diet) %>% select(GO_ID, GO_name)

#### Plot heatmap for these BPs ####

#Get abundance values
  GO_BP_phyloseq %>% microbiome::transform("Z") %>% otu_table %>%
  as.data.frame %>%
  rownames_to_column %>%
  #Get separate columns with the GO ID and the GO human-readable name
  mutate(GO_ID=str_remove(rowname, ": .BP..*")) %>% 
  mutate(GO_name=str_remove(rowname, ".*BP.. *")) %>%
  #Keep only dietary-related GOs
  filter(GO_ID %in% significant_BP_diet) %>%
  melt %>%
  #Plot heatmap
  heat(Xvar="variable", Yvar="GO_name", fill="value", order.rows=FALSE)
  

