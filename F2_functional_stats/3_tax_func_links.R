sink(file="log3_tax_func_links.txt")

#load packages
library(reshape2)
library(phyloseq)
library(plyr)
library(stringr)
library(tidyr)
library(taxize)
library(ggpubr)
library(RColorBrewer)
library(cowplot)
library(tidyverse)

load(".RData")

#Functional analysis - Script 3
#Get taxa that contribute to differentially abundant functions


#### Manipulate the data ####

#Import stratified feature table
#Import stratified GO table
GO_str <- read_tsv("/proj/sllstore2017021/nobackup/MARKELLA/F1_HUMAnN2/all_GOterms_cpm_renamed.tsv")
GO_str <- as.data.frame(GO_str)

#Rename columns
colnames(GO_str) <- str_remove(colnames(GO_str), "_Abundance-RPKs")

#Turn first column to row names
row.names(GO_str) <- GO_str$`# Gene Family`
GO_str <- GO_str[,-1]

#Remove unmapped and ungrouped
GO_str <- GO_str[which(!(grepl("UNGROUPED", rownames(GO_str)))),]
GO_str <- GO_str[which(!(grepl("UNMAPPED", rownames(GO_str)))),]

GO_str_diffabund <- matrix(nrow=0, ncol=ncol(GO_str))
colnames(GO_str_diffabund) <- colnames(GO_str)

#Subset to keep only differentially abundant BP
list <- substring(most_signif_BP$feature, 0, 10)
list <- sub(".", ":", list, fixed=TRUE)
for (i in 1:nrow(significant_BP)) {
  subset <- GO_str[which(grepl(list[i], rownames(GO_str))),]
  GO_str_diffabund <- rbind(GO_str_diffabund, subset)
}

#Also drop unstratified features
GO_str_diffabund <- GO_str_diffabund[which(grepl("|", rownames(GO_str_diffabund), 
                                                 fixed=TRUE)),]

#Remove blanks
GO_str_diffabund <- GO_str_diffabund[, -which(grepl("[BE][ELR]", colnames(GO_str_diffabund)))]


#### Summarize data by genus and BP ####
#Get a table showing the BP, the genus and the average abundance in copies per million
BP_genera <- data.frame(abundance_mean = rowMeans(GO_str_diffabund))
BP_genera$BP <- substring(rownames(BP_genera), 18) %>% 
                          str_remove(pattern=".g__.*|.uncl.*")
#could be a genus classification or empty, if there was no classification
BP_genera$genus <- str_remove(rownames(BP_genera), ".*g__|.*uncl.*") %>%
                          str_remove(pattern=".s__.*")
#Turn these empty values to "unclassified"
BP_genera$genus[BP_genera$genus==""] <- NA

rownames(BP_genera) <- NULL
BP_genera <- BP_genera[,c(2,3,1)] #reorder columns

#Paste BP and genus (will split again after after summarizing)
BP_genera$group <- paste(BP_genera$BP, BP_genera$genus, sep="--")

#Sum the abundances coming from the same genus and BP and then split columns back
BP_genera <- ddply(BP_genera, .(group), summarise, abundance_mean=sum(abundance_mean)) %>% 
  separate(col = group, into = c("BP", "genus"), sep = "--")

#Manually fix some names (keep only the first part, it's enough to get the phylum info)
BP_genera$genus[which(grepl("_", BP_genera$genus))] <-
  str_remove(BP_genera$genus[which(grepl("_", BP_genera$genus))], "_.*")

#Get proportion of contributions by dividing by the sum of each BP
BP_genera <- BP_genera %>% group_by(BP) %>% 
  mutate(proportion=abundance_mean/sum(abundance_mean)) %>% ungroup

#Get phylum info

#Create a separate table, to search each genus only once
#Read from files
small_taxonomy <- read_rds(file="small_taxonomy")
#small_taxonomy <- data.frame(genus = unique(BP_genera$genus))

#small_taxonomy$order <- NA
#for (i in 1:nrow(small_taxonomy)) {
#  if (small_taxonomy$genus[i]=="NA") {
#    small_taxonomy$order[i] <- "unclassified"
#  } else {
#    classification <- classification(sci_id = small_taxonomy$genus[i], db = "ncbi")[[1]]
#    order=classification$name[which(classification$rank=="order")]
#    #In case the "order" field is empty, use 
#    small_taxonomy$order[i] <- ifelse(identical(order, character(0)), "unclassified", order)
#  }
#}

#Use small taxonomy to fill in the data frame
BP_genera$order <- small_taxonomy$order[match(BP_genera$genus, small_taxonomy$genus)]

#Order data by order
BP_genera <- BP_genera[order(BP_genera$order),]

write_rds(small_taxonomy, "small_taxonomy")

####Plot contribution of each order to BP ####

#Keep only BP which the biggest effect
signif_BP_genera <- BP_genera[which(BP_genera$BP %in% str_remove(rownames(GO_unstr_signif), ".*] ")),]

#Keep only orders presented in the taxonomic ancom plot + the 2 orders contributing the most abundance and change everything else to 'Other'
highlighted_genera <- levels %>%
  append(signif_BP_genera %>% group_by(order) %>% summarise(abundance=mean(abundance_mean)) %>% filter(abundance>0.5 & order!="unclassified") %>% pull(order))

signif_BP_genera <- 
  signif_BP_genera %>% mutate(order=ifelse(order %in% append("unclassified", highlighted_genera), order, "Other"))

#Turn order into factor to change legend order
signif_BP_genera$order <- factor(signif_BP_genera$order, 
                                  levels = rev(append("unclassified", highlighted_genera)))
                                  
#Order to match heatmap y-axis
signif_BP_genera$BP <- factor(signif_BP_genera$BP, levels(go_heatmap$data$YYYY))
signif_BP_genera <- signif_BP_genera[order(signif_BP_genera$BP),]

tax_fun_comp <- ggbarplot(signif_BP_genera, x="BP", y="proportion", fill="order", color="order")


#Manually fix palette
tax_func_palette <- ancom_spe_palette
#Add order names not included in diff abundant taxa heatmap
names(tax_func_palette)[is.na(names(tax_func_palette))] <- setdiff(highlighted_genera, names(tax_func_palette))
tax_func_palette["unclassified"] <- "grey34"

tax_fun_comp <- 
  tax_fun_comp + scale_fill_manual(values=tax_func_palette, name="Microbial order") + 
  scale_color_manual(values=tax_func_palette, name="Microbial order") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right", legend.text = element_text(size=8), legend.title = element_text(size=10),
        axis.text.x = element_text(angle=90)) +
  coord_flip() 
tax_fun_comp

tax_func_grid <- plot_grid(go_heatmap + 
                              theme(legend.position="top") +
                              scale_fill_gradient2(low="blue", mid="white", high="red",
  guide = guide_colorbar(barheight = 2, title = "CLR-normalized\nabundance"),
    limits=c(min(clr.pseudozeros.f), NA), na.value="grey"),
                          tax_fun_comp, align = "h", axis="tb",
          ncol = 2, rel_widths = c(25, 10))

ggsave(tax_func_grid, file="/proj/sllstore2017021/nobackup/MARKELLA/F2_functional_stats/tax_func_grid.png",
       height=8, width=15)
       
sessionInfo()

sink(file=NULL)
save.image()
