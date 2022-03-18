# sink(file="log2_maaslin_GO.txt")

#load packages
library(Maaslin2)
library(tidyverse)
library(microbiome)
library(ggpubr)
library(dplyr)

# load(".RData")

#Functional analysis - Script 2
#Differential abundance analysis with Maaslin2
GO_BP_phyloseq<-readRDS("F2_functional_stats/GO_BP_phyloseq")

# make pairs
combs<-c("gorilla","graueri","beringei")
results<-NULL
for (i in combs) {
  print(i)
  #Run Maaslin2 for gene ontology terms
  Maaslin2(input_data = t(otu_table(GO_BP_phyloseq)),normalization = "NONE", analysis_method = "LM",
           input_metadata = as.data.frame(as.matrix(sample_data(GO_BP_phyloseq))),
           fixed_effects = "Spec.subspecies", reference = paste0("Spec.subspecies,",i), plot_heatmap = TRUE,
           output = paste0("F2_functional_stats/maaslin/GO_maaslin_output-",i,"ref"))
  
  #Load data
  significant_BP <- read_tsv(paste0("F2_functional_stats/maaslin/GO_maaslin_output-",i,"ref/significant_results.tsv"))
  
  sig.BP<-length(unique(significant_BP$feature))
  
  #Only keep associations with a corrected p-value (qval) above  0.05
  most_signif_BP <- significant_BP[significant_BP$qval<0.05,]
  
  #How many unique biological processes is this?
  sig.BP.q<-length(unique(most_signif_BP$feature))

  results<-rbind(results,cbind(sig.BP,sig.BP.q))
}

results<-cbind(combs,results)
colnames(results)<-c("subspecies as reference","How many associated biological processes with qval < 0.05?","How many associated biological processes?")
write.table(x = results,file = "F2_functional_stats/comparison-of-subspecies-reference.txt",quote = F,row.names = F)

#Get the IDs of the most significant GO terms
most_signif_BP_ids <- unique(sort(substr(most_signif_BP$feature, 0, 10)))

#Use these IDs to subset the original table
#Create empty matrix
GO_unstr_signif <- as.data.frame(matrix(nrow = length(most_signif_BP_ids), ncol = ncol(GO_unstr)))
for (i in 1:length(most_signif_BP_ids)){
  subset <- GO_unstr[which(grepl(most_signif_BP_ids[i], rownames(GO_unstr))),]
  GO_unstr_signif[i,] <- subset
  rownames(GO_unstr_signif)[i] <- rownames(subset)
}
colnames(GO_unstr_signif) <- colnames(GO_unstr)

#Transform
GO_unstr_signif_z <- microbiome::transform(GO_unstr_signif, "clr")

#Set 0 values (the most negative value of CLR-normalized abundances) to min(clr.pseudocounts)*2 instead of NA (which is what Jaelle did). This is because the heatmap function that I use doesn't recognize NAs
clr.pseudozeros.f <- sapply(colnames(GO_unstr_signif_z), function(x){min(GO_unstr_signif_z[,x])})
for (s in colnames(GO_unstr_signif_z)) {
  GO_unstr_signif_z[which(GO_unstr_signif_z[,s] == clr.pseudozeros.f[s]),s] <- min(clr.pseudozeros.f)*2
}

#Order by subspecies
GO_unstr_signif_z <- GO_unstr_signif_z[,
                             order(sample_data(spe_data_final)$Spec.subspecies[match(colnames(GO_unstr_signif_z), sample_names(spe_data_final))])]

#Melt table
GO_unstr_signif_m <- reshape2::melt(t(GO_unstr_signif_z))
colnames(GO_unstr_signif_m) <- c("Sample", "GO term", "Copies/million")

#Group features based on association to subspecies
GO_codes <- sapply(GO_unstr_signif_m$`GO term`, function(x) {substr(x, 0, 10) %>% gsub(pattern=":", replacement=".")}) 

#Use original maaslin output to find features that are positively/negatively associated to Gbb, Gbg or both
GO_unstr_signif_m$association <- 
  ifelse(GO_codes %in% (filter(significant_BP, coef>1 & value=="beringei") %>% pull(feature) %>% substr(0, 10)), "pos.mountain",
         ifelse(GO_codes %in% (filter(significant_BP, coef<1 & value=="beringei") %>% pull(feature) %>% substr(0, 10)), "neg.mountain", 
         ifelse(GO_codes %in% (filter(significant_BP, coef>1 & value=="graueri") %>% pull(feature) %>% substr(0, 10)), "pos.grauers",
         ifelse(GO_codes %in% (filter(significant_BP, coef<1 & value=="graueri") %>% pull(feature) %>% substr(0, 10)), "neg.grauers", "else")))) %>% factor(level=c("pos.mountain", "pos.grauers", "neg.mountain", "neg.grauers"))

#Order by associations
GO_unstr_signif_m <- GO_unstr_signif_m[order(GO_unstr_signif_m$association),]

#Add host subspecies info
GO_unstr_signif_m$host_subspecies <- metadata$Spec.subspecies[match(GO_unstr_signif_m$Sample, rownames(metadata))]
GO_unstr_signif_m$host_subspecies <- factor(GO_unstr_signif_m$host_subspecies, levels = c("gorilla", "graueri", "beringei"))


#Keep only the human readable names
GO_unstr_signif_m$`GO term` <- substring(GO_unstr_signif_m$`GO term`, 18)

#Labeller for the facet labels
facet.labels <- list("gorilla"="western lowland", "graueri"="Grauer's", "beringei"="mountain")
facet.labeller <- function(variable,value){
  return(facet.labels[value])
}

#Plot heatmap
go_heatmap <- 
  heat(GO_unstr_signif_m, "Sample", "GO term","Copies/million", order.cols = FALSE, order.rows=TRUE) +
  scale_fill_gradient2(low="blue", mid="white", high="red",
  guide = guide_colorbar(barheight = 10, title = "CLR-normalized\nabundance"),
    limits=c(min(clr.pseudozeros.f), NA), na.value="grey") +
  facet_grid(cols=vars(host_subspecies), scales="free", labeller=facet.labeller) +
  theme(plot.title = element_text(size = 20, hjust = 0.8), axis.text.x = element_text(angle = 90, vjust = -0.2), 
        axis.title = element_text(size=15), axis.text.y = element_text(size=10), legend.position="right")
go_heatmap

ggsave(go_heatmap, file="F2_functional_stats/go_heatmap.png",
       height=6, width=9)

write.table(rev(unique(sort(go_heatmap$data$YYYY))), row.names = FALSE, quote = FALSE, file = "F2_functional_stats/GO_signif_BP.txt")


#### Compare abundances of differentially abundant and total functions ####
#Get taxa sums for all taxa in the final dataset

func_ancomVStotal <- taxa_sums(GO_BP_phyloseq)
func_ancomVStotal <- as.data.frame(cbind(func_ancomVStotal, 
                                         ifelse(substring(names(func_ancomVStotal), first=0, last=10) %in%
                                                  #Make the names match
                                                  (significant_BP$feature %>% sub(pattern=".", replacement=":", fixed=TRUE) %>%
                                                  substring(first=0, last=10)),
                                                  "Yes", "No")))
colnames(func_ancomVStotal) <- c("log_abundance", "subspecies_association")
func_ancomVStotal$log_abundance <- log(as.numeric(func_ancomVStotal$log_abundance))

#Plot
ggsave(
  ggboxplot(func_ancomVStotal, x="subspecies_association", y="log_abundance", fill="grey") +
    xlab("Association with host subspecies") + ylab("log-transformed abundance") +
    annotate("text", x=0.8, y=5, label=paste0("n = ", sum(func_ancomVStotal$subspecies_association=="No"))) +
    annotate("text", x=1.8, y=10, label=paste0("n = ", sum(func_ancomVStotal$subspecies_association=="Yes"))),
    file="func_ancomVStotal_boxplot.png",
    device="png")

ggsave(
  gghistogram(func_ancomVStotal, x="log_abundance", fill="subspecies_association"),
  file="func_ancomVStotal_histogram.png",
  device="png")

#### Check enriched BPs per function ####
print("Which subspecies has the most BPs associated with it?")
significant_BP$value %>% table

print("Mountain gorillas: Are these BPs mostly positively or negatively associated?")
significant_BP %>% filter(value=="beringei") %>% pull(coef) %>% summary

print("Grauer's gorillas: Which are these associated BPs?")
significant_BP %>% filter(value=="graueri") %>% as.data.frame

sessionInfo()

sink(file=NULL)
save.image()
