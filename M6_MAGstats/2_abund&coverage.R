#load packages
library(stringr)
library(dplyr)
library(ggpubr)
library(funrar)
library(reshape2)
library(gridExtra)
library(RColorBrewer)

load(".RData")

#MAG relative abundance and coverage
mag_map <- read.csv(file="/proj/sllstore2017021/nobackup/MARKELLA/M5_samplestoMAGs/log_map_readct_210724.txt")

#Read readcount table
readcounts <- read.table(file="/proj/sllstore2017021/nobackup/MARKELLA/plots/readcounts.tsv", sep="\t", header=TRUE)

#Read metadata object
metadata <- readRDS(file="/proj/sllstore2017021/nobackup/MARKELLA/T3_community-level/metadata.RDS")

print("How many NAs in the table?")
table(is.na(mag_map)) #All the NAs are due to 0 depth of coverage for the given MAG
mag_map[is.na(mag_map)] <- 0

#Fix sample names
mag_map$sample <- str_remove(mag_map$sample, "./") %>% str_remove("_m_decontam")

#Include total read count per sample after decontam
mag_map <- mag_map %>% left_join(readcounts[,c(1,12)]) %>% mutate(rel_abund=read_count/decontam_map, .keep="all", decontam_map=NULL)

#Add mag classification
mag_map$MAG <- str_trim(mag_map$MAG) #strip whitespaces
mag_map$MAG_classif <- hq_mq_mags$tree_names[match(mag_map$MAG, hq_mq_mags$user_genome)]
mag_map <- mag_map[,c(1,2,6,3,4,5)] #reorder columns

#Make relative abundance a percentage
mag_map$rel_abund <- mag_map$rel_abund * 100

#Add host subspecies
mag_map$host_subspecies <- metadata$plot.label[match(mag_map$sample, rownames(metadata))]

#Label env controls and blanks as "Blank/Control"
mag_map <- mag_map %>%
  mutate(host_subspecies = ifelse(host_subspecies %in% c("Blank", "Museum_control"), "Blank/control", host_subspecies))
#Make factor
mag_map$host_subspecies <- factor(mag_map$host_subspecies, levels = c("Ggg", "Gbg", "Gbb", "Blank/control"))

#Order rows by samples
mag_map <- mag_map[order(mag_map$sample),]

#Add column to mark the sample from which a MAG was assembled
mag_map$assembled_in_sample <- "No"
mag_map <- mag_map %>% mutate(assembled_in_sample=ifelse(str_remove(MAG, "[:punct:][0-9]+")==sample, "Yes", "No"))


#ggscatter(mag_map, x="sample", y="average_coverage",
#          color="assembled_in_sample",
#          shape=) +
#  theme(axis.text.x = element_text(angle=90)) 

ggsave(
ggscatter(mag_map, x="sample", y="coverage_breadth",
          color="assembled_in_sample") +
  theme(axis.text.x = element_text(angle=90)),
  file="mag_coverage_per_sample_scatterplot.png",
  device="png",
  width=14,
  height=7)

#### Abundances of specific MAGs of interest ####

#mag_group: a grepable string to select a set of MAG 
#e.g mag_group = "Rothia gets "unknown.Rothia"   "unknown.Rothia.1" "unknown.Rothia.2" "unknown.Rothia.3"
tree_abund <- function(mag_group) {
  large_plot_list <- list()
  for (i in hq_mq_mags$tree_names[grepl(mag_group, hq_mq_mags$tree_names)]) {
    small_plot_list = list()
    #Plot relative abundance, average coverage and breath of coverage
    for (j in c("rel_abund", "coverage_breadth")) {
      #Get appropriate y label
      ylabel <- ifelse(j=="rel_abund", "Relative abundance",
                       ifelse(j=="average_coverage", "Average coverage",
                              "Coverage breadth"))
      #Which subspecies was it isolated from
      sample <- str_remove(hq_mq_mags$user_genome[hq_mq_mags$tree_names==i], "[:punct:][0-9]+")
      host_subspecies <- metadata$plot.label[which(rownames(metadata)==sample)]
      #Get the mean and SD per subspecies
      means <- mag_map  %>% filter(MAG_classif==i, host_subspecies!="Blank/control") %>%
        group_by(host_subspecies) %>% 
        #y indicates the position for the annotation: a bit above the largest value
        summarise(mean=round(mean(get(j)), digits = 2), sd=round(sd(get(j)), digits = 2),
                  y=max(get(j))*1.1) %>%
        #add x coordinate
        ungroup %>% mutate(x=rank(host_subspecies)) %>% mutate(y=max(y))

      plot <- mag_map  %>% filter(MAG_classif==i, host_subspecies!="Blank/control") %>%
        ggboxplot(y=j, x="host_subspecies", fill="host_subspecies", colour="assembled_in_sample") +
        theme_bw() + theme(legend.text = element_text(size=12), legend.title = element_text(size=15),
                           axis.text.x = element_text(size=24), legend.position = "none",
                           axis.title.x = element_blank(), axis.title.y = element_text(size=24)) +
        ylab(ylabel) + xlab("Host subspecies") + ggtitle(paste(gsub(".", " ", fixed=TRUE, i), "assembled in", host_subspecies, sep="\n")) +
        scale_x_discrete(labels=c("Ggg" = "Western", "Gbg" = "Grauer's", "Gbb" = "Mountain")) +
        scale_fill_brewer(palette="Set2", labels=c("Ggg" = "Western", "Gbg" = "Grauer's", "Gbb" = "Mountain")) +
        ylim(0, (mag_map %>% filter(MAG_classif==i, host_subspecies!="Blank/control") %>% pull(j) %>% max)*1.2) +
        annotate("text", label= paste0("mean = ", means$mean, "\nSD = ", means$sd), y=means$y, x=means$x, size=7)
      small_plot_list[[paste(i, j, sep="_")]] <- plot
    }
    large_plot_list[[i]] <- small_plot_list
  }
  return(large_plot_list)
}

#Save trees and relative abundance stats for Rothia and Streptococcus

#Rothia
pdf(file = "/crex/proj/sllstore2017021/nobackup/MARKELLA/M6_MAGstats/Rothia_MAGs.pdf", width=13, height=7)
rothia_ggtree
plots <- tree_abund("Rothia")
for (i in 1:length(plots)) {
  grid.arrange(ncol=2, grobs=plots[[i]])
}
dev.off()

#Streptococcus
pdf(file = "/crex/proj/sllstore2017021/nobackup/MARKELLA/M6_MAGstats/Streptococcus_MAGs.pdf", width=13, height=7)
strept_ggtree
plots <- tree_abund("Streptococcus")
for (i in 1:length(plots)) {
  grid.arrange(ncol=2, grobs=plots[[i]])
}
dev.off()

mag_map %>% filter(grepl("Rothia", MAG_classif) & host_subspecies!="Blank/control") %>% group_by(host_subspecies, MAG_classif) %>% summarise(median_rel_abund=median(rel_abund)) %>%
  ggballoonplot(y="MAG_classif", x="host_subspecies", size="median_rel_abund", fill="host_subspecies", size.range = c(1, 20)) +
  scale_fill_manual(values=c(RColorBrewer::brewer.pal(3, "Set2")))

