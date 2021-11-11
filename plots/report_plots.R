#load packages
library(ggpubr)
library("gridExtra")
library(phyloseq)
library(cowplot)
library(microbiome)
library(dplyr)
library(ape)
library(vegan)
library(ggtree)
library(tibble)
library(MutationalPatterns)
load("/proj/sllstore2017021/nobackup/MARKELLA/T3_community-level/.RData")

#### Figure 2 - Dataset summary ####
ggsave(
sample_data(spe_data_final) %>% as.matrix %>% as.data.frame %>%
  select(Spec.subspecies, Seq.centre) %>% table %>% as.data.frame %>% 
  #rearrange order
  mutate(Spec.subspecies=factor(Spec.subspecies, levels=c("gorilla", "graueri", "beringei"))) %>% arrange(Spec.subspecies) %>%
  ggbarplot(x="Spec.subspecies", y="Freq", fill="Seq.centre", color="Seq.centre") +
  scale_x_discrete(labels=c("gorilla" = "Western lowland", "graueri" = "Grauer's", "beringei" = "Mountain")) +
  xlab("Gorilla subspecies") + ylab("Number of samples") +
  scale_fill_brewer(palette="Pastel1", name="Dataset", labels = c("Fellows Yates et al. 2021", "Uppsala")) +
  scale_colour_brewer(palette="Pastel1", name="Dataset", labels = c("Fellows Yates et al. 2021", "Uppsala")) + theme_bw(),
  file="metadata_summary_barplot.png",
  device="png")  

#### Figure 3 - Alpha diversity ####
chao_alpha <- ggboxplot(alpha, y = "Chao1", x = "Spec.subspecies", fill = "Spec.subspecies", title="a)") +
  theme_light() + theme(legend.text = element_text(size=12), legend.title = element_text(size=12),
                     legend.position = "none", axis.text = element_text(size=12), axis.title=element_text(size=12)) +
  ylab("Chao1 estimator") + xlab("Host subspecies") +
  scale_x_discrete(labels=c("gorilla" = "Western", "graueri" = "Grauer's", "beringei" = "Mountain")) +
  scale_fill_brewer(palette="Set2")

ggsave(chao_alpha,
  file="chao_boxplot.png",
  device="png",
  height=5,
  width=7)
  
shan_alpha <- ggboxplot(alpha, y = "Shannon", x = "Spec.subspecies", fill = "Spec.subspecies", title="b)") +
  theme_light() + theme(legend.text = element_text(size=12), legend.title = element_text(size=12),
                     legend.position = "none", axis.text = element_text(size=12), axis.title=element_text(size=12)) +
  ylab("Shannon index") + xlab("Host subspecies") +
  scale_x_discrete(labels=c("gorilla" = "Western", "graueri" = "Grauer's", "beringei" = "Mountain")) +
  scale_fill_brewer(palette="Set2")


ggsave(shan_alpha,
  file="shannon_boxplot.png",
  device="png",
  height=5,
  width=7)

ggsave(
  plot_grid(chao_alpha, shan_alpha, nrow=2, align="v"),
  file="alpha_boxplots.png",
  device="png",
  height=10,
  width=7)


#### Figure 4 - Taxonomic ordinations ####
sample_data(spe_data_final)$Spec.subspecies <- factor(sample_data(spe_data_final)$Spec.subspecies, levels=c("gorilla", "graueri", "beringei"))
sample_data(spe_data_final_norm)$Spec.subspecies <- factor(sample_data(spe_data_final_norm)$Spec.subspecies, levels=c("gorilla", "graueri", "beringei"))


tax_aitchison <- plot_ordination(spe_data_final_norm, clr_5, color="Spec.subspecies", shape = "Seq.centre", title="b) PCoA on Aitchison distances") + 
  theme_light() + geom_point(size=4) + theme(legend.text = element_text(size=14), legend.title = element_text(size=14),
                                           legend.position = "right", plot.title = element_text(size=15)) +
  scale_colour_brewer(palette="Set2", name="Host subspecies", labels = c("beringei"="Mountain", "gorilla"="Western lowland", "graueri"="Grauer's")) +
  scale_shape_discrete(name="Dataset", labels = c("Jena"="Fellows Yates et al. 2021"))

tax_jaccard <- plot_ordination(spe_data_final, jaccard_5, color="Spec.subspecies", shape = "Seq.centre", title="a) PCoA on Jaccard distances") + 
  theme_light() + geom_point(size=4) + theme(legend.text = element_text(size=14), legend.title = element_text(size=14),
                                        legend.position = "none", plot.title = element_text(size=15))  +
  scale_colour_brewer(palette="Set2", name="Host subspecies", labels = c("beringei"="Mountain", "gorilla"="Western lowland", "graueri"="Grauer's")) +
  scale_shape_discrete(name="Dataset", labels = c("Jena"="Fellows Yates et al. 2021"))


#Final plot
tax_ordinations <- grid.arrange(tax_jaccard, tax_aitchison, ncol=2, widths=c(2.3, 3.2))

ggsave(tax_ordinations, file="tax_ordinations.png", device="png", height=7, width=14)
 
#### Figure 5 - Differentially abundant taxa ####
diff_abund_heat <-
  heat(ancom_species_table, Yvar="taxon_name", fill="clr-abundance", Xvar = "sample", order.cols = FALSE, order.rows = TRUE)+
  scale_fill_gradient2(low="blue", mid="white", high="red",
    guide = guide_colorbar(barheight = 10, title = "clr-normalized\nabundance", draw.ulim=FALSE),
    limits=c(min(clr.pseudozeros), NA), na.value="grey") +
  facet_grid(~subspecies, scales="free", labeller=facet.labeller) +
  theme(plot.title = element_text(size = 20, hjust = 0.8), axis.text.x = element_text(angle = 90, vjust = -0.2),
        axis.title = element_text(size=15), axis.text.y = element_text())

levels <- ancom_species_table %>% pull(taxon_order) %>% unique %>% as.character
#Bring 'Other' to be the first in order
levels <- append("Other", levels[-which(levels=="Other")])

#### Colour the labels by order ####
heat_ancom_ylabs <- levels(diff_abund_heat$data$YYYY)
heat_ancom_ylabs <- cbind(heat_ancom_ylabs, as.character(ancom_species_table$taxon_order[match(heat_ancom_ylabs, ancom_species_table$taxon_name)]))
heat_ancom_ylabs <- as.data.frame(heat_ancom_ylabs)
colnames(heat_ancom_ylabs) <- c("species", "order")
heat_ancom_ylabs$order <- factor(heat_ancom_ylabs$order, levels=levels)

#Create a palette
ancom_spe_palette <- c("#999999", "cadetblue1", "navy", "chartreuse3", "darkolivegreen", "gold","lightsalmon", "firebrick1",
                       "darkorchid3", "hotpink", "plum2", "dodgerblue", "khaki", "darkorange", "violetred")
names(ancom_spe_palette) <- levels(heat_ancom_ylabs$order)

#Get species-level palette
ancom_spe_palette.s <- heat_ancom_ylabs
ancom_spe_palette.s$colour <- ancom_spe_palette[match(ancom_spe_palette.s$order, names(ancom_spe_palette))]
ancom_spe_palette.s <- t(ancom_spe_palette.s[, -which(colnames(ancom_spe_palette.s)=="order")])

#Set taxon names as palette names
colnames(ancom_spe_palette.s) <- ancom_spe_palette.s[1,]
ancom_spe_palette.s <- ancom_spe_palette.s[-1,]

#Update plot
diff_abund_heat <-
  diff_abund_heat +  theme(axis.ticks.y = element_line(colour = ancom_spe_palette.s, size=2.5),
                           axis.ticks.length.y = unit(0.5, "cm"),
                           axis.text.y = element_text()) +
                     ylab("Differentially abundant microbial species")

ggsave(diff_abund_heat, file="diff_abund_heat.png", device="png", height=9, width=7)

#Create legend separately
ancom_legend <- ancom_spe_palette[!is.na(names(ancom_spe_palette))]
png("diff_abund_heat_legend.png")
  plot.new()
  legend("center", rev(names(ancom_legend)), col = rev(unname(ancom_legend)), pch = 20, cex = 1, title="Microbial order")
dev.off()

#### Figure 6 - Subspecies-associated biological processes ####
load("/proj/sllstore2017021/nobackup/MARKELLA/F2_functional_stats/.RData")
go_heatmap <-
  heat(GO_unstr_signif_m, "Sample", "GO term","Copies/million", order.cols = FALSE, order.rows=TRUE) +
  scale_fill_gradient2(low="blue", mid="white", high="red",
  guide = guide_colorbar(barheight = 10, title = "CLR-normalized\nabundance"),
    limits=c(min(clr.pseudozeros.f), NA), na.value="grey") +
  facet_grid(cols=vars(host_subspecies), scales="free", labeller=facet.labeller) +
  theme(plot.title = element_text(size = 20, hjust = 0.8), axis.text.x = element_text(angle = 90, vjust = -0.2),
        axis.title = element_text(size=15), axis.text.y = element_text(size=10), legend.position="right")

#Plot order contributions
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
tax_func_grid <- plot_grid(go_heatmap +
                              theme(legend.position="top") +
                              scale_fill_gradient2(low="blue", mid="white", high="red",
  guide = guide_colorbar(barheight = 2, title = "CLR-normalized\nabundance"),
    limits=c(min(clr.pseudozeros.f), NA), na.value="grey"),
                          tax_fun_comp, align = "h", axis="tb",
          ncol = 2, rel_widths = c(25, 10))

ggsave(tax_func_grid, file="BP_heatmap_with_tax_contributions.png",
       height=8, width=15)

#### Figure 7 - Dietary heatmap ####
load("/proj/sllstore2017021/nobackup/MARKELLA/D3_diet_stats/.RData")

#Labeller for the facet labels
facet.labels <- list("gorilla"="western lowland", "graueri"="Grauer's", "beringei"="mountain")
facet.labeller <- function(variable,value){
  return(facet.labels[value])
}

heat_diet <- heat(ancom_euk_abund, Yvar="genus", fill="clr-abundance", Xvar = "sample", order.cols = FALSE, order.rows = FALSE) +
  theme(legend.position = "left", plot.title=element_text(hjust = 0.5)) +
  scale_fill_gradient2(low="blue", mid="white", high="red",
    guide = guide_colorbar(barheight = 2, title = "clr-normalized\nabundance", draw.ulim=FALSE),
    limits=c(min(clr.pseudozeros.dd), NA), na.value="grey") +
    facet_grid(~host_subspecies, scales="free", labeller=facet.labeller) +
  ggtitle("Abundances of differentially abundant gorilla dietary genera across samples")

#### Colour the labels by phylum ####
heat_diet_ylabs <- levels(heat_diet$data$YYYY)
heat_diet_ylabs <- cbind(heat_diet_ylabs, tax_table(euk_genus_diet)[match(heat_diet_ylabs, taxa_names(euk_genus_diet)),3])
heat_diet_ylabs <- as.data.frame(heat_diet_ylabs)
colnames(heat_diet_ylabs) <- c("genus", "phylum")
heat_diet_ylabs$phylum <- factor(heat_diet_ylabs$phylum)

#Also add family infomation (will be needed later)
heat_diet_ylabs$family <- tax_table(euk_genus_diet)[match(heat_diet_ylabs$genus, taxa_names(euk_genus_diet)),6]

#Save numbers of genera per family
heat_diet_fam_freqs <- table(heat_diet_ylabs$family)
heat_diet_fam_freqs <- as.data.frame(heat_diet_fam_freqs)
heat_diet_fam_freqs <- heat_diet_fam_freqs[match(unique(heat_diet_ylabs$family), heat_diet_fam_freqs$Var1),]

#Subset diet_palette.g to get colours for the taxa in this heatmap
diet_palette.g <- diet_palette_full.g[match(heat_diet_ylabs$genus, colnames(diet_palette_full.g))]
names(diet_palette.g) <- heat_diet_ylabs$genus
diet_palette.g <- t(diet_palette.g)

#Also draw lines to show genera that belong to the same family
#From the sorted table with the number of genera per family, get the cumulative sum
#to determine where every line will be along the y axis
heat_diet_fam_freqs$Freq <- cumsum(heat_diet_fam_freqs$Freq)

#Update plot
heat_diet <-
  heat_diet + theme(axis.text.y = element_text(colour = diet_palette.g),
                    legend.position = "top") +
  geom_hline(yintercept=heat_diet_fam_freqs$Freq + 0.5, size=0.5, colour="white")

#### Create a sidebar ####
#based on the families in the studies
ref_sidebar <- heat_diet_ylabs[,c("genus", "family")]
ref_sidebar <- cbind(ref_sidebar, rownames_to_column(diet_ref_fams)[match(ref_sidebar$family, rownames_to_column(diet_ref_fams)$rowname),])
rownames(ref_sidebar) <- NULL
#Add info about genus mentions in literature
ref_sidebar <- cbind(ref_sidebar, rownames_to_column(diet_ref_gen)[match(ref_sidebar$genus, rownames_to_column(diet_ref_gen)$rowname),])
#A bunch of NAs because of the genera that are not mentioned in the references
#turn them into FALSE (because the species is absent)
ref_sidebar[is.na(ref_sidebar)] <- FALSE

#Remove unnecessary colums
ref_sidebar <- ref_sidebar[,-which(colnames(ref_sidebar) %in% c("family", "rowname")),]

#Replace FALSE -> 0, TRUE -> 1
ref_sidebar[ref_sidebar==FALSE] <- 0

#Add columns referring to the same diet.
#This way if the exact genus is mentioned in the literature it takes a value of 2
#But if it is only the family, it takes a value of 1

ref_sidebar$western <- ref_sidebar$western + ref_sidebar$western.1
ref_sidebar$grauers <- ref_sidebar$grauers + ref_sidebar$grauers.1
ref_sidebar$mountain <- ref_sidebar$mountain + ref_sidebar$mountain.1

#Remove extra columns
ref_sidebar <- ref_sidebar[, -which(grepl(".1", colnames(ref_sidebar)))]
rownames(ref_sidebar) <- NULL

#Melt
ref_sidebar <- reshape::melt(ref_sidebar)

colnames(ref_sidebar) <- c("diet_gen", "host_subspecies", "presence/absence")

#Set presence/absence as numeric
ref_sidebar$`presence/absence` <- as.numeric(ref_sidebar$`presence/absence`)

#Plot heatmap
heat_sidebar <- heat(ref_sidebar, Yvar="diet_gen", fill="presence/absence", Xvar = "host_subspecies", order.cols = FALSE, order.rows = FALSE) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.title=element_text(size=10, hjust = 0.5)) +
  guides(fill=FALSE) + scale_fill_gradient2(low = "grey", mid="palegreen2", high = "palegreen4", midpoint = 1) +
  geom_hline(yintercept=heat_diet_fam_freqs$Freq + 0.5, size=0.5, colour="white") +
  ggtitle("Mentioned\nin literature")

#Combine the two plots
heat_diet_complete <- plot_grid(heat_diet, heat_sidebar, align = "h", ncol = 2, rel_widths = c(45, 5),
                                axis="tb")

ggsave(heat_diet_complete, file="diet_heatmap.png", device="png")

#### Figure S1 - read count histogram ####
load("/proj/sllstore2017021/nobackup/MARKELLA/T3_community-level/.RData")

spe_data@sam_data$readcount_log <- log(spe_data@sam_data$readcount.m.before.Kraken)
#Check distribution of read
ggsave(file="read_count_raw_hist.png",
            gghistogram(spe_data@sam_data, "readcount_log", fill="Sample.type", position="stack") +
            geom_vline(xintercept= log(300000), size=1, colour="red"),
            device="png")

#### Figure S2 - compositions of samples + oral proportion histogram ####

#Compositions of samples
source_contribution <- plot_contribution(t(feast_output_env)) +
  theme(axis.text.x=element_text(angle = +90, hjust = 0)) +
  scale_fill_brewer(palette = "Spectral")

png(file = "source_contribution.png", width = 1000, height = 480)
source_contribution
dev.off()

#Oral proportion histogram
ggsave(file="oral_proporton_hist.png", 
        gghistogram(oral_proportion, "oral_proportion_log", fill="Sample.type", position="stack") +
        #Include vertical line with the cutoff of 0.03 (3%)
        geom_vline(xintercept=log(0.03), size=1, colour="red"))

#### Figure S3 - oral/contaminant proportion ####
decontam_prop_core_micr <- 
  heat(decontam_test, Xvar="uppsala_threshold", Yvar="jena_threshold", fill="prop_core_micr",
       order.rows=FALSE, order.cols=FALSE) +
    scale_fill_gradient2(guide = guide_colorbar(barheight = 10, title="core microbiome\ntaxa proportion")) +
    ylab("jena_threshold") + xlab("uppsala_threshold")

ggsave(decontam_prop_core_micr, file="decontam_prop_core_micr.png", device="png")

#### Figure S4 - taxonomic ordinations with duplicate samples ####
jaccard_with_duplicates <- jaccard_with_duplicates + theme_light() + geom_point(size=3) + theme(legend.position="none") + scale_colour_brewer(palette="Pastel1", name="Dataset", labels=c("Uppsala"="Uppsala", "Jena"="Fellows Yates et al. 2021"))

clr_with_duplicates <- clr_with_duplicates + theme_light() + geom_point(size=3) + scale_colour_brewer(palette="Pastel1", name="Dataset", labels=c("Uppsala"="Uppsala", "Jena"="Fellows Yates et al. 2021"))

ordinations_with_duplicates <- grid.arrange(jaccard_with_duplicates, clr_with_duplicates, ncol=2, widths=c(2.4, 3.2))

ggsave(ordinations_with_duplicates, file="ordinations_with_duplicates.png", device="png", height=7, width=14)

#### Figure S5 - taxonomic ordinations with altitude ####
jaccard_with_altitude <- plot_ordination(spe_data_final, jaccard_5, color="Spec.subspecies", shape = "subspecies_locality", title="a) Jaccard distances") +
  theme_light() + geom_point(size=3) + theme(legend.text = element_text(size=12), legend.title = element_text(size=12),
                                             legend.position = "none", plot.title = element_text(size=15))  +
  scale_colour_brewer(palette="Set2", name="Host subspecies", labels = c("beringei"="Mountain", "gorilla"="Western", "graueri"="Grauer's")) + 
  scale_shape_manual(values=c(16, 15, 17, 16), name="Grauer's altitude", labels = c("beringei"="Mountain", "gorilla"="Western", "graueri <1500"="Grauer's <1500m", "graueri >1500"="Grauer's >1500m"))

clr_with_altitude <- plot_ordination(spe_data_final_norm, clr_5, color="Spec.subspecies", shape = "subspecies_locality", title="b) Aitchison distances") + 
  theme_light() + geom_point(size=3) + theme(legend.text = element_text(size=12), legend.title = element_text(size=12),
                                             legend.position = "right", plot.title = element_text(size=15)) +
  scale_colour_brewer(palette="Set2", name="Host subspecies", labels = c("beringei"="Mountain", "gorilla"="Western", "graueri"="Grauer's")) + 
  scale_shape_manual(values=c(16, 15, 17, 16), name="Grauer's altitude", labels = c("beringei"="Mountain", "gorilla"="Western", "graueri <1500"="Grauer's <1500m", "graueri >1500"="Grauer's >1500m"))
  
ordinations_with_altitude <- grid.arrange(jaccard_with_altitude, clr_with_altitude, ncol=2, widths=c(2.4, 3.2))

ggsave(ordinations_with_altitude, file="ordinations_with_altitude.png", device="png", height=7, width=14)
  
  
# PLOTS NOT INCLUDED IN REPORT #
#### Dendrogram ####

#Aitchison
dendr_aitchison <- spe_data_final_norm@otu_table %>% t %>% vegdist("euclidean") %>% hclust("average") %>%
  as.phylo 

#Jaccard
dendr_jaccard <- spe_data_final@otu_table %>% t %>% vegdist("jaccard") %>% hclust("average") %>%
  as.phylo 

#Colour tips by subspecies
tipsubsp <- data.frame(name=dendr_jaccard$tip.label, subspecies=spe_data_final@sam_data$Spec.subspecies[match(dendr_jaccard$tip.label,
                                                                                                              sample_names(spe_data_final))])
ggsave(
grid.arrange(
  ggtree(dendr_aitchison) %<+% tipsubsp + geom_tippoint(aes(color=subspecies), size=4) +
  scale_colour_brewer(palette="Set2") + 
    theme(legend.position = "none") +
    ggtitle("Aitchison distances"),
  ggtree(dendr_jaccard) %<+% tipsubsp + geom_tippoint(aes(color=subspecies), size=4) +
    scale_colour_brewer(palette="Set2", 
                        name="Host subspecies", 
                        labels=c("beringei"="Mountain", "gorilla"="Western", "graueri"="Grauer's")) +
    ggtitle("Jaccard distances"),
  ncol=2, widths=c(2.3, 3.2)),
  file="dendrograms.png",
  device="png")

#### Decontam mapping database ####
#Load contaminant and noncontaminant genome lists
contaminant_genomes <- read.table("/proj/sllstore2017021/nobackup/MARKELLA/RD2_mapping/contaminant_genomes_list_recent.txt", sep="\t", comment.char="", header=F)
noncontaminant_genomes <- read.table("/proj/sllstore2017021/nobackup/MARKELLA/RD2_mapping/noncontaminant_genomes_list_recent.txt", sep="\t", comment.char="", header=F)

colnames(contaminant_genomes) <- c("taxon_ID", "RefSeq", "Assembly_level", "Genome_repres", "Date", "ftp" )
colnames(noncontaminant_genomes) <- c("taxon_ID", "RefSeq", "Assembly_level", "Genome_repres", "Date", "ftp" )

print("How many contaminant genomes?")
nrow(contaminant_genomes)

print("How many non-contaminant genomes?")
nrow(noncontaminant_genomes)

ggsave(
plot_grid(ggpie(as.data.frame(table(contaminant_genomes$Assembly_level)), x="Freq", label="Var1", fill = "Var1", color = "Var1",
                lab.pos="in", lab.font ="black", lab.adjust=10) +
            theme(legend.position = "none"),
          ggpie(as.data.frame(table(contaminant_genomes$RefSeq)), x="Freq", label="Var1", fill = "Var1", color = "Var1",
                lab.pos="in", lab.font ="black")+
            theme(legend.position = "none"),
          ggpie(as.data.frame(table(contaminant_genomes$Genome_repres)), x="Freq", label="Var1", fill = "Var1", color = "Var1",
                lab.pos="in", lab.font ="black")+
            theme(legend.position = "none"),
          ggpie(as.data.frame(table(noncontaminant_genomes$Assembly_level)), x="Freq", label="Var1", fill = "Var1", color = "Var1",
                lab.pos="in", lab.font ="black")+
            theme(legend.position = "none"),
          ggpie(as.data.frame(table(noncontaminant_genomes$RefSeq)), x="Freq", label="Var1", fill = "Var1", color = "Var1",
                lab.pos="in", lab.font ="black")+
            theme(legend.position = "none"),
          ggpie(as.data.frame(table(noncontaminant_genomes$Genome_repres)), x="Freq", label="Var1", fill = "Var1", color = "Var1",
                lab.pos="in", lab.font ="black")+
            theme(legend.position = "none"),
          nrow=2),
          file="decontam_db_references.png",
          device="png")
          

#### Functional ordinations ####

load("/proj/sllstore2017021/nobackup/MARKELLA/F2_functional_stats/.RData")

sample_data(GO_BP_phyloseq)$Spec.subspecies <- factor(sample_data(GO_BP_phyloseq)$Spec.subspecies, levels=c("gorilla", "graueri", "beringei"))

func_ordination <- plot_ordination(GO_BP_phyloseq, GO_eucl, color="Spec.subspecies", shape = "Seq.centre", title="PCoA on euclidean distances of functional profiles") + 
  theme_light() + geom_point(size=4) + theme(legend.text = element_text(size=14), legend.title = element_text(size=14),
                                           legend.position = "right", plot.title = element_text(size=15)) +
  scale_colour_brewer(palette="Set2", name="Host subspecies", labels = c("beringei"="Mountain", "gorilla"="Western lowland", "graueri"="Grauer's")) +
  scale_shape_discrete(name="Dataset", labels = c("Jena"="Fellows Yates et al. 2021"))

ggsave(func_ordination, file="func_ordination.png", device="png", height=7, width=9)
