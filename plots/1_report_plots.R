#load packages
library(tidyverse)
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
library(ggh4x)
load("/home/adrian/PostDocWork/Gorilla_dental_calculus_zenodo/T3_community-level/.RData")

alt.fix<-read_tsv("../Gorilla_dental_calculus_zenodo/T3_community-level/new_altitude.csv") %>% 
  right_join(sample_data(spe_data_final) %>% data.frame() %>% rownames_to_column("Sample_ID")) %>%
  mutate(subspecies_locality=ifelse(Spec.subspecies=="graueri",
                                    ifelse(`Approximate altitude`>1000,"graueri >1000","graueri <1000"),
                                    as.character(subspecies_locality))) %>% 
  pull(subspecies_locality)

sample_data(spe_data_final)$subspecies_locality<-alt.fix

euk_genus_diet<-readRDS("../Gorilla_dental_calculus_zenodo/D3_diet_stats/euk_genus_diet.rds")
euk_genus_diet_norm<-readRDS("../Gorilla_dental_calculus_zenodo/D3_diet_stats/euk_genus_diet_norm.rds")

#### Figure 2 - Dataset summary ####
ggplot2::ggsave(
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

ggplot2::ggsave(chao_alpha,
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


ggplot2::ggsave(shan_alpha,
                file="shannon_boxplot.png",
                device="png",
                height=5,
                width=7)

ggplot2::ggsave(
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

ggplot2::ggsave(tax_ordinations, file="tax_ordinations.png", device="png", height=7, width=14)

#### Figure 5 - Differentially abundant taxa ####

pal<-c(get_palette(palette = "Set2",k=3)[1],
       get_palette(palette = "Set2",k=3)[2],
       get_palette(palette = "Set2",k=3)[2],
       get_palette(palette = "Set2",k=3)[3])

# create vector for ordering by high/low altitude
# luckily these get sorted in the right order
alt.list<-sample_data(spe_data_final) %>% 
  data.frame() %>% 
  rownames_to_column("sample") %>% 
  distinct(plot.label,sample,subspecies_locality) %>%
  arrange(plot.label,subspecies_locality) %>% 
  pull(sample)

facet.labels<-NULL
facet.labels$label.col <- "Microbial Order"
facet.labels$gorilla <- "Western Lowland"
facet.labels$`graueri >1000`<-"LA Grauer's"
facet.labels$`graueri <1000`<-"HA Grauer's"
facet.labels$beringei <- "Mountain"

# add genus identifier
ancom_species_table$genus <- ancom_species_table %>% separate(taxon_name,into = "genus",sep = " ") %>% pull(genus)

# fix problems with genus names
ancom_species_table[ancom_species_table$taxon_name=="Lactobacillus gastricus","taxon_name"]<-"Limosilactobacillus gastricus"
ancom_species_table[ancom_species_table$taxon_name=="Limosilactobacillus gastricus","taxon_order"]<-"Lactobacillales"
ancom_species_table[ancom_species_table$taxon_name=="Limosilactobacillus gastricus","genus"]<-"Limosilactobacillus"

ancom_species_table[ancom_species_table$genus=="Agrobacterium","taxon_order"]<-"Rhizobiales"
ancom_species_table[ancom_species_table$genus=="Rhizobium","taxon_order"]<-"Rhizobiales"

ancom_species_table <-
  ancom_species_table %>% 
  mutate(taxon_order=fct_relevel(taxon_order,"Bacillales",
                                 "Corynebacteriales",
                                 "Enterobacterales",
                                 "Hyphomicrobiales",
                                 "Lactobacillales",
                                 "Limosilactobacillus",
                                 "Micrococcales",
                                 "Pseudonocardiales",
                                 "Rhizobiales",
                                 "Rhodobacterales",
                                 "Other"))

# create a palette for taxon groups
taxa.palette<-ancom_species_table %>% 
  distinct(taxon_order) %>% 
  pull(taxon_order)
taxa.palette<-c(get_palette("Paired",n_distinct(taxa.palette)-1),"gray")

# create vector for ordering rows
# actually the order of genera that we care about here
# so need to produce a list of the genera in order, then expand it to species

row.list<-ancom_species_table %>%
  group_by(genus) %>% 
  summarise(mean=mean(`clr-abundance`)) %>% 
  arrange(mean) %>% pull(genus)

# differential abundance heatmap
diff_abund_heat <-
  ancom_species_table %>% 
  mutate(taxon=as.character(taxon)) %>% 
  left_join(sample_data(spe_data_final) %>% data.frame() %>% rownames_to_column("sample")) %>%
  # right_join(spe_data_final %>% tax_glom(., taxrank="genus") %>% microbiome::transform(.,"clr") %>% otu_table() %>%  data.frame() %>% rownames_to_column("taxon") %>% pivot_longer(cols=-"taxon",names_to="sample",values_to="genus-clr-abundance")) %>% 
  mutate(subspecies_locality=fct_recode(subspecies_locality,
                                        "Western Lowland"="gorilla",
                                        "HA Grauer's"="graueri >1000",
                                        "LA Grauer's"="graueri <1000",
                                        "Mountain"="beringei")) %>% 
  as.data.frame() %>%
  heat(., Yvar="taxon_name", fill="clr-abundance", Xvar = "sample",plot.values = F, order.cols = alt.list, order.rows = F)+
  scale_fill_gradient2(low="blue", mid="white", high="red",
                       guide = guide_colorbar(title = "CLR-normalized\nabundance", draw.ulim=FALSE,direction = "horizontal"),
                       limits=c(min(clr.pseudozeros), NA), na.value="grey") +
  # coord_cartesian(clip = "off")+
  geom_rect(data=. %>% select(genus,taxon_name,taxon_order) %>% mutate(subspecies_locality="Microbial Order"),
            aes(xmin = -0.55, xmax = -0.85, ymin = taxon_name, ymax = taxon_name,color=taxon_order),size=5) +
  scale_color_manual(values=taxa.palette)+
  # facet_grid(factor(genus,row.list)~factor(subspecies_locality,c("Microbial Order","Western Lowland","LA Grauer's","HA Grauer's","Mountain")), 
  #            scales="free",space = "free", margins = F,
  #            switch = "both") +
  facet_grid2(factor(genus,row.list)~factor(subspecies_locality,c("Microbial Order","Western Lowland","LA Grauer's","HA Grauer's","Mountain")), 
              scales="free",space = "free", margins = F,strip = strip_vanilla(size = "constant",clip = "off"),
              switch = "both") +
  guides(color = guide_legend(title = "Microbial Order",override.aes = list(fill=taxa.palette)))+
  theme_bw()+
  theme(strip.placement.y = "outside",                      # Place facet labels outside x axis labels.
        strip.background=element_rect(color = NA), 
        strip.background.y=element_rect(fill = NA), 
        plot.title = element_text(size = 20, hjust = 0.8),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle=90,size=18,hjust=1,vjust = 0.5),
        axis.ticks.x = element_line(),
        line = element_blank(),
        panel.border = element_rect(size=0),
        axis.title = element_text(size=16),
        legend.position='top',
        legend.justification = c(1,1),
        # legend.justification='left',
        legend.direction='horizontal',
        legend.text = element_text(size=17),
        legend.title = element_text(size=19),
        panel.spacing.y=unit(0,'npc'),
        panel.spacing.x=unit(c(0,0.25,0,0.25),'lines'),
        # panel.spacing.x = unit(c(-0.5,0.1,0.1,0.1),"lines"),
        # panel.spacing.y = unit(-0.5, "lines"),
        strip.text.y.left = element_text(angle = 0,size=20,hjust = 1),
        strip.text.x = element_text(color = "white",face="bold",size=20),
        plot.margin = margin(0,0,0,0, "cm"),
        text = element_text(family="calibri"))

# fill the facets with the right colour
diff_abund_heat_filled <- ggplot_gtable(ggplot_build(diff_abund_heat))
stripr <- which(grepl('strip-b', diff_abund_heat_filled$layout$name))
fills <- c("white",pal)
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', diff_abund_heat_filled$grobs[[i]]$grobs[[1]]$childrenOrder))
  diff_abund_heat_filled$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

u_taxa_order<-ancom_species_table %>% 
  mutate(taxon_order=ifelse(taxon_order=="Other","z",taxon_order)) %>% 
  distinct(taxon_order) %>% pull(taxon_order)

ancom_species_table %>% 
  left_join(sample_data(spe_data_final) %>% data.frame() %>% 
              rownames_to_column("sample")) %>% 
  mutate(taxon_order=ifelse(taxon_order=="Other","z",taxon_order)) %>% 
  ggplot()+
  labs(y="")+
  # facet_grid(~facet_label)+
  scale_color_manual(values=c(get_palette(palette = "Paired",k = 9),"gray"),
                     labels=c(u_taxa_order[1:6],u_taxa_order[8:9],"Other"))+
  # guides(color=guide_colorsteps(title="Microbial Order"))+
  theme(strip.text.x = element_text(angle = 90),
        strip.background.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(0,0,0,0, "cm"))

legend<-get_legend(diff_abund_order)

#Get taxa sums for all taxa in the final dataset
comp_ancomVStotal <- taxa_sums(spe_data_final)
comp_ancomVStotal <- as.data.frame(cbind(comp_ancomVStotal, ifelse(names(comp_ancomVStotal) %in% ancom_species$taxa_id, "Diff. abundant", "Not diff. abundant")))
colnames(comp_ancomVStotal) <- c("log_abundance", "type")
comp_ancomVStotal$log_abundance <- log(as.numeric(comp_ancomVStotal$log_abundance) + 1)

#Plot
box.p<-ggboxplot(comp_ancomVStotal, x="type", y="log_abundance", fill="grey") +
  annotate("text",size=4, x=0.8, y=10, label=paste0("n = ", sum(comp_ancomVStotal$type=="Not diff. abundant"))) +
  annotate("text",size=4, x=1.8, y=10.2, label=paste0("n = ", sum(comp_ancomVStotal$type=="Diff. abundant"))) +
  ylab("log-transformed abundance")+
  ggtitle("log-transformed abundance")+
  # annotate(geom = "text",x = 1.5,y=16,label="T test, t(123.59) = 1.24, p=0.22")+
  theme(plot.title = element_text(size = 20,hjust = 0.5),
        axis.text.y = element_text(size=16),
        axis.text.x = element_text(size=16,angle=90),
        axis.ticks.x = element_line(),
        axis.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size=14),
        plot.margin = margin(0,0,0,0, "cm"),
        text = element_text(family="calibri"))

ggarrange(as_ggplot(diff_abund_heat_filled),
          box.p+theme(plot.margin = margin(2,2,20,1, "cm")),
          widths = c(1,0.25),
          heights = c(1,0.25),
          # scale=c(2,0.1),
          # legend,
          nrow = 1)

ggplot2::ggsave(plot = as_ggplot(diff_abund_heat_filled),filename = "T3_community-level/diff_abund_heat.png",dpi=300,height = 18,width = 18)
ggplot2::ggsave(plot = last_plot(),filename = "T3_community-level/diff_abund_heat_with_box.png",dpi=300,height = 16,width = 22)

#### Figure 6 - Subspecies-associated biological processes ####
load("../Gorilla_dental_calculus_zenodo/F2_functional_stats/.RData")
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

ggplot2::ggsave(tax_func_grid, file="BP_heatmap_with_tax_contributions.png",
                height=8, width=15)

#### Plot heatmap ####
#Plot abundances
ancom_euk_abund <- otu_table(subset_taxa(euk_genus_diet_norm, genus %in% ancom_euk$genus))

# rownames(ancom_euk_abund)<-data.frame(tax_table(euk_genus_diet_norm))[data.frame(tax_table(euk_genus_diet_norm))$genus %in% ancom_euk$genus,"genus"]

#Set 0 values (the most negative value of CLR-normalized abundances) to min(clr.pseudocounts)*2 instead of NA (which is what Jaelle did). This is because the heatmap function that I use doesn't recognize NAs
clr.pseudozeros.dd <- sapply(colnames(ancom_euk_abund), function(x){min(ancom_euk_abund[,x])})
for (s in colnames(ancom_euk_abund)) {
  ancom_euk_abund[which(ancom_euk_abund[,s] == clr.pseudozeros.dd[s]),s] <- min(clr.pseudozeros.dd)*2
}

#Melt
ancom_euk_abund <- reshape2::melt(ancom_euk_abund)
colnames(ancom_euk_abund) <- c("taxa_id", "sample", "CLR-abundance")

# just update with the phyloseq metadata
ancom_euk_abund <- sample_data(spe_data_final) %>% data.frame() %>% rownames_to_column("sample") %>% right_join(ancom_euk_abund)

# add genus column
ancom_euk_abund$genus <- taxonomy_euk[,"genus"][match(ancom_euk_abund$taxa_id, rownames(taxonomy_euk))]

#Add phylum column
ancom_euk_abund$phylum_name <- taxonomy_euk[,"phylum"][match(ancom_euk_abund$genus, taxonomy_euk[,"genus"])]

#Add family columns
ancom_euk_abund$family_name <- taxonomy_euk[,"family"][match(ancom_euk_abund$genus, taxonomy_euk[,"genus"])]

# add in subspecies
ancom_euk_abund$Spec.subspecies <- factor(ancom_euk_abund$Spec.subspecies, levels = c("gorilla", "graueri", "beringei"))

#Order table based on host subspecies, phylum and family
ancom_euk_abund <- ancom_euk_abund[order(ancom_euk_abund$Spec.subspecies, ancom_euk_abund$phylum_name, ancom_euk_abund$family_name),]

#Plot
#Labeller for the facet labels
facet.labels <- list("gorilla"="western lowland", "graueri"="Grauer's", "beringei"="mountain")
facet.labeller <- function(variable,value){
  return(facet.labels[value])
}

heat_diet <- 
  ancom_euk_abund %>% 
  as.data.frame() %>%
  mutate(subspecies_locality=fct_recode(subspecies_locality,
                                        # "Microbial Order"=label.col,
                                        "Western Lowland"="gorilla",
                                        "HA Grauer's"="graueri >1000",
                                        "LA Grauer's"="graueri <1000",
                                        "Mountain"="beringei")) %>% 
  heat(., Yvar="taxa_id", fill="CLR-abundance", Xvar = "sample", order.cols = FALSE, order.rows = FALSE) +
  scale_fill_gradient2(low="blue", mid="white", high="red",
                       guide = guide_colorbar(barheight = 2, title = "CLR-normalized\nabundance", draw.ulim=FALSE),
                       limits=c(min(clr.pseudozeros.dd), NA), na.value="grey") +
  facet_grid(~factor(subspecies_locality,c("Western Lowland","LA Grauer's","HA Grauer's","Mountain")), scales="free",space="free",switch = "x")+
  # theme(legend.position = "left", plot.title=element_text(hjust = 0.5)) +
  # ggtitle("Abundances of differentially abundant gorilla dietary genera across samples")
  theme(axis.text.x = element_text(vjust = 0.5,hjust = 1),
        #     axis.ticks.x = element_line(),
        #     line = element_blank(),
        #     axis.title = element_text(size=16),
        #     legend.position='top',
        #     legend.text = element_text(size=17),
        #     legend.title = element_text(size=19),
        panel.spacing.x=unit(c(0.25,0,0.25),'lines'),
        strip.text.x = element_text(color = "white",face="bold",size=14),
        text = element_text(family="calibri"))

#### Colour the labels by phylum ####
heat_diet_ylabs <- unique(heat_diet$data$taxa_id)
heat_diet_ylabs <- cbind(heat_diet_ylabs, tax_table(euk_genus_diet)[match(heat_diet_ylabs, taxa_names(euk_genus_diet)),c("phylum","genus")])
heat_diet_ylabs <- as.data.frame(heat_diet_ylabs)
heat_diet_ylabs$phylum <- factor(heat_diet_ylabs$phylum)

#Also add family infomation (will be needed later)
heat_diet_ylabs$family <- tax_table(euk_genus_diet)[match(rownames(heat_diet_ylabs), taxa_names(euk_genus_diet)),"family"]

#Save numbers of genera per family
heat_diet_fam_freqs <- table(heat_diet_ylabs[,4])
heat_diet_fam_freqs <- as.data.frame(heat_diet_fam_freqs)
heat_diet_fam_freqs <- heat_diet_fam_freqs[match(unique(heat_diet_ylabs[,4]), heat_diet_fam_freqs$Var1),]

#Subset diet_palette.g to get colours for the taxa in this heatmap
diet_palette_full <- c("steelblue1", "sienna3", "khaki3", "olivedrab3")
names(diet_palette_full) <- levels(heat_diet_ylabs$genus) # Give every color an appropriate name

diet_palette_full.g <- heat_diet_ylabs
diet_palette_full.g$colour <- diet_palette_full[match(diet_palette_full.g$phylum, names(diet_palette_full))]


diet_palette.g <- diet_palette_full.g[match(rownames(heat_diet_ylabs), rownames(diet_palette_full.g)),]

#Also draw lines to show genera that belong to the same family
#From the sorted table with the number of genera per family, get the cumulative sum
#to determine where every line will be along the y axis
heat_diet_fam_freqs$Freq <- cumsum(heat_diet_fam_freqs$Freq)

#Update plot
heat_diet <-
  heat_diet + theme(axis.text.y = element_text(colour = diet_palette.g,size=12,hjust=1,vjust = 0.5),
                    legend.position = "top") +
  geom_hline(yintercept=heat_diet_fam_freqs$Freq + 0.5, size=0.5, colour="white")

# fill the facets with the right colour
pal<-c(get_palette(palette = "Set2",k=3)[1],
       get_palette(palette = "Set2",k=3)[2],
       get_palette(palette = "Set2",k=3)[2],
       get_palette(palette = "Set2",k=3)[3])
heat_diet_filled <- ggplot_gtable(ggplot_build(heat_diet))
stripr <- which(grepl('strip-b', heat_diet_filled$layout$name))
fills <- pal
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', heat_diet_filled$grobs[[i]]$grobs[[1]]$childrenOrder))
  heat_diet_filled$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

#### Create a sidebar ####
#based on the families in the studies
ref_sidebar <- heat_diet_ylabs[,c("genus", "family")]
ref_sidebar <- cbind(ref_sidebar, rownames_to_column(diet_ref_fams)[match(ref_sidebar$family, rownames_to_column(diet_ref_fams)$rowname),])

rownames() <- NULL

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
  theme(axis.text.y = element_blank(), 
        axis.text.x = element_text(size=10,hjust=1,vjust = 0.5),
        axis.ticks.y = element_blank(), 
        plot.margin = margin(1.6,0,1.2,0, "cm"),
        plot.title=element_text(size=10, hjust = 0.5)) +
  guides(fill=FALSE) + 
  scale_fill_gradient2(low = "grey", mid="palegreen2", high = "palegreen4", midpoint = 1) +
  geom_hline(yintercept=heat_diet_fam_freqs$Freq + 0.5, size=0.5, colour="white") +
  ggtitle("Mentioned\nin literature")

#Combine the two plots
heat_diet_complete <- plot_grid(as_ggplot(heat_diet_filled), heat_sidebar, 
                                ncol = 2, 
                                align = "v",
                                axis="tb",
                                rel_widths = c(1, 0.1))
heat_diet_complete

ggsave(heat_diet_complete, file="D3_diet_stats/heat_diet_complete.png", device="png", height=7, width=10)

#### Figure S1 - read count histogram ####
load("T3_community-level/.RData")

spe_data@sam_data$readcount_log <- log(spe_data@sam_data$readcount.m.before.Kraken)
#Check distribution of read
ggplot2::ggsave(file="read_count_raw_hist.png",
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
ggplot2::ggsave(file="oral_proporton_hist.png", 
                gghistogram(oral_proportion, "oral_proportion_log", fill="Sample.type", position="stack") +
                  #Include vertical line with the cutoff of 0.03 (3%)
                  geom_vline(xintercept=log(0.03), size=1, colour="red"))

#### Figure S3 - oral/contaminant proportion ####
decontam_prop_core_micr <- 
  heat(decontam_test, Xvar="uppsala_threshold", Yvar="jena_threshold", fill="prop_core_micr",
       order.rows=FALSE, order.cols=FALSE) +
  scale_fill_gradient2(guide = guide_colorbar(barheight = 10, title="core microbiome\ntaxa proportion")) +
  ylab("jena_threshold") + xlab("uppsala_threshold")

ggplot2::ggsave(decontam_prop_core_micr, file="decontam_prop_core_micr.png", device="png")

#### Figure S4 - taxonomic ordinations with duplicate samples ####
jaccard_with_duplicates <- jaccard_with_duplicates + theme_light() + geom_point(size=3) + theme(legend.position="none") + scale_colour_brewer(palette="Pastel1", name="Dataset", labels=c("Uppsala"="Uppsala", "Jena"="Fellows Yates et al. 2021"))

clr_with_duplicates <- clr_with_duplicates + theme_light() + geom_point(size=3) + scale_colour_brewer(palette="Pastel1", name="Dataset", labels=c("Uppsala"="Uppsala", "Jena"="Fellows Yates et al. 2021"))

ordinations_with_duplicates <- grid.arrange(jaccard_with_duplicates, clr_with_duplicates, ncol=2, widths=c(2.4, 3.2))

ggplot2::ggsave(ordinations_with_duplicates, file="ordinations_with_duplicates.png", device="png", height=7, width=14)

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

ggplot2::ggsave(ordinations_with_altitude, file="ordinations_with_altitude.png", device="png", height=7, width=14)


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
ggplot2::ggsave(
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
contaminant_genomes <- read.table("RD2_mapping/contaminant_genomes_list_recent.txt", sep="\t", comment.char="", header=F)
noncontaminant_genomes <- read.table("RD2_mapping/noncontaminant_genomes_list_recent.txt", sep="\t", comment.char="", header=F)

colnames(contaminant_genomes) <- c("taxon_ID", "RefSeq", "Assembly_level", "Genome_repres", "Date", "ftp" )
colnames(noncontaminant_genomes) <- c("taxon_ID", "RefSeq", "Assembly_level", "Genome_repres", "Date", "ftp" )

print("How many contaminant genomes?")
nrow(contaminant_genomes)

print("How many non-contaminant genomes?")
nrow(noncontaminant_genomes)

ggplot2::ggsave(
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

load(".RData")

sample_data(GO_BP_phyloseq)$Spec.subspecies <- factor(sample_data(GO_BP_phyloseq)$Spec.subspecies, levels=c("gorilla", "graueri", "beringei"))

func_ordination <- plot_ordination(GO_BP_phyloseq, GO_eucl, color="Spec.subspecies", shape = "Seq.centre", title="PCoA on euclidean distances of functional profiles") + 
  theme_light() + geom_point(size=4) + theme(legend.text = element_text(size=14), legend.title = element_text(size=14),
                                             legend.position = "right", plot.title = element_text(size=15)) +
  scale_colour_brewer(palette="Set2", name="Host subspecies", labels = c("beringei"="Mountain", "gorilla"="Western lowland", "graueri"="Grauer's")) +
  scale_shape_discrete(name="Dataset", labels = c("Jena"="Fellows Yates et al. 2021"))

ggplot2::ggsave(func_ordination, file="func_ordination.png", device="png", height=7, width=9)