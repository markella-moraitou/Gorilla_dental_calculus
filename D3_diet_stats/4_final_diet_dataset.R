#Start logging
sink(file = "log4_final_diet_dataset.txt")

#load packages
library(phyloseq)
library(microbiome)
library(reshape2)
library(tibble)
library(ggpubr)
library(cowplot)

#Diet analysis - Script 4
#Get final diet dataset. Keep only a priori know dietary taxa

load(".RData")

euk_genus_diet<-readRDS("D3_diet_stats/euk_genus_diet.rds")
euk_genus_diet_norm<-readRDS("D3_diet_stats/euk_genus_diet_norm.rds")
diet_ref_fams<-readRDS("D3_diet_stats/diet_ref_fams.rds")
euk_genus_decontam2<-readRDS("D3_diet_stats/euk_genus_decontam2.rds")
diet_ref_gen <- readRDS("D3_diet_stats/diet_ref_gen.rds")

print("Number of families in the reference db")
nrow(diet_ref_fams) 
print("How many families are found in each gorilla subspecies?")
colSums(diet_ref_fams)

print("Which families ae found in all gorilla subspecies?")
rownames(diet_ref_fams)[which(rowSums(diet_ref_fams)==3)]

print("Which are found in only one subspecies?")
rownames(diet_ref_fams)[which(rowSums(diet_ref_fams)==1)]
print("How many families are unique to each subspecies?")
colSums(diet_ref_fams[which(rowSums(diet_ref_fams)==1),])

#### Extract genera based on a priori knowledge of diet ####
#Extract taxa that belong to these 85 families (also keep only true samples)
euk_genus_diet <- subset_taxa(euk_genus_decontam2, tax_table(euk_genus_decontam2)[,6] %in% rownames(diet_ref_fams))

#Also, retain only samples that were retained in the taxonomic analysis
euk_genus_diet <- subset_samples(euk_genus_diet, sample_names(euk_genus_diet) %in% sample_names(spe_data_final))

#Normalize
euk_genus_diet_norm <- euk_genus_diet
otu_table(euk_genus_diet_norm) <- otu_table(microbiome::transform(otu_table(euk_genus_diet_norm), transform = "clr"),
                                            taxa_are_rows = TRUE)

#### Plot heatmap ####

#Get table with normalized abundances
euk_genus_diet_abund <- otu_table(euk_genus_diet_norm)

#Set 0 values (the most negative value of CLR-normalized abundances) to min(clr.pseudocounts)*2 instead of NA (which is what Jaelle did). This is because the heatmap function that I use doesn't recognize NAs
clr.pseudozeros.d <- sapply(colnames(euk_genus_diet_abund), function(x){min(euk_genus_diet_abund[,x])})
for (s in colnames(euk_genus_diet_abund)) {
  euk_genus_diet_abund[which(euk_genus_diet_abund[,s] == clr.pseudozeros.d[s]),s] <- min(clr.pseudozeros.d)*2
}


#Melt
euk_genus_diet_abund <- reshape2::melt(euk_genus_diet_abund)
colnames(euk_genus_diet_abund) <- c("genus", "sample", "CLR-abundance")

#Add host subspecies column
euk_genus_diet_abund$host_subspecies <- metadata$Spec.subspecies[match(euk_genus_diet_abund$sample, rownames(metadata))]

#Add phylum column
euk_genus_diet_abund$phylum_name <- taxonomy_euk[,3][match(euk_genus_diet_abund$genus, taxonomy_euk[,7])]

#Add family columns
euk_genus_diet_abund$family_name <- taxonomy_euk[,6][match(euk_genus_diet_abund$genus, taxonomy_euk[,7])]

#Order table based on host subspecies, phylum and family
euk_genus_diet_abund$host_subspecies <- factor(euk_genus_diet_abund$host_subspecies, levels = c("gorilla", "graueri", "beringei"))
euk_genus_diet_abund <- euk_genus_diet_abund[order(euk_genus_diet_abund$host_subspecies, euk_genus_diet_abund$phylum_name, euk_genus_diet_abund$family_name),]

heat_diet_full <- heat(euk_genus_diet_abund, Yvar="genus", fill="CLR-abundance", Xvar = "sample", order.cols = FALSE, order.rows=FALSE) +
  scale_fill_gradient2(low="blue", mid="white", high="red",
    guide = guide_colorbar(barheight = 2, title = "CLR-normalized\nabundance", draw.ulim=FALSE),
    limits=c(min(clr.pseudozeros.d), NA), na.value="grey") +
  facet_grid(~host_subspecies, scales="free", labeller=facet.labeller) +
  theme(legend.position = "left", axis.text.y = element_text()) +
  ggtitle("Abundances of known gorilla dietary genera across samples")

#### Colour the axis ticks by phylum ####
heat_diet_full_ylabs <- levels(heat_diet_full$data$YYYY)
heat_diet_full_ylabs <- cbind(heat_diet_full_ylabs, tax_table(euk_genus_diet)[match(heat_diet_full_ylabs, taxa_names(euk_genus_diet)),3])
heat_diet_full_ylabs <- as.data.frame(heat_diet_full_ylabs)
colnames(heat_diet_full_ylabs) <- c("genus", "phylum")
heat_diet_full_ylabs$phylum <- factor(heat_diet_full_ylabs$phylum)

#Also add family infomation (will be needed later)
heat_diet_full_ylabs$family <- tax_table(euk_genus_diet)[match(heat_diet_full_ylabs$genus, taxa_names(euk_genus_diet)),6]

diet_palette_full <- rev(c("#009E93","#56B4E9","#D55E00","#CC79A7"))

names(diet_palette_full) <- levels(heat_diet_full_ylabs$phylum) # Give every color an appropriate name

diet_palette_full.g <- heat_diet_full_ylabs
diet_palette_full.g$colour <- diet_palette_full[match(diet_palette_full.g$phylum, names(diet_palette_full))]

rownames(diet_palette_full.g) <- diet_palette_full.g[,1]
diet_palette_full.g <- as.data.frame(diet_palette_full.g[,-c(1,2,3)])
rownames(diet_palette_full.g) <- heat_diet_full_ylabs$genus
colnames(diet_palette_full.g) <- NULL
diet_palette_full.g <- t(as.matrix(diet_palette_full.g))

#Update plot
heat_diet_full <-
  heat_diet_full + theme(axis.ticks.y = element_line(colour = diet_palette_full.g),
                       axis.ticks.length.y = unit(0.5, "cm"),
                       axis.text.y = element_blank(), legend.position = "top")

#### Create a sidebar ####
#based on a priori knowledge of these taxa
ref_sidebar_full <- heat_diet_full_ylabs[,c("genus", "family")]
#Add info about family mentions in literature
ref_sidebar_full <- cbind(ref_sidebar_full, rownames_to_column(diet_ref_fams)[match(ref_sidebar_full$family, rownames_to_column(diet_ref_fams)$rowname),])
rownames(ref_sidebar_full) <- NULL
#Add info about genus mentions in literature
ref_sidebar_full <- cbind(ref_sidebar_full, rownames_to_column(diet_ref_gen)[match(ref_sidebar_full$genus, rownames_to_column(diet_ref_gen)$rowname),])
#A bunch of NAs because of the genera that are not mentioned in the references
#turn them into FALSE (because the species is absent)
ref_sidebar_full[is.na(ref_sidebar_full)] <- FALSE

#Remove unnecessary colums
ref_sidebar_full <- ref_sidebar_full[,-which(colnames(ref_sidebar_full) %in% c("family", "rowname")),]

#Replace FALSE -> 0, TRUE -> 1
ref_sidebar_full[ref_sidebar_full==FALSE] <- 0

#Add columns referring to the same diet. 
#This way if the exact genus is mentioned in the literature it takes a value of 2
#But if it is only the family, it takes a value of 1

ref_sidebar_full$western <- ref_sidebar_full$western + ref_sidebar_full$western.1
ref_sidebar_full$grauers <- ref_sidebar_full$grauers + ref_sidebar_full$grauers.1
ref_sidebar_full$mountain <- ref_sidebar_full$mountain + ref_sidebar_full$mountain.1

#Remove extra columns
ref_sidebar_full <- ref_sidebar_full[, -which(grepl(".1", colnames(ref_sidebar_full)))]
rownames(ref_sidebar_full) <- NULL

#Melt
ref_sidebar_full <- reshape::melt(ref_sidebar_full)

colnames(ref_sidebar_full) <- c("diet_genus", "host_subspecies", "presence/absence")

#Set presence/absence as numeric
ref_sidebar_full$`presence/absence` <- as.numeric(ref_sidebar_full$`presence/absence`)

#Plot heatmap
heat_sidebar_full <- heat(ref_sidebar_full, Yvar="diet_genus", fill="presence/absence", Xvar = "host_subspecies", order.cols = FALSE, order.rows = FALSE) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_text(angle = 90, vjust = -0.2, size=10),
                                                                                                                                                                  axis.title=element_text(size=14), plot.title=element_text(size=10, hjust = 0.5)) +
  guides(fill=FALSE) + scale_fill_gradient2(low = "grey", mid = "palegreen2", high = "palegreen4", midpoint=1) +
  ggtitle("Mentioned\nin literature")

#Combine the two plots
heat_diet_full_complete <- plot_grid(heat_diet_full, heat_sidebar_full, align = "hv", ncol = 2, rel_widths = c(44, 6),
                                     axis="tb")

#Save plot
ggsave(heat_diet_full_complete, device="png", width=9, height=7,
       filename = "D3_diet_stats/heat_diet_full_complete.png")


#### Summarize at the family level ####

#First remove 'clade' info from taxonomy table (it produces wrong results)
euk_genus_diet@tax_table[,2] <- NA

euk_fam_diet <- tax_glom(euk_genus_diet, taxrank = "family", NArm = FALSE, bad_empty = NA)

#Rename taxa from genera to families
taxa_names(euk_fam_diet) <- make.names(tax_table(euk_fam_diet)[,6], unique=TRUE)

#Just keep non-NA classifications
euk_fam_diet <- subset_taxa(euk_fam_diet, !(is.na(tax_table(euk_fam_diet)[,6])))

#get some metadata
sample_data(euk_fam_diet)$family_richness <- sapply(row.names(sample_data(euk_fam_diet)), function(x) { #Families per sample
  estimate_richness(euk_fam_diet, measures = "Observed")[x,]})

#normalize abundances
euk_fam_diet_norm <- euk_fam_diet
otu_table(euk_fam_diet_norm) <- otu_table(microbiome::transform(otu_table(euk_fam_diet_norm), "clr"), taxa_are_rows = TRUE)

#### Plot heatmap ####

#Get table with normalized abundances
euk_fam_diet_abund <- otu_table(euk_fam_diet_norm)

#Melt
euk_fam_diet_abund <- reshape2::melt(euk_fam_diet_abund)
colnames(euk_fam_diet_abund) <- c("family", "sample", "CLR-abundance")

#Add host subspecies column
euk_fam_diet_abund$host_subspecies <- metadata$Spec.subspecies[match(euk_fam_diet_abund$sample, rownames(metadata))]

#Add phylum column
euk_fam_diet_abund$phylum_name <- taxonomy_euk[,3][match(euk_fam_diet_abund$family, taxonomy_euk[,6])]

#Order table based on host subspecies and phylum
euk_fam_diet_abund$host_subspecies <- factor(euk_fam_diet_abund$host_subspecies, levels = c("gorilla", "graueri", "beringei"))
euk_fam_diet_abund <- euk_fam_diet_abund[order(euk_fam_diet_abund$host_subspecies, euk_fam_diet_abund$phylum_name),]

# update altitude in metadata
alt.fix<-read_tsv("T3_community-level/new_altitude.csv") %>% 
  right_join(sample_data(spe_data_final) %>% data.frame() %>% rownames_to_column("Sample_ID")) %>%
  mutate(subspecies_locality=ifelse(Spec.subspecies=="graueri",
                                    ifelse(`Approximate altitude`>1000,"graueri >1000","graueri <1000"),
                                    as.character(subspecies_locality))) %>% 
  pull(subspecies_locality)

sample_data(spe_data_final)$subspecies_locality<-alt.fix

heat_diet_fam_full <- 
  euk_fam_diet_abund %>% 
  rownames_to_column("taxa_id") %>% 
  mutate(taxa_id=as.character(taxa_id)) %>%
  left_join(sample_data(spe_data_final) %>% data.frame() %>% rownames_to_column("sample")) %>%
  mutate(subspecies_locality=fct_recode(subspecies_locality,
                                        "Western Lowland"="gorilla",
                                        "HA Grauer's"="graueri >1000",
                                        "LA Grauer's"="graueri <1000",
                                        "Mountain"="beringei")) %>% 
  as.data.frame() %>%
  heat(., Yvar="family", fill="CLR-abundance", Xvar = "sample", order.cols = FALSE, order.rows=FALSE) +
  scale_fill_gradient2(low="blue", mid="white", high="red",
                       guide = guide_colorbar(barheight = 2, title = "CLR-normalized \nabundance", draw.ulim=FALSE),
                       limits=c(min(clr.pseudozeros.d), NA), na.value="grey") +
  facet_grid(~factor(subspecies_locality,c("Western Lowland","LA Grauer's","HA Grauer's","Mountain")), scales="free",space="free",switch = "both")+
  theme(
    axis.text.x = element_text(hjust=1,vjust = 0.5),
        #     axis.ticks.x = element_line(),
        #     line = element_blank(),
        #     axis.title = element_text(size=16),
        #     legend.position='top',
        #     legend.text = element_text(size=17),
        #     legend.title = element_text(size=19),
        panel.spacing.x=unit(c(0.25,0,0.25),'lines'),
        strip.text.x = element_text(color = "white",face="bold",size=12),
        text = element_text(family="calibri"))

#### Colour the axis ticks by phylum ####
heat_diet_fam_full_ylabs <- levels(heat_diet_fam_full$data$YYYY)
heat_diet_fam_full_ylabs <- cbind(heat_diet_fam_full_ylabs, tax_table(euk_fam_diet)[match(heat_diet_fam_full_ylabs, taxa_names(euk_fam_diet)),3])
heat_diet_fam_full_ylabs <- as.data.frame(heat_diet_fam_full_ylabs)
colnames(heat_diet_fam_full_ylabs) <- c("family", "phylum")
heat_diet_fam_full_ylabs$phylum <- factor(heat_diet_fam_full_ylabs$phylum)

diet_palette_full.f <- heat_diet_fam_full_ylabs

diet_palette_full.f$colour <- diet_palette_full[match(diet_palette_full.f$phylum, names(diet_palette_full))]

rownames(diet_palette_full.f) <- diet_palette_full.f[,1]
diet_palette_full.f <- as.data.frame(diet_palette_full.f[,-c(1,2)])
rownames(diet_palette_full.f) <- heat_diet_fam_full_ylabs$family
colnames(diet_palette_full.f) <- NULL
diet_palette_full.f <- t(as.matrix(diet_palette_full.f))

#Update plot
heat_diet_fam_full <-
  heat_diet_fam_full + theme(axis.text.y = element_text(colour = diet_palette_full.f, size=10), 
                                    legend.position = "top",
                                    axis.text.x = element_text(angle = 90, vjust = -0.2, size=10), 
                                    strip.text.x=element_text(colour="white", face="bold", size=14)) +
  guides(color=guide_legend(title = "test"))

# fill the facets with the right colour
pal<-c(get_palette(palette = "Set2",k=3)[1],
       get_palette(palette = "Set2",k=3)[2],
       get_palette(palette = "Set2",k=3)[2],
       get_palette(palette = "Set2",k=3)[3])
heat_diet_fam_full_filled <- ggplot_gtable(ggplot_build(heat_diet_fam_full))
stripr <- which(grepl('strip-b', heat_diet_fam_full_filled$layout$name))
fills <- pal
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', heat_diet_fam_full_filled$grobs[[i]]$grobs[[1]]$childrenOrder))
  heat_diet_fam_full_filled$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}


#### Create a sidebar ####
#based on a priori knowledge of these taxa
ref_sidebar_full.f <- heat_diet_fam_full_ylabs[,"family"]
#Add info about family mentions in literature
ref_sidebar_full.f <- cbind(ref_sidebar_full.f, rownames_to_column(diet_ref_fams)[match(ref_sidebar_full.f, rownames_to_column(diet_ref_fams)$rowname),])
rownames(ref_sidebar_full.f) <- NULL

#Remove unnecessary colums
ref_sidebar_full.f <- ref_sidebar_full.f[,-which(colnames(ref_sidebar_full.f)=="rowname")]

colnames(ref_sidebar_full.f)[1] <- "family"

#Remove NAs
ref_sidebar_full.f <- ref_sidebar_full.f[which(!(is.na(ref_sidebar_full.f$western))),]

#Replace FALSE -> 0, TRUE -> 1
ref_sidebar_full.f[ref_sidebar_full.f==FALSE] <- 0

#Melt
ref_sidebar_full.f <- reshape::melt(ref_sidebar_full.f)

colnames(ref_sidebar_full.f) <- c("diet_family", "host_subspecies", "presence/absence")

#Set presence/absence as numeric
ref_sidebar_full.f$`presence/absence` <- as.numeric(ref_sidebar_full.f$`presence/absence`)

#Plot heatmap
heat_sidebar_full.f <- heat(ref_sidebar_full.f, Yvar="diet_family", fill="presence/absence", Xvar = "host_subspecies",
                            order.cols = FALSE, order.rows = FALSE) +
  theme(axis.text.y = element_blank(), 
        axis.text.x = element_text(size=12,hjust=1,vjust = 0.5),
        axis.ticks.y = element_blank(), 
        plot.margin = margin(1.58,0,0.75,0, "cm"),
        plot.title=element_text(size=10, hjust = 0.5)) +
  guides(fill=FALSE) + scale_fill_gradient2(low = "grey", mid = "palegreen2", high = "palegreen4", midpoint=1) +
  ggtitle("Mentioned\nin literature")

#Combine the two plots
heat_diet_fam_full_complete <- plot_grid(as_ggplot(heat_diet_fam_full_filled), heat_sidebar_full.f, 
                                         ncol = 2, 
                                         align = "v",
                                         axis="tb",
                                         rel_widths = c(1, 0.1))

#Save plot
ggsave(heat_diet_fam_full_complete, device="png", height=10, width=10,
       filename = "D3_diet_stats/heat_diet_fam_full_complete_family.png")


#Compare average abundance of all taxa vs diet taxa

#Get the taxa sums of all plants
comp_euk_abund <- taxa_sums(subset_taxa(euk_genus_decontam2, tax_table(euk_genus_decontam2)[,1] == "Eukaryota" &
                                          !(tax_table(euk_genus_decontam2)[,3] %in% "Mammalia")))

comp_euk_abund <- cbind(comp_euk_abund, rep("Total", length(comp_euk_abund)))

#Get the taxa sums of dietary plants
diet_euk_abund <- taxa_sums(euk_genus_diet)
diet_euk_abund <- cbind(diet_euk_abund, rep("Dietary", length(diet_euk_abund)))

comp_euk_abund <- as.data.frame(rbind(comp_euk_abund, diet_euk_abund))

rm(diet_euk_abund)

colnames(comp_euk_abund) <- c("taxon_sum_log", "taxon_type")
comp_euk_abund$taxon_sum <- log(as.numeric(comp_euk_abund$taxon_sum))

gghistogram(comp_euk_abund, x="taxon_sum", fill="taxon_type")
ggboxplot(comp_euk_abund, x="taxon_type", y = "taxon_sum")

# save your work!
saveRDS(euk_genus_diet,"D3_diet_stats/euk_genus_diet.rds")
saveRDS(euk_genus_diet_norm,"D3_diet_stats/euk_genus_diet_norm.rds")

sessionInfo()

#Start logging
sink(file = NULL)

save.image()