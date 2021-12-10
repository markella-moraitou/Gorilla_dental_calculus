#Start logging
sink(file = "F2_functional_stats/log11_altitude.txt")

#load packages
library(phyloseq)
library(stringr)
library(dplyr)
library(readxl)
library(reshape2)
library(tibble)
library(ggplot2)
library(tidyr)
library(vegan)
library(EcolUtils)
library(ade4)
library(cowplot)

load("F2_functional_stats/.RData")

#Check effect of altitude - Running the full script also requires the output of HUMAnN2

#Load location info

#Get functional analysis phyloseq
GO_BP_phyloseq <- readRDS("F2_functional_stats/GO_BP_phyloseq")
GO_eucl <- readRDS("F2_functional_stats/GO_eucl")

#### Grauer's grouped altitude ####
spe_data_final@sam_data$subspecies_locality <-
  #For Western Lowland and Mountain gorillas, use the subspecies classification without any grouping
  ifelse(spe_data_final@sam_data$Spec.subspecies %in% c("gorilla", "beringei"), as.character(spe_data_final@sam_data$Spec.subspecies),
    #For Grauer's group by altitude (less vs greater than 1000)
    ifelse(as.numeric(spe_data_final@sam_data$Approximate.altitude) < 1000, "graueri <1000", "graueri >1000"))

#Turn into factor
sample_data(spe_data_final)$subspecies_locality <-
  factor(sample_data(spe_data_final)$subspecies_locality, 
         levels=c("gorilla", "graueri <1000", "graueri >1000", "beringei"))

sample_data(spe_data_final)$Spec.subspecies <- factor(sample_data(spe_data_final)$Spec.subspecies, levels = c("gorilla", "graueri", "beringei"))

# read in new altitude data
# for some reason my csv's are tab delim
new.alt<-read.csv("T3_community-level/new_altitude.csv",header=T,sep = "\t") %>% column_to_rownames("Sample_ID")

# which samples have changed?
alt.table<-sample_data(spe_data_final) %>% data.frame() %>% 
  select(Spec.subspecies,Approximate.altitude,subspecies_locality) %>% 
  mutate(Approximate.altitude=as.integer(Approximate.altitude)) %>% 
  rownames_to_column("Sample_ID") %>% 
  left_join(new.alt %>% rownames_to_column("Sample_ID"),by="Sample_ID",suffix = c(".old",".new")) %>% 
  rename("subspecies_locality.old"="subspecies_locality") %>% 
  mutate(subspecies_locality.new = ifelse(Spec.subspecies %in% c("gorilla", "beringei"), as.character(Spec.subspecies),
                                      ifelse(as.numeric(Approximate.altitude.new) < 1000, "graueri <1000", "graueri >1000")),
         subspecies_locality.new = factor(subspecies_locality.new, levels=c("gorilla", "graueri <1000", "graueri >1000", "beringei")))
alt.table %>% 
  filter(Approximate.altitude.old != Approximate.altitude.new) %>% 
  select(Sample_ID,Spec.subspecies,Approximate.altitude.old,Approximate.altitude.new,subspecies_locality.old,subspecies_locality.new) %>% 
  arrange(Spec.subspecies) 
alt.table %>% knitr::kable()

# replace old altitude data with new data in phyloseq object 
spe_data_final@sam_data$Approximate.altitude <- new.alt[rownames(spe_data_final@sam_data),"Approximate.altitude"]
spe_data_final@sam_data$Approximate.latitude <- new.alt[rownames(spe_data_final@sam_data),"Approximate.latitude"]
spe_data_final@sam_data$Approximate.longtitude <- new.alt[rownames(spe_data_final@sam_data),"Approximate.longtitude"]
spe_data_final@sam_data$subspecies_locality <- alt.table[alt.table$Sample_ID %in% rownames(spe_data_final@sam_data),"subspecies_locality.new"]
  
GO_BP_phyloseq@sam_data$Approximate.altitude <- new.alt[rownames(GO_BP_phyloseq@sam_data),"Approximate.altitude"]
GO_BP_phyloseq@sam_data$Approximate.latitude <- new.alt[rownames(GO_BP_phyloseq@sam_data),"Approximate.latitude"]
GO_BP_phyloseq@sam_data$Approximate.longtitude <- new.alt[rownames(GO_BP_phyloseq@sam_data),"Approximate.longtitude"]
GO_BP_phyloseq@sam_data$subspecies_locality <- alt.table[alt.table$Sample_ID %in% rownames(GO_BP_phyloseq@sam_data),"subspecies_locality.new"]

#Get the same info in the normalized phyloseq and in the functional phyloseq
sample_data(spe_data_final_norm) <- sample_data(spe_data_final)
sample_data(GO_BP_phyloseq) <- sample_data(spe_data_final)

#### Plot ordination ####
aitch<-plot_ordination(spe_data_final_norm, clr_5, color="Spec.subspecies", shape = "subspecies_locality", title="Aitchison") + 
  geom_point(size=3) +
  scale_colour_brewer(palette="Set2", name="Host subspecies", labels = c("beringei"="Mountain", "gorilla"="Western", "graueri"="Grauer's")) + 
  scale_shape_manual(values=c(16, 15, 17, 16), name="Grauer's altitude", labels = c("beringei"="Mountain", "gorilla"="Western", "graueri <1000"="Grauer's <1000m", "graueri >1000"="Grauer's >1000m")) +
  theme_light() + theme(axis.text = element_text(size=12),axis.title = element_text(size = 14),legend.text = element_text(size=14), legend.title = element_text(size=14),
                        legend.position = "right", plot.title = element_text(size=15,hjust = 0.5))
  
jacc<-plot_ordination(spe_data_final, jaccard_5, color="Spec.subspecies", shape = "subspecies_locality", title="Jaccard") + 
  geom_point(size=3) +
  scale_colour_brewer(palette="Set2", name="Host subspecies", labels = c("beringei"="Mountain", "gorilla"="Western", "graueri"="Grauer's")) + 
  scale_shape_manual(values=c(16, 15, 17, 16), name="Grauer's altitude", labels = c("beringei"="Mountain", "gorilla"="Western", "graueri <1000"="Grauer's <1000m", "graueri >1000"="Grauer's >1000m"))+
  theme_light() + 
  theme(axis.text = element_text(size=12),axis.title = element_text(size = 14),legend.text = element_text(size=14), legend.title = element_text(size=14),
      legend.position = "none", plot.title = element_text(size=15,hjust = 0.5))

ggpubr::ggarrange(jacc,aitch+theme(legend.position = "none"),nrow = 1,labels = c("a)","b)"),legend.grob = ggpubr::get_legend(aitch),legend = "right") 
ggsave(plot = last_plot(),
       file="T3_community-level/comb_ord_altitude.png",
       device="png",width = 12,height = 8,units = "in",dpi = 300)


plot_ordination(GO_BP_phyloseq, GO_eucl, color="Spec.subspecies", shape = "subspecies_locality", title="Jaccard") + 
  theme_light() + geom_point(size=3) + theme(legend.text = element_text(size=12), legend.title = element_text(size=12),
                                             legend.position = "right", plot.title = element_text(size=15))  +
  scale_colour_brewer(palette="Set2", name="Host subspecies", labels = c("beringei"="Mountain", "gorilla"="Western", "graueri"="Grauer's")) + 
  scale_shape_manual(values=c(16, 15, 17, 16), name="Grauer's altitude", labels = c("beringei"="Mountain", "gorilla"="Western", "graueri < 1000"="Grauer's <1000m", "graueri > 1000"="Grauer's >1000m"))
ggsave(plot = last_plot(),
  file="T3_community-level/GO_abund_with_altitude.png", device="png")

#### PERMANOVA ####
print("PERMANOVA, including subspecies locality")
read_count <- sample_data(spe_data_final)$readcount.m.before.Kraken
seq_centre <- sample_data(spe_data_final)$Seq.centre
spec.subspecies <- sample_data(spe_data_final)$Spec.subspecies
subspecies_locality <- sample_data(spe_data_final)$subspecies_locality

model2_jaccard_w_locality <- adonis(t(otu_table(spe_data_final)) ~ read_count + seq_centre + spec.subspecies + subspecies_locality, permutations = 10000, method = "jaccard")
model2_jaccard_w_locality$aov.tab %>%  data.frame() %>% mutate(across(where(is.numeric)),round(.,digits = 2)) %>% knitr::kable(format = "markdown")

model2_clr_w_locality <- adonis(t(otu_table(spe_data_final_norm)) ~ read_count + seq_centre + spec.subspecies + subspecies_locality, permutations = 10000, method = "euclidean")
model2_clr_w_locality$aov.tab %>%  data.frame() %>% mutate(across(where(is.numeric)),round(.,digits = 2)) %>% knitr::kable(format = "markdown")

read_count <- sample_data(GO_BP_phyloseq)$readcount.m.before.Kraken
seq_centre <- sample_data(GO_BP_phyloseq)$Seq.centre
spec.subspecies <- sample_data(GO_BP_phyloseq)$Spec.subspecies
subspecies_locality <- sample_data(GO_BP_phyloseq)$subspecies_locality

GO_BP_model_w_locality <- adonis(t(otu_table(GO_BP_phyloseq)) ~ read_count + seq_centre + spec.subspecies + subspecies_locality, permutations = 10000, method = "euclidean")
GO_BP_model_w_locality

#Pairwise PERMANOVA
print("Pairwise PERMANOVA")

vegdist(t(otu_table(spe_data_final_norm)), method="euclidean") %>%
  adonis.pair(Factor=subspecies_locality)

#### Get distance matrices ####
#Turn altitude, latitude, longitude into numeric values
spe_data_final@sam_data$Approximate.altitude <- as.numeric(spe_data_final@sam_data$Approximate.altitude)
spe_data_final@sam_data$Approximate.latitude <- as.numeric(spe_data_final@sam_data$Approximate.latitude)
spe_data_final@sam_data$Approximate.longitude <- as.numeric(spe_data_final@sam_data$Approximate.longitude)

spe_data_final_norm@sam_data <- spe_data_final@sam_data

GO_BP_phyloseq@sam_data$Approximate.altitude <- as.numeric(spe_data_final@sam_data$Approximate.altitude)
GO_BP_phyloseq@sam_data$Approximate.latitude <- as.numeric(spe_data_final@sam_data$Approximate.latitude)
GO_BP_phyloseq@sam_data$Approximate.longitude <- as.numeric(spe_data_final@sam_data$Approximate.longitude)

#microbial composition
jaccard_dist <- vegdist(t(spe_data_final@otu_table), "jaccard")
aitchison_dist <- vegdist(t(spe_data_final_norm@otu_table), "euclidean")

#Functional composition
func_dist <- vegdist(t(GO_BP_phyloseq@otu_table), "euclidean")

#altitude
altitude_dist <- vegdist(spe_data_final@sam_data$Approximate.altitude, na.rm=TRUE, "euclidean")

#geographic distance
geo_dist <- vegdist(cbind(spe_data_final@sam_data$Approximate.latitude, 
                          spe_data_final@sam_data$Approximate.longitude), "euclidean")


#Perform Mantel test for presence-absence of taxa
print("Mantel test on Jaccard")
jaccard_mantel <- mantel(jaccard_dist, altitude_dist, perm = 10000)
jaccard_mantel
plot(jaccard_dist, altitude_dist)

#Perform Mantel test for relative abundance of taxa
print("Mantel test on Aitchison")
aitchison_mantel <- mantel(aitchison_dist, altitude_dist, perm = 10000, method="spear")
aitchison_mantel
plot(aitchison_dist, altitude_dist)

#Perform Mantel test for functions
print("Mantel test on function euclidean")
func_mantel <- mantel(func_dist, altitude_dist, perm = 10000, method="spear")
func_mantel
plot(func_dist, altitude_dist)

#Partial Mantel tests using geographic
print("Mantel test on Jaccard after accounting for geographic distance")
jaccard_mantel_partial <- mantel.partial(jaccard_dist, altitude_dist, log(geo_dist+1), perm = 10000)
jaccard_mantel_partial

print("Mantel test on Aitchison after accounting for geographic distance")
aitchison_mantel_partial <-mantel.partial(aitchison_dist, altitude_dist, log(geo_dist+1), perm = 10000, method="spear")
aitchison_mantel_partial

print("Mantel test on function euclidead after accounting for geographic distance")
func_mantel_partial <-mantel.partial(func_dist, altitude_dist, log(geo_dist+1), perm = 10000, method="spear")
func_mantel_partial

print("How correlated are altitudinal and geographic?")
mantel(log(geo_dist+1), altitude_dist, perm = 10000) #Very correlated
plot(log(geo_dist+1), altitude_dist)

### Plot metric - altitude - geo relationships ####
distances <- data.frame(jaccard=as.vector(jaccard_dist),
                        aitchison=as.vector(aitchison_dist), 
                        functional=as.vector(func_dist), 
                        altitude=as.vector(altitude_dist),
                        log_geo=as.vector(log(geo_dist+1)))

jacc_plot <-
 ggplot(distances, aes(y=altitude, x=jaccard)) + 
    geom_point(aes(colour=log_geo), size=3) +
    scale_colour_gradient2(low="#8DA0CB", mid="#FFD92F", high="#FC8D62",
                           guide_colorbar(barheight = 2, title = "log-transformed geographical distance (degrees)", midpoint=mean(distances$log_geo), draw.ulim=FALSE)) +
    theme_bw() + theme(legend.position="none", axis.title=element_text(size=15), plot.title=element_text(size=15)) + ggtitle("a)") +
    ylab("Altitudinal distance (km)") + xlab("Jaccard") 


aitch_plot <-
 ggplot(distances, aes(y=altitude, x=aitchison)) + 
    geom_point(aes(colour=log_geo), size=3) +
    scale_colour_gradient2(low="#8DA0CB", mid="#FFD92F", high="#FC8D62",
                           guide_colorbar(barheight = 2, title = "log-transformed geographical distance (degrees)", midpoint=mean(distances$log_geo), draw.ulim=FALSE, label.position="bottom")) +
    theme_bw() + theme(legend.position="bottom", axis.title=element_text(size=15), axis.title.y=element_blank(), plot.title=element_text(size=15)) + ggtitle("b)") +
    ylab("Altitudinal distance (km)") + xlab("Aitchison")
    
        
func_plot <- 
 ggplot(distances, aes(y=altitude, x=functional)) + 
    geom_point(aes(colour=log_geo), size=3) +
    scale_colour_gradient2(low="#8DA0CB", mid="#FFD92F", high="#FC8D62",
                           guide_colorbar(barheight = 2, title = "log-transformed geographical distance (degrees)", midpoint=mean(distances$log_geo), draw.ulim=FALSE)) +
    theme_bw()+ theme(legend.position="none", axis.title=element_text(size=15), axis.title.y=element_blank(), plot.title=element_text(size=15)) + ggtitle("c)") +
    ylab("Altitudinal distance (km)") + xlab("Euclidean functional")
    
#Plot grid
ggsave(plot_grid(jacc_plot, aitch_plot, func_plot,
       align = "h", axis="tb", ncol = 3),
       file="T3_community-level/altitude_vs_geo_scatteplots.png", device="png", height=6, width=15)
    
#### Mantel test only for Grauer's ####
#microbial composition
jaccard_dist_grauers <- vegdist(t(subset_samples(spe_data_final, Spec.subspecies=="graueri")@otu_table), "jaccard")
aitchison_dist_grauers <- vegdist(t(subset_samples(spe_data_final_norm, Spec.subspecies=="graueri")@otu_table), "euclidean")

#Functional composition
func_dist_grauers <- vegdist(t(subset_samples(GO_BP_phyloseq, Spec.subspecies=="graueri")@otu_table), "euclidean")

#altitude
altitude_dist_grauers <- vegdist(subset_samples(spe_data_final, Spec.subspecies=="graueri")@sam_data$Approximate.altitude, na.rm=TRUE, "euclidean")

#geographic distance
geo_dist_grauers <- vegdist(cbind(subset_samples(spe_data_final, Spec.subspecies=="graueri")@sam_data$Approximate.latitude, 
                          subset_samples(spe_data_final, Spec.subspecies=="graueri")@sam_data$Approximate.longitude), "euclidean")
                          
print("Mantel test on Jaccard after accounting for geographic distance - Only Grauer's gorillas")
jaccard_mantel_partial_Gbg <- mantel.partial(jaccard_dist_grauers, altitude_dist_grauers, log(geo_dist_grauers+1), perm = 10000)
jaccard_mantel_partial_Gbg

print("Mantel test on Aitchison after accounting for geographic distance - Only Grauer's gorillas")
aitchison_mantel_partial_Gbg <-mantel.partial(aitchison_dist_grauers, altitude_dist_grauers, log(geo_dist_grauers+1), perm = 10000, method="spear")
aitchison_mantel_partial_Gbg

print("Mantel test on function euclidead after accounting for geographic distance - Only Grauer's gorillas")
func_mantel_partial_Gbg <-mantel.partial(func_dist_grauers, altitude_dist_grauers, log(geo_dist_grauers+1), perm = 10000, method="spear")
func_mantel_partial_Gbg

distances_grauers <- data.frame(jaccard=as.vector(jaccard_dist_grauers),
                        aitchison=as.vector(aitchison_dist_grauers), 
                        functional=as.vector(func_dist_grauers), 
                        altitude=as.vector(altitude_dist_grauers),
                        log_geo=as.vector(log(geo_dist_grauers+1)))

ggplot(distances_grauers, aes(y=altitude, x=log_geo)) + 
    geom_point(aes(colour=jaccard, size=3)) +
    theme_bw() + ylab("Altitudinal distance") + xlab("Geographical distance (log-transformed)") +
    scale_colour_gradient2(low="blue", mid="yellow", high="red", midpoint=mean(distances$jaccard))

                          
#### PERMANOVA ####
print("PERMANOVA, including altitude")
read_count <- spe_data_final@sam_data$readcount.m.before.Kraken %>% as.numeric
seq_centre <- spe_data_final@sam_data$Seq.centre %>% factor
spec.subspecies <- spe_data_final@sam_data$Spec.subspecies %>% factor(level=c("gorilla", "graueri", "beringei"))
altitude <- spe_data_final@sam_data$Approximate.altitude %>% as.numeric

model2_jaccard_w_altitude <- adonis(t(otu_table(spe_data_final)) ~ read_count + seq_centre + spec.subspecies + altitude, permutations = 10000, method = "jaccard")

model2_jaccard_w_altitude

model2_clr_w_altitude <- adonis(t(otu_table(spe_data_final_norm)) ~ read_count + seq_centre + spec.subspecies + altitude, permutations = 10000, method = "euclidean")

model2_clr_w_altitude

read_count <- GO_BP_phyloseq@sam_data$readcount.m.before.Kraken %>% as.numeric
seq_centre <- GO_BP_phyloseq@sam_data$Seq.centre %>% factor
spec.subspecies <- GO_BP_phyloseq@sam_data$Spec.subspecies %>% factor(level=c("gorilla", "graueri", "beringei"))
altitude <- GO_BP_phyloseq@sam_data$Approximate.altitude %>% as.numeric

GO_BP_model_w_altitude <- adonis(t(otu_table(GO_BP_phyloseq)) ~ read_count + seq_centre + spec.subspecies + altitude, permutations = 10000, method = "euclidean")

GO_BP_model_w_altitude

save.image()
#Stop logging
sink(file = NULL)
