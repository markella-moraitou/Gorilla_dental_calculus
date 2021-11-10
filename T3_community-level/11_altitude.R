#Start logging
sink(file = "log11_altitude.txt")

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

load(".RData")

#Check effect of altitude - Running the full script also requires the output of HUMAnN2

#Load location info

#Get functional analysis phyloseq
GO_BP_phyloseq <- readRDS("/proj/sllstore2017021/nobackup/MARKELLA/F2_functional_stats/GO_BP_phyloseq")
GO_eucl <- readRDS("/proj/sllstore2017021/nobackup/MARKELLA/F2_functional_stats/GO_eucl")

#### Grauer's grouped altitude ####

spe_data_final@sam_data$subspecies_locality <-
  #For Western Lowland and Mountain gorillas, use the subspecies classification without any grouping
  ifelse(spe_data_final@sam_data$Spec.subspecies %in% c("gorilla", "beringei"), as.character(spe_data_final@sam_data$Spec.subspecies),
    #For Grauer's group by altitude (less vs greater than 1500)
    ifelse(as.numeric(spe_data_final@sam_data$Approximate.altitude) < 1500, "graueri <1500", "graueri >1500"))

#Turn into factor
sample_data(spe_data_final)$subspecies_locality <-
  factor(sample_data(spe_data_final)$subspecies_locality, 
         levels=c("gorilla", "graueri <1500", "graueri >1500", "beringei"))

sample_data(spe_data_final)$Spec.subspecies <- factor(sample_data(spe_data_final)$Spec.subspecies, levels = c("gorilla", "graueri", "beringei"))

#Get the same info in the normalized phyloseq and in the functional phyloseq
sample_data(spe_data_final_norm) <- sample_data(spe_data_final)
sample_data(GO_BP_phyloseq) <- sample_data(spe_data_final)

#### Plot ordination ####
ggsave(
  plot_ordination(spe_data_final_norm, clr_5, color="Spec.subspecies", shape = "subspecies_locality", title="Aitchison distances") + 
  theme_light() + geom_point(size=3) + theme(legend.text = element_text(size=12), legend.title = element_text(size=12),
                                             legend.position = "right", plot.title = element_text(size=15)) +
  scale_colour_brewer(palette="Set2", name="Host subspecies", labels = c("beringei"="Mountain", "gorilla"="Western", "graueri"="Grauer's")) + 
  scale_shape_manual(values=c(16, 15, 17, 16), name="Grauer's altitude", labels = c("beringei"="Mountain", "gorilla"="Western", "graueri <1500"="Grauer's <1500m", "graueri >1500"="Grauer's >1500m")),
  file="/proj/sllstore2017021/nobackup/MARKELLA/T3_community-level/clr_with_altitude.png",
  device="png")
  
ggsave(
  plot_ordination(spe_data_final, jaccard_5, color="Spec.subspecies", shape = "subspecies_locality", title="Jaccard distances") + 
  theme_light() + geom_point(size=3) + theme(legend.text = element_text(size=12), legend.title = element_text(size=12),
                                             legend.position = "right", plot.title = element_text(size=15))  +
  scale_colour_brewer(palette="Set2", name="Host subspecies", labels = c("beringei"="Mountain", "gorilla"="Western", "graueri"="Grauer's")) + 
  scale_shape_manual(values=c(16, 15, 17, 16), name="Grauer's altitude", labels = c("beringei"="Mountain", "gorilla"="Western", "graueri <1500"="Grauer's <1500m", "graueri >1500"="Grauer's >1500m")),
  file="/proj/sllstore2017021/nobackup/MARKELLA/T3_community-level/jaccard_with_altitude.png", device="png")
  
ggsave(
  plot_ordination(GO_BP_phyloseq, GO_eucl, color="Spec.subspecies", shape = "subspecies_locality", title="Jaccard distances") + 
  theme_light() + geom_point(size=3) + theme(legend.text = element_text(size=12), legend.title = element_text(size=12),
                                             legend.position = "right", plot.title = element_text(size=15))  +
  scale_colour_brewer(palette="Set2", name="Host subspecies", labels = c("beringei"="Mountain", "gorilla"="Western", "graueri"="Grauer's")) + 
  scale_shape_manual(values=c(16, 15, 17, 16), name="Grauer's altitude", labels = c("beringei"="Mountain", "gorilla"="Western", "graueri < 1500"="Grauer's <1500m", "graueri > 1500"="Grauer's >1500m")),
  file="/proj/sllstore2017021/nobackup/MARKELLA/T3_community-level/GO_abund_with_altitude.png", device="png")


#### PERMANOVA ####
print("PERMANOVA, including subspecies locality")
read_count <- sample_data(spe_data_final)$readcount.m.before.Kraken
seq_centre <- sample_data(spe_data_final)$Seq.centre
spec.subspecies <- sample_data(spe_data_final)$Spec.subspecies
subspecies_locality <- sample_data(spe_data_final)$subspecies_locality

model2_jaccard_w_locality <- adonis(t(otu_table(spe_data_final)) ~ read_count + seq_centre + spec.subspecies + subspecies_locality, permutations = 10000, method = "jaccard")
model2_jaccard_w_locality

model2_clr_w_locality <- adonis(t(otu_table(spe_data_final_norm)) ~ read_count + seq_centre + spec.subspecies + subspecies_locality, permutations = 10000, method = "euclidean")
model2_clr_w_locality

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
print("Mantel test on Jaccard distances")
jaccard_mantel <- mantel(jaccard_dist, altitude_dist, perm = 10000)
jaccard_mantel
plot(jaccard_dist, altitude_dist)

#Perform Mantel test for relative abundance of taxa
print("Mantel test on Aitchison distances")
aitchison_mantel <- mantel(aitchison_dist, altitude_dist, perm = 10000, method="spear")
aitchison_mantel
plot(aitchison_dist, altitude_dist)

#Perform Mantel test for functions
print("Mantel test on function euclidean distances")
func_mantel <- mantel(func_dist, altitude_dist, perm = 10000, method="spear")
func_mantel
plot(func_dist, altitude_dist)

#Partial Mantel tests using geographic distances
print("Mantel test on Jaccard distances after accounting for geographic distance")
jaccard_mantel_partial <- mantel.partial(jaccard_dist, altitude_dist, log(geo_dist+1), perm = 10000)
jaccard_mantel_partial

print("Mantel test on Aitchison distances after accounting for geographic distance")
aitchison_mantel_partial <-mantel.partial(aitchison_dist, altitude_dist, log(geo_dist+1), perm = 10000, method="spear")
aitchison_mantel_partial

print("Mantel test on function euclidead distances after accounting for geographic distance")
func_mantel_partial <-mantel.partial(func_dist, altitude_dist, log(geo_dist+1), perm = 10000, method="spear")
func_mantel_partial

print("How correlated are altitudinal and geographic distances?")
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
    ylab("Altitudinal distance (km)") + xlab("Jaccard distances") 


aitch_plot <-
 ggplot(distances, aes(y=altitude, x=aitchison)) + 
    geom_point(aes(colour=log_geo), size=3) +
    scale_colour_gradient2(low="#8DA0CB", mid="#FFD92F", high="#FC8D62",
                           guide_colorbar(barheight = 2, title = "log-transformed geographical distance (degrees)", midpoint=mean(distances$log_geo), draw.ulim=FALSE, label.position="bottom")) +
    theme_bw() + theme(legend.position="bottom", axis.title=element_text(size=15), axis.title.y=element_blank(), plot.title=element_text(size=15)) + ggtitle("b)") +
    ylab("Altitudinal distance (km)") + xlab("Aitchison distances")
    
        
func_plot <- 
 ggplot(distances, aes(y=altitude, x=functional)) + 
    geom_point(aes(colour=log_geo), size=3) +
    scale_colour_gradient2(low="#8DA0CB", mid="#FFD92F", high="#FC8D62",
                           guide_colorbar(barheight = 2, title = "log-transformed geographical distance (degrees)", midpoint=mean(distances$log_geo), draw.ulim=FALSE)) +
    theme_bw()+ theme(legend.position="none", axis.title=element_text(size=15), axis.title.y=element_blank(), plot.title=element_text(size=15)) + ggtitle("c)") +
    ylab("Altitudinal distance (km)") + xlab("Euclidean functional distances")
    
#Plot grid
ggsave(plot_grid(jacc_plot, aitch_plot, func_plot,
       align = "h", axis="tb", ncol = 3),
       file="/proj/sllstore2017021/nobackup/MARKELLA/T3_community-level/altitude_vs_geo_scatteplots.png", device="png", height=6, width=15)
    
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
                          
print("Mantel test on Jaccard distances after accounting for geographic distance - Only Grauer's gorillas")
jaccard_mantel_partial_Gbg <- mantel.partial(jaccard_dist_grauers, altitude_dist_grauers, log(geo_dist_grauers+1), perm = 10000)
jaccard_mantel_partial_Gbg

print("Mantel test on Aitchison distances after accounting for geographic distance - Only Grauer's gorillas")
aitchison_mantel_partial_Gbg <-mantel.partial(aitchison_dist_grauers, altitude_dist_grauers, log(geo_dist_grauers+1), perm = 10000, method="spear")
aitchison_mantel_partial_Gbg

print("Mantel test on function euclidead distances after accounting for geographic distance - Only Grauer's gorillas")
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
