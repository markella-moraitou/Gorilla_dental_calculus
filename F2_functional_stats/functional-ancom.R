library(readr)
library(nlme)
library(tidyverse)
library(ggplot2)
library(compositions)
library(phyloseq)
library(microbiome)
source("T3_community-level/ancom-functions.R")
setwd("~/Lab-Notes/")

# set the palette
pal<-c(ggpubr::get_palette(palette = "Set2",k=3)[2],ggpubr::get_palette(palette = "Set2",k=3)[3],ggpubr::get_palette(palette = "Set2",k=3)[1])

###########
GO_BP_phyloseq<-readRDS("F2_functional_stats/GO_BP_phyloseq")
out<-as.data.frame(otu_table(GO_BP_phyloseq)) %>% mutate_all(as.integer)
meta<-as.data.frame(as.matrix(sample_data(GO_BP_phyloseq))) %>% 
  mutate(Spec.subspecies=recode(Spec.subspecies,
                                "gorilla"="Western",
                                "graueri"="Grauer's",
                                "beringei"="Mountain"))

samples<-bind_cols("SampleID"=colnames(out),"Spec.subspecies"=meta[colnames(out),"Spec.subspecies"],"readcount.m.before.Kraken"=as.numeric(meta[colnames(out),"readcount.m.before.Kraken"]))

feature_table = out
sample_var = "SampleID"
group_var = "Spec.subspecies"
out_cut = 0.05
zero_cut = 0.90
lib_cut = 0
neg_lb = FALSE

prepro = feature_table_pre_process(feature_table, samples, sample_var, group_var,out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

###### Step 2: ANCOM
main_var = "Spec.subspecies" # taxa that are uniquely found in one disease state will be labelled as structural zeros
p_adj_method = "BH"
alpha = 0.05
adj_formula = "readcount.m.before.Kraken"
rand_formula = NULL

# run ANCOM
res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method,alpha, adj_formula, rand_formula)

# save results
write_csv(cbind(res$out,struc_zero), "F2_functional_stats/BP_ancom-read-depth.csv")

# kw test
res$out

###### Step 3: Volcano Plot
# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1), 0.5 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6", "detected_0.5")

# Annotation data
dat_ann = data.frame(x = min(res$fig$data$x), y = cut_off["detected_0.7"], label = "W[0.7]")

volcano_p<-GO_BP_phyloseq %>% microbiome::transform(transform = "clr") %>% 
  otu_table() %>% data.frame(.) %>% 
  rownames_to_column("taxa_id") %>% 
  pivot_longer(cols=-taxa_id,names_to = "SampleID",values_to = "clr") %>% 
  right_join(read_csv("F2_functional_stats/BP_ancom-read-depth.csv")) %>%
  left_join(samples) %>% 
  group_by(taxa_id,W) %>% 
  summarise(mean_difference_mountain_western=mean(clr[Spec.subspecies=="Mountain"]-clr[Spec.subspecies=="Western"],na.rm=T),
            mean_difference_mountain_grauers=mean(clr[Spec.subspecies=="Mountain"]-clr[Spec.subspecies=="Grauer's"],na.rm=T),
            mean_difference_grauers_mountain=mean(clr[Spec.subspecies=="Grauer's"]-clr[Spec.subspecies=="Mountain"],na.rm=T)) %>% 
  pivot_longer(cols = c(-taxa_id,-W),names_to = "contrast",values_to = "clr.mean.diff") %>% 
  ggplot(aes(x=clr.mean.diff,y=W,color=contrast))+geom_point()+
  geom_hline(yintercept = cut_off["detected_0.9"], linetype = "dashed") +
  annotate(geom = "text",x=-3.5,y = cut_off["detected_0.9"]+20,label="W-statistic 0.9 quantile") +
  geom_hline(yintercept = cut_off["detected_0.8"], linetype = "dashed") +
  geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") +
  geom_hline(yintercept = cut_off["detected_0.6"], linetype = "dashed") +
  geom_hline(yintercept = cut_off["detected_0.5"], linetype = "dashed")

ggsave(volcano_p,filename = "BP_ancom-read-depth.png",dpi=300,width=8,height = 8,units="in")

# kruskal wallace
kw<-GO_BP_phyloseq %>% microbiome::transform(transform = "clr") %>% 
  otu_table() %>% data.frame(.) %>% 
  rownames_to_column("taxa_id") %>% 
  pivot_longer(cols=-taxa_id,names_to = "SampleID",values_to = "clr") %>% 
  right_join(res$out) %>%
  filter(W >= cut_off["detected_0.8"]) %>% 
  left_join(samples) %>% 
  group_by(taxa_id,W) %>% 
  summarise(mean_difference_mountain_western=mean(clr[Spec.subspecies=="Mountain"]-clr[Spec.subspecies=="Western"],na.rm=T),
            mean_difference_mountain_grauers=mean(clr[Spec.subspecies=="Mountain"]-clr[Spec.subspecies=="Grauer's"],na.rm=T),
            mean_difference_grauers_mountain=mean(clr[Spec.subspecies=="Grauer's"]-clr[Spec.subspecies=="Mountain"],na.rm=T)) %>% 
  pivot_longer(cols = c(-taxa_id,-W),names_to = "contrast",values_to = "clr.mean.diff")

library(foreach)
foreach(i=unique(kw$taxa_id),.combine = bind_rows) %do% {
  kruskal.test(.)$p.value
}

#### ANCOM PLOT
library(ggpubr)
ancom_res<-read_csv("F2_functional_stats/BP_ancom-read-depth.csv")
out.clr<-out %>% otu_table(taxa_are_rows = T) %>% 
  microbiome::transform(.,transform = "clr") %>% 
  as.data.frame() %>% 
  rownames_to_column("GOID")
ft<- ancom_res %>% 
  filter(detected_0.9==TRUE & !is.infinite(W)) %>%
  left_join(out.clr,by=c("taxa_id"="GOID")) %>% 
  pivot_longer(cols=c(-taxa_id,-W,-detected_0.9,-detected_0.8,-detected_0.7,-detected_0.6),names_to = "SampleID",values_to = "abundance") %>% #na.omit() %>% 
  left_join(samples,by="SampleID") %>%
  filter(!is.na(abundance) & !is.na(Spec.subspecies))

comps <- microbiomeutilities::make_pairs(unique(ft$Spec.subspecies))

ancom_p<-ft %>% 
  mutate(taxa_id=as.factor(trimws(gsub(".*] ","",taxa_id))),
         Spec.subspecies=fct_relevel(Spec.subspecies,"Western","Grauer's","Mountain"),
         taxa_id=fct_relevel(taxa_id,"cellular respiration",
                             "response to tellurium ion",
                             "transition metal ion transport",
                             "pyrimidine nucleotide biosynthetic process",
                             "L-threonine catabolic process to glycine",
                             "protein phosphorylation")) %>% 
  ggplot(aes(x=Spec.subspecies,y=abundance,fill=Spec.subspecies))+
  geom_boxplot(color="black",outlier.shape = 21)+
  facet_grid(.~stringr::str_wrap(taxa_id,width = 30))+
  # geom_signif(comparisons=list(c("Mountain","Grauer's"),c("Mountain","Western"),c("Grauer's","Western")),test = "wilcox.test",map_signif_level = T,tip_length = 0.01,step_increase = 0.04,vjust = 0.5,textsize = 5)+
  labs(y = "CLR Normalized Abundances",x="Subspecies")+
  scale_fill_manual(values = c(pal[3],pal[1],pal[2]))+
  guides(fill=guide_legend("Subspecies"))+
  scale_y_continuous(limits = c(0,6.5),expand = c(0, 0)) +
  theme_bw()+
  theme(axis.text.x = element_blank(),legend.position = "bottom",
        axis.text.y = element_text(size=14),
        axis.title.y = element_text(size=16),
        strip.text = element_blank(),
        axis.ticks.x = element_blank(),axis.title.x = element_blank(),
        plot.background = element_rect(fill="white",color="white"),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16),
        legend.spacing.x = unit(0.75, 'cm'),
        text = element_text(family="calibri"))
ggsave(ancom_p,filename = "F2_functional_stats/ancom_0.9.png",dpi = 300,width = 16,height = 8,units = "in")
