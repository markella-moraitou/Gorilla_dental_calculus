library(tidyverse)
library(vegan)
library(ecodist)
library(phyloseq)
library(microbiome)
library(ggpubr)

setwd("Gorilla_dental_calculus_zenodo")
# load in objects for microbial dataset
load(".RData")

# determine groups for later
groups<-spe_data_final@sam_data %>% data.frame() %>% 
  select(Spec.subspecies,Spec.locality)

# load in dietary object
euk_genus_diet<-readRDS("D3_diet_stats/euk_genus_diet.rds")
euk_genus_diet_norm<-readRDS("D3_diet_stats/euk_genus_diet_norm.rds")

# euclidean
# diet distance
diet_euc<-euk_genus_diet_norm %>% 
  aggregate_taxa(.,"genus") %>% 
  abundances() %>% data.frame() %>% 
  select(which(colSums(.)!=0)) %>% 
  t() %>% 
  vegdist("euclidean")

# microbial distance
micro_euc<-spe_data_final_norm %>% 
  aggregate_taxa(.,"genus") %>%
  abundances() %>% data.frame() %>% 
  select(colnames(as.matrix(diet_euc))) %>% 
  t() %>% 
  vegdist("euclidean")

mrm_euc<-MRM(micro_euc ~ diet_euc,nperm = 10000)

# convert to long format
diet_euc_long<-diet_euc %>% as.matrix() %>% data.frame() %>% rownames_to_column("Sample 1") %>% 
  pivot_longer(cols = -`Sample 1`,
               names_to = "Sample 2",values_to = "diet")

micro_euc_long<-micro_euc %>% as.matrix() %>% data.frame() %>% rownames_to_column("Sample 1") %>% 
  pivot_longer(cols = -`Sample 1`,
               names_to = "Sample 2",values_to = "micro")

# and plot
diet_euc_p<-diet_euc_long %>% 
  left_join(micro_euc_long) %>% 
  filter((diet+micro) > 0) %>% 
  mutate(group=ifelse(groups[`Sample 1`,"Spec.subspecies"]==groups[`Sample 2`,"Spec.subspecies"],"same subspecies","different subspecies")) %>% 
  ggplot(aes(x=diet,y=micro))+
  geom_point(aes(color=group),alpha=0.5)+
  geom_smooth(method = "lm",color="black")+
  labs(y="Microbiome Aitchison dissimilarity",
       x="Dietary Aitchison dissimilarity")+
  annotate(geom = "text",x = 17,y=150,size=4,label=paste0("R^2 == ",round(mrm_euc$r.squared[1], digits = 3)),parse=TRUE)+
  annotate(geom = "text",x = 25,y=149,size=4,label=paste0("p == ",round(mrm_euc$coef[2,2],digits = 3)),parse=TRUE)+
  see::theme_modern()+
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=16),
        legend.title = element_blank(),
        legend.position = c(0.8,0.9),
        legend.box.background = element_rect(fill="gray"))
ggsave(diet_euc_p, file="euc_mrm.png", device="png", height=6, width=6,dpi = 300)

# jaccard
# diet distance
diet_jac<-euk_genus_diet %>% 
  aggregate_taxa("genus") %>%
  abundances() %>% data.frame() %>% 
  mutate(across(where(is.numeric), ~ .x+abs(min(.x, na.rm = TRUE)))) %>%
  select(which(colSums(.)!=0)) %>% 
  t() %>% 
  vegdist("jaccard")

# microbial distance
micro_jac<-spe_data_final %>% 
  aggregate_taxa(.,"genus") %>%
  abundances() %>% data.frame() %>% 
  mutate(across(where(is.numeric), ~ .x+abs(min(.x, na.rm = TRUE)))) %>% 
  select(colnames(as.matrix(diet_jac))) %>% 
  t() %>% 
  vegdist("jaccard")

mrm_jac<-MRM(micro_jac ~ diet_jac,nperm = 10000)

# convert to long format
diet_jac_long<-diet_jac %>% as.matrix() %>% data.frame() %>% rownames_to_column("Sample 1") %>% 
  pivot_longer(cols = -`Sample 1`,
               names_to = "Sample 2",values_to = "diet")

micro_jac_long<-micro_jac %>% as.matrix() %>% data.frame() %>% rownames_to_column("Sample 1") %>% 
  pivot_longer(cols = -`Sample 1`,
               names_to = "Sample 2",values_to = "micro")

# and plot
diet_jac_p<-diet_jac_long %>% 
  left_join(micro_jac_long) %>% 
  filter((diet+micro) > 0) %>% 
  mutate(group=ifelse(groups[`Sample 1`,"Spec.subspecies"]==groups[`Sample 2`,"Spec.subspecies"],"same subspecies","different subspecies")) %>% 
  ggplot(aes(x=diet,y=micro))+
  geom_point(aes(color=group),alpha=0.5)+
  geom_smooth(method = "lm",color="black")+
  labs(y="Microbiome Jaccard dissimilarity",
       x="Dietary Jaccard dissimilarity")+
  annotate(geom = "text",x = 0.5,y=0.35,label=paste0("R2 = ",format(mrm_jac$r.squared[1], digits = 3),", p = ",format(mrm_jac$coef[2,2],digits = 3)))+
  see::theme_modern()+
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=16))
ggsave(diet_jac_p, file="jac_mrm.png", device="png", height=6, width=8,dpi = 300)  

library(cowplot)
legend<-get_legend(diet_euc_p)
ggdraw(plot_grid(plot_grid(diet_euc_p+theme(legend.position = "none"),diet_jac_p+theme(legend.position = "none"),
                 labels = c("a)","b)"),label_size = 16),
       plot_grid(NULL, legend, ncol=1),
       rel_widths=c(1, 0.2)))
ggsave(last_plot(), file="euc_jac_mrm.png", device="png", height=6, width=12,dpi = 300)  
