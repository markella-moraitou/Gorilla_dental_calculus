library(tidyverse)
library(phyloseq)
library(ggpubr)

setwd("Gorilla_dental_calculus_zenodo")
load(".RData")

pal<-c(get_palette(palette = "Set2",k=3)[3],
       get_palette(palette = "Set2",k=3)[1],
       "#ef885d",
       "#f5b499")

alt.fix<-read_tsv("T3_community-level/new_altitude.csv") %>% 
  right_join(sample_data(spe_data_final) %>% data.frame() %>% rownames_to_column("Sample_ID")) %>%
  mutate(subspecies_locality=ifelse(Spec.subspecies=="graueri",
                                    ifelse(`Approximate altitude`>1000,"graueri >1000","graueri ≤1000"),
                                    as.character(Spec.subspecies)))
alt.fix %>% 
  ggplot(aes(fill=subspecies_locality,x=`Approximate altitude`))+
  geom_histogram(color="black",bins = 50)+
  scale_fill_manual(values = pal)+
  see::theme_modern()+
  scale_y_continuous(breaks = c(2,4,6,8,10,12))+
  labs(y="Number of Samples")+
  geom_vline(xintercept = 1000,color=get_palette(palette = "Set2",k=3)[2],linetype="dotted",size=2)+
  theme(axis.text = element_text(size=20),
        axis.ticks = element_line(),
        axis.title = element_text(size=20),
        legend.position = "top",
        legend.title = element_blank())
ggsave(last_plot(), file="altitude_dist.png", device="png", height=6, width=8,dpi = 300)
