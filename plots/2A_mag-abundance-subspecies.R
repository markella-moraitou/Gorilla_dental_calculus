library(tidyverse)
library(ggpubr)
library(ggsignif)
library(plotly)
library(microbiome)
library(pheatmap)

meta<-read.csv("calculus/all-samples/rank-abundance-filtering/reproducibility-tests/Gorilla-Markella/metadata_gorilla_analysis.csv",header=T)

# for now, just care about genera level classifications
tax<-read.delim("calculus/species-specific/DC-gorilla/mag/data/gtdbtk.bac120.summary.tsv",sep = "\t",header=T) %>% 
  separate(col = classification,sep = ";",into=c("domain","phylum","class","order","family","genus","species")) %>% 
  mutate(across(c(genus,species,family),~gsub("s__|g__|f__","",.)),
         # genus=make.unique(genus),
         species=ifelse(genus=="",make.unique(family),ifelse(species=="",make.unique(genus),species))) %>% 
  select(user_genome,family,genus,species)

dat<-read.delim("calculus/species-specific/DC-gorilla/mag/data/bin_abundance_table.tab",header = T) %>% 
  rename_all(.funs = function(x) gsub("_m_decontam","",x)) %>% 
  pivot_longer(cols = -Genomic.bins,names_to = "Seq.label",values_to = "abundance") %>% 
  left_join(meta) %>% 
  left_join(tax,by=c("Genomic.bins"="user_genome")) %>% 
  mutate(Spec.subspecies=recode(Spec.subspecies,
                                "gorilla"="Western",
                                "graueri"="Grauer's",
                                "beringei"="Mountain"),
    Spec.subspecies=ifelse(grepl("ERR|BS",Seq.label),"museum control",Spec.subspecies),
    Spec.subspecies=fct_relevel(Spec.subspecies,"museum control","Grauer's","Mountain","Western"))

# stats
prev.dat<-dat %>% 
  filter(!is.na(Spec.subspecies) & Spec.subspecies != "museum control") %>% 
  group_by(Genomic.bins,family,genus,species,Spec.subspecies) %>% 
  summarise(prev=n_distinct(Seq.label[abundance!=0])/46) %>% 
  arrange(-prev)

prev.dat %>% group_by(Spec.subspecies) %>% summarise(mean=mean(prev))

# set the palette
pal<-ggpubr::get_palette(palette = "Set2",k=3)
# pal<-c(ggpubr::get_palette(palette = "Set2",k=3)[2],ggpubr::get_palette(palette = "Set2",k=3)[3],ggpubr::get_palette(palette = "Set2",k=3)[1])

# which MAGs in which sub spec.?
p<-dat %>% 
  filter(!is.na(Spec.subspecies)) %>% 
  group_by(Genomic.bins,species,Spec.subspecies) %>% 
  summarise(n=n_distinct(Seq.label[abundance!=0])) %>%
  group_by(Genomic.bins) %>% mutate(sum=sum(n),label=paste(species,Genomic.bins)) %>% 
  ggplot(aes(x=reorder(label,-sum),y=n,fill=Spec.subspecies))+geom_bar(stat="identity")+
  scale_fill_manual(values = c("gray",pal))+
  labs(y="Number of Samples with MAG",x="Taxonomic Identity of MAG")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=65,hjust=1),legend.title = element_blank())
ggsave(p,file="calculus/species-specific/DC-gorilla/mag/images/mag-sspecies-membership.png",dpi=300)

# how many samples make up each subspecies?
dat %>% 
  filter(!is.na(Spec.subspecies) & Spec.subspecies!="museum control") %>% 
  group_by(Genomic.bins,Spec.subspecies) %>% 
  summarise(n=n_distinct(Seq.label)) %>% 
  group_by(Spec.subspecies) %>% 
  summarise(mean=mean(n),
            max=max(n),
            min=min(n))
  
# Rothia
dat %>% 
  filter(!is.na(Spec.subspecies) & Spec.subspecies!="museum control") %>% 
  filter(grepl("Rothia",genus)) %>% # only keep MAGs with at least 50% prevalence?
  group_by(genus) %>%
  mutate(sum=sum(abundance,na.rm=T),
         n=n_distinct(Seq.label[abundance!=0]),
         genus=gsub("g__|f__","",genus)) %>% 
  group_by(Seq.label) %>% 
  mutate(rel.abund=abundance/sum,
         Spec.subspecies=fct_relevel(Spec.subspecies,"Western","Grauer's","Mountain")) %>% 
  ungroup() %>%
  ggplot(aes(x=Spec.subspecies,y=rel.abund,fill=Spec.subspecies))+
  geom_boxplot()+
  labs(y="Relative Abundance")+
  scale_x_discrete(limits=c("Western","Grauer's","Mountain"))+
  guides(fill=guide_legend("Subspecies"))+
  scale_y_continuous(expand = expand_scale(mult = c(0, .3)))+
  scale_fill_manual(values = pal)+
  geom_signif(comparisons=list(c("Mountain","Grauer's"),c("Mountain","Western"),c("Grauer's","Western")),map_signif_level = T,step_increase = 0.2)+
  # coord_cartesian(ylim = c(0,0.5))+
  theme_bw()+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),
        strip.text.y = element_text(angle = 0))

anno<-dat %>% 
  filter(!is.na(Spec.subspecies)) %>% 
  distinct(Seq.label,Spec.subspecies) %>% 
  column_to_rownames("Seq.label") %>% 
  as.data.frame() %>% 
  arrange(Spec.subspecies) %>% 
  rename("Subspecies"="Spec.subspecies") %>% 
  mutate(Subspecies=fct_relevel(Subspecies,"Western","Grauer's","Mountain","museum control"))

mat<-dat %>% 
  mutate(species=gsub("Lactobacillus","Limosilactobacillus gorillae",gsub("_H gorillae","",species))) %>% 
  mutate(species=gsub("RUG013","Eggerthellaceae",species)) %>%
  mutate(species=gsub("BM520","Bacteroidales",species)) %>% 
  mutate(species=gsub("F0058","Paludibacteraceae",species)) %>% 
  # filter(!is.na(Spec.subspecies)) %>% 
  select(species,Seq.label,abundance) %>% 
  distinct(species,Seq.label,.keep_all = T) %>% 
  group_by(species) %>% 
  mutate(prev=n_distinct(Seq.label[abundance!=0])/46) %>% 
  arrange(-prev) %>% 
  select(-prev) %>% 
  pivot_wider(id_cols = "species",names_from = "Seq.label",values_from = "abundance") %>% 
  column_to_rownames("species") %>% otu_table(.,taxa_are_rows = T) %>% 
  microbiome::transform(.,transform = "clr",shift = 0.001) %>% 
  data.frame(.) %>% 
  # arrange(-rowSums(.)) %>%
  relocate(rownames(anno)) %>%
  select(which(anno$Subspecies=="Western"),
         which(anno$Subspecies=="Grauer's"),
         which(anno$Subspecies=="Mountain"),
         which(anno$Subspecies=="museum control")) %>% 
  as.matrix()

ann_colors=list(Subspecies= c(Western=pal[1],`Grauer's`=pal[2],Mountain=pal[3],`museum control`="gray"))

# my_palette <- colorRampPalette(c("yellow", "orange", "red")) (n=20)
# breaks <- seq(min(mat, na.rm = T), max(mat, na.rm = T), length.out = 21)
mat[mat < 0] <- NA
pheatmap(mat,legend = F,fontsize_col = 12,fontsize_row = 12,
         # na_col = "gray",
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = anno,
         annotation_colors = ann_colors,
         drop_levels = F,
         legend_breaks = c(2,4,6,8,max(mat,na.rm = T)),
         legend_labels = c("2","4","6","8","CLR-Normalized\nAbundance\n"),fontsize = 8,
         cellwidth = 10,cellheight = 10,scale = "none",width = 12,height = 10,annotation_names_col = F,filename = "calculus/species-specific/DC-gorilla/mag/images/mag-heatmap.png")
 
# Corynebacterium
dat %>% 
  filter(!is.na(Spec.subspecies) & Spec.subspecies!="museum control") %>% 
  filter(grepl("Corynebacterium",genus)) %>% # only keep MAGs with at least 50% prevalence?
  group_by(genus) %>%
  mutate(sum=sum(abundance,na.rm=T),
         n=n_distinct(Seq.label[abundance!=0]),
         genus=gsub("g__|f__","",genus)) %>% 
  group_by(Seq.label) %>% 
  mutate(rel.abund=abundance/sum,
         Spec.subspecies=fct_relevel(Spec.subspecies,"Western","Grauer's","Mountain")) %>% 
  ungroup() %>%
  ggplot(aes(x=Spec.subspecies,y=rel.abund,fill=Spec.subspecies))+
  geom_boxplot()+
  labs(y="Relative Abundance")+
  scale_x_discrete(limits=c("Western","Grauer's","Mountain"))+
  guides(fill=guide_legend("Subspecies"))+
  # scale_y_continuous(expand = expand_scale(mult = c(0, .3)))+
  scale_fill_manual(values = pal)+
  geom_signif(comparisons=list(c("Mountain","Grauer's"),c("Mountain","Western"),c("Grauer's","Western")),map_signif_level = T,step_increase = 0.2)+
  # coord_cartesian(ylim = c(0,0.5))+
  theme_bw()+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),
        strip.text.y = element_text(angle = 0))
