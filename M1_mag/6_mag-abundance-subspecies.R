library(tidyverse)
library(ggpubr)
library(ggsignif)
library(plotly)
library(microbiome)
library(pheatmap)

meta<-readxl::read_excel("metadata_gorilla_analysis.xlsx")

# for now, just care about genera level classifications
tax<-read.delim("M1_mag/gtdbtk.bac120.summary.tsv",sep = "\t",header=T) %>% 
  separate(col = classification,sep = ";",into=c("domain","phylum","class","order","family","genus","species")) %>% 
  mutate(across(c(genus,species,family),~gsub("s__|g__|f__","",.)),
         # genus=make.unique(genus),
         species=ifelse(genus=="",make.unique(family),ifelse(species=="",make.unique(genus),species))) %>% 
  select(user_genome,family,genus,species)

dat<-read.delim("M1_mag/bin_abundance_table.tab",header = T) %>% 
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
  group_by(Genomic.bins,family,genus,species) %>% 
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
ggsave(p,file="M1_mag/mag-sspecies-membership.png",dpi=300)

# how many samples make up each subspecies?
dat %>% 
  filter(!is.na(Spec.subspecies) & Spec.subspecies!="museum control") %>% 
  group_by(Genomic.bins,Spec.subspecies) %>% 
  summarise(n=n_distinct(Seq.label)) %>% 
  group_by(Spec.subspecies) %>% 
  summarise(mean=mean(n),
            max=max(n),
            min=min(n))
  
dat10p<-dat %>% 
  filter(!is.na(Spec.subspecies) & Spec.subspecies!="museum control") %>% 
  group_by(species) %>%
  mutate(sum=sum(abundance,na.rm=T),
         n=n_distinct(Seq.label[abundance!=0])) %>% 
  group_by(Seq.label) %>% 
  mutate(rel.abund=abundance/sum) %>% 
  group_by(Genomic.bins,Spec.subspecies) %>%
  mutate(n_samples=n_distinct(Seq.label)) %>%
  filter(n_samples>=8) # only keep MAGs with at (roughly) 50% prevalence?

dat10p %>% 
  ggplot(aes(x=Spec.subspecies,y=rel.abund,fill=Spec.subspecies))+
  geom_boxplot()+
  labs(y="Relative Abundance")+
  scale_x_discrete(limits=c("Western","Grauer's","Mountain"))+
  facet_wrap(~species,scales = "free_y")+
  guides(fill=guide_legend("Subspecies"))+
  scale_fill_manual(values = pal)+
  # scale_y_log10()+
  geom_signif(comparisons=list(c("Mountain","Grauer's"),c("Mountain","Western"),c("Grauer's","Western")),map_signif_level = T,step_increase = 0.2)+
  theme_bw()+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),
        strip.text.y = element_text(angle = 0),
        legend.position = "none")


`Lactobacillus_H gorillae.p`<-dat10p %>% 
    filter(species=="Lactobacillus_H gorillae") %>% 
    mutate(Spec.subspecies=fct_relevel(Spec.subspecies,"Western","Grauer's","Mountain"),
           species="Limosilactobacillus gorillae") %>% 
    ggplot(aes(x=Spec.subspecies,y=rel.abund,fill=Spec.subspecies))+
    geom_boxplot(outlier.shape = NA)+
    labs(y="Relative Abundance")+
    # scale_x_discrete(limits=c("Western","Grauer's","Mountain"))+
    facet_wrap(~species,scales = "free_y")+
    guides(fill=guide_legend("Subspecies"))+
    # scale_y_continuous(expand = expansion(mult = 0.15))+
    coord_cartesian(ylim=c(0,0.05))+
    scale_fill_manual(values = pal)+
  geom_signif(comparisons=list(c("Mountain","Grauer's"),c("Mountain","Western"),c("Grauer's","Western")),y_position = c(0.01,0.011,0.012),map_signif_level = T,step_increase = 0.01,tip_length = 0.0035)+
  theme_bw()+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size=14),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text.x = element_text(angle = 0,size = 16),legend.position = "bottom",
          legend.title = element_text(size=16),
          legend.text = element_text(size = 14))

Neisseria.p<-dat10p %>% 
  filter(species=="Neisseria") %>% 
  mutate(Spec.subspecies=fct_relevel(Spec.subspecies,"Western","Grauer's","Mountain")) %>% 
  ggplot(aes(x=Spec.subspecies,y=rel.abund,fill=Spec.subspecies))+
  geom_boxplot(outlier.shape = NA)+
  labs(y="Relative Abundance")+
  # scale_x_discrete(limits=c("Western","Grauer's","Mountain"))+
  facet_wrap(~species,scales = "free_y")+
  guides(fill=guide_legend("Subspecies"))+
  # scale_y_continuous(expand = expansion(mult = 0.15))+
  # coord_cartesian(ylim=c(0,0.15))+
  scale_fill_manual(values = pal)+
  geom_signif(comparisons=list(c("Mountain","Grauer's"),c("Mountain","Western"),c("Grauer's","Western")),y_position = c(0.12,0.125,0.13),map_signif_level = T,step_increase = 0.1)+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=14),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(angle = 0,size = 16),
        legend.position = "none")
legend<-get_legend(`Lactobacillus_H gorillae.p`)
tree.legend<-get_legend(n.tree)

Veillonella.p<-dat10p %>% 
  filter(species=="Veillonella") %>% 
  mutate(Spec.subspecies=fct_relevel(Spec.subspecies,"Western","Grauer's","Mountain")) %>% 
  ggplot(aes(x=Spec.subspecies,y=rel.abund,fill=Spec.subspecies))+
  geom_boxplot(outlier.shape = NA)+
  labs(y="Relative Abundance")+
  # scale_x_discrete(limits=c("Western","Grauer's","Mountain"))+
  facet_wrap(~species,scales = "free_y")+
  guides(fill=guide_legend("Subspecies"))+
  # scale_y_continuous(expand = expansion(mult = 0.15))+
  # coord_cartesian(ylim=c(0,0.15))+
  scale_fill_manual(values = pal)+
  geom_signif(comparisons=list(c("Mountain","Grauer's"),c("Mountain","Western"),c("Grauer's","Western")),y_position = c(0.12,0.125,0.13),map_signif_level = T,step_increase = 0.1)+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=14),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(angle = 0,size = 16),
        legend.position = "none")
legend<-get_legend(`Lactobacillus_H gorillae.p`)
ggsave(Veillonella.p,file="M1_mag/Veillonella-abundance-sig.png",dpi=300,height=8,width = 6,units = "in")

tree.legend<-get_legend(v.tree)

#### BUILD NEISSERIA AND L.GORILLAEA TREES FIRST
# with tree plots
gp<-ggarrange(n.tree,Neisseria.p,lb.tree,`Lactobacillus_H gorillae.p`+theme(legend.position = "bottom"),labels = c("a)","b)","c)","d)"),widths = c(1,0.5,1,0.5),font.label = list(size=28),label.x = c(0,-0.05,0,-0.05))
ggsave(gp,file="M1_mag/mag-sspecies-abundance-50p-sig.png",dpi=300,height=16,width = 17,units = "in")

gp<-dat10p %>% 
    filter(species=="Rothia") %>% 
    mutate(Spec.subspecies=fct_relevel(Spec.subspecies,"Grauer's","Mountain","Western")) %>% 
    ggplot(aes(x=Spec.subspecies,y=rel.abund,fill=Spec.subspecies))+
    geom_boxplot(outlier.shape = 21)+
    labs(y="Relative Abundance")+
    facet_wrap(~species,scales = "free_y")+
    guides(fill=guide_legend("Subspecies"))+
    scale_y_continuous(expand = expansion(mult = 0.15))+
    scale_fill_manual(values = pal)+
    geom_signif(comparisons=list(c("Grauer's","Western"),c("Mountain","Grauer's"),c("Mountain","Western")),map_signif_level = T,step_increase = 0.1)+
    theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=14),
        axis.title.y = element_text(size = 16),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(angle = 0,size = 16),
        legend.position = "none")
ggsave(plot=gp,file="M1_mag/mag-rothia-abundance-prop-sig.png",dpi=300,height=5.5,width = 5,units = "in")


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
pheatmap(mat,
         # na_col = "gray",
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = anno,
         # annotation_colors = ann_colors[1],
         drop_levels = F,
         legend_breaks = c(2,4,6,8,max(mat,na.rm = T)),
         legend_labels = c("2","4","6","8","CLR-Normalized\nAbundance\n"),fontsize = 8,
         cellwidth = 10,cellheight = 10,scale = "none",width = 12,height = 10,annotation_names_col = F,filename = "M1_mag/mag-heatmap.png")
 
