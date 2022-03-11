library(tidyverse)
library(phyloseq)
library(microbiome)

cols<-list(Spec.subspecies=c("beringei"=scales::hue_pal()(3)[1],"gorilla"=scales::hue_pal()(3)[2],"graueri"=scales::hue_pal()(3)[3]))

# read in unstratified go terms, in a phyloseq object
GO_BP_phyloseq<-readRDS("F2_functional_statsGO_BP_phyloseq") %>% 
  otu_table(GO_BP_phyloseq) %>% data.frame() %>% 
  mutate_all(as.integer)

# metadata with subspecies info
meta<-as.data.frame(as.matrix(sample_data(GO_BP_phyloseq))) %>% rownames_to_column("SampleID")
samples<-meta %>% select(SampleID,Spec.subspecies) %>% 
  filter(!grepl("BL|BE|LIB|BS|EXB|ERR",SampleID)) %>% 
  mutate(Spec.subspecies=recode(Spec.subspecies,
                                "gorilla"="Western",
                                "graueri"="Grauer's",
                                "beringei"="Mountain"))

# load in ANCOM results
# generated from unstratified dataset
# identified sig GO terms
# w = 0.8
ancom_res<-read_csv("F2_functional_statsBP_ancom-read-depth.csv") %>% 
  filter(detected_0.9==TRUE & !is.infinite(W))

# stratified GO terms
humann2_pathways <- read_tsv("F2_functional_statsall_GOterms_cpm_renamed_wspecies.tsv",col_names = T) %>% 
  rename(id = `# Gene Family`) %>%
  separate(id,into = c("id","genus.species"), sep = "\\|") %>%
  separate(genus.species,into = c("genus","species"), sep = "\\.") %>%
  mutate(genus=replace_na("unidentified")) %>% 
  filter(!is.na(genus) & id!="UNGROUPED" & id %in% c("unidentified",unique(ancom_res$taxa_id))) %>%  
  # pivot_longer(cols = contains("Abundance"),names_to = "SampleID",values_to = "abundance") %>% 
  rename_all(.funs = function(x) gsub("_Abundance-RPKs","",x))

humann2_pathways.otu<-humann2_pathways %>% 
  pivot_longer(cols = c(-id,-genus,-species),names_to = "SampleID",values_to = "abundance") %>% 
  group_by(id,genus,species) %>% 
  summarise(total_abundance=sum(abundance,na.rm=T)) %>% rowwise() %>% 
  mutate(species=ifelse(is.na(species),genus,species)) %>% ungroup() %>% 
  select(-genus) %>% 
  pivot_wider(id_cols = id,names_from = species,values_from = total_abundance,values_fill = 0) %>% 
  column_to_rownames("id") %>% 
  otu_table(.,taxa_are_rows = F)

samps<-humann2_pathways %>% 
  pivot_longer(cols=c(-id,-genus,-species),names_to="SampleID") %>% 
  distinct(SampleID,id)

humann2_pathways.clr<- humann2_pathways.otu %>% 
  microbiome::transform(.,transform = "clr",shift = 0.01) %>% 
  data.frame(otu_table(.)) %>% 
  rownames_to_column("id") %>% 
  left_join(samps) # adding samples back in after transformation

rel.p<-humann2_pathways.otu %>% 
  data.frame(otu_table(.,taxa_are_rows = T)) %>% 
  rownames_to_column("id") %>% 
  pivot_longer(cols = c(-id),names_to = "species",values_to = "abundance") %>% 
  # left_join(samps,by="id") %>%
  # left_join(samples,by="SampleID") %>% 
  filter(!is.na(abundance)) %>% 
  mutate(species=gsub("\\..*","",species),
         species=gsub("s__|g__","",species),
         species=(gsub("_"," ",species)),
         id=gsub(".*] ","",id),
         id=fct_relevel(id,"cellular respiration",
                             "response to tellurium ion",
                             "transition metal ion transport",
                             "pyrimidine nucleotide biosynthetic process",
                             "L-threonine catabolic process to glycine",
                             "protein phosphorylation")) %>%
  separate(col = species,sep = " ",into = "genus",remove = F) %>% 
  group_by(id) %>% 
  summarise(rel_abund=abundance/sum(abundance,na.rn=T),
         thresh=quantile(rel_abund,0.98), # only take the top abundant species in each GO term
         cond.fill=ifelse(genus=="unidentified","zunidentified",ifelse(rel_abund>=thresh,genus,"yLow Abundance Taxa"))) %>%
  arrange(cond.fill) %>%
  ggplot(aes(x=id,y=rel_abund))+
  geom_bar(aes(fill=cond.fill),stat="identity")+ 
  labs(y="Relative Abundance of GO Term")+
  scale_fill_manual(values = c(ggpubr::get_palette(palette = "Paired",k=8),"gray","gray30"))+
  scale_y_reverse(labels=rev(c("0.00","0.25","0.50","0.75","1.00")))+
  guides(fill=guide_legend("Genus"))+
  facet_grid(.~str_wrap(id,30),scales = "free_x")+
  theme_bw()+
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=16),
        axis.text.y = element_text(size=14),
        strip.text.x = element_text(angle = 90,size = 16,hjust = 0),
        legend.text = element_text(size = 14),
        legend.title = element_text(size=16),
        strip.background.x = element_blank())

# run function-ancom first to make plot
comb_p<-cowplot::plot_grid(rel.p,ancom_p,ncol = 1,align = "v")
ggsave(plot = comb_p,filename = paste0("F2_functional_stats/top.sp_ancom_0.9.png"),dpi = 300,width = 14,height = 16,units = "in")

# individual heatmaps
anno<-samples %>% column_to_rownames("SampleID") %>% as.data.frame()

for (i in unique(ft$taxa_id)) {
  p<-ft %>% filter(taxa_id==i) %>% 
    select(species,SampleID,abundance) %>% 
    pivot_wider(id_cols = species,names_from = SampleID,values_from = abundance) %>% 
    filter(!is.na(species)) %>% 
    column_to_rownames("species") %>%
    # as.matrix() %>% 
    pheatmap::pheatmap(annotation_col = anno,annotation_colors = cols,scale = "none")
  n<-i %>% str_split(.,pattern = ":",n = 1)
  ggsave(p,filename = paste0("F2_functional_stats/",n,"_ancom_0.9.png"),dpi = 300,width = 12,height = 8,units = "in")
}

# compare the top n taxa for specific go terms, betweem subspecies (?)
for (i in unique(ft$taxa_id)) {
  p<-ft %>% filter(taxa_id==i) %>% 
    filter(!is.na(species)) %>% 
    group_by(species) %>% 
    mutate(max_abundance=max(abundance,na.rm=T)) %>% ungroup() %>% 
    slice_max(order_by = max_abundance,prop = 0.1) %>% 
    ggplot(aes(x=species,fill=Spec.subspecies,y=abundance))+
    geom_boxplot()+theme_classic()
  ggsave(p,filename = paste0("F2_functional_stats/top10p-",i,"_ancom_0.9.png"),dpi = 300,width = 8,height = 8,units = "in")
}

# instead, why not grab the taxa with the highest average abundance across all samples and GO terms?
top.sp<-ft %>% filter(taxa_id==i & !is.na(species)) %>% 
  group_by(species,Spec.subspecies) %>% 
  summarise(mean_abundance=mean(abundance,na.rm=T),
            var=var(abundance,na.rm=T)) %>% group_by(Spec.subspecies) %>% 
  slice_max(order_by = mean_abundance,n=5) %>% ungroup() %>% distinct(species) %>% pull(species)

for (i in unique(ft$taxa_id)) {
  p<-ft %>% filter(taxa_id==i) %>%
    filter(species%in% top.sp) %>% 
    select(species,SampleID,abundance) %>% 
    pivot_wider(id_cols = species,names_from = SampleID,values_from = abundance) %>% 
    filter(!is.na(species)) %>% 
    column_to_rownames("species") %>%
    # as.matrix() %>% 
    pheatmap::pheatmap(annotation_col = anno,annotation_colors = cols,scale = "none")
  n<-i %>% str_split(.,pattern = ":",n = 1)
  ggsave(p,filename = paste0("F2_functional_stats/top.sp_",n,"_ancom_0.9.png"),dpi = 300,width = 12,height = 8,units = "in")
}