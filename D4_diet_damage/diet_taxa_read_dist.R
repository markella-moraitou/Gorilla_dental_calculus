library(phyloseq)
library(tidyverse)
library(microViz)
library(cowplot)

setwd("Gorilla_dental_calculus_zenodo")

euk_genus_diet<-readRDS("D3_diet_stats/euk_genus_diet.rds")

euk_genus_diet@tax_table<- euk_genus_diet@tax_table %>% data.frame() %>% select(-"clade",-"species") %>% as.matrix() %>%  tax_table()
cat(taxa_names(euk_genus_diet),sep = "\n")

# because we have censored at <= 10 reads, it makes it hard to generate averages of read counts, per taxon, per sample
# use the filtered, decontaminated data to subset the unfiltered, un-decontaminated data

euk_genus_diet_allreads<-prune_samples(sample_names(euk_genus_diet),prune_taxa(taxa_names(euk_genus_diet),euk_genus))

# average number of reads per dietary taxon, per sample
ps_melt(euk_genus_diet_allreads) %>% 
  filter(Abundance>10) %>% 
  group_by(OTU,Sample) %>% 
  summarise(mean=mean(Abundance,na.rm=T)) %>%
  ungroup() %>% 
  summarise(mean=mean(mean,na.rm=T))
  
ps_melt(euk_genus_diet_allreads) %>% 
  filter(Abundance>10) %>% 
  group_by(Sample,OTU) %>% 
  summarise(mean=mean(Abundance,na.rm=T)) %>% 
  ggplot(aes(x=mean))+
  geom_histogram(fill="gray",color="black")+
  geom_vline(xintercept = 191,size=2,color="black")+
  labs(x="Log10 Average reads per sample per taxon")+
  scale_x_log10()+
  see::theme_modern()

# number of reads for the highest abundance dietary taxon in any sample
microbiome::abundances(euk_genus_diet) %>% data.frame() %>% 
  rownames_to_column("Genus") %>% 
  pivot_longer(cols = -Genus,names_to = "SampleID",values_to = "abundance") %>% 
  slice_max(order_by = abundance,n = 10)
  
# across all taxa, abundance values
microbiome::abundances(euk_genus_diet) %>% data.frame() %>% 
  rownames_to_column("Genus") %>% 
  pivot_longer(cols = -Genus,names_to = "SampleID",values_to = "abundance") %>% 
  group_by(SampleID,Genus) %>% 
  summarise(mean=mean(abundance,na.rm=T)) %>% 
  ggplot(aes(x=mean))+
  geom_histogram()+
  scale_x_log10()+
  geom_vline(xintercept=500,size=2,color="black")

# what are the most abundant dietary taxa?
microbiomeutilities::plot_taxa_heatmap(euk_genus_diet,subset.top = 10,transformation = "clr",VariableA = "Spec.subspecies")

# Palisota
palisota<-read_delim("all.gorilla.palisota.kraken.reads.txt",col_names = F,delim = "\t") %>% 
  separate(col = "X1",into = c("Sample","perc_reads_clade"),sep = ":",remove = T) %>% 
  rename("n_reads_clade"=3,"n_reads"=4,"rank"=5,"taxa_id"=6,"name"=7) %>% 
  mutate(Sample=gsub(pattern = "../MARKELLA/D2_kraken2_full_db/",replacement = "",x = gsub(pattern = "_m_bact_arch_vir_removed_kraken2_report.txt",replacement = "",x=Sample))) %>% 
  filter(Sample!="" & Sample %in% sample_names(euk_genus_diet)) %>% 
  group_by(Sample) %>% 
  summarise(n_reads_genus=sum(n_reads,na.rm = T))

# phyllostachys
phyllostachys<-read_delim("all.gorilla.phyllostachys.kraken.reads.txt",col_names = F,delim = "\t") %>% 
  separate(col = "X1",into = c("Sample","perc_reads_clade"),sep = ":",remove = T) %>% 
  rename("n_reads_clade"=3,"n_reads"=4,"rank"=5,"taxa_id"=6,"name"=7) %>% 
  mutate(Sample=gsub(pattern = "../MARKELLA/D2_kraken2_full_db/",replacement = "",x = gsub(pattern = "_m_bact_arch_vir_removed_kraken2_report.txt",replacement = "",x=Sample))) %>% 
  filter(Sample!="" & Sample %in% sample_names(euk_genus_diet)) %>% 
  group_by(Sample) %>% 
  summarise(n_reads_genus=sum(n_reads,na.rm = T)) %>% 
  mutate(genus="Phyllostachys")

# galium
galium<-read_delim("all.gorilla.galium.kraken.reads.txt",col_names = F,delim = "\t") %>% 
  separate(col = "X1",into = c("Sample","perc_reads_clade"),sep = ":",remove = T) %>% 
  rename("n_reads_clade"=3,"n_reads"=4,"rank"=5,"taxa_id"=6,"name"=7) %>% 
  mutate(Sample=gsub(pattern = "../MARKELLA/D2_kraken2_full_db/",replacement = "",x = gsub(pattern = "_m_bact_arch_vir_removed_kraken2_report.txt",replacement = "",x=Sample))) %>% 
  filter(Sample!="" & Sample %in% sample_names(euk_genus_diet)) %>% 
  group_by(Sample) %>% 
  summarise(n_reads_genus=sum(n_reads,na.rm = T)) %>% 
  mutate(genus="Galium")

# rubia
rubia<-read_delim("all.gorilla.rubia.kraken.reads.txt",col_names = F,delim = "\t") %>% 
  separate(col = "X1",into = c("Sample","perc_reads_clade"),sep = ":",remove = T) %>% 
  rename("n_reads_clade"=3,"n_reads"=4,"rank"=5,"taxa_id"=6,"name"=7) %>% 
  mutate(Sample=gsub(pattern = "../MARKELLA/D2_kraken2_full_db/",replacement = "",x = gsub(pattern = "_m_bact_arch_vir_removed_kraken2_report.txt",replacement = "",x=Sample))) %>% 
  filter(Sample!="" & Sample %in% sample_names(euk_genus_diet)) %>% 
  group_by(Sample) %>% 
  summarise(n_reads_genus=sum(n_reads,na.rm = T)) %>% 
  mutate(genus="Rubia")

# galium
pg<-galium %>% 
  ggplot(aes(x=n_reads_genus))+
  geom_histogram(bins=15,fill="gray",color="black")+
  scale_x_log10()+
  geom_vline(xintercept=500,size=2,color="black")+
  # facet_grid(genus~.,scales = "free")+
  labs(y="Number of Samples",x=expression(paste("Log10 reads assigned to",italic(" Galium"))))+
  see::theme_modern()+theme(legend.position = "none",axis.text = element_text(size=20),
                            axis.title = element_text(size=20))
ggsave(plot = pg,filename = "gallium-reads.png",dpi = 300,scale = 2)

# add in mapdamage plot
mapd<-read_delim("misincorporation.txt",skip = 3) %>% 
  filter(Pos %in% 1:25) %>% 
  select(Pos,End,`C>T`,`G>A`) %>% 
  pivot_longer(cols = -c(Pos,End),names_to = "transition",values_to = "freq" ) %>% 
  mutate(Pos=ifelse(End=="3p",Pos*-1,Pos),
         End=factor(End,levels = c("5p","3p")),
         freq=freq/200)

p5p<-mapd %>% 
  filter(End=="5p") %>% 
  ggplot(aes(x=Pos,y=freq,color=transition))+
  geom_line()+
  coord_cartesian(xlim = c(1,25),ylim = c(0,0.15))+
  scale_color_manual(values=c("red","blue"))+
  scale_x_continuous(limits=c(1,25),breaks=1:25)+
  labs(y="Frequency")+
  see::theme_modern()+
  theme(legend.position = "none",
        axis.text.y = element_text(size=16),
        axis.text.x = element_text(size=8),
        axis.title = element_text(size=16),
        axis.title.x = element_blank(),
        axis.ticks = element_line())
  
p3p<-mapd %>% 
  filter(End=="3p") %>% 
  ggplot(aes(x=Pos,y=freq,color=transition))+
  geom_line()+
  coord_cartesian(xlim = c(-25,-1),ylim = c(0,0.15))+
  scale_color_manual(values=c("red","blue"))+
  scale_x_continuous(limits=c(-25,-1),breaks=-25:-1)+
  scale_y_continuous(position="right")+
  see::theme_modern()+
  theme(legend.position = "none",
        axis.text.y = element_text(size=16),
        axis.text.x = element_text(size=8),
        axis.title.x = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks = element_line())

plot_grid(pg+theme(axis.text = element_text(size=16),
                   axis.title = element_text(size=16),
                   axis.ticks = element_line()),
          NULL,
          plot_grid(p5p,p3p,nrow = 1),
          # align = "h",
          nrow = 3,
          scale = c(0.75,0,1),
          rel_heights = c(1.5,-0.2,1),
          rel_widths = c(1,0,2),
          labels = c("a)","","b)"),label_size = 16)
ggsave(plot = last_plot(),filename = "gallium-reads-damage-plot.png",dpi = 300,width = 12,height = 8,units = "in")
