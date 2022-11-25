library(foreach)
library(tidyverse)
library(ggtree)
library(treeio)
library(phangorn)
library(ggsci)
library(RColorBrewer)

# panaroo output:
# gene_presence_absence.csv
# potential pseudo genes will end with the suffix '_stop'
# genes with unusual lengths will end with the suffix '_len'
# fragmented genes will include multiple sequence IDs seperated by a semicolon ';'

# nitrate genes
genes<-c("narX","narG","narY","narK2","norB","modA","moaB","moaE","mog","mopII","moeA","moeB","aniA")

tree<-read.tree("M3_trees/data/phylophlan/input_s__Neisseria_gonorrhoeae/input.tre")
tree<-drop.tip(tree,c("BS005_mappedcontigs.sorted","LIB007.A0117_mappedcontigs.sorted","BS001_mappedcontigs.sorted","LIB013.A0101_mappedcontigs.sorted","LIB014.A0101_mappedcontigs.sorted","LIB012.A0101_mappedcontigs.sorted","ERR2503700_mappedcontigs.sorted","ERR2868193_mappedcontigs.sorted","EXB023.A0101_mappedcontigs.sorted","EXB015.A2101_mappedcontigs.sorted","EXB021.A0101_mappedcontigs.sorted","EXB015.A2501_mappedcontigs.sorted"))
refs<-tree$tip.label[!grepl("sorted",tree$tip.label)]

tree$tip.label<-gsub("_mappedcontigs.sorted","",tree$tip.label)
tree$tip.label<-gsub("(.+?_.+?)_.*" ,"\\1",tree$tip.label)
plot(tree)

panaroo<-foreach(i=genes,.combine = rbind) %do% {
  read_csv("prodigal_training_on_ref/Neisseria_gene_presence_absence.csv",col_types = "c") %>% 
    filter(grepl(i,`Non-unique Gene name`) | grepl(i,Gene) ) %>% 
    mutate(geneid=i,
           issues=ifelse(grepl("_stop",`Non-unique Gene name`) | grepl("_stop",Gene),"pseudo gene",
                         ifelse(grepl("_len",`Non-unique Gene name`) | grepl("_len",Gene),"unusual length",
                                ifelse(grepl(";",`Non-unique Gene name`) | grepl(";",Gene),"unusual length",NA)))) %>% 
    select(-issues) %>% # do something about pseudogenes here?
    pivot_longer(cols = -c(Gene,`Non-unique Gene name`,Annotation,geneid),names_to = "sample",values_to = "sequence_name") %>% 
    mutate(presence.absence=ifelse(!is.na(sequence_name),"present","absent")) %>% 
    select(geneid,sample,presence.absence) %>% 
    filter(!grepl("bin",sample))
}

# any issues with these genes?
table(panaroo$issues)

blast <- read_csv("M4_pangenome/data/tblastn/all-Neisseria.csv",skip=1,col_names=c("file.id","qseqid","sseqid","stitle","pident","qcovs","length","mismatch","gapopen","qstart","qend","sstart","send","qframe","sframe","frames","evalue","bitscore","qseq","sseq","type")) %>%
  mutate(file.id=gsub("Neisseria/|_mappedcontigs.sorted","",gsub("(.+?_.+?)_.*" ,"\\1",file.id))) %>%
  group_by(qseqid) %>%
  mutate(bin.pident=pident[file.id=="genomic_bin"],
         bin.length=length[file.id=="genomic_bin"],
         bin.sseq.length=str_length(str_replace_all(sseq, "[^[:alnum:]]", ""))) %>%
  group_by(file.id,qseqid) %>%
  mutate(presence.absence=ifelse((length/bin.length)>=0.8,"present","absent")) %>%
  filter(!grepl("genomic_bin",file.id)) %>% 
  select(qseqid,file.id,presence.absence)
  
genotype_mat<-panaroo %>% 
  mutate(geneid=fct_recode(geneid,"mogA/moaB"="moaB")) %>% # recode gene names that don't match
  mutate(sample=gsub("Rothia/|_mappedcontigs.sorted.fa","",gsub("(.+?_.+?)_.*" ,"\\1",sample))) %>%
  full_join(blast,by=c("geneid"="qseqid","sample"="file.id","presence.absence")) %>% 
  pivot_wider(id_cols = sample,names_from = geneid,values_from = presence.absence,
              values_fn = function(x) ifelse(any(x=="present"),"present","absent")) %>%
  column_to_rownames("sample") %>% 
  mutate_all(~replace_na(.,"absent")) %>% 
  select(-narY,-mopII) # getting rid of non-essential genes

iso<-read.delim("isolation_sources/neisseria-refs-isolation-sources.csv",header=F,sep = ",")

taxa<-read.table("phylophlan/neisseria-refs-taxaids.txt",colClasses = c("character","character"),header=F)
nrow(taxa) == length(refs)

tree_taxa<-foreach(i=unique(taxa$V2),.combine = rbind,.packages = "taxize") %do% {
  ncbi_get_taxon_summary(id = i)
}

taxa_merge<-tree_taxa %>% 
  right_join(taxa,by=c("uid"="V2")) %>% 
  right_join(iso,by="V1") %>% 
  rename("accession"="V1","isolation.source"="V2") %>% 
  bind_rows(.,bind_cols(uid="unknown",
                        name="",
                        rank="unknown",
                        V3="Gorilla",V4="wild",
                        accession=tree$tip.label[!grepl("_",tree$tip.label)])) %>% 
  mutate(name=gsub("Neisseria","N.",name),
         lab=ifelse(name=="",accession,ifelse(grepl("Eikenella|Snodgrassella",name),name,paste0(name,", ",V4)))) %>% 
  distinct(uid,accession,.keep_all = T) %>% 
  relocate(accession)

n.tree<-ggtree(phangorn::midpoint(drop.tip(tree,tip = c("GCF_000600005.1","GCF_900187105.1","GCF_014054965.1")))) %<+% taxa_merge +
  # geom_nodelab(nudge_x = -0.02,nudge_y = 0.5,size=3)+
  geom_cladelabel(node=58, label="This Study",align=T,fontsize = 5,offset = 0.65,offset.text = 0.02) +
  # geom_text2(aes(subset=!isTip, label=node), hjust=-.5)+
  geom_tiplab(geom="text",aes(label=name),offset = 0.85,align = T,size=6,linesize = 0)+
  # geom_tippoint(aes(color=V3),size=3,position=position_nudge(x = 0.01))+
  vexpand(ratio = 0.2)+
  hexpand(ratio = 0.25)+
  # vexpand(ratio = 0.25)+
  geom_treescale(y = -1)+
  theme(legend.text = element_text(size = 14),legend.title = element_text(size = 16),legend.position = c(0.25,0.75))

n.p<-gheatmap(n.tree,genotype_mat,colnames=TRUE,high="black",low="gray",legend_title="",colnames_angle = 90,colnames_position = "top",hjust = 0,font.size = 4)

ggsave(filename = "M4_pangenome/Neisseria-nitrate-tree-pident.png",plot = n.p,dpi=300,width = 12,height = 12,units = "in")