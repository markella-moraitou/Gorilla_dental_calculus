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
# genes<-c("narG","narY","narH","narI","narJ","narW","nirB","nirD","nrfA","narX1","narX2","narK2","modA","moaB","modA","moaE2","mog","mopll","moeA","moe","nifX","hmp")
genes<-c("narX","narG","narY","narK2","norB","modA","moaB","moaE","mog","mopII","moeA","moeB","aniA")

tree<-read.tree("calculus/species-specific/DC-gorilla/mag/phylophlan/Neisseria/input_s__Neisseria_gonorrhoeae/input.tre")
tree<-drop.tip(tree,c("BS005_mappedcontigs.sorted","LIB007.A0117_mappedcontigs.sorted","BS001_mappedcontigs.sorted","LIB013.A0101_mappedcontigs.sorted","LIB014.A0101_mappedcontigs.sorted","LIB012.A0101_mappedcontigs.sorted","ERR2503700_mappedcontigs.sorted","ERR2868193_mappedcontigs.sorted","EXB023.A0101_mappedcontigs.sorted","EXB015.A2101_mappedcontigs.sorted","EXB021.A0101_mappedcontigs.sorted","EXB015.A2501_mappedcontigs.sorted"))
refs<-tree$tip.label[!grepl("sorted",tree$tip.label)]

tree$tip.label<-gsub("_mappedcontigs.sorted","",tree$tip.label)
tree$tip.label<-gsub("(.+?_.+?)_.*" ,"\\1",tree$tip.label)
plot(tree)

panaroo<-foreach(i=genes,.combine = rbind) %do% {
  read_csv("calculus/species-specific/DC-gorilla/mag/panaroo/prodigal_training_on_ref/Neisseria_gene_presence_absence.csv",col_types = "c") %>% 
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

blast <- read_csv("calculus/species-specific/DC-gorilla/mag/data/all-Neisseria.csv",skip=1,col_names=c("file.id","qseqid","sseqid","stitle","pident","qcovs","length","mismatch","gapopen","qstart","qend","sstart","send","qframe","sframe","frames","evalue","bitscore","qseq","sseq","type")) %>%
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

iso<-read.delim("calculus/species-specific/DC-gorilla/mag/data/neisseria-refs-isolation-sources.csv",header=F,sep = ",")

taxa<-read.table("calculus/species-specific/DC-gorilla/mag/data/neisseria-refs-taxaids.txt",colClasses = c("character","character"),header=F)
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

ggsave(filename = "calculus/species-specific/DC-gorilla/mag/images/Neisseria-nitrate-tree-pident.png",plot = n.p,dpi=300,width = 12,height = 12,units = "in")

###########

tree<-read.tree("calculus/species-specific/DC-gorilla/mag/phylophlan/Rothia/output/input.tre")
tree<-drop.tip(tree,c("BS005_mappedcontigs.sorted","LIB007.A0117_mappedcontigs.sorted","BS001_mappedcontigs.sorted","LIB013.A0101_mappedcontigs.sorted","LIB014.A0101_mappedcontigs.sorted","LIB012.A0101_mappedcontigs.sorted","ERR2503700_mappedcontigs.sorted","ERR2868193_mappedcontigs.sorted","EXB023.A0101_mappedcontigs.sorted","EXB015.A2101_mappedcontigs.sorted","EXB021.A0101_mappedcontigs.sorted","EXB015.A2501_mappedcontigs.sorted"))
refs<-tree$tip.label[!grepl("sorted",tree$tip.label)]

tree$tip.label<-gsub("_mappedcontigs.sorted","",tree$tip.label)
tree$tip.label<-gsub("(.+?_.+?)_.*" ,"\\1",tree$tip.label)
plot(tree)

panaroo<-foreach(i=genes,.combine = rbind) %do% {
  read_csv("calculus/species-specific/DC-gorilla/mag/panaroo/prodigal_training_on_ref/Rothia_gene_presence_absence.csv",col_types = "c") %>% 
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

blast <- read_csv("calculus/species-specific/DC-gorilla/mag/data/all-Rothia.csv",skip=1,col_names=c("qseqid","sseqid","stitle","pident","qcovs","length","mismatch","gapopen","qstart","qend","sstart","send","qframe","sframe","frames","evalue","bitscore","qseq","sseq","file.id","type")) %>%
  mutate(file.id=gsub("Rothia/|_mappedcontigs.sorted","",gsub("(.+?_.+?)_.*" ,"\\1",file.id))) %>%
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

taxa<-read.table("calculus/species-specific/DC-gorilla/mag/data/rothia-refs-taxaids.txt",colClasses = c("character","character"),header=F)

tree_taxa<-foreach(i=unique(taxa$V2),.combine = rbind,.packages = "taxize") %do% {
  ncbi_get_taxon_summary(id = i)
}

taxa_merge<-tree_taxa %>% 
  right_join(taxa,by=c("uid"="V2")) %>% 
  rename("accession"="V1") %>% 
  bind_rows(.,bind_cols(uid="unknown",
                        name="This study",
                        rank="unknown",
                        V3="Gorilla",V4="wild",
                        accession=tree$tip.label[!grepl("_",tree$tip.label)])) %>% 
  mutate(name=gsub("Rothia","R.",name),
         lab=ifelse(name=="This study",accession,name)) %>% 
  distinct(uid,accession,.keep_all = T) %>% 
  relocate(accession)

r.tree<-ggtree(phangorn::midpoint(tree)) %<+% taxa_merge +
  # geom_nodelab(nudge_x = -0.02,nudge_y = 0.5,size=3)+
  # geom_text2(aes(subset=!isTip, label=node), hjust=-.5)+
  # geom_tiplab(geom="text",aes(label=lab),offset = 0.5,align = F,size=4,linesize = 0)+
  geom_cladelabel(node = 293, label="This Study",align=T,fontsize = 7,offset = 0.85,offset.text = 0.02) +
  geom_cladelabel(node = 246,label = "R. mucilaginosa",align = T,fontsize = 7,offset = 0.85)+
  geom_cladelabel(node = 216,label = "R. kristinae",align = T,fontsize = 7,offset = 0.85)+
  geom_cladelabel(node = 234,label = "R. nasimurium",align = T,fontsize = 7,offset = 0.85)+
  geom_cladelabel(node = 239,label = "R. terrae",align = T,fontsize = 7,offset = 0.85)+
  geom_cladelabel(node = 270,label = "R. aeria",align = T,fontsize = 7,offset = 0.85)+
  geom_cladelabel(node = 279,label = "R. dentocariosa",align = T,fontsize = 7,offset = 0.85)+
  geom_cladelabel(node = 243,label = "R. amarae",align = T,fontsize = 7,offset = 0.85)+
  geom_cladelabel(node = 231,label = "R. koreensis",align = T,fontsize = 7,offset = 0.85)+
  # geom_cladelabel(node = 208,label = "Other Micrococcaceae species",align = T,fontsize = 7,offset = 0.85)+
  # geom_tippoint(aes(color=V3),size=3,position=position_nudge(x = 0.01))+
  # vexpand(ratio = 0.2)+
  hexpand(ratio = 0.25)+
  geom_treescale(y = -1)+
  theme(legend.text = element_text(size = 14),legend.title = element_text(size = 16),legend.position = c(0.25,0.75)) 
  
r.tree.sub<-r.tree %>%
  ggtree::collapse(node=208)

r.p<-gheatmap(r.tree.sub,genotype_mat,legend_title="",colnames_angle = 90,colnames_position = "top",hjust = 0,font.size = 8)

ggsave(filename = "calculus/species-specific/DC-gorilla/mag/images/Rothia-nitrate-tree-panaroo.png",plot = r.p,dpi=300,width = 12,height = 12,units = "in")

#####

tree<-read.newick("calculus/species-specific/DC-gorilla/mag/phylophlan/Veillonella/output/input_s__Veillonella_dispar/input.tre")
tree<-drop.tip(tree,c("BS005_mappedcontigs.sorted","LIB007.A0117_mappedcontigs.sorted","BS001_mappedcontigs.sorted","LIB013.A0101_mappedcontigs.sorted","LIB014.A0101_mappedcontigs.sorted","LIB012.A0101_mappedcontigs.sorted","ERR2503700_mappedcontigs.sorted","ERR2868193_mappedcontigs.sorted","EXB023.A0101_mappedcontigs.sorted","EXB015.A2101_mappedcontigs.sorted","EXB021.A0101_mappedcontigs.sorted","EXB015.A2501_mappedcontigs.sorted"))
refs<-tree$tip.label[!grepl("sorted",tree$tip.label)]

tree$tip.label<-gsub("_mappedcontigs.sorted","",tree$tip.label)
tree$tip.label<-gsub("(.+?_.+?)_.*" ,"\\1",tree$tip.label)

tree$node.label<-ifelse(as.numeric(tree$node.label)>=0.75,as.character(round(as.numeric(tree$node.label),digits = 2)),"")

panaroo<-foreach(i=genes,.combine = rbind) %do% {
  read_csv("calculus/species-specific/DC-gorilla/mag/panaroo/prodigal_training_on_ref/Veillonella_gene_presence_absence.csv",col_types = "c") %>% 
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

blast <- read_csv("calculus/species-specific/DC-gorilla/mag/data/all-Veillonella.csv",skip=1,col_names=c("qseqid","sseqid","stitle","pident","qcovs","length","mismatch","gapopen","qstart","qend","sstart","send","qframe","sframe","frames","evalue","bitscore","qseq","sseq","file.id","type")) %>%
  mutate(file.id=gsub("Veillonella/|_mappedcontigs.sorted","",gsub("(.+?_.+?)_.*" ,"\\1",file.id))) %>%
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
  mutate(sample=gsub("Veillonella/|_mappedcontigs.sorted.fa","",gsub("(.+?_.+?)_.*" ,"\\1",sample))) %>%
  full_join(blast,by=c("geneid"="qseqid","sample"="file.id","presence.absence")) %>% 
  pivot_wider(id_cols = sample,names_from = geneid,values_from = presence.absence,
              values_fn = function(x) ifelse(any(x=="present"),"present","absent")) %>%
  column_to_rownames("sample") %>% 
  mutate_all(~replace_na(.,"absent")) %>% 
  select(-narY,-mopII) # getting rid of non-essential genes

refs<-tree$tip.label[grepl("_",tree$tip.label)]
write.table(refs,file = "calculus/species-specific/DC-gorilla/mag/data/veillonella-refs.txt",row.names = F,quote = F,col.names = F)

# system('refs=$(cat M3_trees/veillonella-refs.txt);
#        for line in $refs; do res=$(esearch -query $line -db assembly | elink -target biosample | esummary | xtract -pattern DocumentSummary -first Title -element Accession -group Attribute -if Attribute@harmonized_name -equals "isolation_source" -element Attribute);
#        echo $line $res  >> M3_trees/veillonella-refs-isolation-sources.raw.txt;
#        done;')

# edited by hand here...
taxa<-read.delim("calculus/species-specific/DC-gorilla/mag/data/veillonella-refs-taxaids.txt",header=F,sep = "\t",colClasses = c("factor","character"))
nrow(taxa) == length(refs)

tree_taxa<-foreach(i=unique(taxa$V2),.combine = rbind,.packages = "taxize") %do% {
  ncbi_get_taxon_summary(id = i)
}

iso<-read.csv("calculus/species-specific/DC-gorilla/mag/data/veillonella-refs-isolation-sources.csv",header=F)

taxa_merge<-tree_taxa %>% 
  right_join(taxa,by=c("uid"="V2")) %>% 
  right_join(iso,by="V1") %>% 
  rename("accession"="V1","isolation.source"="V2") %>% 
  bind_rows(.,bind_cols(uid="unknown",
                        name="MAG Constructed\nin this study",
                        rank="unknown",
                        accession=tree$tip.label[!grepl("_",tree$tip.label)],
                        isolation.source="Gorilla")) %>% 
  distinct(uid,accession,.keep_all = T) %>% 
  relocate(accession)

v.tree<-ggtree(phangorn::midpoint(tree)) %<+% taxa_merge + 
  # geom_tiplab(align=TRUE, linetype = "dotted", linesize = .5,geom = "text",size=0)+
  # geom_tiplab(geom="text",aes(label=name)) +
  # geom_text2(aes(subset=!isTip, label=node), hjust=-.5)+
  geom_cladelabel(node = 179,label = "V. parvula",align = T,offset = 2.25,offset.text = 0.1,fontsize = 7)+
  geom_cladelabel(node = 281,label = "V. dispar /\nV. nakazawae",align = T,offset = 2.25,offset.text = 0.1,fontsize = 7)+
  geom_cladelabel(node = 190,label = "V. atypica",align = T,offset = 2.25,offset.text = 0.1,fontsize = 7)+
  geom_cladelabel(node = 270,label = "V. ratti",align = T,offset = 2.25,offset.text = 0.1,fontsize = 7)+
  geom_cladelabel(node = 207,label = "V. sp.",align = T,offset = 2.25,offset.text = 0.1,fontsize = 7)+
  geom_cladelabel(node = 224,label = "Veillonella MAGs\n(this study)",align = T,offset = 2.25,offset.text = 0.1,fontsize = 7)+
  # geom_nodelab(geom = "text",hjust = -0.3)+
  # geom_tippoint(aes(fill=isolation.source),size=3,shape=21)+
  # scale_fill_manual( values = c(ggpubr::get_palette(palette = "Paired",n_distinct(iso$V2)),"gray"))+
  geom_treescale(y = 5)+
  # guides(color=guide_legend("Taxonomy"))+
  vexpand(ratio = 0.2)+
  hexpand(ratio = 0.5)

v.p<-gheatmap(v.tree,genotype_mat,colnames=TRUE,high="black",low="gray",legend_title="",colnames_angle = 90,colnames_position = "top",hjust = 0,font.size = 4)

ggsave(filename = "calculus/species-specific/DC-gorilla/mag/images/Veillonella-nitrate-tree-pident.png",plot = v.p,dpi=300,width = 12,height = 12,units = "in")
