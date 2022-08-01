library(tidyverse)
library(ggpubr)
library(ggsignif)
library(plotly)
library(microbiome)
library(pheatmap)
library(ggtree)
library(treeio)
library(taxize)
library(foreach)
library(phangorn)
library(ggsci)
library(RColorBrewer)

# NEW FIGURE FOR MAG TREES, MAG ABUNDANCE, AND NITRATE GENES
# - 1 FIGURE FOR JUST NEISSERIA
# - THIS CONTAINS TREE, ABUNDANCE, AND NITRATE GENE PANELS
# 2 NEW SUPPLEMENTAL FIGURES
# - 1 FIGURE FOR ROTHIA/VEILLONELLA
# - THESE CONTAINS TREE, ABUNDANCE, AND NITRATE GENE PANELS
# SEPARATE SUPPLEMENTARY FIGURE FOR L. GORILLAE MAG TREE

# start with phylogenies and nitrate genes for all species first ----------------------------

# panaroo output:
# gene_presence_absence.csv
# potential pseudo genes will end with the suffix '_stop'
# genes with unusual lengths will end with the suffix '_len'
# fragmented genes will include multiple sequence IDs seperated by a semicolon ';'

# nitrate genes
# genes<-c("narG","narY","narH","narI","narJ","narW","nirB","nirD","nrfA","narX1","narX2","narK2","modA","moaB","modA","moaE2","mog","mopll","moeA","moe","nifX","hmp")
genes<-c("narX","narG","narY","narK2","norB","modA","moaB","moaE","mog","mopII","moeA","moeB","aniA")

iso.pal<-c(get_palette(palette = "Set1", 5),"gray")

#### NEISSERIA
tree<-read.tree("calculus/species-specific/DC-gorilla/mag/phylophlan/Neisseria/input_s__Neisseria_gonorrhoeae/input.tre")

# clean up tips
tree$tip.label<-gsub("^.*\\|", "", tree$tip.label)
tree$tip.label<-sub("(.*?_.*?)_.*", "\\1", tree$tip.label)

tree$node.label[c(1:18,49:50)]<-""
tree$node.label<-as.character(round(as.numeric(tree$node.label),3))
tree$node.label<-ifelse(tree$node.label=="1","",tree$node.label)

tree<-drop.tip(tree,tree$tip.label[grepl("sorted",tree$tip.label)][-1])
tree$tip.label[grepl("sorted",tree$tip.label)]<-"This Study"
plot(tree)

refs<-tree$tip.label[!grepl("sorted",tree$tip.label)]
write.table(refs,file = "calculus/species-specific/DC-gorilla/mag/data/neisseria-refs.txt",row.names = F,quote = F,col.names = F)

# system('refs=$(cat calculus/species-specific/DC-gorilla/mag/data/neisseria-refs.txt);
#        for line in $refs; do res=$(esearch -query $line -db assembly | elink -target biosample | esummary | xtract -pattern DocumentSummary -first Title -element Accession -group Attribute -if Attribute@harmonized_name -equals "isolation_source" -element Attribute);
#        echo $line $res  >> calculus/species-specific/DC-gorilla/mag/data/neisseria-refs-isolation-sources.raw.txt;
#        done;')


# system('refs=$(cat neisseria-refs.txt);
#        for line in $refs; do res=$(esearch -query $line -db assembly | elink -target biosample | esummary | xtract -pattern DocumentSummary -element Taxonomy);
#        echo $line $res  >> neisseria-refs-taxaids.txt;
#        done;')

iso<-read.delim("calculus/species-specific/DC-gorilla/mag/data/neisseria-refs-isolation-sources.csv",header=F,sep = ",")

taxa<-read.table("calculus/species-specific/DC-gorilla/mag/data/neisseria-refs-taxaids.txt",colClasses = c("character","character"),header=F)
nrow(taxa) == length(refs)

tree_taxa<-foreach(i=unique(taxa$V2),.combine = rbind,.packages = "taxize") %do% {
  print(i)
  ncbi_get_taxon_summary(id = i)
}

taxa_merge<-tree_taxa %>% 
  right_join(taxa,by=c("uid"="V2")) %>% 
  right_join(iso,by="V1") %>% 
  rename("accession"="V1","isolation.source"="V2") %>% 
  bind_rows(.,bind_cols(uid="unknown",
                        name="This Study",
                        rank="unknown",
                        V3="Gorilla",V4="wild",
                        accession="This Study")) %>% # tree$tip.label[grepl("sorted",tree$tip.label)])) 
  mutate(name=gsub("Neisseria","N.",name),
         lab=ifelse(!grepl("Study",name),ifelse(grepl("Eikenella|Snodgrassella",name),name,paste0(name,", ",V4)),name)) %>% 
  distinct(uid,accession,.keep_all = T) %>% 
  relocate(accession) %>% 
  mutate(V3=fct_relevel(V3,c("Gorilla","Human","Avian","Mammalian","Reptile")))

ggtree(phangorn::midpoint(drop.tip(tree,tip = c("GCF_900187105.1","GCF_014054965.1")))) %<+% taxa_merge +
  geom_text2(aes(subset=!isTip, label=node), hjust=-.5)

n.tree<-ggtree(phangorn::midpoint(drop.tip(tree,tip = c("GCF_900187105.1","GCF_014054965.1")))) %<+% taxa_merge +
  geom_nodelab(nudge_x = -0.0375,nudge_y = 0.5,size=4,family = "calibri")+
  # geom_cladelabel(node=59, label="This Study",align=T,fontsize = 5,offset = 0.02,offset.text = 0.02) +
  # geom_text2(aes(subset=!isTip, label=node), hjust=-.5)+
  geom_tiplab(geom="text",aes(label=lab),offset = 0.85,align = T,size=6,family = "calibri")+
  geom_tippoint(aes(color=V3),size=3,position=position_nudge(x = 0.01))+
  geom_tippoint(aes(subset=(isTip & node==1)),color=iso.pal[1],size=6,position=position_nudge(x = 0.01))+
  scale_color_manual( values = iso.pal)+
  guides(color=guide_legend("Host/Isolation Source"))+
  vexpand(ratio = 0.13)+
  hexpand(ratio = 0.25)+
  # vexpand(ratio = 0.25)+
  geom_treescale(y = -1,fontsize = 5,family = "calibri")+
  theme(legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.position = "bottom",
        text = element_text(family = "calibri"))

ggsave(filename = "calculus/species-specific/DC-gorilla/mag/images/neisseria-tree.png",plot = n.tree,dpi=300,width = 12,height = 12,units = "in")

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
  mutate(file.id=gsub("Neisseria/|.fa","",gsub("(.+?_.+?)_.*" ,"\\1",file.id))) %>%
  group_by(qseqid) %>%
  mutate(bin.pident=pident[file.id=="genomic_bin"],
         bin.length=length[file.id=="genomic_bin"],
         bin.sseq.length=str_length(str_replace_all(sseq, "[^[:alnum:]]", ""))) %>%
  group_by(file.id,qseqid) %>%
  mutate(presence.absence=ifelse((length/bin.length)>=0.8,"present","absent")) %>%
  filter(!grepl("genomic_bin",file.id)) %>% 
  select(qseqid,file.id,presence.absence)
# mutate(qseqid=fct_recode(qseqid,"mogA/moaB"="moaB","mogA/moaB"="mogA")) # recode gene names that don't match

genotype_mat<-panaroo %>% 
  mutate(geneid=fct_recode(geneid,"mogA/moaB"="moaB","mogA/moaB"="mogA")) %>% # recode gene names that don't match
  mutate(sample=gsub("Neisseria/|.fa","",gsub("(.+?_.+?)_.*" ,"\\1",sample))) %>%
  full_join(blast,by=c("geneid"="qseqid","sample"="file.id","presence.absence"))

mag_genotype_mat<-genotype_mat %>% 
  filter(!grepl("GC|GF",sample)) %>% 
  group_by(geneid) %>% 
  summarise(presence.absence=ifelse(any(presence.absence=="present"),"present","absent")) %>% 
  mutate(sample="This Study")

genotype_mat_reduced<-genotype_mat %>% 
  filter(grepl("GC|GF",sample)) %>% 
  bind_rows(mag_genotype_mat) %>% 
  pivot_wider(id_cols = sample,names_from = geneid,values_from = presence.absence,
              values_fn = function(x) ifelse(any(x=="present"),"present","absent")) %>%
  mutate_all(~replace_na(.,"absent")) %>% 
  column_to_rownames("sample") %>% 
  relocate(starts_with("nar"),starts_with("mo"))

n.p<-gheatmap(n.tree,
              genotype_mat_reduced,colnames=TRUE,high="black",low="gray",offset = 0.18,
              legend_title="Nitrate Genes",colnames_angle = 90,colnames_position = "top",hjust = 0,width = 0.85,
              family="calibri",colnames_offset_y = -0.4,font.size = 5)+
  # facet_grid()
  vexpand(0.115)+
  theme(legend.position = c(0.15,0.75),
        legend.background = element_blank())

ggsave(filename = "calculus/species-specific/DC-gorilla/mag/images/Neisseria-nitrate-tree.png",plot = n.p,dpi=300,width = 12,height = 12,units = "in")

#### ROTHIA
tree<-read.newick("calculus/species-specific/DC-gorilla/mag/phylophlan/Rothia/output/input.tre")
tree<-drop.tip(tree,tree$tip.label[grepl("sorted",tree$tip.label)][-1])

# clean up tips
tree$tip.label<-gsub("_mappedcontigs.sorted","",tree$tip.label)
#split string at second delimeter
tree$tip.label<-sub("(.*?_.*?)_.*", "\\1", tree$tip.label)
tree$node.label<-ifelse(tree$node.label=="1.000","",as.character(round(as.numeric(tree$node.label),3)))
# tree$node.label<-ifelse(as.numeric(tree$node.label)>=0.75,as.character(round(as.numeric(tree$node.label),digits = 2)),"")

# rescale the branch lengths
# tree$edge.length <- (log(tree$edge.length+0.01)+(min(log(tree$edge.length+0.01))*-1))

refs<-tree$tip.label[grepl("_",tree$tip.label)]
write.table(refs,file = "calculus/species-specific/DC-gorilla/mag/data/rothia-refs.txt",row.names = F,quote = F,col.names = F)

# system('refs=$(cat rothia-refs.txt);
#        for line in $refs; do res=$(esearch -query $line -db assembly | elink -target biosample | esummary | xtract -pattern DocumentSummary -element Taxonomy);
#        echo $line $res  >> rothia-refs-taxaids.txt;
#        done;')

taxa<-read.table("calculus/species-specific/DC-gorilla/mag/data/rothia-refs-taxaids.txt",header=F,colClasses = c("character","character")) %>% distinct(V1,.keep_all = T)
nrow(taxa) == length(refs)

tree_taxa<-foreach(i=unique(taxa$V2),.combine = rbind,.packages = "taxize") %do% {
  print(i)
  ncbi_get_taxon_summary(id = i)
}

taxa_merge<-tree_taxa %>% 
  right_join(taxa,by=c("uid"="V2")) %>% 
  bind_rows(.,bind_cols(uid="unknown",
                        name="This study",
                        rank="unknown",
                        V1=tree$tip.label[!grepl("_",tree$tip.label)])) %>%
  mutate(name=gsub("Rothia","R.",name),
         V3=ifelse(name=="This study","Gorilla","Human")) %>% 
  distinct(uid,V1,.keep_all = T) %>% 
  relocate(V1) %>% 
  mutate(V3=fct_relevel(V3,c("Gorilla","Human")))

# taxa_merge[taxa_merge$uid=="43675","name"]<-"R. mucilaginosa"
# taxa_merge[taxa_merge$uid=="37928","name"]<-"Arthrobacter crystallopoietes"
# taxa_merge[taxa_merge$uid=="37921","name"]<-"Arthrobacter agilis"
# taxa_merge[taxa_merge$uid=="2419771","name"]<-"Gryllotalpicola protaetiae"

# mutate(name=ifelse(!is.na(name),make.unique(name),name),
#   label=ifelse(grepl(paste(c("This study","mucilaginosa",  "kristinae",  "nasimurium",  "terrae",  "aeria",  "dentocariosa",  "amarae",  "koreensis"),collapse = "|"),name),
#                           NA,name))

# to.drop<-taxa_merge %>% filter(!grepl("Rothia|This study|Psychromicrobium|Reinbacterium|Arthrobacter|Pseudarthrobacter|Paenarthrobacter",name)) %>% pull(V1)
to.drop<-taxa_merge %>% 
  filter(V1!="G0021") %>% 
  separate(name,into = c("genus","species"),sep = " ",remove = F) %>% 
  group_by(species) %>% 
  slice(2:n()) %>% 
  pull(V1)

ggtree(drop.tip(tree,tip = to.drop)) %<+% taxa_merge +
  geom_tiplab(geom="text",align = T,size=6,family="calibri",offset = 0.1)+
  geom_text2(aes(label=node), hjust=-.5)+
  geom_treescale(width = 0.1)+
  hexpand(0.3)

r.tree<-ggtree(drop.tip(tree,c(to.drop,"GCF_004135285.1"))) %<+% taxa_merge + 
  # geom_hilight(node=204, fill=pal[1], alpha=.5) +
  geom_tiplab(geom="text",aes(label=name),align = T,offset = 3.7,size=6,family="calibri")+
  geom_tippoint(aes(color=V3),size=3,position=position_nudge(x = 0.01))+
  geom_tippoint(aes(subset=(isTip & node==10)),color=iso.pal[1],size=6,position=position_nudge(x = 0.01))+
  scale_color_manual( values = iso.pal[1:2])+
  guides(color=guide_legend("Host/Isolation Source"))+
  hexpand(ratio = 0.3)+
  geom_treescale(y=-1,width = 0.5,offset.label = 0.5,fontsize = 5,family = "calibri")+ # hide the old scale label behind the heatmap
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.position = c(0.25,0.75),
        text = element_text(family = "calibri"))

ggsave(plot = r.tree,filename = "calculus/species-specific/DC-gorilla/mag/images/rothia-tree.png",dpi=300,width = 18,height = 10,units = "in")

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
# table(panaroo$issues)

blast <- read_csv("calculus/species-specific/DC-gorilla/mag/data/all-Rothia.csv",skip=1,col_names=c("qseqid","sseqid","stitle","pident","qcovs","length","mismatch","gapopen","qstart","qend","sstart","send","qframe","sframe","frames","evalue","bitscore","qseq","sseq","file.id","type")) %>%
  mutate(file.id=gsub("Rothia/|.fa","",gsub("(.+?_.+?)_.*" ,"\\1",file.id))) %>%
  group_by(qseqid) %>%
  mutate(bin.pident=pident[file.id=="genomic_bin"],
         bin.length=length[file.id=="genomic_bin"],
         bin.sseq.length=str_length(str_replace_all(sseq, "[^[:alnum:]]", ""))) %>%
  group_by(file.id,qseqid) %>%
  mutate(presence.absence=ifelse((length/bin.length)>=0.8,"present","absent")) %>%
  filter(!grepl("genomic_bin",file.id)) %>% 
  select(qseqid,file.id,presence.absence)

genotype_mat<-panaroo %>% 
  mutate(geneid=fct_recode(geneid,"mogA/moaB"="moaB",
                           "mogA/moaB"="mogA")) %>% # recode gene names that don't match
  mutate(sample=gsub("Rothia/|.fa","",gsub("(.+?_.+?)_.*" ,"\\1",sample))) %>%
  full_join(blast,by=c("geneid"="qseqid","sample"="file.id","presence.absence"))

mag_genotype_mat<-genotype_mat %>% 
  filter(!grepl("GC|GF",sample)) %>% 
  group_by(geneid) %>% 
  summarise(presence.absence=ifelse(any(presence.absence=="present"),"present","absent")) %>% 
  mutate(sample="G0021")

genotype_mat_reduced<-genotype_mat %>% 
  filter(grepl("GC|GF",sample)) %>% 
  bind_rows(mag_genotype_mat) %>% 
  pivot_wider(id_cols = sample,names_from = geneid,values_from = presence.absence,
              values_fn = function(x) ifelse(any(x=="present"),"present","absent")) %>%
  mutate_all(~replace_na(.,"absent")) %>% 
  column_to_rownames("sample") %>% 
  relocate(starts_with("nar"),starts_with("mo"))

r.p<-gheatmap(r.tree,genotype_mat_reduced,colnames=TRUE,high="black",low="gray",offset = 0.1,
              legend_title="Nitrate Genes",colnames_angle = 90,colnames_position = "top",hjust = 0,width = 5,
              family="calibri",colnames_offset_y = -0.4,font.size = 5)+ 
  vexpand(0.5) +
  theme(legend.position = "none")

ggsave(filename = "calculus/species-specific/DC-gorilla/mag/images/Rothia-nitrate-tree-panaroo.png",plot = r.p,dpi=300,width = 8,height = 8,units = "in")

##### VEILLONELLA #####
tree<-read.newick("calculus/species-specific/DC-gorilla/mag/phylophlan/Veillonella/output/input_s__Veillonella_dispar/input.tre")
tree<-drop.tip(tree,tree$tip.label[grepl("sorted",tree$tip.label)][-1])
tree$tip.label[grepl("sorted",tree$tip.label)]<-"This Study"

# clean up tips
tree$tip.label<-gsub("^.*\\|", "", tree$tip.label)
#split string at second delimeter
tree$tip.label<-sub("(.*?_.*?)_.*", "\\1", tree$tip.label)

tree$node.label<-ifelse(tree$node.label=="1.000","",as.character(round(as.numeric(tree$node.label),3)))

# rescale branch lengths
tree$edge.length<-log(tree$edge.length+1)
# plot(rescale_tree(tree,"edge.length"))
# tree$edge.length <- (log(tree$edge.length+0.01)+(min(log(tree$edge.length+0.01))*-1))

plot(tree)

refs<-tree$tip.label[!grepl("sorted",tree$tip.label)]
refs<-tree$tip.label[grepl("_",tree$tip.label)]
write.table(refs,file = "calculus/species-specific/DC-gorilla/mag/data/veillonella-refs.txt",row.names = F,quote = F,col.names = F)

# system('refs=$(cat calculus/species-specific/DC-gorilla/mag/data/veillonella-refs.txt);
#        for line in $refs; do res=$(esearch -query $line -db assembly | elink -target biosample | esummary | xtract -pattern DocumentSummary -first Title -element Accession -group Attribute -if Attribute@harmonized_name -equals "isolation_source" -element Attribute);
#        echo $line $res  >> calculus/species-specific/DC-gorilla/mag/data/veillonella-refs-isolation-sources.raw.txt;
#        done;')

# edited by hand here...
taxa<-read.delim("calculus/species-specific/DC-gorilla/mag/data/veillonella-refs-taxaids.txt",header=F,sep = "\t",colClasses = c("factor","character"))
nrow(taxa) == length(refs)

tree_taxa<-foreach(i=unique(taxa$V2),.combine = rbind,.packages = "taxize") %do% {
  print(i)
  ncbi_get_taxon_summary(id = i)
}

iso<-read.csv("calculus/species-specific/DC-gorilla/mag/data/veillonella-refs-isolation-sources.csv",header=F)

taxa_merge<-tree_taxa %>% 
  right_join(taxa,by=c("uid"="V2")) %>% 
  right_join(iso,by="V1") %>% 
  rename("accession"="V1","isolation.source"="V2") %>% 
  bind_rows(.,bind_cols(uid="unknown",
                        name="This Study",
                        rank="unknown",
                        isolation.source="Gorilla",
                        short.iso="Gorilla",
                        accession="This Study")) %>% # tree$tip.label[grepl("sorted",tree$tip.label)])) 
  mutate(name=gsub("Veillonella","V.",name),
         lab=ifelse(!grepl("Study",name),paste0(name,", ",isolation.source),name)) %>% 
  distinct(uid,accession,.keep_all = T) %>% 
  relocate(accession) %>% 
  separate(isolation.source,into = "short.iso",sep = " ",remove = F) %>% 
  mutate(short.iso=ifelse(isolation.source %in% c("unknown",""),NA,ifelse(isolation.source=="rat","Mammalian",str_to_title(short.iso))))

to.drop<-taxa_merge %>% 
  filter(name != "This Study") %>% 
  separate(name,into = c("genus","species"),sep = " ",remove = F) %>% 
  group_by(species) %>% 
  slice(2:n()) %>% 
  pull(accession)

to.drop<-c(to.drop,"GCA_001553335","GCF_003006415.1","GCA_002959855","GCA_000215025")

ggtree(drop.tip(tree,to.drop)) %<+% taxa_merge + 
  geom_tiplab(geom="text",aes(label=name)) +
  geom_text2(aes(subset=isTip, label=node), hjust=-.5)

v.tree<-ggtree(drop.tip(tree,to.drop)) %<+% taxa_merge + 
  geom_tiplab(geom="text",aes(label=name),align = T,offset = 2.66,size=6,family="calibri")+
  geom_tippoint(aes(color=short.iso),size=3,position=position_nudge(x = 0.01))+
  geom_tippoint(aes(subset=(isTip & node==5)),color=iso.pal[1],size=6,position=position_nudge(x = 0.01))+
  scale_color_manual( values = iso.pal[c(1,2,4)])+
  guides(color=guide_legend("Host/Isolation Source"))+
  hexpand(ratio = 0.25)+
  geom_treescale(y=-1,fontsize = 5,family = "calibri",width = 0.5)+
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.position = c(0.25,0.75),
        text = element_text(family = "calibri"))

ggsave(filename = "calculus/species-specific/DC-gorilla/mag/images/veillonella-tree.png",plot = v.tree,dpi=300,width = 14,height = 14,units = "in")

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
  ungroup() %>% 
  filter(!grepl("genomic_bin",file.id)) %>%
  dplyr::select(qseqid,file.id,presence.absence) %>% 
  mutate(qseqid=fct_recode(qseqid,`nirB/nirD`="nirD"))# recode gene names that are identical

genotype_mat<-panaroo %>% 
  mutate(geneid=fct_recode(geneid,"mogA/moaB"="moaB")) %>% # recode gene names that don't match
  mutate(sample=gsub("Veillonella/|_mappedcontigs.sorted.fa|.fna","",gsub("(.+?_.+?)_.*" ,"\\1",sample))) %>%
  full_join(blast,by=c("geneid"="qseqid","sample"="file.id","presence.absence"))

mag_genotype_mat<-genotype_mat %>% 
  filter(!grepl("GC|GF",sample)) %>% 
  group_by(geneid) %>% 
  summarise(presence.absence=ifelse(any(presence.absence=="present"),"present","absent")) %>% 
  mutate(sample="This Study")

genotype_mat_reduced<-genotype_mat %>% 
  filter(grepl("GC|GF",sample)) %>% 
  bind_rows(mag_genotype_mat) %>% 
  pivot_wider(id_cols = sample,names_from = geneid,values_from = presence.absence,
              values_fn = function(x) ifelse(any(x=="present"),"present","absent")) %>%
  mutate_all(~replace_na(.,"absent")) %>% 
  column_to_rownames("sample") %>% 
  relocate(starts_with("nar"),starts_with("mo"))

v.p<-gheatmap(v.tree,
              genotype_mat_reduced,width = 1.25,
              colnames=TRUE,legend_title="Nitrate Genes",family = "calibri",offset = 0.1,
              colnames_angle = 90,colnames_position = "top",colnames_offset_y = -0.25,hjust = 0,font.size = 5)+
  vexpand(0.35)+
  theme(legend.position = "bottom",
        legend.background = element_blank())

ggsave(filename = "calculus/species-specific/DC-gorilla/mag/images/Veillonella-nitrate-tree-pident.png",plot = v.p,dpi=300,width = 12,height = 12,units = "in")

##### L gorillae
tree<-read.tree("calculus/species-specific/DC-gorilla/mag/phylophlan/Lactobacillus/output/input.tre")
tree<-drop.tip(tree,tree$tip.label[grepl("sorted",tree$tip.label)][-1])
# clean up tips
tree$tip.label<-gsub("^.*\\|", "", tree$tip.label)

# tree$node.label<-ifelse(as.integer(tree$node.label)<0.5,"< 0.5","")
tree$node.label[c(2:10,13:21)]<-""
tree$node.label[1]<-"1.000"

refs<-tree$tip.label[!grepl("sorted",tree$tip.label)]
write.table(refs,file = "calculus/species-specific/DC-gorilla/mag/data/core-lacto-refs.txt",row.names = F,quote = F,col.names = F)

# system('refs=$(cat core-lacto-refs.txt);
#        for line in $refs; do res=$(esearch -query $line -db assembly | elink -target biosample | esummary | xtract -pattern DocumentSummary -element Taxonomy);
#        echo $line $res  >> core-lacto-refs-taxaids.txt;
#        done;')

taxa<-read.delim("calculus/species-specific/DC-gorilla/mag/data/core-lacto-refs-taxaids.txt",header=F,na.strings = "NA",colClasses = c("character","character","character","character"))
nrow(taxa) == length(refs)

tree_taxa<-foreach(i=unique(taxa$V2),.combine = rbind,.packages = "taxize") %do% {
  ncbi_get_taxon_summary(id = i)
}

taxa_merge<-tree_taxa %>%
  right_join(taxa,by=c("uid"="V2")) %>% 
  bind_rows(.,bind_cols(uid="1450649",rank="unknown",V1=tree$tip.label[grepl("sorted",tree$tip.label)],V3="Gorilla",name="This study",V4="oral",V5="museum specimen")) %>% 
  distinct(uid,V1,.keep_all = T) %>% 
  mutate(lab=ifelse(name=="This study",name,paste0(name,", ",V4))) %>% 
  relocate(V1) %>% 
  mutate(lab=gsub("Limosilactobacillus","L.",lab))

ggtree(phangorn::midpoint(tree)) %<+% taxa_merge + 
  # geom_tiplab(geom="text",aes(label=name)) +
  geom_text2(aes(subset=!isTip, label=node), hjust=-.5)

# x<-phangorn::midpoint(tree)
# m <- MRCA(x, 5)
x <- groupClade(phangorn::midpoint(tree), 5)
y <- groupClade(phangorn::midpoint(tree), 6)
z <- groupClade(phangorn::midpoint(tree), 7)
lb.tree<-ggtree(groupClade(phangorn::midpoint(tree), 7),aes(linetype=group),ladderize = T) %<+% taxa_merge +
  # geom_nodelab(nudge_x = -0.0075,nudge_y = 0.4,size=3)+
  scale_linetype_manual(values = c(1,1))+
  # geom_text2(aes(subset=!isTip, label=node), hjust=-.5)+
  geom_tiplab(geom="text",aes(label=lab),align = T,offset = 0.05,size=6,family="calibri")+
  guides(linetype = "none")+
  geom_tippoint(aes(color=V3),size=3,position = position_nudge(x=0.0025))+
  geom_tippoint(aes(subset=(isTip & node==1)),color=iso.pal[1],size=6,position=position_nudge(x = 0.0025))+
  scale_color_manual( values = c(iso.pal[1],get_palette("Paired",7)[7],iso.pal[2]))+
  guides(color=guide_legend("Host/Isolation Source"))+
  geom_treescale(x=1,offset = 0.1,fontsize = 5)+
  # coord_cartesian(xlim = c(0.325,0.45))+
  vexpand(ratio = 0.1)+
  hexpand(ratio = 0.2)+
  theme(text = element_text(family = "calibri"))

lb.tree$data[lb.tree$data$node == 5,"x"] <- mean(lb.tree$data$x)-0.1
lb.tree$data[lb.tree$data$node == 6,"x"] <- mean(lb.tree$data$x)-0.1
lb.tree$data[lb.tree$data$node == 7,"x"] <- mean(lb.tree$data$x)-0.1

lb.tree<-lb.tree + theme(legend.text = element_text(size = 14),legend.title = element_text(size = 16),legend.position = "bottom")

ggsave(filename = "calculus/species-specific/DC-gorilla/mag/images/lacto-core-tree.png",plot = lb.tree,dpi=300,width = 8,height = 6,units = "in")

# abundances --------------------------------------------------------------
meta<-read.csv("calculus/all-samples/decontamination/rank-abundance-filtering/reproducibility-tests/Gorilla-Markella/metadata_gorilla_analysis.csv",header=T)

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
  group_by(Genomic.bins,family,genus,species) %>% 
  summarise(prev=n_distinct(Seq.label[abundance!=0])/46) %>% 
  arrange(-prev)

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
  scale_color_manual(values = c("gray",pal))+
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

# L. gorillae
`Lactobacillus_H gorillae.p`<-dat10p %>% 
  filter(species=="Lactobacillus_H gorillae") %>% 
  mutate(Spec.subspecies=fct_relevel(Spec.subspecies,"Western","Grauer's","Mountain"),
         species="Limosilactobacillus gorillae") %>% 
  ggplot(aes(x=Spec.subspecies,y=rel.abund,fill=Spec.subspecies))+
  geom_boxplot(outlier.shape = NA)+
  labs(y="Relative Abundance")+
  # scale_x_discrete(limits=c("Western","Grauer's","Mountain"))+
  guides(fill=guide_legend("Subspecies"))+
  # scale_y_continuous(expand = expansion(mult = 0.15))+
  # facet_wrap(~species,scales = "free_y")+
  coord_cartesian(ylim=c(0,0.05))+
  scale_fill_manual(values = pal)+
  geom_signif(comparisons=list(c("Mountain","Grauer's"),c("Mountain","Western"),c("Grauer's","Western")),textsize = 6,
              y_position = c(0.01,0.011,0.012),map_signif_level = T,step_increase = 0.01,tip_length = 0.0035)+
  theme_bw()+
  theme(axis.text = element_text(size=16),
        axis.title.y = element_text(size = 18),
        strip.text = element_text(size=18),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        text = element_text(family = "calibri"))

# Neisseria
Neisseria.p<-dat10p %>% 
  filter(species=="Neisseria") %>% 
  mutate(Spec.subspecies=fct_relevel(Spec.subspecies,"Western","Grauer's","Mountain")) %>% 
  ggplot(aes(x=Spec.subspecies,y=rel.abund,fill=Spec.subspecies))+
  geom_boxplot(outlier.shape = NA)+
  labs(y="Relative Abundance",x="Gorilla Subspecies")+
  # scale_x_discrete(limits=c("Western","Grauer's","Mountain"))+
  # facet_wrap(~species,scales = "free_y")+
  guides(fill=guide_legend("Subspecies"))+
  # scale_y_continuous(expand = expansion(mult = 0.15))+
  # coord_cartesian(ylim=c(0,0.15))+
  scale_fill_manual(values = pal)+
  geom_signif(comparisons=list(c("Mountain","Grauer's"),c("Mountain","Western"),c("Grauer's","Western")),textsize = 8, 
              y_position = c(0.12,0.125,0.13),map_signif_level = T,step_increase = 0.1)+
  theme_bw()+
  theme(axis.text = element_text(size=20),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        text = element_text(family = "calibri"))
n.tree.legend<-get_legend(n.p)

# Veillonella
Veillonella.p<-dat10p %>% 
  filter(species=="Veillonella") %>% 
  mutate(Spec.subspecies=fct_relevel(Spec.subspecies,"Western","Grauer's","Mountain")) %>% 
  ggplot(aes(x=Spec.subspecies,y=rel.abund,fill=Spec.subspecies))+
  geom_boxplot(outlier.shape = NA)+
  labs(y="Relative Abundance")+
  # scale_x_discrete(limits=c("Western","Grauer's","Mountain"))+
  guides(fill=guide_legend("Subspecies"))+
  # scale_y_continuous(expand = expansion(mult = 0.15))+
  facet_wrap(~species,scales = "free_y")+
  coord_cartesian(ylim=c(0,0.3))+
  scale_fill_manual(values = pal)+
  geom_signif(comparisons=list(c("Mountain","Grauer's"),c("Mountain","Western"),c("Grauer's","Western")),textsize = 6,
              y_position = c(0.12,0.125,0.13)*1.5,map_signif_level = T,step_increase = 0.1)+
  theme_bw()+
  theme(axis.text = element_text(size=18),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size=18),
        legend.position = "none",
        text = element_text(family = "calibri"))
legend<-get_legend(`Lactobacillus_H gorillae.p`)
ggsave(Veillonella.p,file="calculus/species-specific/DC-gorilla/mag/images/Veillonella-abundance-sig.png",dpi=300,height=8,width = 6,units = "in")

tree.legend<-get_legend(v.tree)

##### Rothia
Rothia.p<-dat10p %>% 
  filter(species=="Rothia") %>% 
  mutate(Spec.subspecies=fct_relevel(Spec.subspecies,"Western","Grauer's","Mountain")) %>% 
  ggplot(aes(x=Spec.subspecies,y=rel.abund,fill=Spec.subspecies))+
  geom_boxplot(outlier.shape = NA)+
  labs(y="Relative Abundance")+
  # scale_x_discrete(limits=c("Western","Grauer's","Mountain"))+
  facet_wrap(~species,scales = "free_y")+
  guides(fill=guide_legend("Subspecies"))+
  # scale_y_continuous(expand = expansion(mult = 0.15))+
  coord_cartesian(ylim=c(0,0.1))+
  scale_fill_manual(values = pal)+
  geom_signif(comparisons=list(c("Mountain","Grauer's"),c("Mountain","Western"),c("Grauer's","Western")),textsize = 6,
              y_position = c(0.12,0.125,0.13)*0.5,map_signif_level = T,step_increase = 0.1)+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text = element_text(size=18),
        strip.text = element_text(size=18),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        text = element_text(family = "calibri"))


# now combine into multipanel figures -------------------------------------
# Neisseria figure
ngp<-ggarrange(n.p,Neisseria.p,
               labels = c("a)","b)"),nrow = 1,
               widths = c(2,1),
               font.label = list(size=28),
               label.x = c(0,-0.1)) + theme(text = element_text(family="calibri"))
ggsave(ngp,file="calculus/species-specific/DC-gorilla/mag/images/neisseria-comb-figure.png",dpi=300,height=8,width = 20.5,units = "in")

lgp<-ggarrange(lb.tree,
               `Lactobacillus_H gorillae.p`,labels = c("a)","b)"),
               widths = c(1.5,1),
               font.label = list(size=28),label.x = c(0,-0.1),
               ncol=2)
ggsave(lgp,file="calculus/species-specific/DC-gorilla/mag/images/lacto-comb-figure.png",dpi=300,height=6,width = 12,units = "in")

rvgp<-ggarrange(
  annotate_figure(ggarrange(r.p,Rothia.p,
                            v.p+theme(legend.position = "none"),Veillonella.p,
                            labels=c("a)","b)","c)","d)"),
                            font.label = list(size=28),
                            ncol=2,
                            nrow=2,
                            widths = c(1.25,0.65),
                            # widths = c(1.5,1,1.5,0.5,0.5,0.5),
                            # heights = c(1,0.5,1,1,1,1),
                            align = "hv",
                            label.x = c(0,-0.1)),
                  left = grid::textGrob("Relative Abundance", rot = 90, vjust = 38.5, gp = grid::gpar(cex = 2))) +
    theme(text = element_text(family="calibri")),
  get_legend(v.p),
  heights = c(1,0.05),ncol = 1)
ggsave(rvgp,file="calculus/species-specific/DC-gorilla/mag/images/rothia-veillonella-comb-figure.png",dpi=300,height=10,width = 14,units = "in")
