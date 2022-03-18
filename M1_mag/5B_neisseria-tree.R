# neisseria. gorillae tree
# based on pheS sequences
library(foreach)
library(ggtree)
library(treeio)
library(tidyverse)

tree<-read.tree("calculus/DC2/DC-gorilla/mag/phylophlan/Neisseria/input_s__Neisseria_gonorrhoeae/input.tre")
plot(tree)

# clean up tips
tree$tip.label<-gsub("^.*\\|", "", tree$tip.label)
tree$tip.label<-sub("(.*?_.*?)_.*", "\\1", tree$tip.label)

# tree$node.label<-ifelse(tree$node.label>0.8 | tree$node.label==0.169 ,tree$node.label,"")
tree$node.label[c(1:18,49:50)]<-""

refs<-tree$tip.label[!grepl("sorted",tree$tip.label)]
write.table(refs,file = "calculus/DC2/DC-gorilla/mag/data/neisseria-refs.txt",row.names = F,quote = F,col.names = F)

# system('refs=$(cat calculus/DC2/DC-gorilla/mag/data/neisseria-refs.txt);
#        for line in $refs; do res=$(esearch -query $line -db assembly | elink -target biosample | esummary | xtract -pattern DocumentSummary -first Title -element Accession -group Attribute -if Attribute@harmonized_name -equals "isolation_source" -element Attribute);
#        echo $line $res  >> calculus/DC2/DC-gorilla/mag/data/neisseria-refs-isolation-sources.raw.txt;
#        done;')


# system('refs=$(cat neisseria-refs.txt);
#        for line in $refs; do res=$(esearch -query $line -db assembly | elink -target biosample | esummary | xtract -pattern DocumentSummary -element Taxonomy);
#        echo $line $res  >> neisseria-refs-taxaids.txt;
#        done;')

iso<-read.delim("calculus/DC2/DC-gorilla/mag/data/neisseria-refs-isolation-sources.csv",header=F,sep = ",")

taxa<-read.table("calculus/DC2/DC-gorilla/mag/data/neisseria-refs-taxaids.txt",colClasses = c("character","character"),header=F)
nrow(taxa) == length(refs)

tree_taxa<-foreach(i=unique(taxa$V2),.combine = rbind,.packages = "taxize") %do% {
  ncbi_get_taxon_summary(id = i)
}

taxa_merge<-tree_taxa %>% 
  right_join(taxa,by=c("uid"="V2")) %>% 
  right_join(iso,by="V1") %>% 
  rename("accession"="V1","isolation.source"="V2") %>% 
  bind_rows(.,bind_cols(uid="unknown",
                        name="This study",
                        rank="unknown",
                        V3="Gorilla",V4="wild",
                        accession=tree$tip.label[grepl("sorted",tree$tip.label)])) %>% 
  mutate(name=gsub("Neisseria","N.",name),
         lab=ifelse(name=="This study",NA,ifelse(grepl("Eikenella|Snodgrassella",name),name,paste0(name,", ",V4)))) %>% 
  distinct(uid,accession,.keep_all = T) %>% 
  relocate(accession)

ggtree(phangorn::midpoint(drop.tip(tree,tip = c("GCF_900187105.1","GCF_014054965.1")))) %<+% taxa_merge +
  geom_text2(aes(subset=!isTip, label=node), hjust=-.5)

n.tree<-ggtree(phangorn::midpoint(drop.tip(tree,tip = c("GCF_900187105.1","GCF_014054965.1")))) %<+% taxa_merge +
  geom_nodelab(nudge_x = -0.02,nudge_y = 0.5,size=3)+
  geom_cladelabel(node=72, label="This Study",align=T,fontsize = 5,offset = 0.02,offset.text = 0.02) +
  # geom_text2(aes(subset=!isTip, label=node), hjust=-.5)+
  geom_tiplab(geom="text",aes(label=lab),offset = 0.02,align = T,size=4)+
  geom_tippoint(aes(fill=V3),shape=21,size=3,position=position_nudge(x = 0.01))+
  scale_fill_manual( values = c(ggpubr::get_palette(palette = "Set1", 5),"gray"))+
  guides(fill=guide_legend("Host/Isolation Source"))+
  # vexpand(ratio = 0.1)+
  hexpand(ratio = 0.25)+
  # vexpand(ratio = 0.25)+
  geom_treescale(y = -1)+
  theme(legend.text = element_text(size = 14),legend.title = element_text(size = 16),legend.position = c(0.25,0.75))

ggsave(filename = "calculus/DC2/DC-gorilla/mag/images/neisseria-tree.png",plot = n.tree,dpi=300,width = 12,height = 12,units = "in")

# p$data[p$data$node ==7, "x"] <- mean(p$data$x)
