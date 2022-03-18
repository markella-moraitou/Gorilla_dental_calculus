# lacto. gorillae tree
library(tidyverse)
library(foreach)
library(ggtree)
library(treeio)

######## core genes ############
tree<-read.tree("calculus/DC2/DC-gorilla/mag/phylophlan/Lactobacillus/output/input.tre")

# clean up tips
tree$tip.label<-gsub("^.*\\|", "", tree$tip.label)

# tree$node.label<-ifelse(as.integer(tree$node.label)<0.5,"< 0.5","")

tree$node.label[c(2:10,13:21)]<-""
tree$node.label[1]<-"1.000"

refs<-tree$tip.label[!grepl("sorted",tree$tip.label)]
write.table(refs,file = "calculus/DC2/DC-gorilla/mag/data/core-lacto-refs.txt",row.names = F,quote = F,col.names = F)

# system('refs=$(cat core-lacto-refs.txt);
#        for line in $refs; do res=$(esearch -query $line -db assembly | elink -target biosample | esummary | xtract -pattern DocumentSummary -element Taxonomy);
#        echo $line $res  >> core-lacto-refs-taxaids.txt;
#        done;')

taxa<-read.delim("calculus/DC2/DC-gorilla/mag/data/core-lacto-refs-taxaids.txt",header=F,na.strings = "NA",colClasses = c("character","character","character","character"))
nrow(taxa) == length(refs)

tree_taxa<-foreach(i=unique(taxa$V2),.combine = rbind,.packages = "taxize") %do% {
  ncbi_get_taxon_summary(id = i)
}

taxa_merge<-tree_taxa %>%
  right_join(taxa,by=c("uid"="V2")) %>% 
  bind_rows(.,bind_cols(uid="1450649",rank="unknown",V1=tree$tip.label[grepl("sorted",tree$tip.label)],V3="Gorilla",name="This study",V4="oral",V5="museum specimen")) %>% 
  distinct(uid,V1,.keep_all = T) %>% 
  mutate(lab=ifelse(name=="This study",NA,paste0(name,", ",V4))) %>% 
  relocate(V1) %>% mutate(lab=gsub("Limosilactobacillus","L.",lab))

ggtree(phangorn::midpoint(tree))+
  geom_text2(aes(subset=!isTip, label=node), hjust=-.5)

# x<-phangorn::midpoint(tree)
# m <- MRCA(x, 5)
y <- groupClade(phangorn::midpoint(tree), 51)
z <- groupClade(phangorn::midpoint(tree), 30)
lb.tree<-ggtree(groupClade(phangorn::midpoint(tree), 30),aes(linetype=group),ladderize = T) %<+% taxa_merge +#%<+% d2 + aes(linetype=I(lty))
  geom_nodelab(nudge_x = -0.0075,nudge_y = 0.4,size=3)+
  scale_linetype_manual(values = c(1,2))+
  # geom_text2(aes(subset=!isTip, label=node), hjust=-.5)+
  geom_tiplab(geom="text",aes(label=lab),offset = 0.05,align = T,size=5)+
  guides(linetype = "none")+
  geom_tippoint(aes(fill=V3),size=3,shape=21,position = position_nudge(x=0.0025))+
  scale_fill_manual( values = ggpubr::get_palette(palette = "Set1",6)[c(2,6,3)])+
  guides(fill=guide_legend("Host/Isolation Source"))+
  geom_treescale(x=1.25,offset = 0.2,y=-0.75)+
  # coord_cartesian(xlim = c(0.325,0.45))+
  vexpand(ratio = 0.1)+
  hexpand(ratio = 0.001)

lb.tree$data[lb.tree$data$node == 51,"x"] <- mean(lb.tree$data$x)-0.1
lb.tree$data[lb.tree$data$node == 30,"x"] <- mean(lb.tree$data$x)-0.1

lb.tree<-lb.tree + geom_cladelabel(node=49,extend = 12, label="This Study",align=F,offset = 0.035,offset.text = 0,fontsize = 5)+theme(legend.text = element_text(size = 14),legend.title = element_text(size = 16),legend.position = c(0.25,0.75))

ggsave(filename = "calculus/DC2/DC-gorilla/mag/images/lacto-core-tree.png",plot = lb.tree,dpi=300,width = 8,height = 6,units = "in")

######## pheS ##########
# based on pheS sequences
tree<-read.tree("~/Downloads/pheS_tree/pheS.tre")
# clean up tips
tree$tip.label<-gsub("^.*\\|", "", tree$tip.label)
tree$node.label<-ifelse(as.numeric(tree$node.label)>=0.75,as.character(round(as.numeric(tree$node.label),digits = 2)),"")

refs<-tree$tip.label[!grepl("length",tree$tip.label)]
write.table(refs,file = "calculus/DC2/DC-gorilla/mag/data/lacto-refs.txt",row.names = F,quote = F,col.names = F)

# system('refs=$(cat lacto-refs.txt);
#        for line in $refs; do res=$(esearch -query $line -db assembly | elink -target biosample | esummary | xtract -pattern DocumentSummary -element Taxonomy);
#        echo $line $res  >> lacto-refs-taxaids.txt;
#        done;')

taxa<-read.table("calculus/DC2/DC-gorilla/mag/data/lacto-refs-taxaids.txt",colClasses = c("character","character","character"),header=F)
nrow(taxa) == length(refs)

tree_taxa<-foreach(i=unique(taxa$V2),.combine = rbind,.packages = "taxize") %do% {
  ncbi_get_taxon_summary(id = i)
}

taxa_merge<-tree_taxa %>% right_join(taxa,by=c("uid"="V2")) %>% 
  bind_rows(.,bind_cols(uid="unknown",name="L. gorillae MAG\nConstructed in this study",rank="unknown",V1=tree$tip.label[grepl("length",tree$tip.label)],V3="museum specimen")) %>% 
  distinct(uid,V1,.keep_all = T) %>% 
  relocate(V1) %>% mutate(name=gsub("Limosilactobacillus","L.",name))

p<-ggtree(phangorn::midpoint(tree))+
  geom_text2(aes(subset=!isTip, label=node), hjust=-.5)

p<-ggtree(groupClade(phangorn::midpoint(tree), 11),aes(linetype=group),ladderize = T) %<+% taxa_merge +#%<+% d2 + aes(linetype=I(lty))
  geom_nodelab(nudge_x = -0.01,nudge_y = 0.15)+
  scale_linetype_manual(values = c(2,1))+
  # geom_text2(aes(subset=!isTip, label=node), hjust=-.5)+
  geom_tiplab(geom="text",aes(label=name),offset = 0.075,align = T)+
  geom_tippoint(aes(color=V3),size=3)+
  guides(linetype="none")+
  scale_color_manual( values = ggpubr::get_palette(palette = "Set1",4))+
  guides(color=guide_legend("Host/Isolation Source"))+
  geom_treescale(x=0.35,offset = 0.1)+
  vexpand(ratio = 0.1)+
  hexpand(ratio = 0.3)+
  theme(legend.text = element_text(size = 14),legend.title = element_text(size = 16))
p$data[p$data$node == 11,"x"] <- mean(p$data$x)+0.01

ggsave(filename = "calculus/DC2/DC-gorilla/mag/images/lacto-pheS-tree.png",plot = p,dpi=300,width = 10,height = 6,units = "in")
