library(tidyverse)
library(treedataverse)

#collapse_nodes() v0.1
#---------------------------------------------------------------------------------------
#Taking a phylogenetic tree, a data.frame of taxonomy strings and a taxonomic level
#Select tips in each taxonomic group, then calculate most recent common ancestor 
#(MRCA) node and # of tips in each node. Return data.frame.

collapse_nodes <- function(level, tree, taxonomy) {
  require(ape)  #requires tree structure from ape
  require(phytools)  #requires the findMRCA function from phytools
  nodelist<-list() #create empty vectors
  nodedist<-list()
  grouplist<-list()
  groupcount<-list()
  i <- 1
  #Select unique taxonomic levels
  #Ensure data.frame entries only include tip names
  taxonomy <- taxonomy[which(taxonomy[,1] %in% tree$tip.label), ]
  groups <- as.character(unique(taxonomy[,level]))
  
  for(group in groups) {
    print(group)
    #For each group, find all tips
    group_tips <- as.character(taxonomy[which(taxonomy[,level] == group),]$V1)
    #Stuff will break if you call MRCA on groups with only 1 tip.  
    if(length(group_tips) > 1) {
      node <- findMRCA(tree, group_tips, type="node")
      provtree <- keep.tip(tree, group_tips)
      justthetip <- provtree$edge.length
      groupcount[i] <- length(group_tips)
      grouplist[i] <- group
      nodelist[i] <- node
      nodedist[i] <- max(justthetip)
      i <- i + 1
    }
  }
  #Unlisting a list of list is a hacky way to mimic the python "append" function in R. 
  return(data.frame(Group = unlist(grouplist), TipCount = unlist(groupcount), Node = unlist(nodelist), Distance = unlist(nodedist)))
}

# read in full tree
fulltree<-phytools::read.newick("calculus/DC2/DC-gorilla/mag/data/gtdbtk.bac120.classify.clean.tree")
tips <- fulltree$tip.label

# read in full taxonomy
gtdb.taxa<-read.delim("calculus/DC2/DC-gorilla/mag/data/bac120_taxonomy_r202.tsv",header = F) %>% separate(col = V2,sep = ";",into = c("Domain","Phylum","Class","Order","Family","Genus","Species"))

# read in bin taxa
bin.taxa<-read.delim("calculus/DC2/DC-gorilla/mag/data/gtdbtk-refs-taxaids.txt",
                 header=F,colClasses = c("character","character")) %>% 
  rename("classification"="V2") %>% 
  mutate(type=ifelse(grepl("bin",V1),"MAG Constructed\nin this study","GTDB Reference Genome"))

bin.taxa %>% filter(!grepl("Desulfovibrionaceae|s__Exiguobacterium|Propionibacterium|f__Planococcaceae|Erysipelothrix|sp002140795|s__Pseudogracilibacillus",classification))

contams <- c("bin.33","bin.28","bin.12","bin.32","bin.8","bin.11","bin.23","bin.1","bin.3","bin.14")

classif <- read.delim("calculus/DC2/DC-gorilla/mag/data/gtdbtk.bac120.summary.tsv", na.strings = '') %>% separate(col = classification,sep = ";",into = c("Domain","Phylum","Class","Order","Family","Genus","Species"))

############
subtree<-phytools::read.newick("calculus/DC2/DC-gorilla/mag/data/gtdbtk.bac120.classify.clean.sub.tre")

bins<-subtree$tip.label[grepl("bin",subtree$tip.label)]

treedata<-subtree %>% as.treedata %>% as_tibble()

refs<-subtree$tip.label[grepl("RS_",subtree$tip.label)]

ref.tree<-ggtree(subtree,layout = "fan")  %<+% bin.taxa +
  geom_label(aes(label=node))+
  geom_tiplab(geom="text",aes(label=classification),align = T,size=3)
ggsave(filename = "calculus/DC2/DC-gorilla/mag/images/gtdbtk-subtree-ref.png",ref.tree,width = 16,height = 16,units = "in")

tips <- subtree$tip.label
filtered.taxonomy<-gtdb.taxa[gtdb.taxa$V1 %in% tips,]

Phyla.nodes <- collapse_nodes("Phylum", subtree, filtered.taxonomy)
Class.nodes <- collapse_nodes("Class", subtree, filtered.taxonomy)
Order.nodes <- collapse_nodes("Order", subtree, filtered.taxonomy)
Family.nodes <- collapse_nodes("Family", subtree, filtered.taxonomy)
Genus.nodes <- collapse_nodes("Genus", subtree, filtered.taxonomy)

# acido_node <- Phyla.nodes[which(Phyla.nodes$Group == "p__Acidobacteriota"),]$Node
proteo_node <- Phyla.nodes[which(Phyla.nodes$Group == "p__Proteobacteria"),]$Node
actino_node <- 103
firm_node <- Phyla.nodes[which(Phyla.nodes$Group == "p__Firmicutes"),]$Node
bact_node <- Phyla.nodes[which(Phyla.nodes$Group == "p__Bacteroidota"),]$Node

desulf_node <- Class.nodes[which(Class.nodes$Group == "c__Desulfovibrionia"),]$Node
corio_node <- 113
bacill_node <- Class.nodes[which(Class.nodes$Group=="c__Bacilli"),]$Node
gamma_node <- Class.nodes[which(Class.nodes$Group == "c__Gammaproteobacteria"),]$Node
clost_node <- 75
negat_node <- 76
sacc_node <- 118
plant_node <- 119

pal<-ggpubr::get_palette(palette = "Paired",k = 10)

sub.tree<-ggtree(drop.tip(subtree,tip = contams),layout = "fan") %<+% bin.taxa +
  # geom_text2(aes(subset=!isTip, label=node), hjust=-.5)+
  # scale_colour_manual(values = c("black","red"))+
  geom_cladelabel(fontsize=6,node=corio_node,barsize=1,label="Coriobacteriia",offset=0.1,offset.text = 1,align=TRUE,angle = -47.5) +
  geom_cladelabel(fontsize=6,node=sacc_node,barsize=1,label="Saccharimonadia",offset=0.1,offset.text=1.25,align=TRUE,angle=-67.5)+
  geom_cladelabel(fontsize=6,node=plant_node,barsize=1,label="Planctomycetes",offset=0.1,offset.text=1.2,align=TRUE,angle = -77)+
  geom_cladelabel(fontsize=6,node=bact_node,barsize=0,label="Bacteroidota",offset=0.1,offset.text=0,align=TRUE,angle = 80) +
  geom_cladelabel(fontsize=6,node=gamma_node,barsize=1,label="Gammaproteobacteria",offset=0.1,offset.text=0.1,align=TRUE,angle = 45)+
  geom_cladelabel(fontsize=6,node=desulf_node,barsize=1,label="Desulfovibrionia",offset=0.1,offset.text=0.1,align=TRUE,angle = 10) +
  geom_cladelabel(fontsize=6,node=bacill_node,barsize=0,label="Bacilli",offset=0.1,offset.text=-0.5,align=TRUE,angle = -70) +
  geom_cladelabel(fontsize=6,node=negat_node,barsize=1,label="Negativicutes",offset=0.1,offset.text=0.95,align=TRUE,angle = 40) +
  geom_cladelabel(fontsize=6,node=clost_node,barsize=1,label="Clostridia",offset=0.1,offset.text=0.66,align=TRUE,angle = 32.5) +
  geom_cladelabel(fontsize=6,node=actino_node,barsize=0,label="Actinobacteriota",offset=0.1,offset.text=1,align=TRUE,angle = -5) +
  geom_hilight(node=actino_node,fill=pal[2],alpha=0.5) +
  geom_hilight(node=corio_node,fill=pal[3],alpha=0.5) +
  geom_hilight(node=bact_node,fill=pal[1],alpha=0.5) +
  geom_hilight(node=desulf_node,fill=pal[4],alpha=0.5) +
  geom_hilight(node=clost_node,fill=pal[5],alpha=0.5) +
  geom_hilight(node=bacill_node,fill=pal[6],alpha=0.5) +
  geom_hilight(node=sacc_node,fill=pal[7],alpha=0.5)+
  geom_hilight(node=plant_node,fill=pal[8],alpha=0.5)+
  geom_hilight(node=gamma_node,fill=pal[9],alpha=0.5)+
  geom_hilight(node=negat_node,fill=pal[10],alpha=0.5)+
  geom_tippoint(aes(fill=type),shape=21,color="black",size=3)+
  scale_fill_manual(values = c("#7f7f7f","black"))+
  geom_treescale(x = 2)+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),legend.position = c(0.5,0.2))
ggsave(filename = "calculus/DC2/DC-gorilla/mag/images/gtdbtk-subtree.png",sub.tree,width = 12,height = 12,dpi = 600,units = "in")
