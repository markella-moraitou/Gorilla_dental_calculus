#Start logging
sink(file = "log1_mag_stats.txt")

#load packages
library(stringr)
library(microbiome)
library(phyloseq)
library(ape)
library(ggtree)
library(phangorn)
library(phylotools)
library(RColorBrewer)
library(cowplot)
library(dplyr)

#MAG taxonomy and trees

load("/proj/sllstore2017021/nobackup/MARKELLA/T3_community-level/.RData")

#### Load GTDB-Tk summaries ####

hq_mq_mags <- read.table("/proj/sllstore2017021/nobackup/MARKELLA/M4_MAGtaxonomy/hq_mq.bac120.summary.tsv", sep="\t", comment.char="", header=T)

#Add quality info
hq_mq_mags$draft_quality <- ifelse(hq_mq_mags$user_genome %in% 
                                      read.table(file="/crex/proj/sllstore2017021/nobackup/MARKELLA/M3_binning/high_quality_MAGs.txt")$V2, "HQ",
                                        ifelse(hq_mq_mags$user_genome %in% 
                                      read.table(file="/crex/proj/sllstore2017021/nobackup/MARKELLA/M3_binning/medium_quality_MAGs.txt")$V2, "MQ", NA))

#Species
print("Species-level HQ MAGs")
hq_mq_mags %>% filter(draft_quality=="HQ") %>% pull(classification) %>% str_remove("d__.*;s__") 

print("Species-level MQ MAGs")
hq_mq_mags %>% filter(draft_quality=="MQ") %>% pull(classification) %>% str_remove("d__.*;s__") 

#Break down classification to phylum, family, genus and species
hq_mq_mags$phylum <- str_remove(str_remove(hq_mq_mags$classification, "d__.*;p__"), ";c_.*")
hq_mq_mags$family <- str_remove(str_remove(hq_mq_mags$classification, "d__.*;f__"), ";g_.*")
hq_mq_mags$genus <- str_remove(str_remove(hq_mq_mags$classification, "d__.*;g__"), ";s_.*")
hq_mq_mags$species <- str_remove(hq_mq_mags$classification, "d__.*;s__")

#Add info of the host subspecies each genome was isolated from
#Get metadata table
metadata <- readRDS("../T3_community-level/metadata.RDS")

hq_mq_mags$host_subspecies <- metadata$Spec.subspecies[match(str_remove(hq_mq_mags$user_genome, ".[0-9]+$"), rownames(metadata))]

#Write out table
write.table(species_mags[,c("user_genome", "classification", "fastani_reference", "note", "warnings",
                                "phylum", "family", "genus", "species", "host_subspecies", "draft_quality")],
                                sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE,
            file = "/proj/sllstore2017021/nobackup/MARKELLA/M6_MAGstats/species_mags.txt")

#Change names to use in the tree
make.unique.2 = function(x, sep='.'){
    ave(x, x, FUN=function(a){if(length(a) > 1){paste(a, 1:length(a), sep=sep)} else {a}})
}

hq_mq_mags$tree_names <- make.unique.2(ifelse(str_remove(hq_mq_mags$classification, ".*s__")!="", str_remove(hq_mq_mags$classification, ".*s__"),
                                        ifelse(str_remove(str_remove(hq_mq_mags$classification, ".*g__"), ";s__")!="", paste("undescribed", str_remove(str_remove(hq_mq_mags$classification, ".*g__"), ";s__"), sep=" "), 
                                               paste("undescribed", str_remove(str_remove(hq_mq_mags$classification, ".*f__"), ";g__.*;s__")))))
                                 
#### Trees ####
hq_mq_tree <- read.tree("/crex/proj/sllstore2017021/nobackup/MARKELLA/M4_MAGtaxonomy/hq_mq.bac120.classify.tree")

#Extract all access number - they will be used to the the species name
mag_tree_accnum <- hq_mq_tree$tip.label
mag_tree_accnum <- mag_tree_accnum[which(grepl("[A-Z][A-Z]_", mag_tree_accnum))]
mag_tree_accnum <- str_sub(mag_tree_accnum, 4)

write.table(mag_tree_accnum, file = "/crex/proj/sllstore2017021/nobackup/MARKELLA/M6_MAGstats/mag_tree_accnum.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE)

#system("sbatch slurm_accnum_to_species.sh -M snowy")

#Read output
mag_tree_names <- read.table(file = "/crex/proj/sllstore2017021/nobackup/MARKELLA/M6_MAGstats/species_names.txt",
                              header = FALSE, sep="\t")


#Create function to change the tree labels to a more readable version
better_names <- function(tree, mag_table) {
  #get tip labels from tree
  tiplab <- tree$tip.label
  #empty vector to fill in
  tiplab <- as.data.frame(tiplab)
  tiplab$newlab <- NA
  for (i in 1:nrow(tiplab)) {
    #if it's an accession number
    if (grepl("[A-Z][A-Z]_", tiplab$tiplab[i])) {
      #get taxon name
      name <- mag_tree_names$V2[which(mag_tree_names$V1==str_sub(tiplab$tiplab[i], 4))]
      #if there is no replacement keep old name
      name <- ifelse(identical(name, character(0)), tiplab$tiplab[i], name)
    } else {
      #Get classification from GTDB classification.
      classif <- mag_table$classification[mag_table$user_genome==tiplab$tiplab[i]]
      species <- str_remove(classif, ".*s__")
      genus <- str_remove(str_remove(classif, ".*g__"), ";s__")
      family <- str_remove(str_remove(classif, ".*f__"), ";g__.*;s__")
      #If available, use species name, otherwise "undescribed genusname" or "undescribed familyname"
      name <- ifelse(species!="", species,
                     ifelse(genus!="", paste("undescribed", genus, sep=" "), paste("undescribed", family)))
    }
    #Store new name in a different column
    tiplab$newlab[i] <- name
  }
  #Make unique new names
  tiplab$newlab <- make.unique.2(tiplab$newlab)
  tree <- sub.taxa.label(tree, tiplab)
  return(tree)
}

#Rename tree
hq_mq_tree <- better_names(hq_mq_tree, hq_mq_mags)
save.image()


#This function takes a tree and a list of tips to highlight 
#and collapses all other nodes to aid visualization 
##Arguments
#tree: phylo tree to be simplified
#tips_to_keep: vector of tips to be highlighted
simpler_tree <- function(tree, tips_to_keep) {
  nodmin <- length(tree$tip.label)+1 #lower node number
  nodmax <- nodmin + tree$Nnode - 1 #higher node number
  
  tree_plot <- ggtree(tree) + theme_tree2(plot.margin=margin(6, 150, 6, 6)) + 
    coord_cartesian(clip = 'off') + geom_tiplab()  #Create ggtree object
  
  #Colour tip labels based on host of genome (which gorilla subspecies, or if it is a reference)
  tipcateg <- data.frame(name=tree$tip.label, host_subspecies=hq_mq_mags$host_subspecies[match(tree$tip.label, hq_mq_mags$tree_names)])
  tipcateg$host_subspecies <- ifelse(is.na(tipcateg$host_subspecies), "reference", as.character(tipcateg$host_subspecies)) %>%
                              factor(levels=c("gorilla", "graueri", "beringei", "reference"))
  
  #Edit palette 
  my_palette <- data.frame(category=levels(tipcateg$host_subspecies), colour=append(RColorBrewer::brewer.pal(3, "Set2"), "#999999"))
  tipcateg$colour <- my_palette$colour[match(tipcateg$host_subspecies, my_palette$category)]  

  tree_plot <- tree_plot %<+% tipcateg + geom_tiplab(aes(color=host_subspecies)) +
  scale_colour_manual(values=c(reference="#999999", gorilla="#66C2A5", graueri="#FC8D62", beringei="#8DA0CB"), name="Host subspecies", labels=c("Western gorilla", "Grauer's gorilla", "Mountain gorilla", "Reference"))
  
  #get numbers of tips to keep
  tiplist <- which(tree$tip.label %in% tips_to_keep) 
  for (i in nodmin:((nodmin+nodmax)/2)) {
    #For every node, if no tips are on the list and if it's not a final node, collapse
    if (as.logical(prod(!(Descendants(tree, i, type="all") %in% tiplist))) &
        length(Descendants(tree, i, type="all")) != length(Descendants(tree, i, type="tips"))) {
      tree_plot <- tree_plot %>% collapse(i)
    }
  }
  
  #Get node names to be highlighted (orders)
  highlights <- tree$node.label[which(grepl(":o__", tree$node.label))]
  #remove higher level classifications
  highlights <- highlights[which(!grepl("f__", highlights))]
  #Get node numbers
  highlights <- which(tree$node.label %in% highlights) + nodmin
  #Add highlights one by one
  for (i in 1:length(highlights)) {
    tree_plot <- tree_plot + geom_hilight(node=highlights[i], fill="lemonchiffon2")
    #Also add labels
    nodlab <- tree$node.label[highlights[i]-nodmin]
    nodlab <- str_remove(nodlab, ".*o__")
    tree_plot <- tree_plot + geom_cladelabel(node=highlights[i], offset = 1.5,
                                             label=nodlab, color="coral4")
  }
   
  #Add some more graphics
  tree_plot <- tree_plot + geom_nodepoint(colour="red", size=0.5)+
    theme(legend.title = element_blank(), legend.position = c(0.9, 0.1))
  #Change tree labels to replace . with space
  tree_plot$data$label <- gsub(".", " ", fixed=TRUE, tree_plot$data$label)
  return(tree_plot)
}

#Split tree based on the three different phyla in which the MAGs are classified
table(hq_mq_mags$phylum)

# Actinobacteriota #
table(hq_mq_mags$family[hq_mq_mags$phylum=="Actinobacteriota"])

#Subset the clade with the lowest common ancestor of the Actinobacteriota MAGs
hq_actin <- extract.clade(hq_mq_tree, getMRCA(hq_mq_tree, 
                                           hq_mq_mags$tree_names[hq_mq_mags$phylum=="Actinobacteriota"]))

#Get a simpler tree
hq_ggactin <- simpler_tree(hq_actin, hq_mq_mags$tree_names[hq_mq_mags$phylum=="Actinobacteriota"])

hq_ggactin

#Save
ggsave(hq_ggactin, filename = "/proj/sllstore2017021/nobackup/MARKELLA/M6_MAGstats/hq_mq_mags_actinobacteria.png",
       device = "png", width = 8, height=10)

# Firmicutes #
table(hq_mq_mags$family[hq_mq_mags$phylum=="Firmicutes"])

hq_firm <- extract.clade(hq_mq_tree, getMRCA(hq_mq_tree, hq_mq_mags$tree_names
                                          [hq_mq_mags$phylum=="Firmicutes"]))
hq_ggfirm <- simpler_tree(hq_firm, hq_mq_mags$tree_names[hq_mq_mags$phylum=="Firmicutes"])

hq_ggfirm

#Save
ggsave(hq_ggfirm, filename = "/proj/sllstore2017021/nobackup/MARKELLA/M6_MAGstats/hq_mq_mags_firmicutes.png",
       device = "png", width = 8, height=10)

# Proteobacteria #
table(hq_mq_mags$family[hq_mq_mags$phylum=="Proteobacteria"])

hq_prot <- extract.clade(hq_mq_tree, getMRCA(hq_mq_tree, hq_mq_mags$tree_names[hq_mq_mags$phylum=="Proteobacteria"]))
hq_ggprot <- simpler_tree(hq_prot, hq_mq_mags$tree_names[hq_mq_mags$phylum=="Proteobacteria"])

hq_ggprot

#Save
ggsave(hq_ggprot, filename = "/proj/sllstore2017021/nobackup/MARKELLA/M6_MAGstats/M6_MAGstats/hq_mq_mags_proteobacteria.png",
       device = "png", width = 8, height=10)


#### Smaller trees ####

#Create function to format smaller trees in a similar way, but maintaining all labels
format_tree <- function(tree) {
  nodmin <- length(tree$tip.label)+1 #lower node number
  nodmax <- nodmin + tree$Nnode - 1 #higher node number
  
  tree_plot <- ggtree(tree) + theme_tree2(plot.margin=margin(6, 350, 6, 6)) + 
    coord_cartesian(clip = 'off') + geom_tiplab(size=8)  #Create ggtree object
  
  #Colour tip labels based on host of genome (which gorilla subspecies, or if it is a reference)
  tipcateg <- data.frame(name=tree$tip.label, host_subspecies=hq_mq_mags$host_subspecies[match(tree$tip.label, hq_mq_mags$tree_names)])
  tipcateg$host_subspecies <- ifelse(is.na(tipcateg$host_subspecies), "reference", as.character(tipcateg$host_subspecies)) %>%
                              factor(levels=c("gorilla", "graueri", "beringei", "reference"))
  
  #Edit palette 
  my_palette <- append(RColorBrewer::brewer.pal(3, "Set2"), "#999999") #get colours for the number of subspecies + grey for refs
  names(my_palette) <- c("gorilla","graueri","beringei","reference")
  #Subset palette
  my_palette <- my_palette[which(names(my_palette) %in% unique(tipcateg$host_subspecies))]
  
  tree_plot <- tree_plot %<+% tipcateg + geom_tiplab(aes(color=host_subspecies), size=8) +
  scale_colour_manual(values=my_palette, name="Host subspecies", labels=list("gorilla"="Western gorilla", "graueri"="Grauer's gorilla", "beringei"="Mountain gorilla", "reference"="Reference"))
  
  #Add some more graphics
  tree_plot <- tree_plot + geom_nodepoint(colour="red", size=0.5)+
    theme(legend.title = element_blank(), legend.position = c(2, 0.1))
  
  #Change labels to remove dots
  tree_plot$data$label <- gsub(".", " ", fixed=TRUE, tree_plot$data$label)
  
  return(tree_plot)
}


#Rothia
rothia_tree <- extract.clade(hq_mq_tree, 
                             Ancestors(hq_mq_tree, getMRCA(hq_mq_tree, c("undescribed Rothia.1", "undescribed Rothia.15")), type="parent"))
rothia_ggtree <- format_tree(rothia_tree) + theme(legend.position = "bottom", legend.text=element_text(size=12))
rothia_ggtree

ggsave(rothia_ggtree, filename="rothia_ggtree.png", device="png", width=8, height=11)

#Streptococcus
strept_tree <- extract.clade(hq_mq_tree, 
                             getMRCA(hq_mq_tree, c("undescribed Streptococcus.1", "undescribed Streptococcus.4")))

strept_ggtree <- format_tree(strept_tree) + theme(legend.position = "bottom", legend.text=element_text(size=12))

#collapse large clade
node_to_collapse <- getMRCA(strept_tree, c("Streptococcus criceti HS-6", "Streptococcus porcinus"))
strept_ggtree <- strept_ggtree + geom_point2(aes(subset=(node==node_to_collapse)), shape=5, size=10)

strept_ggtree <- strept_ggtree %>% collapse(node_to_collapse) + annotate("text", x=0.3, y=18, label="Clade containing S. porcinus and S. orisratti", size=6)
#strept_ggtree

ggsave(strept_ggtree, filename="strept_ggtree.png", device="png", width=8, height=9)

mag_tree_grid <-
          plot_grid(strept_ggtree + theme(legend.position = "none"), rothia_ggtree + theme(),
          nrow=2, rel_widths = 9/11)
          
ggsave(mag_tree_grid, file="mag_tree_grid.png", device="png", width=8, height=18)

sessionInfo()

sink(file=NULL)