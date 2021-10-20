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

#MAG taxonomy and trees


#### Load GTDB-Tk summaries ####

hq_mags <- read.table("/proj/sllstore2017021/nobackup/MARKELLA/M4_MAGtaxonomy/hq.bac120.summary.tsv", sep="\t", comment.char="", header=T)
mq_mags <- read.table("/proj/sllstore2017021/nobackup/MARKELLA/M4_MAGtaxonomy/mq.bac120.summary.tsv", sep="\t", comment.char="", header=T)


#Species
print("Species-level HQ MAGs")
str_remove(hq_mags$classification, "d__.*;s__") 

print("Species-level MQ MAGs")
str_remove(mq_mags$classification, "d__.*;s__")

#Break down classification to phylum, family, genus and species
hq_mags$phylum <- str_remove(str_remove(hq_mags$classification, "d__.*;p__"), ";c_.*")
hq_mags$family <- str_remove(str_remove(hq_mags$classification, "d__.*;f__"), ";g_.*")
hq_mags$genus <- str_remove(str_remove(hq_mags$classification, "d__.*;g__"), ";s_.*")
hq_mags$species <- str_remove(hq_mags$classification, "d__.*;s__")

mq_mags$phylum <- str_remove(str_remove(mq_mags$classification, "d__.*;p__"), ";c_.*")
mq_mags$family <- str_remove(str_remove(mq_mags$classification, "d__.*;f__"), ";g_.*")
mq_mags$genus <- str_remove(str_remove(mq_mags$classification, "d__.*;g__"), ";s_.*")
mq_mags$species <- str_remove(mq_mags$classification, "d__.*;s__")

#Add info of the host subspecies each genome was isolated from
#Get metadata table
metadata <- readRDS("../T3_community-level/metadata.RDS")

hq_mags$host_subspecies <- metadata$Spec.subspecies[match(str_remove(hq_mags$user_genome, ".[1-9]$"), rownames(metadata))]
mq_mags$host_subspecies <- metadata$Spec.subspecies[match(str_remove(mq_mags$user_genome, ".[1-9]$"), rownames(metadata))]

#Species-level MAGs
species_mags <- rbind(hq_mags[,
                              c("user_genome", "classification", "fastani_reference", "note", "warnings",
                                "phylum", "family", "genus", "species", "host_subspecies")],
                      mq_mags[,
                              c("user_genome", "classification", "fastani_reference", "note", "warnings",
                                "phylum", "family", "genus", "species", "host_subspecies")])

#Add info about quality
species_mags$draft_quality <- c(rep("HQ", nrow(hq_mags)),
                                rep("MQ", nrow(mq_mags)))

#Change names to use in the tree
hq_mags$tree_names <- make.names(ifelse(str_remove(hq_mags$classification, ".*s__")!="", str_remove(hq_mags$classification, ".*s__"),
                                        ifelse(str_remove(str_remove(hq_mags$classification, ".*g__"), ";s__")!="", paste("undescribed", str_remove(str_remove(hq_mags$classification, ".*g__"), ";s__"), sep=" "), 
                                               paste("undescribed", str_remove(str_remove(hq_mags$classification, ".*f__"), ";g__.*;s__")))),
                                 unique = TRUE)
                                 
mq_mags$tree_names <- make.names(ifelse(str_remove(mq_mags$classification, ".*s__")!="", str_remove(mq_mags$classification, ".*s__"),
                                        ifelse(str_remove(str_remove(mq_mags$classification, ".*g__"), ";s__")!="", paste("undescribed", str_remove(str_remove(mq_mags$classification, ".*g__"), ";s__"), sep=" "), 
                                               paste("undescribed", str_remove(str_remove(mq_mags$classification, ".*f__"), ";g__.*;s__")))),
                                 unique = TRUE)

#Write out table
write.table(species_mags, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE,
            file = "/proj/sllstore2017021/nobackup/MARKELLA/M6_MAGstats/species_mags.txt")

#Small table that gives classification by user genomes
hq_mags[, c("user_genome","host_subspecies", "family", "genus", "species")] %>%
  write.table(sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE,
              file = "/proj/sllstore2017021/nobackup/MARKELLA/M6_MAGstats/hq_MAG_classification.txt")

mq_mags[, c("user_genome","host_subspecies", "family", "genus", "species")] %>%
  write.table(sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE,
              file = "/proj/sllstore2017021/nobackup/MARKELLA/M6_MAGstats/mq_MAG_classification.txt")

#### Trees ####
hq_tree <- read.tree("/crex/proj/sllstore2017021/nobackup/MARKELLA/M4_MAGtaxonomy/hq.bac120.classify.tree")
mq_tree <- read.tree("/crex/proj/sllstore2017021/nobackup/MARKELLA/M4_MAGtaxonomy/mq.bac120.classify.tree")

#Extract all access number - they will be used to the the species name
mag_tree_accnum <- unique(append(hq_tree$tip.label, mq_tree$tip.label))
mag_tree_accnum <- mag_tree_accnum[which(grepl("[A-Z][A-Z]_", mag_tree_accnum))]
mag_tree_accnum <- str_sub(mag_tree_accnum, 4)

write.table(mag_tree_accnum, file = "/crex/proj/sllstore2017021/nobackup/MARKELLA/M6_MAGstats/mag_tree_accnum.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE)

system("sbatch slurm_accnum_to_species.sh -M snowy")

#Read output
mag_tree_names <- read.table( file = "/crex/proj/sllstore2017021/nobackup/MARKELLA/M6_MAGstats/species_names.txt",
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
      classif <- mag_table$classification[mq_mags$user_genome==tiplab$tiplab[i]]
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
  tiplab$newlab <- make.names(tiplab$newlab, unique = TRUE)
  tree <- sub.taxa.label(tree, tiplab)
  return(tree)
}

#Rename tree
hq_tree <- better_names(hq_tree, hq_mags)
mq_tree <- better_names(mq_tree, mq_mags)
save.image()

##Focus on HQ MAGs##

#This function takes a tree and a list of tips to highlight 
#and collapses all other nodes to aid visualization 
##Arguments
#tree: phylo tree to be simplified
#tips_to_keep: vector of tips to be highlighted
simpler_tree <- function(tree, tips_to_keep) {
  nodmin <- length(tree$tip.label)+1 #lower node number
  nodmax <- nodmin + tree$Nnode - 1 #higher node number
  
  tree_plot <- ggtree(tree) + theme_tree2(plot.margin=margin(6, 100, 6, 6)) + 
    coord_cartesian(clip = 'off') + geom_tiplab()  #Create ggtree object
  
  #Colour tip labels based on host of genome (which gorilla subspecies, or if it is a reference)
  tipcateg <- data.frame(name=tree$tip.label, host_subspecies=hq_mags$host_subspecies[match(tree$tip.label, hq_mags$tree_names)])
  tipcateg$host_subspecies <- ifelse(is.na(tipcateg$host_subspecies), "reference", as.character(tipcateg$host_subspecies)) %>%
                              factor(levels=c("gorilla", "graueri", "beringei", "reference"))
  
  #Edit palette 
  my_palette <- append(RColorBrewer::brewer.pal(3, "Set2"), "#999999")

  tree_plot <- tree_plot %<+% tipcateg + geom_tiplab(aes(color=host_subspecies)) +
  scale_colour_manual(values=my_palette, name="Host subspecies", labels=c("Western gorilla", "Grauer's gorilla", "Mountain gorilla", "Reference"))
  
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
  
  #get numbers of tips to keep
  tiplist <- which(tree$tip.label %in% tips_to_keep) 
  for (i in nodmin:((nodmin+nodmax)/2)) {
    #For every node, if no tips are on the list and if it's not a final node, collapse
    if (as.logical(prod(!(Descendants(tree, i, type="all") %in% tiplist))) &
        length(Descendants(tree, i, type="all")) != length(Descendants(tree, i, type="tips"))) {
      tree_plot <- tree_plot %>% collapse(i)
    }
  }
  
  #Add some more graphics
  tree_plot <- tree_plot + geom_nodepoint(colour="red", size=0.5)+
    theme(legend.title = element_blank(), legend.position = c(0.9, 0.1))
  #Change tree labels to replace . with space
  tree_plot$data$label <- gsub(".", " ", fixed=TRUE, tree_plot$data$label)
  return(tree_plot)
}

save.image()


#Split tree based on the three different phyla in which the MAGs are classified
table(hq_mags$phylum)

# Actinobacteriota #
table(hq_mags$family[hq_mags$phylum=="Actinobacteriota"])

#Subset the clade with the lowest common ancestor of the Actinobacteriota MAGs
hq_actin <- extract.clade(hq_tree, getMRCA(hq_tree, 
                                           hq_mags$tree_names[hq_mags$phylum=="Actinobacteriota"]))

#Get a simpler tree
hq_ggactin <- simpler_tree(hq_actin, hq_mags$tree_names[hq_mags$phylum=="Actinobacteriota"])

hq_ggactin

#Save
ggsave(hq_ggactin, filename = "/proj/sllstore2017021/nobackup/MARKELLA/M6_MAGstats/HQ_MAGs_actinobacteria.png",
       device = "png", width = 8, height=10)

# Firmicutes #
table(hq_mags$family[hq_mags$phylum=="Firmicutes"])

hq_firm <- extract.clade(hq_tree, getMRCA(hq_tree, hq_mags$tree_names
                                          [hq_mags$phylum=="Firmicutes"]))
hq_ggfirm <- simpler_tree(hq_firm, hq_mags$tree_names[hq_mags$phylum=="Firmicutes"])

hq_ggfirm

#Save
ggsave(hq_ggfirm, filename = "C://Users//MARKELLA//OneDrive - Uppsala universitet//Degree Project//Bioinformatics//M6_MAGstats//HQ_MAGs_firmicutes.png",
       device = "png", width = 8, height=10)

# Proteobacteria #
table(hq_mags$family[hq_mags$phylum=="Proteobacteria"])

hq_prot <- extract.clade(hq_tree, getMRCA(hq_tree, hq_mags$tree_names
                                          [hq_mags$phylum=="Proteobacteria"]))
hq_ggprot <- simpler_tree(hq_prot, hq_mags$tree_names[hq_mags$phylum=="Proteobacteria"])

hq_ggprot

#Save
ggsave(hq_ggprot, filename = "C://Users//MARKELLA//OneDrive - Uppsala universitet//Degree Project//Bioinformatics//M6_MAGstats//HQ_MAGs_proteobacteria.png",
       device = "png", width = 8, height=10)


#### Smaller trees ####

#Create function to format smaller trees in a similar way, but maintaining all labels
format_tree <- function(tree) {
  nodmin <- length(tree$tip.label)+1 #lower node number
  nodmax <- nodmin + tree$Nnode - 1 #higher node number
  
  tree_plot <- ggtree(tree) + theme_tree2(plot.margin=margin(6, 300, 6, 6)) + 
    coord_cartesian(clip = 'off') + geom_tiplab(size=8)  #Create ggtree object
  
  #Colour tip labels based on host of genome (which gorilla subspecies, or if it is a reference)
  tipcateg <- data.frame(name=tree$tip.label, host_subspecies=hq_mags$host_subspecies[match(tree$tip.label, hq_mags$tree_names)])
  tipcateg$host_subspecies <- ifelse(is.na(tipcateg$host_subspecies), "reference", as.character(tipcateg$host_subspecies)) %>%
                              factor(levels=c("gorilla", "graueri", "beringei", "reference"))
  
  #Edit palette 
  my_palette <- append(RColorBrewer::brewer.pal(3, "Set2"), "#999999")

  tree_plot <- tree_plot %<+% tipcateg + geom_tiplab(aes(color=host_subspecies), size=8) +
  scale_colour_manual(values=my_palette, name="Host subspecies", labels=c("Western gorilla", "Grauer's gorilla", "Mountain gorilla", "Reference"))
  
  #Add some more graphics
  tree_plot <- tree_plot + geom_nodepoint(colour="red", size=0.5)+
    theme(legend.title = element_blank(), legend.position = c(2, 0.1))
  
  #Change labels to remove dots
  tree_plot$data$label <- gsub(".", " ", fixed=TRUE, tree_plot$data$label)
  
  return(tree_plot)
}


#Rothia
rothia_tree <- extract.clade(hq_tree, 
                             getMRCA(hq_tree, c("undescribed.Rothia", "undescribed.Rothia.3")))
rothia_ggtree <- format_tree(rothia_tree) + theme(legend.position = "none")
rothia_ggtree

#Olsenella
olsenella_tree <- extract.clade(hq_tree, 
                                length(hq_tree$tip.label) + which(grepl("0.873:g__Olsenella", hq_tree$node.label)))
olsenella_ggtree <- format_tree(olsenella_tree) + theme(legend.position = "none")
olsenella_ggtree

#Eggerthellaceae
eggerth_tree <- extract.clade(hq_tree, 
                              getMRCA(hq_tree, c("Denitrobacterium.detoxificans", "Eggerthella.sp..CAG.1427")))
eggerth_ggtree <- format_tree(eggerth_tree) + theme(legend.position = "none")
eggerth_ggtree + theme_tree2(plot.margin=margin(6, 400, 6, 6), legend.position = "none")
