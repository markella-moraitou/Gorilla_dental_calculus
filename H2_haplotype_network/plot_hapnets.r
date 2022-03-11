library(pegas)
library(data.table)
library(dplyr)
library(tibble)
library(stringr)
library(readxl)

#Based on Zach's script /proj/sllstore3017021/nobackup/ZACH/ORANGS/haplotype_networks/plot_hapnets.r

#Set working directory
setwd(dir="/crex")

#Get subspecies info for DC genomes from metadata
metadata <- readRDS("T3_community-level/metadata.RDS")
subspecies.table <- metadata %>% select(Spec.subspecies)

#Get subspecies info for downloaded genomes
metadata.dl <- fread('downloaded_mtgenomes/gorilla_mt_genomes_to_download.txt', sep = "\t", header = TRUE)
subspecies.table.dl <- metadata.dl %>% select(accession, organism) %>% column_to_rownames("accession")
colnames(subspecies.table.dl) <- "Spec.subspecies"

#Bind the two tables
subspecies.table <- rbind(subspecies.table, subspecies.table.dl) %>% filter(!is.na(Spec.subspecies))

#Change names e.g beringei -> mountain gorilla
subspecies.table[,1] <- subspecies.table[,1] %>% gsub(pattern="gorilla", replacement="Ggg") %>%
                                         gsub(pattern="beringei", replacement="Gbb") %>%
                                         gsub(pattern="graueri", replacement="Gbg")

#Read in alignment, convert individual names to locality names
mt <- read.FASTA(file = "H2_haplotype_network/mt_genomes_10.fa")

#Remove region info from name
names(mt) <- names(mt) %>% str_remove(pattern=":.*")

gethaplos <- haplotype(mt)
haplotable <- stack(setNames(attr(gethaplos, "index"), rownames(gethaplos)))
haplotable <- haplotable[order(haplotable$values),]
haplotable$ID <- names(mt)
haplotable$values <- NULL
colnames(haplotable) <- c('haplotype','sample_id')

#Add previous info about haplotypes, when applicable
downloaded_mt_metadata <- 
  read.table("/crexdownloaded_mtgenomes/duplicate_genomes.txt", sep='\t', header=TRUE) %>%
  rbind(read.table("/crexdownloaded_mtgenomes/sample_per_haplotype.txt", sep='\t', header=TRUE)) %>% unique
  
haplotable$species <- subspecies.table[match(haplotable$sample_id, rownames(subspecies.table)),] %>% unname
haplotable$vdV_haplotype <- downloaded_mt_metadata$mtDNA.haplotype[match(haplotable$sample_id, downloaded_mt_metadata$accession)]

#The only discrepancy:
#XXVII which corresponds to HT_MG_2 and HT_MG_unk

#After assigning haplotypes, record to table so they can be referenced later
write.table(haplotable, file = 'H2_haplotype_network/haplotypes.txt', quote = FALSE, sep = ',', row.names = FALSE)

temp.names <- names(mt)
names(mt) <- subspecies.table[match(temp.names, rownames(subspecies.table)),]

#Calculate haplotypes, create key for coloring in haplotypes by locality, and build minimum spanning network
ind.haplo<-with(stack(setNames(attr(gethaplos, "index"), rownames(gethaplos))), table(hap=ind, individuals=names(mt)[values]))
net <- haploNet(gethaplos)

#Get colours
plotcolours <- ind.haplo %>% rownames
#Get sample ID
plotcolours <- haplotable$sample_id[match(plotcolours, haplotable$haplotype)]
#Get subspecies
plotcolours <- ifelse(plotcolours %in% rownames(metadata),
                      metadata$Spec.subspecies[match(plotcolours, rownames(metadata))] %>% as.character,
                      metadata.dl$organism[match(plotcolours, metadata.dl$accession)])
                      
#Substitute subspecies label with colours
plotcolours[plotcolours %in% c("Ggg", "gorilla")] <- "#63c9ff"
plotcolours[plotcolours %in% c("Gbg", "graueri")] <- "#efbd11"
plotcolours[plotcolours %in% c("Gbb", "beringei")] <- "grey40"

#Draw haplotype network and save to file
png(file = "H2_haplotype_network/haplonet.mito.png")
#setHaploNetOptions(pie.colors.function = plotcolours,
#   haplotype.outer.color = "#CCCC4D",
#   show.mutation = 3, labels = FALSE)
plot(net, size = attr(net, "freq"), scale.ratio = 1, threshold = 0, cex = 0.5, show.mutation = 0, labels = FALSE, pie=ind.haplo, fast = FALSE)
legend("topright", colnames(ind.haplo), col = rainbow(ncol(ind.haplo)), pch = 20)
dev.off()
