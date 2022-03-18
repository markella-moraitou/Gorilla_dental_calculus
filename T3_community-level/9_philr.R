#Start logging
sink(file = "log9_philr.txt")

#load packages
library(vegan)
library(taxize)
library(phyloseq)
library(ape)
library(philr)
library(ggplot2)
library(glmnet)

load(".RData")

#Community-level taxonomic analysis - Script 9
#Philr normalization - ordination and permanova

#### Build classification object and tree ####
#Create classification object
taxonomy_for_philr <- classification(sci_id = taxa_names(spe_data_final), db='ncbi')

#Construct tree based on taxonomy
species_tree <-  class2tree(taxonomy_for_philr)

#Some checks
print("Is the tree rooted?")
is.rooted(species_tree$phylo)
print("Is the tree binary?")
is.binary(species_tree$phylo)

#Save an appropriate class object for the tree
tree <- multi2di(species_tree$phylo)

#Name the internal nodes
tree <- makeNodeLabel(tree, method="number", prefix='n')

#Rename tips so they are labelled with the taxonomic ID instead of species names (for consistency with the phyloseq object)
tree$tip.label <- species_tree$names

#Copy classification matrix
tree$classification <- species_tree$classification

#
print("Which names still don't match?")
which(!(taxa_names(spe_data_final) %in% tree$tip.label))

#### Run philr ####
#Add tree to final phyloseq object
phy_tree(spe_data_final) <- tree
phy_tree(spe_data_final_norm) <- tree

#Get transposed OTU table from final phyloseq
species_table_final <- t(otu_table(spe_data_final))
species_table_final <- as.matrix(species_table_final)

#Philr normalization, after introducing a pseudocount of 0.5
species_philr <- philr((species_table_final + 0.5), tree = tree)

#Calculate philr distances
species_philr_dist <- dist(species_philr, method="euclidean")

#### Ordination and permanova ####
#Plot ordination
philr <- ordinate(spe_data_final, 'PCoA', distance=species_philr_dist)

pdf(file = "T3_community-level/philr.pdf")
plot_ordination(spe_data_final, philr, color="Spec.subspecies", shape="Seq.centre", title="PCoA on philr distances based on species level assignments")
dev.off()

#
print("PERMANOVA using a philr dissimilarity matrix")
model1_philr <- adonis(species_philr ~ sample_data(spe_data_final)$readcount.m.before.Kraken + sample_data(spe_data_final)$Spec.subspecies + sample_data(spe_data_final)$Seq.centre, permutations = 10000, method = "euclidean")
model1_philr

model2_philr <- adonis(species_philr ~ sample_data(spe_data_final)$readcount.m.before.Kraken + sample_data(spe_data_final)$Seq.centre + sample_data(spe_data_final)$Spec.subspecies, permutations = 10000, method = "euclidean")
model2_philr

#Write model output as files
write.table(as.data.frame(model1_philr[1:1]), "T3_community-level/T3_community-level/model1_philr.txt", sep = "\t", quote = FALSE, row.names = TRUE, dec=",")
write.table(as.data.frame(model2_philr[1:1]), "T3_community-level/T3_community-level/model2_philr.txt", sep = "\t", quote = FALSE, row.names = TRUE, dec=",")

#### Logistic regression ####

#Create a factor that indicates if sample is mountain or not
is.mountain <- as.factor(sample_data(spe_data_final)$Spec.subspecies=="beringei")

#Run logistic regression using the philr transformed otu table
mountain_logreg <- glmnet(species_philr, is.mountain, alpha=1, family="binomial")

#Extract coefficents
mountain_coords <- as.matrix(coefficients(mountain_logreg, s=0.2526))
#Keep the balances (nodes?) with a coefficient different than 0
mountain_coords <- rownames(mountain_coords)[which(mountain_coords != 0)]
mountain_coords <- mountain_coords[2:length(mountain_coords)]
mountain_node_names <- name.balance(phy_tree(spe_data_final), tax_table(spe_data_final), "n1066")
mountain_node_names
"species_Pseudoalteromonas ruthenica/species_Pseudoalteromonas sp. MT33b"

#Same process for Grauer's gorillas
#Create a factor that indicates if sample is grauer or not
is.grauer <- as.factor(sample_data(spe_data_final)$Spec.subspecies=="graueri")

#Run logistic regression using the philr transformed otu table
grauer_logreg <- glmnet(species_philr, is.grauer, alpha=1, family="binomial")

#Extract coefficents
grauer_coords <- as.matrix(coefficients(grauer_logreg, s=0.2526))
#Keep the balances (nodes?) with a coefficient different than 0
grauer_coords <- rownames(grauer_coords)[which(grauer_coords != 0)]
#It is just the intercept

#Same process for Western gorillas
#Create a factor that indicates if sample is grauer or not
is.western <- as.factor(sample_data(spe_data_final)$Spec.subspecies=="gorilla")

#Run logistic regression using the philr transformed otu table
western_logreg <- glmnet(species_philr, is.western, alpha=1, family="binomial")

#Extract coefficents
western_coords <- as.matrix(coefficients(western_logreg, s=0.2526))
#Keep the balances (nodes?) with a coefficient different than 0
western_coords <- rownames(western_coords)[which(western_coords != 0)]
#It is just the intercept

#Stop logging
sink(file = NULL)

save.image()