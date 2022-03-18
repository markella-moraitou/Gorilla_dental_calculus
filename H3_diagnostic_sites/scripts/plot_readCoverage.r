#Adapted from /proj/SNIC_PROJECT/nobackup/SAM/Bear_Lineage_Project/scripts/plot_readCoverage.r

#load libraries
library(ggplot2)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
#read in data and give the columns names
data <- read_table2(args[1], col_names = FALSE)

colnames(data) <- c("pos", "count", "allele")


##Read in reference data with pos, western allele and eastern allele (you have something like this, right?)
##The commands below will work on a table with the columns 'pos','allele_west','allele_east'. You might have to
## modify the code/tables accordingly
diagnosticReference <- read_table2("diagnostic_sites_east_west.txt")
colnames(diagnosticReference)=c("pos","Eastern","Western")

#Vector with diagnostic sites positionss from this table:
diagnosticSites <- diagnosticReference$pos

#Extract diagnostic sites from table:
diagnosticData <- data[which(!is.na(match(data$pos, diagnosticSites))),]

#Add a column specifying whether the allele in question is western or eastern, starting from a column with NA-values
diagnosticData$Population <- "NA"
#Then loop through the table and add eastern or western this
for (i in 1:length(diagnosticData$Population)){
  western_allele <- diagnosticReference[which(diagnosticReference$pos == diagnosticData$pos[i]), 'Western']
  eastern_allele <- diagnosticReference[which(diagnosticReference$pos == diagnosticData$pos[i]), 'Eastern']
  if (diagnosticData$allele[i] == western_allele){diagnosticData$Population[i] <- "Western"}
  else if (diagnosticData$allele[i] == eastern_allele){diagnosticData$Population[i] <- "Eastern"}
  else {diagnosticData$Population[i] <- "Other"}
}

#Simplest possible plot from this, colouring the bars by population.
pdf(file=paste0(str_remove(args[1], "_sort_MT_RG.allele_counts"), ".pdf"))
ggplot(diagnosticData, aes(x=as.character(pos), y=count, fill=Population)) + geom_col() + 
  ggtitle(str_remove(args[1], "_sort_MT_RG.allele_counts")) + 
  theme(axis.text.x=element_text(angle=45, hjust =1, size=6))+
  xlab("Diagnostic positions") + ylab("Read depth")

boxplot(count~Population,data=diagnosticData, main=str_remove(args[1], "_sort_MT_RG.allele_counts"),
        xlab="Lineage",
        ylab="Reads at diagnostic sites")
        
tapply(X=diagnosticData$count, INDEX=list(diagnosticData$Population), FUN=mean)
dev.off()