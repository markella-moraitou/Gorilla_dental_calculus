#!/usr/bin/env Rscript

#Based on /proj/SNIC_PROJECT/nobackup/ADRIAN/scripts/functional/addMetadata.R 

require(tidyverse)


args = commandArgs(trailingOnly=TRUE)

  dat<-read.table(args[1],header=T,comment.char = " ")
  subsp_dat<-read.table("metadata_gorilla_analysis.txt", header=T, sep="\t")
  dat.m<- dat %>%
    separate(col = 1,into=c("Seq.label","file"), sep="_") %>%
    #changed left_join to right_join because some of the original samples in the metadata table aren't being used
    right_join(subsp_dat,by="Seq.label") %>%
    mutate(SampleID=paste0(Seq.label,"_",file)) %>%
    select(plot.label) %>%
    t()
  write.table(cbind("plot.label",dat.m),file = args[2],row.names = F,quote = F,col.names = F,sep = "\t")
