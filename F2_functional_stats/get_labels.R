#!/usr/bin/env Rscript

library(tidyverse)
library(stringr)
library(readr)

args = commandArgs(trailingOnly=TRUE)
  #Load abundance table and metadata
  dat<-read_tsv(args[1])
  metadata<-read_tsv("metadata_gorilla_analysis.txt")
  #Get metadata info for each sample in the abundance table
  labels <- colnames(dat)
  labels <- metadata$plot.label[match(str_remove(labels, "_Abundance-RPKs"), metadata$Seq.label)]
  labels <- labels[-1]
write.table(labels,file=args[2],quote=F,sep="\t", row.names=F, col.names=F)