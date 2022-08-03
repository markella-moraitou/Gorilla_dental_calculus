#!/usr/bin/env Rscript

#Based on SNIC_PROJECT/scripts/functional/sortSigOutput.R

require(tidyverse)

args = commandArgs(trailingOnly=TRUE)

dat<-read.delim(args[1],comment="",header=T) %>%
  #chose to filter by Q-value instead of P-value
  filter(Q.value<=0.05) %>%
  separate(col=`X..Feature`,into=c("feature","description"),sep=": ") %>%
  select(feature,description)

write.table(dat,file=args[2],quote=F,sep="\t")
