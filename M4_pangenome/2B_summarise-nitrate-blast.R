library(tidyverse)
species<-c("Rothia","Neisseria","Veillonella")

for(s in species){
filenames <- Sys.glob(paste0(s,"/*.txt"))

# now read them all in
Data <- lapply(filenames, function(.file){
  df<-read_delim(.file,delim="\t",col_names=F,show_col_types=F) %>% 
    mutate(file.id=as.character(.file))
  })

# combine into a single dataframe
df <- do.call(bind_rows, Data)

df %>% mutate(type=ifelse(grepl("sorted",file.id), "MAGs", "REFs")) %>%
  # keep only the highest scoring hits
  group_by(file.id,type,X1) %>%
  rename("sample"=X1)
  # filter(X6==max(X6)) %>%
  # mutate(p.ident.completeness=X7/(X13[grepl("genomic",sample)]-X12[grepl("genomic",sample)]+1))

  write.csv(paste0("SNIC_PROJECT/ADRIAN/calculus/Gorilla/mag/nitrate_genes/blast/all-",s,".csv"),row.names=FALSE)
}
