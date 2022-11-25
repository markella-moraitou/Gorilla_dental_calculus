library(tidyverse)
library(ggpubr)
library(lubridate)

meta <- read_tsv("tables_and_stats/S1_Sample_metadata.tsv", skip=1)

frag_df<-read.delim("all_gorilla_frag_length_retained_reads.txt",header = T,sep = " ") %>% rename("Sample_ID"="SampleID")

meta_frag<-meta %>% left_join(frag_df) %>% 
  filter(`Sample type`=="sample") %>% 
  mutate(Sample_Age=year(today())-`Collection year`)

summary(lm(Average.Length.Retained.Reads~Sample_Age,data = meta_frag))


ggscatter(meta_frag, x = "Sample_Age", y = "Average.Length.Retained.Reads",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
          ) +
  stat_cor(method = "pearson") + # Add correlation coefficient
  labs(y="Average Fragment Length",x="Sample Age")+
  # stat_regline_equation()+
  see::theme_modern()
ggsave(plot = last_plot(),filename = "frag_length_vs_samp_age.png",dpi = 300,width = 6,height = 4,units = "in")
