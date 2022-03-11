library(tidyverse)
library(taxize)

# keywords
oral.list="^oral|mouth|plaque|calculus|tartar|saliva|periodont|peridont|tooth|caries|tongue|molar|throat|cheek|tongue"
host.list="gut|faeces|feces|human|gorilla|Gastrointestinal|Organ|Urogenital|Abdomen|Abscess|Airways|Ankle|Arm|Bladder|Bone|Brain|Ear|Eye|Foot|Hand|Head|Heart|Joint|Kidney|intestine|Leg|Liver|Lung|Lymph|Mammary|Nephridium|Nose|Ovary|Phylloplane|Rectum|Rumen|Skin|Nail|Hair|Spinal Cord|Spleen|Stomach|Thoracic|Torso|Trachea|Urethra|Vagina|Vascular system|Wound|Amniotic|Bladder|Blood|Bronchial|Caecal|Cerebrospinal|Digestive|Mucus|Nasopharyngeal|Phlegm|Rumen|Semen|Sputum|Synovial|Tracheal|Urine|Vagina|Vomit"
environmental.list="soil|dirt|humus|sand|gravel|mud|sediment|rock|rhizosphere|grass|earth|^dust$|^air$|filter|water|river|^sea$|^river|^stream|^pond|aquifer|Environmental|Air|Aquatic|Biofilm|Microbial community|Terrestrial"
contaminants="sludge|digestor|feed|Aquaculture|Beverage|Biofilter|City|Clean room|Composting|Container|Contaminant|Cooling|Currency|Dairy|^Decompos|Domestic waste|Dust|Engineer|Fermented|Field|Food|Greenhouse|Heavy metal|House|Indoor|Industrial|Landfill|Livestock|Machine|Meadow|Meat|^Mine$|Oil|Paddy|Painting|^Park$|Pesticide|Plantation|Preserved|Seafood|Sewage sludge|Solid waste|Spoiled|Starter culture|Sterilized|Sugary food|Surface swab|^Waste|reservoir"
missing.list="na|n/a|none|not applicable|not available|not collected|not provided|not_applicable|not_collected|missing|^-|unknown"

# taxa
taxa<-read_csv(dir("/home/adrian/taxdump", full.names=T, pattern="ncbi_lineages_*"),na = "NA",
               col_types = cols(tax_id = col_character(),superkingdom = col_character(),phylum = col_character(),class = col_character(),order = col_character(),family = col_character(),genus = col_character(),species = col_character(),biotype = col_character(),clade = col_character(),clade1 = col_character(),
                                clade10 = col_character(),clade11 = col_character(),clade12 = col_character(),clade13 = col_character(),clade14 = col_character(),clade15 = col_character(),clade16 = col_character(),clade17 = col_character(),clade18 = col_character(),clade19 = col_character(),clade2 = col_character(),
                                clade3 = col_character(),clade4 = col_character(),clade5 = col_character(),clade6 = col_character(),clade7 = col_character(),clade8 = col_character(),clade9 = col_character(),cohort = col_character(),forma = col_character(),`forma specialis` = col_character(),`forma specialis1` = col_character(),
                                genotype = col_character(),infraclass = col_character(),infraorder = col_character(),isolate = col_character(),kingdom = col_character(),morph = col_character(),`no rank` = col_character(),`no rank1` = col_character(),`no rank2` = col_character(),`no rank3` = col_character(),`no rank4` = col_character(),
                                parvorder = col_character(),pathogroup = col_character(),section = col_character(),series = col_character(),serogroup = col_character(),serotype = col_character(),`species group` = col_character(),`species subgroup` = col_character(),strain = col_character(),subclass = col_character(),subcohort = col_character(),
                                subfamily = col_character(),subgenus = col_character(),subkingdom = col_character(),suborder = col_character(),subphylum = col_character(),subsection = col_character(),subspecies = col_character(),subtribe = col_character(),subvariety = col_character(),superclass = col_character(),superfamily = col_character(),superorder = col_character(),
                                superphylum = col_character(),tribe = col_character(),varietas = col_character())) 

# taxa lists 
genera<-c("Bacillus","Corynebacterium","Erysipelothrix","Erysipelothrix","Exiguobacterium","Limosilactobacillus","Lancefieldella","Lautropia","Macrococcus","Mailhella","Neisseria","Phocaeicola","Propionibacterium","Propionibacterium","Propionimicrobium","Proteus","Providencia","Pseudogracilibacillus","Rothia","Sporosarcina","Vagococcus","Vagococcus","Veillonella")
families<-c("Desulfovibrionaceae","Planococcaceae","Eggerthellaceae","Eggerthellaceae","Actinomycetaceae","Paludibacteraceae","Ruminococcaceae","Saccharimonadaceae","Bacteroidetes","Thermoguttaceae")

# read in custom bacdive database
bacdive<-read_csv("I1_isolation-sources/bacdive-isolation-sources.csv") %>% 
  group_by(ID) %>% 
  fill(Species:Continent,.direction="down") %>% # fill down the empty cells with the correct metadata
  pivot_longer(cols = contains("Category"),names_to = "level",values_to = "source.info")

bacdive.wgroup<-bacdive %>% group_by(source,classification,ID) %>% 
  mutate(isolation_source_category = ifelse(sum(grepl(pattern = oral.list,ignore.case = T,perl = F,source.info))>0,"oral",
                                            ifelse(sum(grepl(pattern=host.list,ignore.case = T,perl = F,source.info))>0,"host",
                                            ifelse(sum(grepl(pattern = environmental.list, ignore.case = T,perl = F,source.info))>0,"environmental",
                                                   ifelse(sum(grepl(pattern = contaminants,ignore.case = T,perl = F,source.info))>0, "contaminant",
                                                          ifelse(sum(grepl(pattern = missing.list,ignore.case = T,perl = F,source.info))>0,"missing","unknown")))))) %>% 
  pivot_wider(id_cols = c(source,classification,ID,isolation_source_category),names_from = "level",values_from = "source.info") %>% 
  mutate(level=ifelse(classification%in%genera,"genus","family"))

# add in HOMD
homd<-read_tsv("T3_community-level/homd_taxonomy_table.txt") %>% 
  mutate(NCBI_taxon_id=as.character(NCBI_taxon_id)) %>% 
  pull(NCBI_taxon_id) %>% 
  taxize::ncbi_get_taxon_summary(.) %>%
  # filter(rank %in% c("family","genus","species")) %>%
  mutate(classification=ifelse(rank=="species",gsub(" .*", "", name),ifelse(rank=="genus",name,ifelse(rank=="family",name,"")))) %>% 
  mutate(isolation_source_category="oral",source="homd") %>% 
  rename(tax_id=uid)

# and hominid
hominid<-read_csv("I1_isolation-sources/fellowsyates2021.csv") %>% 
  mutate(tax_id=as.character(tax_id)) %>% 
  pull(tax_id) %>% 
  taxize::ncbi_get_taxon_summary(.) %>% 
  mutate(classification=ifelse(rank=="species",gsub(" .*", "", name),ifelse(rank=="genus",name,ifelse(rank=="family",name,"")))) %>% 
  mutate(isolation_source_category="oral",source="hominid")

# and NCBI biosample info
ncbi<-read_tsv("I1_isolation-sources/source_results.tsv") %>% 
  mutate(taxid=as.character(taxid)) %>% 
  filter(attribute=="isolation_source") %>% 
  rename(tax_id=taxid) %>% 
  distinct(value,biosample,.keep_all = T) %>% 
  left_join(taxa,by="tax_id") %>% 
  select(value:genus) %>% 
  pivot_longer(cols = phylum:genus,names_to = "class.level",values_to = "classification") %>% 
  mutate(isolation_source_category = ifelse(grepl(pattern = oral.list,ignore.case = T,perl = F,value),"oral",
                                            ifelse(grepl(pattern=host.list,ignore.case = T,perl = F,value),"host",
                                            ifelse(grepl(pattern = environmental.list, ignore.case = T,perl = F,value),"environmental",
                                                   ifelse(grepl(pattern = contaminants,ignore.case = T,perl = F,value), "contaminant",
                                                          ifelse(grepl(pattern = missing.list,ignore.case = T,perl = F,value),"missing","unknown"))))),
         source="ncbi") 
  
# set palette
pal<-c(ggpubr::get_palette(palette = "Paired",k = 5),"gray",NA)
cols<-c(pal[4],pal[3],pal[2],pal[5],pal[6],NA)

comb.df<-bacdive.wgroup %>% 
  # distinct(source,classification,isolation_source_category) %>% 
  bind_rows(homd) %>% 
  bind_rows(hominid) %>%
  bind_rows(ncbi) %>%
  filter(classification%in%genera|classification%in%families) %>% 
  mutate(count=1) %>% 
  group_by(classification) %>% 
  mutate(total=sum(count)) %>%
  group_by(classification,isolation_source_category,total) %>% 
  summarise(rel=sum(count)/total,
            level=ifelse(classification%in%genera,"Genus-level Classification","Family-level Classification")) %>% 
  distinct(classification,isolation_source_category,.keep_all = T)

comb.df<-bind_rows(comb.df,data.frame(classification="Saccharimonadaceae",level="Family-level Classification",
                                      isolation_source_category=c("oral","environmental","host","contaminant","missing"),
                                      total=62,
                                      rel=c(0.064516129,0.758064516,0.129032258,0.032258065,0.016129032)))

oral.order<-comb.df %>% 
  filter(isolation_source_category=="oral") %>% 
  group_by(classification) %>% 
  summarise(max.oral=max(rel,na.rm = T)) %>% 
  arrange(-max.oral)

labels<-read_tsv("M1_mag/mag-table.tsv",col_names = T)

comb.df %>% 
  ungroup() %>% 
  left_join(oral.order) %>% 
  left_join(labels) %>% 
  mutate(isolation_source_category=forcats::fct_relevel(isolation_source_category,c("oral","host","environmental","contaminant","unknown"))) %>% 
  ggplot(aes(x=reorder(label,-max.oral),y=rel,fill=isolation_source_category))+
  geom_bar(stat = "identity")+
  guides(fill=guide_legend(title="Isolation Source Category"))+
  scale_y_reverse(labels=rev(c(0.00,0.25,0.50,0.75,1.00)))+
  scale_fill_manual(values = cols)+
  facet_grid(.~level,scales = "free",space = "free")+
  labs(y="Ratio of Identified Sources")+
  geom_text(aes(label=round(rel*total,digits = 0)),size = 3, position = position_stack(vjust = 0.5)) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 75,hjust = 1,size=12),
        axis.title.x = element_blank(),
        strip.text = element_text(size=14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size=12),
        legend.position = "top")
ggsave(plot = last_plot(),filename = "I1_isolation-sources/isolation-source-proportions.png",dpi=300,height = 8,width = 12,units = "in")
