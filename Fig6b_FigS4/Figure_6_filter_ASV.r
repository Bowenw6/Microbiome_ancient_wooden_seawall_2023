rm(list=ls())

library(Biostrings)
library(tibble)

# import 16S abundance and tax data 
ASV_all_16S_With_Tax<-read.table("ASV_Abundance_16S_split.txt",sep = "\t", header = TRUE,row.names = 1)
ASV_all_16S<-ASV_all_16S_With_Tax[,1:(ncol(ASV_all_16S_With_Tax)-7)]
ASV_all_16S_Tax<-ASV_all_16S_With_Tax[,(ncol(ASV_all_16S_With_Tax)-6):ncol(ASV_all_16S_With_Tax)]

# import ITS abundance and tax data 
ASV_all_ITS_With_Tax<-read.table("ASV_Abundance_ITS_split.txt",sep = "\t", header = TRUE,row.names = 1)
ASV_all_ITS<-ASV_all_ITS_With_Tax[,1:(ncol(ASV_all_ITS_With_Tax)-7)]
ASV_all_ITS_Tax<-ASV_all_ITS_With_Tax[,(ncol(ASV_all_ITS_With_Tax)-6):ncol(ASV_all_16S_With_Tax)]


ASV_filtered_16S<- ASV_all_16S[rowSums(ASV_all_16S)>ncol(ASV_all_16S),]
ASV_filtered_ITS<- ASV_all_ITS[rowSums(ASV_all_ITS)>ncol(ASV_all_ITS),]

ASV_filtered_16S<-ASV_filtered_16S[,order(colnames(ASV_filtered_16S))]
ASV_filtered_ITS<-ASV_filtered_ITS[,order(colnames(ASV_filtered_ITS))]

colSums(ASV_filtered_16S)
colSums(ASV_filtered_ITS)

ASV_filtered_16S_export<- ASV_filtered_16S %>% 
  tibble::rownames_to_column(var = "OTU")
ASV_filtered_ITS_export<- ASV_filtered_ITS %>% 
  tibble::rownames_to_column(var = "OTU")



######## export filtered ASVs #########
write.table(ASV_filtered_16S_export,"ASV_filtered_16S.tsv",
            sep = '\t',quote = FALSE, row.names = FALSE)
write.table(ASV_filtered_ITS_export,"ASV_filtered_ITS.tsv",
            sep = '\t',quote = FALSE, row.names = FALSE)




######## export filtered fasta #########
## 16S 
all_ASV_16S <- readDNAStringSet("16S_otus.filtered.fa")

filtered_ids_16S <- rownames(ASV_filtered_16S)

filtered_asvs_16S <- all_ASV_16S[names(all_ASV_16S) %in% filtered_ids_16S]

writeXStringSet(filtered_asvs_16S, "filtered_ASVs_16S.fa")




## ITS 
all_ASV_ITS <- readDNAStringSet("ITS_otus.fa")

filtered_ids_ITS <- rownames(ASV_filtered_ITS)

filtered_asvs_ITS <- all_ASV_ITS[names(all_ASV_ITS) %in% filtered_ids_ITS]

writeXStringSet(filtered_asvs_ITS, "filtered_ASVs_ITS.fa")









library(dplyr)
library(tidyr)
library(stringr)

ASV_all_ITS_With_Tax<-read.table("ASV_Abundance_ITS_split.txt",sep = "\t", header = TRUE,row.names = 1)

ITS_for_FUNGuild <-ASV_all_ITS_With_Tax%>% 
  mutate(taxonomy = paste0(
  "k__", Kingdom,
  ";p__", Phylum,
  ";c__", Class,
  ";o__", Order,
  ";f__", Family,
  ";g__", Genus,
  ";s__", Species
)) %>% 
  select(-Kingdom, -Phylum, -Class,-Order,-Family,-Genus,-Species) %>% 
  tibble::rownames_to_column(var = "OTU")


head(ITS_for_FUNGuild$taxonomy)

write.table(ITS_for_FUNGuild, "ITS_for_FUNGuild.tsv",
            sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)



