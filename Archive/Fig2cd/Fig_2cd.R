

rm(list=ls())

library(metacoder)
library(tidyverse)

ASV_all<-read.table("Rarefied_ASV_Abundance_16S_Tax.txt",sep = "\t", header = TRUE,row.names = 1)



infor<-data.frame(Samples=c(rep(NA,27)),Groups=c(rep(NA,27)))
infor$Samples<-colnames(ASV_all)[1:27]
infor$Groups<-ifelse(substr(infor$Samples,1,3)=="Bla","Black",
                     ifelse(substr(infor$Samples,1,3)=="Whi","White",
                     ifelse(substr(infor$Samples,1,3)=="Gre","Green","CK")))


ASV_all$ASV_Sum<-apply(ASV_all[,1:27],1,sum)
ASV_all<-arrange(ASV_all,desc(ASV_Sum))
ASV_filtered<-ASV_all[1:(dim(ASV_all)[1]*0.05),]



ASV_filtered<-separate(ASV_filtered,
         Taxonomy,
         into = c("Taxonomy","Species"),
         sep = ",s:")



obj <- parse_tax_data(ASV_filtered,
                      class_cols = "Taxonomy",
                      class_sep = ",",
                      class_regex = "^(.+):(.+)$",
                      class_key = c(tax_rank = "info",
                                    tax_name = "taxon_name"))
print(obj)





obj$data$tax_abund <- calc_taxon_abund(obj, "tax_data",
                                       cols = infor$Samples)

obj$data$tax_abund$total<-rowSums(obj$data$tax_abund[,-1])
  

obj$data$tax_occ <- calc_n_samples(obj, "tax_abund",
                                     groups = infor$Groups,
                                     cols = infor$Samples)


library(RColorBrewer)
display.brewer.all()
display.brewer.pal(9,"YlOrRd")
brewer.pal(9,"YlOrRd")


set.seed(6666)

heat_tree(obj,
          node_label = obj$taxon_names(),
          node_size = obj$n_obs(),
          node_color = log10(obj$data$tax_abund$total),
          node_size_axis_label = "ASVs count",
          node_color_axis_label = "log10 Total reads",
          node_color_range = c(brewer.pal(9,"YlOrRd")[1:7]),
          layout = "reingold-tilford",
          initial_layout = "fruchterman-reingold")

#export pdf 6*6 inches












rm(list=ls())

library(metacoder)
library(tidyverse)

ASV_all<-read.table("Rarefied_ASV_Abundance_ITS_Tax.txt",sep = "\t", header = TRUE,row.names = 1)


ASV_all<-na.omit(ASV_all)


infor<-data.frame(Samples=c(rep(NA,27)),Groups=c(rep(NA,27)))
infor$Samples<-colnames(ASV_all)[1:27]
infor$Groups<-ifelse(substr(infor$Samples,1,3)=="Bla","Black",
                     ifelse(substr(infor$Samples,1,3)=="Whi","White",
                            ifelse(substr(infor$Samples,1,3)=="Gre","Green","CK")))

ASV_all$ASV_Sum<-apply(ASV_all[,1:27],1,sum)
ASV_all<-arrange(ASV_all,desc(ASV_Sum))
ASV_filtered<-ASV_all[1:(dim(ASV_all)[1]*0.3),]




ASV_filtered<-separate(ASV_filtered,
                       Taxonomy,
                       into = c("Taxonomy","Species"),
                       sep = ",s:")


obj <- parse_tax_data(ASV_filtered,
                      class_cols = "Taxonomy",
                      class_sep = ",",
                      class_regex = "^(.+):(.+)$",
                      class_key = c(tax_rank = "info",
                                    tax_name = "taxon_name"))
print(obj)




obj$data$tax_abund <- calc_taxon_abund(obj, "tax_data",
                                       cols = infor$Samples)
obj$data$tax_abund$total<-rowSums(obj$data$tax_abund[,-1])

obj$data$tax_occ <- calc_n_samples(obj, "tax_abund",
                                   groups = infor$Groups,
                                   cols = infor$Samples)


library(RColorBrewer)
display.brewer.all()
display.brewer.pal(9,"YlOrRd")
display.brewer.pal(9,"YlGnBu")
brewer.pal(9,"YlOrRd")

set.seed(666)
heat_tree(obj,
          node_label = obj$taxon_names(),
          node_size = obj$n_obs(),
          node_color = log10(obj$data$tax_abund$total),
          node_size_axis_label = "ASVs count",
          node_color_axis_label = "log10 Total reads",
          node_color_range = c(brewer.pal(9,"YlGnBu")[1:7]),
          layout = "reingold-tilford",
          initial_layout = "fruchterman-reingold")



#export pdf 6*6 inches


