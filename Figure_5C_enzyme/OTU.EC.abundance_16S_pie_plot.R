rm(list=ls())
library(RColorBrewer)
display.brewer.all()
display.brewer.pal(9,'Blues')


OTU_EC <- read.table('OTU.EC.abundance_16S.txt', sep = '\t', row.names = 1, header = T, quote = "")#加了quote = ""就可以读取大表格

OTU_EC_pie<-OTU_EC[,which(colnames(OTU_EC)%in%c("EC.3.2.1.4","EC.3.2.1.21","EC.3.2.1.74",
                                                            "EC.3.2.1.91","EC.2.4.1.20","EC.3.2.1.86"))]

OTU_EC_pie$sum<-apply(OTU_EC_pie,1,sum)

OTU_EC_pie<-OTU_EC_pie[which(OTU_EC_pie$sum!=0),]




CT_Taxonomy <- read.table('feature_table.taxonomy.even_phylum_class.txt', sep = '\t', header = TRUE, quote = "")
CT_Taxonomy$Phylum<-gsub("p__","",CT_Taxonomy$Phylum)
CT_Taxonomy$Class<-gsub("c__","",CT_Taxonomy$Class)

library(dplyr)
## define a helper function 
empty_as_unknown <- function(x){
  if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
  ifelse(as.character(x)!="", x, "unknown")
}

## transform all columns
CT_Taxonomy<-CT_Taxonomy %>% mutate_each(funs(empty_as_unknown))

CT_Taxonomy_pie<-CT_Taxonomy[which(CT_Taxonomy$X%in%row.names(OTU_EC_pie)),]
OTU_EC_pie_plot<-OTU_EC_pie[which(row.names(OTU_EC_pie)%in%CT_Taxonomy_pie$X),]

CT_Taxonomy_pie<-CT_Taxonomy_pie[order(CT_Taxonomy_pie$X),]
OTU_EC_pie_plot<-OTU_EC_pie_plot[order(row.names(OTU_EC_pie_plot)),]

CT_Taxonomy_pie$X==row.names(OTU_EC_pie_plot)

OTU_EC_pie_plot<-OTU_EC_pie_plot[,-6]
CT_Taxonomy_pie<-CT_Taxonomy_pie[,-1]
EC_pie_plot<-cbind(OTU_EC_pie_plot,CT_Taxonomy_pie)




library(dplyr)
EC24120_plot<-EC_pie_plot[,c(1,6)]
EC24120_plot%>% group_by(Phylum) %>% summarise_each(funs(sum))->EC24120_plot_summary

EC32121_plot<-EC_pie_plot[,c(2,6)]
EC32121_plot%>% group_by(Phylum) %>% summarise_each(funs(sum))->EC32121_plot_summary

EC3214_plot<-EC_pie_plot[,c(3,6)]
EC3214_plot%>% group_by(Phylum) %>% summarise_each(funs(sum))->EC3214_plot_summary

EC32186_plot<-EC_pie_plot[,c(4,6)]
EC32186_plot%>% group_by(Phylum) %>% summarise_each(funs(sum))->EC32186_plot_summary

EC32191_plot<-EC_pie_plot[,c(5,6)]
EC32191_plot%>% group_by(Phylum) %>% summarise_each(funs(sum))->EC32191_plot_summary


plot_all<-cbind(EC24120_plot_summary,EC32121_plot_summary,
                EC3214_plot_summary,EC32186_plot_summary,EC32191_plot_summary)
plot_all<-plot_all[,-c(3,5,7,9)]
plot_all_numeric<-plot_all[,2:6]
plot_all$sum<-apply(plot_all_numeric,1,sum)
plot_all_top<-plot_all[which(plot_all$sum>20),]
plot_all_bottom<-plot_all[which(plot_all$sum<20),]

plot_all_bottom<-plot_all_bottom[,-1]
others<-apply(plot_all_bottom,2,sum)
others<-c("others",others)
Plot_all_21<-rbind(plot_all_top,others)
Plot_all_21$sum<-as.numeric(Plot_all_21$sum)
Plot_all_21<-Plot_all_21[order(Plot_all_21$sum),]
Plot_all_21<-Plot_all_21[,-7]
rownames(Plot_all_21)<-Plot_all_21$Phylum
Plot_all_21<-Plot_all_21[,-1]
Plot_all_21<-Plot_all_21[,c(2,3,4,5,1)]

library(reshape2)
Plot_all_21 %>%
  as_tibble(rownames = "from") %>%
  reshape2::melt(id.vars = "from") ->Plot_all_21_long 

Plot_all_21_long$value<-as.numeric(Plot_all_21_long$value)

library(circlize)



library(RColorBrewer)
display.brewer.all()
mycolor<-c(rev(brewer.pal(11,"PuOr")),rev(brewer.pal(11,"RdYlGn"))[2:11],rev(brewer.pal(5,"Accent")))

set.seed(123)
circos.par(gap.after = c(rep(2, length(unique(Plot_all_21_long [[1]]))-1), 12,
                         rep(2, length(unique(Plot_all_21_long [[2]]))-1), 12))
chordDiagram(Plot_all_21_long, grid.col = mycolor)
circos.clear()


#save image 8*8 inches





