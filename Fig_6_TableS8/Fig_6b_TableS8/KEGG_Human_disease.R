#Plot_HD
rm(list = ls())

#input KEGG annotations
KEGG_all<-read.table("Unigenes_KEGG.txt",sep = "\t",header = TRUE,comment.char = "",quote = "")

#filter the col KEGGLevel1 == "Human Diseases"
HD_all<-KEGG_all %>% filter(KEGGLevel1 == "Human Diseases")

#Keep the Query, Identity, E_value, GeneName...
HD_all<-HD_all %>% select("Query","Identity","E_value","GeneName",
                          "GeneDescription","KEGGLevel2")%>%unique()

#input all Unigenes
all_unigenes<-read.table("4_Unigenes_count.txt", sep = "\t",header = TRUE,comment.char = "",quote = "")



library(tidyverse)
HD_unigenes<-filter(all_unigenes,all_unigenes$Unigene_ID %in% HD_all$Query)

HD_all <- merge(HD_all, HD_unigenes, by.x = "Query", by.y = "Unigene_ID", all.x = TRUE)


HD_all %>% group_by(GeneName) %>% 
  mutate(CK_sum = sum(CK)) %>%
  mutate(Black_sum = sum(Black)) %>%
  mutate(Green_sum = sum(Green)) %>%
  mutate(White_sum = sum(White)) -> HD_all


HD_all%>%select(GeneName,GeneDescription,CK_sum,Black_sum,Green_sum,White_sum,KEGGLevel2)%>%unique() ->HD_gene_abundance


HD_gene_abundance$sum<-apply(HD_gene_abundance[,3:6], 1, sum) 
HD_gene_abundance<-arrange(HD_gene_abundance,KEGGLevel2,-sum)


write.table(HD_gene_abundance, "HD.txt",row.names = FALSE,quote = FALSE,sep = "\t")  







HD_plot<-HD_gene_abundance[,3:6]
rownames(HD_plot)<-HD_gene_abundance$GeneName
colnames(HD_plot)<-c("CK","Black","Green","White")

HD_plot<-log10(HD_plot+10)



HD_annotation<-HD_gene_abundance[,7]
HD_annotation<-as.data.frame(HD_annotation)
row.names(HD_annotation)<-HD_gene_abundance$GeneName







library(RColorBrewer)
display.brewer.all()
display.brewer.pal(9,"YlOrRd")
brewer.pal(9,"YlOrRd")

display.brewer.pal(9,"Blues")
brewer.pal(9,"Blues")

display.brewer.pal(7,"Set2")
brewer.pal(7,"Set2")
display.brewer.pal(5,"Set3")
brewer.pal(5,"Set3")

display.brewer.pal(7,"Accent")
brewer.pal(7,"Accent")



bk = unique(c(seq(1,4.5, length=50)))


unique(HD_annotation$KEGGLevel2)
ann_colors<-list(HD_annotation=c('Drug resistance: antimicrobial'='#7FC97F',
                                 'Infectious disease: bacterial'='#F0027F'))

                 

#绘制图片
library(pheatmap)
p1<- pheatmap(HD_plot,
              cluster_cols = F,
              breaks = bk,
              cluster_row = F,
              show_colnames     = T,
              show_rownames     = T,
              color = colorRampPalette(c("#FFFFCC", 
                                         "#FEB24C",
                                         "#FC4E2A"))(50), 
              legend = T,
              border_color = 'NA',
              fontsize = 3,
              fontsize_col = 5,
              angle_col=90,
              annotation_row=HD_annotation, 
             annotation_colors = ann_colors, 
              annotation_names_row = FALSE) 



#save image 5*4.5 inches



