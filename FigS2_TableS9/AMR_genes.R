#Plot_AMRs
rm(list = ls())

#input AMRs
AMR_all<-read.table("Unigenes_CARD.txt",sep = "\t",header = TRUE,comment.char = "",quote = "")

#input all Unigenes
all_unigenes<-read.table("4_Unigenes_count.txt", sep = "\t",header = TRUE,comment.char = "",quote = "")


library(tidyverse)
AMRs_unigenes<-filter(all_unigenes,all_unigenes$Unigene_ID %in% AMR_all$Query)

AMR_all$CK_abundance<-AMRs_unigenes[match(AMR_all$Query,AMRs_unigenes$Unigene_ID),"CK"]
AMR_all$Black_abundance<-AMRs_unigenes[match(AMR_all$Query,AMRs_unigenes$Unigene_ID),"Black"]
AMR_all$Green_abundance<-AMRs_unigenes[match(AMR_all$Query,AMRs_unigenes$Unigene_ID),"Green"]
AMR_all$White_abundance<-AMRs_unigenes[match(AMR_all$Query,AMRs_unigenes$Unigene_ID),"White"]


AMR_all %>% group_by(AMR_Gene_Family) %>% 
  mutate(CK_sum = sum(CK_abundance)) %>%
  mutate(Black_sum = sum(Black_abundance)) %>%
  mutate(Green_sum = sum(Green_abundance)) %>%
  mutate(White_sum = sum(White_abundance)) -> AMR_all


AMR_all%>%select(AMR_Gene_Family,CK_sum,Black_sum,Green_sum,White_sum)%>%unique() ->AMRs_gene_abundance

AMRs_gene_abundance<-arrange(AMRs_gene_abundance,AMR_Gene_Family)

AMRs_annotation<-AMR_all[match(AMRs_gene_abundance$AMR_Gene_Family,AMR_all$AMR_Gene_Family),c("AMR_Gene_Family","Resistance_Mechanism")]

AMRs_annotation<-arrange(AMRs_annotation,AMR_Gene_Family)

AMRs_annotation_tmp<-AMRs_annotation[,2]
AMRs_annotation_tmp<-as.data.frame(AMRs_annotation_tmp)
rownames(AMRs_annotation_tmp)<-AMRs_annotation$AMR_Gene_Family


AMRs_annotation<-AMRs_annotation_tmp


rm(AMRs_annotation_tmp)



AMRs_plot<-cbind(AMRs_gene_abundance,AMRs_annotation)

AMRs_plot<-AMRs_plot[,c(1,6,2:5)]


AMRs_plot$sum<-apply(AMRs_plot[,3:6], 1, sum) 


AMRs_plot<-arrange(AMRs_plot,Resistance_Mechanism,-sum)


write.csv(AMRs_plot, "AMRs.csv",row.names = FALSE,quote = FALSE)  









AMRs_data<-AMRs_plot[,3:6]
rownames(AMRs_data)<-AMRs_plot$AMR_Gene_Family
colnames(AMRs_data)<-c("CK","Black","Green","White")

AMRs_data<-log10(AMRs_data+10)





library(RColorBrewer)
display.brewer.all()
display.brewer.pal(9,"YlOrRd")
brewer.pal(9,"YlOrRd")

display.brewer.pal(7,"Set2")
brewer.pal(7,"Set2")
display.brewer.pal(5,"Set3")
brewer.pal(5,"Set3")



bk = unique(c(seq(1,5, length=50)))


unique(AMRs_annotation$Resistance_Mechanism)
ann_colors<-list(Resistance_Mechanism=c('N/A'='#8DD3C7',
                                        'antibiotic efflux'='#E78AC3',
                                        'antibiotic efflux;antibiotic target alteration'='#A6D854',
                                        'antibiotic efflux;reduced permeability to antibiotic'='#E5C494',
                                        'antibiotic inactivation'='#80B1D3',
                                        'antibiotic target alteration'='#66C2A5',
                                        'antibiotic target alteration;antibiotic target replacement'='#BEBADA',
                                       'antibiotic target protection'='#8DA0CB',
                                       'antibiotic target replacement'='#FFFFB3',
                                       'reduced permeability to antibiotic'='#FFD92F'))


library(pheatmap)
p1<- pheatmap(AMRs_data, 
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
              angle_col=45,
              annotation_row=AMRs_annotation, 
             annotation_colors = ann_colors, 
              annotation_names_row = FALSE) 



#save image 7*7 inches



