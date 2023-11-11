#Plot_MGEs
rm(list = ls())

#input MGEs annotations
MGE_all<-read.table("Unigenes_MGEs.txt",sep = "\t",header = TRUE,comment.char = "",quote = "")

#input all Unigenes
all_unigenes<-read.table("4_Unigenes_count.txt", sep = "\t",header = TRUE,comment.char = "",quote = "")


library(tidyverse)
MGEs_unigenes<-filter(all_unigenes,all_unigenes$Unigene_ID %in% MGE_all$Query)

MGE_all$CK_abundance<-MGEs_unigenes[match(MGE_all$Query,MGEs_unigenes$Unigene_ID),"CK"]
MGE_all$Black_abundance<-MGEs_unigenes[match(MGE_all$Query,MGEs_unigenes$Unigene_ID),"Black"]
MGE_all$Green_abundance<-MGEs_unigenes[match(MGE_all$Query,MGEs_unigenes$Unigene_ID),"Green"]
MGE_all$White_abundance<-MGEs_unigenes[match(MGE_all$Query,MGEs_unigenes$Unigene_ID),"White"]


MGE_all %>% group_by(MGE_class) %>% 
  mutate(CK_sum = sum(CK_abundance)) %>%
  mutate(Black_sum = sum(Black_abundance)) %>%
  mutate(Green_sum = sum(Green_abundance)) %>%
  mutate(White_sum = sum(White_abundance)) -> MGE_all


MGE_all%>%select(MGE_class,CK_sum,Black_sum,Green_sum,White_sum)%>%unique() ->MGEs_gene_abundance

MGEs_gene_abundance<-arrange(MGEs_gene_abundance,MGE_class)


MGEs_annotation_tmp<-MGE_all[,c("MGE_type","MGE_class")] %>% unique()
rownames(MGEs_annotation_tmp)<-MGEs_annotation_tmp$MGE_class
MGEs_annotation_tmp<-as.data.frame(MGEs_annotation_tmp)

MGEs_annotation_tmp<-arrange(MGEs_annotation_tmp,MGE_class)
MGEs_annotation<-MGEs_annotation_tmp[,1]
MGEs_annotation<-as.data.frame(MGEs_annotation)
rownames(MGEs_annotation)<- MGEs_annotation_tmp$MGE_class

rm(MGEs_annotation_tmp)



MGEs_plot<-cbind(MGEs_gene_abundance,MGEs_annotation)

MGEs_plot<-MGEs_plot[,c(1,6,2:5)]


MGEs_plot$sum<-apply(MGEs_plot[,3:6], 1, sum) 


MGEs_plot<-arrange(MGEs_plot,MGEs_annotation,-sum)

write.csv(MGEs_plot, "MGEs.csv",row.names = FALSE,quote = FALSE)  








MGEs_data<-MGEs_plot[,3:6]
rownames(MGEs_data)<-MGEs_plot$MGE_class
colnames(MGEs_data)<-c("CK","Black","Green","White")

MGEs_data<-log10(MGEs_data+10)





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


bk = unique(c(seq(1,4, length=50)))


MGEs_annotation$MGEs_annotation<-factor(MGEs_annotation$MGEs_annotation,levels = c("IS91","ISBf10","ISCR",    
                                                                       "ISHa1152",  "ISPvsp2",   "ISSfl3", 
                                                                       "Tn916",   "integrase", "istA",   
                                                                       "istB",     "plasmid",   "tniA",    
                                                                       "tniB",    "transposase"))

unique(MGEs_annotation$MGEs_annotation)
ann_colors<-list(MGEs_annotation=c('IS91'='#8DA0CB',
                                   'ISBf10'='#FDC086',
                                   'ISCR'='#FFFFB3',
                                   'ISHa1152'='red',
                                   'ISPvsp2'='#FFD92F',
                                   'ISSfl3'='#386CB0',
                                   'Tn916'='#BEBADA',
                                   'integrase'='#E5C494',
                                   'istA'='#66C2A5',
                                   'istB'='#7FC97F',
                                   'plasmid'='#F0027F',
                                   'tniA'='#A6D854',
                                   'tniB'='#E78AC3',
                                   'transposase'='#80B1D3'))

                 

library(pheatmap)
p1<- pheatmap(MGEs_data,
              cluster_cols = F,
              breaks = bk,
              cluster_row = F,
              show_colnames     = T,
              show_rownames     = T,
              color = colorRampPalette(c("#F7FBFF", 
                                         "#6BAED6",
                                         "#08306B"))(50), 
              legend = T,
              border_color = 'NA',
              fontsize = 3,
              fontsize_col = 5,
              angle_col=45,
              annotation_row=MGEs_annotation, 
             annotation_colors = ann_colors,
              annotation_names_row = FALSE) 



#save image 3*3 inches



