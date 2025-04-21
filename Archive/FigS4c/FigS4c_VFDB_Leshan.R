#Plot_VFDBs
rm(list = ls())

#input VFG_ID_annotation file 
all_VFG_ID_annotation<- read.csv("all_VFG_ID_annotation.csv",header = TRUE,sep = ',')
#input VFDBs annotations
VFDB_all<-read.table("Unigenes_VFDB.txt",sep = "\t",header = TRUE,comment.char = "",quote = "")
# repair VF_ID in VFDBs annotations 
VFDB_all$VF_id<-substr(VFDB_all$VF_id,1,9)




#input all Unigenes
all_unigenes<-read.table("4_Unigenes_count.txt", sep = "\t",header = TRUE,comment.char = "",quote = "")


library(tidyverse)
VFDBs_unigenes<-filter(all_unigenes,all_unigenes$Unigene_ID %in% VFDB_all$Query)

VFDB_all$CK_abundance<-VFDBs_unigenes[match(VFDB_all$Query,VFDBs_unigenes$Unigene_ID),"CK"]
VFDB_all$Black_abundance<-VFDBs_unigenes[match(VFDB_all$Query,VFDBs_unigenes$Unigene_ID),"Black"]
VFDB_all$Green_abundance<-VFDBs_unigenes[match(VFDB_all$Query,VFDBs_unigenes$Unigene_ID),"Green"]
VFDB_all$White_abundance<-VFDBs_unigenes[match(VFDB_all$Query,VFDBs_unigenes$Unigene_ID),"White"]


VFDB_all %>% group_by(VF_id) %>% 
  mutate(CK_sum = sum(CK_abundance)) %>%
  mutate(Black_sum = sum(Black_abundance)) %>%
  mutate(Green_sum = sum(Green_abundance)) %>%
  mutate(White_sum = sum(White_abundance)) -> VFDB_all


VFDB_all%>%select(VF_id,Gene,CK_sum,Black_sum,Green_sum,White_sum)%>%unique() ->VFDBs_gene_abundance

VFDBs_gene_abundance<-arrange(VFDBs_gene_abundance,VF_id)


library(tidyverse)
all_VF_Category<-filter(all_VFG_ID_annotation,all_VFG_ID_annotation$all_VFG_ID %in% VFDB_all$VF_id)


VFDBs_gene_abundance<-left_join(VFDBs_gene_abundance, all_VF_Category, by = c("VF_id" = "all_VFG_ID"))

VFDBs_gene_abundance$all_VF_category <- ifelse(is.na(VFDBs_gene_abundance$all_VF_category), "Others", VFDBs_gene_abundance$all_VF_category)


VFDBs_gene_abundance %>% group_by(all_VF_category) %>% 
  mutate(CK_Sum = sum(CK_sum)) %>%
  mutate(Black_Sum = sum(Black_sum)) %>%
  mutate(Green_Sum = sum(Green_sum)) %>%
  mutate(White_Sum = sum(White_sum)) -> VFDBs_gene_abundance

VFDBs_gene_abundance%>%select(all_VF_category,CK_Sum,Black_Sum,Green_Sum,White_Sum)%>%unique() ->VFDBs_gene_abundance

VFDBs_gene_abundance<-arrange(VFDBs_gene_abundance,all_VF_category)

VFDBs_gene_abundance_tmp<-VFDBs_gene_abundance[,-1]

VFDBs_gene_abundance_tmp$Sum<-apply(VFDBs_gene_abundance_tmp,1,sum)

row.names(VFDBs_gene_abundance_tmp)<-VFDBs_gene_abundance$all_VF_category
VFDBs_gene_abundance<-VFDBs_gene_abundance_tmp
rm(VFDBs_gene_abundance_tmp)

colnames(VFDBs_gene_abundance)<-c("CK","Black","Green","White","Sum")

VFDBs_gene_abundance<-arrange(VFDBs_gene_abundance,-Sum)


write.csv(VFDBs_gene_abundance, "VFDBs.csv",row.names = TRUE,quote = FALSE)  


VFDBs_gene_abundance_plot<-log10(VFDBs_gene_abundance)

VFDBs_gene_abundance_plot<-VFDBs_gene_abundance_plot[,-ncol(VFDBs_gene_abundance_plot)]
VFDBs_gene_abundance_plot$category<-rownames(VFDBs_gene_abundance)



library(reshape2)
VFDBs_gene_abundance_plot<-melt(VFDBs_gene_abundance_plot)


VFDBs_gene_abundance_plot$category <- factor(VFDBs_gene_abundance_plot$category,levels = unique(VFDBs_gene_abundance_plot$category))







library(ggsci)
#scale_color_locuszoom()
library(ggplot2)
library(ggbreak)
p1<-ggplot(VFDBs_gene_abundance_plot,aes(category,value,color=category, fill=category))+
  stat_summary(fun = mean, geom="bar",alpha=0.05,
               width = 0.65,
               show.legend=T)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1/2),
               geom="errorbar", width=0.3,show.legend=T)+
  geom_jitter(aes(shape = variable),alpha=0.8,size=2.5,width = 0.27, show.legend=T)+
  scale_shape_manual(values = c(CK=16,Black=15,Green=17,White=18))+
  scale_color_simpsons()+
  scale_fill_simpsons()+
  theme_classic() +
  theme(legend.title = element_text(face = "bold",size = 12))+
  theme(axis.text.x=element_text(colour="black",size=10,face = "bold",angle = 45,hjust = 1,vjust = 0.97)) + 
  theme(axis.text.y=element_text(colour = "black",size=10)) + 
  theme(axis.title.y=element_text(colour = "black",size=12)) + 
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.title.y = element_text(face = "bold")) +
  xlab(NULL)+
  ylab("Log10 abundance value of reads")+
  scale_y_break(breaks = c(1,3.4),scale=20)
  
p1

#save image 6*12 inches









