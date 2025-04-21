#Plot_Survival_on_stone_GO_term
rm(list = ls())

#input GO_term related with Survival on stone
Survival_GO<-read.table("Survival_in_stone_GO_term.txt",sep = "\t",header = TRUE,row.names = 1)

#input all Unigenes
all_unigenes<-read.table("Unigenes_GO.txt", sep = "\t",header = TRUE,comment.char = "",quote = "")


#input Unigene abundance
all_abundance<-read.table("4_Unigenes_count.txt",header = TRUE)


library(tidyverse)
Survival_unigenes<-filter(all_unigenes,all_unigenes$GO_ID %in% Survival_GO$ID)

Survival_abundace<- filter(all_abundance,Unigene_ID %in% Survival_unigenes$Query) 

Survival_unigenes$CK_abundance<-Survival_abundace[match(Survival_unigenes$Query, Survival_abundace$Unigene_ID),"CK"]
Survival_unigenes$Black_abundance<-Survival_abundace[match(Survival_unigenes$Query, Survival_abundace$Unigene_ID),"Black"]
Survival_unigenes$Green_abundance<-Survival_abundace[match(Survival_unigenes$Query, Survival_abundace$Unigene_ID),"Green"]
Survival_unigenes$White_abundance<-Survival_abundace[match(Survival_unigenes$Query, Survival_abundace$Unigene_ID),"White"]

Survival_unigenes %>% group_by(GO_ID) %>% 
  mutate(CK_sum = sum(CK_abundance)) %>%
  mutate(Black_sum = sum(Black_abundance)) %>%
  mutate(Green_sum = sum(Green_abundance)) %>%
  mutate(White_sum = sum(White_abundance)) -> Survival_GO_abundance

Survival_GO_abundance%>%select(GO_ID,GO_Term,GO_Function,GO_Level,
                               CK_sum,Black_sum,Green_sum,White_sum)%>%unique() ->Survival_GO_abundance

Survival_GO_abundance$Functions<-Survival_GO[match(Survival_GO_abundance$GO_ID,Survival_GO$ID),"Functions"]

Survival_GO_abundance<-Survival_GO_abundance[,c(9,1:8)]


write.csv(Survival_GO_abundance, "Survival_GO_abundance.csv",row.names = FALSE,quote = FALSE)  





Survival_GO_abundance_plot<-cbind(Survival_GO_abundance[,c(1,3)],log10(Survival_GO_abundance[,6:9]))


library(reshape2)
Survival_GO_abundance_plot<-melt(Survival_GO_abundance_plot)


Survival_GO_abundance_plot<-arrange(Survival_GO_abundance_plot,Functions) 

Survival_GO_abundance_plot$GO_Term <- factor(Survival_GO_abundance_plot$GO_Term,levels = unique(Survival_GO_abundance_plot$GO_Term))

Survival_GO_abundance_plot <- Survival_GO_abundance_plot %>%
  mutate(Group = str_extract(variable, "^[^_]+"))


library(ggsci)
library(ggplot2)
p1<-ggplot(Survival_GO_abundance_plot,aes(GO_Term,value,color=Functions, fill=Functions))+
  stat_summary(fun = mean, geom="bar",alpha=0.05,
               width = 0.65,
               show.legend=T)+
 
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1/2),
               geom="errorbar", width=0.3,show.legend=T)+

  geom_jitter(aes(shape=Group),alpha=0.8,size=2.5,width = 0.27, show.legend=T)+
  scale_color_locuszoom()+
  scale_fill_locuszoom()+
  theme_classic() +
  theme(legend.title = element_text(face = "bold",size = 12))+
  theme(axis.text.x=element_text(colour="black",size=10,face = "bold",angle = 45,hjust = 1,vjust = 0.97)) + 
  theme(axis.text.y=element_text(colour = "black",size=10)) + 
  theme(axis.title.y=element_text(colour = "black",size=12)) + 
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.title.y = element_text(face = "bold")) +
  xlab(NULL)+
  geom_vline(aes(xintercept =2.5),linetype="dashed", linewidth=0.6, colour="lightgray")+
  geom_vline(aes(xintercept =7.5),linetype="dashed", linewidth=0.6, colour="lightgray")+
  ylab("Log10 abundance value of reads")
p1

#save image 6*12



