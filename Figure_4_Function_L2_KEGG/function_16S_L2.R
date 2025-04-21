rm(list=ls())
library(RColorBrewer)
display.brewer.all()
brewer.pal(9,'Set1')
display.brewer.pal(9,'Set1')
brewer.pal(11,'RdYlGn')
display.brewer.pal(11,'RdYlGn')


all_16S <- read.table('16S_KEGG.pathway.L2.abundance.txt', sep = '\t',header = T)
all_16S<-all_16S[,-9]
colnames(all_16S)<-c("L2","MZ_1","TR_1","MC_1","MZ_2","TR_2","MC_2","MZ_3")
all_16S<-all_16S[,c(1,2,5,8,4,7,3,6)]
only_number_16S<-all_16S[,-1]
only_number_16S_tmp<-t(only_number_16S)  
only_number_16S<-only_number_16S_tmp/rowSums(only_number_16S_tmp)
only_number_16S<-t(only_number_16S)
all_16S<-cbind(all_16S[,1],only_number_16S)
colnames(only_number_16S)
colnames(all_16S)<-c("L2","MZ_1","MZ_2","MZ_3","MC_1","MC_2","TR_1","TR_2")
all_16S<-as.data.frame(all_16S)
row.names(all_16S)<-all_16S$L2
all_16S_tmp<-all_16S[,-1]


all_16S_tmp<-as.data.frame(lapply(all_16S_tmp,as.numeric))
all_16S_tmp$sum<-apply(all_16S_tmp,1,sum)
all_16S<-cbind(all_16S$L2,all_16S_tmp)
colnames(all_16S)[1]<-c("L2")
all_16S<-all_16S[rev(order(all_16S$sum)),]
all_16S<-all_16S[,-9]

all_16S<-all_16S[1:20,]

#plot
library(ggplot2)
library(reshape2)
all_16S_long<-melt(all_16S, id= "L2")
colnames(all_16S_long)<-c("L2","Sample","Relative_Abundance")




p<-ggplot(all_16S_long,
          aes(x=Sample,y=reorder(L2,Relative_Abundance)))+ 
  geom_point(aes(colour=Sample,size=Relative_Abundance),shape=16,alpha=0.85) +
  scale_color_manual(values = c("#A50026","#D73027","#F46D43",
                                "#1A9850", "#66BD63",
                                "#FDAE61", "#FEE08B")) +theme_bw()+
  ylab("L2")+
  xlab("Relative abundance")
p

#save image 6*8 inches


