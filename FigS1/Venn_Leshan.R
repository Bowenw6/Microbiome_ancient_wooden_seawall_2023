
rm(list=ls())
library(RColorBrewer)
display.brewer.all()
brewer.pal(9,'Set1')
display.brewer.pal(9,'Set1')
display.brewer.pal(8,'Set2')
brewer.pal(8,'Set2')

all_16S <- read.table('ASV_Abundance_16S_split.txt', sep = '\t', row.names = 1,header = T)
all_16S<-all_16S[,1:27]

all_16S<-all_16S[rowSums(all_16S)>0,]


CK_16S<-all_16S[,c(7:12)]
Bla_16S<-all_16S[,c(1:6)]
Gre_16S<-all_16S[,c(19:24,15,18)]
Whi_16S<-all_16S[,c(25:27,13,14,16,17)]
CK_16S_sum<-apply(CK_16S,1,sum)
Bla_16S_sum<-apply(Bla_16S,1,sum)
Gre_16S_sum<-apply(Gre_16S,1,sum)
Whi_16S_sum<-apply(Whi_16S,1,sum)


CK_16S<-cbind(CK_16S,CK_16S_sum)
Bla_16S<-cbind(Bla_16S,Bla_16S_sum)
Gre_16S<-cbind(Gre_16S,Gre_16S_sum)
Whi_16S<-cbind(Whi_16S,Whi_16S_sum)


Leshan_16S_venn<-cbind(CK_16S[,ncol(CK_16S)],
                       Bla_16S[,ncol(Bla_16S)],
                       Gre_16S[,ncol(Gre_16S)],
                       Whi_16S[,ncol(Whi_16S)])

colnames(Leshan_16S_venn)<-c("CK","Black","Green","White")
rownames(Leshan_16S_venn)<-rownames(all_16S)




otus <-rownames(Leshan_16S_venn)
head(otus)
sets =colnames(Leshan_16S_venn)
sets

option ="CK-Black-Green-White"
option

newgroup =unlist(strsplit(option,"-"))
newgroup


x =list()

for(i in newgroup){
  x[[i]]=otus[as.numeric(Leshan_16S_venn[,i])>0]
}
str(x)

library(VennDiagram)
venn.diagram(x,col=c("#A65628", "black","#4DAF4A","#B3B3B3"),fill=c("#A65628", "black","#4DAF4A","#B3B3B3"),lwd=3,filename="Leshan_16S_Venn.tif",cex=0.6,cat.cex=0.7,width=1000,height=1000,alpha=0.4)
















rm(list=ls())

all_ITS <- read.table('ASV_Abundance_ITS_split.txt', sep = '\t', row.names = 1,header = T)
all_ITS<-all_ITS[,1:27]

all_ITS<-all_ITS[rowSums(all_ITS)>0,]

CK_ITS<-all_ITS[,c(7:12)]
Bla_ITS<-all_ITS[,c(1:6)]
Gre_ITS<-all_ITS[,c(19:24,15,18)]
Whi_ITS<-all_ITS[,c(25:27,13,14,16,17)]
CK_ITS_sum<-apply(CK_ITS,1,sum)
Bla_ITS_sum<-apply(Bla_ITS,1,sum)
Gre_ITS_sum<-apply(Gre_ITS,1,sum)
Whi_ITS_sum<-apply(Whi_ITS,1,sum)


CK_ITS<-cbind(CK_ITS,CK_ITS_sum)
Bla_ITS<-cbind(Bla_ITS,Bla_ITS_sum)
Gre_ITS<-cbind(Gre_ITS,Gre_ITS_sum)
Whi_ITS<-cbind(Whi_ITS,Whi_ITS_sum)


Leshan_ITS_venn<-cbind(CK_ITS[,ncol(CK_ITS)],
                       Bla_ITS[,ncol(Bla_ITS)],
                       Gre_ITS[,ncol(Gre_ITS)],
                       Whi_ITS[,ncol(Whi_ITS)])

colnames(Leshan_ITS_venn)<-c("CK","Black","Green","White")
rownames(Leshan_ITS_venn)<-rownames(all_ITS)


otus <-rownames(Leshan_ITS_venn)
head(otus)

sets =colnames(Leshan_ITS_venn)
sets

option ="CK-Black-Green-White"
option


newgroup =unlist(strsplit(option,"-"))
newgroup


x =list()

for(i in newgroup){
  x[[i]]=otus[as.numeric(Leshan_ITS_venn[,i])>0]
}
str(x)

library(VennDiagram)
venn.diagram(x,col=c("#A65628", "black","#4DAF4A","#B3B3B3"),fill=c("#A65628", "black","#4DAF4A","#B3B3B3"),lwd=3,filename="Leshan_ITS_Venn.tif",cex=0.6,cat.cex=0.7,width=1000,height=1000,alpha=0.4)


