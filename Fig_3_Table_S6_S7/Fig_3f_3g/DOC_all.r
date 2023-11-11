
###DOC calculation using the DOC r package
#library(devtools)
#install_github("Russel88/DOC")


rm(list=ls())
library(ggplot2)
library(vegan)
library(Hmisc) 


ASV_all<-read.table("ASV_Abundance_16S_split.txt",sep = "\t", header = TRUE,row.names = 1)
ASV_all<-ASV_all[,1:c(length(colnames(ASV_all))-7)]

ASV_filtered<-ASV_all[which(rowSums(ASV_all)>27),]
ASV_filtered<-ASV_filtered[,order(colnames(ASV_filtered))]



A=ASV_filtered 
C=A/rowSums(A)
ASV_filtered_1<-t(C)

ASV_filtered_1[ASV_filtered_1>0]<-1
ASV_filtered_1<-t(ASV_filtered_1)
ASV_filtered<-ASV_filtered[which(rowSums(ASV_filtered_1)>6),]
rm(A,C,ASV_filtered_1)



summary(colSums(ASV_filtered))



library(DOC)
set.seed(666)
DOC_16S_results <- DOC(ASV_filtered, R = 1000)

p_16S_value<-(sum(DOC_16S_results$LME$Slope >= 0) + 1) / (length(DOC_16S_results$LME$Slope) + 1)
Fns_16S_value<-mean(DOC_16S_results$FNS$Fns)
p1<-plot(DOC_16S_results) + theme_bw() +geom_text(aes(x = 0.25,y = 0.95,
                                                         label = paste("Prokaryotes", '\n',"Fns=",round(Fns_16S_value,4), '\n', "P<", round(p_16S_value, 4))))+ 
  theme(axis.text.x = element_text(colour = 'blue'),
        axis.title.x = element_text(colour = 'blue'),
        axis.text.y = element_text(colour='red'),
        axis.title.y = element_text(colour = 'red'))
p1
#save as pdf 4inches*4inches

rm(ASV_all,ASV_filtered,DOC_16S_results)










#ITS

ASV_all<-read.table("ASV_Abundance_ITS_split.txt",sep = "\t", header = TRUE,row.names = 1)
ASV_all<-ASV_all[,1:c(length(colnames(ASV_all))-7)]


ASV_filtered<-ASV_all[which(rowSums(ASV_all)>27),]

ASV_filtered<-ASV_filtered[,order(colnames(ASV_filtered))]


A=ASV_filtered 
C=A/rowSums(A)
ASV_filtered_1<-t(C)

ASV_filtered_1[ASV_filtered_1>0]<-1

ASV_filtered_1<-t(ASV_filtered_1)
ASV_filtered<-ASV_filtered[which(rowSums(ASV_filtered_1)>6),]
rm(A,C,ASV_filtered_1)



summary(colSums(ASV_filtered))




library(DOC)
set.seed(666)
DOC_ITS_results <- DOC(ASV_filtered, R = 1000)

p_ITS_value<-(sum(DOC_ITS_results$LME$Slope >= 0) + 1) / (length(DOC_ITS_results$LME$Slope) + 1)
Fns_ITS_value<-mean(DOC_ITS_results$FNS$Fns)
p2<-plot(DOC_ITS_results) + theme_bw() +geom_text(aes(x = 0.25,y = 0.95,
                                                      label = paste("Eukaryotes", '\n', "Fns=",round(Fns_ITS_value,4), '\n', "P<", round(p_ITS_value, 4)))) +
  theme(axis.text.x = element_text(colour = 'blue'),
        axis.title.x = element_text(colour = 'blue'),
        axis.text.y = element_text(colour='red'),
        axis.title.y = element_text(colour = 'red'))
p2
#save as pdf 4inches*4inches















library(ggpubr)
p3<-ggarrange(p1, p2, labels = c("(c)", "(d)"), 
              widths = c(0.5,0.5),
              ncol = 2, nrow = 1)
p3
#save as pdf 4inches*8inches




