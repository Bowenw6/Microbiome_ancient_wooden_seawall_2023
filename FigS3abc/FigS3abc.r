###Supp file
#Calculate niche width with Shannon
rm(list=ls())


#Calculate All_16S_Bcom


ASV_all_16S_With_Tax<-read.table("ASV_Abundance_16S_split.txt",sep = "\t", header = TRUE,row.names = 1)
ASV_all_16S<-ASV_all_16S_With_Tax[,1:c(length(colnames(ASV_all_16S_With_Tax))-7)]


ASV_filtered_16S<- ASV_all_16S[rowSums(ASV_all_16S)>27,]

ASV_filtered_16S<-ASV_filtered_16S[,order(colnames(ASV_filtered_16S))]

colSums(ASV_filtered_16S)


library(spaa)
library(tidyverse)
t_ASV_filtered_16S <- t(ASV_filtered_16S)
niche_width_16S <- niche.width(t_ASV_filtered_16S, method = "shannon")

niche_width_16S<-t(niche_width_16S)
niche_width_16S_ASV<-niche_width_16S
rownames(niche_width_16S)<-NULL
niche_width_16S<-as.data.frame(niche_width_16S)
niche_width_16S$Group<-c("Prokaryotes")
colnames(niche_width_16S)[1]<-("Bcom")


write.csv(niche_width_16S,"niche_width_16S.csv",
          quote = FALSE,row.names = FALSE)







#Calculate All_ITS_Bcom


ASV_all_ITS_With_Tax<-read.table("ASV_Abundance_ITS_split.txt",sep = "\t", header = TRUE,row.names = 1)
ASV_all_ITS<-ASV_all_ITS_With_Tax[,1:c(length(colnames(ASV_all_ITS_With_Tax))-7)]


ASV_filtered_ITS<- ASV_all_ITS[rowSums(ASV_all_ITS)>27,]

ASV_filtered_ITS<-ASV_filtered_ITS[,order(colnames(ASV_filtered_ITS))]

colSums(ASV_filtered_ITS)


library(spaa)
library(tidyverse)
t_ASV_filtered_ITS <- t(ASV_filtered_ITS)
niche_width_ITS <- niche.width(t_ASV_filtered_ITS, method = "shannon")

niche_width_ITS<-t(niche_width_ITS)
niche_width_ITS_ASV<-niche_width_ITS
rownames(niche_width_ITS)<-NULL
niche_width_ITS<-as.data.frame(niche_width_ITS)
niche_width_ITS$Group<-c("Eukaryotes")
colnames(niche_width_ITS)[1]<-("Bcom")


write.csv(niche_width_ITS,"niche_width_ITS.csv",
          quote = FALSE,row.names = FALSE)









########## calculate 16S Bcom for each group #######


ASV_filtered_CK_16S<-ASV_filtered_16S[,7:12]
ASV_filtered_Bla_16S<-ASV_filtered_16S[,1:6]
ASV_filtered_Gre_16S<-ASV_filtered_16S[,13:20]
ASV_filtered_Whi_16S<-ASV_filtered_16S[,21:27]

ASV_filtered_CK_16S<- ASV_filtered_CK_16S[rowSums(ASV_filtered_CK_16S) > 6,]
ASV_filtered_Bla_16S<- ASV_filtered_Bla_16S[rowSums(ASV_filtered_Bla_16S) > 6,]
ASV_filtered_Gre_16S<- ASV_filtered_Gre_16S[rowSums(ASV_filtered_Gre_16S) > 8,]
ASV_filtered_Whi_16S<- ASV_filtered_Whi_16S[rowSums(ASV_filtered_Whi_16S) > 7,]

t_ASV_filtered_CK_16S <- t(ASV_filtered_CK_16S)
t_ASV_filtered_Bla_16S <- t(ASV_filtered_Bla_16S)
t_ASV_filtered_Gre_16S <- t(ASV_filtered_Gre_16S)
t_ASV_filtered_Whi_16S <- t(ASV_filtered_Whi_16S)






data_frames <- c("CK","Bla", "Gre", "Whi")

for (df_name in data_frames) {
  df_var <- paste0("t_ASV_filtered_", df_name, "_16S")
  
  niche_width_16S <- niche.width(get(df_var), method = "shannon")
  
  niche_width_16S<-t(niche_width_16S)
  niche_width_16S_ASV<-niche_width_16S
  rownames(niche_width_16S)<-NULL
  niche_width_16S<-as.data.frame(niche_width_16S)
  niche_width_16S$Group<-df_name
  colnames(niche_width_16S)[1]<-("Bcom")
  
  write.csv(niche_width_16S,
            paste0("niche_width_", df_name, "_16S.csv"),
            quote = FALSE, row.names = FALSE)
  
}













########## calculate ITS Bcom for each group #######


ASV_filtered_CK_ITS<-ASV_filtered_ITS[,7:12]
ASV_filtered_Bla_ITS<-ASV_filtered_ITS[,1:6]
ASV_filtered_Gre_ITS<-ASV_filtered_ITS[,13:20]
ASV_filtered_Whi_ITS<-ASV_filtered_ITS[,21:27]

ASV_filtered_CK_ITS<- ASV_filtered_CK_ITS[rowSums(ASV_filtered_CK_ITS) > 6,]
ASV_filtered_Bla_ITS<- ASV_filtered_Bla_ITS[rowSums(ASV_filtered_Bla_ITS) > 6,]
ASV_filtered_Gre_ITS<- ASV_filtered_Gre_ITS[rowSums(ASV_filtered_Gre_ITS) > 8,]
ASV_filtered_Whi_ITS<- ASV_filtered_Whi_ITS[rowSums(ASV_filtered_Whi_ITS) > 7,]

t_ASV_filtered_CK_ITS <- t(ASV_filtered_CK_ITS)
t_ASV_filtered_Bla_ITS <- t(ASV_filtered_Bla_ITS)
t_ASV_filtered_Gre_ITS <- t(ASV_filtered_Gre_ITS)
t_ASV_filtered_Whi_ITS <- t(ASV_filtered_Whi_ITS)






data_frames <- c("CK","Bla", "Gre", "Whi")

for (df_name in data_frames) {
  df_var <- paste0("t_ASV_filtered_", df_name, "_ITS")
  
  niche_width_ITS <- niche.width(get(df_var), method = "shannon")
  
  niche_width_ITS<-t(niche_width_ITS)
  niche_width_ITS_ASV<-niche_width_ITS
  rownames(niche_width_ITS)<-NULL
  niche_width_ITS<-as.data.frame(niche_width_ITS)
  niche_width_ITS$Group<-df_name
  colnames(niche_width_ITS)[1]<-("Bcom")
  
  write.csv(niche_width_ITS,
            paste0("niche_width_", df_name, "_ITS.csv"),
            quote = FALSE, row.names = FALSE)
  
}


























library(tidyverse)

niche_width_16S<-read.csv(file = "niche_width_16S.csv",
                                        #    quote = FALSE,
                                        header = TRUE,
                                        sep = ",")

niche_width_ITS<-read.csv(file = "niche_width_ITS.csv",
                                        #    quote = FALSE,
                                        header = TRUE,
                                        sep = ",")


niche_width_CK_16S<-read.csv(file = "niche_width_CK_16S.csv",
                                           #    quote = FALSE,
                                           header = TRUE,
                                           sep = ",")

niche_width_Bla_16S<-read.csv(file = "niche_width_Bla_16S.csv",
                                            #    quote = FALSE,
                                            header = TRUE,
                                            sep = ",")

niche_width_Gre_16S<-read.csv(file = "niche_width_Gre_16S.csv",
                                            #    quote = FALSE,
                                            header = TRUE,
                                            sep = ",")

niche_width_Whi_16S<-read.csv(file = "niche_width_Whi_16S.csv",
                                            #    quote = FALSE,
                                            header = TRUE,
                                            sep = ",")



niche_width_CK_ITS<-read.csv(file = "niche_width_CK_ITS.csv",
                                           #    quote = FALSE,
                                           header = TRUE,
                                           sep = ",")

niche_width_Bla_ITS<-read.csv(file = "niche_width_Bla_ITS.csv",
                                            #    quote = FALSE,
                                            header = TRUE,
                                            sep = ",")

niche_width_Gre_ITS<-read.csv(file = "niche_width_Gre_ITS.csv",
                                            #    quote = FALSE,
                                            header = TRUE,
                                            sep = ",")

niche_width_Whi_ITS<-read.csv(file = "niche_width_Whi_ITS.csv",
                                            #    quote = FALSE,
                                            header = TRUE,
                                            sep = ",")




######### merge files and plot ###


niche_width_All_Boxplot<-rbind(niche_width_16S,
                                 niche_width_ITS)




#### plot Bcom
library(RColorBrewer)
display.brewer.all()
display.brewer.pal(12,"Set3")
brewer.pal(9,"Set1")
wbwPalette <- brewer.pal(9,"Set1")[c(5,4)]
library(ggplot2)
library(ggpubr)
compaired <- list(c('Prokaryotes', 'Eukaryotes'))

niche_width_All_Boxplot$Group <- factor(niche_width_All_Boxplot$Group,levels=c('Prokaryotes', 'Eukaryotes'))


#plot niche_width_All_Boxplot


p11<- ggplot(niche_width_All_Boxplot, aes(Group, Bcom)) + 
  stat_boxplot(geom = 'errorbar', width = 0.4, size=0.2) + 
  geom_boxplot(fill= wbwPalette, width = 0.6, size=0.3, alpha=0.5,outlier.shape = NA) + 
  stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = -0.1) + 
  theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
  geom_violin(aes(fill = Group), width = 0.9, size = 0.5, alpha = 0.6)+ 
  scale_fill_manual(values=wbwPalette)+
  ylab("Bcom")+xlab(NULL)+ 
  ylim(0,4)+
  theme(legend.position="none")+ 
  ggtitle("Bcom-All")+
  theme(plot.title = element_text(hjust = 0.5))

p11











#Calculate Calculate CK_Black_Green_White_16S_Bcom

niche_width_4_16S_Boxplot<-rbind(niche_width_CK_16S,
                                   niche_width_Bla_16S,
                                   niche_width_Gre_16S,
                                   niche_width_Whi_16S)





library(RColorBrewer)
display.brewer.all()
display.brewer.pal(12,"Set3")
brewer.pal(9,"Set1")

library(ggplot2)
library(ggpubr)
compaired <- list(c('CK', 'Bla'),c('CK', 'Gre'),
                  c('CK', 'Whi'))
wbwPalette_4 <- c("#EEA236FF","#D43F3AFF","#5CB85CFF", "#46B8DAFF")

niche_width_4_16S_Boxplot$Group <- factor(niche_width_4_16S_Boxplot$Group,levels=c('CK','Bla','Gre','Whi'))


p12<- ggplot(niche_width_4_16S_Boxplot, aes(Group, Bcom)) + 
  stat_boxplot(geom = 'errorbar', width = 0.4, size=0.2) + 
  geom_boxplot(fill= wbwPalette_4, width = 0.6, size=0.3, alpha=0.5,outlier.shape = NA) + 
  stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = -0.1) + 
  theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
  geom_violin(aes(fill = Group), width = 0.9, size = 0.5, alpha = 0.6)+ 
  scale_fill_manual(values=wbwPalette_4)+
  ylab("Bcom")+xlab(NULL)+ 
  ylim(0,3)+ 
  theme(legend.position="none")+ 
  ggtitle("Bcom-Prokaryotes")+
  theme(plot.title = element_text(hjust = 0.5))
p12














#Calculate Calculate CK_Black_Green_White_ITS_Bcom

niche_width_4_ITS_Boxplot<-rbind(niche_width_CK_ITS,
                                   niche_width_Bla_ITS,
                                   niche_width_Gre_ITS,
                                   niche_width_Whi_ITS)





library(RColorBrewer)
display.brewer.all()
display.brewer.pal(12,"Set3")
brewer.pal(9,"Set1")

library(ggplot2)
library(ggpubr)
compaired <- list(c('CK', 'Bla'),c('CK', 'Gre'),
                  c('CK', 'Whi'))
wbwPalette_4 <- c("#EEA236FF","#D43F3AFF","#5CB85CFF", "#46B8DAFF")

niche_width_4_ITS_Boxplot$Group <- factor(niche_width_4_ITS_Boxplot$Group,levels=c('CK','Bla','Gre','Whi'))


p13<- ggplot(niche_width_4_ITS_Boxplot, aes(Group, Bcom)) + 
  stat_boxplot(geom = 'errorbar', width = 0.4, size=0.2) + 
  geom_boxplot(fill= wbwPalette_4, width = 0.6, size=0.3, alpha=0.5,outlier.shape = NA) + 
  stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = -0.1) +
  theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
  geom_violin(aes(fill = Group), width = 0.9, size = 0.5, alpha = 0.6)+ 
  scale_fill_manual(values=wbwPalette_4)+
  ylab("Bcom")+xlab(NULL)+ 
  ylim(0,3)+ 
  theme(legend.position="none")+ 
  ggtitle("Bcom-Eukaryotes")+
  theme(plot.title = element_text(hjust = 0.5))

p13






#### export 
library(ggpubr)
p1 <- ggarrange(ggarrange(p12, p13, p11,
                          ncol = 3,
                          widths = c(0.37,0.37,0.26)))
p1
# save image 3inches*11inches




