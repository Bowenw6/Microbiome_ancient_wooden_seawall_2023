rm(list=ls())



library(ggplot2)
library(ggpubr)
library(ggsci)

wbwPalette <- c("#EEA236FF","#D43F3AFF","#5CB85CFF", "#46B8DAFF")

beta_nti_CK <- read.table("beta_nti_CK_melted.txt", sep = '\t', row.name = 1, header = TRUE)
beta_nti_Bla <- read.table("beta_nti_Bla_melted.txt", sep = '\t', row.name = 1, header = TRUE)
beta_nti_Gre <- read.table("beta_nti_Gre_melted.txt", sep = '\t', row.name = 1, header = TRUE)
beta_nti_Whi <- read.table("beta_nti_Whi_melted.txt", sep = '\t', row.name = 1, header = TRUE)

raup_crick_CK<- read.table("raup_crick_CK.txt", sep = '\t', row.name = 1, header = TRUE)
raup_crick_Bla<- read.table("raup_crick_Bla.txt", sep = '\t', row.name = 1, header = TRUE)
raup_crick_Gre<- read.table("raup_crick_Gre.txt", sep = '\t', row.name = 1, header = TRUE)
raup_crick_Whi<- read.table("raup_crick_Whi.txt", sep = '\t', row.name = 1, header = TRUE)



beta_nti_CK$Group<-rep(c("CK"),times=dim(beta_nti_CK)[1])
beta_nti_Bla$Group<-rep(c("Bla"),times=dim(beta_nti_Bla)[1])
beta_nti_Gre$Group<-rep(c("Gre"),times=dim(beta_nti_Gre)[1])
beta_nti_Whi$Group<-rep(c("Whi"),times=dim(beta_nti_Whi)[1])

sum_beta_nti<- rbind(beta_nti_CK[,3:4],beta_nti_Bla[,3:4],beta_nti_Gre[,3:4],beta_nti_Whi[,3:4])
colnames(sum_beta_nti)<-c("betaNTI", "Group")

raup_crick_CK$Group<-rep(c("CK"),times=dim(raup_crick_CK)[1])
raup_crick_Bla$Group<-rep(c("Bla"),times=dim(raup_crick_Bla)[1])
raup_crick_Gre$Group<-rep(c("Gre"),times=dim(raup_crick_Gre)[1])
raup_crick_Whi$Group<-rep(c("Whi"),times=dim(raup_crick_Whi)[1])

sum_raup_crick<- rbind(raup_crick_CK[,3:4],raup_crick_Bla[,3:4],raup_crick_Gre[,3:4],raup_crick_Whi[,3:4])
colnames(sum_raup_crick)<-c("RCbray", "Group")
compaired <- list(c('CK', 'Bla'),c('CK','Whi'),
                  c('CK', 'Gre'))

sum_beta_nti$Group <- factor(sum_beta_nti$Group,levels=c("CK","Bla","Gre","Whi"))
sum_raup_crick$Group <- factor(sum_raup_crick$Group,levels=c("CK","Bla","Gre","Whi"))

#################



p1<- ggplot(sum_beta_nti, aes(Group, betaNTI)) + 
  stat_boxplot(geom = 'errorbar', width = 0.4, size=0.3) + 
  geom_boxplot(fill= wbwPalette, width = 0.6, outlier.colour = NA, size=0.8, alpha=0.3) + 
  stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = -0.1) + 
  theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
  geom_jitter(aes(fill = Group), shape=21, width = 0.2, size = 1.2, alpha = 0.8)+ 
  scale_fill_manual(values=wbwPalette)+
  ylab("betaNTI")+xlab(NULL)+ 
  ylim(-2.5,10)+ 
  geom_hline(aes(yintercept=2), linetype='dashed',colour="red",alpha = 0.8) +
  geom_hline(aes(yintercept=-2), linetype='dashed',colour="red",alpha = 0.8) +
  geom_text(aes(x=1.6,y=-2.5,label="Fungi"))+
  theme(legend.position="none") 

p1

ggsave("BetaNTI_ITS_2.22x4.pdf", plot = p1, width = 2.22, height = 4)


p2<- ggplot(sum_raup_crick, aes(Group, RCbray)) + 
  stat_boxplot(geom = 'errorbar', width = 0.4, size=0.3) + 
  geom_boxplot(fill= wbwPalette, width = 0.6, outlier.colour = NA, size=0.8, alpha=0.3) + 
  stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = -0.2) + 
  theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
  geom_jitter(aes(fill = Group), shape=21, width = 0.2, size = 1.2, alpha = 0.8)+ 
  scale_fill_manual(values=wbwPalette)+
  ylab("RCbray")+xlab(NULL)+ 
  ylim(-2.1,2.1)+ 
  geom_hline(aes(yintercept=0.95), linetype='dashed',colour="red",alpha = 0.8) +
  geom_hline(aes(yintercept=-0.95), linetype='dashed',colour="red",alpha = 0.8) +
  theme(legend.position="none") 

p2



library(ggpubr)
p3<-ggarrange(p1, p2, ncol = 2, nrow = 1)
p3 

#save pdf,4inchs*6inchs
#or save pdf,4inchs*6inchs




