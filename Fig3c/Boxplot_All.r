rm(list=ls())



library(ggplot2)
library(ggpubr)
library(ggsci)
library(RColorBrewer)
display.brewer.all()
display.brewer.pal(12,"Set3")
brewer.pal(9,"Set1")


wbwPalette <- brewer.pal(9,"Set1")[c(5,4)]

beta_nti_All_16S <- read.table("beta_nti_All_melted_16S.txt", sep = '\t', row.name = 1, header = TRUE)
beta_nti_All_ITS<- read.table("beta_nti_All_melted_ITS.txt", sep = '\t', row.name = 1, header = TRUE)

raup_crick_All_16S<- read.table("raup_crick_All_16S.txt", sep = '\t', row.name = 1, header = TRUE)
raup_crick_All_ITS<- read.table("raup_crick_All_ITS.txt", sep = '\t', row.name = 1, header = TRUE)


beta_nti_All_16S$Group<-rep(c("Bacteria"),times=dim(beta_nti_All_16S)[1])
beta_nti_All_ITS$Group<-rep(c("Fungi"),times=dim(beta_nti_All_ITS)[1])


sum_beta_nti<- rbind(beta_nti_All_16S[,3:4],beta_nti_All_ITS[,3:4])
colnames(sum_beta_nti)<-c("betaNTI", "Group")

raup_crick_All_16S$Group<-rep(c("Bacteria"),times=dim(raup_crick_All_16S)[1])
raup_crick_All_ITS$Group<-rep(c("Fungi"),times=dim(raup_crick_All_ITS)[1])


sum_raup_crick<- rbind(raup_crick_All_16S[,3:4],raup_crick_All_ITS[,3:4])
colnames(sum_raup_crick)<-c("RCbray", "Group")
compaired <- list(c('Bacteria', 'Fungi'))

sum_beta_nti$Group <- factor(sum_beta_nti$Group,levels=c('Bacteria', 'Fungi'))
sum_raup_crick$Group <- factor(sum_raup_crick$Group,levels=c('Bacteria', 'Fungi'))

#################


p1<- ggplot(sum_beta_nti, aes(Group, betaNTI)) + 
  stat_boxplot(geom = 'errorbar', width = 0.4, size=0.3) +
  geom_boxplot(fill= wbwPalette, width = 0.6, outlier.colour = NA, size=0.8, alpha=0.3) + 
  stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = -0.1) + 
  theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
  geom_jitter(aes(fill = Group), shape=21, width = 0.2, size = 1.2, alpha = 0.8)+ 
  scale_fill_manual(values=wbwPalette)+
  ylab("betaNTI")+xlab(NULL)+ 
  ylim(-2.5,11)+ 
  geom_hline(aes(yintercept=2), linetype='dashed',colour="red",alpha = 0.8) +
  geom_hline(aes(yintercept=-2), linetype='dashed',colour="red",alpha = 0.8) +
  theme(legend.position="none") 

p1
ggsave("beta_nti_All_melted_2.32x4,4.pdf", width = 2.32, height = 4)



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




