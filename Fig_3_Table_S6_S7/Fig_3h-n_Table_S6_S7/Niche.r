#niche breadth 

rm(list=ls())


ASV_all_16S_With_Tax<-read.table("ASV_Abundance_16S_split.txt",sep = "\t", header = TRUE,row.names = 1)
ASV_all_16S<-ASV_all_16S_With_Tax[,1:c(length(colnames(ASV_all_16S_With_Tax))-7)]


ASV_relative_16S<-ASV_all_16S/colSums(ASV_all_16S)

ASV_filtered_16S<- ASV_relative_16S[rowMeans(ASV_relative_16S) > 2*10^-5,]

ASV_filtered_16S<-ASV_filtered_16S[,order(colnames(ASV_filtered_16S))]



colSums(ASV_filtered_16S)


#Calculate All_16S_Bcom
library(spaa)
t_ASV_filtered_16S <- t(ASV_filtered_16S)
niche_breadth_16S <- niche.width(t_ASV_filtered_16S, method = "levins")

niche_breadth_16S<-t(niche_breadth_16S)
niche_breadth_16S_ASV<-niche_breadth_16S
rownames(niche_breadth_16S)<-NULL
niche_breadth_16S<-as.data.frame(niche_breadth_16S)
niche_breadth_16S$Group<-c("Prokaryotes")
colnames(niche_breadth_16S)[1]<-("Bcom")










ASV_all_ITS_With_Tax<-read.table("ASV_Abundance_ITS_split.txt",sep = "\t", header = TRUE,row.names = 1)
ASV_all_ITS<-ASV_all_ITS_With_Tax[,1:c(length(colnames(ASV_all_ITS_With_Tax))-7)]


ASV_relative_ITS<-ASV_all_ITS/colSums(ASV_all_ITS)

ASV_filtered_ITS<- ASV_relative_ITS[rowMeans(ASV_relative_ITS) > 2*10^-5,]



ASV_filtered_ITS<-ASV_filtered_ITS[,order(colnames(ASV_filtered_ITS))]




colSums(ASV_filtered_ITS)


#Calculate All_ITS_Bcom
library(spaa)
t_ASV_filtered_ITS <- t(ASV_filtered_ITS)
niche_breadth_ITS <- niche.width(t_ASV_filtered_ITS, method = "levins")

niche_breadth_ITS<-t(niche_breadth_ITS)
niche_breadth_ITS_ASV<-niche_breadth_ITS
rownames(niche_breadth_ITS)<-NULL
niche_breadth_ITS<-as.data.frame(niche_breadth_ITS)
niche_breadth_ITS$Group<-c("Eukaryotes")
colnames(niche_breadth_ITS)[1]<-("Bcom")




niche_breadth_all_boxplot<-rbind(niche_breadth_16S,niche_breadth_ITS)


library(RColorBrewer)
display.brewer.all()
display.brewer.pal(12,"Set3")
brewer.pal(9,"Set1")

wbwPalette <- brewer.pal(9,"Set1")[c(5,4)]
library(ggplot2)
library(ggpubr)

compaired <- list(c('Prokaryotes', 'Eukaryotes'))


niche_breadth_all_boxplot$Group <- factor(niche_breadth_all_boxplot$Group,levels=c('Prokaryotes', 'Eukaryotes'))


p1<- ggplot(niche_breadth_all_boxplot, aes(Group, Bcom)) + 
  stat_boxplot(geom = 'errorbar', width = 0.4, size=0.2) + 
  geom_boxplot(fill= wbwPalette, width = 0.6, size=0.3, alpha=0.5, outlier.shape = NA) + 
  stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = -0.1) + 
  theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
  geom_violin(aes(fill = Group), width = 0.9, size = 0.5, alpha = 0.6)+
  scale_fill_manual(values=wbwPalette)+
  ylab("Bcom")+xlab(NULL)+ 
  ylim(0,22)+ 
  geom_hline(aes(yintercept=3), linetype='dashed',colour="red",alpha = 0.8) +
  geom_hline(aes(yintercept=1.5), linetype='dashed',colour="red",alpha = 0.8) +

  theme(legend.position="none") 

p1
#save 6inches*3inches，4inches*3inches，4inches*3inches






















#####Calculate CK_Black_Green_White_16S_Bcom
library(spaa)

ASV_relative_CK_16S<-ASV_relative_16S[,7:12]
ASV_relative_Bla_16S<-ASV_relative_16S[,1:6]
ASV_relative_Gre_16S<-ASV_relative_16S[,13:20]
ASV_relative_Whi_16S<-ASV_relative_16S[,21:27]

ASV_filtered_CK_16S<- ASV_relative_CK_16S[rowMeans(ASV_relative_CK_16S) > 2*10^-5,]
ASV_filtered_Bla_16S<- ASV_relative_Bla_16S[rowMeans(ASV_relative_Bla_16S) > 2*10^-5,]
ASV_filtered_Gre_16S<- ASV_relative_Gre_16S[rowMeans(ASV_relative_Gre_16S) > 2*10^-5,]
ASV_filtered_Whi_16S<- ASV_relative_Whi_16S[rowMeans(ASV_relative_Whi_16S) > 2*10^-5,]


t_ASV_filtered_CK_16S <- t(ASV_filtered_CK_16S)
t_ASV_filtered_Bla_16S <- t(ASV_filtered_Bla_16S)
t_ASV_filtered_Gre_16S <- t(ASV_filtered_Gre_16S)
t_ASV_filtered_Whi_16S <- t(ASV_filtered_Whi_16S)

niche_breadth_CK_16S <- niche.width(t_ASV_filtered_CK_16S, method = "levins")
niche_breadth_Bla_16S <- niche.width(t_ASV_filtered_Bla_16S, method = "levins")
niche_breadth_Gre_16S <- niche.width(t_ASV_filtered_Gre_16S, method = "levins")
niche_breadth_Whi_16S <- niche.width(t_ASV_filtered_Whi_16S, method = "levins")




niche_breadth_CK_16S<-t(niche_breadth_CK_16S)
rownames(niche_breadth_CK_16S)<-NULL
niche_breadth_CK_16S<-as.data.frame(niche_breadth_CK_16S)
niche_breadth_CK_16S$Group<-c("CK")
colnames(niche_breadth_CK_16S)[1]<-("Bcom")


niche_breadth_Bla_16S<-t(niche_breadth_Bla_16S)
rownames(niche_breadth_Bla_16S)<-NULL
niche_breadth_Bla_16S<-as.data.frame(niche_breadth_Bla_16S)
niche_breadth_Bla_16S$Group<-c("Bla")
colnames(niche_breadth_Bla_16S)[1]<-("Bcom")


niche_breadth_Gre_16S<-t(niche_breadth_Gre_16S)
rownames(niche_breadth_Gre_16S)<-NULL
niche_breadth_Gre_16S<-as.data.frame(niche_breadth_Gre_16S)
niche_breadth_Gre_16S$Group<-c("Gre")
colnames(niche_breadth_Gre_16S)[1]<-("Bcom")


niche_breadth_Whi_16S<-t(niche_breadth_Whi_16S)
rownames(niche_breadth_Whi_16S)<-NULL
niche_breadth_Whi_16S<-as.data.frame(niche_breadth_Whi_16S)
niche_breadth_Whi_16S$Group<-c("Whi")
colnames(niche_breadth_Whi_16S)[1]<-("Bcom")




niche_breadth_4_16S_boxplot<-rbind(niche_breadth_CK_16S,niche_breadth_Bla_16S,
                                   niche_breadth_Gre_16S,niche_breadth_Whi_16S)



compaired <- list(c('CK', 'Bla'),c('CK', 'Gre'),
                  c('CK', 'Whi'))

wbwPalette_4 <- c("#EEA236FF","#D43F3AFF","#5CB85CFF", "#46B8DAFF")


niche_breadth_4_16S_boxplot$Group <- factor(niche_breadth_4_16S_boxplot$Group,levels=c('CK','Bla','Gre','Whi'))



p2<- ggplot(niche_breadth_4_16S_boxplot, aes(Group, Bcom)) + 
  stat_boxplot(geom = 'errorbar', width = 0.4, size=0.2) + 
  geom_boxplot(fill= wbwPalette_4, width = 0.6, size=0.3, alpha=0.5,outlier.shape = NA) + 
  stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = -0.1) + 
  theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
  geom_violin(aes(fill = Group), width = 0.9, size = 0.5, alpha = 0.6)+ 
  scale_fill_manual(values=wbwPalette_4)+
  ylab("Bcom")+xlab(NULL)+ 
  ylim(0,9)+ 

  geom_text(aes(x = 1.2,y = 0.2,label = c("Prokaryotes")))+
  theme(legend.position="none") 

p2
#save 3inches*4inches
















library(spaa)

ASV_relative_CK_ITS<-ASV_relative_ITS[,7:12]
ASV_relative_Bla_ITS<-ASV_relative_ITS[,1:6]
ASV_relative_Gre_ITS<-ASV_relative_ITS[,13:20]
ASV_relative_Whi_ITS<-ASV_relative_ITS[,21:27]

ASV_filtered_CK_ITS<- ASV_relative_CK_ITS[rowMeans(ASV_relative_CK_ITS) > 2*10^-5,]
ASV_filtered_Bla_ITS<- ASV_relative_Bla_ITS[rowMeans(ASV_relative_Bla_ITS) > 2*10^-5,]
ASV_filtered_Gre_ITS<- ASV_relative_Gre_ITS[rowMeans(ASV_relative_Gre_ITS) > 2*10^-5,]
ASV_filtered_Whi_ITS<- ASV_relative_Whi_ITS[rowMeans(ASV_relative_Whi_ITS) > 2*10^-5,]


t_ASV_filtered_CK_ITS <- t(ASV_filtered_CK_ITS)
t_ASV_filtered_Bla_ITS <- t(ASV_filtered_Bla_ITS)
t_ASV_filtered_Gre_ITS <- t(ASV_filtered_Gre_ITS)
t_ASV_filtered_Whi_ITS <- t(ASV_filtered_Whi_ITS)

niche_breadth_CK_ITS <- niche.width(t_ASV_filtered_CK_ITS, method = "levins")
niche_breadth_Bla_ITS <- niche.width(t_ASV_filtered_Bla_ITS, method = "levins")
niche_breadth_Gre_ITS <- niche.width(t_ASV_filtered_Gre_ITS, method = "levins")
niche_breadth_Whi_ITS <- niche.width(t_ASV_filtered_Whi_ITS, method = "levins")




niche_breadth_CK_ITS<-t(niche_breadth_CK_ITS)
rownames(niche_breadth_CK_ITS)<-NULL
niche_breadth_CK_ITS<-as.data.frame(niche_breadth_CK_ITS)
niche_breadth_CK_ITS$Group<-c("CK")
colnames(niche_breadth_CK_ITS)[1]<-("Bcom")


niche_breadth_Bla_ITS<-t(niche_breadth_Bla_ITS)
rownames(niche_breadth_Bla_ITS)<-NULL
niche_breadth_Bla_ITS<-as.data.frame(niche_breadth_Bla_ITS)
niche_breadth_Bla_ITS$Group<-c("Bla")
colnames(niche_breadth_Bla_ITS)[1]<-("Bcom")


niche_breadth_Gre_ITS<-t(niche_breadth_Gre_ITS)
rownames(niche_breadth_Gre_ITS)<-NULL
niche_breadth_Gre_ITS<-as.data.frame(niche_breadth_Gre_ITS)
niche_breadth_Gre_ITS$Group<-c("Gre")
colnames(niche_breadth_Gre_ITS)[1]<-("Bcom")


niche_breadth_Whi_ITS<-t(niche_breadth_Whi_ITS)
rownames(niche_breadth_Whi_ITS)<-NULL
niche_breadth_Whi_ITS<-as.data.frame(niche_breadth_Whi_ITS)
niche_breadth_Whi_ITS$Group<-c("Whi")
colnames(niche_breadth_Whi_ITS)[1]<-("Bcom")




niche_breadth_4_ITS_boxplot<-rbind(niche_breadth_CK_ITS,niche_breadth_Bla_ITS,
                                   niche_breadth_Gre_ITS,niche_breadth_Whi_ITS)



compaired <- list(c('CK', 'Bla'),c('CK', 'Gre'),
                  c('CK', 'Whi'))

wbwPalette_4 <- c("#EEA236FF","#D43F3AFF","#5CB85CFF", "#46B8DAFF")


niche_breadth_4_ITS_boxplot$Group <- factor(niche_breadth_4_ITS_boxplot$Group,levels=c('CK','Bla','Gre','Whi'))

p3<- ggplot(niche_breadth_4_ITS_boxplot, aes(Group, Bcom)) + 
  stat_boxplot(geom = 'errorbar', width = 0.4, size=0.2) + 
  geom_boxplot(fill= wbwPalette_4, width = 0.6, size=0.3, alpha=0.5,outlier.shape = NA) + 
  stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = -0.1) + 
  theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
  geom_violin(aes(fill = Group), width = 0.9, size = 0.5, alpha = 0.6)+ 
  scale_fill_manual(values=wbwPalette_4)+
  ylab("Bcom")+xlab(NULL)+ 
  ylim(0,8.5)+ 
  geom_text(aes(x = 1.2,y = 0.2,label = c("Eukaryotes")))+
  theme(legend.position="none") 

p3
#save 3inches*4inches











tem_a_16S<-ASV_all_16S_With_Tax[which(rownames(ASV_all_16S_With_Tax)%in%rownames(niche_breadth_16S_ASV)),]
tmp_b_16S<-print(row.names(tem_a_16S)==row.names(niche_breadth_16S_ASV))
summary(tmp_b_16S)

niche_breadth_16S_ASV_Tax<-cbind(niche_breadth_16S_ASV,tem_a_16S$Phylum)
colnames(niche_breadth_16S_ASV_Tax)<-c("Bcom","Phylum")
niche_breadth_16S_ASV_Tax<-as.data.frame(niche_breadth_16S_ASV_Tax)

niche_breadth_16S_ASV_Tax$Type[niche_breadth_16S_ASV_Tax$Bcom < 1.5] <- "Generalists"
niche_breadth_16S_ASV_Tax$Type[niche_breadth_16S_ASV_Tax$Bcom >= 1.5 &
                                 niche_breadth_16S_ASV_Tax$Bcom <= 3] <- "Middle"
niche_breadth_16S_ASV_Tax$Type[niche_breadth_16S_ASV_Tax$Bcom > 3] <- "Specialists"

Top_Phylum_16S_table<-read.csv(file = "Leshan_16S_Phyla.csv",header=TRUE,sep = ',',row.names = 1 )
Top_Phylum_16S<-row.names(Top_Phylum_16S_table)[1:10]
niche_breadth_16S_ASV_Tax$Phylum<-ifelse(niche_breadth_16S_ASV_Tax$Phylum %in% Top_Phylum_16S, niche_breadth_16S_ASV_Tax$Phylum, c("Others"))

niche_breadth_16S_ASV_Tax$Bcom<-as.numeric(niche_breadth_16S_ASV_Tax$Bcom)


library(tidyverse)
Summarise_niche_breadth_16S_ASV_Tax<- niche_breadth_16S_ASV_Tax %>% group_by(Phylum,Type) %>% summarise(n = n()) %>%
  group_by(Phylum) %>% mutate(Phylum_sum = sum(n)) %>% group_by(Phylum,Type) %>% mutate(Relative=n/Phylum_sum)





library(RColorBrewer)
library(ggplot2)
library(ggalluvial)
display.brewer.all()
display.brewer.pal(11,"Spectral")
Palette_niche<-brewer.pal(11,"Spectral")[c(3,9,10)]

Summarise_niche_breadth_16S_ASV_Tax$Type <- factor(Summarise_niche_breadth_16S_ASV_Tax$Type,levels=c('Specialists','Middle','Generalists'))
Summarise_niche_breadth_16S_ASV_Tax$Phylum<-factor(Summarise_niche_breadth_16S_ASV_Tax$Phylum,levels=c(Top_Phylum_16S,"Others"))

write.table(Summarise_niche_breadth_16S_ASV_Tax,
            file = "niche_breadth_16S_ASV_Tax.txt",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = "\t")

p4 <- ggplot(data = Summarise_niche_breadth_16S_ASV_Tax,aes(x = Phylum, y = Relative, alluvium = Type, stratum = Type))+
  geom_alluvium(aes(fill = Type),alpha = .5,width = 0.6) + 
  geom_stratum(aes(fill = Type),width = 0.6) +
 ylab(label = "ASVs Percentage") + xlab(label = NULL) +
 scale_fill_manual(values = Palette_niche) +
  theme_bw()+ theme(panel.grid=element_blank()) + 
  theme(panel.border = element_blank()) +
  theme(panel.background=element_rect(fill='transparent', color='black'),plot.margin = unit(c(3,5,1,1),"mm")) +
  theme(axis.text.x=element_text(colour="black",size=8,face = "bold",angle = -45,vjust = 0.5,hjust = 0)) + 
  theme(axis.text.y=element_text(colour = "black",size = 8)) + 
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.title.y = element_text(size = 8,face = "bold",margin = unit(c(0,1,0,1),"lines")))+
  scale_y_continuous(expand = c(0,0))+
  theme(legend.text = element_text(colour = "black",size = 15)) + 
  theme(legend.title = element_text(size = 10,colour = "black")) + 
  labs(title = "Prokaryotes")+
  theme(text=element_text(size=10))+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")

p4
#save pdf 3*8inches








Summarise_niche_breadth_16S_ALL_ASV_Tax<- niche_breadth_16S_ASV_Tax %>% group_by(Type) %>% summarise(n = n()) %>% mutate(Phylum_sum = sum(n)) %>% mutate(Relative = n/Phylum_sum ) 

Donut_niche_16S<-Summarise_niche_breadth_16S_ALL_ASV_Tax


library(RColorBrewer)
library(ggplot2)
display.brewer.all()
display.brewer.pal(11,"Spectral")
Palette_niche<-brewer.pal(11,"Spectral")[c(3,9,10)]


Donut_niche_16S$Type <- factor(Donut_niche_16S$Type,levels=c('Specialists','Middle','Generalists'))



Donut_niche_16S$ymax <- cumsum(Donut_niche_16S$Relative)

Donut_niche_16S$ymin <- c(0, head(Donut_niche_16S$ymax, n=-1))

Donut_niche_16S$labelPosition <- (Donut_niche_16S$ymax + Donut_niche_16S$ymin) / 2

Donut_niche_16S$label <- paste0(Donut_niche_16S$Type, "\n ", round(Donut_niche_16S$Relative*100, 2), "%")



p5<- ggplot(Donut_niche_16S, aes(ymax=ymax, ymin=ymin, xmax=3, xmin=1, fill=Type)) +
  geom_rect() +
  geom_text( x=5, aes(y=labelPosition, label=label, color=Type), size=4) + # x here controls label position (inner / outer)
  geom_text( x=-1, y=0, label=paste("Prokaryotes"), size=4) +
  scale_fill_manual(values= c(Palette_niche)) +
  scale_color_manual(values=c(Palette_niche)) +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")
p5
#save pdf 3*3inches











tem_a_ITS<-ASV_all_ITS_With_Tax[which(rownames(ASV_all_ITS_With_Tax)%in%rownames(niche_breadth_ITS_ASV)),]
tmp_b_ITS<-print(row.names(tem_a_ITS)==row.names(niche_breadth_ITS_ASV))
summary(tmp_b_ITS)

niche_breadth_ITS_ASV_Tax<-cbind(niche_breadth_ITS_ASV,tem_a_ITS$Phylum)
colnames(niche_breadth_ITS_ASV_Tax)<-c("Bcom","Phylum")
niche_breadth_ITS_ASV_Tax<-as.data.frame(niche_breadth_ITS_ASV_Tax)

niche_breadth_ITS_ASV_Tax$Type[niche_breadth_ITS_ASV_Tax$Bcom < 1.5] <- "Generalists"
niche_breadth_ITS_ASV_Tax$Type[niche_breadth_ITS_ASV_Tax$Bcom >= 1.5 &
                                 niche_breadth_ITS_ASV_Tax$Bcom <= 3] <- "Middle"
niche_breadth_ITS_ASV_Tax$Type[niche_breadth_ITS_ASV_Tax$Bcom > 3] <- "Specialists"


Top_Phylum_ITS_table<-read.csv(file = "Leshan_ITS_Phyla.csv",header=TRUE,sep = ',',row.names = 1 )
Top_Phylum_ITS<-row.names(Top_Phylum_ITS_table)[1:5]
niche_breadth_ITS_ASV_Tax$Phylum<-ifelse(niche_breadth_ITS_ASV_Tax$Phylum %in% Top_Phylum_ITS, niche_breadth_ITS_ASV_Tax$Phylum, c("Others"))

niche_breadth_ITS_ASV_Tax$Bcom<-as.numeric(niche_breadth_ITS_ASV_Tax$Bcom)


library(tidyverse)
Summarise_niche_breadth_ITS_ASV_Tax<- niche_breadth_ITS_ASV_Tax %>% group_by(Phylum,Type) %>% summarise(n = n()) %>%
  group_by(Phylum) %>% mutate(Phylum_sum = sum(n)) %>% group_by(Phylum,Type) %>% mutate(Relative=n/Phylum_sum)





library(RColorBrewer)
library(ggplot2)
library(ggalluvial)
display.brewer.all()
display.brewer.pal(11,"Spectral")
Palette_niche<-brewer.pal(11,"Spectral")[c(3,9,10)]

Summarise_niche_breadth_ITS_ASV_Tax$Type <- factor(Summarise_niche_breadth_ITS_ASV_Tax$Type,levels=c('Specialists','Middle','Generalists'))
Summarise_niche_breadth_ITS_ASV_Tax$Phylum<-factor(Summarise_niche_breadth_ITS_ASV_Tax$Phylum,levels=c(Top_Phylum_ITS,"Others"))


write.table(Summarise_niche_breadth_ITS_ASV_Tax,
            file = "niche_breadth_ITS_ASV_Tax.txt",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = "\t")

p6 <- ggplot(data = Summarise_niche_breadth_ITS_ASV_Tax,aes(x = Phylum, y = Relative, alluvium = Type, stratum = Type))+
  geom_alluvium(aes(fill = Type),alpha = .5,width = 0.6) + 
  geom_stratum(aes(fill = Type),width = 0.6) +
  ylab(label = "ASVs Percentage") + xlab(label = NULL) +
  scale_fill_manual(values = Palette_niche) +
  theme_bw()+ theme(panel.grid=element_blank()) + 
  theme(panel.border = element_blank()) +
  theme(panel.background=element_rect(fill='transparent', color='black'),plot.margin = unit(c(3,5,1,1),"mm")) +
  theme(axis.text.x=element_text(colour="black",size=8,face = "bold",angle = -45,vjust = 0.5,hjust = 0)) + 
  theme(axis.text.y=element_text(colour = "black",size = 8)) + 
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.title.y = element_text(size = 8,face = "bold",margin = unit(c(0,1,0,1),"lines")))+
  scale_y_continuous(expand = c(0,0))+
  theme(legend.text = element_text(colour = "black",size = 15)) + 
  theme(legend.title = element_text(size = 10,colour = "black")) + 
  labs(title = "Eukaryotes")+
  theme(text=element_text(size=10))+
  theme(plot.title = element_text(hjust = 0.5)) 

p6
#save pdf 3*8inches








Summarise_niche_breadth_ITS_ALL_ASV_Tax<- niche_breadth_ITS_ASV_Tax %>% group_by(Type) %>% summarise(n = n()) %>% mutate(Phylum_sum = sum(n)) %>% mutate(Relative = n/Phylum_sum ) 

Donut_niche_ITS<-Summarise_niche_breadth_ITS_ALL_ASV_Tax


library(RColorBrewer)
library(ggplot2)
display.brewer.all()
display.brewer.pal(11,"Spectral")
Palette_niche<-brewer.pal(11,"Spectral")[c(3,9,10)]


Donut_niche_ITS$Type <- factor(Donut_niche_ITS$Type,levels=c('Specialists','Middle','Generalists'))


Donut_niche_ITS$ymax <- cumsum(Donut_niche_ITS$Relative)

Donut_niche_ITS$ymin <- c(0, head(Donut_niche_ITS$ymax, n=-1))

Donut_niche_ITS$labelPosition <- (Donut_niche_ITS$ymax + Donut_niche_ITS$ymin) / 2

Donut_niche_ITS$label <- paste0(Donut_niche_ITS$Type, "\n ", round(Donut_niche_ITS$Relative*100, 2), "%")



p7<- ggplot(Donut_niche_ITS, aes(ymax=ymax, ymin=ymin, xmax=3, xmin=1, fill=Type)) +
  geom_rect() +
  geom_text( x=5, aes(y=labelPosition, label=label, color=Type), size=4) + # x here controls label position (inner / outer)
  geom_text( x=-1, y=0, label=paste("Eukaryotes"), size=4) +
  scale_fill_manual(values= c(Palette_niche)) +
  scale_color_manual(values=c(Palette_niche)) +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")
p7
#save pdf 3*3inches





