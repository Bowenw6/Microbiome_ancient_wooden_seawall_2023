
rm(list=ls())




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
niche_width_16S$Group<-c("Bacteria")
colnames(niche_width_16S)[1]<-("Bcom")



write.csv(niche_width_16S,"niche_width_16S.csv",
          quote = FALSE,row.names = FALSE)






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
niche_width_ITS$Group<-c("Fungi")
colnames(niche_width_ITS)[1]<-("Bcom")


write.csv(niche_width_ITS,"niche_width_ITS.csv",
          quote = FALSE,row.names = FALSE)







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



niche_width_All_Boxplot<-rbind(niche_width_16S,
                                 niche_width_ITS)




library(RColorBrewer)
display.brewer.all()
display.brewer.pal(12,"Set3")
brewer.pal(9,"Set1")
wbwPalette <- brewer.pal(9,"Set1")[c(5,4)]
library(ggplot2)
library(ggpubr)
compaired <- list(c('Bacteria', 'Fungi'))

niche_width_All_Boxplot$Group <- factor(niche_width_All_Boxplot$Group,levels=c('Bacteria', 'Fungi'))


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
  ggtitle("Bcom-Bacteria")+
  theme(plot.title = element_text(hjust = 0.5))

p12











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
  ggtitle("Bcom-Fungi")+
  theme(plot.title = element_text(hjust = 0.5))

p13





library(ggpubr)
p1 <- ggarrange(ggarrange(p12, p13, p11,
                          ncol = 3,
                          widths = c(0.37,0.37,0.26)))
p1
# save image 3inches*11inches















tem_a_16S<-ASV_all_16S_With_Tax[which(rownames(ASV_all_16S_With_Tax)%in%rownames(niche_width_16S_ASV)),]
tmp_b_16S<-print(row.names(tem_a_16S)==row.names(niche_width_16S_ASV))
summary(tmp_b_16S)

niche_width_16S_ASV_Tax<-cbind(niche_width_16S_ASV,tem_a_16S$Phylum)
colnames(niche_width_16S_ASV_Tax)<-c("Bcom","Phylum")
niche_width_16S_ASV_Tax<-as.data.frame(niche_width_16S_ASV_Tax)

niche_width_16S_ASV_Tax$Type[niche_width_16S_ASV_Tax$Bcom < 1.5] <- "Generalists"
niche_width_16S_ASV_Tax$Type[niche_width_16S_ASV_Tax$Bcom >= 1.5 &
                                 niche_width_16S_ASV_Tax$Bcom <= 3] <- "Middle"
niche_width_16S_ASV_Tax$Type[niche_width_16S_ASV_Tax$Bcom > 3] <- "Specialists"


Top_Phylum_16S_table<-read.csv(file = "Leshan_16S_Phyla.csv",header=TRUE,sep = ',',row.names = 1 )
Top_Phylum_16S<-row.names(Top_Phylum_16S_table)[1:10]
niche_width_16S_ASV_Tax$Phylum<-ifelse(niche_width_16S_ASV_Tax$Phylum %in% Top_Phylum_16S, niche_width_16S_ASV_Tax$Phylum, c("Others"))
niche_width_16S_ASV_Tax$Bcom<-as.numeric(niche_width_16S_ASV_Tax$Bcom)

library(tidyverse)
Summarise_niche_width_16S_ASV_Tax<- niche_width_16S_ASV_Tax %>% group_by(Phylum,Type) %>% summarise(n = n()) %>%
  group_by(Phylum) %>% mutate(Phylum_sum = sum(n)) %>% group_by(Phylum,Type) %>% mutate(Relative=n/Phylum_sum)





library(RColorBrewer)
library(ggplot2)
library(ggalluvial)
display.brewer.all()
display.brewer.pal(11,"Spectral")
Palette_niche<-brewer.pal(11,"Spectral")[c(3,9,10)]
Summarise_niche_width_16S_ASV_Tax$Type <- factor(Summarise_niche_width_16S_ASV_Tax$Type,levels=c('Specialists','Middle','Generalists'))
Summarise_niche_width_16S_ASV_Tax$Phylum<-factor(Summarise_niche_width_16S_ASV_Tax$Phylum,levels=c(Top_Phylum_16S,"Others"))

write.table(Summarise_niche_width_16S_ASV_Tax,
            file = "niche_width_16S_ASV_Tax.txt",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = "\t")

p4 <- ggplot(data = Summarise_niche_width_16S_ASV_Tax,aes(x = Phylum, y = Relative, alluvium = Type, stratum = Type))+
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
  labs(title = "Bacteria")+
  theme(text=element_text(size=10))+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")
#theme(text = element_text(family = "Times"))
p4
#save pdf 3*8inches







Summarise_niche_width_16S_ALL_ASV_Tax<- niche_width_16S_ASV_Tax %>% group_by(Type) %>% summarise(n = n()) %>% mutate(Phylum_sum = sum(n)) %>% mutate(Relative = n/Phylum_sum ) 

Donut_niche_16S<-Summarise_niche_width_16S_ALL_ASV_Tax


library(RColorBrewer)
library(ggplot2)
display.brewer.all()
display.brewer.pal(11,"Spectral")
Palette_niche<-brewer.pal(11,"Spectral")[c(3,9,10)]

Donut_niche_16S$Type <- factor(Donut_niche_16S$Type,levels=c('Specialists','Middle','Generalists'))



# Compute the cumulative percentages (top of each rectangle)
Donut_niche_16S$ymax <- cumsum(Donut_niche_16S$Relative)
# Compute the bottom of each rectangle
Donut_niche_16S$ymin <- c(0, head(Donut_niche_16S$ymax, n=-1))
# Compute label position
Donut_niche_16S$labelPosition <- (Donut_niche_16S$ymax + Donut_niche_16S$ymin) / 2
# Compute a good label
Donut_niche_16S$label <- paste0(Donut_niche_16S$Type, "\n ", round(Donut_niche_16S$Relative*100, 2), "%")


# Make the plot
p5<- ggplot(Donut_niche_16S, aes(ymax=ymax, ymin=ymin, xmax=3, xmin=1, fill=Type)) +
  geom_rect() +
  geom_text( x=5, aes(y=labelPosition, label=label, color=Type), size=4) + # x here controls label position (inner / outer)
  geom_text( x=-1, y=0, label=paste("Bacteria"), size=4) +
  scale_fill_manual(values= c(Palette_niche)) +
  scale_color_manual(values=c(Palette_niche)) +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")
p5
#save pdf 3*3inches










tem_a_ITS<-ASV_all_ITS_With_Tax[which(rownames(ASV_all_ITS_With_Tax)%in%rownames(niche_width_ITS_ASV)),]
tmp_b_ITS<-print(row.names(tem_a_ITS)==row.names(niche_width_ITS_ASV))
summary(tmp_b_ITS)

niche_width_ITS_ASV_Tax<-cbind(niche_width_ITS_ASV,tem_a_ITS$Phylum)
colnames(niche_width_ITS_ASV_Tax)<-c("Bcom","Phylum")
niche_width_ITS_ASV_Tax<-as.data.frame(niche_width_ITS_ASV_Tax)

niche_width_ITS_ASV_Tax$Type[niche_width_ITS_ASV_Tax$Bcom < 1.5] <- "Generalists"
niche_width_ITS_ASV_Tax$Type[niche_width_ITS_ASV_Tax$Bcom >= 1.5 &
                                 niche_width_ITS_ASV_Tax$Bcom <= 3] <- "Middle"
niche_width_ITS_ASV_Tax$Type[niche_width_ITS_ASV_Tax$Bcom > 3] <- "Specialists"


Top_Phylum_ITS_table<-read.csv(file = "Leshan_ITS_Phyla.csv",header=TRUE,sep = ',',row.names = 1 )
Top_Phylum_ITS<-row.names(Top_Phylum_ITS_table)[1:5]
niche_width_ITS_ASV_Tax$Phylum<-ifelse(niche_width_ITS_ASV_Tax$Phylum %in% Top_Phylum_ITS, niche_width_ITS_ASV_Tax$Phylum, c("Others"))

niche_width_ITS_ASV_Tax$Bcom<-as.numeric(niche_width_ITS_ASV_Tax$Bcom)


library(tidyverse)
Summarise_niche_width_ITS_ASV_Tax<- niche_width_ITS_ASV_Tax %>% group_by(Phylum,Type) %>% summarise(n = n()) %>%
  group_by(Phylum) %>% mutate(Phylum_sum = sum(n)) %>% group_by(Phylum,Type) %>% mutate(Relative=n/Phylum_sum)





library(RColorBrewer)
library(ggplot2)
library(ggalluvial)
display.brewer.all()
display.brewer.pal(11,"Spectral")
Palette_niche<-brewer.pal(11,"Spectral")[c(3,9,10)]

Summarise_niche_width_ITS_ASV_Tax$Type <- factor(Summarise_niche_width_ITS_ASV_Tax$Type,levels=c('Specialists','Middle','Generalists'))
Summarise_niche_width_ITS_ASV_Tax$Phylum<-factor(Summarise_niche_width_ITS_ASV_Tax$Phylum,levels=c(Top_Phylum_ITS,"Others"))


write.table(Summarise_niche_width_ITS_ASV_Tax,
            file = "niche_width_ITS_ASV_Tax.txt",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = "\t")

p6 <- ggplot(data = Summarise_niche_width_ITS_ASV_Tax,aes(x = Phylum, y = Relative, alluvium = Type, stratum = Type))+
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
  labs(title = "Fungi")+
  theme(text=element_text(size=10))+
  theme(plot.title = element_text(hjust = 0.5)) 

p6
#save pdf 3*8inches







Summarise_niche_width_ITS_ALL_ASV_Tax<- niche_width_ITS_ASV_Tax %>% group_by(Type) %>% summarise(n = n()) %>% mutate(Phylum_sum = sum(n)) %>% mutate(Relative = n/Phylum_sum ) 

Donut_niche_ITS<-Summarise_niche_width_ITS_ALL_ASV_Tax


library(RColorBrewer)
library(ggplot2)
display.brewer.all()
display.brewer.pal(11,"Spectral")
Palette_niche<-brewer.pal(11,"Spectral")[c(3,9,10)]


Donut_niche_ITS$Type <- factor(Donut_niche_ITS$Type,levels=c('Specialists','Middle','Generalists'))



# Compute the cumulative percentages (top of each rectangle)
Donut_niche_ITS$ymax <- cumsum(Donut_niche_ITS$Relative)
# Compute the bottom of each rectangle
Donut_niche_ITS$ymin <- c(0, head(Donut_niche_ITS$ymax, n=-1))
# Compute label position
Donut_niche_ITS$labelPosition <- (Donut_niche_ITS$ymax + Donut_niche_ITS$ymin) / 2
# Compute a good label
Donut_niche_ITS$label <- paste0(Donut_niche_ITS$Type, "\n ", round(Donut_niche_ITS$Relative*100, 2), "%")


# Make the plot
p7<- ggplot(Donut_niche_ITS, aes(ymax=ymax, ymin=ymin, xmax=3, xmin=1, fill=Type)) +
  geom_rect() +
  geom_text( x=5, aes(y=labelPosition, label=label, color=Type), size=4) + # x here controls label position (inner / outer)
  geom_text( x=-1, y=0, label=paste("Fungi"), size=4) +
  scale_fill_manual(values= c(Palette_niche)) +
  scale_color_manual(values=c(Palette_niche)) +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")
p7
#save pdf 3*3inches







library(ggpubr)
p2 <- ggarrange(ggarrange(p4, p5, 
                          ncol = 2,
                          widths = c(0.74,0.26)),
                ggarrange(p6, p7, 
                          ncol = 2,
                          widths = c(0.74,0.26)),
                nrow = 2, heights = c(0.5, 0.5))
p2

#save image 6*11 inches


