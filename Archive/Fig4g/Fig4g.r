#niche width and niche overlap
rm(list=ls())


########## input files from Fig4abchijk_TableS6_S7 ####




library(tidyverse)
library(ggalluvial)
library(ggplot2)

library(RColorBrewer)

display.brewer.all()
display.brewer.pal(11,"Spectral")
Palette_niche<-c("#4991ba","#ffc357","#9874a7")



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










niche_width_4_16S_Boxplot<-rbind(niche_width_CK_16S,
                                 niche_width_Bla_16S,
                                 niche_width_Gre_16S,
                                 niche_width_Whi_16S)
niche_width_4_16S_Boxplot$Group <- factor(niche_width_4_16S_Boxplot$Group,levels=c('CK', 'Bla','Gre','Whi'))




niche_width_4_ITS_Boxplot<-rbind(niche_width_CK_ITS,
                                 niche_width_Bla_ITS,
                                 niche_width_Gre_ITS,
                                 niche_width_Whi_ITS)

niche_width_4_ITS_Boxplot$Group <- factor(niche_width_4_ITS_Boxplot$Group,levels=c('CK', 'Bla','Gre','Whi'))





niche_width_All_Boxplot<-rbind(niche_width_16S,
                               niche_width_ITS)

niche_width_All_Boxplot$Group <- factor(niche_width_All_Boxplot$Group,levels=c('Prokaryotes', 'Eukaryotes'))




rm(niche_width_16S,niche_width_ITS,
   niche_width_Bla_16S,niche_width_Bla_ITS,
   niche_width_CK_16S,niche_width_CK_ITS,
   niche_width_Gre_16S,niche_width_Gre_ITS,
   niche_width_Whi_16S,niche_width_Whi_ITS)







# all
niche_width_All_Boxplot$Type[niche_width_All_Boxplot$Bcom > 5] <- "Generalists"
niche_width_All_Boxplot$Type[niche_width_All_Boxplot$Bcom >= 1.5 &
                               niche_width_All_Boxplot$Bcom <= 5] <- "Middle"
niche_width_All_Boxplot$Type[niche_width_All_Boxplot$Bcom < 1.5] <- "Specialists"

# 16S
niche_width_4_16S_Boxplot$Type[niche_width_4_16S_Boxplot$Bcom > 5] <- "Generalists"
niche_width_4_16S_Boxplot$Type[niche_width_4_16S_Boxplot$Bcom >= 1.5 &
                                 niche_width_4_16S_Boxplot$Bcom <= 5] <- "Middle"
niche_width_4_16S_Boxplot$Type[niche_width_4_16S_Boxplot$Bcom < 1.5] <- "Specialists"

# ITS
niche_width_4_ITS_Boxplot$Type[niche_width_4_ITS_Boxplot$Bcom > 5] <- "Generalists"
niche_width_4_ITS_Boxplot$Type[niche_width_4_ITS_Boxplot$Bcom >= 1.5 &
                                 niche_width_4_ITS_Boxplot$Bcom <= 5] <- "Middle"
niche_width_4_ITS_Boxplot$Type[niche_width_4_ITS_Boxplot$Bcom < 1.5] <- "Specialists"












Summarise_niche_width_4_16S_Boxplot<- niche_width_4_16S_Boxplot %>% 
  group_by(Group,Type) %>% 
  summarise(n = n()) %>%
  group_by(Group) %>% 
  mutate(Group_sum = sum(n)) %>% 
  group_by(Group,Type) %>% 
  mutate(Relative=n/Group_sum)

Summarise_niche_width_4_ITS_Boxplot<- niche_width_4_ITS_Boxplot %>% 
  group_by(Group,Type) %>% 
  summarise(n = n()) %>%
  group_by(Group) %>% 
  mutate(Group_sum = sum(n)) %>% 
  group_by(Group,Type) %>% 
  mutate(Relative=n/Group_sum)

Summarise_niche_width_All_Boxplot<- niche_width_All_Boxplot %>% 
  group_by(Group,Type) %>% 
  summarise(n = n()) %>%
  group_by(Group) %>% 
  mutate(Group_sum = sum(n)) %>% 
  group_by(Group,Type) %>% 
  mutate(Relative=n/Group_sum)



Summarise_niche_width_4_16S_Boxplot$Type <- factor(Summarise_niche_width_4_16S_Boxplot$Type,levels=c('Specialists','Middle','Generalists'))
Summarise_niche_width_4_ITS_Boxplot$Type <- factor(Summarise_niche_width_4_ITS_Boxplot$Type,levels=c('Specialists','Middle','Generalists'))
Summarise_niche_width_All_Boxplot$Type <- factor(Summarise_niche_width_All_Boxplot$Type,levels=c('Specialists','Middle','Generalists'))







##### plot 16S ###

p1 <- ggplot(data = Summarise_niche_width_4_16S_Boxplot,aes(x = Group, y = Relative, alluvium = Type, stratum = Type))+
  geom_alluvium(aes(fill = Type),alpha = .5,width = 0.6) + 
  geom_stratum(aes(fill = Type),width = 0.6) +
  ylab(label = "ASVs Percentage") + xlab(label = NULL) +
  scale_fill_manual(values = Palette_niche) +
  theme_bw()+ 
  theme(panel.grid=element_blank()) + 
  theme(panel.border = element_blank()) +
  theme(panel.background=element_rect(fill='transparent', color='black'),plot.margin = unit(c(3,5,1,1),"mm")) +
  theme(axis.text.x=element_text(colour="black",size=8,face = "bold",
                                 vjust = 0.5,
                                 hjust = 0.5)) + 
  theme(axis.text.y=element_text(colour = "black",size = 8)) + 
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.title.y = element_text(size = 8,face = "bold"))+
  scale_y_continuous(expand = c(0,0))+
  theme(legend.text = element_text(colour = "black",size = 15)) + 
  theme(legend.title = element_text(size = 10,colour = "black")) + 
  labs(title = "Prokaryotes")+
  theme(text=element_text(size=10))+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")
p1







##### plot ITS ###

p2 <- ggplot(data = Summarise_niche_width_4_ITS_Boxplot,aes(x = Group, y = Relative, alluvium = Type, stratum = Type))+
  geom_alluvium(aes(fill = Type),alpha = .5,width = 0.6) + 
  geom_stratum(aes(fill = Type),width = 0.6) +
  ylab(label = "ASVs Percentage") + xlab(label = NULL) +
  scale_fill_manual(values = Palette_niche) +
  theme_bw()+ 
  theme(panel.grid=element_blank()) + 
  theme(panel.border = element_blank()) +
  theme(panel.background=element_rect(fill='transparent', color='black'),plot.margin = unit(c(3,5,1,1),"mm")) +
  theme(axis.text.x=element_text(colour="black",size=8,face = "bold",
                                 vjust = 0.5,
                                 hjust = 0.5)) + 
  theme(axis.text.y=element_text(colour = "black",size = 8)) + 
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.title.y = element_text(size = 8,face = "bold"))+
  scale_y_continuous(expand = c(0,0))+
  theme(legend.text = element_text(colour = "black",size = 15)) + 
  theme(legend.title = element_text(size = 10,colour = "black")) + 
  labs(title = "Eukaryotes")+
  theme(text=element_text(size=10))+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")
p2








##### plot all ###

p3 <- ggplot(data = Summarise_niche_width_All_Boxplot,aes(x = Group, y = Relative, alluvium = Type, stratum = Type))+
  geom_alluvium(aes(fill = Type),alpha = .5,width = 0.6) + 
  geom_stratum(aes(fill = Type),width = 0.6) +
  ylab(label = "ASVs Percentage") + xlab(label = NULL) +
  scale_fill_manual(values = Palette_niche) +
  theme_bw()+ 
  theme(panel.grid=element_blank()) + 
  theme(panel.border = element_blank()) +
  theme(panel.background=element_rect(fill='transparent', color='black'),plot.margin = unit(c(3,5,1,1),"mm")) +
  theme(axis.text.x=element_text(colour="black",size=8,face = "bold",
                                 vjust = 0.5,
                                 hjust = 0.5)) + 
  theme(axis.text.y=element_text(colour = "black",size = 8)) + 
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.title.y = element_text(size = 8,face = "bold"))+
  scale_y_continuous(expand = c(0,0))+
  theme(legend.text = element_text(colour = "black",size = 15)) + 
  theme(legend.title = element_text(size = 10,colour = "black")) + 
  labs(title = "All")+
  theme(text=element_text(size=10))+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")
p3











##### export 
library(ggpubr)
p_all <- ggarrange(p1, p2, p3, 
                   ncol = 3, nrow = 1,
                   widths = c(0.37,0.37,0.26))

  
p_all

ggsave("Raw_Generalists_and_Specialists_Bcom.pdf", plot = p_all, width = 11, height = 2.5, units = "in")


# save image 3*11 inches
# Raw_Generalists_and_Specialists_Bcom.pdf

