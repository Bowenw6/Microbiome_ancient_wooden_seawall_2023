
rm(list=ls())

ASV_all_16S_With_Tax<-read.table("ASV_Abundance_16S_split.txt",sep = "\t", header = TRUE,row.names = 1)
ASV_all_16S<-ASV_all_16S_With_Tax[,1:c(length(colnames(ASV_all_16S_With_Tax))-7)]


ASV_filtered_16S<- ASV_all_16S[rowSums(ASV_all_16S)>27,]


ASV_filtered_16S<-ASV_filtered_16S[,order(colnames(ASV_filtered_16S))]


colSums(ASV_filtered_16S)


library(spaa)
library(tidyverse)
t_ASV_filtered_16S <- t(ASV_filtered_16S)

niche_overlap_16S <- niche.overlap(t_ASV_filtered_16S, method = "pianka")%>% 
  data.matrix() %>% 
  as_tibble(rownames = "otuid") ->niche_overlap_16S


niche_overlap_16S[upper.tri(niche_overlap_16S,diag=FALSE)]=NA
niche_overlap_16S %>% 
  reshape2::melt(id.vars = "otuid") %>% 
  filter(value!="NA") %>% 
  mutate(paired = str_c(otuid,variable,sep = "_") ) -> niche_overlap_16S_Boxplot

niche_overlap_16S_Boxplot_tmp<-niche_overlap_16S_Boxplot$value
niche_overlap_16S_Boxplot_tmp<-as.data.frame(niche_overlap_16S_Boxplot_tmp)
niche_overlap_16S_Boxplot_tmp$Group<-c("Bacteria")
colnames(niche_overlap_16S_Boxplot_tmp)[1]<-("Ocom")


write.csv(niche_overlap_16S_Boxplot_tmp,"niche_overlap_16S_Boxplot_tmp.csv",
          quote = FALSE,row.names = FALSE)






ASV_all_ITS_With_Tax<-read.table("ASV_Abundance_ITS_split.txt",sep = "\t", header = TRUE,row.names = 1)
ASV_all_ITS<-ASV_all_ITS_With_Tax[,1:c(length(colnames(ASV_all_ITS_With_Tax))-7)]


ASV_filtered_ITS<- ASV_all_ITS[rowSums(ASV_all_ITS)>27,]



ASV_filtered_ITS<-ASV_filtered_ITS[,order(colnames(ASV_filtered_ITS))]


colSums(ASV_filtered_ITS)

library(spaa)
library(tidyverse)
t_ASV_filtered_ITS <- t(ASV_filtered_ITS)

niche_overlap_ITS <- niche.overlap(t_ASV_filtered_ITS, method = "pianka")%>% 
  data.matrix() %>% 
  as_tibble(rownames = "otuid") ->niche_overlap_ITS


niche_overlap_ITS[upper.tri(niche_overlap_ITS,diag=FALSE)]=NA
niche_overlap_ITS %>% 
  reshape2::melt(id.vars = "otuid") %>% 
  filter(value!="NA") %>% 
  mutate(paired = str_c(otuid,variable,sep = "_") ) -> niche_overlap_ITS_Boxplot

niche_overlap_ITS_Boxplot_tmp<-niche_overlap_ITS_Boxplot$value
niche_overlap_ITS_Boxplot_tmp<-as.data.frame(niche_overlap_ITS_Boxplot_tmp)
niche_overlap_ITS_Boxplot_tmp$Group<-c("Fungi")
colnames(niche_overlap_ITS_Boxplot_tmp)[1]<-("Ocom")

write.csv(niche_overlap_ITS_Boxplot_tmp,"niche_overlap_ITS_Boxplot_tmp.csv",
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
  
  
  niche_overlap <- niche.overlap(get(df_var), method = "pianka") %>%
    data.matrix() %>%
    as_tibble(rownames = "otuid")
  
  niche_overlap[upper.tri(niche_overlap, diag = FALSE)] <- NA
  
  niche_overlap_Boxplot <- niche_overlap %>%
    reshape2::melt(id.vars = "otuid") %>%
    filter(value != "NA") %>%
    mutate(paired = str_c(otuid, variable, sep = "_"))
  
  niche_overlap_Boxplot_tmp <- as.data.frame(niche_overlap_Boxplot$value)
  niche_overlap_Boxplot_tmp$Group <- df_name
  colnames(niche_overlap_Boxplot_tmp)[1] <- "Ocom"
  
  write.csv(niche_overlap_Boxplot_tmp,
            paste0("niche_overlap_", df_name, "_16S_Boxplot_tmp.csv"),
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
  
  niche_overlap <- niche.overlap(get(df_var), method = "pianka") %>%
    data.matrix() %>%
    as_tibble(rownames = "otuid")
  
  niche_overlap[upper.tri(niche_overlap, diag = FALSE)] <- NA
  
  niche_overlap_Boxplot <- niche_overlap %>%
    reshape2::melt(id.vars = "otuid") %>%
    filter(value != "NA") %>%
    mutate(paired = str_c(otuid, variable, sep = "_"))
  
  niche_overlap_Boxplot_tmp <- as.data.frame(niche_overlap_Boxplot$value)
  niche_overlap_Boxplot_tmp$Group <- df_name
  colnames(niche_overlap_Boxplot_tmp)[1] <- "Ocom"
  
  write.csv(niche_overlap_Boxplot_tmp,
            paste0("niche_overlap_", df_name, "_ITS_Boxplot_tmp.csv"),
            quote = FALSE, row.names = FALSE)
}























library(tidyverse)


niche_overlap_16S_Boxplot_tmp<-read.csv(file = "niche_overlap_16S_Boxplot_tmp.csv",
                                    #    quote = FALSE,
                                        header = TRUE,
                                        sep = ",")


niche_overlap_ITS_Boxplot_tmp<-read.csv(file = "niche_overlap_ITS_Boxplot_tmp.csv",
                                        #    quote = FALSE,
                                        header = TRUE,
                                        sep = ",")



niche_overlap_CK_16S_Boxplot_tmp<-read.csv(file = "niche_overlap_CK_16S_Boxplot_tmp.csv",
                                        #    quote = FALSE,
                                        header = TRUE,
                                        sep = ",")


niche_overlap_Bla_16S_Boxplot_tmp<-read.csv(file = "niche_overlap_Bla_16S_Boxplot_tmp.csv",
                                           #    quote = FALSE,
                                           header = TRUE,
                                           sep = ",")


niche_overlap_Gre_16S_Boxplot_tmp<-read.csv(file = "niche_overlap_Gre_16S_Boxplot_tmp.csv",
                                            #    quote = FALSE,
                                            header = TRUE,
                                            sep = ",")

niche_overlap_Whi_16S_Boxplot_tmp<-read.csv(file = "niche_overlap_Whi_16S_Boxplot_tmp.csv",
                                            #    quote = FALSE,
                                            header = TRUE,
                                            sep = ",")



niche_overlap_CK_ITS_Boxplot_tmp<-read.csv(file = "niche_overlap_CK_ITS_Boxplot_tmp.csv",
                                           #    quote = FALSE,
                                           header = TRUE,
                                           sep = ",")

niche_overlap_Bla_ITS_Boxplot_tmp<-read.csv(file = "niche_overlap_Bla_ITS_Boxplot_tmp.csv",
                                            #    quote = FALSE,
                                            header = TRUE,
                                            sep = ",")


niche_overlap_Gre_ITS_Boxplot_tmp<-read.csv(file = "niche_overlap_Gre_ITS_Boxplot_tmp.csv",
                                            #    quote = FALSE,
                                            header = TRUE,
                                            sep = ",")
niche_overlap_Whi_ITS_Boxplot_tmp<-read.csv(file = "niche_overlap_Whi_ITS_Boxplot_tmp.csv",
                                            #    quote = FALSE,
                                            header = TRUE,
                                            sep = ",")




niche_overlap_All_Boxplot<-rbind(niche_overlap_16S_Boxplot_tmp,
                                 niche_overlap_ITS_Boxplot_tmp)





library(RColorBrewer)
display.brewer.all()
display.brewer.pal(12,"Set3")
brewer.pal(9,"Set1")

wbwPalette <- brewer.pal(9,"Set1")[c(5,4)]
library(ggplot2)
library(ggpubr)

compaired <- list(c('Bacteria', 'Fungi'))


niche_overlap_All_Boxplot$Group <- factor(niche_overlap_All_Boxplot$Group,levels=c('Bacteria', 'Fungi'))


p11<- ggplot(niche_overlap_All_Boxplot, aes(Group, Ocom)) + 
  stat_boxplot(geom = 'errorbar', width = 0.4, size=0.2) + 
  geom_boxplot(fill= wbwPalette, width = 0.6, size=0.3, alpha=0.5,outlier.shape = NA) + 
  stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = -0.1) + 
  theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
  geom_violin(aes(fill = Group), width = 0.9, size = 0.5, alpha = 0.6)+ 
  scale_fill_manual(values=wbwPalette)+
  ylab("Ocom")+xlab(NULL)+ 
  ylim(-0.1,1.2)+ 
  
  theme(legend.position="none")+ 
  ggtitle("Ocom-All")+
  theme(plot.title = element_text(hjust = 0.5))

p11










niche_overlap_4_16S_Boxplot<-rbind(niche_overlap_CK_16S_Boxplot_tmp,
                                   niche_overlap_Bla_16S_Boxplot_tmp,
                                   niche_overlap_Gre_16S_Boxplot_tmp,
                                   niche_overlap_Whi_16S_Boxplot_tmp)





library(RColorBrewer)
display.brewer.all()
display.brewer.pal(12,"Set3")
brewer.pal(9,"Set1")

library(ggplot2)
library(ggpubr)

compaired <- list(c('CK', 'Bla'),c('CK', 'Gre'),
                  c('CK', 'Whi'))

wbwPalette_4 <- c("#EEA236FF","#D43F3AFF","#5CB85CFF", "#46B8DAFF")


niche_overlap_4_16S_Boxplot$Group <- factor(niche_overlap_4_16S_Boxplot$Group,levels=c('CK','Bla','Gre','Whi'))


p12<- ggplot(niche_overlap_4_16S_Boxplot, aes(Group, Ocom)) + 
  stat_boxplot(geom = 'errorbar', width = 0.4, size=0.2) + 
  geom_boxplot(fill= wbwPalette_4, width = 0.6, size=0.3, alpha=0.5,outlier.shape = NA) + 
  stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = -0.1) + 
  theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
  geom_violin(aes(fill = Group), width = 0.9, size = 0.5, alpha = 0.6)+ 
  scale_fill_manual(values=wbwPalette_4)+
  ylab("Ocom")+xlab(NULL)+ 
  ylim(-0.1,1.4)+ 
 
  theme(legend.position="none")+ 
  ggtitle("Ocom-Bacteria")+
  theme(plot.title = element_text(hjust = 0.5))

p12















niche_overlap_4_ITS_Boxplot<-rbind(niche_overlap_CK_ITS_Boxplot_tmp,
                                   niche_overlap_Bla_ITS_Boxplot_tmp,
                                   niche_overlap_Gre_ITS_Boxplot_tmp,
                                   niche_overlap_Whi_ITS_Boxplot_tmp)


library(RColorBrewer)
display.brewer.all()
display.brewer.pal(12,"Set3")
brewer.pal(9,"Set1")

library(ggplot2)
library(ggpubr)

compaired <- list(c('CK', 'Bla'),c('CK', 'Gre'),
                  c('CK', 'Whi'))

wbwPalette_4 <- c("#EEA236FF","#D43F3AFF","#5CB85CFF", "#46B8DAFF")


niche_overlap_4_ITS_Boxplot$Group <- factor(niche_overlap_4_ITS_Boxplot$Group,levels=c('CK','Bla','Gre','Whi'))


p13<- ggplot(niche_overlap_4_ITS_Boxplot, aes(Group, Ocom)) + 
  stat_boxplot(geom = 'errorbar', width = 0.4, size=0.2) + 
  geom_boxplot(fill= wbwPalette_4, width = 0.6, size=0.3, alpha=0.5,outlier.shape = NA) + 
  stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = -0.1) + 
  theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
  geom_violin(aes(fill = Group), width = 0.9, size = 0.5, alpha = 0.6)+ 
  scale_fill_manual(values=wbwPalette_4)+
  ylab("Ocom")+xlab(NULL)+ 
  ylim(-0.1,1.4)+ 
  theme(legend.position="none")+ 
  ggtitle("Ocom-Fungi")+
  theme(plot.title = element_text(hjust = 0.5))

p13





library(ggpubr)
p1 <- ggarrange(ggarrange(p12, p13, p11,
                          ncol = 3,
                          widths = c(0.37,0.37,0.26)))
p1

#save image 3*11 inches

