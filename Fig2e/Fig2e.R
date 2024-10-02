rm(list=ls())
library(RColorBrewer)
display.brewer.all()
display.brewer.pal(12,'Set3')
brewer.pal(12,'Set3')

Palette <- c(brewer.pal(12,'Set3'))
Palette<-Palette[c(1:8,10:11,9)]
Palette<-Palette[c(4:6,1:3,7:11)]


rm(list=ls())

library(vegan)
library(ggplot2)
library(ggsci)

Group_infor<-read.table("Group_Leshan.txt",header = TRUE,sep = "\t",row.names = 1)

ASV_all<-read.table("ASV_Abundance_split.txt",sep = "\t", header = TRUE,row.names = 1)
ASV_all<-ASV_all[,1:c(length(colnames(ASV_all))-7)]


ASV_filtered<-ASV_all[which(rowSums(ASV_all)>27),]

ASV_filtered<-t(ASV_filtered)
ASV_filtered<-as.data.frame(ASV_filtered)


distance <- vegdist(ASV_filtered, method = 'bray')
nmds <- metaMDS(ASV_filtered, k = 2)



summary(nmds)

stress <- nmds$stress

df <- as.data.frame(nmds$points)

rownames(df)==rownames(Group_infor)
df <- cbind(df, Group_infor)

df<-df[order(rownames(df)),]
df$Order <- c(1, 2, 6, 3, 4, 5, 
              1, 3, 4, 5, 2, 6,
              1, 2, 3, 6, 5, 4, 8, 7,
              1, 2, 7, 3, 4, 6, 5)

df$Order_1<-paste(df$Group, df$Order, sep = "_")
df <- df[order(df$Order_1), ]




set.seed(123)
# PERMANOVA
adonis <- adonis2(ASV_filtered ~ Group, data = Group_infor, permutations = 999, method = "bray")
# anosim
anosim = anosim(ASV_filtered, Group_infor$Group, permutations = 999, distance = "bray")



stress_text <- paste("Stress  =", round(stress, 4))
adonis_text <- paste(paste("Adonis  =", round(adonis$R2, 2)), "***")[1] # *** P < 0.001
anosim_text <- paste(paste("Anosim  =", round(anosim$statistic, 2)), "***") # *** P < 0.001


p <- ggplot(df, aes(MDS1, MDS2))+
  geom_point(aes(color = Group), size = 4,)+
  geom_polygon(aes(x = MDS1, y = MDS2, fill = Group, group = Group, color = Group),
               alpha = 0.3, linetype = "longdash")+
  theme(plot.margin = unit(rep(1, 4), 'lines'), 
        panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"), 
        panel.grid = element_blank(), 
        panel.background = element_rect(fill = 'white'),
        legend.position="right")+
  guides(color = "none", fill = "none")+
  ggtitle(paste(paste(stress_text, adonis_text), anosim_text)) + 
  scale_color_locuszoom()+
  scale_fill_locuszoom()   

p

# save pdf 5*5 inches
#save pdf 5*8 inches



