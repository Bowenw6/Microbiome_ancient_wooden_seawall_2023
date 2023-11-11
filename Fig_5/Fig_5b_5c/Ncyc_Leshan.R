
rm(list=ls())


Ncyc_subfamily<-read.table("Ncyc_Genes_subfamily.txt", sep = "\t",header = TRUE)
Ncyc_colnames<-colnames(Ncyc_subfamily)
for (i in 1:length(Ncyc_colnames)) {
  Ncyc_subfamily[,i] <-substr(Ncyc_subfamily[,i],1,nchar(Ncyc_subfamily[,i])-1)
}
colnames(Ncyc_subfamily)<- c("Pathways","Gene_sub_families","Annotation")




Ncyc_abundance<- read.table("NCyc_abund.txt", sep = "\t", header = TRUE)

Ncyc_abundance<-Ncyc_abundance[-1,]


Ncyc_abundance$Sum<-apply(Ncyc_abundance[,2:4], MARGIN = 1, FUN = sum)

Ncyc_abundance<-subset(Ncyc_abundance, Ncyc_abundance$Sum !=0)


Ncyc_abundance_1<-t(Ncyc_abundance[,2:5])
Ncyc_abundance_1<-as.data.frame(Ncyc_abundance_1)
Ncyc_abundance_1$Sum<-apply(Ncyc_abundance_1, MARGIN = 1, FUN = sum)
Ncyc_abundance_1<-Ncyc_abundance_1[,1:ncol(Ncyc_abundance_1-1)]/Ncyc_abundance_1$Sum
Ncyc_abundance_1<-round(Ncyc_abundance_1*100,digits = 2) 
Ncyc_abundance_1<-t(Ncyc_abundance_1)
Ncyc_abundance_1<-as.data.frame(Ncyc_abundance_1)

Ncyc_abundance_1<-Ncyc_abundance_1[-nrow(Ncyc_abundance_1),]
Ncyc_abundance<-cbind(Ncyc_abundance[,1],Ncyc_abundance_1)
colnames(Ncyc_abundance)<-c("Gene","Black","CK","Green","White")
rm(Ncyc_abundance_1)





Ncyc_subfamily_filterd<-Ncyc_subfamily[which(Ncyc_subfamily$Gene_sub_families%in%Ncyc_abundance$Gene),]

write.table(Ncyc_subfamily_filterd, file = "Ncyc_subfamily_filter.txt", quote = FALSE,sep = "\t",
            col.names = TRUE,row.names = FALSE)

Ncyc_subfamily_duplicated<-Ncyc_subfamily_filterd[duplicated(Ncyc_subfamily_filterd$Gene_sub_families),]
Ncyc_subfamily_unique<-Ncyc_subfamily_filterd[!c(duplicated(Ncyc_subfamily_filterd$Gene_sub_families)),]




Ncyc_subfamily_unique<-Ncyc_subfamily_unique[order(Ncyc_subfamily_unique$Gene_sub_families),]
Ncyc_abundance<-Ncyc_abundance[order(Ncyc_abundance$Gene),]
Ncyc_all_unique<-cbind(Ncyc_abundance,Ncyc_subfamily_unique)


Ncyc_abundance_duplicated<-Ncyc_abundance[which(Ncyc_abundance$Gene%in%Ncyc_subfamily_duplicated$Gene_sub_families),]
Ncyc_subfamily_duplicated<-Ncyc_subfamily_duplicated[order(Ncyc_subfamily_duplicated$Gene_sub_families),]
Ncyc_abundance_duplicated<-Ncyc_abundance_duplicated[order(Ncyc_abundance_duplicated$Gene),]
Ncyc_all_duplicated<-cbind(Ncyc_abundance_duplicated,Ncyc_subfamily_duplicated)


Ncyc_all<-rbind(Ncyc_all_unique,Ncyc_all_duplicated)

write.table(Ncyc_all,file = "Ncyc_all.txt",quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)







library(dplyr)
Ncyc_all_pathway_abundance<-Ncyc_all[,2:6]

Ncyc_all_pathway_abundance %>% group_by(Pathways) %>% mutate(Black_sum = sum(Black)) %>%
  group_by(Pathways) %>% mutate(CK_sum = sum(CK)) %>% 
  group_by(Pathways) %>% mutate(Green_sum = sum(Green)) %>%
  group_by(Pathways) %>% mutate(White_sum = sum(White)) -> Ncyc_all_pathway_abundance

Ncyc_all_pathway_abundance<-Ncyc_all_pathway_abundance[,5:9]
Ncyc_all_pathway_abundance<-unique(Ncyc_all_pathway_abundance)



rownames(Ncyc_all_pathway_abundance)<-Ncyc_all_pathway_abundance$Pathways
Ncyc_all_pathway_abundance<-t(Ncyc_all_pathway_abundance)
Ncyc_all_pathway_abundance<-Ncyc_all_pathway_abundance[-1,]
colnames(Ncyc_all_pathway_abundance)
Ncyc_all_pathway_abundance<-as.data.frame(Ncyc_all_pathway_abundance)

Ncyc_all_pathway_abundance_tmp<-lapply(Ncyc_all_pathway_abundance,as.numeric)
Ncyc_all_pathway_abundance_tmp<-as.data.frame(Ncyc_all_pathway_abundance_tmp)
rownames(Ncyc_all_pathway_abundance_tmp)<-rownames(Ncyc_all_pathway_abundance)
Ncyc_all_pathway_abundance<-Ncyc_all_pathway_abundance_tmp
rm(Ncyc_all_pathway_abundance_tmp)

Ncyc_all_pathway_abundance<-t(Ncyc_all_pathway_abundance)
colnames(Ncyc_all_pathway_abundance)<-c("Black","CK","Green", "White")

Ncyc_all_pathway_abundance<-as.data.frame(Ncyc_all_pathway_abundance)
Ncyc_all_pathway_abundance$sum<-apply(Ncyc_all_pathway_abundance,1,sum)
Ncyc_all_pathway_abundance<-arrange(Ncyc_all_pathway_abundance,-sum)

Ncyc_all_pathway_abundance<-Ncyc_all_pathway_abundance[,1:ncol(Ncyc_all_pathway_abundance)-1]
Ncyc_all_pathway_abundance<-as.matrix(Ncyc_all_pathway_abundance)

Ncyc_all_pathway_abundance

Ncyc_all_pathway_abundance<-(colnames(Ncyc_all_pathway_abundance,levels=c("CK","Black","Green", "White")))







library(RColorBrewer)
display.brewer.all()
display.brewer.pal(9,'Set1')
brewer.pal(9,'Set1')
display.brewer.pal(9,'Set2')
brewer.pal(8,'Set2')

Palette<-c(brewer.pal(8,'Set2')[c(1,2,4,3)],
           brewer.pal(4,'Set1')[4],
           brewer.pal(8,'Set2')[c(7,5,6)],
           brewer.pal(8,'Set2')[8])





library(reshape2)
taxon <- melt(Ncyc_all_pathway_abundance)
colnames(taxon) <- c("Pathways","Samples","value")

taxon$Samples<-factor(taxon$Samples,levels=c("CK","Black","Green","White"))

library(ggplot2)
library(ggalluvial)

p <- ggplot(data = taxon,aes(x = Samples, y = value, alluvium = Pathways, stratum = Pathways))
p1 <- p + geom_alluvium(aes(fill = Pathways),alpha = .5,width = 0.6) + 
  geom_stratum(aes(fill = Pathways),width = 0.6)
p2 <- p1 + ylab(label = "Relative Abundance") + xlab(label = "")
p3 <- p2 + scale_fill_manual(values = Palette) +
  theme_bw()+ theme(panel.grid=element_blank()) + 
  theme(panel.border = element_blank()) +
  theme(panel.background=element_rect(fill='transparent', color='black'),plot.margin = unit(c(3,5,1,1),"mm")) +

  theme(axis.text.x=element_text(colour="black",size=15,face = "bold",angle = -45,vjust = 0.5,hjust = 0)) + 
  theme(axis.text.y=element_text(colour = "black",size = 15)) + 
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.title.y = element_text(size = 15,face = "bold",margin = unit(c(0,1,0,1),"lines")))+
  scale_y_continuous(expand = c(0,0))
p5 <- p3 + theme(legend.text = element_text(colour = "black",size = 15)) + 
  theme(legend.title = element_text(size = 16,colour = "black")) + 

  labs(title = "Gene abundance within each process")+
  theme(text=element_text(size=20))+
  theme(plot.title = element_text(hjust = 0.5)) 

p5
#save PDF_6*6inches





Mean_Ncyc_all_pathway_abundance<-rowMeans(Ncyc_all_pathway_abundance)

write.table(Mean_Ncyc_all_pathway_abundance,file = "Mean_Ncyc_all_pathway_abundance.txt",
            quote = FALSE)






library(dplyr)

Ncyc_heatmap<-select(Ncyc_all,Gene,Black,CK,Green,White)

Ncyc_heatmap <- Ncyc_heatmap[!duplicated(Ncyc_heatmap$Gene),]

rownames(Ncyc_heatmap)<-Ncyc_heatmap[,1]
Ncyc_heatmap<-Ncyc_heatmap[,-1]





library(RColorBrewer)
display.brewer.all()
display.brewer.pal(9,'Set1')
brewer.pal(9,'Set1')

display.brewer.pal(9,'Greys')
brewer.pal(9,'Greys')

display.brewer.pal(9,'Blues')
brewer.pal(9,'Blues')


display.brewer.pal(12,'Set3')
brewer.pal(12,'Set3')


bk = unique(c(seq(-1.8,1.8, length=50)))


            
library(pheatmap)
p6<- pheatmap(Ncyc_heatmap, 
              scale = 'row',
              cluster_cols = T,
              breaks = bk,
              cluster_row = T,
              show_colnames     = T,
              show_rownames     = T,
              color = colorRampPalette(c("#80B1D3", 
                                         "white",
                                         "#FDB462"))(50), 
              legend = T,
              border_color = '#F0F0F0',
              treeheight_row=20,
              treeheight_col=20,
              fontsize = 12.5,
              fontsize_col = 12.5,
              angle_col=45)


#save pdf 10*6 inches






all_9_frames<-Ncyc_all[,1:6]

all_9_frames<-all_9_frames[,c(1,3,2,4:6)]

all_9_frames<-split(all_9_frames,all_9_frames$Pathways)


Anammox_frames<-all_9_frames[which(names(all_9_frames)%in%c("Anammox"))]
Anammox_frames<-as.data.frame(Anammox_frames) 
row.names(Anammox_frames)<-Anammox_frames[,1]
Anammox_frames<-Anammox_frames[,2:5]

all_9_frames[which(names(all_9_frames)%in%c("Anammox"))]<-NULL



library(pheatmap)
pdf('collated_pathways.pdf')

for (i in 1:length(all_9_frames)){
  a<-as.data.frame(all_9_frames[i]) 
  row.names(a)<-a[,1]
  a<-a[,2:5]
  
    theplot<-pheatmap(as.data.frame(a),
             scale = 'row',
             cluster_cols = F,
             breaks = bk,
             cluster_row = T,
             show_colnames     = T,
             show_rownames     = T,
             color = colorRampPalette(c("#80B1D3", 
                                        "white",
                                        "#FDB462"))(50), 
             legend = T,
             border_color = '#F0F0F0',
             treeheight_row=20,
             treeheight_col=20,
             cellwidth = 80,
             cellheight = 14,
             fontsize = 12.5,
             fontsize_col = 12.5,
             angle_col=45)
  print(theplot)
}








p7<-pheatmap(Anammox_frames,
         scale = 'row',
         cluster_cols = F,
         breaks = bk,
         cluster_row = F,
         show_colnames     = T,
         show_rownames     = T,
         color = colorRampPalette(c("#80B1D3", 
                                    "white",
                                    "#FDB462"))(50), 
         legend = T,
         border_color = '#F0F0F0',
         treeheight_row=20,
         treeheight_col=20,
         cellwidth = 80,
         cellheight = 14,
         fontsize = 12.5,
         fontsize_col = 12.5,
         angle_col=45)

p7

dev.off()
