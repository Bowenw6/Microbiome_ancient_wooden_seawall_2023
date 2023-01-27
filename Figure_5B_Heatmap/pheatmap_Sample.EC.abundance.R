rm(list=ls())
library(RColorBrewer)
display.brewer.all()
brewer.pal(9,'YlGnBu')
brewer.pal(12,'Paired')
brewer.pal(9,'Set1')
brewer.pal(9,'Pastel1')
display.brewer.pal(9,'Blues')
display.brewer.pal(9,'YlGnBu')

pallette<- brewer.pal(8,'YlGnBu')


sample_EC <- read.table('Sample.EC.abundance_16S.txt', sep = '\t', row.names = 1, header = T, quote = "")

colnames(sample_EC)<-c("MZ_1","TR_1","MC_1","MZ_2","TR_2","MC_2","MZ_3")
sample_EC<-sample_EC[,1:7]
sample_EC_heatmap<-sample_EC[which(rownames(sample_EC)%in%c("EC:3.2.1.4","EC:3.2.1.21","EC:3.2.1.74",
                                                            "EC:3.2.1.91","EC:2.4.1.20","EC:3.2.1.86")),]

#There is no EC:3.2.1.74





#######plot


sample_EC_heatmap<-log10(sample_EC_heatmap)

sample_EC_heatmap<-t(sample_EC_heatmap)


library(pheatmap)
bk1 = unique(c(seq(0,5, length=50)))
p1<- pheatmap(sample_EC_heatmap, 
              scale = 'none',
              cluster_cols = T,
              breaks = bk1,
              display_numbers = T,
              number_color="#FF7F00",
              fontsize_number=9,
              cluster_row = T,
              show_colnames     = T,
              show_rownames     = T,
              color = colorRampPalette(pallette)(50), 
              legend = T,
              border_color = 'lightgray',
              treeheight_row=20,
              treeheight_col=20,
              fontsize = 10,
              fontsize_col = 10,
              angle_col=45)


#save pdf 4*5 inches








