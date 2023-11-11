##########ITS Phyla

rm(list=ls())
library(RColorBrewer)
display.brewer.all()
display.brewer.pal(12,'Set3')
brewer.pal(12,'Set3')

Palette <- c(brewer.pal(12,'Set3'))
Palette<-Palette[c(1:8,10:11,9)]
Palette<-Palette[c(4:6,1:3,7:11)]


abundance <- read.table("Taxonomy_ITS_2Phylum.txt",header = TRUE, sep = "\t")


library(tidyverse)
abundance_tmp<-separate(data=abundance, col=Phylum, 
                              into = c("Kingdom","Phylum"),
                              sep = "\\|")
abundance_tmp<-abundance_tmp[,-1]


unassigned_tmp<-abundance_tmp[abundance_tmp$Phylum=="unassigned",]
unassigned<-apply(unassigned_tmp[,2:length(colnames(unassigned_tmp))],2,sum)
unassigned<-as.data.frame(unassigned)
unassigned<-t(unassigned)
unassigned<-as.data.frame(unassigned)
unassigned$Phylum<-c("unassigned")

uncultured_tmp<-abundance_tmp[abundance_tmp$Phylum=="uncultured",]
uncultured<-apply(uncultured_tmp[,2:length(colnames(uncultured_tmp))],2,sum)
uncultured<-as.data.frame(uncultured)
uncultured<-t(uncultured)
uncultured<-as.data.frame(uncultured)
uncultured$Phylum<-c("uncultured")
uncultured<-as.data.frame(uncultured)


abundance_tmp_1<-filter(abundance_tmp,abundance_tmp$Phylum!="unassigned") 
abundance_tmp_2<-filter(abundance_tmp_1,abundance_tmp_1$Phylum!="uncultured") 
abundance_tmp<-rbind(abundance_tmp_2,uncultured,unassigned)


abundance<-abundance_tmp
rm(abundance_tmp,abundance_tmp_1,abundance_tmp_2,
   unassigned,unassigned_tmp,uncultured,uncultured_tmp )

f.abundance<-abundance[,-1]
rownames(f.abundance)<-abundance$Phylum

f.abundance<-f.abundance/100

Leshan_ITS_abundance<-f.abundance


sum <- apply(Leshan_ITS_abundance,1,sum) 
Leshan_ITS_abundance <- cbind(Leshan_ITS_abundance,sum)
Leshan_ITS_abundance <- as.data.frame(Leshan_ITS_abundance)
Leshan_ITS_abundance <- Leshan_ITS_abundance[order(Leshan_ITS_abundance[,"sum"],decreasing = T),]
Leshan_ITS_abundance <- subset(Leshan_ITS_abundance, select = -sum)

Leshan_ITS_abundance <- Leshan_ITS_abundance[1:7,]
Leshan_ITS_abundance <- t(Leshan_ITS_abundance)
sum <- apply(Leshan_ITS_abundance,1,sum) 
Others <- 1-sum
Leshan_ITS_abundance <- cbind(Leshan_ITS_abundance,Others)
Leshan_ITS_abundance <- t(Leshan_ITS_abundance)

Leshan_ITS_abundance<-t(Leshan_ITS_abundance)
Leshan_ITS_abundance<-as.data.frame(Leshan_ITS_abundance)
Leshan_ITS_abundance$X<-c(rownames(Leshan_ITS_abundance))
Leshan_ITS_abundance<-Leshan_ITS_abundance[order(Leshan_ITS_abundance$X),]
Leshan_ITS_abundance<-subset(Leshan_ITS_abundance, select = -X)
Leshan_ITS_abundance<-t(Leshan_ITS_abundance)

Leshan_ITS_abundance<-Leshan_ITS_abundance[,c(7:12,1:6,13:length(colnames(Leshan_ITS_abundance)))]

file_Leshan_ITS_abundance<-round(Leshan_ITS_abundance*100,digits=2)
write.csv(file_Leshan_ITS_abundance, file="Leshan_ITS_Phyla.csv", quote = FALSE)

library(reshape2)
taxon <- melt(Leshan_ITS_abundance)
colnames(taxon) <- c("Taxon","variable","value")
library(ggplot2)
library(ggalluvial)

p <- ggplot(data = taxon,aes(x = variable, y = value, alluvium = Taxon, stratum = Taxon))
p1 <- p + geom_alluvium(aes(fill = Taxon),alpha = .5,width = 0.6) + 
  geom_stratum(aes(fill = Taxon),width = 0.6)
p2 <- p1 + ylab(label = "Relative Abundance") + xlab(label = "")
p3 <- p2 + scale_fill_manual(values = Palette) +
  theme_bw()+ theme(panel.grid=element_blank()) + 
  theme(panel.border = element_blank()) +
  theme(panel.background=element_rect(fill='transparent', color='black'),plot.margin = unit(c(3,5,1,1),"mm")) +

 theme(axis.text.x=element_text(colour="black",size=12,face = "bold",angle = -45,vjust = 0.5,hjust = 0)) + 
 theme(axis.text.y=element_text(colour = "black",size = 18)) + 
  theme(axis.line = element_line(colour = "black"))+
 theme(axis.title.y = element_text(size = 18,face = "bold",margin = unit(c(0,1,0,1),"lines")))+
  scale_y_continuous(expand = c(0,0))
p5 <- p3 + theme(legend.text = element_text(colour = "black",size = 15)) + 
  theme(legend.title = element_text(size = 18,colour = "black")) + 
  labs(title = "Top 7 Phyla of Eukaryotes")+
  theme(text=element_text(size=20))+
  theme(plot.title = element_text(hjust = 0.5)) 
  
p5
#export PDF_4*12inches

















##########ITS Genera

rm(list=ls())
library(RColorBrewer)
display.brewer.all()
display.brewer.pal(12,'Set3')
display.brewer.pal(9,'Pastel1')
brewer.pal(12,'Set3')
brewer.pal(9,'Pastel1')

Palette_Set3 <- c(brewer.pal(12,'Set3'))
Palette_Pastel1<-c(brewer.pal(8,'Pastel1'))
Palette<-c(Palette_Set3,Palette_Pastel1,"lightgray")



abundance <- read.table("Taxonomy_ITS_6Genus.txt",header = TRUE, sep = "\t")


library(tidyverse)
abundance_tmp<-abundance
abundance_tmp<-separate(data=abundance, col=Genus, 
                        into = c("Kingdom","Phylum", "Class", "Order", "Family", "Genus"),
                        sep = "\\|")
abundance_tmp<-abundance_tmp[,-c(1:5)]


unassigned_tmp<-abundance_tmp[abundance_tmp$Genus=="unassigned",]
unassigned<-apply(unassigned_tmp[,2:length(colnames(unassigned_tmp))],2,sum)
unassigned<-as.data.frame(unassigned)
unassigned<-t(unassigned)
unassigned<-as.data.frame(unassigned)
unassigned$Genus<-c("unassigned")

uncultured_tmp<-abundance_tmp[abundance_tmp$Genus=="uncultured",]
uncultured<-apply(uncultured_tmp[,2:length(colnames(uncultured_tmp))],2,sum)
uncultured<-as.data.frame(uncultured)
uncultured<-t(uncultured)
uncultured<-as.data.frame(uncultured)
uncultured$Genus<-c("uncultured")
uncultured<-as.data.frame(uncultured)


abundance_tmp_1<-filter(abundance_tmp,abundance_tmp$Genus!="unassigned") 
abundance_tmp_2<-filter(abundance_tmp_1,abundance_tmp_1$Genus!="uncultured") 
abundance_tmp<-rbind(abundance_tmp_2,uncultured,unassigned)


abundance<-abundance_tmp
rm(abundance_tmp,abundance_tmp_1,abundance_tmp_2,
   unassigned,unassigned_tmp,uncultured,uncultured_tmp )

f.abundance<-abundance[,-1]
rownames(f.abundance)<-abundance$Genus

f.abundance<-f.abundance/100

Leshan_ITS_abundance<-f.abundance


sum <- apply(Leshan_ITS_abundance,1,sum) 
Leshan_ITS_abundance <- cbind(Leshan_ITS_abundance,sum)
Leshan_ITS_abundance <- as.data.frame(Leshan_ITS_abundance)
Leshan_ITS_abundance <- Leshan_ITS_abundance[order(Leshan_ITS_abundance[,"sum"],decreasing = T),]
Leshan_ITS_abundance <- subset(Leshan_ITS_abundance, select = -sum)

Leshan_ITS_abundance <- Leshan_ITS_abundance[1:20,]
Leshan_ITS_abundance <- t(Leshan_ITS_abundance)
sum <- apply(Leshan_ITS_abundance,1,sum) 
Others <- 1-sum
Leshan_ITS_abundance <- cbind(Leshan_ITS_abundance,Others)
Leshan_ITS_abundance <- t(Leshan_ITS_abundance)

Leshan_ITS_abundance<-t(Leshan_ITS_abundance)
Leshan_ITS_abundance<-as.data.frame(Leshan_ITS_abundance)
Leshan_ITS_abundance$X<-c(rownames(Leshan_ITS_abundance))
Leshan_ITS_abundance<-Leshan_ITS_abundance[order(Leshan_ITS_abundance$X),]
Leshan_ITS_abundance<-subset(Leshan_ITS_abundance, select = -X)
Leshan_ITS_abundance<-t(Leshan_ITS_abundance)

Leshan_ITS_abundance<-Leshan_ITS_abundance[,c(7:12,1:6,13:length(colnames(Leshan_ITS_abundance)))]

file_Leshan_ITS_abundance<-round(Leshan_ITS_abundance*100,digits=2)
write.csv(file_Leshan_ITS_abundance, file="Leshan_ITS_Genus.csv", quote = FALSE)

library(reshape2)
taxon <- melt(Leshan_ITS_abundance)
colnames(taxon) <- c("Taxon","variable","value")
library(ggplot2)
library(ggalluvial)

p <- ggplot(data = taxon,aes(x = variable, y = value, alluvium = Taxon, stratum = Taxon))
p1 <- p + geom_alluvium(aes(fill = Taxon),alpha = .5,width = 0.6) + 
  geom_stratum(aes(fill = Taxon),width = 0.6)
p2 <- p1 + ylab(label = "Relative Abundance") + xlab(label = "")
p3 <- p2 + scale_fill_manual(values = Palette) +
  theme_bw()+ theme(panel.grid=element_blank()) + 
  theme(panel.border = element_blank()) +
  theme(panel.background=element_rect(fill='transparent', color='black'),plot.margin = unit(c(3,5,1,1),"mm")) +
  
  theme(axis.text.x=element_text(colour="black",size=12,face = "bold",angle = -45,vjust = 0.5,hjust = 0)) + 
  theme(axis.text.y=element_text(colour = "black",size = 18)) + 
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.title.y = element_text(size = 18,face = "bold",margin = unit(c(0,1,0,1),"lines")))+
  scale_y_continuous(expand = c(0,0))
p5 <- p3 + theme(legend.text = element_text(colour = "black",size = 10)) + 
  theme(legend.title = element_text(size = 12,colour = "black")) + 
  labs(title = "Top 20 Genus of Eukaryotes")+
  theme(text=element_text(size=20))+
  theme(plot.title = element_text(hjust = 0.5)) 

p5
#export 4*12inches


















