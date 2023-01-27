##########16S_Phylum

rm(list=ls())
library(RColorBrewer)
display.brewer.all()
display.brewer.pal(12,'Set3')
brewer.pal(12,'Set3')

Palette <- c(brewer.pal(12,'Set3'))
Palette<-Palette[c(1:8,10:11,9)]



abundance <- read.table("16S_otu_table.Phylum.even.txt",header = TRUE, row.names = 1, sep = "\t")
abundance <- as.matrix(abundance)

abundance<-abundance[,-8]

colnames(abundance)<-c("MZ_1","TR_1","MC_1","MZ_2","TR_2","MC_2","MZ_3")

f.abundance <- matrix(as.numeric(abundance),nrow = nrow(abundance))
f.abundance<-t(f.abundance)

f.abundance<-f.abundance/rowSums(f.abundance)

rownames(f.abundance) <- colnames(abundance)
colnames(f.abundance) <- rownames(abundance)
f.abundance <- t(f.abundance)



CT_16S_f.abundance<-f.abundance


sum <- apply(CT_16S_f.abundance,1,sum) 
CT_16S_f.abundance <- cbind(CT_16S_f.abundance,sum)
CT_16S_f.abundance <- as.data.frame(CT_16S_f.abundance)
CT_16S_f.abundance <- CT_16S_f.abundance[order(CT_16S_f.abundance[,"sum"],decreasing = T),]
CT_16S_f.abundance <- subset(CT_16S_f.abundance, select = -sum)

CT_16S_f.abundance <- CT_16S_f.abundance[1:10,]
CT_16S_f.abundance <- t(CT_16S_f.abundance)
sum <- apply(CT_16S_f.abundance,1,sum) 
Others <- 1-sum
CT_16S_f.abundance <- cbind(CT_16S_f.abundance,Others)
CT_16S_f.abundance <- t(CT_16S_f.abundance)

CT_16S_f.abundance<-t(CT_16S_f.abundance)
CT_16S_f.abundance<-as.data.frame(CT_16S_f.abundance)
CT_16S_f.abundance$X<-c(rownames(CT_16S_f.abundance))
CT_16S_f.abundance<-CT_16S_f.abundance[order(CT_16S_f.abundance$X),]
CT_16S_f.abundance<-CT_16S_f.abundance[,-12]
CT_16S_f.abundance<-t(CT_16S_f.abundance)


file_CT_16S_f.abundance<-round(CT_16S_f.abundance*100,digits=2)
write.csv(file_CT_16S_f.abundance, file="CT_16S_Phylum.csv", quote = FALSE)

library(reshape2)
taxon <- melt(CT_16S_f.abundance)
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
#p4 <- p3 
 theme(axis.text.x=element_text(colour="black",size=12,face = "bold",angle = -45,vjust = 0.5,hjust = 0)) + 
 theme(axis.text.y=element_text(colour = "black",size = 18)) + 
  theme(axis.line = element_line(colour = "black"))+
 theme(axis.title.y = element_text(size = 24,face = "bold",margin = unit(c(0,1,0,1),"lines")))+
  scale_y_continuous(expand = c(0,0))
p5 <- p3 + theme(legend.text = element_text(colour = "black",size = 15)) + 
  theme(legend.title = element_text(size = 18,colour = "black")) + 
  labs(title = "Top 10 phyla of prokaryotes")+
  theme(text=element_text(size=20))+
  theme(plot.title = element_text(hjust = 0.5)) 
  #theme(text = element_text(family = "Times"))
p5
#Save_PDF_5*8inches


















########16S_Genus

rm(list=ls())
library(RColorBrewer)
display.brewer.all()
display.brewer.pal(12,'Set3')
display.brewer.pal(8,'Set2')
display.brewer.pal(9,'Set1')
brewer.pal(12,'Set3')
brewer.pal(8,'Set2')
brewer.pal(9,'Set1')

Palette <- c(brewer.pal(12,'Set3'),brewer.pal(8,'Set2'),brewer.pal(9,'Set1'))
Palette<-Palette[c(1:19,21,20)]



abundance <- read.table("16S_otu_table.Genus.even.txt",header = TRUE, row.names = 1, sep = "\t")
abundance <- as.matrix(abundance)

abundance<-abundance[,-8]

colnames(abundance)<-c("MZ_1","TR_1","MC_1","MZ_2","TR_2","MC_2","MZ_3")

f.abundance <- matrix(as.numeric(abundance),nrow = nrow(abundance))
f.abundance<-t(f.abundance)

f.abundance<-f.abundance/rowSums(f.abundance)

rownames(f.abundance) <- colnames(abundance)
colnames(f.abundance) <- rownames(abundance)
f.abundance <- t(f.abundance)



CT_16S_f.abundance<-f.abundance


sum <- apply(CT_16S_f.abundance,1,sum) 
CT_16S_f.abundance <- cbind(CT_16S_f.abundance,sum)
CT_16S_f.abundance <- as.data.frame(CT_16S_f.abundance)
CT_16S_f.abundance <- CT_16S_f.abundance[order(CT_16S_f.abundance[,"sum"],decreasing = T),]
CT_16S_f.abundance <- subset(CT_16S_f.abundance, select = -sum)

CT_16S_f.abundance <- CT_16S_f.abundance[1:20,]
CT_16S_f.abundance <- t(CT_16S_f.abundance)
sum <- apply(CT_16S_f.abundance,1,sum) 
Others <- 1-sum
CT_16S_f.abundance <- cbind(CT_16S_f.abundance,Others)
CT_16S_f.abundance <- t(CT_16S_f.abundance)


CT_16S_f.abundance<-t(CT_16S_f.abundance)
CT_16S_f.abundance<-as.data.frame(CT_16S_f.abundance)
CT_16S_f.abundance$X<-c(rownames(CT_16S_f.abundance))
CT_16S_f.abundance<-CT_16S_f.abundance[order(CT_16S_f.abundance$X),]
CT_16S_f.abundance<-CT_16S_f.abundance[,-22]
CT_16S_f.abundance<-t(CT_16S_f.abundance)


file_CT_16S_f.abundance<-round(CT_16S_f.abundance*100,digits=2)
write.csv(file_CT_16S_f.abundance, file="CT_16S_Genus.csv", quote = FALSE)


library(reshape2)
taxon <- melt(CT_16S_f.abundance)
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
  #p4 <- p3 
  theme(axis.text.x=element_text(colour="black",size=12,face = "bold",angle = -45,vjust = 0.5,hjust = 0)) + 
  theme(axis.text.y=element_text(colour = "black",size = 18)) + 
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.title.y = element_text(size = 24,face = "bold",margin = unit(c(0,1,0,1),"lines")))+
  scale_y_continuous(expand = c(0,0))
p5 <- p3 + theme(legend.text = element_text(colour = "black",size = 15)) + 
  theme(legend.title = element_text(size = 18,colour = "black")) + 
  labs(title = "Top 20 genera of prokaryotes")+
  theme(text=element_text(size=20))+
  theme(plot.title = element_text(hjust = 0.5)) 
#theme(text = element_text(family = "Times"))
p5
#Save_PDF_5*10inches












##########ITS_Phylum

rm(list=ls())
library(RColorBrewer)
display.brewer.all()
display.brewer.pal(12,'Set3')
brewer.pal(12,'Set3')

Palette <- c(brewer.pal(12,'Set3'))
Palette<-Palette[c(1:8,10:11,9)]


abundance <- read.table("ITS_otu_table.Phylum.even.txt",header = TRUE, row.names = 1, sep = "\t")
abundance <- as.matrix(abundance)

abundance<-abundance[,-6]

colnames(abundance)<-c("MZ_1","TR_1","MZ_2","TR_2","MC_2")

f.abundance <- matrix(as.numeric(abundance),nrow = nrow(abundance))
f.abundance<-t(f.abundance)

f.abundance<-f.abundance/rowSums(f.abundance)

rownames(f.abundance) <- colnames(abundance)
colnames(f.abundance) <- rownames(abundance)
f.abundance <- t(f.abundance)



CT_ITS_f.abundance<-f.abundance


sum <- apply(CT_ITS_f.abundance,1,sum) 
CT_ITS_f.abundance <- cbind(CT_ITS_f.abundance,sum)
CT_ITS_f.abundance <- as.data.frame(CT_ITS_f.abundance)
CT_ITS_f.abundance <- CT_ITS_f.abundance[order(CT_ITS_f.abundance[,"sum"],decreasing = T),]
CT_ITS_f.abundance <- subset(CT_ITS_f.abundance, select = -sum)

CT_ITS_f.abundance <- CT_ITS_f.abundance[1:10,]
CT_ITS_f.abundance <- t(CT_ITS_f.abundance)
sum <- apply(CT_ITS_f.abundance,1,sum) 
Others <- 1-sum
CT_ITS_f.abundance <- cbind(CT_ITS_f.abundance,Others)
CT_ITS_f.abundance <- t(CT_ITS_f.abundance)



CT_ITS_f.abundance<-t(CT_ITS_f.abundance)
CT_ITS_f.abundance<-as.data.frame(CT_ITS_f.abundance)
CT_ITS_f.abundance$X<-c(rownames(CT_ITS_f.abundance))
CT_ITS_f.abundance<-CT_ITS_f.abundance[order(CT_ITS_f.abundance$X),]
CT_ITS_f.abundance<-CT_ITS_f.abundance[,-12]
CT_ITS_f.abundance<-t(CT_ITS_f.abundance)


file_CT_ITS_f.abundance<-round(CT_ITS_f.abundance*100,digits=2)
write.csv(file_CT_ITS_f.abundance, file="CT_ITS_Phylum.csv", quote = FALSE)


library(reshape2)
taxon <- melt(CT_ITS_f.abundance)
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
  #p4 <- p3 
  theme(axis.text.x=element_text(colour="black",size=12,face = "bold",angle = -45,vjust = 0.5,hjust = 0)) + 
  theme(axis.text.y=element_text(colour = "black",size = 18)) + 
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.title.y = element_text(size = 24,face = "bold",margin = unit(c(0,1,0,1),"lines")))+
  scale_y_continuous(expand = c(0,0))
p5 <- p3 + theme(legend.text = element_text(colour = "black",size = 15)) + 
  theme(legend.title = element_text(size = 18,colour = "black")) + 
  labs(title = "Top 10 phyla of eukaryotes")+
  theme(text=element_text(size=20))+
  theme(plot.title = element_text(hjust = 0.5)) 
#theme(text = element_text(family = "Times"))
p5
#Save_PDF_5*9inches



















########ITS_Genus

rm(list=ls())
library(RColorBrewer)
display.brewer.all()
display.brewer.pal(12,'Set3')
display.brewer.pal(8,'Set2')
display.brewer.pal(9,'Set1')
brewer.pal(12,'Set3')
brewer.pal(8,'Set2')
brewer.pal(9,'Set1')

Palette <- c(brewer.pal(12,'Set3'),brewer.pal(8,'Set2'),brewer.pal(9,'Set1'))
Palette<-Palette[c(1:19,21,20)]



abundance <- read.table("ITS_otu_table.Genus.even.txt",header = TRUE, row.names = 1, sep = "\t")
abundance <- as.matrix(abundance)

abundance<-abundance[,-6]

colnames(abundance)<-c("MZ_1","TR_1","MZ_2","TR_2","MC_2")

f.abundance <- matrix(as.numeric(abundance),nrow = nrow(abundance))
f.abundance<-t(f.abundance)

f.abundance<-f.abundance/rowSums(f.abundance)

rownames(f.abundance) <- colnames(abundance)
colnames(f.abundance) <- rownames(abundance)
f.abundance <- t(f.abundance)



CT_ITS_f.abundance<-f.abundance


sum <- apply(CT_ITS_f.abundance,1,sum) 
CT_ITS_f.abundance <- cbind(CT_ITS_f.abundance,sum)
CT_ITS_f.abundance <- as.data.frame(CT_ITS_f.abundance)
CT_ITS_f.abundance <- CT_ITS_f.abundance[order(CT_ITS_f.abundance[,"sum"],decreasing = T),]
CT_ITS_f.abundance <- subset(CT_ITS_f.abundance, select = -sum)

CT_ITS_f.abundance <- CT_ITS_f.abundance[1:20,]
CT_ITS_f.abundance <- t(CT_ITS_f.abundance)
sum <- apply(CT_ITS_f.abundance,1,sum) 
Others <- 1-sum
CT_ITS_f.abundance <- cbind(CT_ITS_f.abundance,Others)
CT_ITS_f.abundance <- t(CT_ITS_f.abundance)


CT_ITS_f.abundance<-t(CT_ITS_f.abundance)
CT_ITS_f.abundance<-as.data.frame(CT_ITS_f.abundance)
CT_ITS_f.abundance$X<-c(rownames(CT_ITS_f.abundance))
CT_ITS_f.abundance<-CT_ITS_f.abundance[order(CT_ITS_f.abundance$X),]
CT_ITS_f.abundance<-CT_ITS_f.abundance[,-22]
CT_ITS_f.abundance<-t(CT_ITS_f.abundance)


file_CT_ITS_f.abundance<-round(CT_ITS_f.abundance*100,digits=2)
write.csv(file_CT_ITS_f.abundance, file="CT_ITS_Genus.csv", quote = FALSE)


library(reshape2)
taxon <- melt(CT_ITS_f.abundance)
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
  #p4 <- p3 
  theme(axis.text.x=element_text(colour="black",size=12,face = "bold",angle = -45,vjust = 0.5,hjust = 0)) + 
  theme(axis.text.y=element_text(colour = "black",size = 18)) + 
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.title.y = element_text(size = 24,face = "bold",margin = unit(c(0,1,0,1),"lines")))+
  scale_y_continuous(expand = c(0,0))
p5 <- p3 + theme(legend.text = element_text(colour = "black",size = 15)) + 
  theme(legend.title = element_text(size = 18,colour = "black")) + 
  labs(title = "Top 20 genera of eukaryotes")+
  theme(text=element_text(size=20))+
  theme(plot.title = element_text(hjust = 0.5)) 
#theme(text = element_text(family = "Times"))
p5
#Save_PDF_5*9inches




