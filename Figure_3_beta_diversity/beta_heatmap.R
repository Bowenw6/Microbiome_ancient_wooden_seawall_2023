rm(list=ls())

library(RColorBrewer)
display.brewer.all()
brewer.pal(9,'YlOrRd')
display.brewer.pal(9,'YlOrRd')
brewer.pal(9,'RdYlBu')
display.brewer.pal(9,'RdYlBu')





data_16S <- read.table('16S_feature_table.taxonomy.even.txt', sep = '\t', row.names = 1, header = T)

data_16S<- data_16S[,-8]

colnames(data_16S)<-c("MZ_1","TR_1","MC_1","MZ_2","TR_2","MC_2","MZ_3")

data_16S<-data_16S[which(rowSums(data_16S)>=7),]

data_16S<-t(data_16S)

#Bray_Curtis
library(vegan)
CT_bray_16S <- vegdist(data_16S,method = "bray")
CT_bray_16S<-as.matrix(CT_bray_16S)







library(pheatmap)

bk1 = unique(c(seq(0,1, length=50)))
p1<- pheatmap(CT_bray_16S, 
              scale = 'none',
              cluster_cols = T,
              breaks = bk1,
              cluster_row = T,
              show_colnames     = T,
              show_rownames     = T,
              color = colorRampPalette(rev(brewer.pal(9,'RdYlBu')
                                       ))(50), 
              legend = T,
              border_color = NA,
              treeheight_row=30,
              treeheight_col=30,
              fontsize = 12.5,
              fontsize_col = 12.5,
              na_col = 'white',
              angle_col=45)

p1
#save pdf 3*4 inches






#Unifrac




x<- scan('dna_sequences_16S.fasta',what = "c")
x[grep('>',x)]<-paste0('&',x[grep('>',x)],'@')
x<-unlist(strsplit(unlist(strsplit(paste(x,collapse =''),"&")),"@"))
write.table(x,'dna_sequences_16S_revised.fasta',sep = '\n',row.names = F,col.names =FALSE,quote =FALSE)


Fasta_16S<-read.table("dna_sequences_16S_revised.fasta")
re_Fasta_16S<-Fasta_16S[seq(1,nrow(Fasta_16S),2),]
re_Fasta_16S<-as.data.frame(re_Fasta_16S)
re_Fasta_16S$V2<-Fasta_16S[seq(0,nrow(Fasta_16S),2),]
colnames(re_Fasta_16S)<- c("OTUs_names","Sequences")
re_Fasta_16S$OTUs_names<- gsub(">*","",re_Fasta_16S$OTUs_names) 
row.names(re_Fasta_16S)<- re_Fasta_16S$OTUs_names


data_16S<-t(data_16S)
data_16S<-as.data.frame(data_16S)

OTU_filtered_fasta_16S<-re_Fasta_16S[match(row.names(data_16S), row.names(re_Fasta_16S)),]

print(row.names(OTU_filtered_fasta_16S)==row.names(data_16S))

OTU_filtered_fasta_16S$OTUs_names<- paste(">",OTU_filtered_fasta_16S$OTUs_names, sep = "")

write.table(OTU_filtered_fasta_16S, "OTU_filtered_fasta_16S.fasta",sep =",", row.names =FALSE, col.names =FALSE, quote =FALSE)

filtered_fasta_16S<-readLines("OTU_filtered_fasta_16S.fasta",encoding = "UTF-8")

filtered_fasta_16S<-gsub(",","\n",filtered_fasta_16S)

writeLines(filtered_fasta_16S,"filtered_fasta_16S.fasta")




#in MacOS command line, perform alignment with muscle
#nohup muscle -in filtered_fasta_16S.fasta -out aligned_filtered_fasta_16S.fasta &
# in MacOS command line, construct tree with iqtree
#nohup iqtree -s aligned_filtered_fasta_16S.fasta -m MFP -B 1000 --bnni -cptime 6000 -T AUTO &




data_16S<-t(data_16S)
data_16S<-as.data.frame(data_16S)

library(phyloseq)
library(GUniFrac)
library(treeio)


#
tree <- read.tree("aligned_filtered_fasta_16S.fasta.treefile")
unifrac <- phyloseq(otu_table(data_16S,taxa_are_rows = F),phy_tree(tree))
#W_Unifrac
W.unifrac <- distance(unifrac,method = "wunifrac")
W.unifrac <- as.matrix(W.unifrac)
#unW_Unifrac
U.unifrac <- distance(unifrac,method = "unifrac")
U.unifrac <- as.matrix(U.unifrac)


#W_Unifrac_16S

library(pheatmap)

bk2 = unique(c(seq(0,1, length=50)))
p2<- pheatmap(W.unifrac, 
              scale = 'none',
              cluster_cols = T,
              breaks = bk2,
              cluster_row = T,
              show_colnames     = T,
              show_rownames     = T,
              color = colorRampPalette(rev(brewer.pal(9,'RdYlBu')
              ))(50), 
              legend = T,
              border_color = NA,
              treeheight_row=30,
              treeheight_col=30,
              fontsize = 12.5,
              fontsize_col = 12.5,
              na_col = 'white',
              angle_col=45)

p2
#save pdf 3*4 inches





#unW_Unifrac_16S

library(pheatmap)

bk3 = unique(c(seq(0,1, length=50)))
p3<- pheatmap(U.unifrac, 
              scale = 'none',
              cluster_cols = T,
              breaks = bk3,
              cluster_row = T,
              show_colnames     = T,
              show_rownames     = T,
              color = colorRampPalette(rev(brewer.pal(9,'RdYlBu')
              ))(50), 
              legend = T,
              border_color = NA,
              treeheight_row=30,
              treeheight_col=30,
              fontsize = 12.5,
              fontsize_col = 12.5,
              na_col = 'white',#设置NA值的颜色
              angle_col=45)

p3
#save pdf 3*4 inches


















rm(list = ls())
#ITS

data_ITS <- read.table('ITS_feature_table.taxonomy.even.txt', sep = '\t', row.names = 1, header = T)

data_ITS<- data_ITS[,-6]

colnames(data_ITS)<-c("MZ_1","TR_1","MZ_2","TR_2","MC_2")

data_ITS<-data_ITS[which(rowSums(data_ITS)>=5),]

data_ITS<-t(data_ITS)

#Bray_Curtis
library(vegan)
CT_bray_ITS <- vegdist(data_ITS,method = "bray")
CT_bray_ITS<-as.matrix(CT_bray_ITS)


#bray_ITS

library(pheatmap)

bk4 = unique(c(seq(0,1, length=50)))
p4<- pheatmap(CT_bray_ITS, 
              scale = 'none',
              cluster_cols = T,
              breaks = bk4,
              cluster_row = T,
              show_colnames     = T,
              show_rownames     = T,
              color = colorRampPalette(rev(brewer.pal(9,'RdYlBu')
              ))(50), 
              legend = T,
              border_color = NA,
              treeheight_row=30,
              treeheight_col=30,
              fontsize = 12.5,
              fontsize_col = 12.5,
              na_col = 'white',
              angle_col=45)

p4
#save pdf 3*4 inches






#ITS_Unifrac

x<- scan('dna_sequences_ITS.fasta',what = "c")
x[grep('>',x)]<-paste0('&',x[grep('>',x)],'@')
x<-unlist(strsplit(unlist(strsplit(paste(x,collapse =''),"&")),"@"))
write.table(x,'dna_sequences_ITS_revised.fasta',sep = '\n',row.names = F,col.names =FALSE,quote =FALSE)


Fasta_ITS<-read.table("dna_sequences_ITS_revised.fasta")
re_Fasta_ITS<-Fasta_ITS[seq(1,nrow(Fasta_ITS),2),]#提取奇数行放在新表格的V1
re_Fasta_ITS<-as.data.frame(re_Fasta_ITS)
re_Fasta_ITS$V2<-Fasta_ITS[seq(0,nrow(Fasta_ITS),2),]#提取偶数行放在新表格的V2
colnames(re_Fasta_ITS)<- c("OTUs_names","Sequences")#改个列名
re_Fasta_ITS$OTUs_names<- gsub(">*","",re_Fasta_ITS$OTUs_names) #把OTU编号前的>去掉
row.names(re_Fasta_ITS)<- re_Fasta_ITS$OTUs_names

data_ITS<-t(data_ITS)
data_ITS<-as.data.frame(data_ITS)
OTU_filtered_fasta_ITS<-re_Fasta_ITS[match(row.names(data_ITS), row.names(re_Fasta_ITS)),]
print(row.names(OTU_filtered_fasta_ITS)==row.names(data_ITS))
OTU_filtered_fasta_ITS$OTUs_names<- paste(">",OTU_filtered_fasta_ITS$OTUs_names, sep = "")
write.table(OTU_filtered_fasta_ITS, "OTU_filtered_fasta_ITS.fasta",sep =",", row.names =FALSE, col.names =FALSE, quote =FALSE)
filtered_fasta_ITS<-readLines("OTU_filtered_fasta_ITS.fasta",encoding = "UTF-8")
filtered_fasta_ITS<-gsub(",","\n",filtered_fasta_ITS)
writeLines(filtered_fasta_ITS,"filtered_fasta_ITS.fasta")



#in MacOS command line, perform alignment with muscle
#nohup muscle -in filtered_fasta_ITS.fasta -out aligned_filtered_fasta_ITS.fasta &
# in MacOS command line, construct tree with iqtree
#nohup iqtree -s aligned_filtered_fasta_ITS.fasta -m MFP -B 1000 --bnni -cptime 6000 -T AUTO &






data_ITS<-t(data_ITS)
data_ITS<-as.data.frame(data_ITS)

library(phyloseq)
library(GUniFrac)
library(treeio)


tree <- read.tree("aligned_filtered_fasta_ITS.fasta.treefile")
unifrac <- phyloseq(otu_table(data_ITS,taxa_are_rows = F),phy_tree(tree))
#W_Unifrac
W.unifrac <- distance(unifrac,method = "wunifrac")
W.unifrac <- as.matrix(W.unifrac)
#unW_Unifrac
U.unifrac <- distance(unifrac,method = "unifrac")
U.unifrac <- as.matrix(U.unifrac)


#W_Unifrac_ITS

library(pheatmap)
bk5 = unique(c(seq(0,1, length=50)))
p5<- pheatmap(W.unifrac, 
              scale = 'none',
              cluster_cols = T,
              breaks = bk5,
              cluster_row = T,
              show_colnames     = T,
              show_rownames     = T,
              color = colorRampPalette(rev(brewer.pal(9,'RdYlBu')
              ))(50), 
              legend = T,
              border_color = NA,
              treeheight_row=30,
              treeheight_col=30,
              fontsize = 12.5,
              fontsize_col = 12.5,
              na_col = 'white',
              angle_col=45)

p5
#save pdf 3*4 inches





#unW_Unifrac_ITS

library(pheatmap)
bk6 = unique(c(seq(0,1, length=50)))
p6<- pheatmap(U.unifrac, 
              scale = 'none',
              cluster_cols = T,
              breaks = bk6,
              cluster_row = T,
              show_colnames     = T,
              show_rownames     = T,
              color = colorRampPalette(rev(brewer.pal(9,'RdYlBu')
              ))(50), 
              legend = T,
              border_color = NA,
              treeheight_row=30,
              treeheight_col=30,
              fontsize = 12.5,
              fontsize_col = 12.5,
              na_col = 'white',
              angle_col=45)

p6
#save pdf 3*4 inches









