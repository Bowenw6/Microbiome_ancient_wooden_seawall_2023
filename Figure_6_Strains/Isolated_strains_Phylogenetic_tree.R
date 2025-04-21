
#in MacOS command line, perform alignment with muscle
#muscle -in isolated_strains_tree.fasta -out aligned_isolated_strains_tree.fasta 
# in server, construct tree with iqtree2
#nohup iqtree -s aligned_isolated_strains_tree.fasta -m MFP -B 1000 --bnni -T AUTO > aligned_isolated_strains_tree_nohup.out & 

rm(list=ls())

library(RColorBrewer)
display.brewer.all()
display.brewer.pal(9,'Set1')
brewer.pal(9,'Set1')[c(2,3,4)]

treepallette<-brewer.pal(9,'Set1')[c(2,3,4)]



library(ggtree)
library(ggplot2)
library(ggtreeExtra)
library(treeio)
tree_16S<-read.tree("aligned_isolated_strains_tree.fasta.treefile")
map_16S<- read.table("Taxonomy_16S_tree.txt", header = TRUE, sep = "\t")


p<-ggtree(tree_16S, layout="fan", size=0.8,open.angle=165) %<+% map_16S + 
  geom_tippoint(aes(color=Phylum, shape=Sample), size=3, alpha =0.9)+

  scale_colour_manual(values = treepallette)+
  scale_shape_manual(values=c("MZ_2" = 1, 
                              "MZ_3" = 16,
                              "MC_1" = 2,
                              "MC_2" = 17,
                              "TR_1" = 5,
                              "TR_2" = 18))+
  geom_tiplab(aes(color=Phylum),align=T, linetype=3, linesize=0.5,size=2,hjust=-0.07)+
  theme(legend.title=element_text(face="bold"), legend.position="bottom", legend.box="horizontal", legend.text=element_text(size=rel(0.8)))+
  geom_text2(aes(subset=!isTip, label=label), hjust=-0.05, size=2, color="red") + 
  geom_nodepoint(color="orange", alpha=1/2, size=1.5)+
  xlim(0,0.8) 

p<-rotate_tree(p,-16)
p


#save image 9 * 16 inches





