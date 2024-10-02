# test_ggClusterNet
rm(list=ls())

# devtools::install_github("taowenmicro/ggClusterNet")
# install.packages("ggClusterNet-master/", repos = NULL, type = "source")
# BiocManager::install("GO.db")
# BiocManager::install("impute")
# BiocManager::install("preprocessCore")
# BiocManager::install("ggraph")
# BiocManager::install("tidyfst")

library(ggClusterNet)
library(phyloseq)
library(dplyr)
library(WGCNA)
library(sna)
library(igraph)
library(tidyverse)
library(network)
library(ggrepel)
library(ggplot2)


##### 1. import data #####
# import 16S abundance and tax data 
ASV_all_16S_With_Tax<-read.table("ASV_Abundance_16S_split.txt",sep = "\t", header = TRUE,row.names = 1)
ASV_all_16S<-ASV_all_16S_With_Tax[,1:(ncol(ASV_all_16S_With_Tax)-7)]
ASV_all_16S_Tax<-ASV_all_16S_With_Tax[,(ncol(ASV_all_16S_With_Tax)-6):ncol(ASV_all_16S_With_Tax)]

# import ITS abundance and tax data 
ASV_all_ITS_With_Tax<-read.table("ASV_Abundance_ITS_split.txt",sep = "\t", header = TRUE,row.names = 1)
ASV_all_ITS<-ASV_all_ITS_With_Tax[,1:(ncol(ASV_all_ITS_With_Tax)-7)]
ASV_all_ITS_Tax<-ASV_all_ITS_With_Tax[,(ncol(ASV_all_ITS_With_Tax)-6):ncol(ASV_all_16S_With_Tax)]


infor<-data.frame(SampleID=c(rep(NA,27)),Group=c(rep(NA,27)))
infor$SampleID<-colnames(ASV_all_16S)[1:27]
infor$Group<-ifelse(substr(infor$SampleID,1,3)=="Bla","Black",
                     ifelse(substr(infor$SampleID,1,3)=="Whi","White",
                            ifelse(substr(infor$SampleID,1,3)=="Gre","Green","CK")))
row.names(infor)<- infor$SampleID


####### 2. rearrange abundance data #####
ASV_filtered_16S<- ASV_all_16S[rowSums(ASV_all_16S)>ncol(ASV_all_16S),]
ASV_filtered_ITS<- ASV_all_ITS[rowSums(ASV_all_ITS)>ncol(ASV_all_ITS),]

ASV_filtered_16S<-ASV_filtered_16S[,order(colnames(ASV_filtered_16S))]
ASV_filtered_ITS<-ASV_filtered_ITS[,order(colnames(ASV_filtered_ITS))]

colSums(ASV_filtered_16S)
colSums(ASV_filtered_ITS)




####### 3. rearrange tax data #####
common_rows_16S <- intersect(rownames(ASV_all_16S_Tax), rownames(ASV_filtered_16S))
ASV_filtered_16S_Tax<- ASV_all_16S_Tax[common_rows_16S, ]
ASV_filtered_16S_Tax <- ASV_filtered_16S_Tax %>%
  dplyr::rename(Kingdom = Domain)

common_rows_ITS <- intersect(rownames(ASV_all_ITS_Tax), rownames(ASV_filtered_ITS))
ASV_filtered_ITS_Tax<- ASV_all_ITS_Tax[common_rows_ITS, ]







###### 4. import taxonomy info #####

# import data
ps_16S<-phyloseq(otu_table(as.matrix(ASV_filtered_16S),taxa_are_rows = TRUE),
                      tax_table(as.matrix(ASV_filtered_16S_Tax)),
                      sample_data(infor)) 
ps_ITS<-phyloseq(otu_table(as.matrix(ASV_filtered_ITS),taxa_are_rows = TRUE),
                      tax_table(as.matrix(ASV_filtered_ITS_Tax)),
                      sample_data(infor)) 




# 
# ###### 5. merge files ######
# skip...



##### 6. calculate network #####

folder <- "./16S_network"
if (!dir.exists(folder)) {
  dir.create(folder)
}

tab.r_16S <- network.pip(
  ps = ps_16S,
  # ra = 0.05,
  big = TRUE,
  select_layout = FALSE,
  layout_net = "model_maptree2",
  r.threshold = 0.6,
  p.threshold = 0.05,
  maxnode = 2,
  method = "spearman",
  label = TRUE,
  lab = "elements",
  group = "Group",
  fill = "Phylum",
  size = "igraph.degree",
  zipi = TRUE,
  ram.net = TRUE,
  clu_method = "cluster_fast_greedy",
  N = 400, 
  step = 100, 
  R=10, 
  ncpus = 1
)


# save the results
saveRDS(tab.r_16S,paste0(folder,"/","16S_network.pip.rds"))

# Compute the 0/1 matrix and save it
dat_16S <- tab.r_16S [[2]]
cortab_16S <- dat_16S$net.cor.matrix$cortab
saveRDS(cortab_16S,paste0(folder,"/","16S_cor.matrix.all.group.rds"))



#######################################################
#### read rds
tab.r_16S <- readRDS(paste0(folder,"/","16S_network.pip.rds"))
cortab_16S <- readRDS(paste0(folder,"/","16S_cor.matrix.all.group.rds"))







####### 7. plot draft #######

plot <- tab.r_16S[[1]]
p1_tmp <- plot[[1]]
ggsave(paste0(folder,"/","plot.network2.pdf"),p1_tmp,width = 7.5,height = 5)







######## 8. plot #######
dat_16S <- tab.r_16S [[2]]
node_16S <- dat_16S$net.cor.matrix$node
edge_16S <- dat_16S$net.cor.matrix$edge
head(node_16S)
head(edge_16S)



node_plot_16S<- node_16S %>%
  mutate(Phylum_plot = ifelse(Phylum %in% c("Proteobacteria", "Actinobacteriota", "Cyanobacteria", "Chloroflexi",
                                            "Acidobacteriota", "Bacteroidota", "Planctomycetota", "Gemmatimonadota",
                                            "Myxococcota", "Deinococcota"), 
                              Phylum, "Others"))
  


edge_plot_16S<-edge_16S %>% 
  group_by(group) %>%
  mutate(links = n()) %>%
  ungroup() %>%
  mutate(label_plot = paste0(substr(label, 1, nchar(as.character(label)) - 2),
                             " ",
                             links,
                             ")"))

node_plot_16S <- node_plot_16S %>%
  mutate(label_plot = edge_plot_16S$label_plot[match(group, sapply(strsplit(edge_plot_16S$label_plot, ":"), `[`, 1))])



order_of_group<-c("CK", "Black", "Green", "White")
node_plot_16S <- node_plot_16S %>%
  mutate(label_plot = factor(label_plot, levels = unique(label_plot[order(match(sapply(strsplit(label_plot, ":"), `[`, 1), order_of_group))])))
edge_plot_16S <- edge_plot_16S %>%
  mutate(label_plot = factor(label_plot, levels = unique(label_plot[order(match(sapply(strsplit(label_plot, ":"), `[`, 1), order_of_group))])))




Palette_16S<-c("Proteobacteria"="#FB8072", "Actinobacteriota" ="#80B1D3",
           "Cyanobacteria"="#FDB462", "Chloroflexi"="#8DD3C7",
           "Acidobacteriota"="#FFFFB3","Bacteroidota"="#BEBADA",
           "Planctomycetota"="#B3DE69", "Gemmatimonadota"="#FCCDE5",
           "Myxococcota"="#BC80BD","Deinococcota"="#CCEBC5",
           "Others"="grey80")


p1 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = cor
),
data = edge_plot_16S, linewidth= 0.3,alpha = 0.6) +
  geom_point(aes(X1, X2,
                 fill = Phylum_plot,
                 size = igraph.degree),
             pch = 21, data = node_plot_16S,color = "gray40") +
  facet_wrap(.~ label_plot,scales="free_y",nrow = 1) +
  theme(strip.background = element_rect(fill = NA, color = NA))+
  theme(strip.text = element_text(size = 14, hjust = 0.5, face = "bold") ) + 

  scale_colour_manual(values = c("#6D98B5","#D48852"), name = "Correlation") +
  scale_fill_manual(values = Palette_16S, name = "Phylum") +
  
  scale_size(range = c(0.8, 4),name = "Degree") +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5)
  ) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()
  ) +
  theme(legend.position = "bottom",
        legend.key.height = unit(0.01, "cm"))+
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

p1
ggsave(paste0(folder,"/","network_16S.pdf"),p1,width = 15,height = 4)

## Notes: The degree of a vertex is its most basic structural property, the number of its adjacent edges
## from https://igraph.org/r/doc/degree.html










####### 9. plot #######

Samples_Palette <- c("CK"="#EEA236FF",
                     "Black"="#D43F3AFF",
                     "Green"="#5CB85CFF", 
                     "White"="#46B8DAFF")

library(patchwork)
Robustness_list_1<- Robustness.Targeted.removal(ps = ps_16S,
                                  corg = cortab_16S,
                                  degree = TRUE,
                                  zipi = FALSE
)



p2_tmp <- Robustness_list_1[[1]]
p2_tmp

Robustness_dat1 <- Robustness_list_1[[2]]
Robustness_dat1$Group <- factor(Robustness_dat1$Group, levels = c("CK", "Black", "Green", "White"))


# plot
library(ggplot2)

p2 <- ggplot(Robustness_dat1, aes(Number.hub.removed, remain.mean, color = Group)) +
  geom_point(shape=19, alpha=0.6) +  
  geom_errorbar(aes(ymin = remain.mean - remain.sd, ymax = remain.mean + remain.sd), width = 0.2, alpha=0.6) +  # 添加误差条
  geom_smooth(se = FALSE) +
  facet_wrap(~weighted, ncol = 2,
             labeller = labeller(weighted = c("unweighted" = "unweighted", 
                                              "weighted" = " weighted")))+ 
  theme_bw()+
  theme(strip.background = element_rect(fill = NA, color = NA),
        strip.text = element_text(size = 12,face = "bold"),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank())+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))+
  scale_color_manual(values = Samples_Palette)+
  theme(legend.position = "bottom",
        legend.key.height = unit(0.01, "inches"))+
  labs(title = "Targeted removal")+
  theme(plot.title = element_text(size =14, face = "bold", hjust = 0.5,vjust = -0.5))+
  xlab("Number of module hubs removed")+
  ylab("Proportion of species remained (%)")
p2  
ggsave(paste0(folder,"/","Robustness_Targeted removal.pdf"),p2,width = 4.5,height = 4)





####### 10. plot #######

Robustness_list_2<- Robustness.Random.removal(ps = ps_16S,
                                                corg = cortab_16S,
                                              )


p3_tmp <- Robustness_list_2[[1]]
p3_tmp

Robustness_dat2 <- Robustness_list_2[[2]]
Robustness_dat2$Group <- factor(Robustness_dat2$Group, levels = c("CK", "Black", "Green", "White"))



# plot
library(ggplot2)

p3 <- ggplot(Robustness_dat2, aes(Proportion.removed, remain.mean, color = Group)) +
  geom_point(shape=16,alpha=0.6) +  
  geom_errorbar(aes(ymin = remain.mean - remain.sd, ymax = remain.mean + remain.sd), width = 0.02, alpha=0.6) +  # 添加误差条
  geom_smooth(se = FALSE, linewidth = 0.7) +
  facet_wrap(~weighted, ncol = 2,
             labeller = labeller(weighted = c("unweighted" = "unweighted", 
                                              "weighted" = "weighted")))+ 
  theme_bw()+
  theme(strip.background = element_rect(fill = NA, color = NA),
        strip.text = element_text(size = 12,face = "bold"),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank())+
  theme(axis.text.x = element_text(angle = 30,hjust = 1))+
  scale_color_manual(values = Samples_Palette)+
  theme(legend.position = "bottom",
        legend.key.height = unit(0.01, "inches"))+
  labs(title = "Randomly removal")+
  theme(plot.title = element_text(size =14, face = "bold", hjust = 0.5,vjust = -0.5))+
  xlab("Proportion of species removed (%)")+
  ylab("Proportion of species remained (%)")
p3
ggsave(paste0(folder,"/","Robustness_Randomly_removal.pdf"),p3,width = 4.5,height = 4)





####### 11. plot #######


# natural connectivity 
library("pulsar")
Robustness_list_3<- natural.con.microp(ps = ps_16S,
                                              corg = cortab_16S,
                                              # norm = TRUE,
                                              method = "spearman",
                                              end = 280,
                                              start = 0
)


p4_tmp <- Robustness_list_3[[1]]
p4_tmp

Robustness_dat3 <- Robustness_list_3[[2]]
Robustness_dat3$Group <- factor(Robustness_dat3$Group, levels = c("CK", "Black", "Green", "White"))



p4 <- ggplot(Robustness_dat3, aes(Num.of.remove.nodes, Natural.connectivity, color = Group)) +
  geom_point(shape=16,alpha=0.3,size=0.7) +  
  geom_smooth(se = FALSE, linewidth = 0.7) +
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  theme(axis.text.x = element_text(angle = 30,hjust = 1))+
  scale_color_manual(values = Samples_Palette)+
  theme(legend.position = "bottom",
        legend.key.height = unit(0.01, "inches"),
        legend.key.width = unit(0.04, "inches"))+
  labs(title = "Natural connectivity") +
  theme(plot.title = element_text(size =14, face = "bold", hjust = 0.5,vjust = -0.5))+
  xlab("Number of removed nodes")+
  ylab("Natural connectivity")
p4
ggsave(paste0(folder,"/","Robustness_natural_connectivity.pdf"),p4,width = 3,height = 4)








####### 12. plot #######


# negative.correlation.ratio
Robustness_list_4<- negative.correlation.ratio(ps = ps_16S,
                                       corg = cortab_16S,
                                       # norm = TRUE,
                                       method = "spearman",
                                       
)


p5_tmp <- Robustness_list_4[[1]]
p5_tmp

Robustness_dat4 <- Robustness_list_4[[2]]
Robustness_dat4$Group <- factor(Robustness_dat4$ID, levels = c("CK", "Black", "Green", "White"))



p5 <- ggplot(Robustness_dat4, aes(x = Group, y = ratio, color=Group, fill =Group)) +
  geom_bar(stat = "identity",alpha = 0.2) +
  geom_text(aes(label = paste0(sprintf("%.2f", ratio)," %")), vjust = -0.5) +
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  theme(axis.text.x = element_text(angle = 30,hjust = 1))+
  scale_color_manual(values = Samples_Palette)+
  scale_fill_manual(values = Samples_Palette)+
  theme(legend.position = "none")+
  labs(title = "Negative correlation ratio") +
  theme(plot.title = element_text(size =14, face = "bold", hjust = 0.5, vjust = -0.5))+
  ylim(0,30)+
  xlab(NULL)+
  ylab("Negative correlation ratio (%)")
p5
ggsave(paste0(folder,"/","Negative_correlation_ratio.pdf"),p5,width = 3,height = 4)






###### 13 . merge figures #####

library(ggpubr)
p_all <- ggarrange(
  p1,
  ggarrange(p2, p3, p4, p5, ncol = 4, widths = c(4.5, 4.5, 3, 3)),
  ncol = 1,
  heights = c(8, 7) 
)
p_all


ggsave(paste0(folder,"/","Network_all_16S.pdf"),p_all,width = 15,height = 8)


