
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
# import ITS abundance and tax data 
ASV_all_ITS_With_Tax<-read.table("ASV_Abundance_ITS_split.txt",sep = "\t", header = TRUE,row.names = 1)
ASV_all_ITS<-ASV_all_ITS_With_Tax[,1:(ncol(ASV_all_ITS_With_Tax)-7)]
ASV_all_ITS_Tax<-ASV_all_ITS_With_Tax[,(ncol(ASV_all_ITS_With_Tax)-6):ncol(ASV_all_ITS_With_Tax)]

# import ITS abundance and tax data 
ASV_all_ITS_With_Tax<-read.table("ASV_Abundance_ITS_split.txt",sep = "\t", header = TRUE,row.names = 1)
ASV_all_ITS<-ASV_all_ITS_With_Tax[,1:(ncol(ASV_all_ITS_With_Tax)-7)]
ASV_all_ITS_Tax<-ASV_all_ITS_With_Tax[,(ncol(ASV_all_ITS_With_Tax)-6):ncol(ASV_all_ITS_With_Tax)]


infor<-data.frame(SampleID=c(rep(NA,27)),Group=c(rep(NA,27)))
infor$SampleID<-colnames(ASV_all_ITS)[1:27]
infor$Group<-ifelse(substr(infor$SampleID,1,3)=="Bla","Black",
                     ifelse(substr(infor$SampleID,1,3)=="Whi","White",
                            ifelse(substr(infor$SampleID,1,3)=="Gre","Green","CK")))
row.names(infor)<- infor$SampleID


####### 2. rearrange abundance data #####
ASV_filtered_ITS<- ASV_all_ITS[rowSums(ASV_all_ITS)>ncol(ASV_all_ITS),]

ASV_filtered_ITS<-ASV_filtered_ITS[,order(colnames(ASV_filtered_ITS))]

colSums(ASV_filtered_ITS)




####### 3. rearrange tax data #####
common_rows_ITS <- intersect(rownames(ASV_all_ITS_Tax), rownames(ASV_filtered_ITS))
ASV_filtered_ITS_Tax<- ASV_all_ITS_Tax[common_rows_ITS, ]

common_rows_ITS <- intersect(rownames(ASV_all_ITS_Tax), rownames(ASV_filtered_ITS))
ASV_filtered_ITS_Tax<- ASV_all_ITS_Tax[common_rows_ITS, ]




###### 4. import taxonomy info #####
## skip...

# import data
ps_ITS<-phyloseq(otu_table(as.matrix(ASV_filtered_ITS),taxa_are_rows = TRUE),
                      tax_table(as.matrix(ASV_filtered_ITS_Tax)),
                      sample_data(infor)) 





# ###### 5. merge files ######
# skip...



##### 6. calculate network #####

folder <- "./ITS_network"
if (!dir.exists(folder)) {
  dir.create(folder)
}

tab.r_ITS <- network.pip(
  ps = ps_ITS,
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
saveRDS(tab.r_ITS,paste0(folder,"/","ITS_network.pip.rds"))

# Compute the 0/1 matrix and save it
dat_ITS <- tab.r_ITS [[2]]
cortab_ITS <- dat_ITS$net.cor.matrix$cortab
saveRDS(cortab_ITS,paste0(folder,"/","ITS_cor.matrix.all.group.rds"))



#######################################################
#### read rds
tab.r_ITS <- readRDS(paste0(folder,"/","ITS_network.pip.rds"))
cortab_ITS <- readRDS(paste0(folder,"/","ITS_cor.matrix.all.group.rds"))







####### 7. plot draft #######
plot <- tab.r_ITS[[1]]
p1_tmp <- plot[[1]]
ggsave(paste0(folder,"/","plot.network2.pdf"),p1_tmp,width = 7.5,height = 5)







######## 8. plot #######
dat_ITS <- tab.r_ITS [[2]]
node_ITS <- dat_ITS$net.cor.matrix$node
edge_ITS <- dat_ITS$net.cor.matrix$edge
head(node_ITS)
head(edge_ITS)



node_plot_ITS<- node_ITS %>%
  mutate(Phylum_plot = ifelse(Phylum %in% c("Ascomycota", "unassigned", "Basidiomycota", "Fungi_phy_Incertae_sedis",
                                            "Chytridiomycota", "Mortierellomycota", "Aphelidiomycota"), 
                              Phylum, "Others"))
  

edge_plot_ITS<-edge_ITS %>% 
  group_by(group) %>%
  mutate(links = n()) %>%
  ungroup() %>%
  mutate(label_plot = paste0(substr(label, 1, nchar(as.character(label)) - 2),
                             " ",
                             links,
                             ")"))


node_plot_ITS <- node_plot_ITS %>%
  mutate(label_plot = edge_plot_ITS$label_plot[match(group, sapply(strsplit(edge_plot_ITS$label_plot, ":"), `[`, 1))])



order_of_group<-c("CK", "Black", "Green", "White")
node_plot_ITS <- node_plot_ITS %>%
  mutate(label_plot = factor(label_plot, levels = unique(label_plot[order(match(sapply(strsplit(label_plot, ":"), `[`, 1), order_of_group))])))
edge_plot_ITS <- edge_plot_ITS %>%
  mutate(label_plot = factor(label_plot, levels = unique(label_plot[order(match(sapply(strsplit(label_plot, ":"), `[`, 1), order_of_group))])))





Palette_ITS<-c("Ascomycota"="#a6cee3","unassigned" ="#fb9a99",
                "Basidiomycota"="#e31a1c","Fungi_phy_Incertae_sedis"="#ff7f00",
                "Chytridiomycota" = "#6a3d9a", "Mortierellomycota" ="#1f78b4",
                "Aphelidiomycota" = "#33a02c", "Others"= "gray90")


p1 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = cor
),
data = edge_plot_ITS, linewidth= 0.3,alpha = 0.6) +
  geom_point(aes(X1, X2,
                 fill = Phylum_plot,
                 size = igraph.degree),
             pch = 21, data = node_plot_ITS,color = "gray40") +
  facet_wrap(.~ label_plot,scales="free_y",nrow = 1) +
  theme(strip.background = element_rect(fill = NA, color = NA))+
  theme(strip.text = element_text(size = 14, hjust = 0.5, face = "bold") ) + 

  scale_colour_manual(values = c("#6D98B5","#D48852"), name = "Correlation") +
  scale_fill_manual(values = Palette_ITS, name = "Phylum") +
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
ggsave(paste0(folder,"/","network_ITS.pdf"),p1,width = 15,height = 4)

## Notes: The degree of a vertex is its most basic structural property, the number of its adjacent edges
## from https://igraph.org/r/doc/degree.html










####### 9. plot #######

Samples_Palette <- c("CK"="#EEA236FF",
                     "Black"="#D43F3AFF",
                     "Green"="#5CB85CFF", 
                     "White"="#46B8DAFF")

library(patchwork)
Robustness_list_1<- Robustness.Targeted.removal(ps = ps_ITS,
                                  corg = cortab_ITS,
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
                                              "weighted" = "weighted")))+ 
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

Robustness_list_2<- Robustness.Random.removal(ps = ps_ITS,
                                                corg = cortab_ITS,
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
Robustness_list_3<- natural.con.microp(ps = ps_ITS,
                                              corg = cortab_ITS,
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
Robustness_list_4<- negative.correlation.ratio(ps = ps_ITS,
                                       corg = cortab_ITS,
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
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5, vjust = -0.5))+
  ylim(0,5)+
  xlab(NULL)+
  ylab("Negative correlation ratio (%)")
p5
ggsave(paste0(folder,"/","Negative_correlation_ratio.pdf"),p5,width = 3,height = 4)






###### 13 merge figures #####

library(ggpubr)
p_all <- ggarrange(
  p1,
  ggarrange(p2, p3, p4, p5, ncol = 4, widths = c(4.5, 4.5, 3, 3)),
  ncol = 1,
  heights = c(8, 7) 
)
p_all


ggsave(paste0(folder,"/","Network_all_ITS.pdf"),p_all,width = 15,height = 8)


