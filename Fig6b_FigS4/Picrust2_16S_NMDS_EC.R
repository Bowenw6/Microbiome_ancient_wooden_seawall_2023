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

picrust_16S_EC<-read.table("picrust2_out_pipeline_16S/EC_metagenome_out/pred_metagenome_unstrat.tsv",
                           sep = "\t", header = TRUE,row.names = 1)



library(vegan)
library(ggplot2)
library(ggrepel)

ec_filter <- picrust_16S_EC[rowSums(picrust_16S_EC > 0) >= 0.1 * ncol(picrust_16S_EC), ]

ec_rel <- decostand(ec_filter, method = "total", MARGIN = 2)  

ec_hel <- decostand(ec_rel, method = "hellinger", MARGIN = 2)

ec_hel_t <- t(ec_hel)

set.seed(123)
nmds <- metaMDS(ec_hel_t, distance = "bray", k = 2, trymax = 100)
stressplot(nmds)


nmds_points <- as.data.frame(scores(nmds, display = "sites"))
nmds_points$Group <- Group_infor[rownames(nmds_points), "Group"]  



set.seed(123)
adonis_result <- adonis2(ec_hel_t ~ Group, 
                         data = Group_infor, 
                         permutations = 999, 
                         method = "bray")


print(adonis_result)


set.seed(123)
bray_dist <- vegdist(ec_hel_t, method = "bray")  
anosim_result <- anosim(bray_dist, grouping = Group_infor$Group, permutations = 999)


print(anosim_result)



ordi <- ordiellipse(nmds, Group_infor$Group, kind = "se", conf = 0.95)


group_colors <- c("Black" = "#D43F3AFF", "CK" = "#EEA236FF", 
                  "White" = "#46B8DAFF", "Green" = "#5CB85CFF")




p <- ggplot(nmds_points, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Group), shape = 19, size = 4) +  
  scale_color_manual(values = group_colors,
                     breaks = names(group_colors)) +  
  scale_fill_manual(values = group_colors,
                    breaks = names(group_colors)) +   
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.95, 0.95),  
    legend.justification = c(1, 1),    
    legend.background = element_blank(),  
    legend.box.background = element_blank(),  
    legend.key = element_blank(),      
    legend.margin = margin(2, 2, 2, 2) 
  ) 


df_ell <- data.frame()
for(g in levels(factor(nmds_points$Group))){
  df_ell <- rbind(df_ell, 
                  cbind(as.data.frame(with(nmds_points[nmds_points$Group == g,],
                                           vegan:::veganCovEllipse(ordi[[g]]$cov,
                                                                   ordi[[g]]$center,
                                                                   ordi[[g]]$scale))),
                        Group = g))
}

p <- p + 
  geom_polygon(data = df_ell, 
               aes(x = NMDS1, y = NMDS2, 
                   fill = Group, color = Group),  
               alpha = 0.1,      
               linewidth = 0.2,   
               show.legend = FALSE)  

permanova_label <- paste("PERMANOVA: R²=", round(adonis_result$R2[1], 3),
                         " p<", round(adonis_result$`Pr(>F)`[1], 3),  
                         sep = "")
anosim_label <- paste("ANOSIM: R=", round(anosim_result$statistic, 3),
                      " p<", anosim_result$signif, sep = "")

p <- p + annotate("text", x = min(nmds_points$NMDS1), y = max(nmds_points$NMDS2),
             label = paste(permanova_label, "\n", anosim_label),
             hjust = 0, vjust = 1, size = 4)
p

ggsave(plot=p, filename ="EC_NMDS_plot.pdf", 
       width = 7, 
       height = 5,
       units = "in")



