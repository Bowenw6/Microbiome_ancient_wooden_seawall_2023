#Robustness
##核心计算思路来自yongqi shao老师的文章From surviving to thriving, the assembly processes of microbial communities in stone biodeterioration- A case study of the West Lake UNESCO World Heritage area in China
###针对本项目进行了一些个性化改动
rm(list = ls())


nc <- function(adj_matrix) {
  
  adj_matrix <- as.matrix(adj_matrix)
  adj_matrix[abs(adj_matrix) != 0] <- 1
  

  lambda <- eigen(adj_matrix, only.values = TRUE)$values
  lambda <- sort(lambda, decreasing = TRUE)
  

  lambda_sum <- 0
  N = length(lambda)
  for (i in 1:N) lambda_sum = lambda_sum + exp(lambda[i])
  lambda_average <- log(lambda_sum/N, base = exp(1))
  lambda_average
}




library(igraph)
library(ggplot2)

adj_matrix_16S<-read.table("Whole Network matrix 0.77 RMT_result_16S.txt",header = TRUE,
                           sep = '\t',row.names = 1)
adj_matrix_16S<-as.matrix(adj_matrix_16S)
adj_matrix_ITS<-read.table("Whole Network matrix 0.73 RMT result_ITS.txt",header = TRUE,
                           sep = '\t',row.names = 1)
adj_matrix_ITS<-as.matrix(adj_matrix_ITS)


net_attack<-function(adj_matrix){
  
  g <- graph_from_adjacency_matrix(adj_matrix, mode = 'undirected', diag = FALSE)
  
  order_degree<-order(degree(g),decreasing = T)
  order_betwe<-order(betweenness(g),decreasing = T)
  order_random<-sample(1:dim(adj_matrix)[1],dim(adj_matrix)[1])
  
  
  degree_dist <- table(degree(g))
  degree_num <- as.numeric(names(degree_dist))
  degree_count <- as.numeric(degree_dist)
  names(degree_count) <- degree_num
  degs <- rep(degree_num, degree_count)
  
  adj_matrix_rand1 <- adj_matrix
  adj_matrix_rand2 <- adj_matrix
  adj_matrix_rand3 <- adj_matrix
  
  natural_connectivity_rand1 <- nc(adj_matrix_rand1)
  natural_connectivity_rand2 <- nc(adj_matrix_rand2)
  natural_connectivity_rand3 <- nc(adj_matrix_rand3)
  
  for (i in 1:(dim(adj_matrix)[1]-1)) {
    
   
    adj_matrix_rand1_remove <- adj_matrix_rand1[-order_degree[1:i],-order_degree[1:i]]
    adj_matrix_rand2_remove <- adj_matrix_rand2[-order_betwe[1:i],-order_betwe[1:i]]
    adj_matrix_rand3_remove <- adj_matrix_rand3[-order_random[1:i],-order_random[1:i]]
    
    natural_connectivity_rand1 <- c(natural_connectivity_rand1, nc(adj_matrix_rand1_remove))
    natural_connectivity_rand2 <- c(natural_connectivity_rand2, nc(adj_matrix_rand2_remove))
    natural_connectivity_rand3 <- c(natural_connectivity_rand3, nc(adj_matrix_rand3_remove))
    if(i%%10==0) {print(i)}
  }
  
  dat <- data.frame(remove_node = rep(1:dim(adj_matrix)[1],3),
                    natural_connectivity = c(natural_connectivity_rand1, natural_connectivity_rand2, natural_connectivity_rand3),
                    network = c( rep('Degree Attack', dim(adj_matrix)[1]), rep('Betweenness Attack', dim(adj_matrix)[1]), rep('Random Attack', dim(adj_matrix)[1])))
  

  dat$nc_per<-dat$natural_connectivity/max(dat$natural_connectivity)
  dat$remove_node_per<-dat$remove_node/max(dat$remove_node)
  return(dat)
  
}




Robustness_16S<-net_attack(adj_matrix_16S)
Robustness_ITS<-net_attack(adj_matrix_ITS)


Robustness_16S$group<-c("Prokaryotes")
Robustness_ITS$group<-c("Eukaryotes")

Robustness<-rbind(Robustness_16S,Robustness_ITS)


library(ggplot2)
#plot
ggplot(Robustness, aes(nc_per, remove_node_per, color = group)) +
  geom_smooth(se = FALSE) +
  facet_wrap(~network, ncol = 3)+
  theme_bw()+
  scale_color_manual(values = c("#B2182B","#E69F00"))+
  theme(axis.text.x = element_text(angle = 30))
  
  #save image 6*8 inches


