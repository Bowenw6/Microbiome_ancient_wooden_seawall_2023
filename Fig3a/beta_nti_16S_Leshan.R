
rm(list=ls())
library(vegan)
library(Hmisc) 

Group_infor<-read.table("Group_Leshan.txt",header = TRUE,sep = "\t",row.names = 1)

ASV_all<-read.table("ASV_Abundance_16S_split.txt",sep = "\t", header = TRUE,row.names = 1)
ASV_all<-ASV_all[,1:c(length(colnames(ASV_all))-7)]


ASV_filtered<-ASV_all[which(rowSums(ASV_all)>27),]
ASV_filtered<-ASV_filtered[,order(colnames(ASV_filtered))]
ASV_filtered_CK<-ASV_filtered[,c(7:12)]
ASV_filtered_Bla<-ASV_filtered[,c(1:6)]
ASV_filtered_Gre<-ASV_filtered[,c(13:20)]
ASV_filtered_Whi<-ASV_filtered[,c(21:27)]
ASV_filtered_All<-ASV_filtered


A=ASV_filtered_CK 
C=A/rowSums(A)
ASV_filtered_CK_1<-t(C)

ASV_filtered_CK_1[ASV_filtered_CK_1>0]<-1

ASV_filtered_CK_1<-t(ASV_filtered_CK_1)
ASV_filtered_CK<-ASV_filtered_CK[which(rowSums(ASV_filtered_CK_1)>=3),]
rm(A,C,ASV_filtered_CK_1)


A=ASV_filtered_Bla 
C=A/rowSums(A)
ASV_filtered_Bla_1<-t(C)

ASV_filtered_Bla_1[ASV_filtered_Bla_1>0]<-1
ASV_filtered_Bla_1<-t(ASV_filtered_Bla_1)
ASV_filtered_Bla<-ASV_filtered_Bla[which(rowSums(ASV_filtered_Bla_1)>=3),]
rm(A,C,ASV_filtered_Bla_1)


A=ASV_filtered_Gre 
C=A/rowSums(A)
ASV_filtered_Gre_1<-t(C)
ASV_filtered_Gre_1[ASV_filtered_Gre_1>0]<-1
ASV_filtered_Gre_1<-t(ASV_filtered_Gre_1)
ASV_filtered_Gre<-ASV_filtered_Gre[which(rowSums(ASV_filtered_Gre_1)>=4),]
rm(A,C,ASV_filtered_Gre_1)


A=ASV_filtered_Whi 
C=A/rowSums(A)
ASV_filtered_Whi_1<-t(C)

ASV_filtered_Whi_1[ASV_filtered_Whi_1>0]<-1
ASV_filtered_Whi_1<-t(ASV_filtered_Whi_1)
ASV_filtered_Whi<-ASV_filtered_Whi[which(rowSums(ASV_filtered_Whi_1)>=3.5),]
rm(A,C,ASV_filtered_Whi_1)


A=ASV_filtered_All 
C=A/rowSums(A)
ASV_filtered_All_1<-t(C)

ASV_filtered_All_1[ASV_filtered_All_1>0]<-1
ASV_filtered_All_1<-t(ASV_filtered_All_1)
ASV_filtered_All<-ASV_filtered_All[which(rowSums(ASV_filtered_All_1)>6),]
rm(A,C,ASV_filtered_All_1)




summary(colSums(ASV_filtered_CK))
summary(colSums(ASV_filtered_Bla))
summary(colSums(ASV_filtered_Gre))
summary(colSums(ASV_filtered_Whi))
summary(colSums(ASV_filtered_All))







raw_Fasta_16S<-read.table("otus.filtered_16S.fa")
library(tidyverse)

format_fasta <- function(input_file, output_file) {
  lines <- readLines(input_file)
  
  output <- c()
  current_seq <- ""
  
  for (line in lines) {
    line <- trimws(line)
    if (nchar(line) == 0) next  
    
    if (grepl("^>", line)) {
      if (nchar(current_seq) > 0) {
        output <- c(output, current_seq)
        current_seq <- ""
      }
      output <- c(output, line)
    } else {
      current_seq <- paste0(current_seq, line)
    }
  }
  
  if (nchar(current_seq) > 0) {
    output <- c(output, current_seq)
  }
  
  writeLines(output, output_file)
}
format_fasta("otus.filtered_16S.fa", "fixed_otus.filtered_16S.fa")

Fasta_16S<-read.table("fixed_otus.filtered_16S.fa")

re_Fasta_16S<-Fasta_16S[seq(1,nrow(Fasta_16S),2),]
re_Fasta_16S<-as.data.frame(re_Fasta_16S)
re_Fasta_16S$V2<-Fasta_16S[seq(0,nrow(Fasta_16S),2),]
colnames(re_Fasta_16S)<- c("ASVs_names","Sequences")
re_Fasta_16S$ASVs_names<- gsub(">A*","A",re_Fasta_16S$ASVs_names) 
row.names(re_Fasta_16S)<- re_Fasta_16S$ASVs_names




ASV_filtered_CK_fasta_16S<-re_Fasta_16S[match(row.names(ASV_filtered_CK), row.names(re_Fasta_16S)),]
print(row.names(ASV_filtered_CK_fasta_16S)==row.names(ASV_filtered_CK))
ASV_filtered_CK_fasta_16S$ASVs_names<- gsub("ASV",">ASV",ASV_filtered_CK_fasta_16S$ASVs_names)
write.table(ASV_filtered_CK_fasta_16S, "ASV_filtered_CK_fasta_16S.fasta",sep =",", row.names =FALSE, col.names =FALSE, quote =FALSE)
filtered_CK_fasta_16S<-readLines("ASV_filtered_CK_fasta_16S.fasta",encoding = "UTF-8")
filtered_CK_fasta_16S<-gsub(",","\n",filtered_CK_fasta_16S)
writeLines(filtered_CK_fasta_16S,"filtered_CK_fasta_16S.fasta")
rm(ASV_filtered_CK_fasta_16S,filtered_CK_fasta_16S)




ASV_filtered_Bla_fasta_16S<-re_Fasta_16S[match(row.names(ASV_filtered_Bla), row.names(re_Fasta_16S)),]
print(row.names(ASV_filtered_Bla_fasta_16S)==row.names(ASV_filtered_Bla))
ASV_filtered_Bla_fasta_16S$ASVs_names<- gsub("ASV",">ASV",ASV_filtered_Bla_fasta_16S$ASVs_names)
write.table(ASV_filtered_Bla_fasta_16S, "ASV_filtered_Bla_fasta_16S.fasta",sep =",", row.names =FALSE, col.names =FALSE, quote =FALSE)
filtered_Bla_fasta_16S<-readLines("ASV_filtered_Bla_fasta_16S.fasta",encoding = "UTF-8")
filtered_Bla_fasta_16S<-gsub(",","\n",filtered_Bla_fasta_16S)
writeLines(filtered_Bla_fasta_16S,"filtered_Bla_fasta_16S.fasta")
rm(ASV_filtered_Bla_fasta_16S,filtered_Bla_fasta_16S)







ASV_filtered_Gre_fasta_16S<-re_Fasta_16S[match(row.names(ASV_filtered_Gre), row.names(re_Fasta_16S)),]

print(row.names(ASV_filtered_Gre_fasta_16S)==row.names(ASV_filtered_Gre))

ASV_filtered_Gre_fasta_16S$ASVs_names<- gsub("ASV",">ASV",ASV_filtered_Gre_fasta_16S$ASVs_names)

write.table(ASV_filtered_Gre_fasta_16S, "ASV_filtered_Gre_fasta_16S.fasta",sep =",", row.names =FALSE, col.names =FALSE, quote =FALSE)

filtered_Gre_fasta_16S<-readLines("ASV_filtered_Gre_fasta_16S.fasta",encoding = "UTF-8")

filtered_Gre_fasta_16S<-gsub(",","\n",filtered_Gre_fasta_16S)
writeLines(filtered_Gre_fasta_16S,"filtered_Gre_fasta_16S.fasta")
rm(ASV_filtered_Gre_fasta_16S,filtered_Gre_fasta_16S)





ASV_filtered_Whi_fasta_16S<-re_Fasta_16S[match(row.names(ASV_filtered_Whi), row.names(re_Fasta_16S)),]
print(row.names(ASV_filtered_Whi_fasta_16S)==row.names(ASV_filtered_Whi))
ASV_filtered_Whi_fasta_16S$ASVs_names<- gsub("ASV",">ASV",ASV_filtered_Whi_fasta_16S$ASVs_names)
write.table(ASV_filtered_Whi_fasta_16S, "ASV_filtered_Whi_fasta_16S.fasta",sep =",", row.names =FALSE, col.names =FALSE, quote =FALSE)
filtered_Whi_fasta_16S<-readLines("ASV_filtered_Whi_fasta_16S.fasta",encoding = "UTF-8")
filtered_Whi_fasta_16S<-gsub(",","\n",filtered_Whi_fasta_16S)
writeLines(filtered_Whi_fasta_16S,"filtered_Whi_fasta_16S.fasta")
rm(ASV_filtered_Whi_fasta_16S,filtered_Whi_fasta_16S)


ASV_filtered_All_fasta_16S<-re_Fasta_16S[match(row.names(ASV_filtered_All), row.names(re_Fasta_16S)),]
print(row.names(ASV_filtered_All_fasta_16S)==row.names(ASV_filtered_All))
ASV_filtered_All_fasta_16S$ASVs_names<- gsub("ASV",">ASV",ASV_filtered_All_fasta_16S$ASVs_names)
write.table(ASV_filtered_All_fasta_16S, "ASV_filtered_All_fasta_16S.fasta",sep =",", row.names =FALSE, col.names =FALSE, quote =FALSE)
filtered_All_fasta_16S<-readLines("ASV_filtered_All_fasta_16S.fasta",encoding = "UTF-8")
filtered_All_fasta_16S<-gsub(",","\n",filtered_All_fasta_16S)
writeLines(filtered_All_fasta_16S,"filtered_All_fasta_16S.fasta")
rm(ASV_filtered_All_fasta_16S,filtered_All_fasta_16S)



############ Run this script in Mac: alignment_tree.sh ###############
############ or run the following commands respectively ##############
### in MacOS command line, perform alignment with muscle
#nohup muscle -align filtered_CK_fasta_16S.fasta -output aligned_filtered_CK_fasta_16S.fasta &
###  in MacOS command line, construct tree with iqtree
#nohup iqtree -s aligned_filtered_CK_fasta_16S.fasta -m MFP -B 1000 --bnni -cptime 6000 -T AUTO &

### in MacOS command line, perform alignment with muscle
#nohup muscle -align filtered_Bla_fasta_16S.fasta -output aligned_filtered_Bla_fasta_16S.fasta &
###  in MacOS command line, construct tree with iqtree
#nohup iqtree -s aligned_filtered_Bla_fasta_16S.fasta -m MFP -B 1000 --bnni -cptime 6000 -T AUTO &

### in MacOS command line, perform alignment with muscle
#nohup muscle -align filtered_Gre_fasta_16S.fasta -output aligned_filtered_Gre_fasta_16S.fasta &
###  in MacOS command line, construct tree with iqtree
#nohup iqtree -s aligned_filtered_Gre_fasta_16S.fasta -m MFP -B 1000 --bnni -cptime 6000 -T AUTO &

### in MacOS command line, perform alignment with muscle
#nohup muscle -align filtered_Whi_fasta_16S.fasta -output aligned_filtered_Whi_fasta_16S.fasta &
###  in MacOS command line, construct tree with iqtree
#nohup iqtree -s aligned_filtered_Whi_fasta_16S.fasta -m MFP -B 1000 --bnni -cptime 6000 -T AUTO &

### in MacOS command line, perform alignment with muscle
#nohup muscle -align filtered_All_fasta_16S.fasta -output aligned_filtered_All_fasta_16S.fasta &
###  in MacOS command line, construct tree with iqtree
#nohup iqtree -s aligned_filtered_All_fasta_16S.fasta -m MFP -B 1000 --bnni -cptime 6000 -T AUTO &

### To remove window returns if needed
# dos2unix alignment_tree.sh









library(treeio)
tree_16S_CK<-read.tree("aligned_filtered_CK_fasta_16S.fasta.treefile")
tree_16S_Bla<-read.tree("aligned_filtered_Bla_fasta_16S.fasta.treefile")
tree_16S_Gre<-read.tree("aligned_filtered_Gre_fasta_16S.fasta.treefile")
tree_16S_Whi<-read.tree("aligned_filtered_Whi_fasta_16S.fasta.treefile")
tree_16S_All<-read.tree("aligned_filtered_All_fasta_16S.fasta.treefile")




ASV_filtered_CK<-t(ASV_filtered_CK)
head(rowSums(ASV_filtered_CK))
ASV_filtered_CK<-as.data.frame(ASV_filtered_CK)

ASV_filtered_Bla<-t(ASV_filtered_Bla)
head(rowSums(ASV_filtered_Bla))
ASV_filtered_Bla<-as.data.frame(ASV_filtered_Bla)

ASV_filtered_Gre<-t(ASV_filtered_Gre)
head(rowSums(ASV_filtered_Gre))
ASV_filtered_Gre<-as.data.frame(ASV_filtered_Gre)

ASV_filtered_Whi<-t(ASV_filtered_Whi)
head(rowSums(ASV_filtered_Whi))
ASV_filtered_Whi<-as.data.frame(ASV_filtered_Whi)

ASV_filtered_All<-t(ASV_filtered_All)
head(rowSums(ASV_filtered_All))
ASV_filtered_All<-as.data.frame(ASV_filtered_All)








library(picante)
beta_nti <- function(otu_niche,otu_tree,reps,threads){
  library(picante)
  library(ape)
  prune_tree<-prune.sample(otu_niche,otu_tree) 
 
  otu_phydist <- cophenetic(prune_tree)
  match.phylo.otu = match.phylo.comm(prune_tree, otu_niche)
  
  beta.mntd.weighted = as.matrix(comdistnt(match.phylo.otu$comm,cophenetic(match.phylo.otu$phy),abundance.weighted=T));
  
  
  beta.reps = reps; 
  rand.weighted.bMNTD.comp = NULL
  dim(rand.weighted.bMNTD.comp)
  library(abind)
  arraybind <- function(...){
    abind(...,along = 3,force.array=TRUE)
  }
  library(foreach)
  library(doParallel)
  registerDoParallel(cores = threads)
  rand.weighted.bMNTD.comp <- foreach (rep = 1:beta.reps, .combine = "arraybind") %dopar%{
    as.matrix(comdistnt(match.phylo.otu$comm,taxaShuffle(cophenetic(match.phylo.otu$phy)),abundance.weighted=T,exclude.conspecifics = F))
  }
  weighted.bNTI = matrix(c(NA),nrow=nrow(match.phylo.otu$comm),ncol=nrow(match.phylo.otu$comm));
  dim(weighted.bNTI);
  for (columns in 1:(nrow(match.phylo.otu$comm)-1)) {
    for (rows in (columns+1):nrow(match.phylo.otu$comm)) {
      rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
      weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
      rm("rand.vals")
    }
  }
  rownames(weighted.bNTI) = rownames(match.phylo.otu$comm)
  colnames(weighted.bNTI) = rownames(match.phylo.otu$comm)
  return(weighted.bNTI)
}



raup_crick_abundance = function(spXsite,
                                plot_names_in_col1 = TRUE,
                                classic_metric = FALSE,
                                split_ties = TRUE,
                                reps = 999,
                                set_all_species_equal = FALSE,
                                as.distance.matrix = TRUE,
                                report_similarity = FALSE,
                                threads = 1) {
  if (plot_names_in_col1) {
    row.names(spXsite) <- spXsite[, 1]
    spXsite <- spXsite[, -1]
  }
  library(doParallel)
  library(foreach)
  n_sites <- nrow(spXsite)
  gamma <- ncol(spXsite)
  ##build a site by site matrix for the results, with the names of the sites in the row and col names:
  results <-
    matrix(
      data = NA,
      nrow = n_sites,
      ncol = n_sites,
      dimnames = list(row.names(spXsite), row.names(spXsite))
    )
  ##make the spXsite matrix into a new, pres/abs. matrix:
  ceiling(spXsite / max(spXsite)) -> spXsite.inc
  ##create an occurrence vector- used to give more weight to widely distributed species in the null model:
  occur <- apply(spXsite.inc, MARGIN = 2, FUN = sum)
  ##create an abundance vector- used to give more weight to abundant species in the second step of the null model:
  abundance <- apply(spXsite, MARGIN = 2, FUN = sum)
  ##make_null:
  ##looping over each pairwise community combination:
  for (null.one in 1:(nrow(spXsite) - 1)) {
    for (null.two in (null.one + 1):nrow(spXsite)) {
      null_bray_curtis <- NULL
      registerDoParallel(cores = threads)
      null_bray_curtis <- foreach (i = 1:reps,.combine = "c") %dopar%{
        ##two empty null communities of size gamma:
        com1 <- rep(0, gamma)
        com2 <- rep(0, gamma)
        ##add observed number of species to com1, weighting by species occurrence frequencies:
        com1[sample(1:gamma,
                    sum(spXsite.inc[null.one, ]),
                    replace = FALSE,
                    prob = occur)] <- 1
        com1.samp.sp = sample(which(com1 > 0),
                              (sum(spXsite[null.one, ]) - sum(com1)),
                              replace = TRUE,
                              prob = abundance[which(com1 > 0)])
        com1.samp.sp = cbind(com1.samp.sp, 1)
        # head(com1.samp.sp);
        com1.sp.counts = as.data.frame(tapply(com1.samp.sp[, 2], com1.samp.sp[, 1], FUN = sum))
        colnames(com1.sp.counts) = 'counts'
        # head(com1.sp.counts);
        com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts))
        # head(com1.sp.counts);
        com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts
        # com1;
        #sum(com1) - sum(spXsite[null.one,]); ## this should be zero if everything work properly
        rm('com1.samp.sp', 'com1.sp.counts')
        ##same for com2:
        com2[sample(1:gamma,
                    sum(spXsite.inc[null.two, ]),
                    replace = FALSE,
                    prob = occur)] <- 1
        com2.samp.sp = sample(which(com2 > 0),
                              (sum(spXsite[null.two, ]) - sum(com2)),
                              replace = TRUE,
                              prob = abundance[which(com2 > 0)])
        com2.samp.sp = cbind(com2.samp.sp, 1)
        # head(com2.samp.sp);
        com2.sp.counts = as.data.frame(tapply(com2.samp.sp[, 2], com2.samp.sp[, 1], FUN =
                                                sum))
        colnames(com2.sp.counts) = 'counts'
        # head(com2.sp.counts);
        com2.sp.counts$sp = as.numeric(rownames(com2.sp.counts))
        # head(com2.sp.counts);
        com2[com2.sp.counts$sp] = com2[com2.sp.counts$sp] + com2.sp.counts$counts
        # com2;
        # sum(com2) - sum(spXsite[null.two,]); ## this should be zero if everything work properly
        rm('com2.samp.sp', 'com2.sp.counts')
        null.spXsite = rbind(com1, com2)
        # null.spXsite;
        ##calculate null bray curtis
        null_bray_curtis = vegdist(null.spXsite, method = 'bray')
        #print(c(i,date()))
      }
      # end reps loop
      ## empirically observed bray curtis
      obs.bray = vegdist(spXsite[c(null.one, null.two), ], method = 'bray')
      ##how many null observations is the observed value tied with?
      num_exact_matching_in_null = sum(null_bray_curtis == obs.bray)
      ##how many null values are smaller than the observed *dissimilarity*?
      num_less_than_in_null = sum(null_bray_curtis < obs.bray)
      rc = (num_less_than_in_null) / reps
      # rc;
      if (split_ties) {
        rc = ((
          num_less_than_in_null + (num_exact_matching_in_null) / 2
        ) / reps)
      }
      if (!classic_metric) {
        ##our modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
        rc = (rc - .5) * 2
      }
      results[null.two, null.one] = round(rc, digits = 2)
      ##store the metric in the results matrix
      print(c(null.one, null.two, date()))
    }
    ## end null.two loop
  }
  ## end null.one loop
  if (as.distance.matrix) {
    ## return as distance matrix if so desired
    results <- as.dist(results)
  }
  return(results)
}
## end function






paired_beta_nti_CK <- beta_nti(otu_niche = ASV_filtered_CK, otu_tree = tree_16S_CK,1000,8)
write.table(paired_beta_nti_CK,"paired_beta_nti_CK.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)

paired_beta_nti_Bla <- beta_nti(otu_niche = ASV_filtered_Bla, otu_tree = tree_16S_Bla,1000,8)
write.table(paired_beta_nti_Bla,"paired_beta_nti_Bla.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)

paired_beta_nti_Gre <- beta_nti(otu_niche = ASV_filtered_Gre, otu_tree = tree_16S_Gre,1000,8)
write.table(paired_beta_nti_Gre,"paired_beta_nti_Gre.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)

paired_beta_nti_Whi <- beta_nti(otu_niche = ASV_filtered_Whi, otu_tree = tree_16S_Whi,1000,8)
write.table(paired_beta_nti_Whi,"paired_beta_nti_Whi.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)

paired_beta_nti_All <- beta_nti(otu_niche = ASV_filtered_All, otu_tree = tree_16S_All,1000,8)
write.table(paired_beta_nti_All,"paired_beta_nti_All.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)










paired_raup_crick_CK<-raup_crick_abundance(ASV_filtered_CK,plot_names_in_col1 = FALSE,threads = 8)
paired_raup_crick_CK<-as.matrix(paired_raup_crick_CK)
paired_raup_crick_CK<-as.data.frame(paired_raup_crick_CK)
write.table(paired_raup_crick_CK, "paired_raup_crick_CK.txt", 
            sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE)

paired_raup_crick_Bla<-raup_crick_abundance(ASV_filtered_Bla,plot_names_in_col1 = FALSE,threads = 8)
paired_raup_crick_Bla<-as.matrix(paired_raup_crick_Bla)
paired_raup_crick_Bla<-as.data.frame(paired_raup_crick_Bla)
write.table(paired_raup_crick_Bla, "paired_raup_crick_Bla.txt", 
            sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE)

paired_raup_crick_Gre<-raup_crick_abundance(ASV_filtered_Gre,plot_names_in_col1 = FALSE,threads = 8)
paired_raup_crick_Gre<-as.matrix(paired_raup_crick_Gre)
paired_raup_crick_Gre<-as.data.frame(paired_raup_crick_Gre)
write.table(paired_raup_crick_Gre, "paired_raup_crick_Gre.txt", 
            sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE)

paired_raup_crick_Whi<-raup_crick_abundance(ASV_filtered_Whi,plot_names_in_col1 = FALSE,threads = 8)
paired_raup_crick_Whi<-as.matrix(paired_raup_crick_Whi)
paired_raup_crick_Whi<-as.data.frame(paired_raup_crick_Whi)
write.table(paired_raup_crick_Whi, "paired_raup_crick_Whi.txt", 
            sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE)

paired_raup_crick_All<-raup_crick_abundance(ASV_filtered_All,plot_names_in_col1 = FALSE,threads = 8)
paired_raup_crick_All<-as.matrix(paired_raup_crick_All)
paired_raup_crick_All<-as.data.frame(paired_raup_crick_All)
write.table(paired_raup_crick_All, "paired_raup_crick_All.txt", 
            sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE)











library(reshape2)
library(dplyr)
paired_beta_nti_CK<-read.table(file = "paired_beta_nti_CK.txt", sep = '\t', header = TRUE, row.names = 1)
paired_beta_nti_Bla<-read.table(file = "paired_beta_nti_Bla.txt", sep = '\t', header = TRUE, row.names = 1)
paired_beta_nti_Gre<-read.table(file = "paired_beta_nti_Gre.txt", sep = '\t', header = TRUE, row.names = 1)
paired_beta_nti_Whi<-read.table(file = "paired_beta_nti_Whi.txt", sep = '\t', header = TRUE, row.names = 1)
paired_beta_nti_All<-read.table(file = "paired_beta_nti_All.txt", sep = '\t', header = TRUE, row.names = 1)


paired_raup_crick_CK<- read.table(file = "paired_raup_crick_CK.txt", sep = '\t', header = TRUE, row.names = 1)
paired_raup_crick_Bla<- read.table(file = "paired_raup_crick_Bla.txt", sep = '\t', header = TRUE, row.names = 1)
paired_raup_crick_Gre<- read.table(file = "paired_raup_crick_Gre.txt", sep = '\t', header = TRUE, row.names = 1)
paired_raup_crick_Whi<- read.table(file = "paired_raup_crick_Whi.txt", sep = '\t', header = TRUE, row.names = 1)
paired_raup_crick_All<- read.table(file = "paired_raup_crick_All.txt", sep = '\t', header = TRUE, row.names = 1)


dist_melt <- function(dist_a) {
  a <- dist_a 
  a <- as.matrix(a)
  a[upper.tri(a)] <- 10000
  diag(a) <- 10000
  betamat <- melt(a)
  betamat <- betamat[!(betamat$value %in% c(10000)), ]
  return(betamat)
}



beta_nti_CK<- dist_melt(paired_beta_nti_CK)
write.table(beta_nti_CK,"beta_nti_CK_melted.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)

beta_nti_Bla<- dist_melt(paired_beta_nti_Bla)
write.table(beta_nti_Bla,"beta_nti_Bla_melted.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)

beta_nti_Gre<- dist_melt(paired_beta_nti_Gre)
write.table(beta_nti_Gre,"beta_nti_Gre_melted.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)

beta_nti_Whi<- dist_melt(paired_beta_nti_Whi)
write.table(beta_nti_Whi,"beta_nti_Whi_melted.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)

beta_nti_All<- dist_melt(paired_beta_nti_All)
write.table(beta_nti_All,"beta_nti_All_melted.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)




raup_crick_CK<- dist_melt(paired_raup_crick_CK)
write.table(raup_crick_CK,"raup_crick_CK.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)

raup_crick_Bla<- dist_melt(paired_raup_crick_Bla)
write.table(raup_crick_Bla,"raup_crick_Bla.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)

raup_crick_Gre<- dist_melt(paired_raup_crick_Gre)
write.table(raup_crick_Gre,"raup_crick_Gre.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)

raup_crick_Whi<- dist_melt(paired_raup_crick_Whi)
write.table(raup_crick_Whi,"raup_crick_Whi.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)

raup_crick_All<- dist_melt(paired_raup_crick_All)
write.table(raup_crick_All,"raup_crick_All.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)









CK_null_mat <- cbind(beta_nti_CK, raup_crick_CK[, c(3, 2)])[, 1:4]

CK_null_mat$assembly[CK_null_mat$value <= -2] <- "Homo_Selection"
CK_null_mat$assembly[CK_null_mat$value >= 2] <- "Hetero_Selection"
CK_null_mat$assembly[(abs(CK_null_mat$value) < 2) &
                         CK_null_mat$value.1 >= 0.95] <-
  "Dispersal_Limitation"
CK_null_mat$assembly[(abs(CK_null_mat$value) < 2) &
                         CK_null_mat$value.1 <= -0.95] <-
  "Homogenizing_Dispersal"
CK_null_mat$assembly[(abs(CK_null_mat$value) < 2) &
                         CK_null_mat$value.1 > -0.95 &
                         CK_null_mat$value.1 < 0.95] <- "Drift"
colnames(CK_null_mat)=c("Var1","Var2","beta_NTI","RCbray","assembly")
write.table(CK_null_mat,"CK_null_mat.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)


Summarise_CK_null_mat<- CK_null_mat %>% group_by(assembly) %>% summarise(n = n())












Bla_null_mat <- cbind(beta_nti_Bla, raup_crick_Bla[, c(3, 2)])[, 1:4]

Bla_null_mat$assembly[Bla_null_mat$value <= -2] <- "Homo_Selection"
Bla_null_mat$assembly[Bla_null_mat$value >= 2] <- "Hetero_Selection"
Bla_null_mat$assembly[(abs(Bla_null_mat$value) < 2) &
                        Bla_null_mat$value.1 >= 0.95] <-
  "Dispersal_Limitation"
Bla_null_mat$assembly[(abs(Bla_null_mat$value) < 2) &
                       Bla_null_mat$value.1 <= -0.95] <-
  "Homogenizing_Dispersal"
Bla_null_mat$assembly[(abs(Bla_null_mat$value) < 2) &
                        Bla_null_mat$value.1 > -0.95 &
                        Bla_null_mat$value.1 < 0.95] <- "Drift"
colnames(Bla_null_mat)=c("Var1","Var2","beta_NTI","RCbray","assembly")
write.table(Bla_null_mat,"Bla_null_mat.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)


Summarise_Bla_null_mat<- Bla_null_mat %>% group_by(assembly) %>% summarise(n = n())










Gre_null_mat <- cbind(beta_nti_Gre, raup_crick_Gre[, c(3, 2)])[, 1:4]

Gre_null_mat$assembly[Gre_null_mat$value <= -2] <- "Homo_Selection"
Gre_null_mat$assembly[Gre_null_mat$value >= 2] <- "Hetero_Selection"
Gre_null_mat$assembly[(abs(Gre_null_mat$value) < 2) &
                        Gre_null_mat$value.1 >= 0.95] <-
  "Dispersal_Limitation"
Gre_null_mat$assembly[(abs(Gre_null_mat$value) < 2) &
                        Gre_null_mat$value.1 <= -0.95] <-
  "Homogenizing_Dispersal"
Gre_null_mat$assembly[(abs(Gre_null_mat$value) < 2) &
                        Gre_null_mat$value.1 > -0.95 &
                        Gre_null_mat$value.1 < 0.95] <- "Drift"
colnames(Gre_null_mat)=c("Var1","Var2","beta_NTI","RCbray","assembly")
write.table(Gre_null_mat,"Gre_null_mat.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)


Summarise_Gre_null_mat<- Gre_null_mat %>% group_by(assembly) %>% summarise(n = n())










Whi_null_mat <- cbind(beta_nti_Whi, raup_crick_Whi[, c(3, 2)])[, 1:4]

Whi_null_mat$assembly[Whi_null_mat$value <= -2] <- "Homo_Selection"
Whi_null_mat$assembly[Whi_null_mat$value >= 2] <- "Hetero_Selection"
Whi_null_mat$assembly[(abs(Whi_null_mat$value) < 2) &
                        Whi_null_mat$value.1 >= 0.95] <-
  "Dispersal_Limitation"
Whi_null_mat$assembly[(abs(Whi_null_mat$value) < 2) &
                        Whi_null_mat$value.1 <= -0.95] <-
  "Homogenizing_Dispersal"
Whi_null_mat$assembly[(abs(Whi_null_mat$value) < 2) &
                        Whi_null_mat$value.1 > -0.95 &
                        Whi_null_mat$value.1 < 0.95] <- "Drift"
colnames(Whi_null_mat)=c("Var1","Var2","beta_NTI","RCbray","assembly")
write.table(Whi_null_mat,"Whi_null_mat.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)


Summarise_Whi_null_mat<- Whi_null_mat %>% group_by(assembly) %>% summarise(n = n())







All_null_mat <- cbind(beta_nti_All, raup_crick_All[, c(3, 2)])[, 1:4]

All_null_mat$assembly[All_null_mat$value <= -2] <- "Homo_Selection"
All_null_mat$assembly[All_null_mat$value >= 2] <- "Hetero_Selection"
All_null_mat$assembly[(abs(All_null_mat$value) < 2) &
                        All_null_mat$value.1 >= 0.95] <-
  "Dispersal_Limitation"
All_null_mat$assembly[(abs(All_null_mat$value) < 2) &
                        All_null_mat$value.1 <= -0.95] <-
  "Homogenizing_Dispersal"
All_null_mat$assembly[(abs(All_null_mat$value) < 2) &
                        All_null_mat$value.1 > -0.95 &
                        All_null_mat$value.1 < 0.95] <- "Drift"
colnames(All_null_mat)=c("Var1","Var2","beta_NTI","RCbray","assembly")
write.table(All_null_mat,"All_null_mat.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)


Summarise_All_null_mat<- All_null_mat %>% group_by(assembly) %>% summarise(n = n())

















library(ggplot2)
library(RColorBrewer)
display.brewer.all()
display.brewer.pal(11,'RdYlBu')
brewer.pal(11,'RdYlBu')

# Create data
CK_Donut_data <- as.data.frame(Summarise_CK_null_mat)
# Compute percentages
CK_Donut_data$fraction <- CK_Donut_data$n / sum(CK_Donut_data$n)
# Compute the cumulative percentages (top of each rectangle)
CK_Donut_data$ymax <- cumsum(CK_Donut_data$fraction)
# Compute the bottom of each rectangle
CK_Donut_data$ymin <- c(0, head(CK_Donut_data$ymax, n=-1))
# Compute label position
CK_Donut_data$labelPosition <- (CK_Donut_data$ymax + CK_Donut_data$ymin) / 2
# Compute a good label
CK_Donut_data$label <- paste0(CK_Donut_data$assembly, "\n ", round(CK_Donut_data$n / sum(CK_Donut_data$n) *100, 2), "%")


# Make the plot
p1<- ggplot(CK_Donut_data, aes(ymax=ymax, ymin=ymin, xmax=3, xmin=1, fill=assembly)) +
  geom_rect() +
  geom_text( x=5, aes(y=labelPosition, label=label, color=assembly), size=6) + # x here controls label position (inner / outer)
  geom_text( x=-1, y=0, label=paste("CK"), size=8) +
  scale_fill_manual(values=(c(Hetero_Selection="#D73027", Homo_Selection="#313695", 
                    Dispersal_Limitation="#FDAE61", Drift="gray", Homogenizing_Dispersal="#74ADD1"))) +
  scale_color_manual(values=(c(Hetero_Selection="#D73027", Homo_Selection="#313695", 
                               Dispersal_Limitation="#FDAE61", Drift="gray", Homogenizing_Dispersal="#74ADD1"))) +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")
p1
#save image 5 inches * 5 inches
#then rearrange the positions of texts with AI









library(ggplot2)
library(RColorBrewer)
display.brewer.all()
display.brewer.pal(11,'RdYlBu')
brewer.pal(11,'RdYlBu')

# Create data
Black_Donut_data <- as.data.frame(Summarise_Bla_null_mat)
# Compute percentages
Black_Donut_data$fraction <- Black_Donut_data$n / sum(Black_Donut_data$n)
# Compute the cumulative percentages (top of each rectangle)
Black_Donut_data$ymax <- cumsum(Black_Donut_data$fraction)
# Compute the bottom of each rectangle
Black_Donut_data$ymin <- c(0, head(Black_Donut_data$ymax, n=-1))
# Compute label position
Black_Donut_data$labelPosition <- (Black_Donut_data$ymax + Black_Donut_data$ymin) / 2
# Compute a good label
Black_Donut_data$label <- paste0(Black_Donut_data$assembly, "\n ", round(Black_Donut_data$n / sum(Black_Donut_data$n) *100, 2), "%")


# Make the plot
p2<- ggplot(Black_Donut_data, aes(ymax=ymax, ymin=ymin, xmax=3, xmin=1, fill=assembly)) +
  geom_rect() +
  geom_text( x=5, aes(y=labelPosition, label=label, color=assembly), size=6) + # x here controls label position (inner / outer)
  geom_text( x=-1, y=0, label=paste("Black"), size=8) +
  scale_fill_manual(values=(c(Hetero_Selection="#D73027", Homo_Selection="#313695", 
                              Dispersal_Limitation="#FDAE61", Drift="gray", Homogenizing_Dispersal="#74ADD1"))) +
  scale_color_manual(values=(c(Hetero_Selection="#D73027", Homo_Selection="#313695", 
                               Dispersal_Limitation="#FDAE61", Drift="gray", Homogenizing_Dispersal="#74ADD1"))) +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")
p2
#save image 5 inches * 5 inches
#then rearrange the positions of texts with AI











library(ggplot2)
library(RColorBrewer)
display.brewer.all()
display.brewer.pal(11,'RdYlBu')
brewer.pal(11,'RdYlBu')

# Create data
Green_Donut_data <- as.data.frame(Summarise_Gre_null_mat)
# Compute percentages
Green_Donut_data$fraction <- Green_Donut_data$n / sum(Green_Donut_data$n)
# Compute the cumulative percentages (top of each rectangle)
Green_Donut_data$ymax <- cumsum(Green_Donut_data$fraction)
# Compute the bottom of each rectangle
Green_Donut_data$ymin <- c(0, head(Green_Donut_data$ymax, n=-1))
# Compute label position
Green_Donut_data$labelPosition <- (Green_Donut_data$ymax + Green_Donut_data$ymin) / 2
# Compute a good label
Green_Donut_data$label <- paste0(Green_Donut_data$assembly, "\n ", round(Green_Donut_data$n / sum(Green_Donut_data$n) *100, 2), "%")


# Make the plot
p3<- ggplot(Green_Donut_data, aes(ymax=ymax, ymin=ymin, xmax=3, xmin=1, fill=assembly)) +
  geom_rect() +
  geom_text( x=5, aes(y=labelPosition, label=label, color=assembly), size=6) + # x here controls label position (inner / outer)
  geom_text( x=-1, y=0, label=paste("Green"), size=8) +
  scale_fill_manual(values=(c(Hetero_Selection="#D73027", Homo_Selection="#313695", 
                              Dispersal_Limitation="#FDAE61", Drift="gray", Homogenizing_Dispersal="#74ADD1"))) +
  scale_color_manual(values=(c(Hetero_Selection="#D73027", Homo_Selection="#313695", 
                               Dispersal_Limitation="#FDAE61", Drift="gray", Homogenizing_Dispersal="#74ADD1"))) +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")
p3
#save image 5 inches * 5 inches
#then rearrange the positions of texts with AI











library(ggplot2)
library(RColorBrewer)
display.brewer.all()
display.brewer.pal(11,'RdYlBu')
brewer.pal(11,'RdYlBu')

# Create data
White_Donut_data <- as.data.frame(Summarise_Whi_null_mat)
# Compute percentages
White_Donut_data$fraction <- White_Donut_data$n / sum(White_Donut_data$n)
# Compute the cumulative percentages (top of each rectangle)
White_Donut_data$ymax <- cumsum(White_Donut_data$fraction)
# Compute the bottom of each rectangle
White_Donut_data$ymin <- c(0, head(White_Donut_data$ymax, n=-1))
# Compute label position
White_Donut_data$labelPosition <- (White_Donut_data$ymax + White_Donut_data$ymin) / 2
# Compute a good label
White_Donut_data$label <- paste0(White_Donut_data$assembly, "\n ", round(White_Donut_data$n / sum(White_Donut_data$n) *100, 2), "%")


# Make the plot
p4<- ggplot(White_Donut_data, aes(ymax=ymax, ymin=ymin, xmax=3, xmin=1, fill=assembly)) +
  geom_rect() +
  geom_text( x=5, aes(y=labelPosition, label=label, color=assembly), size=6) + # x here controls label position (inner / outer)
  geom_text( x=-1, y=0, label=paste("White"), size=8) +
  scale_fill_manual(values=(c(Hetero_Selection="#D73027", Homo_Selection="#313695", 
                              Dispersal_Limitation="#FDAE61", Drift="gray", Homogenizing_Dispersal="#74ADD1"))) +
  scale_color_manual(values=(c(Hetero_Selection="#D73027", Homo_Selection="#313695", 
                               Dispersal_Limitation="#FDAE61", Drift="gray", Homogenizing_Dispersal="#74ADD1"))) +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")
p4
#save image 5 inches * 5 inches
#then rearrange the positions of texts with AI









library(ggplot2)
library(RColorBrewer)
display.brewer.all()
display.brewer.pal(11,'RdYlBu')
brewer.pal(11,'RdYlBu')

# Create data
All_Donut_data <- as.data.frame(Summarise_All_null_mat)
# Compute percentages
All_Donut_data$fraction <- All_Donut_data$n / sum(All_Donut_data$n)
# Compute the cumulative percentages (top of each rectangle)
All_Donut_data$ymax <- cumsum(All_Donut_data$fraction)
# Compute the bottom of each rectangle
All_Donut_data$ymin <- c(0, head(All_Donut_data$ymax, n=-1))
# Compute label position
All_Donut_data$labelPosition <- (All_Donut_data$ymax + All_Donut_data$ymin) / 2
# Compute a good label
All_Donut_data$label <- paste0(All_Donut_data$assembly, "\n ", round(All_Donut_data$n / sum(All_Donut_data$n) *100, 2), "%")


# Make the plot
p5<- ggplot(All_Donut_data, aes(ymax=ymax, ymin=ymin, xmax=3, xmin=1, fill=assembly)) +
  geom_rect() +
  geom_text( x=5, aes(y=labelPosition, label=label, color=assembly), size=6) + # x here controls label position (inner / outer)
  geom_text( x=-1, y=0, label=paste("All"), size=8) +
  scale_fill_manual(values=(c(Hetero_Selection="#D73027", Homo_Selection="#313695", 
                              Dispersal_Limitation="#FDAE61", Drift="gray", Homogenizing_Dispersal="#74ADD1"))) +
  scale_color_manual(values=(c(Hetero_Selection="#D73027", Homo_Selection="#313695", 
                               Dispersal_Limitation="#FDAE61", Drift="gray", Homogenizing_Dispersal="#74ADD1"))) +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")
p5
#save image 5 inches * 5 inches
#then rearrange the positions of texts with AI













library(ggpubr)
p6<-ggarrange(p1,p2,p3,p4,p5, ncol = 5, nrow = 1,common.legend = TRUE,legend="bottom",align = "h")
p6

#save image 7 inches * 25 inches
#then rearrange the positions of texts with AI

