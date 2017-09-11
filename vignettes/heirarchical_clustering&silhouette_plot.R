require(glmnet)
require(SNFtool)
require(iClusterPlus)
require(treeOptAndComTools)
require(epiChoose)
require(tidyverse)
require(GenomicRanges)
require(biomaRt)
require(SummarizedExperiment)
require(devtools)
require(Rtsne)
require(factoextra)
require(Hmisc)
require(dendextend)
require(SNFtool)
require(clusternomics)
require(pheatmap)
require(heatmap.plus)
require(ggplot2)
require(gplots)
require(limma)
require(minfi)
require(wesanderson)
require(viridis)
load_all()

#######################################################

#' @author Nam Janjumratsang
#' @description hierarchical clustering of all_data and
#' the Silhouette plot for both first and second set
#' of data

#######################################################

###############################
####all data BLUEPRINT only####
###############################

load("/r_data/all_data.RData")

#' picking only the BLUEPRINT samples
dat_slice_1 = all_data

sample_ix = c(17:50)

for(i in 1:length(dat_slice_1)) {
  dat_slice_1[[i]]$res = dat_slice_1[[i]]$res[sample_ix,]
  dat_slice_1[[i]]$annot = dat_slice_1[[i]]$annot[sample_ix,]
}

#######################################################

marks=c("rna","H3K27ac", "H3K4me3", "H3K27me3")

for(i in 2:length(dat_slice_1)) {

  dat = dat_slice_1[[i]]$res

  remove_rows = which(apply(dat, 1, function(x) all(is.na(x)))) # remove samples with no data na = no data epigenetics with no data
  if(length(remove_rows)) {
    dat_na_rm = dat[-remove_rows,]
  } else {
    dat_na_rm = dat
  }

  dat_na_rm = dat_na_rm[,!apply(dat_na_rm, 2, function(x) sd(x)==0)] # remove regions with no variance (rows with 0 variance)
  #dim(dat_na_rm)
  dat_na_rm = dat_na_rm[,!apply(is.na(dat_na_rm), 2, all)] # remove regions with no data (remove columns i.e. genes with no data)
  #dim(dat_na_rm)

  png(filename=paste("all_data_hclust_",marks[i],".png", collapse=",", sep=""), height=1000, width=1500) ; par(mar=c(30, 4, 4, 2))
  clusters <- as.dendrogram(hclust(dist(dat_na_rm), method = "complete"))
  plot(clusters, main=marks[i])
  #plot( varclus(dat_na_rm, similarity="pearson") )
  rect.dendrogram(clusters, k=7)
  dev.off()

  #Silhouette plot
  png(filename=paste("silhouette_",marks[i],".png", collapse=",", sep=""), height=800, width=3200)
  silhouette <- fviz_silhouette(clusters)
  plot(silhouette, main=marks[i])
  dev.off()

  #' write tree cut clustering to a file
  current <- as.matrix(cutree(clusters, k = 7))
  write.csv(current, file = paste0("alldata_TreeCut",marks[[i]],".csv"))
  
}

#######################################################

#' after investigating into which clustering method
#' gives the best results.

#######################################################

#########################
####fine grained data####
#########################

marks=c("H3K27ac", "H3K4me3", "H3K27me3")
load("/r_data/new_blueprint_chip.RData")
fine_grained_data = prep_gsk_chip_filter(blueprint_chip)

for(i in 1:length(fine_grained_data)) {

  dat = fine_grained_data[[i]]$res

  remove_rows = which(apply(dat, 1, function(x) all(is.na(x)))) # remove samples with no data na = no data epigenetics with no data
  if(length(remove_rows)) {
    dat_na_rm = dat[-remove_rows,]
  } else {
    dat_na_rm = dat
  }

  dat_na_rm = dat_na_rm[,!apply(dat_na_rm, 2, function(x) sd(x)==0)] # remove regions with no variance (rows with 0 variance)
  dat_na_rm = dat_na_rm[,!apply(is.na(dat_na_rm), 2, all)] # remove regions with no data (remove columns i.e. genes with no data)

  png(filename=paste("finegrained_data_varclus_",marks[i],".png", collapse=",", sep=""), height=1000, width=1500) ; par(mar=c(30, 4, 4, 2))
  clusters <- as.dendrogram(varclus(t(dat_na_rm), similarity="spearman", method = "ward.D2"))
  plot(clusters)
  #plot( varclus(dat_na_rm, similarity="pearson") )
  rect.dendrogram(clusters, k=4)
  dev.off()

  #Silhouette plot
  png(filename=paste("silhouette_",marks[i],".png", collapse=",", sep=""), height=800, width=3200)
  silhouette <- fviz_silhouette(clusters)
  plot(silhouette, main=marks[i])
  dev.off()

  #' clusterings are evaluated in intergrationMethods.R
  
}

