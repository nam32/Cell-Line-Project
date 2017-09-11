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
#' @description   heatmap for the original data set and 
#' fine grained data set are drawn. Log10 is applied to
#'  the data in the order to better visualise the data. 
#' This script saves the images to png files which will
#' appear in the working directory

#######################################################

load("/r_data/all_data.RData")
marks=c("rna","H3K27ac", "H3K4me3", "H3K27me3")
sample_ix = c(17:50)

for(i in 1:length(all_data)) {
  
  all_data[[i]]$res = all_data[[i]]$res[sample_ix,]
  all_data[[i]]$annot = all_data[[i]]$annot[sample_ix,]
    
  dat = t(all_data[[i]]$res)
  
  remove_rows = which(apply(dat, 1, function(x) all(is.na(x)))) # remove samples with no data na = no data epigenetics with no data
  if(length(remove_rows)) {
    dat_na_rm = dat[-remove_rows,]
  } else {
    dat_na_rm = dat
  }
  
  dat_na_rm = dat_na_rm[,!apply(dat_na_rm, 2, function(x) sd(x)==0)] # remove regions with no variance (rows with 0 variance)
  dim(dat_na_rm)
  
  dat_na_rm = dat_na_rm[,!apply(is.na(dat_na_rm), 2, all)] # remove regions with no data (remove columns i.e. genes with no data)
  dim(dat_na_rm)
  
  a <- standardNormalization(dat_na_rm)
  #a <- dat_na_rm
  png(filename=paste("heatmap_all_data_standardNorm_",marks[i],".png", collapse=",", sep=""), height=4000, width=7000)
  heatmap.2(a, Rowv = TRUE, col=redgreen(100),Colv = TRUE, dendrogram = 'both', main = paste(marks[i]), tracecol=NA, margins = c(120, 6), breaks=seq(-1, 1, length.out=101), cexCol = 4)
  dev.off()

}

#####################################
####heatmap for fine grained data####
#####################################

########################
####data preparation####
########################

load("/r_data/new_blueprint_chip.RData")
prepped_data = prep_gsk_chip_filter(blueprint_chip)
marks <-c("H3K27ac", "H3K4me3", "H3K27me3")

for(i in 1:length(marks)) {
  dat = t(prepped_data[[i]]$res)
  
  remove_rows = which(apply(dat, 1, function(x) all(is.na(x)))) # remove samples with no data na = no data epigenetics with no data
  if(length(remove_rows)) {
    dat_na_rm = dat[-remove_rows,]
  } else {
    dat_na_rm = dat
  }
  
  dat_na_rm = dat_na_rm[,!apply(dat_na_rm, 2, function(x) sd(x)==0)] # remove regions with no variance (rows with 0 variance)
  dim(dat_na_rm)
  
  dat_na_rm = dat_na_rm[,!apply(is.na(dat_na_rm), 2, all)] # remove regions with no data (remove columns i.e. genes with no data)
  dim(dat_na_rm)
  
  dat = eval(as.name(paste0("new_cluster_data_", i)))
  a <- standardNormalization(dat)
  png(filename=paste0("heatmap_finegrained_standardNorm_",marks[i],".png"), height=4000, width=7000)
  heatmap.2(a, Rowv = TRUE, col=redgreen(100),Colv = TRUE, dendrogram = 'both', main = paste(marks[i]), tracecol=NA, margins = c(120, 6), breaks=seq(-1, 1, length.out=101), cexCol = 4)
  dev.off()

}
