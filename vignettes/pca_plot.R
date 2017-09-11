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

#' data taken from GSK and BLUEPRINT databases

load("/r_data/all_data.RData")
load("/r_data/new_blueprint_chip.RData")

#######################################################

#' extended from epiChoose vinettes by Dr. Aiden
#' MacNamara
#'
#' @author Nam Janjumratsang
#' @description This script plots the pca of all_data
#' and fine grained data sets and save them as .png
#' files. The functions and parts of the code are taken
#'  from Dr. Aiden McNamara's epiChoose package.
#'
#' Each cluster is labelled according to biology
#' literature such as Haematopoiesis studies.

#######################################################

#########################
####all data pca plot####
#########################

#' taking only the BLUEPRINT data
for(i in 1:4) {

  all_data[[i]]$annot <- all_data[[i]]$annot[(17:dim(all_data[[i]]$annot)[1]),]
  all_data[[i]]$res <- all_data[[i]]$res[(17:dim(all_data[[i]]$res)[1]),]

}

#' group labels from biology literatures (lineage trees)
group_labels = as.character(c(3,1,3,3,1,1,5,2,2,2,6,2,2,6,2,4,2,2,7,6,6,7,6,7,6,4,6,2,2,2,2,6,3,3))
single_labels = rownames(all_data[[2]]$res)

pca_data = prep_for_plot(all_data, annot_1=group_labels, annot_2=single_labels, marks=c("rna","H3K27ac", "H3K4me3", "H3K27me3"))

png(filename="all_data_pca.png", height=800, width=3200)
pcaAllData <- ggplot(pca_data, aes(x=x, y=y, color=annot_1)) + geom_point(size=5, shape=17) + theme_thesis() + facet_wrap(~mark, nrow=1)
pcaAllData +scale_fill_manual(values=wes_palette(n=7, name="GrandBudapest"))
dev.off()

#########################################
####finegrained_gsk_chip_filtered PCA####
#########################################

#' group labels from biology literatures (lineage trees)
fingrained_group_labels = as.character(c(4,4,4,4,4,4,4,4,3,4,4,4,1,1,3,1,3,3,2,2,1,3,2,2,3,2,1,2,1,3,2,1))
fingrained_group_single_labels = rownames(finegrained_gsk_chip_filtered[[2]]$res)
fingrained_pca_data = prep_for_plot(finegrained_gsk_chip_filtered, annot_1=fingrained_group_labels, annot_2=fingrained_group_single_labels, marks=c("H3K27ac", "H3K4me3", "H3K27me3"))

#' fine grained data plot
png(filename="finegrained_pca.png", height=800, width=3200)
pcaFineGrainedData <- ggplot(fingrained_pca_data, aes(x=x, y=y, color=annot_1)) + geom_point(size=5, shape=17) + theme_thesis() + facet_wrap(~mark, nrow=1)
pcaFineGrainedData + scale_color_manual(breaks = c("1", "2", "3", "4"), values=c("red", "blue", "green", "black"))
dev.off()
