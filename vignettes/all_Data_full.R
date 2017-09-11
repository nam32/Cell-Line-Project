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
#' @description  hierarchical clustering of all_data 
#' and compare them with the full synthetic lineage 
#' tree

#######################################################

marks=c("rna","H3K27ac", "H3K4me3", "H3K27me3")
load("/r_data/all_data.RData")

###############################
####full lineage dendrogram####
###############################

lineage_tree <- list()
lineage_tree$merge <- matrix(c(-1,-2,
                     1, -3,
                     -4,-5,
                     3, -6,
                     4,-7,
                     5,-8,
                     6,-9,
                     7,-10,
                     8,-11,
                     9,-12,
                     10,-13,
                     11,-14,
                     12,-15,
                     13,2,
                     
                     -16,-17,
                     15,-18,
                     16,-19,
                     17,-20,

                     18,14,
                     
                     -21,-22,
                     
                     20,19,
                     
                     21,-23,
                     
                     -24,-25,
                     23,-26,
                     24,-27,
                     25,-28,
                     26,-29,
                     27,-30,
                     28,-31,
                     
                     29,22,
                     
                     -32,-33,
                     31,-34,
                     32,30

), nc=2, byrow=TRUE ) 
lineage_tree$height <- c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.6,0.4,0.6,0.6,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.6,0.8,0.8,1)    # define merge heights
lineage_tree$order <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34)            # order of leaves(trivial if hand-entered)
lineage_tree$labels <- c("C004SQ_CD14-positive, CD16-negative classical monocyte", "C005PS_CD14-positive, CD16-negative classical monocyte", "S000RD_CD14-positive, CD16-negative classical monocyte", "C005VG_macrophage", "S00390_macrophage", "S001MJ_inflammatory macrophage", "S001S7_macrophage", "S0022I_macrophage", "S00622_alternatively activated macrophage", "S0022I_inflammatory macrophage", "S001S7_inflammatory macrophage", "S00BS4_alternatively activated macrophage", "S00BHQ_macrophage", "S00C1H_alternatively activated macrophage", "S00DVR_macrophage", "C00184_mature neutrophil", "C12012_mature neutrophil", "C004GD_mature neutrophil", "BM060814_neutrophilic metamyelocyte", "BM060814_segmented neutrophil of bone marrow", "S002S3_erythroblast", "S002R5_erythroblast", "S004BT_CD34-negative, CD41-positive, CD42-positive megakaryocyte cell", "C0066P_CD8-positive, alpha-beta T cell", "S000RD_CD4-positive, alpha-beta T cell", "S007G7_CD4-positive, alpha-beta T cell", "S007DD_CD4-positive, alpha-beta T cell", "S008H1_CD4-positive, alpha-beta T cell", "S00C2F_CD8-positive, alpha-beta T cell", "S009W4_CD4-positive, alpha-beta T cell", "S0018A_CD4-positive, alpha-beta T cell","S00BJM_endothelial cell of umbilical vein (resting)","S00DCS_endothelial cell of umbilical vein (proliferating)","S00BJM_endothelial cell of umbilical vein (proliferating)")
class(lineage_tree) <- "hclust"
png(filename="not_simplified_all_data_dendrogram.png", height=1800, width=4800, res = 300)
plot(lineage_tree, hang = -1, cex = 0.4, main = "lineage tree")
lineage_tree_dend <- as.dendrogram(lineage_tree)
dev.off()

#######################################################
#' cluster all_data

sample_ix = c(17:50)

for(i in 1:length(all_data)) {
  all_data[[i]]$res = all_data[[i]]$res[sample_ix,]
  all_data[[i]]$annot = all_data[[i]]$annot[sample_ix,]
}

#par(mfrow=c(2,2))

for(i in 1:length(all_data)) {
  
  dat = all_data[[i]]$res
  
  remove_rows = which(apply(dat, 1, function(x) all(is.na(x)))) # remove samples with no data na = no data epigenetics with no data
  if(length(remove_rows)) {
    dat_na_rm = dat[-remove_rows,]
  } else {
    dat_na_rm = dat
  }
  
  dat_na_rm = dat_na_rm[,!apply(dat_na_rm, 2, function(x) sd(x)==0)] # remove regions with no variance (rows with 0 variance)
  dat_na_rm = dat_na_rm[,!apply(is.na(dat_na_rm), 2, all)] # remove regions with no data (remove columns i.e. genes with no data)
  #png(filename=paste("all_data_hclust_",marks[i],".png", collapse=",", sep=""), height=800, width=3200)
  clusters <- hclust(dist(dat_na_rm), method = "complete")
  assign(paste0("cluster_", i), clusters)
  #plot(clusters, main=marks[i])
  #dev.off()
}

#######################################################

#' compare the lineage tree with the clusters using both cophenetic and bakers gamma correlation
cor_cophenetic(lineage_tree_dend, cluster_2)
#0.4940695
cor_cophenetic(lineage_tree_dend, cluster_3)
#0.5443037
cor_cophenetic(lineage_tree_dend, cluster_4)
#0.8868998
cor_bakers_gamma(lineage_tree_dend, cluster_2)
#0.2729104
cor_bakers_gamma(lineage_tree_dend, cluster_3)
#0.5068829
cor_bakers_gamma(lineage_tree_dend, cluster_4)
#0.9720476

#' creating tanglegrams using Dendextend package to compare them visually
png(filename="not_simplified_pear_lineage_&_H3K27ac_tanglegram.png", height=1800, width=4800, res = 300)
tanglegram(lineage_tree_dend, cluster_2, main="lineage tree & H3K27ac")
dev.off()

png(filename="not_simplified_pear_lineage_&_H3K4me3_tanglegram.png", height=1800, width=4800, res = 300)
tanglegram(lineage_tree_dend, cluster_3, main="lineage tree & H3K4me3")
dev.off()

png(filename="not_simplified_pear_lineage_&_H3K27me3_tanglegram.png", height=1800, width=4800, res = 300)
tanglegram(lineage_tree_dend, cluster_4, main="lineage tree & H3K27me3")
dev.off()
