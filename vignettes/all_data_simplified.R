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

########################################################

#' @author Nam Janjumratsang
#' @description A test on simplified lineage dendrogram
#' and simplified heirarchical clustering marks.

########################################################

#####################################
####simplified lineage dendrogram####
#####################################

simLineage <- list()
simLineage$merge <- matrix(c(-1, -2,
                             -3, -4,
                             1, 2,
                             3, -5,
                             4, -6,
                             5, -7,
                             6, -8
                             
), nc=2, byrow=TRUE ) 
simLineage$height <- c(1, 1, 2, 3, 3, 4, 5) 
simLineage$order <- c(1,2,3,4,5,6,7,8)     
simLineage$labels <- c("monocyte", "macrophage", "neutrophil", "metamyelocyte neutrophil", "erythroyte", "megakaryocyte", "T lymphocyte", "endothelial cell of umbilical vein")
class(simLineage) <- "hclust"      
#plot(simLineage)                  
png(filename="simplified_lineage_tree.png", height=1500, width=3200, res = 300)
plot(simLineage, hang = -1, cex = 0.6, main = "lineage tree")
simLineageDend <- as.dendrogram(simLineage)
dev.off()

########################################################

#####################################
####simplified H3K4me3 dendrogram####
#####################################

simH3K4me3 <- list()  
simH3K4me3$merge <- matrix(c(-1, -2,
                             -3, -4,
                             1, 2,
                             3, -5,
                             4, -6,
                             -7,-8,
                             5, 6
), nc=2, byrow=TRUE ) 
simH3K4me3$height <- c(1, 1, 2, 3, 4, 4, 5)   
simH3K4me3$order <- c(8,7,6,5,4,3,2,1)     
simH3K4me3$labels <- c('monocyte', 'neutrophil', 'erythroyte', 'megakaryocyte', "metamyelocyte neutrophil", 'T lymphocyte', 'macrophage', 'endothelial cell of umbilical vein')    # labels of leaves
class(simH3K4me3) <- "hclust"        # make it an hclust object
#plot(simH3K4me3, center=TRUE, horiz=TRUE)             

#png(filename="simplified_H3K4me3_dendrogram.png", height=1500, width=3200, res = 300)
#plot(simH3K4me3, hang = -1, cex = 0.6, main="simplified H3K4me3 dendrogram")
simH3K4me3Dend <- as.dendrogram(simH3K4me3)
#plot(simH3K4me3, hang = -1, cex = 0.6)
#dev.off()

########################################################

#####################################
####simplified H3K27ac dendrogram####
#####################################

simH3K27ac <- list()
simH3K27ac$merge <- matrix(c(-1, -2,
                             1, -3,
                             -4,-5,
                             3, -6,
                             2,4,
                             5, -7,
                             6,-8
                             
), nc=2, byrow=TRUE ) 
simH3K27ac$height <- c(1, 2, 1, 2, 3, 4, 5)   
simH3K27ac$order <- c(8,7,6,5,4,3,2,1)         
simH3K27ac$labels <- c("erythroyte", "megakaryocyte", "endothelial cell of umbilical vein", "neutrophil", "metamyelocyte neutrophil", "monocyte", "T lymphocyte", "macrophage")    # labels of leaves
class(simH3K27ac) <- "hclust" 
#plot(simH3K27ac)
#png(filename="simplified_H3K27ac_dendrogram.png", height=1500, width=3200, res = 300)
#plot(simH3K27ac, hang = -1, cex = 0.6, main = "simplified H3K27ac dendrogram")

simH3K27acDend <- as.dendrogram(simH3K27ac)
#dev.off()

########################################################

######################################
####simplified H3K27me3 dendrogram####
######################################

simH3K27me3 <- list()  # initialize empty object
# define merging pattern: 
#    negative numbers are leaves, 
#    positive are merged clusters (defined by row number in $merge)
simH3K27me3$merge <- matrix(c(-1, -2,
                              1, -3,
                              -4, -5,
                              2, 3,
                              4, -6,
                              5, -7,
                              6, -8
), nc=2, byrow=TRUE ) 
simH3K27me3$height <- c(1, 2, 2, 3, 4, 5, 6)  
simH3K27me3$order <- c(8,7,6,5,4,3,2,1) 
simH3K27me3$labels <- c("neutrophil", "metamyelocyte neutrophil", "monocyte", "erythroyte", "megakaryocyte", "macrophage", "endothelial cell of umbilical vein", "T lymphocyte")    # labels of leaves
class(simH3K27me3) <- "hclust"        

#plot(simH3K27me3)                

#png(filename="simplified_H3K27me3_dendrogram.png", height=1500, width=3200, res = 300)
#plot(simH3K27me3, hang = -1, cex = 0.6, main = "simplified H3K27me3 dendrogram")
simH3K27me3Dend <- as.dendrogram(simH3K27me3)
#dev.off()

########################################################

#' using different method to clarify the similarities
#' between the lineage dendrogram and the simplified
#' clustered dendrograms

########################################################

##################################################
####compare lineage tree to H3K4me3 dendrogram####
##################################################

png(filename="lineage_tree_&_H3K4me3_tanglegram.png", height=600, width=1600, res = 300)
tanglegram(simLineageDend, simH3K4me3Dend, main="lineage tree & H3K4me3", "labels_cex"=10)
dev.off()
#Cophenetic correlation
cor_cophenetic(simLineageDend, simH3K4me3Dend)
#0.3938627
#Baker’s Gamma Index
cor_bakers_gamma(simLineageDend, simH3K4me3Dend)
#0.3799525


##################################################
####compare lineage tree to H3K27ac dendrogram####
##################################################

png(filename="lineage_tree_&_H3K27ac_tanglegram.png", height=600, width=1600, res = 300)
tanglegram(simLineageDend, simH3K27acDend, main="lineage tree & H3K27ac")
dev.off()
#Cophenetic correlation
cor_cophenetic(simLineageDend, simH3K27acDend)
#0.04545455
#Baker’s Gamma Index
cor_bakers_gamma(simLineageDend, simH3K27acDend)
#0.03144962

###################################################
####compare lineage tree to H3K27me3 dendrogram####
###################################################

png(filename="lineage_tree_&_H3K27me3_tanglegram.png", height=600, width=1600, res = 300)
tanglegram(simLineageDend, simH3K27me3Dend, main="lineage tree & H3K27me3")
dev.off()
#Cophenetic correlation
cor_cophenetic(simLineageDend, simH3K27me3Dend)
#0.7235668
#Baker’s Gamma Index
cor_bakers_gamma(simLineageDend, simH3K27me3Dend)
#0.7128041

##############################################
####compare H3K27ac to H3K27me3 dendrogram####
##############################################

png(filename="H3K27ac_&_H3K27me3_tanglegram.png", height=600, width=1600, res = 300)
tanglegram(simLineageDend, simH3K4me3Dend, main="H3K27ac & H3K27me3")
dev.off()
#Cophenetic correlation
cor_cophenetic(simH3K27acDend, simH3K27me3Dend)
#0.5621336
cor_bakers_gamma(simH3K27acDend, simH3K27me3Dend)
#0.5406116

#############################################
####compare H3K4me3 to H3K27ac dendrogram####
#############################################

png(filename="H3K4me3_&_H3K27ac_tanglegram.png", height=600, width=1600, res = 300)
tanglegram(simLineageDend, simH3K4me3Dend, main="H3K4me3 & H3K27ac")
dev.off()
#Cophenetic correlation
cor_cophenetic(simH3K4me3Dend, simH3K27acDend)
#0.5731824
cor_bakers_gamma(simH3K4me3Dend, simH3K27acDend)
#0.5298419

##############################################
####compare H3K4me3 to H3K27me3 dendrogram####
##############################################

png(filename="H3K4me3_&_H3K27me3_tanglegram.png", height=600, width=1600, res = 300)
tanglegram(simLineageDend, simH3K4me3Dend, main="H3K4me3 & H3K27me3")
dev.off()
#Cophenetic correlation
cor_cophenetic(simH3K4me3Dend, simH3K27me3Dend)
#0.6850908
cor(cophenetic(simH3K4me3Dend), cophenetic(simH3K27me3Dend))
#0.888374
cor_bakers_gamma(simH3K4me3Dend, simH3K27me3Dend)
#0.611768
