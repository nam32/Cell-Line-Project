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
#' @description Two different methods of integrations
#' are done:
#' 1. combining the matrices of data together and
#'    performing heirarchical clustering on the
#'    resulting matrix
#' 2. using SNF (similarity network fusion) method
#'    provided by SNFtool package. The parameters of the
#'    method is then varied in the order to optimise the
#'    resulting data into the most similar tree to the
#'    synthetic tree as possible.
#' all different combinations of the 3 data matrices
#' were tried for all different integration methods and
#' the best results accumulates in a data frame called
#' integrationResults

########################################################

##################
####tree prep#####
##################
newdata <-list()
newdata$merge <- matrix(c(-1, -2,#1
                          1, -3,#2
                          2, -4,#3
                          3, -5,#4
                          4, -6,#5
                          5, -7,#6
                          6, -8,#7
                          7, -9,#8
                          8,-10,#9
                          9,-11,#10

                          #inflammatory
                          -12,-13,#11
                          11,-14,#12
                          12,-15,#13
                          13,-16,#14
                          14,-17,#15
                          15,-18,#16

                          #alternatively activated macrophage
                          -19,-20,#17
                          17,-21,#18
                          18,-22,#19
                          19,-23,#20
                          20,-24,#21
                          21,-25,#22

                          #macrophage
                          -26,-27,#23
                          23,-28,#24
                          24,-29,#25
                          25,-30,#26
                          26,-31,#27
                          27,-32,#28

                          #merge
                          16, 22,#29
                          29, 28,#30
                          30, 10#31

), nc=2, byrow=TRUE )

newdata$height <- myTreePrep(c(28,4,2,6,1,18))
newdata$order <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32)            # order of leaves(trivial if hand-entered)
newdata$labels <- c("C000S5_CD14-positive, CD16-negative classical monocyte", "C0010K_CD14-positive, CD16-negative classical monocyte", "C0011I_CD14-positive, CD16-negative classical monocyte", "C001UY_CD14-positive, CD16-negative classical monocyte", "C00264_CD14-positive, CD16-negative classical monocyte", "C00280_CD14-positive, CD16-negative classical monocyte", "C004SQ_CD14-positive, CD16-negative classical monocyte", "C005PS_CD14-positive, CD16-negative classical monocyte", "S000RD_CD14-positive, CD16-negative classical monocyte", "Primary_monocytes_Donor_E_monocyte", "Primary_monocytes_Donor_F_monocyte", "S001MJ_inflammatory macrophage", "S001S7_inflammatory macrophage", "S0022I_inflammatory macrophage","S007SK_inflammatory macrophage", "S00H6O_inflammatory macrophage", "S01F8K_inflammatory macrophage","S01H5I_inflammatory macrophage", "S00622_alternatively activated macrophage", "S006VI_alternatively activated macrophage", "S00BS4_alternatively activated macrophage", "S00C1H_alternatively activated macrophage", "S00FTN_alternatively activated macrophage", "S00T2L_alternatively activated macrophage", "S01FW9_alternatively activated macrophage", "C005VG_macrophage", "S001S7_macrophage", "S0022I_macrophage", "S00390_macrophage", "S00BHQ_macrophage", "S00DVR_macrophage", "S01F8K_macrophage")
class(newdata) <- "hclust"
newdata <- as.dendrogram(newdata)

png(filename=paste("fined_grained_dendrogram.png", collapse=",", sep=""), height=800, width=3200) ; par(mar=c(24, 4, 4, 2))
plot(newdata)
dev.off()

#################
####data prep####
#################

load("/r_data/new_blueprint_chip.RData")
prepped_data = prep_gsk_chip_filter(blueprint_chip)
marks <-c("H3K27ac", "H3K4me3", "H3K27me3")

#' adding data to 'integrationResults' dataframe for comparison
integrationResults <- as.data.frame(setNames(replicate(3,numeric(0), simplify = F), c("name", "cophenetic index", "baker's gamma index")))

dat_slice_1 = prepped_data
for(i in 1:length(marks)) {

  currentDataslice <- dat_slice_1[[i]]
  currentMarks <- marks[[i]]
  assign(paste0("nolog_new_cluster_data_", i), dataPrep(currentDataslice, currentMarks, i))

}

#' data heirarchical clustering without integration
nolog_raw_data_dend1 <- as.dendrogram(varclus(nolog_new_cluster_data_1, similarity="pearson", method = "ward.D2"))
nolog_raw_data_dend2 <- as.dendrogram(varclus(nolog_new_cluster_data_2, similarity="pearson", method = "ward.D2"))
nolog_raw_data_dend3 <- as.dendrogram(varclus(nolog_new_cluster_data_3, similarity="pearson", method = "ward.D2"))

clusterVec <- c("nolog_raw_data_dend1", "nolog_raw_data_dend2", "nolog_raw_data_dend3")

for(i in 1:length(clusterVec)){

  current <- eval(as.name(paste0(clusterVec[[i]])))
  cophen <- cor_cophenetic(newdata, current)
  bakers <- cor_bakers_gamma(newdata, current)
  integrationResults[nrow(integrationResults) + 1, ] <- c(paste0(clusterVec[[i]]),cophen,bakers)

  #' save cut tree data to .csv file
  current <- as.matrix(cutree(current, k = 4))
  write.csv(current, file = paste0(clusterVec[[i]],"_hierarchicalClustering_TreeCut_4clusters.csv"))
  current <- as.matrix(cutree(current, k = 5))
  write.csv(current, file = paste0(clusterVec[[i]],"_hierarchicalClustering_TreeCut_5clusters.csv"))

}

##################################################
####integration method (dataframe combination)####
##################################################
##################################################

#' hierarchical clustering of matrices combination
#' cut tree to form 4 clusters

##################################################

combine12 <- rbind(nolog_new_cluster_data_1, nolog_new_cluster_data_2)
combine12dend <- as.dendrogram(varclus(combine12, similarity="pearson", method = "ward.D2"))

combine13 <- rbind(nolog_new_cluster_data_1, nolog_new_cluster_data_3)
combine13dend <- as.dendrogram(varclus(combine13, similarity="pearson", method = "ward.D2"))

combine23 <- rbind(nolog_new_cluster_data_2, nolog_new_cluster_data_3)
combine23dend <- as.dendrogram(varclus(combine23, similarity="pearson", method = "ward.D2"))

combine123 <- rbind(nolog_new_cluster_data_1, nolog_new_cluster_data_2, nolog_new_cluster_data_3)
combine123dend <- as.dendrogram(varclus(combine123, similarity="pearson", method = "ward.D2"))

#' save dendrogram of combined matrix
png(filename="combine123dend.png", height=2500, width=4500, res = 300); par(mar=c(24,6,1,1))
plot(combine123dend)
rect.dendrogram(combine123dend, k=4)
dev.off()

combineVec <- c("combine12dend", "combine13dend", "combine23dend", "combine123dend")

#' adding data to the dataframe

#integrationResults <- as.data.frame(setNames(replicate(3,numeric(0), simplify = F), c("name", "cophenetic index", "baker's gamma index")))

for(i in 1:length(combineVec)){

  current <- eval(as.name(paste0(combineVec[[i]])))
  cophen <- cor_cophenetic(newdata, current)
  bakers <- cor_bakers_gamma(newdata, current)
  integrationResults[nrow(integrationResults) + 1, ] <- c(paste0(combineVec[[i]]),cophen,bakers)

  current <- as.matrix(cutree(current, k = 4))
  write.csv(current, file = paste0("hierarchical_clustering_treecut_4_clusters_",combineVec[[i]], ".csv"))

  current <- as.matrix(cutree(current, k = 5))
  write.csv(current, file = paste0("hierarchical_clustering_treecut_5_clusters_",combineVec[[i]], ".csv"))


}

##############################
####SNF integration method####
##############################

#' data preperation

new_nolog_new_cluster_data_1 <- t(nolog_new_cluster_data_1)
new_nolog_new_cluster_data_2 <- t(nolog_new_cluster_data_2)
new_nolog_new_cluster_data_3 <- t(nolog_new_cluster_data_3)

truelabel = c(matrix(1,32,1),matrix(2,32,1))

Data1 = standardNormalization(new_nolog_new_cluster_data_1)
Data2 = standardNormalization(new_nolog_new_cluster_data_2)
Data3 = standardNormalization(new_nolog_new_cluster_data_3)

Dist1 = dist2(as.matrix(Data1),as.matrix(Data1))
Dist2 = dist2(as.matrix(Data2),as.matrix(Data2))
Dist3 = dist2(as.matrix(Data3),as.matrix(Data3))

#'data integration

#' K          number of neighbors, usually (10~30)
#' alpha 	    hyperparameter, usually (0.3~0.8)
#T = 10       Number of Iterations, usually (10~20)
#C       			number of clusters

alphaVec = seq(0.3,0.8,+0.1)
C = 4

integrationTestResults <- as.data.frame(setNames(replicate(3,numeric(0), simplify = F), c("name", "cophenetic index", "baker's gamma index")))

#' testing to find the best parameters combinations by varying K,
#' alpha and T within the recommended range

for(k in 10:30){

  for(j in length(alphaVec)){

    for(t in 10:20){

      W1 = affinityMatrix(Dist1, k, alphaVec[[j]])
      W2 = affinityMatrix(Dist2, k, alphaVec[[j]])
      W3 = affinityMatrix(Dist3, k, alphaVec[[j]])

      cellTypes <- row.names(new_nolog_new_cluster_data_2)

      W23 <- SNF(list(W2,W3), k, t)
      W12 <- SNF(list(W1,W2), k, t)
      W13 <- SNF(list(W1,W3), k, t)
      W123 <- SNF(list(W13,W2), k, t)

      group23 = spectralClustering(W23, C) 	  # the final subtypes information
      group12 = spectralClustering(W12, C) 	  # the final subtypes information
      group13 = spectralClustering(W13, C) 	  # the final subtypes information
      group123 = spectralClustering(W123, C) 	# the final subtypes information
      #displayClusters(W, group)

      row.names(W23) <- cellTypes
      colnames(W23) <- c(1:32)
      row.names(W12) <- cellTypes
      colnames(W12) <- c(1:32)
      row.names(W13) <- cellTypes
      colnames(W13) <- c(1:32)
      row.names(W123) <- cellTypes
      colnames(W123) <- c(1:32)

      resultsW23 <- cbind(cellTypes, group23)
      resultsW23 <-(data.frame(resultsW23))
      resultsW23 <- resultsW23[with(resultsW23, order(resultsW23$group23)), ]
      resultsW12 <- cbind(cellTypes, group12)
      resultsW12 <-(data.frame(resultsW12))
      resultsW12 <- resultsW12[with(resultsW12, order(resultsW12$group12)), ]
      resultsW13 <- cbind(cellTypes, group13)
      resultsW13 <-(data.frame(resultsW13))
      resultsW13 <- resultsW13[with(resultsW13, order(resultsW13$group13)), ]
      resultsW123 <- cbind(cellTypes, group123)
      resultsW123 <-(data.frame(resultsW123))
      resultsW123 <- resultsW123[with(resultsW123, order(resultsW123$group123)), ]

      W_dend23 <- hclust(dist(W23), method = "ward.D2")
      W_dend12 <- hclust(dist(W12), method = "ward.D2")
      W_dend13 <- hclust(dist(W13), method = "ward.D2")
      W_dend123 <- hclust(dist(W123), method = "ward.D2")

      plot(W_dend23)
      plot(W_dend12)
      plot(W_dend13)
      plot(W_dend123)

      Ws = c("23", "12", "13", "123")

      kjt <- paste0("k",k, "j", j, "t", t)

      for(a in 1:length(Ws)){

        currentW <- eval(as.name(paste0("W_dend", Ws[[a]])))
        cophen <- cor_cophenetic(newdata, currentW)
        bakers <- cor_bakers_gamma(newdata, currentW)
        integrationTestResults[nrow(integrationTestResults) + 1, ] <- c(paste0("W_dend", Ws[[a]], ".", kjt),cophen,bakers)

      }

    }

  }

}

maxintegrationResults <- max(integrationTestResults$cophenetic.index, na.rm = TRUE)
maxintegrationResultsdataframe <- integrationTestResults[grep(paste0(maxintegrationResults), integrationTestResults$cophenetic.index), ]
integrationResults <- rbind(integrationResults, maxintegrationResultsdataframe)

save(integrationResults, file="integrationResults.RData")

#' spectral clustering on combination that gives the best scores after evaluating the data
#' k=10
#' j=6
#' t=10

C = 4
k=10
j=6
t=10

W1 = affinityMatrix(Dist1, k, j)
W2 = affinityMatrix(Dist2, k, j)
W3 = affinityMatrix(Dist3, k, j)

cellTypes <- row.names(new_nolog_new_cluster_data_2)

W23 <- SNF(list(W2,W3), k, t)
W12 <- SNF(list(W1,W2), k, t)
W13 <- SNF(list(W1,W3), k, t)
W123 <- SNF(list(W13,W2), k, t)

#' spectral clustering on the resulting data
group23 = spectralClustering(W23, C) 	  # the final subtypes information
group12 = spectralClustering(W12, C) 	  # the final subtypes information
group13 = spectralClustering(W13, C) 	  # the final subtypes information
group123 = spectralClustering(W123, C) 	# the final subtypes information
#displayClusters(W, group)

row.names(W23) <- cellTypes
colnames(W23) <- c(1:32)
row.names(W12) <- cellTypes
colnames(W12) <- c(1:32)
row.names(W13) <- cellTypes
colnames(W13) <- c(1:32)
row.names(W123) <- cellTypes
colnames(W123) <- c(1:32)

spectral_resultsW23 <- cbind(cellTypes, group23)
spectral_resultsW23 <-(data.frame(spectral_resultsW23))
spectral_resultsW23 <- spectral_resultsW23[with(spectral_resultsW23, order(spectral_resultsW23$group23)), ]
spectral_resultsW12 <- cbind(cellTypes, group12)
spectral_resultsW12 <-(data.frame(spectral_resultsW12))
spectral_resultsW12 <- spectral_resultsW12[with(spectral_resultsW12, order(spectral_resultsW12$group12)), ]
spectral_resultsW13 <- cbind(cellTypes, group13)
spectral_resultsW13 <-(data.frame(spectral_resultsW13))
spectral_resultsW13 <- spectral_resultsW13[with(spectral_resultsW13, order(spectral_resultsW13$group13)), ]
spectral_resultsW123 <- cbind(cellTypes, group123)
spectral_resultsW123 <-(data.frame(spectral_resultsW123))
spectral_resultsW123 <- spectral_resultsW123[with(spectral_resultsW123, order(spectral_resultsW123$group123)), ]

data <- as.dist(W123)
#row.names(data) <- cellTypes
#colnames(data) <- c(1:32)
png(filename=paste("SNF123dendrogram_4clusters.png", collapse=",", sep=""), height=2500, width=4500) ; par(mar=c(24, 4, 4, 2))
snfDend <- as.dendrogram(varclus(t(W123), method = "ward.D2", members = NULL))
plot(snfDend)
rect.dendrogram(snfDend, k=4)
dev.off()

png(filename=paste("SNF123dendrogram_5clusters.png", collapse=",", sep=""), height=2500, width=4500) ; par(mar=c(24, 4, 4, 2))
snfDend <- as.dendrogram(varclus(t(W123), method = "ward.D2", members = NULL))
plot(snfDend)
rect.dendrogram(snfDend, k=5)
dev.off()

#' save cut tree data to .csv file
snf123HirarchicaltreeCut_4clusters <- as.matrix(cutree(snfDend, k = 4))
write.csv(snf123HirarchicaltreeCut_4clusters, file = "snf123HirarchicaltreeCut_4clusters.csv")

snf123HirarchicaltreeCut_5clusters <- as.matrix(cutree(snfDend, k = 5))
write.csv(snf123HirarchicaltreeCut_5clusters, file = "snf123HirarchicaltreeCut_5clusters.csv")

####################################################
####spectral clustering for matrices combination####
####################################################

####4 clusters
##mark 1
nolog_raw_data_1_spectral_4clusters <- spectralClusterSummary(new_nolog_new_cluster_data_1, c(1:32), 4)

##mark 2
nolog_raw_data_2_spectral_4clusters <- spectralClusterSummary(new_nolog_new_cluster_data_2, c(1:32), 4)

##mark 3
nolog_raw_data_3_spectral_4clusters <- spectralClusterSummary(new_nolog_new_cluster_data_3, c(1:32), 4)

cellTypes <- row.names(new_nolog_new_cluster_data_2)
combineVec <- c("combine12", "combine13", "combine23", "combine123")
for(i in 1:length(combineVec)){

  current <- t(eval(as.name(paste0(combineVec[[i]]))))
  assign(paste0("nolog_spectral_4clusters_", combineVec[[i]]), spectralClusterSummary(current, c(1:32), 4))

}

####5 clusters
##mark 1
nolog_raw_data_1_spectral_5clusters <- spectralClusterSummary(new_nolog_new_cluster_data_1, c(1:32), 5)

##mark 2
nolog_raw_data_2_spectral_5clusters <- spectralClusterSummary(new_nolog_new_cluster_data_2, c(1:32), 5)

##mark 3
nolog_raw_data_3_spectral_5clusters <- spectralClusterSummary(new_nolog_new_cluster_data_3, c(1:32), 5)

cellTypes <- row.names(new_nolog_new_cluster_data_2)
combineVec <- c("combine12", "combine13", "combine23", "combine123")
for(i in 1:length(combineVec)){
  
  current <- t(eval(as.name(paste0(combineVec[[i]]))))
  assign(paste0("nolog_spectral_5clusters_", combineVec[[i]]), spectralClusterSummary(current, c(1:32), 5))
  
}

#' writing all spectral results to .csv files
namesVec = c("nolog_raw_data_1_spectral_4clusters","nolog_raw_data_2_spectral_4clusters", 
             "nolog_raw_data_3_spectral_4clusters", "nolog_spectral_4clusters_combine12", 
             "nolog_spectral_4clusters_combine13", "nolog_spectral_4clusters_combine23", 
             "nolog_spectral_4clusters_combine123", "spectral_resultsW23","spectral_resultsW12",
             "spectral_resultsW13", "spectral_resultsW123", "nolog_raw_data_1_spectral_5clusters",
             "nolog_raw_data_2_spectral_5clusters", "nolog_raw_data_3_spectral_5clusters",
             "nolog_spectral_5clusters_combine12", "nolog_spectral_5clusters_combine13", 
             "nolog_spectral_5clusters_combine23", "nolog_spectral_5clusters_combine123")


for(i in 1:length(namesVec)){

  write.csv(eval(as.name(paste0(namesVec[[i]]))), file = paste0(namesVec[[i]],".csv"))

}

########################
####icluster package####
########################
########################################################

#' although the project does not include whole data
#' analysis of icluster results due to limitations in
#' computer resources (memory space too small), we did
#' do a shorter cluster analysis of iCluster as shown
#' below.

########################################################

load("/r_data/new_blueprint_chip.RData")
marks <-c("H3K27ac", "H3K4me3", "H3K27me3")
data_slice_1 = prep_gsk_chip_filter(blueprint_chip)

for(i in 1:length(marks)) {

  dat = dat_slice_1[[i]]$res

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

  dat_slice_1[[i]]$res <- dat_na_rm
}

########################################################

#' full data iCluster analysis. lambda (lasso penalty
#' term) was found using cv.glmnet() function from glmnet
#' package

########################################################

trainingData <- c(4,4,4,4,4,4,4,4,3,4,4,4,1,1,3,1,3,3,2,2,1,3,2,2,3,2,1,2,1,3,2,1)

cvFullFit1 <- cv.glmnet(x = dat_slice_1[[1]]$res, y = trainingData)
response1 <- cvFit1$lambda.1se

cvFullFit2 <- cv.glmnet(x = dat_slice_1[[2]]$res, y = trainingData)
response2 <- cvFit2$lambda.1se

cvFullFit3 <- cv.glmnet(x = dat_slice_1[[3]]$res, y = trainingData)
response3 <- cvFit3$lambda.1se

#' using the lambda aquired in iCluster
iClusterFull <- iCluster(list(dat_slice_1[[1]]$res, dat_slice_1[[2]]$res, dat_slice_1[[2]]$res), k=4, lambda=c(response1,response2,response3))

cellTypes <- row.names(dat_slice_1[[1]]$res)
iClusterFullResults <- cbind(cellTypes, iClusterFull$clusters)
iClusterFullResults <- iClusterFullResults[with(iClusterFull, order(iClusterFull$clusters)), ]
save(iClusterFullResults, file="ClusterResultsiClusterFull.RData")

########################################################

#' iCluster analysis on a shorter version of data taking
#'  those with the largest differences in the the column
#' First we take 1000 values with the highest
#' differences for each data slice

########################################################

load("/r_data/new_blueprint_chip.RData")

min1 <-apply(dat_slice_1[[1]]$res, 1, min)
max1 <- apply(dat_slice_1[[1]]$res, 1, max)
dat_slice_1[[1]]$res <- dat_slice_1[[1]]$res[order(max1-min1, decreasing=TRUE),]
nolog_new_cluster_data_1short <- dat_slice_1[[1]]$res[, 1:1000]

min2 <-apply(dat_slice_1[[2]]$res, 1, min)
max2 <- apply(dat_slice_1[[2]]$res, 1, max)
dat_slice_1[[2]]$res <- dat_slice_1[[2]]$res[order(max2-min2, decreasing=TRUE),]
nolog_new_cluster_data_2short <- dat_slice_1[[2]]$res[, 1:1000]

min3 <-apply(dat_slice_1[[3]]$res, 1, min)
max3 <- apply(dat_slice_1[[3]]$res, 1, max)
dat_slice_1[[3]]$res <- dat_slice_1[[3]]$res[order(max3-min3, decreasing=TRUE),]
nolog_new_cluster_data_3short <- dat_slice_1[[3]]$res[, 1:1000]

#' using variances

variance1 <- colVar(dat_slice_1[[1]]$res)
dat_slice_1[[1]]$res <- dat_slice_1[[1]]$res[order(variance1, decreasing=TRUE),]
nolog_new_cluster_data_1short <- dat_slice_1[[1]]$res[, 1:1000]

variance2 <- colVar(dat_slice_1[[2]]$res)
dat_slice_1[[2]]$res <- dat_slice_1[[2]]$res[order(variance2, decreasing=TRUE),]
nolog_new_cluster_data_2short <- dat_slice_1[[2]]$res[, 1:1000]

variance3 <- colVar(dat_slice_1[[3]]$res)
dat_slice_1[[3]]$res <- dat_slice_1[[3]]$res[order(variance3, decreasing=TRUE),]
nolog_new_cluster_data_3short <- dat_slice_1[[3]]$res[, 1:1000]

#' we then find the best lambdas (lasso penalty terms)
trainingData <- c(4,4,4,4,4,4,4,4,3,4,4,4,1,1,3,1,3,3,2,2,1,3,2,2,3,2,1,2,1,3,2,1)
cvFit1 <- cv.glmnet(x = nolog_new_cluster_data_1short, y = trainingData)
response1 <- cvFit1$lambda.1se

cvFit2 <- cv.glmnet(x = nolog_new_cluster_data_2short, y = trainingData)
response2 <- cvFit2$lambda.1se

cvFit3 <- cv.glmnet(x = nolog_new_cluster_data_3short, y = trainingData)
response3 <- cvFit3$lambda.1se

#' we then run iCluster on the three data matrices using the lamda values from glmnet
#' 4 clusters
shortCluster_4clusters <- iCluster(list(nolog_new_cluster_data_1short, nolog_new_cluster_data_2short, nolog_new_cluster_data_3short), k=4, lambda=c(response1,response2,response3))
shortClusterResults_4clusters <- as.matrix(cbind(cellTypes, shortCluster_4clusters$clusters))
shortClusterResults_4clusters <- shortClusterResults_4clusters[with(shortCluster_4clusters, order(shortCluster_4clusters$clusters)), ]
save(shortClusterResults_4clusters, file="shortClusterResultsiCluster.RData")
write.csv(shortClusterResults_4clusters, file="shortClusterResultsiCluster_4cluster.csv")

#' 5 clusters
shortCluster_5clusters <- iCluster(list(nolog_new_cluster_data_1short, nolog_new_cluster_data_2short, nolog_new_cluster_data_3short), k=5, lambda=c(response1,response2,response3))
cellTypes <- row.names(dat_slice_1[[1]]$res)
shortClusterResults_5clusters <- as.matrix(cbind(cellTypes, shortCluster_5clusters$clusters))
shortClusterResults_5clusters <- shortClusterResults_5clusters[with(shortCluster_5clusters, order(shortCluster_5clusters$clusters)), ]
#save(shortClusterResults, file="shortClusterResultsiCluster.RData")
write.csv(shortClusterResults_5clusters, file="shortClusterResultsiCluster_5cluster.csv")
