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
#' @description  loops through the whole synthetic tree
#' of the fine grained data with varying heights for 
#' each level in the order to find optimal height for 
#' this particular data set. This is validated using 
#' H3K27ac3 which was hypothesized to give the best 
#' clusters. 3 differently preprocessing of data are 
#' used (log2(), log10() and the raw data). The results 
#' are then plotted and the heights that give the best 
#' scores for each type are stored in a dataframe.

########################################################

load("/r_data/new_blueprint_chip.RData")
marksVec <- c("H3K27ac", "H3K4me3", "H3K27me3")
correlationsVec <- c("spearman", "pearson")
linkagesVec <- c("single", "complete", "average", "mcquitty", "median", "centroid", "ward.D2")
numbers=c(1: length(marks))

######################################################
########################no log########################
######################################################

#' data preperation

dat_slice_1 = prep_gsk_chip_filter(blueprint_chip)

for(i in 1:length(marks)) {
  
  currentDataslice <- dat_slice_1[[i]]
  currentMarks <- marks[[i]]
  assign(paste0("nolog_new_cluster_data_", i), dataPrep(currentDataslice, currentMarks, i))
  
}

heightsTestNolog <- as.data.frame(setNames(replicate(8,numeric(0), simplify = F), c("mark", "corr", "linkage", "i", "j", "k", "cophenetic value", "baker's gamma value")))
for(x in 1:length(marksVec)){
  
  for (y in 1:length(correlationsVec)){
    
    for (z in 1:length(linkagesVec)){
      
      for (i in 1:5){
        
        for (j in 6:11){
          
          #for (k in 47:60){
          for (k in 17:22){
            
            newdata <- list()
            
            newdata$height <- myTreePrep(c(28,i,2,j,1,k))
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
            
            newdata$order <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32)            # order of leaves(trivial if hand-entered)
            newdata$labels <- c("C000S5_CD14-positive, CD16-negative classical monocyte", "C0010K_CD14-positive, CD16-negative classical monocyte", "C0011I_CD14-positive, CD16-negative classical monocyte", "C001UY_CD14-positive, CD16-negative classical monocyte", "C00264_CD14-positive, CD16-negative classical monocyte", "C00280_CD14-positive, CD16-negative classical monocyte", "C004SQ_CD14-positive, CD16-negative classical monocyte", "C005PS_CD14-positive, CD16-negative classical monocyte", "S000RD_CD14-positive, CD16-negative classical monocyte", "Primary_monocytes_Donor_E_monocyte", "Primary_monocytes_Donor_F_monocyte", "S001MJ_inflammatory macrophage", "S001S7_inflammatory macrophage", "S0022I_inflammatory macrophage","S007SK_inflammatory macrophage", "S00H6O_inflammatory macrophage", "S01F8K_inflammatory macrophage","S01H5I_inflammatory macrophage", "S00622_alternatively activated macrophage", "S006VI_alternatively activated macrophage", "S00BS4_alternatively activated macrophage", "S00C1H_alternatively activated macrophage", "S00FTN_alternatively activated macrophage", "S00T2L_alternatively activated macrophage", "S01FW9_alternatively activated macrophage", "C005VG_macrophage", "S001S7_macrophage", "S0022I_macrophage", "S00390_macrophage", "S00BHQ_macrophage", "S00DVR_macrophage", "S01F8K_macrophage")
            class(newdata) <- "hclust"
            newdata <- as.dendrogram(newdata)
            
            compareToDend <- as.dendrogram(varclus(eval(as.name(paste0("nolog_new_cluster_data_",x))), similarity=paste0(corVec[[y]]), method=paste0(linkagesVec[[z]])))
            
            cophen <- cor_cophenetic(compareToDend, newdata)
            
            bakers <- cor_bakers_gamma(compareToDend, newdata)
            
            heightsTestNolog[nrow(heightsTestNolog) + 1, ] <- c(paste0(marks[[x]]), corVec[[y]], linkagesVec[[z]], i, j, k, paste0(cophen), paste0(bakers))
            
          }
          
        }
        
      }
      
    }
    
  }
  
}

save(heightsTestNolog, file="heightsTestNolog.RData")

#' plotting graphs and comparing data

cophenetictest <- as.numeric(heightsTestNolog$cophenetic.value)
bakertest <- as.numeric(heightsTestNolog$baker.s.gamma.value)

png(filename=paste("copheneticvsbakerheightsTestnolog.png", collapse=",", sep=""), height=1000, width=2000)
palette(c("red","blue"))
plot(cophenetictest, bakertest, col = factor(heightsTestNolog$corr), cex=2, pch=16)
text(cophenetictest, bakertest, heightsTestNolog$corr, cex=1, pos = 3, col="red")
dev.off()

#heightsTestNologOnlycerfeatures <- heightsTestNolog[grep("spearman",heightsTestNolog$corr, ""), ]
heightsTestNologOnlycerfeatures <- heightsTestNolog[grep("",heightsTestNolog$corr, ""), ]

nologbefore23 <- heightsTestNologOnlycerfeatures[as.numeric(heightsTestNologOnlycerfeatures$k) <23,]
nologafter23 <- heightsTestNologOnlycerfeatures[as.numeric(heightsTestNologOnlycerfeatures$k) >23,]
sixtynolog <- vectorNum("60", "nologafter23", "")
nologafter23 <- cbind(nologafter23, sixtynolog)
png(filename=paste("heightsTestNolog(k60).png", collapse=",", sep=""), height=800, width=5000)
palette(c("green", "orange"))
plot(nologafter23$cophenetic.value, col = factor(nologafter23$X60), cex=2, pch=16)
#text(nologafter23$cophenetic.value, nologafter23$corr, cex=1, pos = 3, col="red")
grid(nx = NULL, ny = NULL)
dev.off()

#' tagging a specific number with a different color
twentytwonolog <- vectorNum("22", "nologbefore23", "")
nologbefore23 <- cbind(nologbefore23, twentytwonolog)
png(filename=paste("heightsTestNolog(k22).png", collapse=",", sep=""), height=800, width=5000)
palette(c("green", "orange"))
plot(nologbefore23$cophenetic.value, col = factor(nologbefore23$X22), cex=2, pch=16)
text(nologbefore23$cophenetic.value, nologbefore23$corr, cex=1, pos = 3, col="red")
grid(nx = NULL, ny = NULL)
dev.off()

######################################################
#########################LOG2#########################
######################################################

#' data preperation

dat_slice_1 = prep_gsk_chip_filter(blueprint_chip)

for(i in 1:length(marks)) {
  
  currentDataslice <- dat_slice_1[[i]]
  currentMarks <- marks[[i]]
  assign(paste0("log2_new_cluster_data_", i), dataPrepLog2(currentDataslice, currentMarks, i))
  
}

heightsTestLog2 <- as.data.frame(setNames(replicate(8,numeric(0), simplify = F), c("mark", "corr", "linkage", "i", "j", "k", "cophenetic value", "baker's gamma value")))

for(x in 1:length(marksVec)){
  
  for (y in 1:length(correlationsVec)){
    
    for (z in 1:length(linkagesVec)){
      
      for (i in 1:5){
        
        for (j in 6:11){
          
          #for (k in 47:60){
            
          for (k in 17:22){
            
            
            #tree with varying for loop heights
            newdata <- list()
            
            newdata$height <- myTreePrep(c(28,i,2,j,1,k))
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
            
            newdata$order <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32)            # order of leaves(trivial if hand-entered)
            newdata$labels <- c("C000S5_CD14-positive, CD16-negative classical monocyte", "C0010K_CD14-positive, CD16-negative classical monocyte", "C0011I_CD14-positive, CD16-negative classical monocyte", "C001UY_CD14-positive, CD16-negative classical monocyte", "C00264_CD14-positive, CD16-negative classical monocyte", "C00280_CD14-positive, CD16-negative classical monocyte", "C004SQ_CD14-positive, CD16-negative classical monocyte", "C005PS_CD14-positive, CD16-negative classical monocyte", "S000RD_CD14-positive, CD16-negative classical monocyte", "Primary_monocytes_Donor_E_monocyte", "Primary_monocytes_Donor_F_monocyte", "S001MJ_inflammatory macrophage", "S001S7_inflammatory macrophage", "S0022I_inflammatory macrophage","S007SK_inflammatory macrophage", "S00H6O_inflammatory macrophage", "S01F8K_inflammatory macrophage","S01H5I_inflammatory macrophage", "S00622_alternatively activated macrophage", "S006VI_alternatively activated macrophage", "S00BS4_alternatively activated macrophage", "S00C1H_alternatively activated macrophage", "S00FTN_alternatively activated macrophage", "S00T2L_alternatively activated macrophage", "S01FW9_alternatively activated macrophage", "C005VG_macrophage", "S001S7_macrophage", "S0022I_macrophage", "S00390_macrophage", "S00BHQ_macrophage", "S00DVR_macrophage", "S01F8K_macrophage")
            class(newdata) <- "hclust"
            newdata <- as.dendrogram(newdata)
            
            compareToDend <- as.dendrogram(varclus(eval(as.name(paste0("log2_new_cluster_data_",x))), similarity=paste0(corVec[[y]]), method=paste0(linkagesVec[[z]])))
            
            cophen <- cor_cophenetic(compareToDend, newdata)
            
            bakers <- cor_bakers_gamma(compareToDend, newdata)
            
            heightsTestLog2[nrow(heightsTestLog2) + 1, ] <- c(paste0(marks[[x]]), corVec[[y]], linkagesVec[[z]], i, j, k, paste0(cophen), paste0(bakers))
            
          }
          
        }
        
      }
      
    }
    
  }
  
}

save(heightsTestLog2, file="heightsTestLog2.RData")

cophenetictest <- as.numeric(heightsTestLog2$cophenetic.value)
bakertest <- as.numeric(heightsTestLog2$baker.s.gamma.value)

png(filename=paste("copheneticvsbakerheightsTestLog2.png", collapse=",", sep=""), height=1000, width=2000)
palette(c("red","blue"))
plot(cophenetictest, bakertest, col = factor(heightsTestLog2$corr), cex=2, pch=16)
text(cophenetictest, bakertest, heightsTestLog2$corr, cex=1, pos = 3, col="red")
dev.off()

heightsTestLog2Onlycerfeatures <- heightsTestLog2[grep("",heightsTestLog2$corr, ""), ]

log2before23 <- heightsTestLog2Onlycerfeatures[as.numeric(heightsTestLog2Onlycerfeatures$k) <23,]
log2after23 <- heightsTestLog2Onlycerfeatures[as.numeric(heightsTestLog2Onlycerfeatures$k) >23,]
sixtylog2 <- vectorNum("60", "log2after23", "")
log2after23 <- cbind(log2after23, sixtylog2)
png(filename=paste("heightsTestLog2(k60).png", collapse=",", sep=""), height=800, width=5000)
palette(c("green", "orange"))
plot(log2after23$cophenetic.value, col = factor(log2after23$X60), cex=2, pch=16)
#text(log2after23$cophenetic.value, log2after23$corr, cex=1, pos = 3, col="red")
grid(nx = NULL, ny = NULL)
dev.off()

twentytwolog2 <- vectorNum("22", "log2before23", "")
log2before23 <- cbind(log2before23, twentytwolog2)
png(filename=paste("heightsTestLog2(k22).png", collapse=",", sep=""), height=800, width=5000)
palette(c("green", "orange"))
plot(log2before23$cophenetic.value, col = factor(log2before23$X22), cex=2, pch=16)
text(log2before23$cophenetic.value, log2before23$corr, cex=1, pos = 3, col="red")
grid(nx = NULL, ny = NULL)
dev.off()

#####################################################
########################LOG10########################
#####################################################

#' data preperation

dat_slice_1 = prep_gsk_chip_filter(blueprint_chip)

for(i in 1:length(marks)) {
  
  currentDataslice <- dat_slice_1[[i]]
  currentMarks <- marks[[i]]
  assign(paste0("log10_new_cluster_data_", i), dataPrepLog10(currentDataslice, currentMarks, i))
  
}

heightsTestLog10 <- as.data.frame(setNames(replicate(8,numeric(0), simplify = F), c("mark", "corr", "linkage", "i", "j", "k", "cophenetic value", "baker's gamma value")))

for(x in 1:length(marksVec)){
  
  for (y in 1:length(correlationsVec)){
    
    for (z in 1:length(linkagesVec)){
      
      for (i in 1:5){
        
        for (j in 6:11){
          
          #for (k in 47:60){
          
          for (k in 17:22){
            
            #tree with varying for loop heights
            newdata <- list()
            
            newdata$height <- myTreePrep(c(28,i,2,j,1,k))
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
            
            newdata$order <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32)            # order of leaves(trivial if hand-entered)
            newdata$labels <- c("C000S5_CD14-positive, CD16-negative classical monocyte", "C0010K_CD14-positive, CD16-negative classical monocyte", "C0011I_CD14-positive, CD16-negative classical monocyte", "C001UY_CD14-positive, CD16-negative classical monocyte", "C00264_CD14-positive, CD16-negative classical monocyte", "C00280_CD14-positive, CD16-negative classical monocyte", "C004SQ_CD14-positive, CD16-negative classical monocyte", "C005PS_CD14-positive, CD16-negative classical monocyte", "S000RD_CD14-positive, CD16-negative classical monocyte", "Primary_monocytes_Donor_E_monocyte", "Primary_monocytes_Donor_F_monocyte", "S001MJ_inflammatory macrophage", "S001S7_inflammatory macrophage", "S0022I_inflammatory macrophage","S007SK_inflammatory macrophage", "S00H6O_inflammatory macrophage", "S01F8K_inflammatory macrophage","S01H5I_inflammatory macrophage", "S00622_alternatively activated macrophage", "S006VI_alternatively activated macrophage", "S00BS4_alternatively activated macrophage", "S00C1H_alternatively activated macrophage", "S00FTN_alternatively activated macrophage", "S00T2L_alternatively activated macrophage", "S01FW9_alternatively activated macrophage", "C005VG_macrophage", "S001S7_macrophage", "S0022I_macrophage", "S00390_macrophage", "S00BHQ_macrophage", "S00DVR_macrophage", "S01F8K_macrophage")
            class(newdata) <- "hclust"
            newdata <- as.dendrogram(newdata)
            
            compareToDend <- as.dendrogram(varclus(eval(as.name(paste0("log10_new_cluster_data_",x))), similarity=paste0(corVec[[y]]), method=paste0(linkagesVec[[z]])))
            
            cophen <- cor_cophenetic(compareToDend, newdata)
            
            bakers <- cor_bakers_gamma(compareToDend, newdata)
            
            heightsTestLog10[nrow(heightsTestLog10) + 1, ] <- c(paste0(marks[[x]]), corVec[[y]], linkagesVec[[z]], i, j, k, paste0(cophen), paste0(bakers))
            
          }
          
        }
        
      }
      
    }
    
  }
  
}

save(heightsTestLog10, file="heightsTestLog10.RData")

cophenetictest <- as.numeric(heightsTestLog10$cophenetic.value)
bakertest <- as.numeric(heightsTestLog10$baker.s.gamma.value)

png(filename=paste("copheneticvsbakerheightsTestLog10.png", collapse=",", sep=""), height=1000, width=2000)
palette(c("red","blue"))
plot(cophenetictest, bakertest, col = factor(heightsTestLog10$corr), cex=2, pch=16)
text(cophenetictest, bakertest, heightsTestLog10$corr, cex=1, pos = 3, col="red")
dev.off()

heightsTestLog10Onlycerfeatures <- heightsTestLog10[grep("",heightsTestLog10$corr, ""), ]

log10before23 <- heightsTestLog10Onlycerfeatures[as.numeric(heightsTestLog10Onlycerfeatures$k) <23,]
log10after23 <- heightsTestLog10Onlycerfeatures[as.numeric(heightsTestLog10Onlycerfeatures$k) >23,]
sixtylog10 <- vectorNum("60", "log10after23", "")
log10after23 <- cbind(log10after23, sixtylog10)
png(filename=paste("heightsTestLog10(k60).png", collapse=",", sep=""), height=800, width=5000)
palette(c("green", "orange"))
plot(log10after23$cophenetic.value, col = factor(log10after23$X60), cex=2, pch=16)
#text(log10after23$cophenetic.value, log10after23$corr, cex=1, pos = 3, col="red")
grid(nx = NULL, ny = NULL)
dev.off()

twentytwolog10 <- vectorNum("22", "log10before23", "")
log10before23 <- cbind(log10before23, twentytwolog10)
png(filename=paste("heightsTestLog10(k22).png", collapse=",", sep=""), height=800, width=5000)
palette(c("green", "orange"))
plot(log10before23$cophenetic.value, col = factor(log10before23$X22), cex=2, pch=16)
text(log10before23$cophenetic.value, log10before23$corr, cex=1, pos = 3, col="red")
grid(nx = NULL, ny = NULL)
dev.off()

######################################################

#' identifying max value of each data types and store
#' them in a dataframe called optimalVal.

######################################################

#' creating results data frame
optimalVal <- as.data.frame(setNames(replicate(10,numeric(0), simplify = F), c("mark", "corr", "linkage", "i", "j", "k", "cophenetic value", "baker's gamma value", "preprocessing")))

######################################################
########################no log########################
######################################################

maxCorrnolog <- max(heightsTestNolog$cophenetic.value, na.rm = TRUE)
maxCorrnologdataframe <-heightsTestNolog[grep(paste0(maxCorrnolog), heightsTestNolog$cophenetic.value), ]
#pearson 0.945765257746349
#spearsman 0.88877892069262
nameVec <- as.data.frame(setNames(replicate(1,numeric(0), simplify = F), c("preprocessing")))
for(i in 1:length(maxCorrnologdataframe$mark)){
  nameVec[nrow(nameVec) + 1, ] <- c(paste0("not logged"))
}
maxCorrnologdataframe <- cbind(maxCorrnologdataframe, nameVec)
optimalVal<-rbind(optimalVal, maxCorrnologdataframe)

######################################################
#########################LOG2#########################
######################################################

maxCorrLog2 <- max(heightsTestLog2$cophenetic.value, na.rm = TRUE)
maxCorrLog2dataframe <- heightsTestLog2[grep(paste0(maxCorrLog2), heightsTestLog2$cophenetic.value), ]
#pearson 0.89145627407426
#spearman0.909243010430466

nameVec <- as.data.frame(setNames(replicate(1,numeric(0), simplify = F), c("preprocessing")))
for(i in 1:length(maxCorrLog2dataframe$mark)){
  nameVec[nrow(nameVec) + 1, ] <- c(paste0("Log2"))
}
maxCorrLog2dataframe <- cbind(maxCorrLog2dataframe, nameVec)
optimalVal<-rbind(optimalVal, maxCorrLog2dataframe)

#####################################################
########################LOG10########################
#####################################################

maxCorrLog10 <- max(heightsTestLog10$cophenetic.value, na.rm = TRUE)
maxCorrLog10dataframe <- heightsTestLog10[grep(paste0(maxCorrLog10), heightsTestLog10$cophenetic.value), ]

nameVec <- as.data.frame(setNames(replicate(1,numeric(0), simplify = F), c("preprocessing")))
for(i in 1:length(maxCorrLog10dataframe$mark)){
  nameVec[nrow(nameVec) + 1, ] <- c(paste0("Log10"))
}
maxCorrLog10dataframe <- cbind(maxCorrLog10dataframe, nameVec)
optimalVal<-rbind(optimalVal, maxCorrLog10dataframe)

#save & write to csv
save(optimalVal, file="optimalVal.RData")
write.csv(optimalVal, file = "optimalVal.csv")

############################
####separating the marks####
############################
######################################################

nameVec <- as.data.frame(setNames(replicate(1,numeric(0), simplify = F), c("preprocessing")))
for(i in 1:length(heightsTestNolog$mark)){
  nameVec[nrow(nameVec) + 1, ] <- c(paste0("not logged"))
}
heightsTestNolog <- cbind(heightsTestNolog, nameVec)

nameVec <- as.data.frame(setNames(replicate(1,numeric(0), simplify = F), c("preprocessing")))
for(i in 1:length(heightsTestLog2$mark)){
  nameVec[nrow(nameVec) + 1, ] <- c(paste0("log2"))
}
heightsTestLog2 <- cbind(heightsTestLog2, nameVec)

nameVec <- as.data.frame(setNames(replicate(1,numeric(0), simplify = F), c("preprocessing")))
for(i in 1:length(heightsTestLog10$mark)){
  nameVec[nrow(nameVec) + 1, ] <- c(paste0("log10"))
}
heightsTestLog10 <- cbind(heightsTestLog10, nameVec)

######################################################

H3K27acMarkNolog <- heightsTestNolog[grep(paste0("H3K27ac"), heightsTestNolog$mark), ]
H3K27acMarkLog2 <- heightsTestLog2[grep(paste0("H3K27ac"), heightsTestLog2$mark), ]
H3K27acMarkLog10 <- heightsTestLog10[grep(paste0("H3K27ac"), heightsTestLog10$mark), ]

H3K27acMark <- H3K27acMarkNolog
H3K27acMark<-rbind(H3K27acMark, H3K27acMarkLog2)
H3K27acMark<-rbind(H3K27acMark, H3K27acMarkLog10)

maxH3K27acMark <- max(H3K27acMark$cophenetic.value, na.rm = TRUE)
maxForAllMark <- H3K27acMark[grep(paste0(maxH3K27acMark), H3K27acMark$cophenetic.value), ]

########################################################

H3K4me3MarkNolog <- heightsTestNolog[grep(paste0("H3K4me3"), heightsTestNolog$mark), ]
H3K4me3MarkLog2 <- heightsTestLog2[grep(paste0("H3K4me3"), heightsTestLog2$mark), ]
H3K4me3MarkLog10 <- heightsTestLog10[grep(paste0("H3K4me3"), heightsTestLog10$mark), ]

H3K4me3Mark<-H3K4me3MarkNolog
H3K4me3Mark<-rbind(H3K4me3Mark, H3K4me3MarkLog2)
H3K4me3Mark<-rbind(H3K4me3Mark, H3K4me3MarkLog10)

maxH3K4me3Mark <- max(H3K4me3Mark$cophenetic.value, na.rm = TRUE)
temp <- H3K4me3Mark[grep(paste0(maxH3K4me3Mark), H3K4me3Mark$cophenetic.value), ]
maxForAllMark <- rbind(maxForAllMark, temp)

########################################################

H3K27me3MarkNolog <- heightsTestNolog[grep(paste0("H3K27me3"), heightsTestNolog$mark), ]
H3K27me3MarkLog2 <- heightsTestLog2[grep(paste0("H3K27me3"), heightsTestLog2$mark), ]
H3K27me3MarkLog10 <- heightsTestLog10[grep(paste0("H3K27me3"), heightsTestLog10$mark), ]

H3K27me3Mark<-H3K27me3MarkNolog
H3K27me3Mark<-rbind(H3K27me3Mark, H3K27me3MarkLog2)
H3K27me3Mark<-rbind(H3K27me3Mark, H3K27me3MarkLog10)

maxH3K27me3Mark <- max(H3K27me3Mark$cophenetic.value, na.rm = TRUE)
temp <- H3K27me3Mark[grep(paste0(maxH3K27me3Mark), H3K27me3Mark$cophenetic.value), ]
maxForAllMark <- rbind(maxForAllMark, temp)

save(maxForAllMark, file="maxForAllMark.RData")
write.csv(maxForAllMark, file = "maxForAllMark.csv")