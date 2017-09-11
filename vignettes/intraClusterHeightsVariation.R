require(tidyverse)
require(GenomicRanges)
require(biomaRt)
require(devtools)
require(SummarizedExperiment)
library(devtools)
library(Rtsne)
library(factoextra)
library(Hmisc)
library(dendextend)
library(SNFtool)
library(treeOptAndComTools)
library(epiChoose)
load("/r_data/new_blueprint_chip.RData")

########################################################

#' @author Nam Janjumratsang
#' @description  Trees are tested for the best
#' parameters (linkage and distance) for raw data The
#' results are plotted in scatter graph and the top
#' performing combination(s) is/are stored in .RData
#' form.

########################################################

#' data download and preperation

marks <-c("H3K27ac", "H3K4me3", "H3K27me3")
numbers=c(1: length(marks))
dat_slice_1 = prep_gsk_chip_filter(blueprint_chip)

#' tree preperation (tree must have this name in the order to be passed into marksLoop function or else error will occur)
newdataMerge <- matrix(c(-1, -2,#1
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

                         #merge all together
                         16, 22,#29
                         29, 28,#30
                         30, 10#31

), nc=2, byrow=TRUE )
#newdataHeight <- myTreePrep(c(28,2,2,6,1,20))
newdataHeight <- myTreePrep(c(28,4,2,6,1,18))
newdataOrder <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32)
newdataLabels <- c("C000S5_CD14-positive, CD16-negative classical monocyte", "C0010K_CD14-positive, CD16-negative classical monocyte", "C0011I_CD14-positive, CD16-negative classical monocyte", "C001UY_CD14-positive, CD16-negative classical monocyte", "C00264_CD14-positive, CD16-negative classical monocyte", "C00280_CD14-positive, CD16-negative classical monocyte", "C004SQ_CD14-positive, CD16-negative classical monocyte", "C005PS_CD14-positive, CD16-negative classical monocyte", "S000RD_CD14-positive, CD16-negative classical monocyte", "Primary_monocytes_Donor_E_monocyte", "Primary_monocytes_Donor_F_monocyte", "S001MJ_inflammatory macrophage", "S001S7_inflammatory macrophage", "S0022I_inflammatory macrophage","S007SK_inflammatory macrophage", "S00H6O_inflammatory macrophage", "S01F8K_inflammatory macrophage","S01H5I_inflammatory macrophage", "S00622_alternatively activated macrophage", "S006VI_alternatively activated macrophage", "S00BS4_alternatively activated macrophage", "S00C1H_alternatively activated macrophage", "S00FTN_alternatively activated macrophage", "S00T2L_alternatively activated macrophage", "S01FW9_alternatively activated macrophage", "C005VG_macrophage", "S001S7_macrophage", "S0022I_macrophage", "S00390_macrophage", "S00BHQ_macrophage", "S00DVR_macrophage", "S01F8K_macrophage")

##################
#####raw data#####
##################

for(i in 1:length(marks)) {

  currentDataslice <- dat_slice_1[[i]]
  currentMarks <- marks[[i]]
  assign(paste0("new_cluster_data_", i), dataPrep(currentDataslice, currentMarks, i))

}

marks <-c("H3K27ac", "H3K4me3", "H3K27me3")
correlationsVec <- c("spearman")
linkagesVec <- c("ward.D2")
finalDataFrame <- marksLoop(marks,correlationsVec,linkagesVec)
save(finalDataFrame, file="finalDataFrame(raw_data).RData")
write.csv(finalDataFrame, file="finalDataFrame(raw_data).csv")

#' find standard deviation of each mark
mark1cop <- sd(finalDataFrame[1:20,]$cophenetic.value)
mark2cop <- sd(finalDataFrame[21:40,]$cophenetic.value)
mark3cop <- sd(finalDataFrame[31:60,]$cophenetic.value)

mark1bak <- sd(finalDataFrame[1:20,]$baker.s.gamma.value)
mark2bak <- sd(finalDataFrame[21:40,]$baker.s.gamma.value)
mark3bak <- sd(finalDataFrame[31:60,]$baker.s.gamma.value)

#####################
#####log 10 data#####
#####################

for(i in 1:length(marks)) {

  currentDataslice <- dat_slice_1[[i]]
  currentMarks <- marks[[i]]
  assign(paste0("new_cluster_data_", i), dataPrepLog10(currentDataslice, currentMarks, i))
  dev.off()

}

#marks <-c("H3K27ac", "H3K4me3", "H3K27me3")
#correlationsVec <- c("spearman", "pearson")
#linkagesVec <- c("single", "complete", "average", "mcquitty", "median", "centroid", "ward.D2")
finalDataFrame <- marksLoop(marks,correlationsVec,linkagesVec)
save(finalDataFrame, file="finalDataFrame(log10).RData")
write.csv(finalDataFrame, file="finalDataFrame(log10).csv")

#####################
#####plot graphs#####
#####################
#load data
finalDataFrameRaw_data <- get(load("/finalDataFrame(raw_data).RData"))
finalDataFrameLog10 <- get(load("/finalDataFrame(log10).RData"))

#delete row containing certain elements
#finalDataFrameRaw_data <- finalDataFrameRaw_data[!finalDataFrameRaw_data$correlation.used == "ward.D", ]

################
####raw data####
################

cophenRaw <- finalDataFrameRaw_data$cophenetic.value
bakerRaw <- finalDataFrameRaw_data$baker.s.gamma.value
col1Raw <- finalDataFrameRaw_data$correlation.used
col2Raw <- finalDataFrameRaw_data$type.of.linkage.used
col3Raw <- finalDataFrameRaw_data$second.object

temp <- as.data.frame(setNames(replicate(3,numeric(0), simplify = F), c("correlation", "linkage", "mark")))


marks <-c("fine grained dataset H3K27ac cluster", "fine grained dataset H3K4me3 cluster", "fine grained dataset H3K27me3 cluster")
correlation <- c("single", "complete", "average", "mcquitty", "median", "centroid", "ward.D2")
linkage <- c("spearman", "pearson")

for(i in 1:length(finalDataFrameRaw_data[[1]])){

  number1 <- col(col1Raw[[i]], correlation)
  number2 <- col(col2Raw[[i]], linkage)
  number3 <- col(col3Raw[[i]], marks)

  #print(paste0(number1, " ", number2))
  temp[nrow(temp) + 1, ] <- c(number1, number2, number3)

}

finalDataFrameRaw_data <- cbind(finalDataFrameRaw_data, temp)

png(filename=paste("randomisedTestRawDifferentmarksAndLinkage.png", collapse=",", sep=""), height=4000, width=7000); par(mar=c(4, 4, 4, 4))
palette(rainbow(30))
#plot(cophenRaw, bakerRaw, col = palette(rainbow(30))[as.numeric(finalDataFrameRaw_data$correlation)], cex=2, pch = c(18, 17)[as.numeric(finalDataFrameRaw_data$linkage)],xlab="Cophenetic index", ylab="Baker's gamma index", ylim=range(-0.1:1.5))
plot(cophenRaw, bakerRaw, col = palette(rainbow(30))[as.numeric(finalDataFrameRaw_data$correlation)], cex=6, pch = c(18, 17)[as.numeric(finalDataFrameRaw_data$linkage)],xlab="Cophenetic index", ylab="Baker's gamma index", ylim=range(-0.1,0.4))
text(as.numeric(cophenRaw), as.numeric(bakerRaw), finalDataFrameRaw_data$correlation.used, cex=5, pos = 3, col="black")
grid(nx = NULL, ny = NULL)
dev.off()


png(filename=paste("randomisedTestRawMarksColouring.png", collapse=",", sep=""), height=4000, width=7000); par(mar=c(4, 4, 4, 4)+10)
#palette(rainbow(30))
#plot(cophenRaw, bakerRaw, col = palette(rainbow(30))[as.numeric(finalDataFrameRaw_data$correlation)], cex=2, pch = c(18, 17)[as.numeric(finalDataFrameRaw_data$linkage)],xlab="Cophenetic index", ylab="Baker's gamma index", ylim=range(-0.1:1.5))
plot(cophenRaw, bakerRaw, col = palette(c("red", "green", "blue"))[as.numeric(finalDataFrameRaw_data$mark)], cex=6, pch = c(15, 16, 17)[as.numeric(finalDataFrameRaw_data$mark)], ylim=range(0.0,0.4), cex.axis=5, cex.lab=6, ann = FALSE)
text(as.numeric(cophenRaw), as.numeric(bakerRaw), finalDataFrameRaw_data$mark, cex=5, pos = 3, col="black")
title(ylab="Baker's gamma index", cex.lab = 6, line = 7)
title(xlab="Cophenetic index", cex.lab = 6, line = 7)
grid(nx = NULL, ny = NULL, col="black")
dev.off()

#####################
#####log 10 data#####
#####################

cophenLog10 <- finalDataFrameLog10_data$cophenetic.value
bakerLog10 <- finalDataFrameLog10_data$baker.s.gamma.value
col1Log10 <- finalDataFrameLog10_data$correlation.used
col2Log10 <- finalDataFrameLog10_data$type.of.linkage.used
col3Log10 <- finalDataFrameLog10_data$second.object

temp <- as.data.frame(setNames(replicate(3,numeric(0), simplify = F), c("correlation", "linkage", "mark")))


marks <-c("fine grained dataset H3K27ac cluster", "fine grained dataset H3K4me3 cluster", "fine grained dataset H3K27me3 cluster")
correlation <- c("single", "complete", "average", "mcquitty", "median", "centroid", "ward.D2")
linkage <- c("spearman", "pearson")

for(i in 1:length(finalDataFrameLog10_data[[1]])){

  number1 <- col(col1Log10[[i]], correlation)
  number2 <- col(col2Log10[[i]], linkage)
  number3 <- col(col3Log10[[i]], marks)

  #print(paste0(number1, " ", number2))
  temp[nrow(temp) + 1, ] <- c(number1, number2, number3)

}

finalDataFrameLog10_data <- cbind(finalDataFrameLog10_data, temp)

png(filename=paste("randomisedTestLog10DifferentmarksAndLinkage.png", collapse=",", sep=""), height=4000, width=7000)
palette(rainbow(30))
#plot(cophenLog10, bakerLog10, col = palette(rainbow(30))[as.numeric(finalDataFrameLog10_data$correlation)], cex=2, pch = c(18, 17)[as.numeric(finalDataFrameLog10_data$linkage)],xlab="Cophenetic index", ylab="Baker's gamma index", ylim=range(-0.1:1.5))
plot(cophenLog10, bakerLog10, col = palette(rainbow(30))[as.numeric(finalDataFrameLog10_data$correlation)], cex=6, pch = c(18, 17)[as.numeric(finalDataFrameLog10_data$linkage)],xlab="Cophenetic index", ylab="Baker's gamma index", ylim=range(-0.1,0.4))
text(as.numeric(cophenLog10), as.numeric(bakerLog10), finalDataFrameLog10_data$correlation.used, cex=5, pos = 3, col="black")
grid(nx = NULL, ny = NULL)
dev.off()


png(filename=paste("randomisedTestLog10MarksColouring.png", collapse=",", sep=""), height=4000, width=7000)
#palette(rainbow(30))
#plot(cophenLog10, bakerLog10, col = palette(rainbow(30))[as.numeric(finalDataFrameLog10_data$correlation)], cex=2, pch = c(18, 17)[as.numeric(finalDataFrameLog10_data$linkage)],xlab="Cophenetic index", ylab="Baker's gamma index", ylim=range(-0.1:1.5))
plot(cophenLog10, bakerLog10, col = palette(c("red", "green", "blue"))[as.numeric(finalDataFrameLog10_data$mark)], cex=6, pch = c(15, 16, 17)[as.numeric(finalDataFrameLog10_data$mark)],xlab="Cophenetic index", ylab="Baker's gamma index", ylim=range(-0.1,0.4))
text(as.numeric(cophenLog10), as.numeric(bakerLog10), finalDataFrameLog10_data$mark, cex=5, pos = 3, col="black")
grid(nx = NULL, ny = NULL)
dev.off()

########################################################

#' if multiple tests are carried out, the maximum values
#'  can be found and stored using the code below for
#'  both raw and log10 data.

########################################################

maxRandomTest <- as.data.frame(setNames(replicate(6,numeric(0), simplify = F), c("first object", "second object", "type of linkage used", "correlation used", "cophenetic value", "baker's gamma value")))

################
####raw data####
################

maxCorrRaw <- max(finalDataFrameRaw_data$cophenetic.value, na.rm = TRUE)
maxCorrRawdataframe <- finalDataFrameRaw_data[grep(paste0(maxCorrRaw), finalDataFrameRaw_data$cophenetic.value), ]

nameVec <- as.data.frame(setNames(replicate(1,numeric(0), simplify = F), c("preprocessing")))
for(i in 1:length(maxCorrRawdataframe[[1]])){
  nameVec[nrow(nameVec) + 1, ] <- c(paste0("no log"))
}
maxCorrRawdataframe <- cbind(maxCorrRawdataframe, nameVec)
maxRandomTest<-rbind(maxRandomTest, maxCorrRawdataframe)

#####################
#####log 10 data#####
#####################

finalDataFrameLog10 <- get(load("/finalDataFrame(log10).RData"))
finalDataFrameRaw_data <- get(load("/finalDataFrame(raw_data).RData"))

maxCorrLog10 <- max(finalDataFrameLog10$cophenetic.value, na.rm = TRUE)
maxCorrLog10dataframe <- finalDataFrameLog10[grep(paste0(maxCorrLog10), finalDataFrameLog10$cophenetic.value), ]

nameVec <- as.data.frame(setNames(replicate(1,numeric(0), simplify = F), c("preprocessing")))
for(i in 1:length(maxCorrLog10dataframe[[1]])){
  nameVec[nrow(nameVec) + 1, ] <- c(paste0("Log10"))
}
maxCorrLog10dataframe <- cbind(maxCorrLog10dataframe, nameVec)
maxRandomTest <-rbind(maxRandomTest, maxCorrLog10dataframe)

################

save(maxRandomTest, file="maxRandomTest.RData")
write.csv(maxRandomTest, file = "maxRandomTest.csv")

#################
