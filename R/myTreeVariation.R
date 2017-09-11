# myTreeVariation
#
# This is a function called 'myTreeVariation'
#
#
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

myTreeVariation <- function(name_1, name_2, linkage, correlation, newdataMerge, newdataHeight, newdataOrder, newdataLabels, clusterToCompareTo, repeat1){

  resultsDataframe <- as.data.frame(setNames(replicate(6,numeric(0), simplify = F), c("first object", "second object", "type of linkage used", "correlation used", "cophenetic value", "baker's gamma value")))

  for (j in 1:paste0(repeat1))
  {
    factorTime <- runif(paste0(length(newdataHeight)), 0.8, 1.5)
    newdata <- list()
    newdata$merge <- newdataMerge
    newdata$height <- (newdataHeight*factorTime)
    newdata$order <- newdataOrder
    newdata$labels <- newdataLabels
    class(newdata) <- "hclust"

    cluster_synthetic <- as.dendrogram(newdata)

    compareToDend <- as.dendrogram(varclus(clusterToCompareTo, similarity=paste0(linkage), method=paste0(correlation)))

    resultsDataframe[nrow(resultsDataframe) + 1, ] <- c(name_1,name_2, linkage, correlation, cor_cophenetic(cluster_synthetic, compareToDend), cor_bakers_gamma(cluster_synthetic, compareToDend))
  }

  return(resultsDataframe)

}
