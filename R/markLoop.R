# marksLoop
#
# This is a function called 'marksLoop'
#
#
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

marksLoop <- function(marksVec, correlationsVec, linkagesVec){

  finalDataframe <- as.data.frame(setNames(replicate(6,numeric(0), simplify = F), c("first object", "second object", "type of linkage used", "correlation used", "cophenetic value", "baker's gamma value")))

  for (i in 1:length(marksVec)){
    for (j in 1:length(correlationsVec)){
      for (k in 1:length(linkagesVec)){

        currentCorr <- correlationsVec[[j]]
        currentMark <- marksVec[[i]]
        currentLinkage <- linkagesVec[[k]]
        clusterName <- paste0("fine grained dataset ",currentMark," cluster")
        dataMyTreeVariation <- myTreeVariation("fine grained synthetic cluster",paste0(clusterName), paste0(currentCorr), paste0(currentLinkage), newdataMerge, newdataHeight, newdataOrder, newdataLabels,  eval(as.name(paste0("new_cluster_data_", i))), 20)
        finalDataframe <- rbind(finalDataframe, dataMyTreeVariation)

      }
    }
  }

  return(finalDataframe)
}

