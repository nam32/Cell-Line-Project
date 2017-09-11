# spectralClusterSummary
#
# This is a function called 'spectralClusterSummary'
#
#
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

spectralClusterSummary <- function(new_new_cluster_data_1, vec, groupNo){
  new_cluster_data_1_dist <- dist2(new_new_cluster_data_1, new_new_cluster_data_1)
  group_original1 = spectralClustering(new_cluster_data_1_dist, as.numeric(groupNo))

  row.names(new_cluster_data_1_dist) <- cellTypes
  colnames(new_cluster_data_1_dist) <- vec
  new_cluster_data_1_dist_spectral <- cbind(cellTypes, group_original1)
  new_cluster_data_1_dist_spectral<-(data.frame(new_cluster_data_1_dist_spectral))
  new_cluster_data_1_dist_spectral <- new_cluster_data_1_dist_spectral[with(new_cluster_data_1_dist_spectral, order(new_cluster_data_1_dist_spectral$group_original1)), ]

  return(new_cluster_data_1_dist_spectral)
}
