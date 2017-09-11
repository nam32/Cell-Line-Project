# dataPrepLog2
#
# This is a function called 'dataPrepLog2'
#
#
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

dataPrepLog2 <- function(inData, marks, i){

  dat = inData$res

  remove_rows = which(apply(dat, 1, function(x) all(is.na(x))))
  if(length(remove_rows)) {
    dat_na_rm = dat[-remove_rows,]
  } else {
    dat_na_rm = dat
  }

  dat_na_rm = log2(dat_na_rm)
  dat_na_rm = dat_na_rm[,!apply(dat_na_rm, 2, function(x) sd(x)==0)] # remove regions with no variance (rows with 0 variance)
  dim(dat_na_rm)

  dat_na_rm = dat_na_rm[,!apply(is.na(dat_na_rm), 2, all)] # remove regions with no data (remove columns i.e. genes with no data)
  dim(dat_na_rm)

  dat_na_rm = t(dat_na_rm)
  assign(paste0("new_cluster_", paste0(i)), as.dendrogram(varclus(dat_na_rm, similarity="pearson")))
  assign(paste0("new_cluster_data_", paste0(i)), dat_na_rm)
  return(eval(as.name(paste0("new_cluster_data_", i))))

}
