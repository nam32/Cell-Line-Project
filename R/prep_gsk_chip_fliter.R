# prep_gsk_chip_filter (from Dr. Aiden epiChoose package)
#
# This is a function called 'prep_gsk_chip_filter'
#
#
#
# Some useful keyboard shortcuts for package authoring:
#

prep_gsk_chip_filter <- function(gsk_chip) {

  # what are the common labels unlist : list to data frame
  all_labels = sort(unique(unlist(lapply(gsk_chip, function(x) x$annot$Label))))
  all_ix = lapply(gsk_chip, function(x) match(all_labels, x$annot$Label))

  gsk_chip_filtered = gsk_chip

  for(i in 1:length(gsk_chip)) {
    gsk_chip_filtered[[i]]$res = gsk_chip[[i]]$res[all_ix[[i]],]
    gsk_chip_filtered[[i]]$annot = gsk_chip[[i]]$annot[all_ix[[i]],]
  }

  return(gsk_chip_filtered)

}


