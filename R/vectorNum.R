# vectorNum
#
# This is a function called 'vectorNum'
# which it
#
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

vectorNum <- function(target, vector, specificregion){

  vec1<-eval(as.name(paste0(vector)))
  spe <- vec1[grep(paste0(specificregion),vec1$corr), ]

  numberVec <- as.data.frame(setNames(replicate(1,numeric(0), simplify = F), c(target)))
  for(i in 1:length(spe$k)){
    number <- as.numeric(spe$k[[i]])
    if(number == as.numeric(target)) {numberVec[nrow(numberVec) + 1, ] <- c(paste0("y"))} else {numberVec[nrow(numberVec) + 1, ] <- c(paste0("n"))}
  }

  return(numberVec)
}

