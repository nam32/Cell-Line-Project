# myTreePrep
#
# This is a function called 'myTreePrep'
#
#
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

myTreePrep <- function(vector) {
  finalVec <- c()
  test <- split(vector, ceiling(seq_along(vector)/2))
  for (i in 1:length(test)){
    height <- tail(test[[i]], n=1)
    number <- head(test[[i]], n=1)
    Vec <- rnorm(number, height, 0.0)
    finalVec <- append(finalVec, Vec)
  }

  return(finalVec)

}
