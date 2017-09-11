# colVar
#
# This is a function called 'colVar'
#
#


colVar <- function(x) {
  colSums((x - colMeans(x))^2)/(dim(x)[2] - 1)
}
