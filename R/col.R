# col
#
# This is a function called 'col'
#
#
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

col <- function(word, vec){
  for(j in 1:length(vec)){

    if(word == vec[[j]]){

      return(j)
    }

  }

}
