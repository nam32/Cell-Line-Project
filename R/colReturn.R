# colReturn
#
# This is a function called 'colReturn'
#
#
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

colReturn <- function(word, vec){
  for(j in 1:length(vec)){

    if(grepl(vec[[j]] , word)){

      return(vec[[j]])

    }

  }
}
