# colReturn
#
# This is a function called 'colContains'
#
#
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

colContains <- function(word, vec){
  for(j in 1:length(vec)){

    if(grepl(vec[[j]] , word)){

      return(j)

    }

  }
}
