# theme_thesis
#
# This is a function called 'theme_thesis'
# taken from epiChoose package by Dr. Aiden MacNamara
#

theme_thesis <- function(base_size = 30){
  theme_bw(base_size = base_size) + theme(text = element_text(size = base_size, family = ""), axis.title = element_text(size = rel(1)), axis.text = element_text(size = rel(0.75)), axis.title.y = element_text(vjust = 0.3), axis.title.x = element_text(vjust = 0.3), plot.title = element_text(size = rel(1.33), vjust = 2), legend.title = element_blank(), legend.key = element_blank())
}

