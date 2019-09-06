###################################################################################
# Formation_patterns.R 
#
# Auxiliary routines for Bandt-Pompe distribution implementation and variations
#
# Author: Eduarda Chagas
# Date : Aug 2019
# Contact: eduarda.chagas@dcc.ufmg.br
###################################################################################

#################################### Packages #####################################

require(gtools)


################## One dimensional Bandt-Pompe Version ############################

formationPattern<-function(serie,dimension,delay,option){
  n_symbols = i = 1
  n = length(serie)
  p_patterns = elements = index2 = matrix(nrow=n,ncol=dimension)
  index = c(0:(dimension-1)) 
  while(i <= n){    
    first = i
    if((i+dimension-1)<=n){
      index2[n_symbols,] = i:(i+dimension-1)
      elements[n_symbols,] = serie[i:(i+dimension-1)]
      p_patterns[n_symbols,] = index[order(elements[n_symbols,])]
      i = first + delay
      n_symbols = n_symbols + 1
    }else break
  }
  if(option == 0){
    p_patterns = na.omit(p_patterns)
    return(p_patterns[1:dim(p_patterns)[1],])
  }else if(option == 1){
    elements = na.omit(elements)
    return(elements[1:dim(elements)[1],])    
  }else{
    index2 = na.omit(index2)
    return(index2[1:dim(index2)[1],])    
  }
}

################## Bidimensional Horizontal Patterns ############################

formationPatternsTexturesHorizontal <- function(img, dimension, delay, xy){
  initial.element.row = initial.element.col = number.patterns = 1
  matriz <- array(0, dim = c(dimension, dimension, xy))
  while(number.patterns <= xy){
    if(((initial.element.col + dimension - 1) <= dim(img)[1]) && ((initial.element.row + dimension - 1) <= dim(img)[2])){
      for(j in 1:dimension){
        matriz[, j,number.patterns] =  img[initial.element.col:(initial.element.col + dimension - 1), initial.element.row + j - 1]
      }
      number.patterns = number.patterns + 1
      initial.element.col = initial.element.col + 1
    }
    else if(((initial.element.col + dimension - 1) > dim(img)[1]) && ((initial.element.row + delay + dimension - 1) <= dim(img)[2])){
      initial.element.col = 1
      initial.element.row = initial.element.row + delay 
    }else{
      break
    }
  }
  matriz[,,1:(number.patterns-1)]
}

################## Bidimensional Vertical Patterns ############################

formationPatternsTexturesVertical <- function(img, dimension, delay, xy){
  initial.element.row = initial.element.col = number.patterns = 1
  matriz <- array(0, dim = c(dimension, dimension, xy))
  while(number.patterns <= xy){
    if(((initial.element.col + dimension - 1) <= dim(img)[1]) && ((initial.element.row + dimension - 1) <= dim(img)[2])){
      for(j in 1:dimension){
        matriz[j, ,number.patterns] =  img[initial.element.col:(initial.element.col + dimension - 1), initial.element.row + j - 1]
      }
      number.patterns = number.patterns + 1
      initial.element.col = initial.element.col + 1
    }
    else if(((initial.element.col + dimension - 1) > dim(img)[1]) && ((initial.element.row + delay + dimension - 1) <= dim(img)[2])){
      initial.element.col = 1
      initial.element.row = initial.element.row + delay 
    }else{
      break
    }
  }
  matriz[,,1:(number.patterns-1)]
}

################## Bidimensional Hilbert Patterns ############################

formationPatternsTexturesHilbert <- function(img, dimension, delay){
  hilbertcurve <- unlist(read.table("../../Data/Hilbert/HilbertCurves128.txt")) + 1
  matriz <- array(0, dim = c(dim(img)[1]*dim(img)[1], dimension))
  nsymbols = init = 1
  vector.img.hilbert = img[hilbertcurve]
  while(init <= (dim(img)[1]) + dimension){
    matriz[nsymbols,] = vector.img.hilbert[init:(init+dimension-1)]
    init = init + delay
    nsymbols = nsymbols + 1
  }
  matriz = matriz[1:nsymbols-1,]
  matriz
}
