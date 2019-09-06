###################################################################################
# WPE_RMS.R 
#
# Implementation of Weighted-permutation entropy with RMS modification function 
#
# Author: Eduarda Chagas
# Date : Aug 2019
# Contact: eduarda.chagas@dcc.ufmg.br
####################################################################################

require(gtools)

weights.RMS<-function(series,dimension,delay){
  groups = formationPattern(series,dimension,delay,1)
  weight = w = rep(0, dim(groups)[1])
  for(i in 1:dim(groups)[1]){
    weight[i] = sqrt(mean(groups[i,]^2))
  }
  weight
}

WPE.RMS<-function(series,dimension,delay){
  
  w = weights.RMS(series,dimension,delay)
  simbols = formationPattern(series,dimension,delay,0)
  patterns = definePatterns(dimension)
  sw = rep(0,factorial(dimension))
  
  for(i in 1:factorial(dimension)){
    for(j in 1:dim(simbols)[1]){
      if(all(simbols[j,] != patterns[i,]))
        sw[i] = sw[i] + w[j]
    }
  }
  
  pw = sw/sum(sw)
  pw
}

