########################################################################################################
# WATG.R
#
# Analysis of SAR images using Hilbert space-filling curves and WATG
#
# Author: Eduarda Chagas
# Date : Apr 2020
# Contact: eduarda.chagas@dcc.ufmg.br
########################################################################################################

############################################# Packages #################################################

if(!require(gtools)) install.packages("gtools")
if(!require(raster)) install.packages("raster")
source("theory_information.R")

############################### Transition graphs functions ############################################

FP <- function(n, dimension, delay){
  dyn.load("FormationPatterns.so")
  p <- .Call("FormationPatterns", n, dimension, delay)
  p = t(p) + 1
  return(p)
}

formationPattern <- function(serie, dimension, delay, option){
  i = 1
  n = length(serie)
  p_patterns = elements = index2 = matrix(nrow=n,ncol=dimension)
  index = c(0:(dimension-1))
  
  index2 = FP(length(serie), dimension, delay)
  
  while((i + ((dimension-1)*delay)) <= n){ 
    elements[i,] = serie[index2[i,]]
    p_patterns[i,] = index[order(elements[i,])]
    i = i + 1
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


define.symbols <- function(dimension){
  d = c(1:dimension)
  symbol = matrix(unlist(permutations(n=dimension,r=dimension,v=d)),nrow = factorial(dimension),ncol = dimension,byrow = FALSE)
  symbol
}

pattern.wedding <- function(patterns){
  patterns = patterns + 1
  m = dim(patterns)[1]
  D = dim(patterns)[2]
  symbols = define.symbols(D)
  wedding = rep(0, m)
  for(i in 1:m){
    e = 0
    j = 1
    stop = F
    while(j <= factorial(D) && stop == F){
      if(sum(symbols[j,] == patterns[i,]) == D){
        wedding[i] = j
        stop = T
      }
      j = j + 1
    }
  }
  wedding
}

transition.graph.weight <- function(series, dimension, delay){
  
  graph = matrix(0, nrow = factorial(dimension), ncol = factorial(dimension))
  patterns = formationPattern(series, dimension, delay, 0)
  elements = formationPattern(series, dimension, delay, 1)
  wedding = pattern.wedding(patterns)
  m = length(wedding)
  weight.total = 0
  
  for(i in 1:(m-1)){
    weight.i1 = (max(elements[i,]) - min(elements[i,]))
    weight.i2 = (max(elements[i+1,]) - min(elements[i+1,]))
    graph[wedding[i],wedding[i+1]] = graph[wedding[i],wedding[i+1]] + abs(weight.i1 - weight.i2)
    weight.total = weight.total + abs(weight.i1 - weight.i2)
  }
  graph = graph/weight.total
  return(graph)
}


###################################### Image Sample Parameters #######################################

#Defining the number of samples from each region analyzed
ns.guatemala = 40
dimen.guatemala = matrix(nrow = ns.guatemala, ncol = 4)
ns.canaveral.behavior1 = 40
dimen.canaveral.behavior1 = matrix(nrow = ns.canaveral.behavior1, ncol = 4)
ns.canaveral.behavior2 = 40
dimen.canaveral.behavior2 = matrix(nrow = ns.canaveral.behavior2, ncol = 4)
ns.munich = 40
dimen.munich = matrix(nrow = ns.munich, ncol = 4)
n.total = (ns.guatemala + ns.canaveral.behavior1 + ns.canaveral.behavior2 + ns.munich)

#The SAR data is available on https://drive.google.com/file/d/1jtbOcYwQfysfcUp4UhoA7lSl4_tPIqfa/view?usp=sharing and
# correspond to HHHH band of an image taken from the Cape Canaveral (acquired Sep 22, 2016)

#Ocean regions in Cape Canaveral
row1 = c(50, 100, 150, 200, 250, 350, 450, 550, 650, 750)
row2 = c(50, 100, 150, 200, 250, 300, 350, 400, 450, 550)
row3 = c(50, 150, 250, 350, 450, 550, 650, 750, 800, 850)
row4 = c(250, 350, 450, 550, 650, 750, 850, 950, 1050)
row5 = c(50, 150, 250, 350, 450, 550, 650, 750, 850, 950, 1050)
row6 = c(50, 150, 250, 350, 450, 550, 650, 750, 800, 850, 950)
row7 = c(50, 150, 250, 350, 450, 550, 650, 750, 850, 950)
cols = c(1700, 1850, 1550, 1400, 1, 200, 350, 550)
#{Behavior 1}
dimen.canaveral.behavior1[1:9,] = c(row1[1:9], rep(128, 9), rep(cols[1], 9), rep(128, 9))
dimen.canaveral.behavior1[10,] = c(row1[10], 128, 1300, 128)
dimen.canaveral.behavior1[11:19,] = c(row2[1:9], rep(128, 9), rep(cols[2], 9), rep(128, 9))
dimen.canaveral.behavior1[20,] = c(row2[10], 128, 1350, 128)
dimen.canaveral.behavior1[21:30,] = c(row3, rep(128, 10), rep(cols[3], 10), rep(128, 10))
dimen.canaveral.behavior1[31:40,] = c(row3, rep(128, 10), rep(cols[4], 10), rep(128, 10))

#{Behavior 2}
dimen.canaveral.behavior2[1:9,] = c(row4, rep(128, 9), rep(cols[5], 9), rep(128, 9))
dimen.canaveral.behavior2[10:20,] = c(row5, rep(128, 11), rep(cols[6], 11), rep(128, 11))
dimen.canaveral.behavior2[21:30,] = c(row7, rep(128, 10), rep(cols[7], 10), rep(128, 10))
dimen.canaveral.behavior2[31:40,] = c(row7, rep(128, 10), rep(cols[8], 10), rep(128, 10))

#The SAR data is available on https://drive.google.com/file/d/1pO6p_UI9Cgdci9y6jVynAv8SrrAvv7K8/view?usp=sharing and
# correspond to HHHH band of an image taken from the Munich, Germany (acquired Jun 5, 2015) 

#Urban regions in Munich
row1 = seq(3000, 3950, by = 50)
row2 = rep(c(4300, 4350), 5)
row3 = c(rep(seq(2300, 2450, by = 50), 2), 2400, 2450)
cols = 400
cols2 = c(1300, 1300, 1350, 1350, 1400, 1400, 1450, 1450, 1500, 1500)
cols3 = c(rep(500, 4), rep(400, 4), rep(300, 2))
dimen.munich[1:20,] = c(row1, rep(128, 20), rep(cols, 20), rep(128, 20))
dimen.munich[21:30,] = c(row2, rep(128, 10), cols2, rep(128, 10))
dimen.munich[31:40,] = c(row3, rep(128, 10), cols3, rep(128, 10))

#Forest regions in Guatemala
row1 = seq(5150, 6100, by = 50)
row2 = seq(5200, 5650, by = 50)
row3 = seq(4100, 4200, by = 50)
row4 = seq(1000, 1150, by = 50)
cols = c(2700, 2800, 2930, 1930, 1870)
dimen.guatemala[1:20,] = c(row1, rep(128, 20), rep(cols[1], 20), rep(128, 20))
dimen.guatemala[21:30,] = c(row2, rep(128, 10), rep(cols[2], 10), rep(128, 10))
dimen.guatemala[31:33,] = c(row3, rep(128, 3), rep(cols[3], 3), rep(128, 3))
dimen.guatemala[34:37,] = c(row4, rep(128, 4), rep(cols[4], 4), rep(128, 4))
dimen.guatemala[38:40,] = c(row4[1:3], rep(128, 3), rep(cols[5], 3), rep(128, 3))

###################################### Function of Analysis ##########################################

transition.graph.weight.analysis <- function(){
  
  a = b = 0
  n = c(3,4,5,6) #Dimension parameter
  tal = c(1,2,3,4,5) #Delay parameter
  hilbertcurve = unlist(read.table("../Data/Hilbert/HilbertCurves128.txt")) + 1
  Entropy.Complexity.csv = matrix(nrow = n.total*length(n)*length(tal), ncol = 2)
  
  for(i in 1:(length(n)*length(tal))){
    cat("- Plane: ", i, "de 20 ", "\n")
    if(i%%5 == 1){
      a = a + 1
      b = 0
    }
    b = b + 1
    
    Entropy.Complexity <- matrix(nrow = n.total, ncol = 2)
    
    #Guatemala
    sar_data = raster(paste("../Data/", "guatemala", "/HHHH", ".grd", sep = ""))
    for(j in c(1:ns.guatemala)){
      img = getValuesBlock(sar_data, row = dimen.guatemala[j,1], nrows = dimen.guatemala[j,2], col = dimen.guatemala[j,3], ncols = dimen.guatemala[j,4], format = "matrix")
      ts = img[hilbertcurve]/max(img[hilbertcurve])
      g = transition.graph.weight(ts, n[a], tal[b])
      Entropy.Complexity[j, 1] <- shannonNormalized(as.vector(g))
      Entropy.Complexity[j, 2] <- Ccomplexity(as.vector(g))
      cat("Guatemala ", j, "\n")
    }
    #Cape Canaveral - behavior 1
    sar_data = raster(paste("../Data/", "cape", "/HHHH", ".grd", sep = ""))
    for(j in c(1:ns.canaveral.behavior1)){
      img = getValuesBlock(sar_data, row = dimen.canaveral.behavior1[j,1], nrows = dimen.canaveral.behavior1[j,2], col = dimen.canaveral.behavior1[j,3], ncols = dimen.canaveral.behavior1[j,4], format = "matrix")
      ts = img[hilbertcurve]/max(img[hilbertcurve])
      g = transition.graph.weight(ts, n[a], tal[b])
      Entropy.Complexity[ns.guatemala + j, 1] <- shannonNormalized(as.vector(g))
      Entropy.Complexity[ns.guatemala + j, 2] <- Ccomplexity(as.vector(g))
      cat("Cape 1 ", j, "\n")
    }
    #Cape Canaveral - behavior 2
    sar_data = raster(paste("../Data/", "cape", "/HHHH", ".grd", sep = ""))
    for(j in c(1:ns.canaveral.behavior2)){
      img = getValuesBlock(sar_data, row = dimen.canaveral.behavior2[j,1], nrows = dimen.canaveral.behavior2[j,2], col = dimen.canaveral.behavior2[j,3], ncols = dimen.canaveral.behavior2[j,4], format = "matrix")
      ts = img[hilbertcurve]/max(img[hilbertcurve])
      g = transition.graph.weight(ts, n[a], tal[b])
      Entropy.Complexity[(ns.canaveral.behavior1 + ns.guatemala) + j, 1] <- shannonNormalized(as.vector(g))
      Entropy.Complexity[(ns.canaveral.behavior1 + ns.guatemala) + j, 2] <- Ccomplexity(as.vector(g))
      cat("Cape 2 ", j, "\n")
    }
    #Munich
    sar_data = raster(paste("../Data/", "munich", "/HHHH", ".grd", sep = ""))
    for(j in c(1:ns.munich)){
      img = getValuesBlock(sar_data, row = dimen.munich[j,1], nrows = dimen.munich[j,2], col = dimen.munich[j,3], ncols = dimen.munich[j,4], format = "matrix")
      ts = img[hilbertcurve]/max(img[hilbertcurve])
      g = transition.graph.weight(ts, n[a], tal[b])
      Entropy.Complexity[(ns.canaveral.behavior1 + ns.canaveral.behavior2 + ns.guatemala) + j, 1] <- shannonNormalized(as.vector(g))
      Entropy.Complexity[(ns.canaveral.behavior1 + ns.canaveral.behavior2 + ns.guatemala) + j, 2] <- Ccomplexity(as.vector(g))
      cat("Munich ", j, "\n")
    }
    
    write.csv(Entropy.Complexity, paste('../Data/EntropyComplexityWATG', i, '.csv', sep = ""), row.names = FALSE)
    
    Entropy.Complexity.csv[((n.total*(i-1))+1):(n.total*i), 1] = Entropy.Complexity[,1]
    Entropy.Complexity.csv[((n.total*(i-1))+1):(n.total*i), 2] = Entropy.Complexity[,2]
  }
  
  write.csv(Entropy.Complexity.csv, '../Data/EntropyComplexityWATG.csv', row.names = FALSE)
  
}

transition.graph.weight.analysis()
