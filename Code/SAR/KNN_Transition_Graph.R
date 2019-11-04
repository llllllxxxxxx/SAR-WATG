###################################################################################################################
# Transition_graph.R
#
# Generate Plot of Analysis of SAR images using Hilbert space-filling curves and ordinal patterns transition graphs
# (sliding window)
#
# Author: Eduarda Chagas
# Date : Oct 2019
# Contact: eduarda-chagas@ufmg.br
####################################################################################################################

############################################# Packages #################################################

require(gtools)
require(ggplot2)
require(ggthemes)
require(igraph)
require(ggpubr)
require(class)
source("SAR_TimeSerie.R")
source("../Theory_Information.R")
#source("Theory_Information.R")

############################### Transition graphs functions ############################################

sliding.window <- function(series, dimension, delay){
  i = j = 1
  n = length(series)
  elements = matrix(nrow = n-(dimension-1)*delay, ncol = dimension)
  s = seq(0, (0 + (dimension - 1)*delay), by = delay)
  while(i+((dimension-1)*delay) <= length(series)){
    elements[j,] = series[s + i]
    i = i + 1
    j = j + 1
  }
  elements
}

formation.pattern <- function(elements){
  patterns = matrix(nrow = dim(elements)[1], ncol = dim(elements)[2])
  index = c(1:(dim(elements)[2])) 
  for(i in 1:dim(elements)[1]){
    patterns[i,] = index[order(elements[i,])]
  }
  patterns
}


define.symbols<-function(dimension){
  d = c(1:dimension)
  symbol = matrix(unlist(permutations(n=dimension,r=dimension,v=d)),nrow = factorial(dimension),ncol = dimension,byrow = FALSE)
  symbol
}

pattern.wedding <- function(patterns){
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

transition.graph <- function(elements, wedding, D){
  m = length(wedding)
  graph = matrix(0,nrow = factorial(D), ncol = factorial(D))
  weight.total = 0
  
  for(i in 1:(m-1)){
    weight.i1 = (max(elements[i,]) - min(elements[i,]))
    weight.i2 = (max(elements[i+1,]) - min(elements[i+1,]))
    graph[wedding[i],wedding[i+1]] = graph[wedding[i],wedding[i+1]] + abs(weight.i1 - weight.i2)
    weight.total = weight.total + abs(weight.i1 - weight.i2)
  }
  graph = graph/weight.total
  graph
}

###################################### Function of Plot ##############################################

HC.color.shape.signal <- function(color.signal, shape.signal, signal.values){
  
  shape.select <- c(17,18,19,8)
  
  # Paleta montada a partir de https://coolors.co/
  rainbow.colors <- palette(c("#494947", #DarkGreen
                              "#7494EA", #MutedDarkBlue
                              "#B14AED", #Violet
                              "#44CCFF", #BrightLightBlue
                              "#35FF69", #BrightGreen
                              "#ED8438", #Orange
                              "#E7AD99", #Pink
                              "#C18C5D", #LightBrown
                              "#BF6F00", #DarkYellow
                              "#FB4D3D", #BrightRed
                              "#495867", #DarkGray
                              "#98CE00", #SHEEN GREEN
                              "#98838F")) #MOUNTBATTEN PINK
  
  Color = rainbow.colors[color.signal]
  Shape = shape.select[shape.signal]
  Texture = color.signal
  signal.values <- data.frame("H" = signal.values[,1], "C" = signal.values[,2], "Color" = Color, "Shape" = Shape, "Texture" = Texture)
  
  p = ggplot(signal.values, aes(x = signal.values$H, y = signal.values$C)) +
    geom_point(shape = signal.values$Shape, color = signal.values$Color, size = 1) + 
    labs(x="", y="") + 
    theme_few(base_size = 18, base_family = "serif")  + theme(plot.title = element_text(hjust=0.5)) + 
    scale_colour_few("Dark")
  return(p)
} 


###################################### Image Sample Parameters #######################################

#Defining the number of samples from each region analyzed
ns.guatemala = 40
dimen.guatemala <- matrix(nrow = ns.guatemala, ncol = 4)
ns.canaveral.behavior1 = 40
dimen.canaveral.behavior1 <- matrix(nrow = ns.canaveral.behavior1, ncol = 4)
ns.canaveral.behavior2 = 40
dimen.canaveral.behavior2 <- matrix(nrow = ns.canaveral.behavior2, ncol = 4)
ns.munich = 40
dimen.munich <- matrix(nrow = ns.munich, ncol = 4)
n.total = (ns.guatemala + ns.canaveral.behavior1 + ns.canaveral.behavior2 + ns.munich)

#The SAR data is available on https://drive.google.com/file/d/1jtbOcYwQfysfcUp4UhoA7lSl4_tPIqfa/view?usp=sharing and
# correspond to HHHH band of an image taken from the Cape Canaveral (acquired Sep 22, 2016)

#Ocean regions in Cape Canaveral
row1 <- c(50, 100, 150, 200, 250, 350, 450, 550, 650, 750)
row2 <- c(50, 100, 150, 200, 250, 300, 350, 400, 450, 550)
row3 <- c(50, 150, 250, 350, 450, 550, 650, 750, 800, 850)
row4 <- c(250, 350, 450, 550, 650, 750, 850, 950, 1050)
row5 <- c(50, 150, 250, 350, 450, 550, 650, 750, 850, 950, 1050)
row6 <- c(50, 150, 250, 350, 450, 550, 650, 750, 800, 850, 950)
row7 <- c(50, 150, 250, 350, 450, 550, 650, 750, 850, 950)
cols <- c(1700, 1850, 1550, 1400, 1, 200, 350, 550)
#{Behavior 1}
dimen.canaveral.behavior1[1:10,] <- c(row1, rep(128, 10), rep(cols[1], 10), rep(128, 10))
dimen.canaveral.behavior1[11:20,] <- c(row2, rep(128, 10), rep(cols[2], 10), rep(128, 10))
dimen.canaveral.behavior1[21:30,] <- c(row3, rep(128, 10), rep(cols[3], 10), rep(128, 10))
dimen.canaveral.behavior1[31:40,] <- c(row3, rep(128, 10), rep(cols[4], 10), rep(128, 10))
#{Behavior 2}
dimen.canaveral.behavior2[1:9,] <- c(row4, rep(128, 9), rep(cols[5], 9), rep(128, 9))
dimen.canaveral.behavior2[10:20,] <- c(row5, rep(128, 11), rep(cols[6], 11), rep(128, 11))
dimen.canaveral.behavior2[21:30,] <- c(row7, rep(128, 10), rep(cols[7], 10), rep(128, 10))
dimen.canaveral.behavior2[31:40,] <- c(row7, rep(128, 10), rep(cols[8], 10), rep(128, 10))

#The SAR data is available on https://drive.google.com/file/d/1pO6p_UI9Cgdci9y6jVynAv8SrrAvv7K8/view?usp=sharing and
# correspond to HHHH band of an image taken from the Munich, Germany (acquired Jun 5, 2015) 

#Urban regions in Munich
row1 <- seq(3000, 3950, by = 50)
row2 <- rep(c(4300, 4350), 5)
row3 <- c(rep(seq(2300, 2450, by = 50), 2), 2400, 2450)
cols <- 400
cols2 <- c(1300, 1300, 1350, 1350, 1400, 1400, 1450, 1450, 1500, 1500)
cols3 <- c(rep(500, 4), rep(400, 4), rep(300, 2))
dimen.munich[1:20,] <- c(row1, rep(128, 20), rep(cols, 20), rep(128, 20))
dimen.munich[21:30,] <- c(row2, rep(128, 10), cols2, rep(128, 10))
dimen.munich[31:40,] <- c(row3, rep(128, 10), cols3, rep(128, 10))

#Forest regions in Guatemala
row1 <- seq(5150, 6100, by = 50)
row2 <- seq(5200, 5650, by = 50)
row3 <- seq(4100, 4200, by = 50)
row4 <- seq(1000, 1150, by = 50)
cols <- c(2700, 2800, 2930, 1930, 1870)
dimen.guatemala[1:20,] <- c(row1, rep(128, 20), rep(cols[1], 20), rep(128, 20))
dimen.guatemala[21:30,] <- c(row2, rep(128, 10), rep(cols[2], 10), rep(128, 10))
dimen.guatemala[31:33,] <- c(row3, rep(128, 3), rep(cols[3], 3), rep(128, 3))
dimen.guatemala[34:37,] <- c(row4, rep(128, 4), rep(cols[4], 4), rep(128, 4))
dimen.guatemala[38:40,] <- c(row4[1:3], rep(128, 3), rep(cols[5], 3), rep(128, 3))

###################################### Function of Analysis ##########################################

plot.transition.graph.knn <- function(){
  
  n = 3 #Dimension parameter
  tal = 1 #Delay parameter
  hilbertcurve <- unlist(read.table("../../Data/Hilbert/HilbertCurves128.txt")) + 1
  types <- c(rep(1,40), rep(2,40), rep(3,40), rep(4, 40))
  regions <- c(rep(1,40), rep(2,80), rep(4, 40))
  
  
  Entropy.Complexity.csv <-read.csv(file="EntropyComplexity.csv", header=TRUE, sep=",")
  
  Entropy.Complexity <- matrix(nrow = n.total, ncol = 3)
  
  Entropy.Complexity[,1] = Entropy.Complexity.csv[1:n.total,2]
  Entropy.Complexity[,2] = Entropy.Complexity.csv[1:n.total,3]
  
  #Diferent classes
  Entropy.Complexity[1:80, 3] = 1      #Guatemala
  Entropy.Complexity[81:120, 3] = 2     #Canaveral
  Entropy.Complexity[121:160, 3] = 3    #Munich
  
  ##Generate a random number that is 95% of the total number of rows in dataset
  ran <- sample(1:nrow(Entropy.Complexity), 0.95*nrow(Entropy.Complexity))
  HC.train <- Entropy.Complexity[ran,1:2] 
  HC.test <- Entropy.Complexity[-ran,1:2]  
  HC.target.category <- Entropy.Complexity[ran,3]
  HC.test.category <- Entropy.Complexity[-ran,3]
  
  ##run knn function
  pr <- knn(HC.train, HC.test, cl = HC.target.category, k = 3)
  
  ##create confusion matrix
  tab <- table(pr, HC.test.category)
  
  ##this function divides the correct predictions by total number of predictions that tell us how accurate teh model is.
  accuracy <- function(x){sum(diag(x)/(sum(rowSums(x)))) * 100}
  accuracy(tab)
}

plot.transition.graph.knn()