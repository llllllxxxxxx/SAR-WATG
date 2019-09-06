###################################################################################
# TheoryInformation.R 
#
# Implementation of theory information descriptors 
#
# Author: Eduarda Chagas
# Date : Set 2019
# Contact: eduarda.chagas@dcc.ufmg.br
####################################################################################

require(ggplot2)
require(ggthemes)

Entropy.Complexity.Points <- function(probabilities, ns, dimension){
  C <- H <-rep(0,ns)
  fact <- factorial(dimension)
  for(i in c(1:ns)){
    H[i] = shannonNormalized(probabilities[i,1:fact])
    C[i] = Ccomplexity(probabilities[i,1:fact])
  }
  Entropy.Complexity <- data.frame(H, C)
  return(Entropy.Complexity)
}

Entropy.Complexity.Fisher.Points <- function(probabilities, ns, dimension){
  C <- H <- Fs<-rep(0,ns)
  fact <- factorial(dimension)
  for(i in c(1:ns)){
    H[i] = shannonNormalized(probabilities[i,1:fact])
    C[i] = Ccomplexity(probabilities[i,1:fact])
    Fs[i] = fis(probabilities[i,1:fact])
  }
  Entropy.Complexity.Fisher <- data.frame(H, Fs, C)
  return(Entropy.Complexity.Fisher)
}

shannonEntropy <- function(p){
  h <- p * log(p)
  h[is.nan(h)] <- 0
  return(-sum(h))
}

shannonNormalized <- function(p){
  return(shannonEntropy(p)/log(length(p)))
}

jensenDivergence<-function(p){
  cc = rep(1/length(p),length(p))
  s_p = shannonEntropy(p)
  s_q = shannonEntropy(cc)
  s_pq = shannonEntropy((p+cc)/2)
  divergence = sum(s_pq - (s_p/2) - (s_q/2))
  return(divergence)
}

constant <- function(p){
  k = (0.5)/length(p)
  a1 = (0.5 + k) * log(0.5 + k)
  a2 = (length(p) - 1) * k * log(k)
  a3 = (1 - 0.5) * log(length(p))
  b = -1/(a1 + a2 + a3)
  return(b)
}

Ccomplexity<-function(p){
  cc <- jensenDivergence(p) * constant(p) * shannonNormalized(p)
  return(cc)
}

cotas <- function(dimension){
  c1x = readingMPR(dimension,1)
  c1y = readingMPR(dimension,2)
  c2x = readingMPR(dimension,3)
  c2y = readingMPR(dimension,4)
  
  p = qplot(xlab=expression(H), ylab=expression(C)) +
    theme(plot.title = element_text(hjust=0.5)) +
    geom_line(aes(x=c2x, y=c2y), size=0.5, color="gray") +
    geom_line(aes(x=c1x, c1y), size=0.5, color="gray")
  return(p)
}

HCPlane <- function(Entropy.Complexity){
  p = ggplot(Entropy.Complexity, aes(x = H, y = C, label = Texture))+theme(plot.title = element_text(hjust=0.5))+geom_point(size=2)+geom_text(check_overlap = TRUE, hjust = 0, nudge_x = 0.0015, size = 2) 
  return(p)
} 

readingMPR<-function(dimension,option=0){
  if(dimension == 3){ 
    continua = "../../Data/trozos/continuaN6.txt"
    trozo = "../../Data/trozos/trozosN6.txt"
  }
  if(dimension == 4){ 
    continua = "../../Data/trozos/continuaN24.txt"
    trozo = "../../Data/trozos/trozosN24.txt"
  }
  if(dimension == 5){ 
    continua = "../../Data/trozos/continuaN120.txt"
    trozo = "../../Data/trozos/trozosN120.txt"
  }
  if(dimension == 6){ 
    continua = "../../Data/trozos/continuaN720.txt"
    trozo = "../../Data/trozos/trozosN720.txt"
  }
  curva1x = read.table(continua, stringsAsFactors=FALSE, fileEncoding="latin1")[,1]
  if(option==1) return(curva1x)
  curva1y = read.table(continua, stringsAsFactors=FALSE, fileEncoding="latin1")[,2]
  if(option==2) return(curva1y)
  curva2x = read.table(trozo, stringsAsFactors=FALSE, fileEncoding="latin1")[,1]
  if(option==3) return(curva2x)
  curva2y = read.table(trozo, stringsAsFactors=FALSE, fileEncoding="latin1")[,2]
  if(option==4) return(curva2y)
}

HCPlane.Classes <- function(probability, dimension, Entropy.Complexity){
  Color <- Shape <-rep(0,dim(probability)[1])
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
  
  Color = rainbow.colors[probability[,factorial(dimension) + 1]]
  Shape = shape.select[probability[,factorial(dimension) + 2]]
  Texture = probability[,factorial(dimension) + 1]
  Entropy.Complexity <- data.frame("H" = Entropy.Complexity$H, "C" = Entropy.Complexity$C, "Color" = Color, "Shape" = Shape, "Texture" = Texture)
  p = cotas(dimension)
  p = p + 
    geom_point(aes(x = Entropy.Complexity$H, y = Entropy.Complexity$C), shape = Entropy.Complexity$Shape, color = Entropy.Complexity$Color, size = 1) +
    theme(plot.title = element_text(hjust=0.5)) 
  #geom_text(check_overlap = TRUE, hjust = 0, nudge_x = 0.0015, size = 2) 
  return(p)
} 
