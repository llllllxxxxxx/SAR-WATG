require(ggplot2)
require(ggthemes)
require(statcomp)

Bandt.Pompe <- function(elements, dimension, elementsize){
  dyn.load("../BandtPompe.so")
  probability <- .Call("BandtPompe", elements, dimension, elementsize)
  return (probability)
}


Theory.Information.Descriptors <- function(probabilities, ns, dimension){
  Complexity = Fisher = Shannon = Euclidian =rep(0,ns)
  fact = factorial(dimension)
  
  for(i in c(1:ns)){
    Fisher[i] = fis(probabilities[i,1:fact])
    Complexity[i] = Ccomplexity(probabilities[i,1:fact])
    Shannon[i] = shannonNormalized(probabilities[i,1:fact])
    Euclidian[i] = euclidianDistance(probabilities[i,1:fact])
  }
  
  theory.descriptors =  data.frame(Fisher, Complexity, Shannon, Euclidian)
  return(theory.descriptors)
}

HCPlane.Groups <- function(Entropy.Complexity, xpoints){
  Color <-rep(0,dim(Entropy.Complexity)[1])
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
  
  Color = rainbow.colors[Entropy.Complexity$Groups]
  Texture = c(1:dim(Entropy.Complexity)[1])
  Entropy.Complexity <- data.frame("H" = Entropy.Complexity$H, "C" = Entropy.Complexity$C, "Color" = Color, "Texture" = Texture)
  #print(Entropy.Complexity)
  XMAX = max(Entropy.Complexity$H)
  XMIN = min(Entropy.Complexity$H) 
  YMAX = max(Entropy.Complexity$C)
  YMIN = min(Entropy.Complexity$C)
  p = ggplot(Entropy.Complexity, aes(x = Entropy.Complexity$H, y = Entropy.Complexity$C, label = Texture))+geom_point(color = Color, size = 2) + scale_x_continuous(limits=c(XMAX, XMIN)) + scale_y_continuous(limits=c(YMAX, YMIN)) + theme(plot.title = element_text(hjust=0.5))+geom_text(check_overlap = TRUE, hjust = 0, nudge_x = xpoints, size = 2) 
  return(p)
} 

get.number.patterns <- function(nrows, ncols, dimension, delay){
  if(dimension > delay){
    x = nrows - ((dimension - 1)*delay)
    y = ncols - ((dimension - 1)*delay)
  }else if(dimension < delay){
    x = nrows %/% delay
    y = ncols %/% delay
    if(nrows %% delay >= dimension) x = x + 1
    if(ncols %% delay >= dimension) y = y + 1
  }else{
    x = nrows%/%dimension
    y = ncols%/%dimension
  }
  return(x*y)
}