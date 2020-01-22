########################################################################################################
# Transition_graph.R
#
# Analysis of SAR images using Hilbert space-filling curves and ordinal patterns transition graphs
# (sliding window)
#
# Author: Eduarda Chagas
# Date : Oct 2019
# Contact: eduarda-chagas@ufmg.br
########################################################################################################

############################################# Packages #################################################

if(!require(gtools)) install.packages("gtools")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(ggthemes)) install.packages("ggthemes")
if(!require(igraph)) install.packages("igraph")
if(!require(ggpubr)) install.packages("ggpubr")
source("Theory_Information.R")

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
  33
  - `Band&Pompe.R`- 
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

transition.graph.analysis <- function(){
  
  a = b = 0
  n = c(3,4,5,6) #Dimension parameter
  tal = c(1,2,3,4,5) #Delay parameter
  plots = array(list(), 20)
  hilbertcurve <- unlist(read.table("../../Data/Hilbert/HilbertCurves128.txt")) + 1
  types <- c(rep(1,40), rep(2,40), rep(3,40), rep(4, 40))
  regions <- c(rep(1,40), rep(2,80), rep(4, 40))
  
  for(i in 1:(length(n)*length(tal))){
    cat("- Plane: ", i, "de 20 ", "\n")
    if(i%%5 == 1){
      a = a + 1
      b = 0
    }
    b = b + 1
    
    Entropy.Complexity <- matrix(nrow = n.total, ncol = 2)
    
    #Guatemala
    sar_data <- raster(paste("../../Data/", "guatemala", "/HHHH", ".grd", sep = ""))
    for(j in c(1:ns.guatemala)){
      img <- getValuesBlock(sar_data, row = dimen.guatemala[j,1], nrows = dimen.guatemala[j,2], col = dimen.guatemala[j,3], ncols = dimen.guatemala[j,4], format = "matrix")
      ts = img[hilbertcurve]/max(img[hilbertcurve])
      e = sliding.window(ts, n[a], tal[b])
      p = formation.pattern(e)
      w = pattern.wedding(p)
      g = transition.graph(e, w, n[a])
      Entropy.Complexity[j, 1] <- shannonNormalized(as.vector(g))
      Entropy.Complexity[j, 2] <- Ccomplexity(as.vector(g))
      cat("Guatemala ", j, "\n")
    }
    #Cape Canaveral - behavior 1
    sar_data <- raster(paste("../../Data/", "cape", "/HHHH", ".grd", sep = ""))
    for(j in c(1:ns.canaveral.behavior1)){
      img <- getValuesBlock(sar_data, row = dimen.canaveral.behavior1[j,1], nrows = dimen.canaveral.behavior1[j,2], col = dimen.canaveral.behavior1[j,3], ncols = dimen.canaveral.behavior1[j,4], format = "matrix")
      ts = img[hilbertcurve]/max(img[hilbertcurve])
      e = sliding.window(ts, n[a], tal[b])
      p = formation.pattern(e)
      w = pattern.wedding(p)
      g = transition.graph(e, w, n[a])
      Entropy.Complexity[ns.guatemala + j, 1] <- shannonNormalized(as.vector(g))
      Entropy.Complexity[ns.guatemala + j, 2] <- Ccomplexity(as.vector(g))
      cat("Cape 1 ", j, "\n")
    }
    #Cape Canaveral - behavior 2
    sar_data <- raster(paste("../../Data/", "cape", "/HHHH", ".grd", sep = ""))
    for(j in c(1:ns.canaveral.behavior2)){
      img <- getValuesBlock(sar_data, row = dimen.canaveral.behavior2[j,1], nrows = dimen.canaveral.behavior2[j,2], col = dimen.canaveral.behavior2[j,3], ncols = dimen.canaveral.behavior2[j,4], format = "matrix")
      ts = img[hilbertcurve]/max(img[hilbertcurve])
      e = sliding.window(ts, n[a], tal[b])
      p = formation.pattern(e)
      w = pattern.wedding(p)
      g = transition.graph(e, w, n[a])
      Entropy.Complexity[(ns.canaveral.behavior1 + ns.guatemala) + j, 1] <- shannonNormalized(as.vector(g))
      Entropy.Complexity[(ns.canaveral.behavior1 + ns.guatemala) + j, 2] <- Ccomplexity(as.vector(g))
      cat("Cape 2 ", j, "\n")
    }
    #Munich
    sar_data <- raster(paste("../../Data/", "munich", "/HHHH", ".grd", sep = ""))
    for(j in c(1:ns.munich)){
      img <- getValuesBlock(sar_data, row = dimen.munich[j,1], nrows = dimen.munich[j,2], col = dimen.munich[j,3], ncols = dimen.munich[j,4], format = "matrix")
      ts = img[hilbertcurve]/max(img[hilbertcurve])
      e = sliding.window(ts, n[a], tal[b])
      p = formation.pattern(e)
      w = pattern.wedding(p)
      g = transition.graph(e, w, n[a])
      Entropy.Complexity[(ns.canaveral.behavior1 + ns.canaveral.behavior2 + ns.guatemala) + j, 1] <- shannonNormalized(as.vector(g))
      Entropy.Complexity[(ns.canaveral.behavior1 + ns.canaveral.behavior2 + ns.guatemala) + j, 2] <- Ccomplexity(as.vector(g))
      cat("Munich ", j, "\n")
    }
    
    XMIN = min(Entropy.Complexity[,1], XMIN) 
    YMAX = max(Entropy.Complexity[,2], YMAX)
    
    name = paste("transitionGraphD",n[a],"t",tal[b], ".png", sep="")
    png(name, width = 1500, height = 850)
    
    plots[[i]] = HC.color.shape.signal(regions, types, Entropy.Complexity)
    
    print(plots[[i]])
    dev.off()
    
    #V(net)$color <- "#C9C5C5"  
    #V(net)$size <- 40
    #E(net)$width <- 2
    #E(net)$arrow.size <- .1
    #E(net)$edge.color <- "#838181"
    #E(net)$width <- 1
    #E(net)$label <- floor(100*E(net)$weight)/100
    #plot(net, vertex.label = V(net)$name)
  }
  
  for(i in 1:(length(n)*length(tal))){
    plots[[i]] = plots[[i]] + xlim(limits=c(XMIN, XMAX)) + ylim(limits=c(YMIN, YMAX)) 
  }
  
  
  png("transitionGraphHilbert.png", width = 1500, height = 850)
  
  p = ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]],
            plots[[6]], plots[[7]], plots[[8]], plots[[9]], plots[[10]],
            plots[[11]], plots[[12]], plots[[13]], plots[[14]], plots[[15]],
            plots[[16]], plots[[17]], plots[[18]], plots[[19]], plots[[20]],
            ncol=5, nrow=4, common.legend = TRUE, legend = "right") + 
    ggtitle(expression(italic("SAR Images - Transition Graph - Sliding Window"))) +
    xlab(expression(italic(H))) + ylab(expression(italic(C))) + labs(colour=expression(italic(Regions))) +
    theme_igray() + theme(text=element_text(size=14, family="Times New Roman"), axis.text.x=element_blank(), axis.text.y=element_blank(),plot.title = element_text(hjust=0.5)) + 
    guides(colour = guide_legend(override.aes = list(size=3)))
  
  print(p)
  dev.off()
  
}

transition.graph.analysis()