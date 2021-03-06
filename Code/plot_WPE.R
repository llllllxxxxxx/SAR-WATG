########################################################################################################
# plot_WPE.R
#
# Generate Plot of Analysis of SAR images using WPE
#
# Author: Eduarda Chagas
# Date : May 9, 2020
# Contact: eduarda.chagas@dcc.ufmg.br
########################################################################################################

############################################# Packages #################################################

setwd("/home/eduarda/Desktop/Research/Repositories/PolSARfromITQualitative/Code")
if(!require(gtools)) install.packages("gtools")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(ggthemes)) install.packages("ggthemes")
if(!require(ggpubr)) install.packages("ggpubr")
source("theory_information.R")

###################################### Function of Plot ################################################

HC.Plane.no.cota <- function(dimension, color.signal, shape.signal, signal.values){
  
  shape.select = c(17,18,19,8)
  XMIN = min(signal.values[,1]) + 0.00005
  XMAX = min(max(signal.values[,1]) + 0.0005, 1)
  YMIN = max(0,min(signal.values[,2]) - 0.0005)
  YMAX = max(signal.values[,2]) + 0.0005
  
  # Paleta montada a partir de https://coolors.co/
  rainbow.colors = palette(c("#3F84E5",
                             "#B20D30", 
                             "#3F784C",
                             "#EFCD4A"))
  
  Color = rainbow.colors[color.signal]
  Shape = shape.select[shape.signal]
  Regions =  c("Forest", "Ocean", "", "Urban")[color.signal]
  signal.values = data.frame("H" = signal.values[,1], "C" = signal.values[,2], "Color" = Color, "Shape" = Shape, "Regions" = Regions)
  
  p = ggplot(signal.values, aes(H, C, color = Regions, shape = Shape)) + geom_point(size = 2) +
    scale_shape_identity() +
    xlim(limits=c(XMIN, XMAX)) + ylim(limits=c(YMIN, YMAX)) + 
    theme_few(base_size = 18, base_family = "serif")  + 
    theme(plot.title = element_text(hjust=0.5)) + 
    scale_colour_few("Dark")
  print(p)
  return(p)
}


###################################### Function of Analysis ##########################################

plot.WPE.d3t1 <- function(){
  
  n = 3 #Dimension parameter
  tal = 1 #Delay parameter
  plots = array(list(), 20)
  hilbertcurve = unlist(read.table("../Data/Hilbert/HilbertCurves128.txt")) + 1
  types = c(rep(1,40), rep(2,80), rep(4, 40))
  regions = c(rep(1,40), rep(2,80), rep(4, 40))
  n.total = 160
  
  Entropy.Complexity.csv = read.csv(file="../Data/EntropyComplexityWPED3T1.csv", header=TRUE, sep=",")
  Entropy.Complexity = matrix(nrow = n.total, ncol = 2)
  
  Entropy.Complexity[,1] = Entropy.Complexity.csv[, 1]
  Entropy.Complexity[,2] = Entropy.Complexity.csv[, 2]
  plot.WPE = HC.Plane.no.cota(n, regions, types, Entropy.Complexity) + ggtitle(expression(italic("WPE"))) +
    labs(x="", y="")
  return(plot.WPE)
}

pdf("WPE.pdf", width = 10, height = 8) 
plot.WPE = plot.WPE.d3t1()
dev.off() 