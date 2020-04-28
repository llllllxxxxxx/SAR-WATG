###################################################################################
# Plot_sar_timeSeries.R 
#
# Generate Plot of SAR data as one-dimensional data
#
# Author: Eduarda Chagas
# Date : Apr 2020
# Contact: eduarda.chagas@dcc.ufmg.br
####################################################################################
#setwd("/home/eduarda/Desktop/Research/Repositories/PolSARfromITQualitative/Code")

if(!require(ggpubr)) install.packages("ggpubr")
if(!require(raster)) install.packages("raster")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(devtools)) install.packages("devtools")
if(!require(ggthemes)) install.packages("ggthemes")


ns.munich = 4
ns.guatemala = 4
ns.canaveral = 8
dimen.munich = matrix(nrow = ns.munich, ncol = 4)
dimen.guatemala = matrix(nrow = ns.guatemala, ncol = 4)
dimen.canaveral = matrix(nrow = ns.canaveral, ncol = 4)

#Ocean regions in Cape Canaveral
#{Behavior 1}
dimen.canaveral[1,] = c(100, 128, 1700, 128)
dimen.canaveral[2,] = c(300, 128, 1700, 128)
dimen.canaveral[3,] = c(100, 128, 1900, 128)
dimen.canaveral[4,] = c(300, 128, 1900, 128)
#{Behavior 2}
dimen.canaveral[5,] = c(100, 128, 100, 128)
dimen.canaveral[6,] = c(400, 128, 100, 128)
dimen.canaveral[7,] = c(700, 128, 100, 128)
dimen.canaveral[8,] = c(1100, 128, 100, 128)

#Urban regions in Munich
dimen.munich[1,] = c(3000, 128, 400, 128)
dimen.munich[2,] = c(3000, 128, 600, 128)
dimen.munich[3,] = c(3200, 128, 400, 128)
dimen.munich[4,] = c(3200, 128, 600, 128)

#Forest regions in Guatemala
dimen.guatemala[1,] = c(5600, 128, 2700, 128) #region 1
dimen.guatemala[2,] = c(5200, 128, 2800, 128) #region 2
dimen.guatemala[3,] = c(4100, 128, 2930, 128) #region 3
dimen.guatemala[4,] = c(1075, 128, 1930, 128) #region 4

SAR.timeSeries <- function(analysis, image.region, j){
  if(image.region == 1){
    sar_data = raster(paste("../Data/", "guatemala", "/HHHH", ".grd", sep = ""))
    img = getValuesBlock(sar_data, row = dimen.guatemala[j,1], nrows = dimen.guatemala[j,2], col = dimen.guatemala[j,3], ncols = dimen.guatemala[j,4], format = "matrix")
  }else if(image.region == 2){
    sar_data = raster(paste("../Data/", "cape", "/HHHH", ".grd", sep = ""))
    img = getValuesBlock(sar_data, row = dimen.canaveral[j,1], nrows = dimen.canaveral[j,2], col = dimen.canaveral[j,3], ncols = dimen.canaveral[j,4], format = "matrix")
  }else{
    sar_data = raster(paste("../Data/", "munich", "/HHHH", ".grd", sep = ""))
    img = getValuesBlock(sar_data, row = dimen.munich[j,1], nrows = dimen.munich[j,2], col = dimen.munich[j,3], ncols = dimen.munich[j,4], format = "matrix")
  }
  
  if(analysis == 1){ # Raster 1
    time.series = as.vector(t(img))
  }else if(analysis == 2){ # Raster 2
    time.series = as.vector(img)
  }else{ # Hilbert
    hilbertcurve = unlist(read.table("../Data/Hilbert/HilbertCurves128.txt")) + 1
    time.series = img[hilbertcurve]
  }
  return(time.series)
}

n = 3
tal = 1
analysis = 3

#################################################################################################################

Guatemala.signal = SAR.timeSeries(analysis, 1, j = 1)
Guatemala.signal = (Guatemala.signal - min(Guatemala.signal))/(max(Guatemala.signal)- min(Guatemala.signal))

#watg.guatemala = weight.transition.graph(Guatemala.signal, n, tal)
#round(watg.guatemala, digits = 3)

#tg.guatemala = transition.graphs(Guatemala.signal, n, tal)
#round(tg.guatemala, digits = 3)

guatemala = qplot(x = c(1:length(Guatemala.signal)), y = Guatemala.signal, geom = "line", main = "Guatemala Forest", xlab = "Observation", ylab = "Serie") +
  theme_few(base_size = 10, base_family = "serif")  + theme(plot.title = element_text(hjust=0.5)) + 
  scale_colour_few("Dark")

#################################################################################################################

Canaveral.signal.1 = SAR.timeSeries(analysis = 3, image.region = 2, j = 1)
Canaveral.signal.1 = (Canaveral.signal.1 - min(Canaveral.signal.1))/(max(Canaveral.signal.1)- min(Canaveral.signal.1))

#watg.cape = weight.transition.graph(Canaveral.signal, n, tal)
#round(watg.cape, digits = 3)

#tg.cape = transition.graphs(Canaveral.signal, n, tal)
#round(tg.cape, digits = 3)

canaveral.1 = qplot(x = c(1:length(Canaveral.signal.1)), y = Canaveral.signal.1, geom = "line", main = "Canaveral Ocean - Contrast type 1", xlab = "Observation", ylab = "Serie") +
  theme_few(base_size = 10, base_family = "serif")  + theme(plot.title = element_text(hjust=0.5)) + 
  scale_colour_few("Dark")

#################################################################################################################

Canaveral.signal = SAR.timeSeries(analysis = 3, image.region = 2, j = 7)
Canaveral.signal = (Canaveral.signal - min(Canaveral.signal))/(max(Canaveral.signal)- min(Canaveral.signal))

#watg.cape = weight.transition.graph(Canaveral.signal, n, tal)
#round(watg.cape, digits = 3)

#tg.cape = transition.graphs(Canaveral.signal, n, tal)
#round(tg.cape, digits = 3)

canaveral.2 = qplot(x = c(1:length(Canaveral.signal)), y = Canaveral.signal, geom = "line", main = "Canaveral Ocean - Contrast type 2", xlab = "Observation", ylab = "Serie") +
  theme_few(base_size = 10, base_family = "serif")  + theme(plot.title = element_text(hjust=0.5)) + 
  scale_colour_few("Dark")

#################################################################################################################

Munich.signal = SAR.timeSeries(analysis = 3, image.region = 3, j = 1)
Munich.signal = (Munich.signal - min(Munich.signal))/(max(Munich.signal)- min(Munich.signal))

#watg.munich = weight.transition.graph(Munich.signal, n, tal)
#round(watg.munich, digits = 3)

#tg.munich = transition.graphs(Munich.signal, n, tal)
#round(tg.munich, digits = 3)

munich = qplot(x = c(1:length(Munich.signal)), y = Munich.signal, geom = "line", main = "Munich Urban Area", xlab = "Observation", ylab = "Serie") +
  theme_few(base_size = 10, base_family = "serif")  + theme(plot.title = element_text(hjust=0.5)) + 
  scale_colour_few("Dark")

#################################################################################################################

ggarrange(guatemala, canaveral.1, canaveral.2, munich + rremove("x.text"), ncol = 4, nrow = 1)

########################################### Plot Weigthed Graphs ################################################

if(!require(igraph)) install.packages("igraph")

nodes.names = c("123", "132", "213", "231", "312", "321")
#nodes.names.pi = c(expression(paste(pi, "1")), expression(paste(pi, "2")), expression(paste(pi, "3")), expression(paste(pi, "4")), expression(paste(pi, "5")), expression(paste(pi, "6")))

g1 = graph_from_adjacency_matrix(g.cape, mode = "directed", weighted = TRUE)

plot(g1, layout = layout_with_dh(g1),
     edge.arrow.size = .4, edge.curved = 0.1, edge.color = "#657D81", edge.label = round(E(g1)$weight, 3), edge.label.cex = 0.8,
     label.x = 1.5, label.y = 5, edge.label.font = 2,
     vertex.color = "#C1D1E1", vertex.frame.color = "#C1D1E1",
     vertex.label = nodes.names, vertex.label.color = "#294C60",
     vertex.label.cex=.7, vertex.size = 30) 
