#######################################################################################################
# RMS_R2_Analysis.R
#
# Analysis of the SAR images using RMS and R-squared 
#
# Author: Eduarda Chagas
# Date : Set 2019
# Contact: eduarda.chagas@dcc.ufmg.br
########################################################################################################

############################################# Packages #################################################
source("SARTimeSerie.R")
require(ggplot2)
require(ggthemes)

###################################### Image Sample Parameters #########################################

#Defining the number of samples from each region analyzed
ns.guatemala = 40
dimen.guatemala <- matrix(nrow = ns.guatemala, ncol = 4)
ns.canaveral.behavior1 = 40
dimen.canaveral.behavior1 <- matrix(nrow = ns.canaveral.behavior1, ncol = 4)
ns.canaveral.behavior2 = 40
dimen.canaveral.behavior2 <- matrix(nrow = ns.canaveral.behavior2, ncol = 4)
ns.munich = 40
dimen.munich <- matrix(nrow = ns.munich, ncol = 4)

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

###################################### Function of Plot ##############################################

RMS.R2.Analysis.Classes <- function(color.signal, shape.signal, signal.values){
  
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
  signal.values <- data.frame("RSM" = signal.values$RSM, "R2" = signal.values$R2, "Color" = Color, "Shape" = Shape, "Texture" = Texture)
  
  p = ggplot(signal.values, aes(x = signal.values$R2, y = signal.values$RSM)) +
    geom_point(shape = signal.values$Shape, color = signal.values$Color, size = 1) + 
    theme(plot.title = element_text(hjust=0.5))
  return(p)
} 

###################################### Function of Analysis ##########################################

SAR.RMS.R2.Analysis <- function(){
  
  hilbertcurve <- unlist(read.table("../../Data/Hilbert/HilbertCurves128.txt")) + 1
  RSM = R2 = pearson  = color = shape = rep(0, (ns.guatemala + ns.canaveral.behavior1 + ns.canaveral.behavior2 + ns.munich))
  
  #Guatemala
  sar_data <- raster(paste("../../Data/", "guatemala", "/HHHH", ".grd", sep = ""))
  for(j in c(1:ns.guatemala)){
    img <- getValuesBlock(sar_data, row = dimen.guatemala[j,1], nrows = dimen.guatemala[j,2], col = dimen.guatemala[j,3], ncols = dimen.guatemala[j,4], format = "matrix")
    ts = (img[hilbertcurve] - min(img[hilbertcurve]))/(max(img[hilbertcurve])- min(img[hilbertcurve]))
    x = hist.default(ts)$breaks[-1]
    y = hist.default(ts)$counts + 0.5
    lm.model <- lm(log(y) ~ log(x))
    R2[j] = summary(lm(log(y) ~ log(x)))$r.squared
    RSM[j] = sqrt(mean(ts^2))
    pearson[j] = cor(log(x), log(y))
    color[j] = 1
    shape[j] = 1
    cat("Guatemala ", j, "\n")
  }
  
  #Cape Canaveral - behavior 1
  sar_data <- raster(paste("../../Data/", "cape", "/HHHH", ".grd", sep = ""))
  for(j in c(1:ns.canaveral.behavior1)){
    img <- getValuesBlock(sar_data, row = dimen.canaveral.behavior1[j,1], nrows = dimen.canaveral.behavior1[j,2], col = dimen.canaveral.behavior1[j,3], ncols = dimen.canaveral.behavior1[j,4], format = "matrix")
    ts = (img[hilbertcurve] - min(img[hilbertcurve]))/(max(img[hilbertcurve])- min(img[hilbertcurve]))
    x = hist.default(ts)$breaks[-1]
    y = hist.default(ts)$counts + 0.5
    lm.model <- lm(log(y) ~ log(x))
    R2[ns.guatemala + j] = summary(lm(log(y) ~ log(x)))$r.squared
    RSM[ns.guatemala + j] = sqrt(mean(ts^2))
    pearson[ns.guatemala + j] = cor(log(x), log(y))
    color[ns.guatemala + j] = 2
    shape[ns.guatemala + j] = 2
    cat("Canaveral 1 ", j, "\n")
  }
  
  #Cape Canaveral - behavior 2
  sar_data <- raster(paste("../../Data/", "cape", "/HHHH", ".grd", sep = ""))
  for(j in c(1:ns.canaveral.behavior2)){
    img <- getValuesBlock(sar_data, row = dimen.canaveral.behavior2[j,1], nrows = dimen.canaveral.behavior2[j,2], col = dimen.canaveral.behavior2[j,3], ncols = dimen.canaveral.behavior2[j,4], format = "matrix")
    ts = (img[hilbertcurve] - min(img[hilbertcurve]))/(max(img[hilbertcurve])- min(img[hilbertcurve]))
    x = hist.default(ts)$breaks[-1]
    y = hist.default(ts)$counts + 0.5
    lm.model <- lm(log(y) ~ log(x))
    R2[(ns.canaveral.behavior1 + ns.guatemala) + j] = summary(lm(log(y) ~ log(x)))$r.squared
    RSM[(ns.canaveral.behavior1 + ns.guatemala) + j] = sqrt(mean(ts^2))
    pearson[(ns.canaveral.behavior1 + ns.guatemala) + j] = cor(log(x), log(y))
    color[(ns.canaveral.behavior1 + ns.guatemala) + j] = 3
    shape[(ns.canaveral.behavior1 + ns.guatemala) + j] = 2
    cat("Canaveral 2 ", j, "\n")
  }
  
  #Munich
  sar_data <- raster(paste("../../Data/", "munich", "/HHHH", ".grd", sep = ""))
  for(j in c(1:ns.munich)){
    img <- getValuesBlock(sar_data, row = dimen.munich[j,1], nrows = dimen.munich[j,2], col = dimen.munich[j,3], ncols = dimen.munich[j,4], format = "matrix")
    ts = (img[hilbertcurve] - min(img[hilbertcurve]))/(max(img[hilbertcurve])- min(img[hilbertcurve]))
    x = hist.default(ts)$breaks[-1]
    y = hist.default(ts)$counts + 0.5
    lm.model <- lm(log(y) ~ log(x))
    R2[(ns.canaveral.behavior1 + ns.canaveral.behavior2 + ns.guatemala) + j] = summary(lm(log(y) ~ log(x)))$r.squared
    RSM[(ns.canaveral.behavior1 + ns.canaveral.behavior2 + ns.guatemala) + j] = sqrt(mean(ts^2))
    pearson[(ns.canaveral.behavior1 + ns.canaveral.behavior2 + ns.guatemala) + j] = cor(log(x), log(y))
    color[(ns.canaveral.behavior1 + ns.canaveral.behavior2 + ns.guatemala) + j] = 4
    shape[(ns.canaveral.behavior1 + ns.canaveral.behavior2 + ns.guatemala) + j] = 3
    cat("Munich ", j, "\n")
  }
  
  signal.values <- data.frame(RSM, R2)
  
  p <- RMS.R2.Analysis.Classes(color, shape, signal.values) + theme(axis.title.x=element_blank(),
                                                                    axis.title.y=element_blank(), 
                                                                    legend.position = "none")
  plot(p) + ggtitle(expression(italic("SAR Images - R2 and RMS"))) +
    xlab(expression(italic(R2))) + ylab(expression(italic(RMS))) + labs(colour=expression(italic(Regions))) +
    theme_igray() + theme(text=element_text(size=14, family="Times New Roman"), axis.text.x=element_blank(), axis.text.y=element_blank(),plot.title = element_text(hjust=0.5)) + 
    guides(colour = guide_legend(override.aes = list(size=3)))
}
