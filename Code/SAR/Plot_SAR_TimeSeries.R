###################################################################################
# Textures_hilbert.R 
#
# Plot SAR data as one-dimensional data
#
# Author: Eduarda Chagas
# Date : Jun 2019
# Contact: eduarda-chagas@ufmg.br
####################################################################################

source("SARTimeSerie.R")
source("../Band&Pompe.R")
source("../Textures/TextureTimeSeries.R")
require(seewave)
require(plm)

ns.guatemala = 8
dimen.guatemala <- matrix(nrow = ns.guatemala, ncol = 4)
ns.canaveral = 8
dimen.canaveral <- matrix(nrow = ns.canaveral, ncol = 4)
ns.munich = 4
dimen.munich <- matrix(nrow = ns.munich, ncol = 4)

#The SAR data is available on https://drive.google.com/file/d/1jtbOcYwQfysfcUp4UhoA7lSl4_tPIqfa/view?usp=sharing and
# correspond to HHHH band of an image taken from the Cape Canaveral (acquired Sep 22, 2016)

#Ocean regions in Cape Canaveral
#{Behavior 1}
dimen.canaveral[1,] <- c(100, 128, 1700, 128)
dimen.canaveral[2,] <- c(300, 128, 1700, 128)
dimen.canaveral[3,] <- c(100, 128, 1900, 128)
dimen.canaveral[4,] <- c(300, 128, 1900, 128)
#{Behavior 2}
dimen.canaveral[5,] <- c(100, 128, 100, 128)
dimen.canaveral[6,] <- c(400, 128, 100, 128)
dimen.canaveral[7,] <- c(700, 128, 100, 128)
dimen.canaveral[8,] <- c(1100, 128, 100, 128)

#The SAR data is available on https://drive.google.com/file/d/1pO6p_UI9Cgdci9y6jVynAv8SrrAvv7K8/view?usp=sharing and
# correspond to HHHH band of an image taken from the Munich, Germany (acquired Jun 5, 2015) 

#Urban regions in Munich
dimen.munich[1,] <- c(3000, 128, 400, 128)
dimen.munich[2,] <- c(3000, 128, 600, 128)
dimen.munich[3,] <- c(3200, 128, 400, 128)
dimen.munich[4,] <- c(3200, 128, 600, 128)

#Forest regions in Guatemala
dimen.guatemala[1,] <- c(5600, 128, 2700, 128) #region 1
dimen.guatemala[2,] <- c(5200, 128, 2800, 128) #region 2
dimen.guatemala[3,] <- c(4100, 128, 2930, 128) #region 3
dimen.guatemala[4,] <- c(1075, 128, 1930, 128) #region 4

#Crop regions in Guatemala
dimen.guatemala[5,] <- c(400, 128, 1500, 128) #region 5

#Ground regions in Guatemala
#---------------------------------------------------
#Possibilly you will have problems with this regions
#It isn't uniform
dimen.guatemala[6,] <- c(250, 128, 1100, 128) #region 6
dimen.guatemala[7,] <- c(250, 128, 1850, 128) #region 7
dimen.guatemala[8,] <- c(1675, 128, 1930, 128) #region 8

SAR.timeSeries <- function(analysis, image.region, j){
  if(image.region == 1){
    sar_data <- raster(paste("../../Data/", "guatemala", "/HHHH", ".grd", sep = ""))
    img <- getValuesBlock(sar_data, row = dimen.guatemala[j,1], nrows = dimen.guatemala[j,2], col = dimen.guatemala[j,3], ncols = dimen.guatemala[j,4], format = "matrix")
  }else if(image.region == 2){
    sar_data <- raster(paste("../../Data/", "cape", "/HHHH", ".grd", sep = ""))
    img <- getValuesBlock(sar_data, row = dimen.canaveral[j,1], nrows = dimen.canaveral[j,2], col = dimen.canaveral[j,3], ncols = dimen.canaveral[j,4], format = "matrix")
  }else{
    sar_data <- raster(paste("../../Data/", "munich", "/HHHH", ".grd", sep = ""))
    img <- getValuesBlock(sar_data, row = dimen.munich[j,1], nrows = dimen.munich[j,2], col = dimen.munich[j,3], ncols = dimen.munich[j,4], format = "matrix")
  }
  
  if(analysis == 1){ # Raster 1
    time.series = as.vector(t(img))
  }else if(analysis == 2){ # Raster 2
    time.series = as.vector(img)
  }else{ # Hilbert
    hilbertcurve <- unlist(read.table("../../Data/Hilbert/HilbertCurves128.txt")) + 1
    time.series = img[hilbertcurve]
  }
  
  #qplot(x=c(1:length(time.series)),y=time.series,geom="line",xlab="Time",ylab="Serie") +
  #  theme(plot.title = element_text(hjust=0.5))
  return(time.series)
}

#Guatemala
j = 1
analysis = 3


#Domain setup
T <- 163.83 
dt <- 0.01 #s
n <- T/dt
F <-1/dt
df <- 1/T
freq<-5 #Hz
t <- seq(0,T,by=dt) #also try ts function
f <- 1:length(t)/T

#################################################################################################################

Guatemala.signal <- SAR.timeSeries(analysis, 1, j)
ts = (Guatemala.signal - min(Guatemala.signal))/(max(Guatemala.signal)- min(Guatemala.signal))

wpe.rms.guatemala = WPE.RMS(Guatemala.signal, 3, 1)
h.guatemala = shannonNormalized(wpe.rms.guatemala)
c.guatemala = Ccomplexity(wpe.guatemala)
HC.wpe.rms.guatemala = data.frame(H = h.guatemala, C = c.guatemala) 
#cotas(3) + geom_point(aes(x = HC.wpe.rms.guatemala$H, y = HC.wpe.rms.guatemala$C), size = 1) + theme(plot.title = element_text(hjust=0.5)) 

wpe.guatemala = WPE(Guatemala.signal, 3, 1)
h.guatemala = shannonNormalized(wpe.guatemala)
c.guatemala = Ccomplexity(wpe.guatemala)
HC.wpe.guatemala = data.frame(H = h.guatemala, C = c.guatemala) 
#cotas(3) + geom_point(aes(x = HC.wpe.guatemala$H, y = HC.wpe.guatemala$C), size = 1) + theme(plot.title = element_text(hjust=0.5)) 

#################################################################################################################

Canaveral.signal <- SAR.timeSeries(analysis, 2, j)
ts = (Canaveral.signal - min(Canaveral.signal))/(max(Canaveral.signal)- min(Canaveral.signal))

wpe.rms.cape = WPE.RMS(Canaveral.signal, 3, 1)
h.cape = shannonNormalized(wpe.rms.cape)
c.cape = Ccomplexity(wpe.rms.cape)
HC.wpe.rms.cape = data.frame(H = h.cape, C = c.cape) 
#cotas(3) + geom_point(aes(x = HC.wpe.rms.cape$H, y = HC.wpe.rms.cape$C), size = 1) + theme(plot.title = element_text(hjust=0.5)) 

wpe.cape = WPE(Canaveral.signal, 3, 1)
h.cape = shannonNormalized(wpe.cape)
c.cape = Ccomplexity(wpe.cape)
HC.wpe.cape = data.frame(H = h.cape, C = c.cape) 
#cotas(3) + geom_point(aes(x = HC.wpe.cape$H, y = HC.wpe.cape$C), size = 1) + theme(plot.title = element_text(hjust=0.5)) 

#################################################################################################################

Munich.signal <- SAR.timeSeries(analysis, 3, j)
ts = (Munich.signal - min(Munich.signal))/(max(Munich.signal)- min(Munich.signal))

wpe.rms.munich = WPE.RMS(Munich.signal, 3, 1)
h.munich = shannonNormalized(wpe.rms.munich)
c.munich = Ccomplexity(wpe.rms.munich)
HC.wpe.rms.munich = data.frame(H = h.munich, C = c.munich) 
#cotas(3) + geom_point(aes(x = HC.wpe.rms.munich$H, y = HC.wpe.rms.munich$C), size = 1) + theme(plot.title = element_text(hjust=0.5)) 

wpe.munich = WPE(Munich.signal, 3, 1)
h.munich = shannonNormalized(wpe.munich)
c.munich = Ccomplexity(wpe.munich)
HC.wpe.munich = data.frame(H = h.munich, C = c.munich) 
#cotas(3) + geom_point(aes(x = HC.wpe.munich$H, y = HC.wpe.munich$C), size = 1) + theme(plot.title = element_text(hjust=0.5)) 