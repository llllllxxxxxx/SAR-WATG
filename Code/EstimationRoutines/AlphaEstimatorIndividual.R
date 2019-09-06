source("../SAR/SARTimeSerie.R")
source("../Textures/TextureTimeSeries.R")
source("../Band&Pompe.R")
require(maxLik)
require(stats4)

GI0.Estimator.m1m2 <- function(z, L){
  m1 <- mean(z)
  m2 <- mean(z^2)
  m212 <- m2/m1^2
  
  a <- -2 - (L+1) / (L * m212)
  g <- m1 * (2 + (L+1) / (L * m212))
  
  return(list("alpha"=a, "gamma"=g))
}

GI0.Estimator.LogLikelihoodLknown <- function(params){
  p_alpha <- -abs(params[1])
  p_gamma <- abs(params[2])
  p_L <- abs(params[3])
  
  n <- length(z)
  
  return(
    n*(lgamma(p_L-p_alpha) - p_alpha*log(p_gamma) - lgamma(-p_alpha)) + 
      (p_alpha-p_L)*sum(log(p_gamma + z*p_L)) 
  )
}

densGI0 <- function(x, alpha, gamma, Looks) {
  return(((Looks^Looks * gamma(Looks-alpha))/ (gamma^alpha * gamma(Looks) * gamma(-alpha))) *
           (x^(Looks-1) / (gamma + Looks * x)^(Looks-alpha)))
}


looks = 36

SAR.Analysis.Regions.individual <- function(option){
  ns = 1
  if(option == 1){ #Forest regions in Guatemala
    dimen.guatemala <- c(5150, 128, 2700, 128)
    sar_data <- raster(paste("../../Data/", "guatemala", "/HHHH", ".grd", sep = ""))
    img <- getValuesBlock(sar_data, row = dimen.guatemala[1], nrows = dimen.guatemala[2], col = dimen.guatemala[3], ncols = dimen.guatemala[4], format = "matrix")
    
  }else if(option == 2){ #Ocean regions in Cape Canaveral
    dimen.canaveral <- c(250, 128, 1, 128)
    sar_data <- raster(paste("../../Data/", "cape", "/HHHH", ".grd", sep = ""))
    img <- getValuesBlock(sar_data, row = dimen.canaveral[1], nrows = dimen.canaveral[2], col = dimen.canaveral[3], ncols = dimen.canaveral[4], format = "matrix")
    
  }else if(option == 3){ #Urban regions in Munich
    dimen.munich <- c(3000, 128, 400, 128)
    sar_data <- raster(paste("../../Data/", "munich", "/HHHH", ".grd", sep = ""))
    img <- getValuesBlock(sar_data, row = dimen.munich[1], nrows = dimen.munich[2], col = dimen.munich[3], ncols = dimen.munich[4], format = "matrix")
    
  }
  img
}

img <- SAR.Analysis.Regions.individual(2)
z <- OneDimensionTextureHilbert(img)
#w <- 2*ns.guatemala^(-1/3)*IQR(z)
#z <- data.frame(z)
#ggplot(z, aes(x=z)) + geom_histogram()
#ggplot(z, aes(x=z)) + geom_histogram(binwidth = w)
estimators.moments <- GI0.Estimator.m1m2(z,looks)
estimators.moments
estimators.ML <- maxNR(GI0.Estimator.LogLikelihoodLknown, start=c(estimators.moments$alpha, estimators.moments$gamma, looks), activePar=c(TRUE,TRUE,FALSE))$estimate[1:2]
estimators.ML


dens1 <- densGI0(z, estimators.moments$alpha, estimators.moments$gamma, looks)
dens1 <- densGI0(z, estimators.ML[1], estimators.ML[2], looks)

qplot (xlab="Intensity", ylab="GI0 Densities") +
  ggtitle("Probability density function of a GI0 distribution") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_line(aes(z,dens1, colour="black"), size=1.3)