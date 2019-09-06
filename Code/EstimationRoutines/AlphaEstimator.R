########################################################################################################
# AlphaEstimator.R 
#
# Estimating the alpha parameter for textures of SAR images
#
# Author: Eduarda Chagas
# Date : Aug 2019
# Contact: eduarda-chagas@ufmg.br
########################################################################################################

source("../SAR/SARTimeSerie.R")
source("../Textures/TextureTimeSeries.R")
source("../Band&Pompe.R")
require(maxLik)

########################### Parameter estimation functions ############################################ 

GI0.Estimator.m1m2 <- function(z, L){
  m1 <- mean(z)
  m2 <- mean(z^2)
  m212 <- m2/m1^2
  
  a <- -2 - (L+1) / (L * m212)
  g <- m1 * (2 + (L+1) / (L * m212))
  
  return(list("alpha"=a, "gamma"=g, "looks" = L))
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

########################### Defining samples for each region analyzed #####################################
ns = 40
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


########################### Alpha parameter estimation #####################################
get.alpha <- function(img){
  z <- OneDimensionTextureHilbert(img)
  estimators.moments <- GI0.Estimator.m1m2(z,looks)
  cat("Alpha: ", estimators.moments$alpha, " Gamma: ", estimators.moments$gamma, " Looks: ", estimators.moments$looks, "\n")
  estimators <- maxNR(GI0.Estimator.LogLikelihoodLknown,
                      start=c(estimators.moments$alpha, estimators.moments$gamma, estimators.moments$looks),
                      activePar=c(TRUE,TRUE,FALSE))$estimate[1:2]
  print(estimators)
  estimators
}
SAR.Analysis.Regions <- function(option){
  
  estimators.ML = matrix(ncol = 2, nrow = ns)
  estimators.moments = array(list(), ns)
  looks = 36
  
  if(option == 1){ #Forest regions in Guatemala
    sar_data <- raster(paste("../../Data/", "guatemala", "/HHHH", ".grd", sep = ""))
    for(j in c(1:ns.guatemala)){
      img <- getValuesBlock(sar_data, row = dimen.guatemala[j,1], nrows = dimen.guatemala[j,2], col = dimen.guatemala[j,3], ncols = dimen.guatemala[j,4], format = "matrix")
      z <- OneDimensionTextureHilbert(img)
      estimators.moments[[j]] <- GI0.Estimator.m1m2(z,looks)
      estimators.ML[j,] <- maxNR(GI0.Estimator.LogLikelihoodLknown, start=c(estimators.moments[[j]]$alpha, estimators.moments[[j]]$gamma, looks), activePar=c(TRUE,TRUE,FALSE))$estimate[1:2]
    }
  }
  else if(option == 2){ #Cape Canaveral - behavior 1
    sar_data <- raster(paste("../../Data/", "cape", "/HHHH", ".grd", sep = ""))
    for(j in c(1:ns.canaveral.behavior1)){
      img <- getValuesBlock(sar_data, row = dimen.canaveral.behavior1[j,1], nrows = dimen.canaveral.behavior1[j,2], col = dimen.canaveral.behavior1[j,3], ncols = dimen.canaveral.behavior1[j,4], format = "matrix")
      z <- OneDimensionTextureHilbert(img)
      estimators.moments[[j]] <- GI0.Estimator.m1m2(z,looks)
      estimators.ML[j,] <- maxNR(GI0.Estimator.LogLikelihoodLknown, start=c(estimators.moments[[j]]$alpha, estimators.moments[[j]]$gamma, looks), activePar=c(TRUE,TRUE,FALSE))$estimate[1:2]
    }
  }
  else if(option == 3){ #Cape Canaveral - behavior 2
    sar_data <- raster(paste("../../Data/", "cape", "/HHHH", ".grd", sep = ""))
    for(j in c(1:ns.canaveral.behavior2)){
      img <- getValuesBlock(sar_data, row = dimen.canaveral.behavior2[j,1], nrows = dimen.canaveral.behavior2[j,2], col = dimen.canaveral.behavior2[j,3], ncols = dimen.canaveral.behavior2[j,4], format = "matrix")
      z <- OneDimensionTextureHilbert(img)
      estimators.moments[[j]] <- GI0.Estimator.m1m2(z,looks)
      estimators.ML[j,] <- maxNR(GI0.Estimator.LogLikelihoodLknown, start=c(estimators.moments[[j]]$alpha, estimators.moments[[j]]$gamma, looks), activePar=c(TRUE,TRUE,FALSE))$estimate[1:2]
    }
  }
  else{ #Urban regions in Munich
    sar_data <- raster(paste("../../Data/", "munich", "/HHHH", ".grd", sep = ""))
    for(j in c(1:ns.munich)){
      img <- getValuesBlock(sar_data, row = dimen.munich[j,1], nrows = dimen.munich[j,2], col = dimen.munich[j,3], ncols = dimen.munich[j,4], format = "matrix")
      z <- OneDimensionTextureHilbert(img)
      estimators.moments[[j]] <- GI0.Estimator.m1m2(z,looks)
      estimators.ML[j,] <- maxNR(GI0.Estimator.LogLikelihoodLknown, start=c(estimators.moments[[j]]$alpha, estimators.moments[[j]]$gamma, looks), activePar=c(TRUE,TRUE,FALSE))$estimate[1:2]
    }
  }
  parameters <- data.frame(alpha.moments = estimators.moments[[j]]$alpha, gamma.moments = estimators.moments[[j]]$gamma, alpha.ML = estimators.ML[,1], gamma.ML = estimators.ML[,2])
  parameters
}

alpha.forest.region = SAR.Analysis.Regions(1)
#alpha.cape.behavior1.region = SAR.Analysis.Regions(2)
#alpha.cape.behavior2.region = SAR.Analysis.Regions(3)
#alpha.urban.region = SAR.Analysis.Regions(4)
