########################################################################################################
# SAR_Hilbert.R 
#
# Analysis of SAR images using Hilbert space-filling curves
#
# Author: Eduarda Chagas
# Date : Aug 2019
# Contact: eduarda-chagas@ufmg.br
########################################################################################################

source("SARTimeSerie.R")
source("../Textures/TextureTimeSeries.R")
source("../Band&Pompe.R")
require(grid)
require(gridExtra)
require(ggsci)
require(extrafont)
require(ggpubr)
require(dbscan)
require(DataDriven)
require(gtools)

#Defining the number of samples from each region analyzed
ns.guatemala = 40
dimen.guatemala <- matrix(nrow = ns.guatemala, ncol = 4)
ns.canaveral.behavior1 = 40
dimen.canaveral.behavior1 <- matrix(nrow = ns.canaveral.behavior1, ncol = 4)
ns.canaveral.behavior2 = 40
dimen.canaveral.behavior2 <- matrix(nrow = ns.canaveral.behavior2, ncol = 4)
ns.munich = 40
dimen.munich <- matrix(nrow = ns.munich, ncol = 4)

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

#Crop regions in Guatemala
#dimen.guatemala[5,] <- c(400, 128, 1500, 128) #region 5

#Ground regions in Guatemala
#---------------------------------------------------
#Possibilly you will have problems with this regions
#It isn't uniform
#dimen.guatemala[6,] <- c(250, 128, 1100, 128) #region 6
#dimen.guatemala[7,] <- c(250, 128, 1850, 128) #region 7
#dimen.guatemala[8,] <- c(1675, 128, 1930, 128) #region 8


a = b = 0
n = c(3,4,5,6) #Dimension parameter
tal = c(1,2,3,4,5) #Delay parameter
plots = array(list(), 20)

SAR.Hilbert.Analysis <- function(Analysis){
  
  if(Analysis == 1){ #Analysis of SAR characterization
    
    XMAX = 1
    XMIN = 1
    YMAX = 0
    YMIN = 0
    
    for(i in 1:(length(n)*length(tal))){
      cat("- Plane: ", i, "de 20 ", "\n")
      if(i%%5 == 1){
        a = a + 1
        b = 0
      }
      b = b + 1
      symbols <- definePatterns(n[a])
      probability <- matrix(nrow = (ns.guatemala + ns.canaveral.behavior1 + ns.canaveral.behavior2 + ns.munich), ncol = factorial(n[a]) + 2)
      #GUatemala
      sar_data <- raster(paste("../../Data/", "guatemala", "/HHHH", ".grd", sep = ""))
      for(j in c(1:ns.guatemala)){
        img <- getValuesBlock(sar_data, row = dimen.guatemala[j,1], nrows = dimen.guatemala[j,2], col = dimen.guatemala[j,3], ncols = dimen.guatemala[j,4], format = "matrix")
        get.serie <- formationPatternsTexturesHilbert(img, n[a], tal[b])
        probability[j, 1:factorial(n[a])] <- Bandt.Pompe(get.serie, n[a], dim(get.serie)[1])
      }
      #Cape Canaveral - behavior 1
      sar_data <- raster(paste("../../Data/", "cape", "/HHHH", ".grd", sep = ""))
      for(j in c(1:ns.canaveral.behavior1)){
        img <- getValuesBlock(sar_data, row = dimen.canaveral.behavior1[j,1], nrows = dimen.canaveral.behavior1[j,2], col = dimen.canaveral.behavior1[j,3], ncols = dimen.canaveral.behavior1[j,4], format = "matrix")
        get.serie <- formationPatternsTexturesHilbert(img, n[a], tal[b])
        probability[ns.guatemala + j, 1:factorial(n[a])] <- Bandt.Pompe(get.serie, n[a], dim(get.serie)[1])
      }
      #Cape Canaveral - behavior 2
      sar_data <- raster(paste("../../Data/", "cape", "/HHHH", ".grd", sep = ""))
      for(j in c(1:ns.canaveral.behavior2)){
        img <- getValuesBlock(sar_data, row = dimen.canaveral.behavior2[j,1], nrows = dimen.canaveral.behavior2[j,2], col = dimen.canaveral.behavior2[j,3], ncols = dimen.canaveral.behavior2[j,4], format = "matrix")
        get.serie <- formationPatternsTexturesHilbert(img, n[a], tal[b])
        probability[(ns.canaveral.behavior1 + ns.guatemala) + j, 1:factorial(n[a])] <- Bandt.Pompe(get.serie, n[a], dim(get.serie)[1])
      }
      #Munich
      sar_data <- raster(paste("../../Data/", "munich", "/HHHH", ".grd", sep = ""))
      for(j in c(1:ns.munich)){
        img <- getValuesBlock(sar_data, row = dimen.munich[j,1], nrows = dimen.munich[j,2], col = dimen.munich[j,3], ncols = dimen.munich[j,4], format = "matrix")
        get.serie <- formationPatternsTexturesHilbert(img, n[a], tal[b])
        probability[(ns.canaveral.behavior1 + ns.canaveral.behavior2 + ns.guatemala) + j, 1:factorial(n[a])] <- Bandt.Pompe(get.serie, n[a], dim(get.serie)[1])
      }
      
      #Diferent regions 
      probability[1:40, factorial(n[a]) + 1] = 1
      probability[41:80, factorial(n[a]) + 1] = 2
      probability[81:120, factorial(n[a]) + 1] = 3
      probability[121:160, factorial(n[a]) + 1] = 4
      
      #Diferent places
      probability[1:80, factorial(n[a]) + 2] = 1      #Guatemala
      probability[81:120, factorial(n[a]) + 2] = 2     #Canaveral
      probability[121:160, factorial(n[a]) + 2] = 3    #Munich
      
      Entropy.Complexity <- Entropy.Complexity.Points(probability, dim(probability)[1], n[a])
      
      XMIN = min(Entropy.Complexity$H, XMIN) 
      YMAX = max(Entropy.Complexity$C, YMAX)
      
      plots[[i]] <- HCPlane.Classes(probability, n[a], Entropy.Complexity) + theme(axis.title.x=element_blank(),
                                                               axis.title.y=element_blank(), 
                                                               legend.position = "none")
    }
    
    for(i in 1:(length(n)*length(tal))){
      plots[[i]] = plots[[i]] + xlim(limits=c(XMIN, XMAX)) + ylim(limits=c(YMIN, YMAX)) 
    }
    
    ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]],
              plots[[6]], plots[[7]], plots[[8]], plots[[9]], plots[[10]],
              plots[[11]], plots[[12]], plots[[13]], plots[[14]], plots[[15]],
              plots[[16]], plots[[17]], plots[[18]], plots[[19]], plots[[20]],
              #plots[[21]], plots[[22]], plots[[23]], plots[[24]], plots[[25]],
              ncol=5, nrow=4, common.legend = TRUE, legend = "right") + 
      ggtitle(expression(italic("SAR Images - Hilbert Curve Analysis"))) +
      xlab(expression(italic(H))) + ylab(expression(italic(C))) + labs(colour=expression(italic(Regions))) +
      theme_igray() + theme(text=element_text(size=14, family="Times New Roman"), axis.text.x=element_blank(), axis.text.y=element_blank(),plot.title = element_text(hjust=0.5)) + 
      guides(colour = guide_legend(override.aes = list(size=3)))
  }
  else if(Analysis == 2){ #Analysis of Guatemala region
    
    XMAX = 1
    XMIN = 1
    YMAX = 0
    YMIN = 0
    
    for(i in 1:(length(n)*length(tal))){
      cat("- Plane: ", i, "de 20 ", "\n")
      if(i%%5 == 1){
        a = a + 1
        b = 0
      }
      b = b + 1
      symbols <- definePatterns(n[a])
      probability <- matrix(nrow = (ns.guatemala), ncol = factorial(n[a]) + 2)
      #GUatemala
      sar_data <- raster(paste("../../Data/", "guatemala", "/HHHH", ".grd", sep = ""))
      for(j in c(1:ns.guatemala)){
        img <- getValuesBlock(sar_data, row = dimen.guatemala[j,1], nrows = dimen.guatemala[j,2], col = dimen.guatemala[j,3], ncols = dimen.guatemala[j,4], format = "matrix")
        get.serie <- formationPatternsTexturesHilbert(img, n[a], tal[b])
        probability[j, 1:factorial(n[a])] <- Bandt.Pompe(get.serie, n[a], dim(get.serie)[1])
      }
      
      #Diferent regions 
      probability[, factorial(n[a]) + 1] = 1
      
      #Diferent places
      probability[, factorial(n[a]) + 2] = 1      #Guatemala
      
      Entropy.Complexity <- Entropy.Complexity.Points(probability, dim(probability)[1], n[a])
      
      XMIN = min(Entropy.Complexity$H, XMIN) 
      YMAX = max(Entropy.Complexity$C, YMAX)
      
      plots[[i]] <- HCPlane.Classes(probability, n[a], Entropy.Complexity) + theme(axis.title.x=element_blank(),
                                                               axis.title.y=element_blank(), legend.position = "none")
    }
    
    for(i in 1:(length(n)*length(tal))){
      plots[[i]] = plots[[i]] + xlim(limits=c(XMIN, XMAX)) + ylim(limits=c(YMIN, YMAX)) 
    }
    
    ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]],
              plots[[6]], plots[[7]], plots[[8]], plots[[9]], plots[[10]],
              plots[[11]], plots[[12]], plots[[13]], plots[[14]], plots[[15]],
              plots[[16]], plots[[17]], plots[[18]], plots[[19]], plots[[20]],
              #plots[[21]], plots[[22]], plots[[23]], plots[[24]], plots[[25]],
              ncol=5, nrow=4, common.legend = TRUE, legend = "right") + 
      ggtitle(expression(italic("Guatemala Region - Hilbert Curve Analysis"))) +
      xlab(expression(italic(H))) + ylab(expression(italic(C))) + labs(colour=expression(italic(Regions))) +
      theme_igray() + theme(text=element_text(size=14, family="Times New Roman"), axis.text.x=element_blank(), axis.text.y=element_blank(),plot.title = element_text(hjust=0.5)) + 
      guides(colour = guide_legend(override.aes = list(size=3)))
  }
  else if(Analysis == 3){ #Analysis of Canaveral region
    
    XMAX = 1
    XMIN = 1
    YMAX = 0
    YMIN = 0
    
    for(i in 1:(length(n)*length(tal))){
      cat("- Plane: ", i, "de 20 ", "\n")
      if(i%%5 == 1){
        a = a + 1
        b = 0
      }
      b = b + 1
      symbols <- definePatterns(n[a])
      probability <- matrix(nrow = (ns.canaveral.behavior1 + ns.canaveral.behavior2), ncol = factorial(n[a]) + 2)
      #Cape Canaveral - behavior 1
      sar_data <- raster(paste("../../Data/", "cape", "/HHHH", ".grd", sep = ""))
      for(j in c(1:ns.canaveral.behavior1)){
        img <- getValuesBlock(sar_data, row = dimen.canaveral.behavior1[j,1], nrows = dimen.canaveral.behavior1[j,2], col = dimen.canaveral.behavior1[j,3], ncols = dimen.canaveral.behavior1[j,4], format = "matrix")
        img <- getValuesBlock(sar_data, row = dimen.canaveral.behavior1[j,1], nrows = dimen.canaveral.behavior1[j,2], col = dimen.canaveral.behavior1[j,3], ncols = dimen.canaveral.behavior1[j,4], format = "matrix")
        get.serie <- formationPatternsTexturesHilbert(img, n[a], tal[b])
        probability[j, 1:factorial(n[a])] <- Bandt.Pompe(get.serie, n[a], dim(get.serie)[1])
      }
      #Cape Canaveral - behavior 2
      sar_data <- raster(paste("../../Data/", "cape", "/HHHH", ".grd", sep = ""))
      for(j in c(1:ns.canaveral.behavior2)){
        img <- getValuesBlock(sar_data, row = dimen.canaveral.behavior2[j,1], nrows = dimen.canaveral.behavior2[j,2], col = dimen.canaveral.behavior2[j,3], ncols = dimen.canaveral.behavior2[j,4], format = "matrix")
        get.serie <- formationPatternsTexturesHilbert(img, n[a], tal[b])
        probability[(ns.canaveral.behavior1) + j, 1:factorial(n[a])] <- dataDriven(get.serie, n[a], tal[b])
      }
      
      #Diferent regions 
      probability[1:40, factorial(n[a]) + 1] = 1
      probability[41:80, factorial(n[a]) + 1] = 2
      
      #Diferent places
      probability[, factorial(n[a]) + 2] = 1      #Canaveral
      
      Entropy.Complexity <- Entropy.Complexity.Points(probability, dim(probability)[1], n[a])
      
      XMIN = min(Entropy.Complexity$H, XMIN) 
      YMAX = max(Entropy.Complexity$C, YMAX)
      
      plots[[i]] <- HCPlane.Classes(probability, n[a], Entropy.Complexity) + theme(axis.title.x=element_blank(),
                                                               axis.title.y=element_blank(), legend.position = "none")
    }
    
    for(i in 1:(length(n)*length(tal))){
      plots[[i]] = plots[[i]] + xlim(limits=c(XMIN, XMAX)) + ylim(limits=c(YMIN, YMAX)) 
    }
    
    ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]],
              plots[[6]], plots[[7]], plots[[8]], plots[[9]], plots[[10]],
              plots[[11]], plots[[12]], plots[[13]], plots[[14]], plots[[15]],
              plots[[16]], plots[[17]], plots[[18]], plots[[19]], plots[[20]],
              #plots[[21]], plots[[22]], plots[[23]], plots[[24]], plots[[25]],
              ncol=5, nrow=4, common.legend = TRUE, legend = "right") + 
      ggtitle(expression(italic("Cape Canaveral Region - Hilbert Curve Analysis"))) +
      xlab(expression(italic(H))) + ylab(expression(italic(C))) + labs(colour=expression(italic(Regions))) +
      theme_igray() + theme(text=element_text(size=14, family="Times New Roman"), axis.text.x=element_blank(), axis.text.y=element_blank(),plot.title = element_text(hjust=0.5)) + 
      guides(colour = guide_legend(override.aes = list(size=3)))
  }
  else if(Analysis == 4){ #Analysis of Munich region
    
    XMAX = 1
    XMIN = 1
    YMAX = 0
    YMIN = 0
    
    for(i in 1:(length(n)*length(tal))){
      cat("- Plane: ", i, "de 20 ", "\n")
      if(i%%5 == 1){
        a = a + 1
        b = 0
      }
      b = b + 1
      symbols <- definePatterns(n[a])
      probability <- matrix(nrow = (ns.munich), ncol = factorial(n[a]) + 2)
      #Munich
      sar_data <- raster(paste("../../Data/", "munich", "/HHHH", ".grd", sep = ""))
      for(j in c(1:ns.munich)){
        img <- getValuesBlock(sar_data, row = dimen.munich[j,1], nrows = dimen.munich[j,2], col = dimen.munich[j,3], ncols = dimen.munich[j,4], format = "matrix")
        get.serie <- formationPatternsTexturesHilbert(img, n[a], tal[b])
        probability[j, 1:factorial(n[a])] <- Bandt.Pompe(get.serie, n[a], dim(get.serie)[1])
      }
      
      #Diferent regions 
      probability[, factorial(n[a]) + 1] = 1
      
      #Diferent places
      probability[, factorial(n[a]) + 2] = 1      #Munich
      
      Entropy.Complexity <- Entropy.Complexity.Points(probability, dim(probability)[1], n[a])
      
      XMIN = min(Entropy.Complexity$H, XMIN) 
      YMAX = max(Entropy.Complexity$C, YMAX)
      
      plots[[i]] <- HCPlane.Classes(probability, n[a], Entropy.Complexity) + theme(axis.title.x=element_blank(),
                                                               axis.title.y=element_blank(), legend.position = "none")
    }
    
    for(i in 1:(length(n)*length(tal))){
      plots[[i]] = plots[[i]] + xlim(limits=c(XMIN, XMAX)) + ylim(limits=c(YMIN, YMAX)) 
    }
    
    ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]],
              plots[[6]], plots[[7]], plots[[8]], plots[[9]], plots[[10]],
              plots[[11]], plots[[12]], plots[[13]], plots[[14]], plots[[15]],
              plots[[16]], plots[[17]], plots[[18]], plots[[19]], plots[[20]],
              #plots[[21]], plots[[22]], plots[[23]], plots[[24]], plots[[25]],
              ncol=5, nrow=4, common.legend = TRUE, legend = "right") + 
      ggtitle(expression(italic("Munich Region - Hilbert Curve Analysis"))) +
      xlab(expression(italic(H))) + ylab(expression(italic(C))) + labs(colour=expression(italic(Regions))) +
      theme_igray() + theme(text=element_text(size=14, family="Times New Roman"), axis.text.x=element_blank(), axis.text.y=element_blank(),plot.title = element_text(hjust=0.5)) + 
      guides(colour = guide_legend(override.aes = list(size=3)))
  }
}

SAR.Hilbert.Analysis(1)
SAR.Hilbert.Analysis(2)
SAR.Hilbert.Analysis(3)
SAR.Hilbert.Analysis(4)