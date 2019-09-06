###################################################################################
# Textures_raster1.R 
#
# Analysis of the brodatz textures using Raster-1 curves
#
# Author: Eduarda Chagas
# Date : May 2019
# Contact: eduarda-chagas@ufmg.br
####################################################################################

source("TextureTimeSeries.R")
source("Band&Pompe.R")
require(grid)
require(gridExtra)
require(ggsci)
require(extrafont)
require(ggpubr)
require(png)
require(dbscan)

textures.names = c("1.1.01", "1.1.02", "1.1.03", "1.1.04", "1.1.05",
                    "1.1.06", "1.1.07", "1.1.08", "1.1.09", "1.1.10",
                    "1.1.11", "1.1.12", "1.1.13", "1.3.01", "1.3.02",
                    "1.3.03", "1.3.04", "1.3.05", "1.3.06", "1.3.07",
                    "1.3.08", "1.3.09", "1.3.10", "1.3.11", "1.3.12",
                    "1.3.13", "1.4.01", "1.4.02", "1.4.03", "1.4.04", 
                    "1.4.05", "1.4.06", "1.4.07", "1.4.08", "1.4.09", 
                    "1.4.10", "1.4.11", "1.4.12", "1.5.01", "1.5.02", 
                    "1.5.03", "1.5.04", "1.5.05", "1.5.06", "1.5.07", 
                    "1102", "1201", "1202", "1203", "1204", "1205", "1206",
                    "1207", "1208", "1209", "1210", "1211", "1212","1213")

textures.no.contrast = c("1.1.01", "1.1.02", "1.1.03", "1.1.04", "1.1.05",
                          "1.1.06", "1.1.07", "1.1.08", "1.1.09", "1.1.10",
                          "1.1.11", "1.1.12", "1.1.13")
textures.contrast = c("1201", "1202", "1203", "1204", "1205", 
                       "1206", "1207", "1208", "1209", "1210", 
                       "1211", "1212","1213")

set.seed(123)
ntex = 45
ntex.contrast = 13
n = c(2,3,4,5,6)
tal = c(1,2,3,4,5)
a = b = 0
plots = array(list(), 25)
dim_image = 128

Textures.Raster1.Analysis <- function(Analysis){
  
  if(Analysis == 1){ #Analysis of textures characterization
    for(i in 1:(length(n)*length(tal))){
      cat("- Plane: ", i, "de 25 ", "\n")
      if(i%%length(tal) == 1){
        a = a + 1
        b = 0
      }
      b = b + 1
      probability <- matrix(nrow = ntex, ncol = factorial(n[a]))
      for(j in 1:ntex){
        cat("- Texture: ", j, "de 45 ", "\n")
        name.file <- paste("../Images/Brodatz/", textures.names[j], ".tiff", sep = "")
        img <- Read_tiff(name.file)
        img <- img[1:dim_image, 1:dim_image]
        xy <- get.number.patterns(dim(img)[1], dim(img)[2], n[a], tal[b])
        my.elements <- matrix(nrow = xy, ncol = n[a]*n[a]) 
        
        get.serie <- formationPatternsTexturesHorizontal(img, n[a], tal[b], xy)
        for(w in c(1:dim(get.serie)[3])){
          my.elements[w,] <- as.vector(t(get.serie[,,w]))
        }
        probability[j, ] <- Bandt.Pompe(as.vector(t(my.elements)), n[a], dim(get.serie)[3])
      }
      
      Entropy.Complexity <- Entropy.Complexity.Points(probability, ntex, n[a])
      Entropy.Complexity = data.frame("H" = Entropy.Complexity$H, "C" = Entropy.Complexity$C, "Texture" = 1:ntex)
      plots[[i]] <- HCPlane(Entropy.Complexity) + theme(axis.title.x=element_blank(),
                                                        axis.title.y=element_blank(), 
                                                        legend.position = "none")
    }
    
    ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]],
              plots[[6]], plots[[7]], plots[[8]], plots[[9]], plots[[10]],
              plots[[11]], plots[[12]], plots[[13]], plots[[14]], plots[[15]],
              plots[[16]], plots[[17]], plots[[18]], plots[[19]], plots[[20]],
              plots[[21]], plots[[22]], plots[[23]], plots[[24]], plots[[25]],
              ncol=5, nrow=5, common.legend = TRUE, legend = "right") + 
      ggtitle(expression(italic("Textures Raster-1 Analysis"))) + xlab(expression(italic(H))) + ylab(expression(italic(C))) + theme_igray() +
      theme(text=element_text(size=14, family="Times New Roman"), axis.text.x=element_blank(), axis.text.y=element_blank(),plot.title = element_text(hjust=0.5)) + 
      guides(colour = guide_legend(override.aes = list(size=2)))
  }
  else if(Analysis == 2){ #Analysis of textures characterization with or without contrast
    for(i in 1:(length(n)*length(tal))){
      cat("- Plane: ", i, "de 25 ", "\n")
      if(i%%length(tal) == 1){
        a = a + 1
        b = 0
      }
      b = b + 1
      probability <- matrix(nrow = (ntex.contrast*2), ncol = factorial(n[a])+2)
      for(j in 1:(ntex.contrast*2)){
        cat("- Texture: ", j, "de 26", "\n")
        if(j <= ntex.contrast){
          name.file <- paste("../Images/Brodatz/", textures.no.contrast[j], ".tiff", sep = "")
          img <- Read_tiff(name.file)
        }else{
          name.file <- paste("../Images/Brodatz/", textures.contrast[j - ntex.contrast], "png.png", sep = "")
          img <- readPNG(name.file)
        }
        img <- img[1:dim_image, 1:dim_image]
        xy <- get.number.patterns(dim(img)[1], dim(img)[2], n[a], tal[b])
        get.serie <- formationPatternsTexturesHorizontal(img, n[a], tal[b], xy)
        my.elements <- matrix(nrow = dim(get.serie)[3], ncol = n[a]*n[a]) 
        for(w in c(1:dim(get.serie)[3])){
          my.elements[w,] <- as.vector(t(get.serie[,,w]))
        }
        probability[j, 1:factorial(n[a])] <- Bandt.Pompe(as.vector(t(my.elements)), n[a], dim(get.serie)[3])
      } 
      
      #Diferent textures
      probability[, factorial(n[a])+1] = rep(1:ntex.contrast, length = (ntex.contrast*2))
      #Contrast or not
      probability[, factorial(n[a])+2] = rep(c(1,4), each = ntex.contrast)
      
      plots[[i]] <- HCPlane.Classes(probability, n[a]) + theme(axis.title.x=element_blank(),
                                                               axis.title.y=element_blank(), 
                                                               legend.position = "none")
    }
    
    ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]],
              plots[[6]], plots[[7]], plots[[8]], plots[[9]], plots[[10]],
              plots[[11]], plots[[12]], plots[[13]], plots[[14]], plots[[15]],
              plots[[16]], plots[[17]], plots[[18]], plots[[19]], plots[[20]],
              plots[[21]], plots[[22]], plots[[23]], plots[[24]], plots[[25]],
              ncol=5, nrow=5, common.legend = TRUE, legend = "right") + 
      ggtitle(expression(italic("Textures with and without contrast - Raster-1 Analysis"))) + xlab(expression(italic(H))) + ylab(expression(italic(C))) + theme_igray() +
      theme(text=element_text(size=14, family="Times New Roman"), axis.text.x=element_blank(), axis.text.y=element_blank(),plot.title = element_text(hjust=0.5)) + 
      guides(colour = guide_legend(override.aes = list(size=2)))
  }
  else{ #Analysis of Clusters in HC Plane
    for(i in 1:(length(n)*length(tal))){
      cat("- Plane: ", i, "de 25 ", "\n")
      if(i%%length(tal) == 1){
        a = a + 1
        b = 0
      }
      b = b + 1
      probability <- matrix(nrow = ntex, ncol = factorial(n[a]))
      for(j in 1:ntex){
        cat("- Texture: ", j, "de 45 ", "\n")
        name.file <- paste("../Images/Brodatz/", textures.names[j], ".tiff", sep = "")
        img <- Read_tiff(name.file)
        img <- img[1:dim_image, 1:dim_image]
        xy <- get.number.patterns(dim(img)[1], dim(img)[2], n[a], tal[b])
        get.serie <- formationPatternsTexturesHorizontal(img, n[a], tal[b], xy)
        my.elements <- matrix(nrow = dim(get.serie)[3], ncol = n[a]*n[a]) 
        for(w in c(1:dim(get.serie)[3])){
          my.elements[w,] <- as.vector(t(get.serie[,,w]))
        }
        probability[j, ] <- Bandt.Pompe(as.vector(t(my.elements)), n[a], dim(get.serie)[1])
      }
      
      Entropy.Complexity <- Entropy.Complexity.Points(probability, ntex, n[a])
      db <- dbscan::dbscan(Entropy.Complexity,  eps = .025, minPts = 2)
      groups <- db$cluster+1
      Entropy.Complexity <- data.frame("H" = Entropy.Complexity$H, "C" = Entropy.Complexity$C, "Groups" = groups)
      
      plots[[i]] <- HCPlane.Groups(Entropy.Complexity) + theme(axis.title.x=element_blank(),
                                                               axis.title.y=element_blank(), 
                                                               legend.position = "none")
    }
    
    ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]],
              plots[[6]], plots[[7]], plots[[8]], plots[[9]], plots[[10]],
              plots[[11]], plots[[12]], plots[[13]], plots[[14]], plots[[15]],
              plots[[16]], plots[[17]], plots[[18]], plots[[19]], plots[[20]],
              plots[[21]], plots[[22]], plots[[23]], plots[[24]], plots[[25]],
              ncol=5, nrow=5, common.legend = TRUE, legend = "right") + 
      ggtitle(expression(italic("DBSCAN applied in Textures - Raster-1 Analysis"))) + xlab(expression(italic(H))) + ylab(expression(italic(C))) + theme_igray() +
      theme(text=element_text(size=14, family="Times New Roman"), axis.text.x=element_blank(), axis.text.y=element_blank(),plot.title = element_text(hjust=0.5)) + 
      guides(colour = guide_legend(override.aes = list(size=2)))
  }
}
