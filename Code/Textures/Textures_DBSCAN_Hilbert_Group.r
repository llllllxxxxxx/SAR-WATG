###################################################################################
# Textures_hilbert.R 
#
# Individual analysis of the brodatz textures using Hilbert space filling curves
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

set.seed(123)
ntex = 45
n = c(2,3,4,5)
tal = c(5,2,3,5)
dim_image = 128
plots = array(list(), 4)
xpoints = c(0.0035, 0.0025, 0.0035, 0.0015)
subtitle = c("D = 2 and t = 5", "D = 3 and t = 2", "D = 4 and t = 3", "D = 5 and t = 5")

Group.Analysis.Textures.Hilbert <- function(){
  for(i in 1:4){
    probability <- matrix(nrow = ntex, ncol = factorial(n[i]))
    
    for(j in 1:ntex){
      cat("- Texture: ", j, "de 45 ", "\n")
      name.file <- paste("../Images/Brodatz/", textures.names[j], ".tiff", sep = "")
      img <- Read_tiff(name.file)
      img <- img[1:dim_image, 1:dim_image]
      
      get.serie <- formationPatternsTexturesHilbert(img, n[i], tal[i])
      probability[j, ] <- Bandt.Pompe(get.serie, n[i], dim(get.serie)[1])
    }
    
    Entropy.Complexity <- Entropy.Complexity.Points(probability, ntex, n[i])
    db <- dbscan::dbscan(Entropy.Complexity,  eps = .01, minPts = 2)
    groups <- db$cluster+1
    Entropy.Complexity <- data.frame("H" = Entropy.Complexity$H, "C" = Entropy.Complexity$C, "Groups" = groups)
    
    plots[[i]] = HCPlane.Groups(Entropy.Complexity, xpoints[i]) + xlab(subtitle[i]) + theme(axis.title.y=element_blank(), 
                                                                                            legend.position = "none")
  }
  
  ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]],
            ncol=4, nrow=1, common.legend = TRUE, legend = "right") + 
            ggtitle(expression(italic("DBSCAN applied in Textures - Hilbert Curve Analysis"))) + xlab(expression(italic(H))) + ylab(expression(italic(C))) + theme_igray() + 
            theme(text=element_text(size=14, family="Times New Roman"), axis.text.x=element_blank(), axis.text.y=element_blank(),plot.title = element_text(hjust=0.5)) + 
            guides(colour = guide_legend(override.aes = list(size=2)))
}