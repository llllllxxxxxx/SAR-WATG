#######################################################################################################
# PCA_Analysis.R
#
# SAR image analysis using parameter PCA: gamma distribution parameters and information theory 
# descriptors applied in BPW  
#
# Author: Eduarda Chagas
# Date : Set 2019
# Contact: eduarda.chagas@dcc.ufmg.br
########################################################################################################

############################################# Packages #################################################
require(devtools)
require(ggbiplot)
require(factoextra)
require(ggfortify)
require(plotly)

############################### gamma distribution parameters ##########################################
header = paste("../../Data/Gamma/")

file.guatemala = paste(header, "gamma.guatemala",".csv", sep = "")
gamma.guatemala = read.csv(file.guatemala)
gamma.guatemala = data.frame(gamma.guatemala$Estimate.shape, gamma.guatemala$Estimate.rate, gamma.guatemala$Sd.shape, gamma.guatemala$Sd.rate)

file.cape1 = paste(header, "gamma.cape1",".csv", sep = "")
gamma.cape1 = read.csv(file.cape1)
gamma.cape1 = data.frame(gamma.cape1$Estimate.shape, gamma.cape1$Estimate.rate, gamma.cape1$Sd.shape, gamma.cape1$Sd.rate)

file.cape2 = paste(header, "gamma.cape2",".csv", sep = "")
gamma.cape2 = read.csv(file.cape2)
gamma.cape2 = data.frame(gamma.cape2$Estimate.shape, gamma.cape2$Estimate.rate, gamma.cape2$Sd.shape, gamma.cape2$Sd.rate)

file.munich = paste(header, "gamma.munich",".csv", sep = "")
gamma.munich = read.csv(file.munich)
gamma.munich = data.frame(gamma.munich$Estimate.shape, gamma.munich$Estimate.rate, gamma.munich$Sd.shape, gamma.munich$Sd.rate)

gamma.parameters = data.frame("Estimate.shape" = numeric(160), "Estimate.rate" = numeric(160), "Sd.shape" = numeric(160), "Sd.rate" = numeric(160))
gamma.parameters[1:40, ] = gamma.guatemala
gamma.parameters[41:80, ] = gamma.cape1
gamma.parameters[81:120, ] = gamma.cape2
gamma.parameters[121:160, ] = gamma.munich

############################### Theory Information Descriptors ##########################################
header = paste("../../Data/WPE/SAR.Entropy.Complexity.Fisher")

a = b = 0
n = c(3,4,5,6) #Dimension parameter
tal = c(1,2,3,4,5) #Delay parameter
pca = res.pca = array(list(), 20)
res = array(dim = c(20, 160, 4))

for(i in 1:(length(n)*length(tal))){
  cat("- Plane: ", i, "de 20 ", "\n")
  if(i%%5 == 1){
    a = a + 1
    b = 0
  }
  b = b + 1
  
  HCF.file = paste(header, "D", n[a], "t", tal[b], ".csv", sep = "")
  Entropy.complexity.fisher = read.csv(HCF.file)
  Entropy.complexity.fisher = data.frame("H" = Entropy.complexity.fisher$H, "Fs" = Entropy.complexity.fisher$Fs, "C" = Entropy.complexity.fisher$C)
  
  regions = array(dim = 160)
  regions[1:40] = "Guatemala"
  regions[41:80] = "Cape 1"
  regions[81:120] = "Cape 2"
  regions[121:160] = "Munich"
  
  parameters = data.frame(Entropy.complexity.fisher, gamma.parameters)
  
  pca[[i]] = prcomp(parameters, center = TRUE, scale. = TRUE, rank. = 3)
  res.pca[[i]] = data.frame(pca[[i]]$x, regions)
  res[i, , ] =  matrix(unlist(res.pca[[i]]), nrow = 160, ncol = 4)
}

############################################# Plotting #################################################

plot.pca.2d <- function(res, j){
  
  shape.select <- c(17,18,19,8)
  
  # Paleta montada a partir de https://coolors.co/
  rainbow.colors <- palette(c("#494947", #DarkGreen
                              "#7494EA", #MutedDarkBlue
                              "#B14AED", #Violet
                              "#44CCFF" #BrightLightBlue
  )) 
  
  Color = rainbow.colors[res[j,,3]]
  Shape = shape.select[res[j,,3]]
  PCA <- data.frame("PCA1" = res[j,,1], "PCA2" = res[j,,2], "Color" = Color, "Shape" = Shape)
  p = qplot(xlab=expression(PCA1), ylab=expression(pCA2)) +
      theme(plot.title = element_text(hjust=0.5)) + 
      geom_point(aes(x = PCA$PCA1, y = PCA$PCA2), shape = PCA$Shape, color = PCA$Color, size = 1)
  print(p)
  
}

plot.pca.3d.individual <- function(res, j){
  p <- plot_ly(mtcars, x = res[j,,1], y = res[j,,2], z = res[j,,3], color = res[j,,4], colors = c('#BF382A', '#0C4B8E', '#B14AED', '#44CCFF')) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(title = 'PCA1'),
                        yaxis = list(title = 'PCA2'),
                        zaxis = list(title = 'PCA3')))
  print(p)
  
}

plot.pca.3d <- function(res){
  
  Sys.setenv("plotly_username"="EduardaChagas")
  Sys.setenv("plotly_api_key"="HPG76tLedL8bBCjcECRr")
  
  aa <- c("PCA_D1_T1", "PCA_D1_T2", "PCA_D1_T3", "PCA_D1_T4", "PCA_D1_T5", "PCA_D2_T1", "PCA_D2_T2", 
          "PCA_D2_T3", "PCA_D2_T4", "PCA_D2_T5", "PCA_D3_T1", "PCA_D3_T2", "PCA_D3_T3", "PCA_D3_T4", 
          "PCA_D3_T5", "PCA_D4_T1", "PCA_D4_T2", "PCA_D4_T3", "PCA_D4_T4", "PCA_D4_T5", "PCA_D5_T1", 
          "PCA_D5_T2", "PCA_D5_T3", "PCA_D5_T4", "PCA_D5_T5")
  
  for(j in 1:dim(res)[1]){
    
    cat("- Plane: ", j, "de 20 ", "\n")
    
    p <- plot_ly(mtcars, x = res[j,,1], y = res[j,,2], z = res[j,,3], color = res[j,,4], colors = c('#BF382A', '#0C4B8E', '#B14AED', '#44CCFF')) %>%
      add_markers() %>%
      layout(scene = list(xaxis = list(title = 'PCA1'),
                          yaxis = list(title = 'PCA2'),
                          zaxis = list(title = 'PCA3')))
    options(browser = 'false')
    api_create(p, filename = aa[j])
  }
}

plot.pca.3d(res)