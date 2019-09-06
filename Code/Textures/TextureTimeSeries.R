require(tiff)

Read_tiff<-function(name.file){
  myImage<-readTIFF(name.file)
  return(myImage)
}

OneDimensionTextureHilbert <- function(img){
  hilbertcurve <- unlist(read.table("../../Data/Hilbert/HilbertCurves128.txt")) + 1
  timeseries <- img[hilbertcurve]
  timeseries
}