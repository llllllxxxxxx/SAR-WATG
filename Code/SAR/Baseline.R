require(igraph)

images.filter <- function(img, n, x, y, dim){
  seq1 = seq(0, n, by = 1)
  for(i in seq1){
    if((x-(n-i)) > 0 && (y+i) < dim && (y-i) > 0 && ((x-(n-i)) < (y-i) || (x-(n-i)) > (y+i)))
      img[x-(n-i), (y-i):(y+i)] = img[x-(n-i), (y-i):(y+i)] + 1
    else if((x-(n-i)) > 0 && (y+i) < dim && (y-i) <= 0 && ((x-(n-i)) > (y+i)))
      img[x-(n-i), 1:(y+i)] = img[x-(n-i), 1:(y+i)] + 1
    else if((x-(n-i)) > 0 && (y+i) > dim && (y-i) > 0 && ((x-(n-i)) < (y-i)))
      img[x-(n-i), (y-i):dim] = img[x-(n-i), (y-i):dim] +1
    
    if((x+(n-i)) < dim && (y+i) < dim && (y-i) > 0 && ((x+(n-i)) < (y-i) || (x+(n-i)) > (y+i)))
      img[x+(n-i), (y-i):(y+i)] = img[x+(n-i), (y-i):(y+i)] + 1
    else if((x+(n-i)) < dim && (y+i) < dim && (y-i) <= 0 && ((x+(n-i)) > (y+i)))
      img[x+(n-i), 1:(y+i)] = img[x+(n-i), 1:(y+i)] + 1
    else if((x+(n-i)) < dim && (y+i) > dim && (y-i) > 0 && ((x+(n-i)) < (y-i)))
      img[x+(n-i), (y-i):dim] = img[x+(n-i), (y-i):dim] + 1
  }
  img[x, y] = img[x, y] + 1
  img
}

images.into.network <- function(r, dim){
  
  texture.adjacency = matrix(0, ncol = dim, nrow = dim)
  for(i in c(1:dim)){
    for(j in c(1:dim)){
      texture.adjacency = images.filter(texture.adjacency, r, i, j, dim)
    }
  }
  texture.adjacency
}

images.into.network(2, 5)