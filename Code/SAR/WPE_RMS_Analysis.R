require(MFDFA)
require(gtools)
source("PBWeight.R")

series_generator_fk <- function(pp, y, n, k){
  Series <- vector(mode="numeric")
  filtro <- (1:n)^-(k/2)
  filtro <- filtro / sum(filtro)
  y1 <- y * filtro    
  x1 <- IFFT(y1, plan=pp)  
  Series <- c(Re(x1)) 
  Series
}

N <- 10^4
D <- c(3,4,5,6)
Tau <- c(1,2,3,4,5)

set.seed(seed = 1234567890, kind = "Mersenne-Twister")
x <- rnorm(N)
x <- x - mean(x)
pp <- planFFT(N)
y <- FFT(x, plan=pp)
fk <- series_generator_fk(pp, y, N, k = 0)

shannon <- function(ts, t, d, tal){
  h = rep(0, length(ts)/t)
  i = j = 1
  while((i + t - 1) <= length(ts)){
    bp = distribution(ts[i:(i + t - 1)], d, tal)
    h[j] = shannonNormalized(bp)
    i = (i + t)
    j = j + 1
  }
  h
}

wpe <- function(ts, t, d, tal){
  h = rep(0, length(ts)/t)
  i = j = 1
  while((i + t - 1) <= length(ts)){
    bp = WPE(ts[i:(i + t - 1)], d, tal)
    h[j] = shannonNormalized(bp)
    i = (i + t)
    j = j + 1
  }
  h
}

wpe.rsm <- function(ts, t, d, tal){
  h = rep(0, length(ts)/t)
  i = j = 1
  while((i + t - 1) <= length(ts)){
    bp = WPE.RMS(ts[i:(i + t - 1)], d, tal)
    h[j] = shannonNormalized(bp)
    i = (i + t)
    j = j + 1
  }
  h
}


ts.10 = MFsim(N = 2^10,a = 0.75)

h.ts.10 = shannon(ts.10, 2^4, 3, 1)
plot(h.ts.10, type = "l")

wpe.ts.10 = wpe(ts.10, 2^4, 3, 1)
plot(wpe.ts.10, type = "l")

wpe.rsm.ts.10 = wpe.rsm(ts.10, 2^4, 3, 1)
plot(wpe.rsm.ts.10, type = "l")



wpe.rms.ts.10 = WPE.RMS(ts.10, 3, 1)
h = shannonNormalized(wpe.rms.ts.10)
c = Ccomplexity(wpe.rms.ts.10)
HC.wpe.rms.ts.10 = data.frame(H = h, C = c) 
print(HC.wpe.rms.ts.10)
#cotas(3) + geom_point(aes(x = HC.wpe.rms.cape$H, y = HC.wpe.rms.cape$C), size = 1) + theme(plot.title = element_text(hjust=0.5)) 

wpe.ts.10 = WPE(ts.10, 3, 1)
h = shannonNormalized(wpe.ts.10)
c = Ccomplexity(wpe.ts.10)
HC.wpe.ts.10 = data.frame(H = h, C = c) 
print(HC.wpe.ts.10)
#cotas(3) + geom_point(aes(x = HC.wpe.cape$H, y = HC.wpe.cape$C), size = 1) + theme(plot.title = element_text(hjust=0.5)) 


bp = distribution(ts.10, 3, 1)
h = shannonNormalized(bp)
c = Ccomplexity(bp)
HC.ts.10 = data.frame(H = h, C = c) 
print(HC.ts.10)

ts.16 = MFsim(N = 2^16,a = 0.75)
plot(ts.16, type = "l")