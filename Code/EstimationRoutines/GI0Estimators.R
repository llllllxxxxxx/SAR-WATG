
### Parameters of GI0
params_alpha = c(-1.5, -3, -8)
params_gamma = - params_alpha - 1
params_Looks = c(1, 3, 5, 8)

### Number of samples
n = c(30, 100, 1000)

### Parameters of Monte Carlo
R = 100000
Rep = R/n

### Generate GI0 by Multiplicative Model 
randGI0 <- function(n, alpha, gamma, Looks) {
  X = rgamma(n, shape = -alpha, rate = gamma)
  Y = rgamma(n, shape = Looks, rate = Looks)
  
  return (Y/X)
}


# Maximum Likelihood Estimators (MLE) #

library(stats4)

densGI0 <- function(x, alpha, gamma, Looks) {
    return(((Looks^Looks * gamma(Looks-alpha))/ (gamma^alpha * gamma(Looks) * gamma(-alpha))) *
           (x^(Looks-1) / (gamma + Looks * x)^(Looks-alpha)))
}

LL_Brent <- function (alpha) {
    D = densGI0(Z, alpha, -alpha-1, params_Looks[4])
    -sum(log(D))
}

vector_alpha = c()
vector_error = c()
for(i in 1:Rep[2]){
  
  Z = randGI0(n[2], params_alpha[3], params_gamma[3], params_Looks[4])
  estimateBrent = mle(minuslogl = LL_Brent, start = list(alpha = -8.0), method = "Brent", lower = c(-10), upper = c(-1))
  value = coef(estimateBrent)
  list = data.frame(value)
  
  value_alpha = list$value[1]
  vector_alpha = c(vector_alpha, value_alpha)
  
  value_error = ((value_alpha - params_alpha[3])^2)
  vector_error = c(vector_error, value_error)
    
}


# Results 
estim_alpha = mean(vector_alpha)
mse = mean(vector_error)
bias = estim_alpha - params_alpha[3]

results = c(estim_alpha, bias, mse)
print(results)

######################################################################################

# Fractional Moments Estimators (Moments-1/2) #

library(rootSolve)

fun_Mom <- function (alpha) {
  return ( (1/n[3])*sum(sqrt(Z)) - sqrt((-alpha-1)/params_Looks[4])*((gamma(-alpha-1/2)*gamma(params_Looks[4]+1/2))/(gamma(-alpha)*gamma(params_Looks[4]))) )
} 

vector_alpha = c()
vector_error = c()

for (i in 1:Rep[3]){
    
  Z = randGI0(n[3], params_alpha[3], params_gamma[3], params_Looks[4])
    
  value_alpha = uniroot(fun_Mom, c(-10.0, -1.1))$root
  vector_alpha = c(vector_alpha, value_alpha)
    
  value_error = ((value_alpha - params_alpha[3])^2)
  vector_error = c(vector_error, value_error)
    
}

print(vector_alpha)

# Results 
estim_alpha <- mean(vector_alpha)
mse = mean(vector_error)
bias = estim_alpha - params_alpha[3]

results = c(estim_alpha, bias, mse)

print(results)


#Test

# print(vector_alpha)

# fun_Mom <- function (alpha) {
#   return ( (1/n[1])*sum(sqrt(Z)) - sqrt((-alpha-1)/params_Looks[4])*((gamma(-alpha-1/2)*gamma(params_Looks[4]+1/2))/(gamma(-alpha)*gamma(params_Looks[4]))) )
# }
# 
# Z = randGI0(n[1], params_alpha[1], params_gamma[1], params_Looks[4])
# x <- seq(-40.0 , -1.1 , 0.01)
# y <- fun_Mom(x)
# plot(x, y, ylim = c(-1,1),type = "l")
# abline(h=0)

######################################################################################

# Log-cumulants Estimators #

library(rootSolve)

fun_LCum <- function (alpha) {
    return ( (1/n[3])*sum(log(Z)) - log((-1-alpha)/params_Looks[4]) - digamma(params_Looks[4]) + digamma(-alpha) )
} 

vector_alpha = c()
vector_error = c()

for (i in 1:Rep[3]){
  
  Z = randGI0(n[3], params_alpha[3], params_gamma[3], params_Looks[4])
  
  value_alpha = uniroot(fun_LCum, c(-15.0, -1.1), maxiter = 100000)$root
  vector_alpha = c(vector_alpha, value_alpha)
  
  value_error = ((value_alpha - params_alpha[3])^2)
  vector_error = c(vector_error, value_error)
  
}

print(vector_alpha)

# Results 
estim_alpha <- mean(vector_alpha)
mse = mean(vector_error)
bias = estim_alpha - params_alpha[3]

results = c(estim_alpha, bias, mse)

print(results)


fun_LCum <- function (alpha) {
  return ( (1/n[2])*sum(log(Z)) - log((-1-alpha)/params_Looks[1]) - digamma(params_Looks[1]) + digamma(-alpha) )
} 

Z = randGI0(n[2], params_alpha[3], params_gamma[3], params_Looks[1])
x <- seq(-7.5 , -2.1 , 0.01)
y <- fun_LCum(x)
plot(x, y, ylim = c(-1,1),type = "l")
abline(h=0)

############################################################################

# Estimation by the Minimization of Stochastic Distances

library(cubature)

vector_alpha = c()
vector_error = c()

for (i in 1:Rep[3]){

  Z <- randGI0(n[3], params_alpha[3], params_gamma[3], params_Looks[4])
  hist <- hist(Z, breaks="FD", plot=FALSE)
  #lines(hist$mids, densGI0(hist$mids, params_alpha[1], params_gamma[1], params_Looks[1]))
  f_estimate <- densGI0(hist$mids, params_alpha[3], params_gamma[3], params_Looks[4])
  
  smaller_dTriangular = 999
  set_alpha = seq(-1, -10, by = -0.01)
  estimateMinDistStoc = -20
  
  for(alpha in set_alpha){
    
    f <- function(Z) { (((densGI0(Z, alpha, params_gamma[3], params_Looks[4])-f_estimate)^2)/(densGI0(Z, alpha, params_gamma[3], params_Looks[4])+f_estimate)) }
    
    integ <- adaptIntegrate(f, lowerLimit = c(0.0), upperLimit = c(0.2))$integral
    
    dTriangular = integ
    
    if(dTriangular < smaller_dTriangular){
      estimateMinDistStoc = alpha
      smaller_dTriangular = dTriangular
    }
    
  }
  
  vector_alpha = c(vector_alpha, estimateMinDistStoc)
  
  value_error = ((estimateMinDistStoc - params_alpha[3])^2)
  vector_error = c(vector_error, value_error)

}

print(vector_alpha)

# Results 
estim_alpha <- mean(vector_alpha)
mse = mean(vector_error)
bias = estim_alpha - params_alpha[3]

results = c(estim_alpha, bias, mse)


print(results)

####################################################################################################

