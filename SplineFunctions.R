#This file contains all the functions needed to run the script

#################################################################################
#################################################################################
#Below is the Gauss function used to specify the spline basis function

Gauss.func <- function(mu, x, bandwidth){
  phi = exp(-bandwidth*(x - mu)^2)
  return(phi)
}

#Below is the function used to build the spline basis matrix
phi.mat.func <- function(x, M, bandwidth){
  phi.mat = c()
  for(i in c(seq(min(x), max(x), dist(range(x))/(M-1)))){
    phi.mat = rbind(phi.mat, Gauss.func(i, x, bandwidth))
  }
  return(phi.mat)
}

#################################################################################
#################################################################################