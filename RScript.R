#Below is the code needed to run all simulation studies
#To change the TV-AR(1) or TV-GARCH(1, 1) functions refer back to "SimulateData.R"

setwd("./StanFiles") #Set your Working Directory
#Download these R packages intall.packages("rstan"); install.packages("loo");
library(rstan); library(loo);

#fit.1 is the model with standard Horseshoe prior for TV-AR(1) spline basis coefficients
#fit.1.GARCH is the model with standard Horseshoe prior for TV-GARCH(1, 1) spline basis coefficients

#fit.1.Matrix is the model with multivariate Horseshoe prior for TV-AR(1) spline basis coefficients
#fit.1.Matrix.GARCH is the model with multivariate Horseshoe prior for TV-GARCH(1, 1) spline basis coefficients

#fit.1.Wishart is the model with a Inverse-Wishart prior on the covariance matrix of a zero-centred multivariate prior for TV-AR(1) spline basis coefficients
#fit.1.Matrix.GARCH is the model with Inverse-Wishart prior on the covariance matrix of a zero-centred multivariate prior for TV-GARCH(1, 1) spline basis coefficients

#################################################################################
#################################################################################
#Below is the code to run a model with standard Horseshoe priors

#Set-up for running RStan
mat.1 <- model.matrix(y1.GARCH[2:1001] ~ 1)
x1 <- x[2:1001]
data.1 <- list(N = 1000, ytm1 = y1.GARCH[1:1000], y = y1.GARCH[2:1001],
               KTVAR = 15, KVar = 15,
               modmat = mat.1, B_TVAR = t(phi.mat.func(x = x1, M = 15, bandwidth = 0.000025)))

#RStan Model
fit.1 <- stan("GammaTVAR1.stan",
              data = data.1,
              pars = c("beta",
                       "beta_TVAR",
                       "tau0",
                       "tau",
                       "log_lik",
                       "y_rep"),
              iter = 1500, verbose = FALSE, chains = 4, warmup = 500,
              control = list(adapt_delta = 0.99),
              init = "0", cores = 2)

print(fit.1)

#LOOIC and WAIC Estimates
fit.1.LL <- extract_log_lik(fit.1)
fit.1.looic <- loo(fit.1.LL)
fit.1.waic <- waic(fit.1.LL)

########################################
########################################
#Below are model predictions

pred <- exp(c(t(as.matrix(fit.1)[,1])) +
              t(log(y1.GARCH[1:1000])*
                  t(phi.mat.func(x = x1, M = 15, bandwidth = 0.000025))%*%t(as.matrix(fit.1)[,c(2:16)])))

pred.quant <- c()
for(j in c(1:1000)){
  pred.quant <- rbind(pred.quant, quantile(pred[,j], c(0.10, 0.5, 0.90)))
}

plot(y1.GARCH[2:1001] ~ x1, type = "l")
points(pred.quant[,2] ~ x1, type = "l", col = "red")
points(pred.quant[,1] ~ x1, type = "l", col = "red", lty = 2)
points(pred.quant[,3] ~ x1, type = "l", col = "red", lty = 2)

####################
####################

pred.TVAR <- t(t(phi.mat.func(x = x1, M = 15, 0.000025))%*%t(as.matrix(fit.1)[,c(2:16)]))

pred.TVAR.quant <- c()
for(j in c(1:1000)){
  pred.TVAR.quant <- rbind(pred.TVAR.quant, quantile(pred.TVAR[,j], c(0.1, 0.5, 0.9)))
}

plot((f.AR1) ~ x, type = "l", lty = 3, ylim = c(-0.1, 1))
points((pred.TVAR.quant[,2]) ~ x1, type = "l", col = "red")
points((pred.TVAR.quant[,1]) ~ x1, type = "l", lty = 2, col = "red")
points((pred.TVAR.quant[,3]) ~ x1, type = "l", lty = 2, col = "red")

#################################################################################
#################################################################################
#Set-up for the TV-GARCH(1, 1) process

x1.GARCH <- x1[2:1000]
fit.1.pred <- pred.quant[,2]
fit.1.resi <- c((y1.GARCH[2:1001]) - (fit.1.pred))
fit.1.sigma <- mean(as.matrix(fit.1)[,17])

data.1.GARCHStan <- list(N = 999, y = y1.GARCH[3:1001], pred = fit.1.pred[2:1000],
                         resi2 = fit.1.resi[2:1000],
                         KARCH = 15, KGARCH = 15,
                         B_ARCH = t(phi.mat.func(x = x1.GARCH, M = 15, bandwidth = 0.000025)),
                         B_GARCH = t(phi.mat.func(x = x1.GARCH, M = 15, bandwidth = 0.000025)))

#RStan Model
fit.1.GARCH <- stan("GammaTVGARCH.stan",
                     data = data.1.GARCHStan,
                     pars = c("tau_0",
                              "tau_ARCH",
                              "tau_GARCH",
                              "log_lik",
                              "y_rep"),
                    iter = 1500, verbose = FALSE, chains = 4, warmup = 500,
                     control = list(adapt_delta = 0.99),
                     init = "0", cores = 2)

print(fit.1.GARCH)

#LOOIC and WAIC Estimates
fit.1.GARCH.LL <- extract_log_lik(fit.1.GARCH)
fit.1.GARCH.looic <- loo(fit.1.GARCH.LL)
fit.1.GARCH.waic <- waic(fit.1.GARCH.LL)

####################
####################
#Predictions of the ARCH and GARCH Process

pred.TVARCH <- t(t(phi.mat.func(x = x1.GARCH, M = 15, bandwidth = 0.000025))%*%t(as.matrix(fit.1.GARCH)[,c(2:16)]))

pred.TVARCH.quant <- c()
for(j in c(1:999)){
  pred.TVARCH.quant <- rbind(pred.TVARCH.quant, quantile(pred.TVARCH[,j], c(0.1, 0.5, 0.9)))
}

plot(f.ARCH1, type = "l", ylim = c(-1, 1))
points(c((pred.TVARCH.quant[,2])), type = "l", col = "red")
points(c((pred.TVARCH.quant[,1])), type = "l", col = "red", lty = 2)
points(c((pred.TVARCH.quant[,3])), type = "l", col = "red", lty = 2)

####################
####################

pred.TVGARCH <- t(t(phi.mat.func(x = x1.GARCH, M = 15, bandwidth = 0.000025))%*%t(as.matrix(fit.1.GARCH)[,c(17:31)]))

pred.TVGARCH.quant <- c()
for(j in c(1:999)){
  pred.TVGARCH.quant <- rbind(pred.TVGARCH.quant, quantile(pred.TVGARCH[,j], c(0.1, 0.5, 0.9)))
}

plot(pred.TVGARCH[3,], type = "l")

plot(f.GARCH1, type = "l", ylim = c(-1, 1))
points((pred.TVGARCH.quant[,2]), type = "l", col = "red")
points((pred.TVGARCH.quant[,1]), type = "l", col = "red", lty = 2)
points((pred.TVGARCH.quant[,3]), type = "l", col = "red", lty = 2)

#################################################################################
#################################################################################
#################################################################################
#################################################################################
#Below is the code to run a model with multivariate Horseshoe priors

#Set-up for running RStan
mat.1 <- model.matrix(y1.GARCH[2:1001] ~ 1) #x[3:52])

x1 <- x[2:1001]
data.1.Matrix <- list(N = 1000, ytm1 = y1.GARCH[1:1000], y = y1.GARCH[2:1001],
                      KTVAR = 15, KVar = 15, vec_null = rep(0, 15),
                      modmat = mat.1, B_TVAR = t(phi.mat.func(x = x1, M = 15, bandwidth = 0.000025)),
                      B_var = t(phi.mat.func(x = x1, M = 15, bandwidth = 0.000025)))

#RStan Model
fit.1.Matrix <- stan("GammaTVAR1Matrix.stan",
                     data = data.1.Matrix,
                     pars = c("beta",
                              "beta_TVAR",
                              "tau0",
                              "eta_beta_TVAR",
                              "log_lik",
                              "y_rep"),
                     iter = 1500, verbose = FALSE, chains = 4, warmup = 500,
                     control = list(adapt_delta = 0.99),
                     init = "0", cores = 2, thin = 1)

print(fit.1.Matrix)

#LOOIC and WAIC estimates
fit.1.Matrix.LL <- extract_log_lik(fit.1.Matrix)
fit.1.Matrix.Looic <- loo(fit.1.Matrix.LL)
fit.1.Matrix.waic <- waic(fit.1.Matrix.LL)

########################################
########################################
#Prediction of the TV-AR(1) process

pred.Matrix <- exp(c(t(as.matrix(fit.1.Matrix)[,1])) +
                     t(log(y1.GARCH[1:1000])*
                         t(phi.mat.func(x = x1, M = 15, bandwidth = 0.000025))%*%t(as.matrix(fit.1.Matrix)[,c(2:16)])))

pred.Matrix.quant <- c()
for(j in c(1:1000)){
  pred.Matrix.quant <- rbind(pred.Matrix.quant, quantile(pred.Matrix[,j], c(0.10, 0.5, 0.90)))
}

plot(y1.GARCH ~ x, type = "l")
points(pred.Matrix.quant[,2] ~ x1, type = "l", col = "blue")
points(pred.Matrix.quant[,1] ~ x1, type = "l", col = "blue", lty = 2)
points(pred.Matrix.quant[,3] ~ x1, type = "l", col = "blue", lty = 2)

####################
####################

pred.Matrix.TVAR <- t(t(phi.mat.func(x = x1, M = 15, bandwidth = 0.000025))%*%t(as.matrix(fit.1.Matrix)[,c(2:16)]))

pred.Matrix.TVAR.quant <- c()
for(j in c(1:1000)){
  pred.Matrix.TVAR.quant <- rbind(pred.Matrix.TVAR.quant, quantile(pred.Matrix.TVAR[,j], c(0.10, 0.5, 0.90)))
}

plot((f.AR1) ~ x, type = "l", ylim = c(-0.1, 1))
points((pred.Matrix.TVAR.quant[,2]) ~ x1, type = "l", col = "blue")
points((pred.Matrix.TVAR.quant[,1]) ~ x1, type = "l", lty = 2, col = "blue")
points((pred.Matrix.TVAR.quant[,3]) ~ x1, type = "l", lty = 2, col = "blue")

#################################################################################
#################################################################################
#Set-up for the TV-GARCH(1, 1) process

x1.GARCH <- x1[2:1000]
fit.1.Matrix.pred <- pred.Matrix.quant[,2]
fit.1.Matrix.resi <- c(y1.GARCH[2:1001] - fit.1.Matrix.pred)
fit.1.Matrix.sigma <- mean(as.matrix(fit.1.Matrix)[,17])

data.1.GARCHStan.matrix <- list(N = 999, y = y1.GARCH[3:1001], pred = fit.1.Matrix.pred[2:1000],
                                resi2 = fit.1.Matrix.resi[2:1000],
                                KARCH = 15, KGARCH = 15, vec_null = rep(0, 15),
                                B_ARCH = t(phi.mat.func(x = x1.GARCH, M = 15, bandwidth = 0.000025)),
                                B_GARCH = t(phi.mat.func(x = x1.GARCH, M = 15, bandwidth = 0.000025)))

#RStan model
fit.1.Matrix.GARCH <- stan("GammaTVGARCHMatrix.stan",
                            data = data.1.GARCHStan.matrix,
                            pars = c("tau_0",
                                     "tau_ARCH",
                                     "tau_GARCH",
                                     "log_lik",
                                     "y_rep"),
                           iter = 1500, verbose = FALSE, chains = 4, warmup = 500,
                            control = list(adapt_delta = 0.99),
                            init = "0", cores = 2)

print(fit.1.Matrix.GARCH)

#LOOIC and WAIC estimates
fit.1.Matrix.GARCH.LL <- extract_log_lik(fit.1.Matrix.GARCH)
fit.1.Matrix.GARCH.looic <- loo(fit.1.Matrix.GARCH.LL)
fit.1.Matrix.GARCH.waic <- waic(fit.1.Matrix.GARCH.LL)

####################
####################
#Estimates for the TV-GARCH(1, 1) process

pred.Matrix.TVARCH <- t(t(phi.mat.func(x = x1.GARCH, M = 15, bandwidth = 0.000025))%*%t(as.matrix(fit.1.Matrix.GARCH)[,c(2:16)]))

pred.Matrix.TVARCH.quant <- c()
for(j in c(1:999)){
  pred.Matrix.TVARCH.quant <- rbind(pred.Matrix.TVARCH.quant, quantile(pred.Matrix.TVARCH[,j], c(0.10, 0.5, 0.90)))
}

plot(f.ARCH1, type = "l", ylim = c(-1, 1))
points(c((pred.Matrix.TVARCH.quant[,2])), col = "blue", type = "l")
points(c((pred.Matrix.TVARCH.quant[,1])), col = "blue", type = "l", lty = 2)
points(c((pred.Matrix.TVARCH.quant[,3])), col = "blue", type = "l", lty = 2)

####################
####################

pred.Matrix.TVGARCH <- t(t(phi.mat.func(x = x1.GARCH, M = 15, bandwidth = 0.000025))%*%t(as.matrix(fit.1.Matrix.GARCH)[,c(17:31)]))

pred.Matrix.TVGARCH.quant <- c()
for(j in c(1:999)){
  pred.Matrix.TVGARCH.quant <- rbind(pred.Matrix.TVGARCH.quant, quantile(pred.Matrix.TVGARCH[,j], c(0.10, 0.5, 0.90)))
}

plot(f.GARCH1, type = "l", ylim = c(-1, 1))
points(c((pred.Matrix.TVGARCH.quant[,2])), col = "blue", type = "l")
points(c((pred.Matrix.TVGARCH.quant[,1])), col = "blue", type = "l", lty = 2)
points(c((pred.Matrix.TVGARCH.quant[,3])), col = "blue", type = "l", lty = 2)

#################################################################################
#################################################################################
#Below is the code to run the model with a Inverse-Wishart prior on the covariance matrix of a zero-centred multivariate prior

#Set-up for running RStan
mat.1 <- model.matrix(y1.GARCH[2:1001] ~ 1)

x1 <- x[2:1001]
data.1.Wishart <- list(N = 1000, ytm1 = y1.GARCH[1:1000], y = y1.GARCH[2:1001],
                       KTVAR = 15, KVar = 15, vec_null = rep(0, 15),
                       modmat = mat.1, B_TVAR = t(phi.mat.func(x = x1, M = 15, bandwidth = 0.000025)),
                       B_var = t(phi.mat.func(x = x1, M = 15, bandwidth = 0.000025)),
                       Sig_Idnt = diag(15))

#RStan Model
fit.1.Wishart <- stan("GammaTVAR1Wishart.stan",
                      data = data.1.Wishart,
                      pars = c("beta",
                               "beta_TVAR",
                               "tau0",
                               "tau",
                               "log_lik",
                               "y_rep"),
                      iter = 1500, verbose = FALSE, chains = 4, warmup = 500,
                      control = list(adapt_delta = 0.99),
                      init = 0, cores = 2)

print(fit.1.Wishart)

fit.1.Wishart.LL <- extract_log_lik(fit.1.Wishart)
fit.1.Wishart.looic <- loo(fit.1.Wishart.LL)
fit.1.Wishart.waic <- waic(fit.1.Wishart.LL)

########################################
########################################
#Estimates of the AR(1) process

pred.Wishart <- exp(c(t(as.matrix(fit.1.Wishart)[,1])) +
                      t(log(y1.GARCH[1:1000])*
                          t(phi.mat.func(x = x1, M = 15, bandwidth = 0.000025))%*%t(as.matrix(fit.1.Wishart)[,c(2:16)])))

pred.Wishart.quant <- c()
for(j in c(1:1000)){
  pred.Wishart.quant <- rbind(pred.Wishart.quant, quantile(pred.Wishart[,j], c(0.10, 0.5, 0.90)))
}

plot(y1.GARCH[2:1001] ~ x1, type = "l")
points(pred.Wishart.quant[,2] ~ x1, col = "green", type = "l")
points(pred.Wishart.quant[,1] ~ x1, col = "green", type = "l", lty = 2)
points(pred.Wishart.quant[,3] ~ x1, col = "green", type = "l", lty = 2)

####################
####################

pred.Wishart.TVAR <- t(t(phi.mat.func(x = x1, M = 15, bandwidth = 0.000025))%*%t(as.matrix(fit.1.Wishart)[,c(2:16)]))

pred.Wishart.TVAR.quant <- c()
for(j in c(1:1000)){
  pred.Wishart.TVAR.quant <- rbind(pred.Wishart.TVAR.quant, quantile(pred.Wishart.TVAR[,j], c(0.10, 0.5, 0.90)))
}

plot((f.AR1) ~ x, type = "l", ylim = c(-0.1, 1))
points((pred.Wishart.TVAR.quant[,2]) ~ x1, col = "green", type = "l")
points((pred.Wishart.TVAR.quant[,1]) ~ x1, col = "green", type = "l", lty = 2)
points((pred.Wishart.TVAR.quant[,3]) ~ x1, col = "red", type = "l", lty = 2)

#################################################################################
#################################################################################
#Set-up for the TV-GARCH(1, 1) process

x1.GARCH <- x1[2:1000]
fit.1.Wishart.pred <- pred.Wishart.quant[,2]
fit.1.Wishart.resi <- c(y1.GARCH[2:1001] - fit.1.Wishart.pred)
fit.1.Wishart.sigma <- mean(as.matrix(fit.1.Wishart)[,17])

data.1.GARCHStan.matrix <- list(N = 999, y = y1.GARCH[3:1001], pred = fit.1.Wishart.pred[2:1000],
                                resi2 = fit.1.Wishart.resi[2:1000],
                                KARCH = 15, KGARCH = 15, vec_null = rep(0, 15),
                                B_ARCH = t(phi.mat.func(x = x1.GARCH, M = 15, bandwidth = 0.000025)),
                                B_GARCH = t(phi.mat.func(x = x1.GARCH, M = 15, bandwidth = 0.000025)),
                                Sig_Idnt = diag(15))

#RStan model
fit.1.Wishart.GARCH <- stan("GammaTVGARCHWishart.stan",
                            data = data.1.GARCHStan.matrix,
                            pars = c("tau_0",
                                     "tau_ARCH",
                                     "tau_GARCH",
                                     "log_lik"),
                            iter = 1500, verbose = FALSE, chains = 4, warmup = 500,
                            init = "0", cores = 2)

print(fit.1.Wishart.GARCH)

#LOOIC and WAIC estimates
fit.1.Wishart.GARCH.LL <- extract_log_lik(fit.1.Wishart.GARCH)
fit.1.Wishart.GARCH.looic <- loo(fit.1.Wishart.GARCH.LL)
fit.1.Wishart.GARCH.waic <- waic(fit.1.Wishart.GARCH.LL)

####################
####################
#Estimates of the TV-GARCH(1, 1) process

pred.Wishart.TVARCH <- t(t(phi.mat.func(x = x1.GARCH, M = 15, bandwidth = 0.000025))%*%t(as.matrix(fit.1.Wishart.GARCH)[,c(2:16)]))

pred.Wishart.TVARCH.quant <- c()
for(j in c(1:999)){
  pred.Wishart.TVARCH.quant <- rbind(pred.Wishart.TVARCH.quant, quantile(pred.Wishart.TVARCH[,j], c(0.10, 0.5, 0.90)))
}


plot(f.ARCH1, type = "l", ylim = c(-1, 1))
points(c((pred.Wishart.TVARCH.quant[,2])), col = "green", type = "l")
points(c((pred.Wishart.TVARCH.quant[,1])), col = "green", type = "l", lty = 2)
points(c((pred.Wishart.TVARCH.quant[,3])), col = "green", type = "l", lty = 2)

####################
####################

pred.Wishart.TVGARCH <- t(t(phi.mat.func(x = x1.GARCH, M = 15, bandwidth = 0.000025))%*%t(as.matrix(fit.1.Wishart.GARCH)[,c(17:31)]))

pred.Wishart.TVGARCH.quant <- c()
for(j in c(1:999)){
  pred.Wishart.TVGARCH.quant <- rbind(pred.Wishart.TVGARCH.quant, quantile(pred.Wishart.TVGARCH[,j], c(0.10, 0.5, 0.90)))
}

plot(f.GARCH1, type = "l", ylim = c(-1, 1), lty = 2)
points((pred.Wishart.TVGARCH.quant[,2]), col = "green", type = "l")
points((pred.Wishart.TVGARCH.quant[,1]), col = "green", type = "l")
points((pred.Wishart.TVGARCH.quant[,3]), col = "green", type = "l")

########################################
########################################
#Printing the LOOIC and WAIC estimates for each model

fit.1.looic
fit.1.Matrix.Looic
fit.1.Wishart.looic

fit.1.GARCH.looic
fit.1.Matrix.GARCH.looic
fit.1.Wishart.GARCH.looic

fit.1.waic
fit.1.Matrix.waic
fit.1.Wishart.waic

fit.1.GARCH.waic
fit.1.Matrix.GARCH.waic
fit.1.Wishart.GARCH.waic

#################################################################################
#################################################################################
