#Below is the code to generate a TV-AR(1) and TV-GARCH(1, 1) process for a Gamma model via standardised residuals
#Given Caution when changing the functions below as this can break the time-series

#################################################################################
#################################################################################
x <- c(0:1000);

#Below are the functions for the simulation study in section 5.1
  #Change the functions to replicate the study for ARCH1 and GARCH1 processes
f.AR1 <- 0.45*sin(x*0.0125) + (0.0005*x); plot(f.AR1, type = "l");
#f.AR1 <- rep(0, length(x)); plot(f.AR1, type = "l");
#f.ARCH1 <- 0.45*sin(x*0.0125) + (0.0005*x); plot(f.ARCH1, type = "l");
f.ARCH1 <- rep(0, length(x)); plot(f.ARCH1, type = "l");
#f.GARCH1 <- 0.45*sin(x*0.0125) + (0.0005*x); plot(f.GARCH1, type = "l");
f.GARCH1 <- rep(0, length(x)); plot(f.GARCH1, type = "l");


#Below are the functions for the simulation study in section 5.2
#f.AR1 <- sin((x + 10)/75)*(50/(x + 100)); plot(f.AR1, type = "l")
#f.ARCH1 <- (exp(-((x - 500)/200)^10))*0.5; plot(f.ARCH1, type = "l");
#f.GARCH1 <- 0.0000000027*((x - 350)^3) + 0.00005*x; plot(f.GARCH1, type = "l")

#################################################################################
#################################################################################

AR1 <- func_AR1[1]; intercept <- 1; y.1 <- log(2.5); var.GARCH1 <- 0.5;
ARCH1 <- func_ARCH1[1]; GARCH1 <- func_GARCH1[1]; var.intercept <- 0.5;

Generate.TVARGARCH <- function(start_value, start_value_var, intercept, var_intercept,
                               func_AR1, func_ARCH1, func_GARCH1){
  alpha1.i <- (exp(start_value))^2/start_value_var; beta1.i = exp(start_value)/start_value_var
  y1.i <- rgamma(1, shape = alpha1.i, rate = beta1.i)
  
  y1.GARCH1 <- c(y1.i); var.GARCH.save <- c(start_value_var)
  for(i in c(1:1000)){
    y.AR1.i <- (func_AR1[i + 1]*(log(y1.GARCH1[i])))
    resi.i <- (y1.GARCH1[i] - exp(c(intercept + y.AR1.i)))
    var.ARCH1.i <- func_ARCH1[i + 1]*((resi.i/sqrt(var.GARCH.save[i]))^2)
    var.GARCH1.i <- func_GARCH1[i + 1]*log(var.GARCH.save[i])
  
    var.i <- exp(log(var_intercept) + var.ARCH1.i + var.GARCH1.i)
    var.GARCH.save <- c(var.GARCH.save, var.i)
    
    alpha.i <- exp(c(intercept + y.AR1.i))^2/var.i; beta.i = exp(c(intercept + y.AR1.i))/var.i
    y1.GARCH1.i <- rgamma(1, shape = alpha.i, rate = beta.i)
    y1.GARCH1 <- c(y1.GARCH1, y1.GARCH1.i)
  }
  return(y1.GARCH1)
}

y1.GARCH <- Generate.TVARGARCH(start_value = log(2.5), start_value_var = 0.5,
                               intercept = 1, var_intercept = 0.5,
                               func_AR1 = f.AR1, func_ARCH1 = f.ARCH1, func_GARCH1 = f.GARCH1);

plot(y1.GARCH, type = "l");

#################################################################################
#################################################################################