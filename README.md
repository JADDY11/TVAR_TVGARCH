# TVAR_TVGARCH
Read Me file for the files to replicate the study for the paper Horseshoe Priors for Time-Varying AR and GARCH processes

Author J W G Addy
Date 10/11/2023

####################################################################################
####################################################################################

“RScript.R” is the R file to rerun the analysis of the simulation study
“SimulateDat.R” contains a function to simulate a time-vary AR(1) and GARCH(1, 1) processes for a Gamma model
“SplineFunctionsR” generates the functions for B-splines with Gaussian functions
“StanFiles” is a folder that contains all Stan files
	“GammaTVAR1.stan” is a Stan file for a Gamma model with a TV-AR(1) process with Horseshoe priors
	“GammaTVAR1Matrix.stan” is a Stan file for a Gamma model with a TV-AR(1) process with Multivariate Horseshoe priors
	“GammaTVAR1Wishart.stan” is a Stan file for a Gamma model with a TV-AR(1) process with an Inverse-Wishart prior
	“GammaTVGARCH.stan” is a Stan file for a Gamma model with a TV-GARCH(1) process with Horseshoe priors
	“GammaTVGARCHMatrix.stan” is a Stan file for a Gamma model with a TV-GARCH(1) process with Multivariate Horseshoe priors
	“GammaTVGARCHWishart.stan” is a Stan file for a Gamma model with a TV-GARCH(1) process with an Inverse-Wishart prior
