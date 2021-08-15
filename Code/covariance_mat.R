library(geoR)
library(gstat)


# geoR --------------------------------------------------------------------
??geoR
cov.spatial(W_scalar, cov.pars=cov(W_scalar)) #what should cov.pars be?

?varcov.spatial
varcov.spatial(dists.lowertri = W_scalar[lower.tri(W_scalar)],
               cov.pars = cov(W_scalar)) 
  #longer object length is not a multiple of shorter object length
  #still a matrix output though

cov_params <- matrix(c(2, 3, 4,
                       0.3, 0.3, 0.3), byrow=FALSE, ncol=2)
varcov.spatial(dists.lowertri = W_scalar[lower.tri(W_scalar)],
               cov.pars = cov_params) 

?globalvar

?varcovBGCCM #bivariate Gaussian common component geostatistical model
varcovBGCCM(W_scalar)
varcovBGCCM(W_scalar[lower.tri(W_scalar)])
  #Error: $ operator is invalid for atomic vectors

?variofit
?variog

stations_gd <- as.geodata(stations_sub)
vario <- variog(stations_gd, data = matrix(c(log(stations_sub$scale), stations_sub$shape), byrow=FALSE, ncol=2) )
summary(vario)
vario

variofit(vario) #object in vario$v is a matrix. This function works for only 1 empirical variogram at once

######
?likfitBGCCM
stations_gd <- as.geodata(stations_sub)
stations_gd$scale <- log(stations_sub$scale)

stations_gd2 <- as.geodata(stations_sub)
stations_gd$shape <- stations_sub$shape

bg_model <- likfitBGCCM(stations_gd, stations_gd2)
summary(bg_model)
bg_model


# gstat -------------------------------------------------------------------
??gstat
