## Playing around with functions for CAR modeling, also has weight matrix normalization

source('./Test/functionV2.R') # read in updated functions, this houses the hausMat function to get 
                              # a matrix of extended Hausdorff distances
library(CARBayes)
#library(care)
library(spatialreg)
library(sp)
library(spdep)

# Reading in data ---------------------------------------------------------
regions <- rgdal::readOGR('~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Test/watershed_region_updated')
proj4string(regions) <- sp::CRS("+init=epsg:2278")
regions <- spTransform(regions, CRS("+init=epsg:2278"))

regions_alt <- readOGR('~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Test/watershed_region_alt')
regions_alt <- spTransform(regions_alt, CRS("+init=epsg:2278"))

regions_ts <- read.csv("~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/max_precip_regions.csv")

ws_reg_avg <- readRDS("~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/ws_reg_avg.rds")
ws_reg_avg_alt <- readRDS("~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/ws_reg_alt_avg.rds")

## function to normalize the distance matrix from hausMat()
normalizeMatrix <- function(mat, method="row") {
  if (method == "row") {
    #row standardized, the weights of a row should add up to 1
    for (i in 1:nrow(mat)) {
      mat[i,] <- mat[i,] / sum(mat[i,])
    }
    return(mat)
  }
  else if (method == "scalar") {
    alpha <- 1 / max(mat)
    mat <- alpha * mat
    return(mat)
  }
  else if (method == "inverse") {
    for (i in 1:nrow(mat)) {
      for (j in 1:ncol(mat)) {
        if (i != j) {
          mat[i, j] <- mat[i, j]^(-1)
        }
      }
    }
    return(mat)
    # for (j in 1:ncol(mat)) {
    #   full_sum <- sum(mat[, j]^(-1)) #will subtract the current value from this
    #   for (i in 1:nrow(mat)) {
    #     mat[i, j] <- mat[i, j] / (full_sum - (mat[i, j]^(-1)))
    #   }
    # }
    # return(mat)
  }
  else if (method == "exp") {
    for (i in 1:nrow(mat)) {
      for (j in 1:ncol(mat)) {
        if (i != j) {
          mat[i, j] <- exp(-0.1*mat[i, j])
        }
      }
    }
    return(mat)
    # for (j in 1:ncol(mat)) {
    #   full_sum <- sum(exp(-1 * mat[, j])) #will subtract the current value from this
    #   for (i in 1:nrow(mat)) {
    #     mat[i, j] <- exp(-1 * mat[i, j]) / (full_sum - exp(-1 * mat[i, j]))
    #   }
    # }
    # return(mat)
  }
  else {
    stop(paste("Method of normalization does not match any of the options."))
  }
}


# CARBayes ----------------------------------------------------------------
## modeling with the original hydrologic region setup
# getting weight matrix
distMat <- hausMat(regions, 0.5)
distMat

# W_row <- normalizeMatrix(distMat, "row") 
# isSymmetric.matrix(W_row) #not symmetric
# 
# W_scalar <- normalizeMatrix(distMat, "scalar")
# isSymmetric.matrix(W_scalar) #symmetric

#take the inverse first, and then scalar/row normalize
(W <- normalizeMatrix(distMat, "inverse"))
isSymmetric.matrix(W) #symmetric
(W_scalar <- normalizeMatrix(W, "scalar"))
isSymmetric.matrix(W_scalar) #symmetric

(W_row <- normalizeMatrix(W, "row"))
isSymmetric.matrix(W_row) #not symmetric

# (W <- normalizeMatrix(distMat, "exp")) #all 0
# isSymmetric.matrix(W) #symmetric

ws_reg_avg

# using CARBayes functions: NONE OF THESE WORK?
?MVS.CARleroux
# MVS.CARleroux(region_1 ~ 1, family="gaussian", data = regions_ts[-1, 2:5], W=W_scalar) 
  #the formula inputted contains an error, e.g the variables may be different lengths
MVS.CARleroux(ws_reg_avg ~ 1, family="gaussian", W=W_scalar, burnin=50, n.sample=100)
  #system is computationally singular: reciprocal condition number = 1.11679e-17

#example had "Turn the data and trials into matrices where each row is an area"
MVS.CARleroux(t(ws_reg_avg)~ 1, family="gaussian", W=W_scalar, burnin=20, n.sample=50)
MVS.CARleroux(scale_avg~1, data=t(ws_reg_avg), family="gaussian", W=W_scalar, burnin=20, n.sample=50)


?S.CARbym 
# S.CARbym(formula = region1~1, family = "binomial", data = regions_ts, trials = 10, W = W_scalar)
S.CARbym(t(ws_reg_avg)~1, family="binomial", W=W_scalar, trials=30) #response variable needs to be integer values
  
?S.CARdissimilarity
S.CARdissimilarity(t(ws_reg_avg)~1, family="gaussian", W=W_scalar, Z=) #need dissimilarity matrix

?S.CARleroux
S.CARleroux(t(ws_reg_avg)~1, family="gaussian", W=W_scalar, burnin=20, n.sample=50)
  #t(X.standardised) %*% beta.offset : non-conformable arguments, NAs produced

?S.CARlocalised #requires integer response variable
# S.CARlocalised(formula = region1~1, family = "binomial", data = regions_ts, trials = 10, W = W_scalar) 
# S.CARlocalised(formula = region1~1, family = "poisson", data = regions_ts, W = W_scalar) 


?S.CARmultilevel
S.CARmultilevel(t(ws_reg_avg)~1, family="gaussian", W=W_scalar, ind.area = c(1,2,3), burnin=20, n.sample=50)
  #t(X.standardised) %*% beta.offset : non-conformable arguments, NAs produced

S.CARmultilevel(scale_avg~1, data=t(ws_reg_avg),family="gaussian",
                W=W_scalar, ind.area = c(1,2,3), burnin=20, n.sample=50)




# care package ------------------------------------------------------------
car.models(ws_reg_avg) #this doesn't work


# CARBayesST (Spatio-Temporal) --------------------------------------------
library(CARBayesST)
ST.CARadaptive(t(ws_reg_avg)~1, family="gaussian", W=W_scalar, burnin=20, n.sample=50)
ST.CARadaptive(regions_ts[3:5] ~1, family="gaussian", W=W_scalar, burnin=20, n.sample=50)


# Short Course ------------------------------------------------------------
# car.ml <- spatialreg::spautolm(region1 ~1, data=regions_ts[3:5], family="CAR", method="eigen",
#                                listw=mat2listw(W_scalar), na.action=na.exclude, zero.policy=TRUE)
# 
# subset1 <- regions_ts[15706:44332, 3:5]
# subset1 <- regions_ts[17000:44332, 3:5]
# 
# car.ml <- spatialreg::spautolm(region1 ~1, data=subset1, family="CAR", method="eigen",
#                                listw=mat2listw(W_scalar), na.action=na.exclude, zero.policy=TRUE)

car.ml <- spatialreg::spautolm(scale_avg ~1, data=data.frame(t(ws_reg_avg)), family="CAR",
                               listw=mat2listw(W_scalar))
summary(car.ml)$Coef

spatialreg::spautolm(shape_avg ~1, data=data.frame(t(ws_reg_avg)), family="CAR",
                     listw=mat2listw(W_scalar))

spatialreg::spautolm(rate_avg ~1, data=data.frame(t(ws_reg_avg)), family="CAR",
                     listw=mat2listw(W_scalar))

# Repeat for alternative region setup -------------------------------------
## get distance matrix
distMatAlt <- hausMat(regions_alt, 0.5)
distMatAlt
distMatAlt - distMat #how different are the two matrices?

## get weight marix: normalize with inverse/exp first, and then scalar/row
# (W_alt <- normalizeMatrix(distMatAlt, "exp")) #all 0
# isSymmetric.matrix(W_alt) #symmetric

(W_alt <- normalizeMatrix(distMatAlt, "inverse"))
isSymmetric.matrix(W_alt) #symmetric
(W_scalar_alt <- normalizeMatrix(W_alt, "scalar"))
isSymmetric.matrix(W_scalar_alt) #symmetric

(W_row_alt <- normalizeMatrix(W_alt, "row"))
isSymmetric.matrix(W_row_alt) #not symmetric

## CAR modeling
spatialreg::spautolm(scale_avg ~1, data=data.frame(t(ws_reg_avg_alt)), family="CAR",
                     listw=mat2listw(W_scalar_alt))

