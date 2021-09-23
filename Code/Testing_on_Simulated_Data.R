#######################
# Testing the models on simulated data
#
# by Carly Fagnant
#######################
library(extRemes)
library(rgeos)
library(rgdal)
library(sp)
library(gstat)
library(dplyr)
library(spdep)
library(spatialreg)
library(foreach)

thresh = 253

### Load data
grid_regs_df <- readRDS(file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/grid_regs_df.rds")


# # Function to project data given...
# #     df: a data frame with data and coordinates
# project_data <- function(df){
#   coordinates(df)=~long+lat
#   proj4string(df) <- "+proj=longlat +datum=WGS84"
#   df <- spTransform(df, CRS("+init=epsg:2278")) # projecting data
#   is.projected(df)
#   return(df)
# }


# Function to format data to stations within the regions, given...
#   sim_data: data frame of simulated rainfall data (each column is the data for a location on the grid)
#   reg_polygons: the SpatialPolygonsDataFrame for the regions of interest
#   grid_regs_df: data frame with the coordinates and Region information of the simulated grid points (pull from saved .rds file above)
format_data_from_sim <- function(sim_data, reg_polygons, grid_regs_df){
  # Accessing Model Fits and saving the parameter values (scale, shape, rate) for each station
  ws_scale <- NULL
  ws_shape <- NULL
  ws_rate  <- NULL
  for(i in 1:ncol(sim_data)){
    fit <- extRemes::fevd(sim_data[,i], threshold=thresh, type="GP", method="MLE")
    if(is.na(fit)){
      ws_scale[i] <- ws_shape[i] <- ws_rate[i]  <- NA
    }else{
      ws_scale[i] <- fit$results$par[1]
      ws_shape[i] <- fit$results$par[2]
      ws_rate[i]  <- fit$rate
    }
  }
  grid_regs_df <- grid_regs_df %>%
    dplyr::mutate(scale = ws_scale, shape = ws_shape, rate = ws_rate) %>%
    dplyr::filter(!is.na(scale))
  
  # "Projecting" grid data... coordinates are already projected, but have to remind it of proj4string
  grid_regs <- grid_regs_df
  coordinates(grid_regs)=~x1+x2
  proj4string(grid_regs) <- "+init=epsg:2278"
  
  # already has Region information! so don't have to do intersect again
  # also no station indices, so don't have to worry about that
  
  grid_regs_df$Reg_fac <- as.factor(grid_regs_df$REGION)
  h0 <- model.matrix( ~ Reg_fac - 1, data=grid_regs_df)
  
  # Using grid_regs_df and h0 is indicator of regions
  grid_by_reg <- grid_regs_df %>% dplyr::mutate(Reg1 = h0[,1], Reg2 = h0[,2], Reg3 = h0[,3])
  sort_dat <- grid_by_reg %>% dplyr::arrange(Reg3, Reg2) # put stations in order of which regions they fall within
  
  return(sort_dat)
}

sim_formatted_data <- format_data_from_sim(test_sim_data, ws_regs, grid_regs_df)



# Model 1 -----------------------------------------------------------------
# Loading in median Hausdorff distance between regions
hMat <- readRDS(file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/hMat_med.rds")
hMat_miles <- hMat/5280


# Function to calculate a vector of return levels by region, given...
#   n_vec: a vector of the number of stations in each region
#   hMat: the extended Hausdorff distance matrix you want to use (make sure to convert to miles first)
#   c: the constant to define the distance between points within the same region (defaults to 1)
create_mat <- function(n_vec, hMat, c = 1){
  n1 = n_vec[1]
  n2 = n_vec[2]
  n3 = n_vec[3]
  # Creating blocks for block matrix
  one_to_one <- matrix(data = rep(c, n1^2), nrow = n1, ncol = n1)
  two_to_two <- matrix(data = rep(c, n2^2), nrow = n2, ncol = n2)
  three_to_three <- matrix(data = rep(c, n3^2), nrow = n3, ncol = n3)
  one_to_two <- matrix(data = rep(hMat[1,2], n1*n2), nrow = n1, ncol = n2)
  one_to_three <- matrix(data = rep(hMat[1,3], n1*n3), nrow = n1, ncol = n3)
  two_to_one <- matrix(data = rep(hMat[1,2], n1*n2), nrow = n2, ncol = n1)
  three_to_one <- matrix(data = rep(hMat[1,3], n1*n3), nrow = n3, ncol = n1)
  two_to_three <- matrix(data = rep(hMat[2,3], n2*n3), nrow = n2, ncol = n3)
  three_to_two <- matrix(data = rep(hMat[2,3], n2*n3), nrow = n3, ncol = n2)
  # putting blocks together
  one <- cbind(one_to_one, one_to_two, one_to_three)
  two <- cbind(two_to_one, two_to_two, two_to_three)
  three <- cbind(three_to_one, three_to_two, three_to_three)
  D_all <- rbind(one, two, three) # n x n matrix, where n = n1 + n2 + n3
  return(D_all)
}


# Function to create a symmetric matrix usable in the CAR model, given...
#   mat: a square (preferably symmetric) matrix
jitter_n_sym <- function(mat){
  # adding rnorm to jitter values
  mat_jit <- mat + matrix(rnorm(nrow(mat) * ncol(mat), sd = 0.1), ncol = ncol(mat))
  # Make symmetric
  mat_jit_sym <- mat_jit
  mat_jit_sym[upper.tri(mat_jit_sym)] <- t(mat_jit_sym)[upper.tri(mat_jit_sym)]
  return(mat_jit_sym)
}


# Shortcut function to call other functions to get to the symmetrical distance matrix 
# (still make sure to take reciprocal of distance matrix for weight matrix - a.k.a. inverse of each entry, NOT inverse of matrix)
get_sym_car_mat <- function(data, hMat, c=1){
  test_num <- summary(data$Reg_fac) # taking summary of factors gives the number of stations in each factor
  test_D <- create_mat(test_num, hMat, c)
  test_D_sym <- jitter_n_sym(test_D)
  return(test_D_sym)
}


# summary(sim_formatted_data$Reg_fac) # taking summary of factors gives the number of stations in each factor
D_sim <- get_sym_car_mat(sim_formatted_data, hMat_miles)

ptm <- proc.time()
car_shape_sim <- spatialreg::spautolm(shape ~ -1 + Reg1 + Reg2 + Reg3, data = sim_formatted_data, family="CAR",
                                     listw=mat2listw(1/D_sim))
proc.time() - ptm

ptm <- proc.time()
car_ln.scale_sim <- spatialreg::spautolm(log(scale) ~ -1 + Reg1 + Reg2 + Reg3, data = sim_formatted_data, family="CAR",
                                        listw=mat2listw(1/D_sim))
proc.time() - ptm

ptm <- proc.time()
car_rate_sim <- spatialreg::spautolm(rate ~ -1 + Reg1 + Reg2 + Reg3, data = sim_formatted_data, family="CAR",
                                    listw=mat2listw(1/D_sim))
proc.time() - ptm

summary(car_shape_sim)
summary(car_ln.scale_sim)
summary(car_rate_sim)

# Should just hold rate constant, it will be estimated correctly with SE of 3.9465e-19 (basically 0)



# Model 2 -----------------------------------------------------------------

# Change to perform on gridded data:

# Function to project data given...
#     sim_df: a data frame of the simulated grid with data and coordinates
project_sim_data <- function(sim_df){
  # "Projecting" grid data... coordinates are already projected, but have to remind it of proj4string
  coordinates(sim_df)=~x1+x2
  proj4string(sim_df) <- "+init=epsg:2278"
  return(sim_df)
}

sim_formatted_data_df <- sim_formatted_data
sim_formatted_data <- project_sim_data(sim_formatted_data)

# Fitting cross-variogram to real data (log(scale) and shape):
g <- gstat(id = "ln.scale", formula = log(scale)~1, data = sim_formatted_data)
g <- gstat(g, id = "shape", formula = shape~1, data = sim_formatted_data)
vg <- variogram(g)
cv.fit = fit.lmc(vg, g, vgm("Mat"))
plot(vg, model = cv.fit)
cv.fit

# cokriging
cv.fit$set=list(nocheck=1)
krig_regs_sim <- predict(cv.fit, newdata = ws_regs)
krig_regs_sim$ln.scale.pred
### estimates:
exp(krig_regs_sim$ln.scale.pred)
krig_regs_sim$shape.pred

# # testing run time
# ptm <- proc.time()
# krig_regs_timetest <- predict(cv.fit, newdata = ws_regs)
# (timer_krig_to_pol <- proc.time() - ptm)

# sp::zerodist(sim_formatted_data)

## Need to bring in rate parameter
# Fitting variogram to real data (rate)
g3 <- gstat(id = "rate", formula = rate~1, data = sim_formatted_data)
vg3 <- variogram(g3)
var.fit3 <- fit.variogram(vg3, vgm("Sph")) # singular model b/c all the same values
var.fit3
plot(vg3, model = var.fit3, main = "rate")

# cv.fit$set=list(nocheck=1)
# krig_regs_rate <- predict(var.fit3, newdata = ws_regs)
# krig_regs$ln.scale.pred

krig_regs_rate_sim <- krige(rate~1, sim_formatted_data, ws_regs, model = var.fit3)
krig_regs_rate_sim$var1.pred
krig_regs_rate_sim$var1.var

### Again, don't need to model rate because we are just holding it constant
# it will be estimated correctly with SE of 2.009807e-31 (basically 0)

rate_sim <- 795/14610
rate_sim1 <- krig_regs_rate_sim$var1.pred[1]
rate_sim2 <- car_rate_sim$fit$coefficients[[1]]

# Function to get parameter estimate matrix from kriging model fits, given...
#   krig_regs: the co-krige fit object for ln.scale and shape
###  rate_sim: constant rate set by simulation (0.05441478)
get_par_krig_sim <- function(krig_regs){    # }, rate_sim){
  krig_fit <- rbind(exp(krig_regs$ln.scale.pred)*(1 + (krig_regs$ln.scale.var/2)),    # updated with transformation/approximation
                    # exp(krig_regs$ln.scale.pred),   # old estimate
                    krig_regs$shape.pred)
                    # , rep(rate_sim, 3))
  rownames(krig_fit) <- c("scale", "shape") #, "rate")
  colnames(krig_fit) <- c("Reg1", "Reg2", "Reg3")
  return(krig_fit)
}

# Function to get parameter estimate matrix from CAR model fits, given...
#   ln.scale.fit: the CAR fit object for ln.scale
#   shape.fit: the CAR fit object for shape
###   rate_sim: constant rate set by simulation (0.05441478)
get_par_car_sim <- function(ln.scale.fit, shape.fit){   # }, rate_sim){
  car_fit <- rbind(exp(ln.scale.fit$fit$coef) * (1 + (diag(ln.scale.fit$fit$imat * ln.scale.fit$fit$s2)/2)),   # updated with transformation/approximation
                   # exp(ln.scale.fit$fit$coef),  # old estimate
                   shape.fit$fit$coef)
                   # , rep(rate_sim, 3))
  rownames(car_fit) <- c("scale", "shape") #, "rate")
  return(car_fit)
}

par_car_sim <- get_par_car_sim(car_ln.scale_sim, car_shape_sim)
par_krig_sim <- get_par_krig_sim(krig_regs_sim)

# Function to get parameter estimate matrix from consolidated data fits, given...
#   fitreg1: GPD fit object for Region 1
#   fitreg2: GPD fit object for Region 2
#   fitreg3: GPD fit object for Region 3
###   rate_sim: constant rate set by simulation (0.05441478)
get_par_consol_sim <- function(fitreg1, fitreg2, fitreg3){
  consol_fit <- rbind(c(fitreg1$results$par[1], fitreg2$results$par[1], fitreg3$results$par[1]), # scale
                      c(fitreg1$results$par[2], fitreg2$results$par[2], fitreg3$results$par[2])) # shape
                      # , c(fitreg1$rate, fitreg2$rate, fitreg3$rate)) # rate
  rownames(consol_fit) <- c("scale", "shape") #, "rate")
  colnames(consol_fit) <- c("Reg1", "Reg2", "Reg3")
  return(consol_fit)
}


# Function to generate a list of varcov matrices for each region, given...
#   ln.scale.fit: the CAR fit object for ln.scale
#   shape.fit: the CAR fit object for shape
make_varcov_mats_car <- function(ln.scale.fit, shape.fit){
  varcov_list <- list()
  # reg1 <- reg2 <- reg3 <- matrix(data=NA, nrow=2, ncol=2)
  shape_vcov <- shape.fit$fit$imat * shape.fit$fit$s2
  shape.var <- diag(shape_vcov)
  ln.scale_vcov <- ln.scale.fit$fit$imat * ln.scale.fit$fit$s2
  ln.scale.var <- diag(ln.scale_vcov)
  ln.scale.pred <- ln.scale.fit$fit$coefficients
  
  for(i in 1:length(shape.fit$fit$coefficients)){
    mat <- matrix(data=NA, nrow=2, ncol=2)
    mat[1,1] <- exp(ln.scale.pred[i])^2 * ln.scale.var[i] # doing transformation
    mat[2,2] <- shape.var[i]
    mat[1,2] <- mat[2,1] <- 0 # setting cov to 0 since we model separately... a simplifying assumption that will make our CIs wider
    varcov_list[[i]] <- mat
  }
  return(varcov_list)
}

# Function to get SE's from CAR model fits
get_se_car_sim <- function(ln.scale.fit, shape.fit){  # still returns log(scale)...
  mat <- NULL
  ln.scale_vcov <- ln.scale.fit$fit$imat * ln.scale.fit$fit$s2
  ln.scale.se <- sqrt(diag(ln.scale_vcov))
  ln.scale.pred <- ln.scale.fit$fit$coefficients
  # exp(ln.scale.pred[i])^2 * ln.scale.var[i] # VAR
  # exp(ln.scale.pred[i]) * ln.scale.se[i] # SE
  mat <- rbind(mat, exp(ln.scale.pred) * ln.scale.se)
  shape_vcov <- shape.fit$fit$imat * shape.fit$fit$s2
  mat <- rbind(mat, sqrt(diag(shape_vcov)))
  rownames(mat) <- c("scale", "shape") #, "rate")
  colnames(mat) <- c("Reg1", "Reg2", "Reg3")
  return(mat)
}

se_car_sim <- get_se_car_sim(car_ln.scale_sim, car_shape_sim)

# # Function to get SE's from CAR model fits
# get_se_krig_sim <- function(ln.scale.var, ln.scale.pred, shape.var, cov.ln.scale.shape=0){
#   get_se_from_varcov(make_varcov_mats_krig(ln.scale.var, ln.scale.pred, shape.var, cov.ln.scale.shape))
# }

get_se_from_varcov(make_varcov_mats_car(car_ln.scale_sim, car_shape_sim))


# Function to generate a list of varcov matrices for each region, given...  (useful for kriging model)
#   ln.scale.var: the vector of log(scale) variances
#   ln.scale.pred: the vector of log(scale) estimates
#   shape.var: the vector of shape variances
#   cov.ln.scale.shape: the vector of covariances between log(scale) and shape
make_varcov_mats_krig <- function(ln.scale.var, ln.scale.pred, shape.var, cov.ln.scale.shape){      # still returns cov log(scale)...
  varcov_list <- list()
  for(i in 1:length(ln.scale.var)){
    mat <- matrix(data=NA, nrow=2, ncol=2)
    mat[1,1] <- exp(ln.scale.pred[i])^2 * ln.scale.var[i] # doing transformation  
    mat[2,2] <- shape.var[i]
    mat[1,2] <- mat[2,1] <- exp(ln.scale.pred[i]) * cov.ln.scale.shape[i] # transforming the covariance
    varcov_list[[i]] <- mat
  }
  return(varcov_list)
}

# renamed from give_se
# Function to get SE's from a list of varcov matrices
get_se_from_varcov <- function(varcov_list){ 
  se_vec <- NULL
  for (i in 1:length(varcov_list)){
    se_vec <- cbind(se_vec, sqrt(diag(varcov_list[[i]])))
  }
  rownames(se_vec) <- c("scale", "shape")
  colnames(se_vec) <- c("Reg1", "Reg2", "Reg3")
  return(se_vec)
}

# Function to get SE's from Kriging model fits
get_se_krig_sim <- function(ln.scale.var, ln.scale.pred, shape.var, cov.ln.scale.shape=0){
  get_se_from_varcov(make_varcov_mats_krig(ln.scale.var, ln.scale.pred, shape.var, cov.ln.scale.shape))
}

se_krig_sim <- get_se_krig_sim(krig_regs_sim$ln.scale.var, krig_regs_sim$ln.scale.pred, krig_regs_sim$shape.var)

scale.means
exp(ln.scale.means)
shape.means



# Repeating Simulation for More Iterations --------------------------------

# test_sims is the list of simulated data 1-100. Take the first 10

car_shape_sim_list <- list()
car_ln.scale_sim_list <- list()
krig_regs_sim_list <- list()
par_car_sim_list <- list()
par_krig_sim_list <- list()
  

# Select an index and loop through:
i = 2
ptm <- proc.time()
sim_formatted_data_i <- format_data_from_sim(test_sims_constant[[i]], ws_regs, grid_regs_df)  # this step takes the longest
# sim_formatted_data_i <- format_data_from_sim(test_sims[[i]], ws_regs, grid_regs_df)  # this step takes the longest
proc.time() - ptm

D_sim_i <- get_sym_car_mat(sim_formatted_data_i, hMat_miles)

ptm <- proc.time()
car_shape_sim_i <- spatialreg::spautolm(shape ~ -1 + Reg1 + Reg2 + Reg3, data = sim_formatted_data_i, family="CAR",
                                      listw=mat2listw(1/D_sim_i))
proc.time() - ptm

ptm <- proc.time()
car_ln.scale_sim_i <- spatialreg::spautolm(log(scale) ~ -1 + Reg1 + Reg2 + Reg3, data = sim_formatted_data_i, family="CAR",
                                         listw=mat2listw(1/D_sim_i))
proc.time() - ptm

summary(car_shape_sim_i)
summary(car_ln.scale_sim_i)


### Model 2
sim_formatted_data_i_df <- sim_formatted_data_i
sim_formatted_data_i <- project_sim_data(sim_formatted_data_i)

# Fitting cross-variogram to real data (log(scale) and shape):
g <- gstat(id = "ln.scale", formula = log(scale)~1, data = sim_formatted_data_i)
g <- gstat(g, id = "shape", formula = shape~1, data = sim_formatted_data_i)
vg <- variogram(g)
cv.fit = fit.lmc(vg, g, vgm("Mat"))
plot(vg, model = cv.fit)
cv.fit

# cokriging
cv.fit$set=list(nocheck=1)
ptm <- proc.time()
krig_regs_sim_i <- predict(cv.fit, newdata = ws_regs)
proc.time() - ptm

## Calculating parameters
par_car_sim_i <- get_par_car_sim(car_ln.scale_sim_i, car_shape_sim_i)
par_krig_sim_i <- get_par_krig_sim(krig_regs_sim_i)

## Saving fits to list
car_shape_sim_list[[i]] <- car_shape_sim_i
car_ln.scale_sim_list[[i]] <- car_ln.scale_sim_i
krig_regs_sim_list[[i]] <- krig_regs_sim_i
par_car_sim_list[[i]] <- par_car_sim_i
par_krig_sim_list[[i]] <- par_krig_sim_i

# par_car_sim_list_old <- par_car_sim_list
# par_krig_sim_list_old <- par_krig_sim_list

### create lists ONLY IF you have not already saved data manually above
car_shape_sim_list <- list()
car_ln.scale_sim_list <- list()
krig_regs_sim_list <- list()
par_car_sim_list <- list()
par_krig_sim_list <- list()


### Set up in a loop - change indicies as needed
ptm <- proc.time()
for (i in 21:50){
  sim_formatted_data_i <- format_data_from_sim(test_sims_constant[[i]], ws_regs, grid_regs_df)  # this step usually takes the longest
  
  ### Model 1
  D_sim_i <- get_sym_car_mat(sim_formatted_data_i, hMat_miles)
  car_shape_sim_i <- spatialreg::spautolm(shape ~ -1 + Reg1 + Reg2 + Reg3, data = sim_formatted_data_i, family="CAR",
                                          listw=mat2listw(1/D_sim_i))
  car_ln.scale_sim_i <- spatialreg::spautolm(log(scale) ~ -1 + Reg1 + Reg2 + Reg3, data = sim_formatted_data_i, family="CAR",
                                             listw=mat2listw(1/D_sim_i))
  
  ### Model 2
  sim_formatted_data_i_df <- sim_formatted_data_i
  sim_formatted_data_i <- project_sim_data(sim_formatted_data_i)
  
  # Fitting cross-variogram to real data (log(scale) and shape):
  g <- gstat(id = "ln.scale", formula = log(scale)~1, data = sim_formatted_data_i)
  g <- gstat(g, id = "shape", formula = shape~1, data = sim_formatted_data_i)
  vg <- variogram(g)
  cv.fit = fit.lmc(vg, g, vgm("Mat"))
  plot(vg, model = cv.fit)
  cv.fit
  
  # cokriging
  cv.fit$set=list(nocheck=1)
  krig_regs_sim_i <- predict(cv.fit, newdata = ws_regs)
  
  
  ### Calculating parameters
  par_car_sim_i <- get_par_car_sim(car_ln.scale_sim_i, car_shape_sim_i)
  par_krig_sim_i <- get_par_krig_sim(krig_regs_sim_i)
  
  ## Saving fits to list
  car_shape_sim_list[[i]] <- car_shape_sim_i
  car_ln.scale_sim_list[[i]] <- car_ln.scale_sim_i
  krig_regs_sim_list[[i]] <- krig_regs_sim_i
  par_car_sim_list[[i]] <- par_car_sim_i
  par_krig_sim_list[[i]] <- par_krig_sim_i
}
proc.time() - ptm


# saveRDS(par_car_sim_list, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/par_car_sim_list.rds")
# saveRDS(par_krig_sim_list, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/par_krig_sim_list.rds")
# saveRDS(round.scale.means, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/true_scale.rds")
# saveRDS(round.shape.means, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/true_shape.rds")
# saveRDS(par_car_sim_list, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/par_car_sim_list_50.rds")
# saveRDS(par_krig_sim_list, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/par_krig_sim_list_50.rds")
# saveRDS(car_ln.scale_sim_list, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/car_ln.scale_sim_list_50.rds")
# saveRDS(car_shape_sim_list, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/car_shape_sim_list_50.rds")
# saveRDS(krig_regs_sim_list, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/krig_regs_sim_list_50.rds")




### Function to calculate MSE
mse <- function(data, truth){
  out = mean((data-truth)^2)
  return(out)
}

### Function to calculate RMSE
rmse <- function(data, truth){
  out = mse(data, truth)
  return(sqrt(out))
}

### Function to calculate MAE
mae <- function(data, truth){
  out = mean(abs(data-truth))
  return(out)
}

mse(c(4,2,4,5,4), 4)
rmse(c(4,2,4,5,4), 4)
mae(c(4,2,4,5,4), 4)

### Function to create simulation output table

create_sim_table_scale <- function(par_car_sim_list, par_krig_sim_list, par_consol_sim_list, length){
  car_reg1 <- car_reg2 <- car_reg3 <- NULL
  krig_reg1 <- krig_reg2 <- krig_reg3 <- NULL
  consol_reg1 <- consol_reg2 <- consol_reg3 <- NULL
  for (i in 1:length){
    car_reg1[i] <- par_car_sim_list[[i]][1,1]
    car_reg2[i] <- par_car_sim_list[[i]][1,2]
    car_reg3[i] <- par_car_sim_list[[i]][1,3]
    krig_reg1[i] <- par_krig_sim_list[[i]][1,1]
    krig_reg2[i] <- par_krig_sim_list[[i]][1,2]
    krig_reg3[i] <- par_krig_sim_list[[i]][1,3]
    consol_reg1[i] <- par_consol_sim_list[[i]][1,1]
    consol_reg2[i] <- par_consol_sim_list[[i]][1,2]
    consol_reg3[i] <- par_consol_sim_list[[i]][1,3]
  }
  # plus the truth as the first column
  scale_car <- rbind( c(round.scale.means[1], mean(car_reg1), rmse(car_reg1, round.scale.means[1]), mae(car_reg1, round.scale.means[1])),
                      c(round.scale.means[2], mean(car_reg2), rmse(car_reg2, round.scale.means[2]), mae(car_reg2, round.scale.means[2])),
                      c(round.scale.means[3], mean(car_reg3), rmse(car_reg3, round.scale.means[3]), mae(car_reg3, round.scale.means[3])) )
  scale_krig <- rbind(c(mean(krig_reg1), rmse(krig_reg1, round.scale.means[1]), mae(krig_reg1, round.scale.means[1])),
                      c(mean(krig_reg2), rmse(krig_reg2, round.scale.means[2]), mae(krig_reg2, round.scale.means[2])),
                      c(mean(krig_reg3), rmse(krig_reg3, round.scale.means[3]), mae(krig_reg3, round.scale.means[3])) )
  scale_consol <- rbind( c(mean(consol_reg1), rmse(consol_reg1, round.scale.means[1]), mae(consol_reg1, round.scale.means[1])),
                         c(mean(consol_reg2), rmse(consol_reg2, round.scale.means[2]), mae(consol_reg2, round.scale.means[2])),
                         c(mean(consol_reg3), rmse(consol_reg3, round.scale.means[3]), mae(consol_reg3, round.scale.means[3])) )
  scale_out <- cbind(scale_car, scale_krig, scale_consol)
  colnames(scale_out) <- c("Truth", rep(c("Mean", "RMSE", "MAE"), 3))
  rownames(scale_out) <- c("Reg1", "Reg2", "Reg3")
  return(scale_out)
}

# Truth
round.scale.means
round.shape.means

# Shape
create_sim_table_shape <- function(par_car_sim_list, par_krig_sim_list, par_consol_sim_list, length){
  car_reg1 <- car_reg2 <- car_reg3 <- NULL
  krig_reg1 <- krig_reg2 <- krig_reg3 <- NULL
  consol_reg1 <- consol_reg2 <- consol_reg3 <- NULL
  for (i in 1:length){
    car_reg1[i] <- par_car_sim_list[[i]][2,1]
    car_reg2[i] <- par_car_sim_list[[i]][2,2]
    car_reg3[i] <- par_car_sim_list[[i]][2,3]
    krig_reg1[i] <- par_krig_sim_list[[i]][2,1]
    krig_reg2[i] <- par_krig_sim_list[[i]][2,2]
    krig_reg3[i] <- par_krig_sim_list[[i]][2,3]
    consol_reg1[i] <- par_consol_sim_list[[i]][2,1]
    consol_reg2[i] <- par_consol_sim_list[[i]][2,2]
    consol_reg3[i] <- par_consol_sim_list[[i]][2,3]
  }
  # plus the truth as the first column
  shape_car <- rbind( c(round.shape.means[1], mean(car_reg1), rmse(car_reg1, round.shape.means[1]), mae(car_reg1, round.shape.means[1])),
                      c(round.shape.means[2], mean(car_reg2), rmse(car_reg2, round.shape.means[2]), mae(car_reg2, round.shape.means[2])),
                      c(round.shape.means[3], mean(car_reg3), rmse(car_reg3, round.shape.means[3]), mae(car_reg3, round.shape.means[3])) )
  shape_krig <- rbind(c(mean(krig_reg1), rmse(krig_reg1, round.shape.means[1]), mae(krig_reg1, round.shape.means[1])),
                      c(mean(krig_reg2), rmse(krig_reg2, round.shape.means[2]), mae(krig_reg2, round.shape.means[2])),
                      c(mean(krig_reg3), rmse(krig_reg3, round.shape.means[3]), mae(krig_reg3, round.shape.means[3])) )
  shape_consol <- rbind( c(mean(consol_reg1), rmse(consol_reg1, round.shape.means[1]), mae(consol_reg1, round.shape.means[1])),
                         c(mean(consol_reg2), rmse(consol_reg2, round.shape.means[2]), mae(consol_reg2, round.shape.means[2])),
                         c(mean(consol_reg3), rmse(consol_reg3, round.shape.means[3]), mae(consol_reg3, round.shape.means[3])) )
  shape_out <- cbind(shape_car, shape_krig, shape_consol)
  colnames(shape_out) <- c("Truth", rep(c("Mean", "RMSE", "MAE"), 3))
  rownames(shape_out) <- c("Reg1", "Reg2", "Reg3")
  return(shape_out)
}

# sim_scale_table_10 <- create_sim_table_scale(par_car_sim_list, par_krig_sim_list, length=10)
# sim_shape_table_10 <- create_sim_table_shape(par_car_sim_list, par_krig_sim_list, length=10)
# sim_scale_table <- create_sim_table_scale(par_car_sim_list, par_krig_sim_list, length=50)
# sim_shape_table <- create_sim_table_shape(par_car_sim_list, par_krig_sim_list, length=50)

sim_scale_table <- create_sim_table_scale(par_car_sim_list, par_krig_sim_list, par_consol_sim_list, length=50)
sim_shape_table <- create_sim_table_shape(par_car_sim_list, par_krig_sim_list, par_consol_sim_list, length=50)


# saveRDS(sim_scale_table, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/sim_scale_table.rds")
# saveRDS(sim_shape_table, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/sim_shape_table.rds")



# Model 3 -----------------------------------------------------------------



# #### Code from Karen Yuan
# ## Finding intersection of points and regions
# # sp::over will let us get what region each station falls within (or NA if none) in order of the 601 stations
# region_intersect <- sp::over(station_info, regions)
# class(region_intersect$REGION)
# 
# ## Aggregate station info into data by each region
# #create function that will retrieve relevant stations' columns for subsetting purposes, given a region
# num_stations <- 601
# get_columns_index <- function(intersect_vector, region) {
#   stations <- seq(1, num_stations, 1)[intersect_vector == region] #get station number for those in "region", there will be NA's
#   col_nums <- stations[!is.na(stations)] + 1 #add 1 to column indices since column 1 is Date
#   col_nums <- c(1, col_nums) #include Date column
#   return(col_nums)
# }
# #make 6 subsets of the data to easily collect a max value from each date for each region subset (region1_alt technically a duplicate)
# rain_region1 <- rain_data[, get_columns_index(region_intersect, 1)]
# rain_region2 <- rain_data[, get_columns_index(region_intersect, 2)]
# rain_region3 <- rain_data[, get_columns_index(region_intersect, 3)]
# rain_region1_alt <- rain_data[, get_columns_index(region_intersect_alt, 1)]
# rain_region2_alt <- rain_data[, get_columns_index(region_intersect_alt, 2)]
# rain_region3_alt <- rain_data[, get_columns_index(region_intersect_alt, 3)]
# 
# ## Populate consolidated time series
# #create final data frame to populate
# max_precip_regions <- data.frame(matrix(ncol=7, nrow=nrow(rain_data)))
# colnames(max_precip_regions) <- c("Date", "region1", "region2", "region3", "region1_alt", "region2_alt", "region3_alt")
# max_precip_regions$Date <- rain_data$Date
# 
# #create function that will take in a row from each region subset and find the max value
# get_max_value <- function(values_vector) {
#   if (all(is.na(values_vector))) {
#     return(NA)
#   }
#   else {
#     return(max(values_vector, na.rm=TRUE))
#   }
# }
# #loop through every date in the data set and get the max value for each region subset for that day
# for (i in 1:nrow(rain_data)) {
#   max_precip_regions$region1[i] <- get_max_value(rain_region1[i, -1]) 
#   max_precip_regions$region2[i] <- get_max_value(rain_region2[i, -1]) 
#   max_precip_regions$region3[i] <- get_max_value(rain_region3[i, -1]) 
#   max_precip_regions$region1_alt[i] <- get_max_value(rain_region1_alt[i, -1]) 
#   max_precip_regions$region2_alt[i] <- get_max_value(rain_region2_alt[i, -1]) 
#   max_precip_regions$region3_alt[i] <- get_max_value(rain_region3_alt[i, -1]) 
# }



#### Fitting Model 3 (Regional Max)

# Edit code below to run Model 3

get_columns_index <- function(grid_regs_df, region) {
  col_nums <- which(grid_regs_df$REGION == region)
  return(col_nums)
}
# These are the column indicies of which stations are in each region
reg1_col <- get_columns_index(grid_regs_df, 1)
reg2_col <- get_columns_index(grid_regs_df, 2)
reg3_col <- get_columns_index(grid_regs_df, 3)
# Run this once before looping beacuse it will be the same for all

#create function that will take in a row from each region subset and find the max value
get_max_value <- function(values_vector) {
  if (all(is.na(values_vector))) {
    return(NA)
  }
  else {
    return(max(values_vector, na.rm=TRUE))
  }
}


### create lists ONLY IF you have not already saved some data to it
consol_data_sim_list <- list()
fitreg1_sim_list <- list()
fitreg2_sim_list <- list()
fitreg3_sim_list <- list()
par_consol_sim_list <- list()

### Set up in a loop - change indicies as needed
ptm <- proc.time()
for (i in 31:50){
  ### Model 3
  data_i <- test_sims_constant_ordered[[i]]
  rain_reg1 <- data_i[, reg1_col]
  rain_reg2 <- data_i[, reg2_col]
  rain_reg3 <- data_i[, reg3_col]
  
  consol_data_sim_i <- data.frame(matrix(ncol=3, nrow=nrow(data_i)))
  colnames(consol_data_sim_i) <- c("region1", "region2", "region3")
  # max_precip_regions <- data.frame(matrix(ncol=3, nrow=nrow(data_i)))
  # colnames(max_precip_regions) <- c("region1", "region2", "region3")
  
  #loop through every date in the data set and get the max value for each region subset for that day
  for (j in 1:nrow(data_i)) {
    consol_data_sim_i$region1[j] <- get_max_value(rain_reg1[j, ]) 
    consol_data_sim_i$region2[j] <- get_max_value(rain_reg2[j, ]) 
    consol_data_sim_i$region3[j] <- get_max_value(rain_reg3[j, ]) 
  }
  
  fitreg1_sim_i <- extRemes::fevd(consol_data_sim_i$region1, threshold=thresh, type="GP", method="MLE")
  fitreg2_sim_i <- extRemes::fevd(consol_data_sim_i$region2, threshold=thresh, type="GP", method="MLE")
  fitreg3_sim_i <- extRemes::fevd(consol_data_sim_i$region3, threshold=thresh, type="GP", method="MLE")
  
  ## Calculating parameters
  par_consol_sim_i <- get_par_consol_sim(fitreg1_sim_i, fitreg2_sim_i, fitreg3_sim_i)
  
  ## Saving fits to list
  consol_data_sim_list[[i]] <- consol_data_sim_i
  fitreg1_sim_list[[i]] <- fitreg1_sim_i
  fitreg2_sim_list[[i]] <- fitreg2_sim_i
  fitreg3_sim_list[[i]] <- fitreg3_sim_i
  par_consol_sim_list[[i]] <- par_consol_sim_i
}
proc.time() - ptm


# saveRDS(par_consol_sim_list, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/par_consol_sim_list_50.rds")
# saveRDS(consol_data_sim_list, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/consol_data_sim_list_50.rds")
# saveRDS(fitreg1_sim_list, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/fitreg1_sim_list_50.rds")
# saveRDS(fitreg2_sim_list, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/fitreg2_sim_list_50.rds")
# saveRDS(fitreg3_sim_list, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/fitreg3_sim_list_50.rds")



## save to compare
# test_par <- par_consol_sim_i
# View(par_consol_sim_list[[1]])
