#######################
# Testing the models on real data
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

thresh <- 253

# Loading our data --------------------------------------------------------

setwd("~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Test")
wd <- getwd()

ws_regs <- rgdal::readOGR(dsn = paste0(wd, "/watershed_region_updated"), layer = "watershed_region_updated")
plot(ws_regs)

# Alternate version of 3 regions
ws_reg_alt <- rgdal::readOGR(dsn = paste0(wd, "/watershed_region_alt"), layer = "watershed_region_alt")
plot(ws_reg_alt)

# load pre-cleaned spatial points data frame
stations_sub <- readRDS(file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/stations_sub.rds")
stations_sub_df <- readRDS(file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/stations_sub_df.rds")

# avg parameter values by region
ws_reg_avg <- readRDS("~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/ws_reg_avg.rds")


# Function to project data given...
#     df: a data frame with data and coordinates
project_data <- function(df){
  coordinates(df)=~long+lat
  proj4string(df) <- "+proj=longlat +datum=WGS84"
  df <- spTransform(df, CRS("+init=epsg:2278")) # projecting data
  is.projected(df)
  return(df)
}

is.projected(stations_sub)
is.projected(ws_regs)
# set regions and points to have the same projection
proj4string(ws_regs) <- proj4string(stations_sub)


## stations_sub is already using data fits from last 40-year window (1981-2020)
## but if want raccess to all fits - load the following window.rds file
## (access "window_1day_dclust_updated.rds" from zipped file in Box - should be 1.67 GB unzipped)

# setwd("~/Documents/HCFCD Watershed Data") # wherever you have the following file saved 
# window <- readRDS("window_1day_dclust_updated.rds")  # list of lists (82 x 601) of GPD fits



# Model 1 - Random Effects Model using CAR and median Hausdorff distance --------------------------------

## this was already run in "Simulations_using_gstat.R", we will bring in the saved estimates here
# car_fit <- readRDS(file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/car_fit.rds")

# Finding return levels given the parameters:
car_fit

# extRemes::return.level - calculates return levels given an fevd object
# extRemes::rlevd() - calculates return levels given a distribution and parameter estimates

# Function to calculate a vector of return levels by region, given...
#   parmat: the parameter matrix for the regions
#   return_period: a return period
pars_to_rl <- function(par_mat, return_period){
  RL_vec <- NULL
  for(reg in 1:3){
    scale <- par_mat[1, reg]
    shape <- par_mat[2, reg]
    rate  <- par_mat[3, reg]
    RL_vec[reg] <- extRemes::rlevd(period=return_period, type="GP", scale=scale, shape=shape, rate=rate, threshold=thresh)
  }
  return(RL_vec)
}

pars_to_rl(car_fit, 25)/254
pars_to_rl(car_fit, 100)/254
pars_to_rl(car_fit, 500)/254

# issue - extRemes::rlevd function cannot calculate confidence intervals



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

hMat <- readRDS(file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/hMat_med.rds")
hMat_miles <- hMat/5280


### Trying CAR-type model on a different window!
# startyear = 1900 + window - 1
# endyear = 1900 + window + 39 - 1 = startyear + 39
# ex: 1900 + 82 - 1 = 1981
# ex: 1900 + 82 + 39 - 1 = 1981 + 39 = 2020

# NEVERMIND, window 2 only has 3 stations
# Let's do windows: 2, 42, and 82
# (1901-1940), (1941-1980), (1981-2020)

# Let's do windows: 22, 52, and 82
# (1921-1960), (1951-1990), (1981-2020)

# 1. load in window and turn into stations_sub format with columns Reg1, Reg2, Reg3
# 2. order data by columns using: 
      # # Using stations_sub_df and h0 is indicator of regions
      # stations_sub_by_reg <- stations_sub_df %>% dplyr::mutate(Reg1 = h0[,1], Reg2 = h0[,2], Reg3 = h0[,3])
      # sort_dat <- stations_sub_by_reg %>% dplyr::arrange(Reg3, Reg2) # put stations in order of which regions they fall within
# 3. Make big distance matrix using functions above

## Testing
# window_test <- window[[82]]
# identical(window_test[[3]], window[[82]][[3]])
stations <- read.csv("~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/station_info.csv")
stat_nos <- stations[,1]

# Function to format data to stations within the regions, given...
#   window_fit: the window list subset to the moving window of choice (e.g. input window[[82]] for the final 40-yr window)
#   reg_polygons: the SpatialPolygonsDataFrame for the regions of interest
format_data <- function(window_fit, reg_polygons){
  # Accessing Model Fits and saving the parameter values (scale, shape, rate) for each station
  ws_scale <- NULL
  ws_shape <- NULL
  ws_rate  <- NULL
  j = 1
  for(i in stat_nos){
    fit <- window_fit[[i]] 
    if(is.na(fit)){
      ws_scale[j] <- ws_shape[j] <- ws_rate[j]  <- NA
    }else{
      ws_scale[j] <- fit$results$par[1]
      ws_shape[j] <- fit$results$par[2]
      ws_rate[j]  <- fit$rate
    }
    j = j+1
  }
  # add new location (long/lat) columns that we will transform to coordinates, as well as parameter values
  stations_df <- stations %>%
    dplyr::mutate(scale = ws_scale, shape = ws_shape, rate = ws_rate, long = LONGITUDE, lat = LATITUDE)
  
  stations <- project_data(stations_df)
  region_int <- sp::over(stations, reg_polygons)
  region_int$REGION <- as.factor(region_int$REGION)
  h <- model.matrix( ~ REGION - 1, data=region_int)
  inds <- as.integer(rownames(h))
  
  stations_sub_df <- stations_df %>%
    dplyr::filter(STAT_NO %in% inds) %>%
    dplyr::filter(!is.na(scale))
  
  stations_sub <- project_data(stations_sub_df)
  
  ## Have to redo intersect, because now NAs are removed. 
  ## If remove them before the intersect, the station indicies get messed up. So have to do intersect both before and after
  region_int <- sp::over(stations_sub, reg_polygons)
  region_int$REGION <- as.factor(region_int$REGION)
  h0 <- model.matrix( ~ REGION - 1, data=region_int)
  
  # Using stations_sub_df and h0 is indicator of regions
  stations_sub_by_reg <- stations_sub_df %>% dplyr::mutate(Reg_fac = region_int$REGION, Reg1 = h0[,1], Reg2 = h0[,2], Reg3 = h0[,3])
  sort_dat <- stations_sub_by_reg %>% dplyr::arrange(Reg3, Reg2) # put stations in order of which regions they fall within
  
  return(sort_dat)
}


test_format_data <- format_data(window[[82]], ws_regs)

test_num <- summary(test_format_data$Reg_fac) # taking summary of factors gives the number of stations in each factor
test_D <- create_mat(test_num, hMat_miles)
test_D_sym <- jitter_n_sym(test_D)

# Shortcut function to call other functions to get to the symmetrical distance matrix 
# (still make sure to take reciprocal of distance matrix for weight matrix - a.k.a. inverse of each entry, NOT inverse of matrix)
get_sym_car_mat <- function(data, hMat, c=1){
  test_num <- summary(data$Reg_fac) # taking summary of factors gives the number of stations in each factor
  test_D <- create_mat(test_num, hMat, c)
  test_D_sym <- jitter_n_sym(test_D)
  return(test_D_sym)
}


car_shape <- spatialreg::spautolm(shape ~ -1 + Reg1 + Reg2 + Reg3, data = test_format_data, family="CAR",
                                       listw=mat2listw(1/test_D_sym))
car_ln.scale <- spatialreg::spautolm(log(scale) ~ -1 + Reg1 + Reg2 + Reg3, data = test_format_data, family="CAR",
                                          listw=mat2listw(1/test_D_sym))
car_rate <- spatialreg::spautolm(rate ~ -1 + Reg1 + Reg2 + Reg3, data = test_format_data, family="CAR",
                                      listw=mat2listw(1/test_D_sym))


# dat_win_2 <- format_data(window[[2]], ws_regs) # only 3 stations...
dat_win_42 <- format_data(window[[42]], ws_regs) # 35
dat_win_62 <- format_data(window[[62]], ws_regs) # 26
dat_win_82 <- format_data(window[[82]], ws_regs) # 149

dat_win_22 <- format_data(window[[22]], ws_regs) # 34
dat_win_52 <- format_data(window[[52]], ws_regs) # 28
dat_win_82 <- format_data(window[[82]], ws_regs) # 149


D_22 <- get_sym_car_mat(dat_win_22, hMat_miles)
D_52 <- get_sym_car_mat(dat_win_52, hMat_miles)
D_82 <- get_sym_car_mat(dat_win_82, hMat_miles)

car_shape_22 <- spatialreg::spautolm(shape ~ -1 + Reg1 + Reg2 + Reg3, data = dat_win_22, family="CAR",
                                  listw=mat2listw(1/D_22))
ptm <- proc.time()
car_ln.scale_22 <- spatialreg::spautolm(log(scale) ~ -1 + Reg1 + Reg2 + Reg3, data = dat_win_22, family="CAR",
                                     listw=mat2listw(1/D_22))
proc.time() - ptm
ptm <- proc.time()
car_rate_22 <- spatialreg::spautolm(rate ~ -1 + Reg1 + Reg2 + Reg3, data = dat_win_22, family="CAR",
                                 listw=mat2listw(1/D_22))
proc.time() - ptm

summary(car_shape_22)




# Model 2 - Kriging and Aggregation ---------------------------------------------------------------

# Fitting cross-variogram to real data (log(scale) and shape):
g <- gstat(id = "ln.scale", formula = log(scale)~1, data = stations_sub)
g <- gstat(g, id = "shape", formula = shape~1, data = stations_sub)
vg <- variogram(g)
cv.fit = fit.lmc(vg, g, vgm("Mat"))
plot(vg, model = cv.fit)
cv.fit

# cokriging
cv.fit$set=list(nocheck=1)
krig_regs <- predict(cv.fit, newdata = ws_regs)
krig_regs$ln.scale.pred
### estimates:
exp(krig_regs$ln.scale.pred)
krig_regs$shape.pred

# testing run time
ptm <- proc.time()
krig_regs_timetest <- predict(cv.fit, newdata = ws_regs)
(timer_krig_to_pol <- proc.time() - ptm)


krig_points <- predict(cv.fit, newdata = stations_sub)

### create a grid (of points)
# library(raster)
grid <- sp::makegrid(ws_regs, cellsize = 2000) # cellsize in map units
grid <- SpatialPoints(grid, proj4string = CRS(proj4string(ws_regs)))
grid <- grid[ws_regs, ] # grid points only within the polygon
plot(ws_regs)
plot(grid, pch = ".", add = T)
# gridded(grid) = TRUE # turn into spatial pixels data frame

# to grid
krig_grid <- predict(cv.fit, newdata = grid)

# from Eduardo - function calculating avgs for regions from grid
compute_means <- function(co_krig, regions) {
  num_points <- length(co_krig)
  num_regions <- length(regions)
  ln.scale_avg <- rep(0, num_regions)
  shape_avg <- rep(0, num_regions)
  for (i in 1:num_regions) {
    print(paste0("Outer iter: ", i))
    region <- regions[i, ]
    counter <- 0
    pct_old <- 0
    for (j in 1:num_points) {
      pct <- 100 * j / num_points
      point <- co_krig[j, ]
      if (is.null(rgeos::gDifference(point, region, byid = T))) {
        ln.scale_avg[i] <- ln.scale_avg[i] + co_krig@data$ln.scale.pred[j]
        shape_avg[i] <- shape_avg[i] + co_krig@data$shape.pred[j]
        counter <- counter + 1
      }
      if (pct - pct_old > 5) {
        print(paste0(pct, "% of outer iter completed"))
        pct_old <- pct
      }
    }
    if (counter != 0) {
      ln.scale_avg[i] <- ln.scale_avg[i] / counter
      shape_avg[i] <- shape_avg[i] / counter
    }
  }
  avg_mat <- rbind(ln.scale_avg, shape_avg)
  return(avg_mat)
}

ptm <- proc.time()
mat <- compute_means(krig_grid, ws_regs)
(timer <- proc.time() - ptm)
mat



# plot(ws_regs)
# plot(skrig, pch = "•", add = T)
# spplot(skrig) 
# spplot(skrigws)
# spplot(skrigws[3]) # ln.scale.pred
# spplot(skrigws[5]) # shape.pred


## Need to bring in rate parameter
# Fitting variogram to real data (rate)
g3 <- gstat(id = "rate", formula = rate~1, data = stations_sub)
vg3 <- variogram(g3)
var.fit3 <- fit.variogram(vg3, vgm("Sph"))
var.fit3
plot(vg3, model = var.fit3, main = "rate")

# cv.fit$set=list(nocheck=1)
# krig_regs_rate <- predict(var.fit3, newdata = ws_regs)
# krig_regs$ln.scale.pred

krig_regs_rate <- krige(rate~1, stations_sub, ws_regs, model = var.fit3)
krig_regs_rate$var1.pred

# to grid
krig_grid_rate <- krige(rate~1, stations_sub, newdata = grid, model = var.fit3)

# krig_regs_rate_test <- gstat::krige(rate~1, stations_sub, newdata = ws_regs, model = var.fit3)
# krig_regs_rate_test$var1.pred

# function calculating rate avgs for regions from grid
compute_means_rate <- function(co_krig, regions) {
  num_points <- length(co_krig)
  num_regions <- length(regions)
  rate_avg <- rep(0, num_regions)
  # shape_avg <- rep(0, num_regions)
  for (i in 1:num_regions) {
    print(paste0("Outer iter: ", i))
    region <- regions[i, ]
    counter <- 0
    pct_old <- 0
    for (j in 1:num_points) {
      pct <- 100 * j / num_points
      point <- co_krig[j, ]
      if (is.null(rgeos::gDifference(point, region, byid = T))) {
        rate_avg[i] <- rate_avg[i] + co_krig@data$var1.pred[j]
        # shape_avg[i] <- shape_avg[i] + co_krig@data$shape.pred[j]
        counter <- counter + 1
      }
      if (pct - pct_old > 5) {
        print(paste0(pct, "% of outer iter completed"))
        pct_old <- pct
      }
    }
    if (counter != 0) {
      rate_avg[i] <- rate_avg[i] / counter
      # shape_avg[i] <- shape_avg[i] / counter
    }
  }
  # avg_mat <- rbind(ln.scale_avg, shape_avg)
  return(rate_avg)
}

ptm <- proc.time()
rate_avg <- compute_means_rate(krig_grid_rate, ws_regs)
(timer_rate <- proc.time() - ptm)
rate_avg


### estimates:
exp(krig_regs$ln.scale.pred)
krig_regs$shape.pred
krig_regs_rate$var1.pred

krig_fit <- rbind(exp(krig_regs$ln.scale.pred),
                  krig_regs$shape.pred,
                  krig_regs_rate$var1.pred)
rownames(krig_fit) <- c("scale", "shape", "rate")
colnames(krig_fit) <- colnames(car_fit)
# saveRDS(krig_fit, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/krig_fit.rds")


# # Function to generate a list of varcov matrices for each region, given...
# #   ln.scale.var: the vector of log(scale) variances
# #   shape.var: the vector of shape variances
# #   cov.ln.scale.shape: the vector of covariances between log(scale) and shape
# make_varcov_mats <- function(ln.scale.var, shape.var, cov.ln.scale.shape){      # still returns log(scale)...
#   varcov_list <- list()
#   # reg1 <- reg2 <- reg3 <- matrix(data=NA, nrow=2, ncol=2)
#   for(i in 1:length(ln.scale.var)){
#     mat <- matrix(data=NA, nrow=2, ncol=2)
#     mat[1,1] <- ln.scale.var[i]
#     mat[2,2] <- shape.var[i]
#     mat[1,2] <- mat[2,1] <- cov.ln.scale.shape[i]
#     varcov_list[[i]] <- mat
#   }
#   return(varcov_list)
# }



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
    mat[1,2] <- mat[2,1] <- cov.ln.scale.shape[i] # still need to transform the covariance!
    varcov_list[[i]] <- mat
  }
  return(varcov_list)
}


# Function to get SE's from a list of varcov matrices (useful for kriging model)
give_se <- function(varcov_list){ 
  se_vec <- NULL
  for (i in 1:length(varcov_list)){
    se_vec <- cbind(se_vec, sqrt(diag(varcov_list[[i]])))
  }
  rownames(se_vec) <- c("scale", "shape")
  colnames(se_vec) <- c("Reg1", "Reg2", "Reg3")
  return(se_vec)
}
# krig_fit <- rbind(exp(krig_regs$ln.scale.pred),
#                   krig_regs$shape.pred,
#                   krig_regs_rate$var1.pred)
# rownames(krig_fit) <- c("scale", "shape", "rate")
# colnames(krig_fit) <- colnames(car_fit)

# cbind(se_vec, sqrt(diag(krig_varcov[[1]])))


# krig_varcov <- make_varcov_mats(krig_regs$ln.scale.var, krig_regs$shape.var, krig_regs$cov.ln.scale.shape)
# give_se(krig_varcov)
krig_varcov <- make_varcov_mats_krig(krig_regs$ln.scale.var, krig_regs$ln.scale.pred, krig_regs$shape.var, krig_regs$cov.ln.scale.shape)
give_se(krig_varcov)
get_se_consol(krig_varcov)


pars_to_rl(krig_fit, 25)/254
pars_to_rl(krig_fit, 100)/254
pars_to_rl(krig_fit, 500)/254


### Fitting Model 2 to other windows 

# Already formatted data using format_data function, but this also changes order by region
# dat_win_22 and dat_win_52

# Function to format data (a stations_sub data frame) to stations within the regions, given...
#   window_fit: the window list subset to the moving window of choice (e.g. input window[[82]] for the final 40-yr window)
#   reg_polygons: the SpatialPolygonsDataFrame for the regions of interest
format_data_sub <- function(window_fit, reg_polygons){
  # Accessing Model Fits and saving the parameter values (scale, shape, rate) for each station
  ws_scale <- NULL
  ws_shape <- NULL
  ws_rate  <- NULL
  j = 1
  for(i in stat_nos){
    fit <- window_fit[[i]] 
    if(is.na(fit)){
      ws_scale[j] <- ws_shape[j] <- ws_rate[j]  <- NA
    }else{
      ws_scale[j] <- fit$results$par[1]
      ws_shape[j] <- fit$results$par[2]
      ws_rate[j]  <- fit$rate
    }
    j = j+1
  }
  # add new location (long/lat) columns that we will transform to coordinates, as well as parameter values
  stations_df <- stations %>%
    dplyr::mutate(scale = ws_scale, shape = ws_shape, rate = ws_rate, long = LONGITUDE, lat = LATITUDE)
  
  stations <- project_data(stations_df)
  region_int <- sp::over(stations, reg_polygons)
  region_int$REGION <- as.factor(region_int$REGION)
  h <- model.matrix( ~ REGION - 1, data=region_int)
  inds <- as.integer(rownames(h))
  
  stations_sub_df <- stations_df %>%
    dplyr::filter(STAT_NO %in% inds) %>%
    dplyr::filter(!is.na(scale))
  
  # stations_sub <- project_data(stations_sub_df)
  
  return(stations_sub_df)
}

# test_format_data <- format_data(window[[82]], ws_regs)
# testing_format_data_82 <- format_data_sub(window[[82]], ws_regs)
# identical(testing_format_data_82, stations_sub_df)

stations_sub_22_df <- format_data_sub(window[[22]], ws_regs)
stations_sub_52_df <- format_data_sub(window[[52]], ws_regs)

stations_sub_22 <- project_data(stations_sub_22_df)
stations_sub_52 <- project_data(stations_sub_52_df)

# # Could also just use data formatted for CAR model - just reorders the data by region
# stations_sub_22 <- project_data(dat_win_22)
# test_df <- stations_sub_22@data
# sp::zerodist(stations_sub_22) # (11, 15, 22)  (13, 16)  (17, 24)  # same hing as before, but out of order
# stations_sub_52 <- project_data(dat_win_52)

# singular results because mulitple observations at same location
sp::zerodist(stations_sub_52) # (8, 15) - remove station 8
stations_sub_52_df <- stations_sub_52_df[-8, ]
stations_sub_52 <- project_data(stations_sub_52_df)

sp::zerodist(stations_sub_22) # (4, 12, 20)  (8, 14)  (15, 22)  # keep 12, 14, 15
# remove 4, 8, 20, 22 
stations_sub_22_df <- stations_sub_22_df[-c(4, 8, 20, 22), ]
stations_sub_22 <- project_data(stations_sub_22_df)

### Window 22
# Fitting cross-variogram to real data (log(scale) and shape):
g <- gstat(id = "ln.scale", formula = log(scale)~1, data = stations_sub_22)
g <- gstat(g, id = "shape", formula = shape~1, data = stations_sub_22)
vg <- variogram(g)
cv.fit_22 = fit.lmc(vg, g, vgm("Mat"))
# cv.fit_22 = fit.lmc(vg, g, vgm(NA, "Mat", nugget = 0))
# cv.fit_22 = fit.lmc(vg, g, vgm("Sph")) # new
plot(vg, model = cv.fit_22)
cv.fit_22

# cokriging
# cv.fit_22$set=list(nocheck=1)
krig_regs_22 <- predict(cv.fit_22, newdata = ws_regs)
krig_regs_22$ln.scale.pred
### estimates:
exp(krig_regs_22$ln.scale.pred)
krig_regs_22$shape.pred
sqrt(krig_regs_22$shape.var)

## Need to bring in rate parameter
# Fitting variogram to real data (rate)
g3 <- gstat(id = "rate", formula = rate~1, data = stations_sub_22)
vg3 <- variogram(g3)
var.fit_22 <- fit.variogram(vg3, vgm("Sph")) # singular model
# var.fit_22 <- fit.variogram(vg3, vgm("Mat")) # no convergence
var.fit_22
plot(vg3, model = var.fit_22, main = "rate")
# kriging
krig_regs_rate_22 <- krige(rate~1, stations_sub_22, ws_regs, model = var.fit_22)
krig_regs_rate_22$var1.pred  # estimate
sqrt(krig_regs_rate_22$var1.var) # se


### Window 52
# Fitting cross-variogram to real data (log(scale) and shape):
g <- gstat(id = "ln.scale", formula = log(scale)~1, data = stations_sub_52)
g <- gstat(g, id = "shape", formula = shape~1, data = stations_sub_52)
vg <- variogram(g)
cv.fit_52 = fit.lmc(vg, g, vgm("Mat"))
# cv.fit_52 = fit.lmc(vg, g, vgm(NA, "Mat", nugget = 0))
# cv.fit_52 = fit.lmc(vg, g, vgm("Sph")) # new
plot(vg, model = cv.fit_52)
cv.fit_52

# cokriging
cv.fit_52$set=list(nocheck=1)
krig_regs_52 <- predict(cv.fit_52, newdata = ws_regs)
### estimates:
exp(krig_regs_52$ln.scale.pred)
krig_regs_52$shape.pred
sqrt(krig_regs_52$shape.var)

## Need to bring in rate parameter
# Fitting variogram to real data (rate)
g3 <- gstat(id = "rate", formula = rate~1, data = stations_sub_52)
vg3 <- variogram(g3)
var.fit_52 <- fit.variogram(vg3, vgm("Sph"))
# var.fit_52 <- fit.variogram(vg3, vgm("Mat"))
var.fit_52
plot(vg3, model = var.fit_52, main = "rate")
# kriging
krig_regs_rate_52 <- krige(rate~1, stations_sub_52, ws_regs, model = var.fit_52)
krig_regs_rate_52$var1.pred  # estimate
sqrt(krig_regs_rate_52$var1.var) # se


# Function to get parameter estimate matrix from kriging model fits, given...
#   krig_regs: the co-krige fit object for ln.scale and shape
#   krig_regs_rate: the krige fit object for rate
get_par_krig <- function(krig_regs, krig_regs_rate){
  krig_fit <- rbind(exp(krig_regs$ln.scale.pred),
                    krig_regs$shape.pred,
                    krig_regs_rate$var1.pred)
  rownames(krig_fit) <- c("scale", "shape", "rate")
  colnames(krig_fit) <- c("Reg1", "Reg2", "Reg3")
  return(krig_fit)
}

par_krig_22 <- get_par_krig(krig_regs_22, krig_regs_rate_22)
par_krig_52 <- get_par_krig(krig_regs_52, krig_regs_rate_52)

pars_to_rl(krig_fit, 25)/254
pars_to_rl(krig_fit, 100)/254
pars_to_rl(krig_fit, 500)/254

pars_to_rl(par_krig_22, 25)/254
pars_to_rl(par_krig_22, 100)/254
pars_to_rl(par_krig_22, 500)/254

pars_to_rl(par_krig_52, 25)/254
pars_to_rl(par_krig_52, 100)/254
pars_to_rl(par_krig_52, 500)/254



krig_varcov_22 <- make_varcov_mats_krig(krig_regs_22$ln.scale.var, krig_regs_22$ln.scale.pred, krig_regs_22$shape.var, krig_regs_22$cov.ln.scale.shape)
give_se(krig_varcov_22)
get_se_consol(krig_varcov_22)

krig_varcov_52 <- make_varcov_mats_krig(krig_regs_52$ln.scale.var, krig_regs_52$ln.scale.pred, krig_regs_52$shape.var, krig_regs_52$cov.ln.scale.shape)
give_se(krig_varcov_52)
get_se_consol(krig_varcov_52)


# Additional Functions and testing on the CAR model -----------------------


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

car_varcov <- make_varcov_mats_car(car_test_ln.scale, car_test_shape)
car_varcov[[3]]

summary(car_test_ln.scale)

car_test_ln.scale$fit$coefficients # param_ests
sqrt(diag((car_test_ln.scale$fit$imat * car_test_ln.scale$fit$s2))) # standard errors
diag((car_test_ln.scale$fit$imat * car_test_ln.scale$fit$s2)) # variances

car_test_ln.scale$fit$coefficients[1]
sqrt(diag((car_test_ln.scale$fit$imat * car_test_ln.scale$fit$s2)))[1]

exp(car_test_ln.scale$fit$coefficients[1])
exp(car_test_ln.scale$fit$coefficients[1]) * sqrt(diag((car_test_ln.scale$fit$imat * car_test_ln.scale$fit$s2)))[1]
# > mean(y)
# [1] 10
# > sd(y)
# [1] 0.03
# > lm=mean(log(y))
# > ls=sd(log(y))
# > exp(lm)*ls
# [1] 0.0300104 

sd(car_test_ln.scale$fit$fitted.values[ind1])
mean(car_test_ln.scale$fit$fitted.values[ind2])
mean(car_test_ln.scale$fit$fitted.values[ind3])
car_test_ln.scale$fit$signal_trend
car_test_ln.scale$fit$signal_stochastic
sd(car_test_ln.scale$fit$signal_stochastic[ind1])
sd(car_test_ln.scale$fit$fitted.values[ind1])
sd(car_test_ln.scale$fit$residuals[ind1])
# sd(test_dat$scale[ind1])
# sd(log(test_dat$scale[ind1]))
mean(test_dat$scale)
mean(log(test_dat$scale))

# Function to get parameter estimate matrix from CAR model fits, given...
#   ln.scale.fit: the CAR fit object for ln.scale
#   shape.fit: the CAR fit object for shape
#   rate.fit: the CAR fit object for rate
get_par_car <- function(ln.scale.fit, shape.fit, rate.fit){
  car_fit <- rbind(exp(summary(ln.scale.fit)$fit$coef),
                   summary(shape.fit)$fit$coef,
                   summary(rate.fit)$fit$coef)
  rownames(car_fit) <- c("scale", "shape", "rate")
  return(car_fit)
}

# Function to get SE's from CAR model fits
get_se <- function(ln.scale.fit, shape.fit, rate.fit){  # still returns log(scale)...
  mat <- NULL
  ln.scale_vcov <- ln.scale.fit$fit$imat * ln.scale.fit$fit$s2
  mat <- rbind(mat, sqrt(diag(ln.scale_vcov)))
  shape_vcov <- shape.fit$fit$imat * shape.fit$fit$s2
  mat <- rbind(mat, sqrt(diag(shape_vcov)))
  rate_vcov <- rate.fit$fit$imat * rate.fit$fit$s2
  mat <- rbind(mat, sqrt(diag(rate_vcov)))
  rownames(mat) <- c("ln.scale", "shape", "rate")
  colnames(mat) <- c("Reg1", "Reg2", "Reg3")
  return(mat)
}

car_fit
car_varcov[[1]]
car_fit[3,1]

# pars_to_rl <- function(par_mat, return_period){
#   RL_vec <- NULL
#   for(reg in 1:3){
#     scale <- par_mat[1, reg]
#     shape <- par_mat[2, reg]
#     rate  <- par_mat[3, reg]
#     RL_vec[reg] <- extRemes::rlevd(period=return_period, type="GP", scale=scale, shape=shape, rate=rate, threshold=thresh)
#   }
#   return(RL_vec)
# }


# Function to calculate return level CIs given parameter estimates and varcov matrix...
# using code pulled from functions in extRemes package
rl_with_ci <- function(par_mat, varcov_list, return_period, type = "ci", alpha = 0.05){
  out_ci <- out_se <- NULL
  for(i in 1:length(varcov_list)){  # for Region i
    cov.theta <- varcov_list[[i]]
    scale <- par_mat[1,i]
    shape <- par_mat[2,i]
    rate  <- par_mat[3,i]
    p <- extRemes::rlevd(period=return_period, type="GP", scale=scale, shape=shape, rate=rate, threshold=thresh)
    # p <- rlevd(period = return.period, loc = loc, scale = scale, 
    #            shape = shape, threshold = x$threshold, type = mod, 
    #            npy = x$npy, rate = lam)
    # p gives the basic return level!
    
    z.alpha <- qnorm(alpha/2, lower.tail = FALSE)

    # grads <- rlgrad.fevd(x, period = return.period)
    lam <- rate
    n = 14610 # while working with 40-year windows
    npy = 365.25
      m <- return_period * npy
      mlam <- m * lam
      grads <- cbind(scale * m^(-shape) * lam^(-shape - 1), 
                   (shape)^(-1) * ((mlam)^(shape) - 1), 
                   -scale * (shape)^(-2) * ((mlam)^(shape) - 1) + (scale/shape) * (mlam)^(shape) * log(mlam))
      grads <- t(grads)
      cov.theta <- rbind(c(lam * (1 - lam)/n, 0, 0), 
                                cbind(0, cov.theta))
      var.theta <- t(grads) %*% cov.theta %*% grads
      
      # return(var.theta)
      # which.par = 1

    # now make it work for all 3 regions!

    if(type=="se"){  # if type="se", give return level and standard error of RL
      rl_se <- sqrt(var.theta)
      out_se <- rbind(out_se, c(p, rl_se))
    }else{          # if type="ci" (default), give return level with CI
      out <- c(p, p - z.alpha * sqrt(var.theta),
                p + z.alpha * sqrt(var.theta))
      out_ci <- rbind(out_ci, out)
    }
  }
  if(type=="se"){  # if type="se", give return level and standard error of RL
    rownames(out_se) <- c("Reg1", "Reg2", "Reg3")
    colnames(out_se) <- c(paste0(return_period, "-yr RL"), "SE")
    return(out_se)
  }else{          # if type="ci" (default), give return level with CI
    rownames(out_ci) <- c("Reg1", "Reg2", "Reg3")
    colnames(out_ci) <- c(paste0(return_period, "-yr RL"), "LB", "UB")
    return(out_ci)
  }
}

rl_with_ci(car_fit, car_varcov, 100, "se")/254
rl_with_ci(car_fit, car_varcov, 100)/254

rl_with_ci(car_fit, car_varcov, 25, "se")/254
rl_with_ci(car_fit, car_varcov, 100, "se")/254
rl_with_ci(car_fit, car_varcov, 500, "se")/254

car_varcov_22 <- make_varcov_mats_car(car_ln.scale_22, car_shape_22)  # se of scale 8.5, 4, 9.25
# test_not_log_22 <- spatialreg::spautolm(scale ~ -1 + Reg1 + Reg2 + Reg3,  # se of scale 9.2, 5, 11.5
#                                         data = dat_win_22, family="CAR", listw=mat2listw(1/D_22)) 
rl_with_ci(par_22, car_varcov_22, 25, "se")/254
rl_with_ci(par_22, car_varcov_22, 100, "se")/254
rl_with_ci(par_22, car_varcov_22, 500, "se")/254

car_varcov_52 <- make_varcov_mats_car(car_ln.scale_52, car_shape_52)  # se of scale 5.5, 3.2, 7.8
# test_not_log_52 <- spatialreg::spautolm(scale ~ -1 + Reg1 + Reg2 + Reg3,  # se of scale 5, 2.5, 5.9
#                                         data = dat_win_52, family="CAR", listw=mat2listw(1/D_52))
rl_with_ci(par_52, car_varcov_52, 25, "se")/254
rl_with_ci(par_52, car_varcov_52, 100, "se")/254
rl_with_ci(par_52, car_varcov_52, 500, "se")/254


# Testing if transformed varcov is on same scale as if modeled scale instead of log(scale) -- for CAR model --
test_not_log_scale <- spatialreg::spautolm(scale ~ -1 + Reg1 + Reg2 + Reg3, data = test_dat, family="CAR",
                                                                listw=mat2listw(1/D_all_jit_sym))
# param ests somewhat different, but SE (of scale) is on the same scale, around 3 or 4
summary(test_not_log_scale)
sqrt(car_varcov[[1]])
sqrt(car_varcov[[2]])
sqrt(car_varcov[[3]])
car_fit

# Testing if transformed varcov is on same scale as if modeled scale instead of log(scale) -- for krige model --
# Fitting cross-variogram to real data (log(scale) and shape):
gt <- gstat(id = "scale", formula = scale~1, data = stations_sub)
gt <- gstat(gt, id = "shape", formula = shape~1, data = stations_sub)
vgt <- variogram(gt)
cv.fit2 = fit.lmc(vgt, gt, vgm("Mat"))
plot(vgt, model = cv.fit2)
cv.fit2
# cokriging
cv.fit2$set=list(nocheck=1)
krig_regs_not_log <- predict(cv.fit2, newdata = ws_regs)
### estimates:
krig_regs_not_log$scale.pred
krig_regs_not_log$shape.pred
krig_regs_not_log$scale.var
krig_regs_not_log$shape.var
krig_regs_not_log$cov.scale.shape



# # Can get varcov matrix between the coefficient estimates... but need varcov between scale and shape
# coef_vcov <- car_shape_22$fit$imat * car_shape_22$fit$s2
# sqrt(diag(coef_vcov))


# Trying jointly:
cbind(TOT, AMI)
lm_combo <- lm(cbind(shape, log(scale)) ~ -1 + Reg1 + Reg2 + Reg3, data = dat_win_22)
# Doesn't work with CAR fn
# car_combo_22 <- spatialreg::spautolm(cbind(shape, scale) ~ -1 + Reg1 + Reg2 + Reg3, data = dat_win_22, family="CAR",
#                                      listw=mat2listw(1/D_22))

# Can get varcov matrix between the coefficient estimates... but need varcov between scale and shape
coef_vcov <- car_test_shape$fit$imat * car_test_shape$fit$s2
#               Reg1          Reg2          Reg3
# Reg1  2.669983e-04 -3.232997e-05 -4.191154e-06
# Reg2 -3.232997e-05  2.946188e-04 -2.505771e-05
# Reg3 -4.191154e-06 -2.505771e-05  1.661473e-04
sqrt(2.669983e-04)  # 0.01634008
sqrt(diag(coef_vcov))

car_shape_22$fdHess # Numerical Hessian-based variance-covariance matrix... 
# is the same as matrix above, but with lambda added
car_shape_22$lambda.se; sqrt(car_shape_22$fdHess[1,1]) # first variable is lambda
sqrt(diag(car_shape_22$fdHess))
# sqrt(car_shape_22$fdHess[2,2]); sqrt(car_shape_22$fdHess[3,3]); sqrt(car_shape_22$fdHess[4,4])

summary(car_test_shape)
# car_test_shape$fit$signal_trend
# car_test_shape$fit$signal_stochastic

identical(car_test_shape$fit$fitted.values, car_test_shape$fit$signal_trend + car_test_shape$fit$signal_stochastic)

ind1 <- which(test_dat$Reg1==1)
ind2 <- which(test_dat$Reg2==1)
ind3 <- which(test_dat$Reg3==1)
sqrt(var(car_test_shape$fit$fitted.values[ind1]))
sqrt(var(car_test_shape$fit$signal_trend[ind1]))
sqrt(var(car_test_shape$fit$signal_stochastic[ind1]))
sqrt(var(car_test_shape$fit$residuals[ind1]))

car_test_shape$fit$fitted.values
car_test_ln.scale$fit$fitted.values
var(car_test_shape$fit$fitted.values)
var(car_test_ln.scale$fit$fitted.values)
cov(car_test_shape$fit$fitted.values,
    car_test_ln.scale$fit$fitted.values)
plot(car_test_shape$fit$fitted.values, car_test_ln.scale$fit$fitted.values)

# plot(car_test_shape$fit$fitted.values, car_test_shape$fit$residuals)
# plot(car_test_ln.scale$fit$fitted.values, car_test_ln.scale$fit$residuals)
plot(car_test_shape$fit$residuals, car_test_ln.scale$fit$residuals)

# get_se(car_ln.scale_22, car_shape_22, car_rate_22)
get_se(car_test_ln.scale, car_test_shape, car_test_rate)
get_par_car(car_test_ln.scale, car_test_shape, car_test_rate) # check, same as car_fit - yes!
# car_fit 


library(matrixcalc)
is.positive.definite(1/D_all_jit_sym)



# Moving Windows ----------------------------------------------------------

# dat_win_22 <- format_data(window[[22]], ws_regs) # 34
# dat_win_52 <- format_data(window[[52]], ws_regs) # 28
# dat_win_82 <- format_data(window[[82]], ws_regs) # 149
# 
# D_22 <- get_sym_car_mat(dat_win_22, hMat_miles)
# D_52 <- get_sym_car_mat(dat_win_52, hMat_miles)
# D_82 <- get_sym_car_mat(dat_win_82, hMat_miles)

car_shape_22 <- spatialreg::spautolm(shape ~ -1 + Reg1 + Reg2 + Reg3, data = dat_win_22, family="CAR",
                                     listw=mat2listw(1/D_22))
ptm <- proc.time()
car_ln.scale_22 <- spatialreg::spautolm(log(scale) ~ -1 + Reg1 + Reg2 + Reg3, data = dat_win_22, family="CAR",
                                        listw=mat2listw(1/D_22))
proc.time() - ptm
ptm <- proc.time()
car_rate_22 <- spatialreg::spautolm(rate ~ -1 + Reg1 + Reg2 + Reg3, data = dat_win_22, family="CAR",
                                    listw=mat2listw(1/D_22))
proc.time() - ptm

par_22 <- get_par_car(car_ln.scale_22, car_shape_22, car_rate_22)
se_22 <- get_se(car_ln.scale_22, car_shape_22, car_rate_22)


ptm <- proc.time()
car_ln.scale_52 <- spatialreg::spautolm(log(scale) ~ -1 + Reg1 + Reg2 + Reg3, data = dat_win_52, family="CAR",
                                        listw=mat2listw(1/D_52))
proc.time() - ptm
ptm <- proc.time()
car_shape_52 <- spatialreg::spautolm(shape ~ -1 + Reg1 + Reg2 + Reg3, data = dat_win_52, family="CAR",
                                     listw=mat2listw(1/D_52))
proc.time() - ptm
ptm <- proc.time()
car_rate_52 <- spatialreg::spautolm(rate ~ -1 + Reg1 + Reg2 + Reg3, data = dat_win_52, family="CAR",
                                    listw=mat2listw(1/D_52))
proc.time() - ptm

par_52 <- get_par_car(car_ln.scale_52, car_shape_52, car_rate_52)
se_52 <- get_se(car_ln.scale_52, car_shape_52, car_rate_52)


## Window 82 is already found through car_fit! (car_test_shape...)
get_se(car_test_ln.scale, car_test_shape, car_test_rate)
get_par_car(car_test_ln.scale, car_test_shape, car_test_rate) # check! same as car_fit
# car_fit 

pars_to_rl(par_22, 25)/254
pars_to_rl(par_22, 100)/254
pars_to_rl(par_22, 500)/254

pars_to_rl(par_52, 25)/254
pars_to_rl(par_52, 100)/254
pars_to_rl(par_52, 500)/254

pars_to_rl(car_fit, 25)/254
pars_to_rl(car_fit, 100)/254
pars_to_rl(car_fit, 500)/254



# Model 3 - Consolidated Series -------------------------------------------------------------------

# load consolidated series and subset to same 40-year window (1981-2020)
setwd("~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data")

thresh <- 253

# Was it declustered when created? NO, so we decluster it here
precip <- read.csv("max_precip_regions.csv")
precip <- precip[,-1]   # delete first column since is repeat of row names (numbers)

### Declustering
precip_dclust = precip
for ( i in seq(from=2, to=length(precip[1,]))){
  station = precip[,i]
  non_na = which(station>=0)
  station[is.na(station)] <- 0
  
  dec = extRemes::decluster(station, threshold=0, clusterfun = "max", replace.with = 0)
  stat_clus = as.numeric(dec)
  
  station = stat_clus
  precip_dclust[non_na,i] = station[non_na]
}

### Save as data and format Date column
consol_data <- precip_dclust
consol_data$Date <- lubridate::ymd(as.character(consol_data$Date)) #put Date column into date format

# # # # # # # # # # # # # # #

# 40-year period
start <- "1981-01-01"
end <- "2020-12-31"

start <- which(consol_data$Date==start)  #finding indexes corresponding to start & end dates
end   <- which(consol_data$Date==end)
sub <- consol_data[start:end, ]  #subset data to those 40 years
fitreg1 <- extRemes::fevd(sub$region1, threshold=thresh, type="GP", method="MLE")
fitreg2 <- extRemes::fevd(sub$region2, threshold=thresh, type="GP", method="MLE")
fitreg3 <- extRemes::fevd(sub$region3, threshold=thresh, type="GP", method="MLE")
fitreg1$rate; fitreg2$rate; fitreg3$rate
mean(c(fitreg1$rate, fitreg2$rate, fitreg3$rate))

fitreg1$results$par[1]
fitreg1$results$par[2]
fitreg1$rate
consol_fit <- rbind(c(fitreg1$results$par[1], fitreg2$results$par[1], fitreg3$results$par[1]), # scale
                    c(fitreg1$results$par[2], fitreg2$results$par[2], fitreg3$results$par[2]), # shape
                    c(fitreg1$rate, fitreg2$rate, fitreg3$rate)) # rate
rownames(consol_fit) <- c("scale", "shape", "rate")
colnames(consol_fit) <- colnames(car_fit)
# saveRDS(consol_fit, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/consol_fit.rds")

pars_to_rl(consol_fit, 25)/254
pars_to_rl(consol_fit, 100)/254
pars_to_rl(consol_fit, 500)/254



# Function to get parameter estimate matrix from consolidated data fits, given...
#   fitreg1: GPD fit object for Region 1
#   fitreg2: GPD fit object for Region 2
#   fitreg3: GPD fit object for Region 3
get_par_consol <- function(fitreg1, fitreg2, fitreg3){
  consol_fit <- rbind(c(fitreg1$results$par[1], fitreg2$results$par[1], fitreg3$results$par[1]), # scale
                      c(fitreg1$results$par[2], fitreg2$results$par[2], fitreg3$results$par[2]), # shape
                      c(fitreg1$rate, fitreg2$rate, fitreg3$rate)) # rate
  rownames(consol_fit) <- c("scale", "shape", "rate")
  colnames(consol_fit) <- colnames(car_fit)
  return(consol_fit)
}


# Function to generate a list of varcov matrices for each region, given...
#   fitreg1: GPD fit object for Region 1
#   fitreg2: GPD fit object for Region 2
#   fitreg3: GPD fit object for Region 3
make_varcov_mats_consol <- function(fitreg1, fitreg2, fitreg3){
  varcov_list <- list()
  varcov_list[[1]] <- extRemes::parcov.fevd(fitreg1)
  varcov_list[[2]] <- extRemes::parcov.fevd(fitreg2)
  varcov_list[[3]] <- extRemes::parcov.fevd(fitreg3)
  return(varcov_list)
}

# Function to get SE's from a list of varcov matrices (useful for model 3)
get_se_consol <- function(varcov_list){
  se_vec <- NULL
  for (i in 1:length(varcov_list)){
    se_vec <- cbind(se_vec, sqrt(diag(varcov_list[[i]])))
  }
  rownames(se_vec) <- c("scale", "shape")
  colnames(se_vec) <- c("Reg1", "Reg2", "Reg3")
  return(se_vec)
}

consol_varcov <- make_varcov_mats_consol(fitreg1, fitreg2, fitreg3)
rl_with_ci(consol_fit, consol_varcov, 25, "se")/254
rl_with_ci(consol_fit, consol_varcov, 100, "se")/254
rl_with_ci(consol_fit, consol_varcov, 500, "se")/254

# See if matches built-in return.level function - YES!
extRemes::return.level(fitreg1, return.period = c(25, 100, 500), do.ci=T)/254

# Need to fit data for this - subset to the years of the windows
# 22: 1921-1960
# 52: 1951-1990
#last 1981-2020

# 40-year period
start <- "1921-01-01"
end <- "1960-12-31"
start <- which(consol_data$Date==start)  #finding indexes corresponding to start & end dates
end   <- which(consol_data$Date==end)
sub <- consol_data[start:end, ]  #subset data to those 40 years
# there are NAs... need to remove
# fitreg1_22_comp <- extRemes::fevd(na.omit(sub$region1), threshold=thresh, type="GP", method="MLE")
fitreg1_22 <- extRemes::fevd(sub$region1, threshold=thresh, type="GP", method="MLE", na.action = na.exclude)
fitreg2_22 <- extRemes::fevd(sub$region2, threshold=thresh, type="GP", method="MLE")
fitreg3_22 <- extRemes::fevd(sub$region3, threshold=thresh, type="GP", method="MLE", na.action = na.exclude)

par_consol_22 <- get_par_consol(fitreg1_22, fitreg2_22, fitreg3_22)
consol_varcov_22 <- make_varcov_mats_consol(fitreg1_22, fitreg2_22, fitreg3_22)
rl_with_ci(par_consol_22, consol_varcov_22, 25, "se")/254
rl_with_ci(par_consol_22, consol_varcov_22, 100, "se")/254
rl_with_ci(par_consol_22, consol_varcov_22, 500, "se")/254

# 40-year period
start <- "1951-01-01"
end <- "1990-12-31"
start <- which(consol_data$Date==start)  #finding indexes corresponding to start & end dates
end   <- which(consol_data$Date==end)
sub <- consol_data[start:end, ]  #subset data to those 40 years
fitreg1_52 <- extRemes::fevd(sub$region1, threshold=thresh, type="GP", method="MLE")
fitreg2_52 <- extRemes::fevd(sub$region2, threshold=thresh, type="GP", method="MLE")
fitreg3_52 <- extRemes::fevd(sub$region3, threshold=thresh, type="GP", method="MLE")

par_consol_52 <- get_par_consol(fitreg1_52, fitreg2_52, fitreg3_52)
consol_varcov_52 <- make_varcov_mats_consol(fitreg1_52, fitreg2_52, fitreg3_52)
rl_with_ci(par_consol_52, consol_varcov_52, 25, "se")/254
rl_with_ci(par_consol_52, consol_varcov_52, 100, "se")/254
rl_with_ci(par_consol_52, consol_varcov_52, 500, "se")/254

### To put in tables
consol_fit
get_se_consol(consol_varcov)
par_consol_22
get_se_consol(consol_varcov_22)
par_consol_52
get_se_consol(consol_varcov_52)


# Plot Moving Window Results ----------------------------------------------
library(ggplot2)

# test_25 <- rl_with_ci(par_consol_52, consol_varcov_52, 25, "ci")/254
# test_100 <- rl_with_ci(par_consol_52, consol_varcov_52, 100, "ci")/254
# test_500 <- rl_with_ci(par_consol_52, consol_varcov_52, 500, "ci")/254
# dim(test_500)[1]
# test_500[1,]

# Function to convert return level CI matrices into vectors for easy plotting

# separate by region! For each region, will have matrix made of diff RL values (25, 100, 500 and bounds)
rl_ci_to_plot_vec <- function(rl_mat_25, rl_mat_100, rl_mat_500){
  reg <- list()
  for ( i in 1:dim(rl_mat_25)[1]){ #i.e. 1:3
    vec_25  <- rl_mat_25[i, ]
    vec_100 <- rl_mat_100[i, ]
    vec_500 <- rl_mat_500[i, ]
    mat <- rbind(vec_25, vec_100, vec_500)
    rownames(mat) <- c("25-yr", "100-yr", "500-yr")
    colnames(mat)[1] <- "RL"
    reg[[i]] <- mat
  }
  return(reg)
}

# win_2 <- rl_ci_to_plot_vec(test_25, test_100, test_500)
# win_2[[3]]


# Function to create dataset for plotting, given...
#   win_1, win_2, win_3: the list of RLs and CIs by region for each window
# currently only for 3 windows
win_rl_to_plot_dat <- function(win_1, win_2, win_3){
  reg <- list()
  for(i in 1:3){ # by region
    # win_1[[i]] # is region i values from window 1
    dat_25 <- rbind(win_1[[i]][1, ], win_2[[i]][1, ], win_3[[i]][1, ])
    dat_100 <- rbind(win_1[[i]][2, ], win_2[[i]][2, ], win_3[[i]][2, ])
    dat_500 <- rbind(win_1[[i]][3, ], win_2[[i]][3, ], win_3[[i]][3, ])
    dat <- cbind(dat_25, dat_100, dat_500)
    colnames(dat) <- c("RL_25", "LB_25", "UB_25", "RL_100", "LB_100", "UB_100","RL_500", "LB_500", "UB_500")
    reg[[i]] <- as.data.frame(dat)
  }
  return(reg)
}


# Function to create dataset for plotting, given...  (useful for RLs with no CI)
#   par_krig_22, par_krig_52, krig_fit: the parameter matrices for each window
#  Note: this does not have CIs, just RLs
pars_to_plot_dat <- function(par_krig_22, par_krig_52, krig_fit){
  mat_25 <- rbind(pars_to_rl(par_krig_22, 25)/254,
                  pars_to_rl(par_krig_52, 25)/254,
                  pars_to_rl(krig_fit, 25)/254)
  mat_100 <-rbind(pars_to_rl(par_krig_22, 100)/254,
                  pars_to_rl(par_krig_52, 100)/254,
                  pars_to_rl(krig_fit, 100)/254)
  mat_500 <-rbind(pars_to_rl(par_krig_22, 500)/254,
                  pars_to_rl(par_krig_52, 500)/254,
                  pars_to_rl(krig_fit, 500)/254)
  colnames(mat_25) <- colnames(mat_100) <- colnames(mat_500) <- c("Reg1", "Reg2", "Reg3")
  reg <- list()
  for(i in 1:dim(mat_25)[1]){ #i.e. 1:3
    dat_25 <- mat_25[, i]
    dat_100 <- mat_100[, i]
    dat_500 <- mat_500[, i]
    dat <- cbind(dat_25, dat_100, dat_500)
    colnames(dat) <- c("RL_25", "RL_100", "RL_500")
    reg[[i]] <- as.data.frame(dat)
  }
  return(reg)
}




### Fitting to Model 1 ###
car_rl_25_22  <- rl_with_ci(par_22, car_varcov_22, 25, "ci")/254
car_rl_100_22 <- rl_with_ci(par_22, car_varcov_22, 100, "ci")/254
car_rl_500_22 <- rl_with_ci(par_22, car_varcov_22, 500, "ci")/254

car_rl_25_52  <- rl_with_ci(par_52, car_varcov_52, 25, "ci")/254
car_rl_100_52 <- rl_with_ci(par_52, car_varcov_52, 100, "ci")/254
car_rl_500_52 <- rl_with_ci(par_52, car_varcov_52, 500, "ci")/254

car_rl_25_82  <- rl_with_ci(car_fit, car_varcov, 25, "ci")/254
car_rl_100_82 <- rl_with_ci(car_fit, car_varcov, 100, "ci")/254
car_rl_500_82 <- rl_with_ci(car_fit, car_varcov, 500, "ci")/254

win_1 <- rl_ci_to_plot_vec(car_rl_25_22, car_rl_100_22, car_rl_500_22)
win_2 <- rl_ci_to_plot_vec(car_rl_25_52, car_rl_100_52, car_rl_500_52)
win_3 <- rl_ci_to_plot_vec(car_rl_25_82, car_rl_100_82, car_rl_500_82)

plot_dat <- win_rl_to_plot_dat(win_1, win_2, win_3)
reg1 <- plot_dat[[1]] # Region 1 data frame for plotting
reg2 <- plot_dat[[2]] 
reg3 <- plot_dat[[3]] 

## Plotting
# Fixed the colors from old Richmond code
# started with ylim = c(5,30)... maybe c(5,33)?
ggplot(data=reg1, aes(x=c(1960, 1990, 2020))) + 
  coord_cartesian(ylim = c(5, 35)) +
  geom_line(aes(y=RL_500, color="500-Year")) + geom_point(aes(y=RL_500, color="500-Year")) + geom_ribbon(aes(ymin=LB_500, ymax=UB_500), linetype=2, alpha=0.07, fill="red") +
  geom_line(aes(y=RL_100, color="100-Year")) + geom_point(aes(y=RL_100, color="100-Year"))+ geom_ribbon(aes(ymin=LB_100, ymax=UB_100), linetype=2, alpha=0.07, fill="blue") +
  geom_line(aes(y=RL_25, color="25-Year")) + geom_point(aes(y=RL_25, color="25-Year")) + geom_ribbon(aes(ymin=LB_25, ymax=UB_25), linetype=2, alpha=0.15, fill="green") +
  labs(x="Last Year of 40-Year Window", y="Return Level (in)", title="Estimated Return Levels - Model 1 - Region 1") +
  scale_colour_manual(name = "Return Period", values = c('500-Year' = "red", '100-Year' = "blue", '25-Year' = "green"))

ggplot(data=reg2, aes(x=c(1960, 1990, 2020))) + 
  coord_cartesian(ylim = c(5, 35)) +
  geom_line(aes(y=RL_500, color="500-Year")) + geom_point(aes(y=RL_500, color="500-Year")) + geom_ribbon(aes(ymin=LB_500, ymax=UB_500), linetype=2, alpha=0.07, fill="red") +
  geom_line(aes(y=RL_100, color="100-Year")) + geom_point(aes(y=RL_100, color="100-Year"))+ geom_ribbon(aes(ymin=LB_100, ymax=UB_100), linetype=2, alpha=0.07, fill="blue") +
  geom_line(aes(y=RL_25, color="25-Year")) + geom_point(aes(y=RL_25, color="25-Year")) + geom_ribbon(aes(ymin=LB_25, ymax=UB_25), linetype=2, alpha=0.15, fill="green") +
  labs(x="Last Year of 40-Year Window", y="Return Level (in)", title="Estimated Return Levels - Model 1 - Region 2") +
  scale_colour_manual(name = "Return Period", values = c('500-Year' = "red", '100-Year' = "blue", '25-Year' = "green"))

ggplot(data=reg3, aes(x=c(1960, 1990, 2020))) + 
  coord_cartesian(ylim = c(5, 35)) +
  geom_line(aes(y=RL_500, color="500-Year")) + geom_point(aes(y=RL_500, color="500-Year")) + geom_ribbon(aes(ymin=LB_500, ymax=UB_500), linetype=2, alpha=0.07, fill="red") +
  geom_line(aes(y=RL_100, color="100-Year")) + geom_point(aes(y=RL_100, color="100-Year"))+ geom_ribbon(aes(ymin=LB_100, ymax=UB_100), linetype=2, alpha=0.07, fill="blue") +
  geom_line(aes(y=RL_25, color="25-Year")) + geom_point(aes(y=RL_25, color="25-Year")) + geom_ribbon(aes(ymin=LB_25, ymax=UB_25), linetype=2, alpha=0.15, fill="green") +
  labs(x="Last Year of 40-Year Window", y="Return Level (in)", title="Estimated Return Levels - Model 1 - Region 3") +
  scale_colour_manual(name = "Return Period", values = c('500-Year' = "red", '100-Year' = "blue", '25-Year' = "green"))


### Fitting to Model 2 ###
krig_plot_dat <- pars_to_plot_dat(par_krig_22, par_krig_52, krig_fit)
reg1_krig <- krig_plot_dat[[1]] # Region 1 data frame for plotting
reg2_krig <- krig_plot_dat[[2]] 
reg3_krig <- krig_plot_dat[[3]] 

## Plotting - no CIs
ggplot(data=reg1_krig, aes(x=c(1960, 1990, 2020))) + 
  coord_cartesian(ylim = c(5, 35)) +
  geom_line(aes(y=RL_500, color="500-Year")) + geom_point(aes(y=RL_500, color="500-Year")) +
  geom_line(aes(y=RL_100, color="100-Year")) + geom_point(aes(y=RL_100, color="100-Year"))+ 
  geom_line(aes(y=RL_25, color="25-Year")) + geom_point(aes(y=RL_25, color="25-Year")) + 
  labs(x="Last Year of 40-Year Window", y="Return Level (in)", title="Estimated Return Levels - Model 2 - Region 1") +
  scale_colour_manual(name = "Return Period", values = c('500-Year' = "red", '100-Year' = "blue", '25-Year' = "green"))

ggplot(data=reg2_krig, aes(x=c(1960, 1990, 2020))) + 
  coord_cartesian(ylim = c(5, 35)) +
  geom_line(aes(y=RL_500, color="500-Year")) + geom_point(aes(y=RL_500, color="500-Year")) +
  geom_line(aes(y=RL_100, color="100-Year")) + geom_point(aes(y=RL_100, color="100-Year"))+ 
  geom_line(aes(y=RL_25, color="25-Year")) + geom_point(aes(y=RL_25, color="25-Year")) + 
  labs(x="Last Year of 40-Year Window", y="Return Level (in)", title="Estimated Return Levels - Model 2 - Region 2") +
  scale_colour_manual(name = "Return Period", values = c('500-Year' = "red", '100-Year' = "blue", '25-Year' = "green"))

ggplot(data=reg3_krig, aes(x=c(1960, 1990, 2020))) + 
  coord_cartesian(ylim = c(5, 35)) +
  geom_line(aes(y=RL_500, color="500-Year")) + geom_point(aes(y=RL_500, color="500-Year")) +
  geom_line(aes(y=RL_100, color="100-Year")) + geom_point(aes(y=RL_100, color="100-Year"))+ 
  geom_line(aes(y=RL_25, color="25-Year")) + geom_point(aes(y=RL_25, color="25-Year")) + 
  labs(x="Last Year of 40-Year Window", y="Return Level (in)", title="Estimated Return Levels - Model 2 - Region 3") +
  scale_colour_manual(name = "Return Period", values = c('500-Year' = "red", '100-Year' = "blue", '25-Year' = "green"))


### Fitting to Model 3 ###
consol_rl_25_22  <- rl_with_ci(par_consol_22, consol_varcov_22, 25, "ci")/254
consol_rl_100_22 <- rl_with_ci(par_consol_22, consol_varcov_22, 100, "ci")/254
consol_rl_500_22 <- rl_with_ci(par_consol_22, consol_varcov_22, 500, "ci")/254

consol_rl_25_52  <- rl_with_ci(par_consol_52, consol_varcov_52, 25, "ci")/254
consol_rl_100_52 <- rl_with_ci(par_consol_52, consol_varcov_52, 100, "ci")/254
consol_rl_500_52 <- rl_with_ci(par_consol_52, consol_varcov_52, 500, "ci")/254

consol_rl_25_82  <- rl_with_ci(consol_fit, consol_varcov, 25, "ci")/254
consol_rl_100_82 <- rl_with_ci(consol_fit, consol_varcov, 100, "ci")/254
consol_rl_500_82 <- rl_with_ci(consol_fit, consol_varcov, 500, "ci")/254

consol_win_1 <- rl_ci_to_plot_vec(consol_rl_25_22, consol_rl_100_22, consol_rl_500_22)
consol_win_2 <- rl_ci_to_plot_vec(consol_rl_25_52, consol_rl_100_52, consol_rl_500_52)
consol_win_3 <- rl_ci_to_plot_vec(consol_rl_25_82, consol_rl_100_82, consol_rl_500_82)

consol_plot_dat <- win_rl_to_plot_dat(consol_win_1, consol_win_2, consol_win_3)
reg1_consol <- consol_plot_dat[[1]] # Region 1 data frame for plotting
reg2_consol <- consol_plot_dat[[2]] 
reg3_consol <- consol_plot_dat[[3]] 

## Plotting
ggplot(data=reg1_consol, aes(x=c(1960, 1990, 2020))) + 
  coord_cartesian(ylim = c(5, 35)) +
  geom_line(aes(y=RL_500, color="500-Year")) + geom_point(aes(y=RL_500, color="500-Year")) + geom_ribbon(aes(ymin=LB_500, ymax=UB_500), linetype=2, alpha=0.07, fill="red") +
  geom_line(aes(y=RL_100, color="100-Year")) + geom_point(aes(y=RL_100, color="100-Year"))+ geom_ribbon(aes(ymin=LB_100, ymax=UB_100), linetype=2, alpha=0.07, fill="blue") +
  geom_line(aes(y=RL_25, color="25-Year")) + geom_point(aes(y=RL_25, color="25-Year")) + geom_ribbon(aes(ymin=LB_25, ymax=UB_25), linetype=2, alpha=0.15, fill="green") +
  labs(x="Last Year of 40-Year Window", y="Return Level (in)", title="Estimated Return Levels - Model 3 - Region 1") +
  scale_colour_manual(name = "Return Period", values = c('500-Year' = "red", '100-Year' = "blue", '25-Year' = "green"))

ggplot(data=reg2_consol, aes(x=c(1960, 1990, 2020))) + 
  coord_cartesian(ylim = c(5, 35)) +
  geom_line(aes(y=RL_500, color="500-Year")) + geom_point(aes(y=RL_500, color="500-Year")) + geom_ribbon(aes(ymin=LB_500, ymax=UB_500), linetype=2, alpha=0.07, fill="red") +
  geom_line(aes(y=RL_100, color="100-Year")) + geom_point(aes(y=RL_100, color="100-Year"))+ geom_ribbon(aes(ymin=LB_100, ymax=UB_100), linetype=2, alpha=0.07, fill="blue") +
  geom_line(aes(y=RL_25, color="25-Year")) + geom_point(aes(y=RL_25, color="25-Year")) + geom_ribbon(aes(ymin=LB_25, ymax=UB_25), linetype=2, alpha=0.15, fill="green") +
  labs(x="Last Year of 40-Year Window", y="Return Level (in)", title="Estimated Return Levels - Model 3 - Region 2") +
  scale_colour_manual(name = "Return Period", values = c('500-Year' = "red", '100-Year' = "blue", '25-Year' = "green"))

ggplot(data=reg3_consol, aes(x=c(1960, 1990, 2020))) + 
  coord_cartesian(ylim = c(5, 35)) +
  geom_line(aes(y=RL_500, color="500-Year")) + geom_point(aes(y=RL_500, color="500-Year")) + geom_ribbon(aes(ymin=LB_500, ymax=UB_500), linetype=2, alpha=0.07, fill="red") +
  geom_line(aes(y=RL_100, color="100-Year")) + geom_point(aes(y=RL_100, color="100-Year"))+ geom_ribbon(aes(ymin=LB_100, ymax=UB_100), linetype=2, alpha=0.07, fill="blue") +
  geom_line(aes(y=RL_25, color="25-Year")) + geom_point(aes(y=RL_25, color="25-Year")) + geom_ribbon(aes(ymin=LB_25, ymax=UB_25), linetype=2, alpha=0.15, fill="green") +
  labs(x="Last Year of 40-Year Window", y="Return Level (in)", title="Estimated Return Levels - Model 3 - Region 3") +
  scale_colour_manual(name = "Return Period", values = c('500-Year' = "red", '100-Year' = "blue", '25-Year' = "green"))


### Fitting to Model 2 ###
# consol_rl_25_22  <- rl_with_ci(par_consol_22, consol_varcov_22, 25, "ci")/254
# consol_rl_100_22 <- rl_with_ci(par_consol_22, consol_varcov_22, 100, "ci")/254
# consol_rl_500_22 <- rl_with_ci(par_consol_22, consol_varcov_22, 500, "ci")/254
# 
# consol_rl_25_52  <- rl_with_ci(par_consol_52, consol_varcov_52, 25, "ci")/254
# consol_rl_100_52 <- rl_with_ci(par_consol_52, consol_varcov_52, 100, "ci")/254
# consol_rl_500_52 <- rl_with_ci(par_consol_52, consol_varcov_52, 500, "ci")/254
# 
# consol_rl_25_82  <- rl_with_ci(consol_fit, consol_varcov, 25, "ci")/254
# consol_rl_100_82 <- rl_with_ci(consol_fit, consol_varcov, 100, "ci")/254
# consol_rl_500_82 <- rl_with_ci(consol_fit, consol_varcov, 500, "ci")/254




# Saving Model Fits -------------------------------------------------------------------------------

# saveRDS(car_test_shape, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/car_test_shape.rds")
# saveRDS(car_test_ln.scale, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/car_test_ln.scale.rds")
# saveRDS(car_test_rate, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/car_test_rate.rds")
# saveRDS(car_fit, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/car_fit.rds")
# 
# saveRDS(krig_regs, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/krig_regs.rds")
# saveRDS(krig_regs_rate, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/krig_regs_rate.rds")
# saveRDS(krig_grid, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/krig_grid.rds")
# saveRDS(krig_grid_rate, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/krig_grid_rate.rds")
# saveRDS(krig_fit, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/krig_fit.rds")
# 
# saveRDS(fitreg1, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/consol_fitreg1.rds")
# saveRDS(fitreg2, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/consol_fitreg2.rds")
# saveRDS(fitreg3, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/consol_fitreg3.rds")
# saveRDS(consol_fit, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/consol_fit.rds")


# saveRDS(car_varcov, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/car_varcov.rds")
# saveRDS(D_all_jit_sym, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/D_82.rds")
# saveRDS(D_all_jit_sym, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/D_all_jit_sym.rds")
# 
# saveRDS(car_shape_22, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/car_shape_22.rds")
# saveRDS(car_ln.scale_22, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/car_ln.scale_22.rds")
# saveRDS(car_rate_22, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/car_rate_22.rds")
# saveRDS(par_22, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/car_par_22.rds")
# # se not quite right yet
# saveRDS(car_varcov_22, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/car_varcov_22.rds")
# saveRDS(D_22, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/D_22.rds")
# 
# saveRDS(car_shape_52, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/car_shape_52.rds")
# saveRDS(car_ln.scale_52, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/car_ln.scale_52.rds")
# saveRDS(car_rate_52, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/car_rate_52.rds")
# saveRDS(par_52, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/car_par_52.rds")
# # se not quite right yet
# saveRDS(car_varcov_52, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/car_varcov_52.rds")
# saveRDS(D_52, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/D_52.rds")

# saveRDS(consol_data, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/consol_data.rds")
# saveRDS(consol_varcov, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/consol_varcov.rds")

# saveRDS(fitreg1_22, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/consol_fitreg1_22.rds")
# saveRDS(fitreg2_22, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/consol_fitreg2_22.rds")
# saveRDS(fitreg3_22, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/consol_fitreg3_22.rds")
# saveRDS(par_consol_22, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/par_consol_22.rds")
# saveRDS(consol_varcov_22, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/consol_varcov_22.rds")
# 
# saveRDS(fitreg1_52, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/consol_fitreg1_52.rds")
# saveRDS(fitreg2_52, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/consol_fitreg2_52.rds")
# saveRDS(fitreg3_52, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/consol_fitreg3_52.rds")
# saveRDS(par_consol_52, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/par_consol_52.rds")
# saveRDS(consol_varcov_52, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/consol_varcov_52.rds")

### Fits from kriging model
# saveRDS(krig_regs_22, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/krig_regs_22.rds")
# saveRDS(krig_regs_rate_22, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/krig_regs_rate_22.rds")
# saveRDS(krig_regs_52, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/krig_regs_52.rds")
# saveRDS(krig_regs_rate_52, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/krig_regs_rate_52.rds")
# saveRDS(par_krig_22, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/par_krig_22.rds")
# saveRDS(par_krig_52, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/Saved_Fit_Objects/par_krig_52.rds")




# Extra Code - Finding CIs ------------------------------------------------------------------------

fit <- window[[79]][[588]]

# confidence intervals
extRemes::return.level(fit, return.period = 100, do.ci = T)  # same as  extRemes::ci.fevd(fit, type="return.level", return.period = 100)
extRemes::ci.fevd(fit, type="parameter")


z = qnorm(0.975) # for 95% CI (default)

# saving estimates
scale_est <- fit$results$par[1]
shape_est <- fit$results$par[2]

#saving standard errors
cov_mat <- extRemes::parcov.fevd(fit)
scale_se <- sqrt(extRemes::parcov.fevd(fit)[1,1])
shape_se <- sqrt(extRemes::parcov.fevd(fit)[2,2])

# Hand-calculated CIs, matches ci.fevd function with type="parameter"
# scale
scale_est - z*scale_se;   scale_est;    scale_est + z*scale_se
# shape
shape_est - z*shape_se;   shape_est;    shape_est + z*shape_se


# Return levels:
# from extRemes::rlevd
calc_rl <- function(threshold, scale, shape, rate, period){
  npy = 365.25
  m <- period * npy * rate
  if (shape == 0) 
    res <- threshold + scale * log(m)
  else res <- threshold + (scale/shape) * (m^shape - 1)
  return(res)
}

# Reference CI:
extRemes::return.level(fit, return.period = 100, do.ci = T)

calc_rl(threshold = 253, scale=scale_est, shape=shape_est, rate=fit$rate, period = 100)
calc_rl(threshold = 253, scale=scale_est - z*scale_se, shape=shape_est - z*shape_se, rate=fit$rate, period = 100)
calc_rl(threshold = 253, scale=scale_est + z*scale_se, shape=shape_est + z*shape_se, rate=fit$rate, period = 100)


### HOW TO: Find standard error (SE) of parameter estimates and return levels
window[[79]][[588]]
extRemes::parcov.fevd(window[[79]][[588]])
sqrt(extRemes::parcov.fevd(window[[79]][[588]])[1,1])
sqrt(extRemes::parcov.fevd(window[[79]][[588]])[2,2])
# Gives same output, but is messier output (still gives all of the summary)
# summary(window[[79]][[588]])$se.theta[1]
# summary(window[[79]][[588]])$se.theta[2]

### CI of parameter estimates
# extRemes::ci.fevd(window[[79]][[588]], type="parameter")
#Trying to figure out CI for return level

# from extRemes::ci.fevd.mle
cov.theta <- parcov.fevd(x)
var.theta <- diag(cov.theta)

if (is.element(x$type, c("PP", "GP", "Beta", "Pareto", 
                         "Exponential"))) 
  lam <- mean(c(datagrabber(x)[, 1]) > x$threshold)
# else lam <- 1
loc <- 0
theta.hat <- x$results$par
scale <- theta.hat["scale"]
shape <- theta.hat["shape"]
# if (x$type == "PP") 
#   mod <- "GEV"
else mod <- x$type
p <- rlevd(period = return.period, loc = loc, scale = scale, 
           shape = shape, threshold = x$threshold, type = mod, 
           npy = x$npy, rate = lam)
# p gives the basic return level!

else if (type == "return.level") {
  grads <- rlgrad.fevd(x, period = return.period)
  grads <- t(grads)
  if (is.element(x$type, c("GP", "Beta", "Pareto", 
                           "Exponential"))) {
    # if (x$type == "Exponential") 
    #   cov.theta <- diag(c(lam * (1 - lam)/x$n, var.theta))
    else cov.theta <- rbind(c(lam * (1 - lam)/x$n, 
                              0, 0), cbind(0, cov.theta))
  }
  # else lam <- 1
  var.theta <- t(grads) %*% cov.theta %*% grads
}
else stop("ci: invalid type argument.  Must be return.level or parameter.")
if (length(p) > 1) {
  if (type == "return.level") 
    se.theta <- sqrt(diag(var.theta))
  out <- cbind(p - z.alpha * se.theta, p, p + z.alpha * 
                 se.theta)
  rownames(out) <- par.name
  conf.level <- paste(round((1 - alpha) * 100, digits = 2), 
                      "%", sep = "")
  colnames(out) <- c(paste(conf.level, " lower CI", 
                           sep = ""), "Estimate", paste(conf.level, " upper CI", 
                                                        sep = ""))
  attr(out, "data.name") <- x$call
  attr(out, "method") <- method.name
  attr(out, "conf.level") <- (1 - alpha) * 100
  class(out) <- "ci"
  return(out)
}
else out <- c(p - z.alpha * sqrt(var.theta[which.par]), 
              p, p + z.alpha * sqrt(var.theta[which.par]))



## Aside - rlgrad.fevd()
p <- x$results$par
if (is.element("log.scale", names(p))) {
  id <- names(p) == "log.scale"
  p[id] <- exp(p[id])
  names(p)[id] <- "scale"
}
# if (p["shape"] == 0)
#   if (is.element(type, c("gp", "exponential"))) 
#     type <- "exponential"

else if (is.element(type, c("gp", "beta", "pareto"))) {
  lam <- mean(c(datagrabber(x)[, 1]) > x$threshold)
  m <- period * x$npy
  mlam <- m * lam
  res <- cbind(p["scale"] * m^(-p["shape"]) * lam^(-p["shape"] - 1), 
               (p["shape"])^(-1) * ((mlam)^(p["shape"]) - 1), 
               -p["scale"] * (p["shape"])^(-2) * ((mlam)^(p["shape"]) - 1) + (p["scale"]/p["shape"]) * (mlam)^(p["shape"]) * log(mlam))
}
# else if (type == "exponential") {
#   lam <- mean(c(datagrabber(x)[, 1]) > x$threshold)
#   m <- period * x$npy
#   mlam <- m * lam
#   res <- cbind(p["scale"]/mlam, log(mlam))
# }

rlgrad.fevd(fit, period=100)
fit$n
fit$rate



## Find SE of Return Level by using Conf Int
z = qnorm(0.975) # for 95% CI (default)
## Same return level CI Output from these 2 functions:
# extRemes::ci.fevd(window[[79]][[588]], alpha = 0.05, type = "return.level", return.period = 100)
# extRemes::return.level(window[[79]][[588]], return.period=100, do.ci=T)
test_rl <- extRemes::return.level(window[[79]][[588]], return.period=100)/10
((extRemes::return.level(window[[79]][[588]], return.period=100, do.ci=T)/10)[3] - test_rl)/z
## Examples from earlier code:
# ((extRemes::return.level(window[[79]][[588]], return.period=100, do.ci=T)/10)[3] - trend_RL[[588]][79])/z
## new_se_100[i] <- ((extRemes::return.level(window[[79]][[i]], return.period=100, do.ci=T)/10)[3] - new_rl_100[i])/z


# Extra Kriging/Variogram Code --------------------------------------------------------------------
#--------
# To real data:
g <- gstat(id = "ln.scale", formula = log(scale)~1, data = stations_sub)
g <- gstat(g, id = "shape", formula = shape~1, data = stations_sub)
# examine variograms and cross variogram:
vg <- variogram(g)
plot(vg)

# Trying variograms separately
g1 <- gstat(id = "ln.scale", formula = log(scale)~1, data = stations_sub)
g2 <- gstat(id = "shape", formula = shape~1, data = stations_sub)
g3 <- gstat(id = "rate", formula = rate~1, data = stations_sub)
vg1 <- variogram(g1)
vg2 <- variogram(g2)
vg3 <- variogram(g3)
var.fit1 <- fit.variogram(vg1, vgm("Mat")) # no convergence
var.fit2 <- fit.variogram(vg2, vgm("Mat")) # no convergence
# var.fit3 <- fit.variogram(vg3, vgm("Mat")) # no convergence, singular model  # range is at 44941.16
# var.fit3 <- fit.lmc(vg3, g3, vgm("Mat")) # won't plot  # using fit.lmc... seems to default to range of 44941.16...

sp::zerodist(stations_sub)

var.fit1
var.fit2
var.fit3

var.fit1 <- fit.variogram(vg1, vgm("Sph")) # singular model # Sph model gives larger range than Mat for these
var.fit2 <- fit.variogram(vg2, vgm("Sph")) # singular model
var.fit3 <- fit.variogram(vg3, vgm("Sph"))

plot(vg1, model = var.fit1, main = "Log(scale)")
plot(vg2, model = var.fit2, main = "shape")
plot(vg3, model = var.fit3, main = "rate")

abline(v=8955.58, col="red")
abline(v=50000, col="red")
lines(x=rep(8955.58, 2), y = c(0, 1), col="red")
lines(x=rep(50000, 2), y = c(0, 1), col="red")


## "Sph", "Gau", "Exp", "Mat"
cv.fit1 = fit.lmc(vg, g, vgm("Sph"))
# cv.fit1 = fit.lmc(vg, g, vgm("Sph"), correct.diagonal = 1.01)
cv.fit = fit.lmc(vg, g, vgm("Mat"))
plot(vg, model = cv.fit)
cv.fit
plot(vg, model = cv.fit1)
cv.fit1

cv.fitr = fit.lmc(vg, g, vgm("Mat"), fit.ranges = TRUE) # no convergence, fitting different ranges
plot(vg, model = cv.fitr)
cv.fitrs = fit.lmc(vg, g, vgm("Sph"), fit.ranges = TRUE) # singular model, fitting different ranges
plot(vg, model = cv.fitrs)

g.test <- gstat(g, model = vgm("Sph", range=37500), fill.all = TRUE)
test.fit1 <- fit.lmc(vg, g.test) # Does not actually set it to my range... but is a constant range of 44941.16
test.fit1 # same as cv.fit1 above
test.fit2 <- fit.lmc(vg, g.test, fit.ranges = TRUE) # singular model, fitting different ranges
test.fit2

test.fit = fit.lmc(vg, g, vgm("Mat"), fit.ranges = TRUE) # no convergence
test.fit = fit.lmc(vg, g, vgm("Sph"), fit.ranges = TRUE) # singular model - "a possible solution MIGHT be to scale semivariances and/or distances"
test.fit = fit.lmc(vg, g, vgm("Sph"), fit.ranges = TRUE, correct.diagonal = 1.01) # still singular, no difference

cv.fit$set=list(nocheck=1)
cv.fit1$set=list(nocheck=1)
pred <- predict(cv.fit1, newdata = stations_sub, nsim = 1)
pred <- predict(cv.fit, newdata = stations_sub, nsim = 1)
pred_orig <- pred
plot(ws_regs)
plot(pred, pch = "•", add = T)
spplot(pred) 

###

# Now try kriging these values
sg <- gstat(id = "ln.scale", formula = log(scale)~1, data = sim_fit)
sg <- gstat(sg, id = "shape", formula = shape~1, data = sim_fit)
# examine variograms and cross variogram:
svg <- variogram(sg)
plot(svg)
## "Sph", "Gau", "Exp", "Mat"
s.cv.fit1 = fit.lmc(svg, sg, vgm("Sph")) #, correct.diagonal = 1.01)
s.cv.fit = fit.lmc(svg, sg, vgm("Mat"))
plot(svg, model = s.cv.fit)
s.cv.fit
plot(svg, model = s.cv.fit1)
s.cv.fit1

s.cv.fit$set=list(nocheck=1)
# cv.fit1$set=list(nocheck=1)
skrigws <- predict(s.cv.fit1, newdata = ws_regs)
skrigws <- predict(s.cv.fit, newdata = ws_regs)

skrig <- predict(s.cv.fit1, newdata = stations_sub)
skrig <- predict(s.cv.fit, newdata = stations_sub)

plot(ws_regs)
plot(skrig, pch = "•", add = T)
spplot(skrig) 
spplot(skrigws)
spplot(skrigws[3]) # ln.scale.pred
spplot(skrigws[5]) # shape.pred

# Need to add in rate parameter somehow

### Decision: shape and scale parameters are negatively correlated, 
#   but the rate parameter can be considered independent 
#   (makes sense because rate is determined by the amount of data above the threshold, which differs for each station)


### Trying to plot cross-variogram
# coordinates(stations_sub)=~long+lat
# let's do some manual fitting of two direct variograms and a cross variogram
g <- gstat(id = "ln.scale", formula = log(scale)~1, data = stations_sub)
g <- gstat(g, id = "shape", formula = shape~1, data = stations_sub)
# examine variograms and cross variogram:
plot(variogram(g))
vg <- variogram(g)
## "Sph", "Gau", "Exp", "Mat"
cv.fit = fit.lmc(vg, g, vgm("Sph"))
cv.fit1 = fit.lmc(vg, g, vgm("Sph"), correct.diagonal = 1.01)
# cv.fit = fit.lmc(vg, g, vgm("Gau"))
# cv.fit = fit.lmc(vg, g, vgm("Exp"))
# cv.fit = fit.lmc(vg, g, vgm("Mat"))
plot(vg, model = cv.fit)
cv.fit

try.cok <- predict(cv.fit, newdata = ws_regs) # cokriging?
plot(try.cok)

plot(vg, model = cv.fit1)
cv.fit1
try.cok <- predict(cv.fit1, newdata = ws_regs) # cokriging?
plot(try.cok)

### Eduardo added the following line to make it work
cv.fit$set=list(nocheck=1)
ws_reg_grid <- predict(cv.fit, newdata = ws_regs, nsim = 1)

### create a grid (of points)
# library(raster)
grid <- sp::makegrid(ws_regs, cellsize = 2000) # cellsize in map units
grid <- SpatialPoints(grid, proj4string = CRS(proj4string(ws_regs)))
grid <- grid[ws_regs, ] # grid points only within the polygon
plot(ws_regs)
plot(grid, pch = ".", add = T)
try.cok <- predict(cv.fit1, newdata = grid) # cokriging?
plot(try.cok)
gridded(grid) = TRUE # turn into spatial pixels data frame
