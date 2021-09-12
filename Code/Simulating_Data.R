#######################
# Simulating Data
#
# by Carly Fagnant
#######################

library(rgdal)
library(sp)
library(gstat)
library(dplyr)
library(spdep)
library(spatialreg)
library(foreach)

# Using the data we have, capture the mean and covariance structure.
# The mean structure will come in at the region level - we will assume the same mean value for all points that fall within a region.
# The covariance structure comes in at the point level - we can retrieve this structure from our raw data.



# Step 1: Regress our data parameter values against constant means for the region, and save the residuals

# Step 2: Perform variogram modeling and spatial simulations on the residuals. Add back in the mean afterwards

# Are the parameter values different for each station? Yes. But only "simulate" them once.

# Step 3: use the extRemes::revd fn to simulate daily rainfall data using the parameters given (be sure to transorm scale back using exp())
# NOTE: Holding the rate constant, we will only generate as many observations as would be above the threshold. (14610 * rate for 40-year)



# Loading our data --------------------------------------------------------

setwd("~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Test")
wd <- getwd()

ws_regs <- rgdal::readOGR(dsn = paste0(wd, "/watershed_region_updated"), layer = "watershed_region_updated")
plot(ws_regs)

# load pre-cleaned spatial points data frame
stations_sub <- readRDS(file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/stations_sub.rds")
stations_sub_df <- readRDS(file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/stations_sub_df.rds")

# put in same projection!
proj4string(ws_regs) <- proj4string(stations_sub)

# avg parameter values by region
ws_reg_avg <- readRDS("~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/ws_reg_avg.rds")




# Step 1 -----------------------------------------------------
# Regress our data parameter values against constant means for the region, and save the residuals

# Finding which region the subset stations fall within
region_int <- sp::over(stations_sub, ws_regs)
stations_regs <- stations_sub_df %>% mutate(REGION = region_int$REGION) # adding region information to the spatially referenced data

# Same, but make REGION a factor for easy regression with categorical variables
region_int$REGION <- as.factor(region_int$REGION)
summary(region_int$REGION)
h0 <- model.matrix( ~ REGION - 1, data=region_int)
stations_regs_cat <- stations_sub_df %>% mutate(REGION = region_int$REGION) # adding region information to the spatially referenced data, as a factor

# shape
model <- lm(shape ~ REGION - 1, data = stations_regs_cat)
summary(model)
shape.means <- c(summary(model)$coef[1],
                 summary(model)$coef[2],
                 summary(model)$coef[3])
plot(model$residuals)
hist(model$residuals)
qqnorm(model$residuals)

# log(scale)
model.sc <- lm(log(scale) ~ REGION - 1, data = stations_regs_cat)
summary(model.sc)
ln.scale.means <- c(summary(model.sc)$coef[1],
                    summary(model.sc)$coef[2],
                    summary(model.sc)$coef[3])
plot(model.sc$residuals)
hist(model.sc$residuals)
qqnorm(model.sc$residuals)

# # scale - Issue, can't take log of negative numbers, so taking log of residuals won't work. Use log of scale parameters, then find residuals
model.scale <- lm(scale ~ REGION - 1, data = stations_regs_cat)
# summary(model.sc)
scale.means <- c(summary(model.scale)$coef[1],
                 summary(model.scale)$coef[2],
                 summary(model.scale)$coef[3])
# plot(model.sc$residuals)
# hist(model.sc$residuals)
# hist(log(model.sc$residuals))
# qqnorm(model.sc$residuals)
# qqnorm(log(model.sc$residuals))
# stations_resids <- stations_regs_cat %>% mutate(ln.scale.resid = log(model.sc$residuals), shape.resid = model$residuals) # taking log after regressing for scale created NaNs

stations_resids_df <- stations_regs_cat %>% mutate(ln.scale.resid = model.sc$residuals, shape.resid = model$residuals)
stations_resids <- project_data(stations_resids_df)


# Testing the means
region1_test <- stations_regs %>% filter(REGION == 1)
region2_test <- stations_regs %>% filter(REGION == 2)
region3_test <- stations_regs %>% filter(REGION == 3)

mean(region1_test$scale)

scale.means
exp(ln.scale.means)
exp(mean(log(region1_test$scale))) * (1 + var(log(region1_test$scale))/2)
exp(mean(log(region2_test$scale))) * (1 + var(log(region2_test$scale))/2)
exp(mean(log(region3_test$scale))) * (1 + var(log(region3_test$scale))/2)

# same as exp(ln.scale.means)
exp(mean(log(region1_test$scale))); exp(mean(log(region2_test$scale))); exp(mean(log(region3_test$scale)))

# Step 2 ------------------------------------------------------------------
# Perform variogram modeling and spatial simulations on the residuals. Add back in the means afterwards

# No issue of repeated locations
sp::zerodist(stations_resids)

# To residuals:
rg <- gstat(id = "ln.scale.resid", formula = ln.scale.resid~1, data = stations_resids)
rg <- gstat(rg, id = "shape.resid", formula = shape.resid~1, data = stations_resids)
# examine variograms and cross variogram:
rvg <- variogram(rg)
plot(rvg)
## "Sph", "Gau", "Exp", "Mat"
r.cv.fit1 = fit.lmc(rvg, rg, vgm("Sph"))
r.cv.fit = fit.lmc(rvg, rg, vgm("Mat"))
plot(rvg, model = r.cv.fit)
r.cv.fit
plot(rvg, model = r.cv.fit1)
r.cv.fit1

r.cv.fit$set=list(nocheck=1)
r.cv.fit1$set=list(nocheck=1)
pred1 <- predict(r.cv.fit1, newdata = stations_sub, nsim = 1)
pred_at_stat <- predict(r.cv.fit, newdata = stations_sub, nsim = 1)
# Because of kriging structure, we are going to get very precise (the same values we already have) unless we move to new points/grid/area

### create a grid (of points)
grid <- sp::makegrid(ws_regs, cellsize = 10000) # cellsize in map units
grid <- SpatialPoints(grid, proj4string = CRS(proj4string(ws_regs)))
grid <- grid[ws_regs, ] # grid points only within the polygon
plot(ws_regs)
plot(grid, pch = "+", add = T, cex = 0.7)
plot(grid, pch = ".", add = T)
plot(pred_at_stat, pch = "•", add = T)
# gridded(grid) = TRUE # turn into spatial pixels data frame


grid3m <- sp::makegrid(ws_regs, cellsize = 3 * 5280) # cellsize in map units
grid3m <- SpatialPoints(grid3m, proj4string = CRS(proj4string(ws_regs)))
grid3m <- grid3m[ws_regs, ] # grid points only within the polygon
plot(ws_regs, main = "3-mile grid")
plot(grid3m, pch = ".", add = T)
plot(pred_at_stat, pch = "•", add = T)

# saveRDS(grid3m, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/grid3m.rds")
plot(ws_regs, main = "3-Mile Grid for Simulations")
plot(grid3m, pch = "+", add = T, cex = 0.7)
plot(pred_at_stat, pch = "•", col = "blue", add = T)

pred <- predict(r.cv.fit, newdata = grid3m, nsim = 1)

spplot(pred) 
spplot(pred[1]) 
spplot(pred[2]) 

spplot(pred_at_stat)
spplot(pred_at_stat[1]) 
spplot(pred_at_stat[2])

range(pred@data$ln.scale.resid.sim1)
range_ln.scale <- range(pred_at_stat@data$ln.scale.resid.sim1)
brk <- seq(from=range_ln.scale[1], to=range_ln.scale[2], length.out = 7)

spplot(pred[1], at=brk) 



# Make sure to add back in the means!!

ln.scale.resid <- pred@data$ln.scale.resid.sim1
shape.resid <- pred@data$shape.resid.sim1

# Need to edit now that we have new data on a grid - find out which region each new point falls within
region_int_grid <- sp::over(pred, ws_regs)
grid_regs_df <- cbind(pred@data, pred@coords)
grid_regs_df <- grid_regs_df %>% mutate(REGION = region_int_grid$REGION) # adding region information to the spatially referenced data

scale <- shape <- NULL
for (i in 1:length(shape.resid)){
  # value = residual + mean of region of station (index i)
  scale[i] <- exp( ln.scale.resid[i] + ln.scale.means[grid_regs_df$REGION[i]] )  # transforming it back from the log scale
  shape[i] <- shape.resid[i] + shape.means[grid_regs_df$REGION[i]]
}

# range(scale)
# range(shape)
# range(stations_sub_df$scale)
# range(stations_sub_df$shape)

# saveRDS(grid_regs_df, file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/grid_regs_df.rds")


### Testing the means and variances
region1_df <- grid_regs_df %>% filter(REGION == 1)
region2_df <- grid_regs_df %>% filter(REGION == 2)
region3_df <- grid_regs_df %>% filter(REGION == 3)

mean(region1_test$scale)
summary(region1_df$ln.scale.resid.sim1); range(region1_df$ln.scale.resid.sim1)[2] - range(region1_df$ln.scale.resid.sim1)[1]
summary(region2_df$ln.scale.resid.sim1); range(region2_df$ln.scale.resid.sim1)[2] - range(region2_df$ln.scale.resid.sim1)[1]
summary(region3_df$ln.scale.resid.sim1); range(region3_df$ln.scale.resid.sim1)[2] - range(region3_df$ln.scale.resid.sim1)[1]

var(region1_df$ln.scale.resid.sim1)
var(region2_df$ln.scale.resid.sim1)
var(region3_df$ln.scale.resid.sim1)

scale.means
exp(ln.scale.means)
exp(mean(log(region1_test$scale))) * (1 + var(log(region1_test$scale))/2)
exp(mean(log(region2_test$scale))) * (1 + var(log(region2_test$scale))/2)
exp(mean(log(region3_test$scale))) * (1 + var(log(region3_test$scale))/2)

exp(mean(log(region1_test$scale)))

## Testing different variances for transformation
scale.var.real <- NULL
for (i in 1:length(ln.scale.resid)){
  # value = residual + mean of region of station (index i)
  reg <- grid_regs_df$REGION[i]
  if(reg==1){
    var.ln.sc = var(log(region1_test$scale))
  }else if(reg==2){
    var.ln.sc = var(log(region2_test$scale))
  }else if(reg==3){
    var.ln.sc = var(log(region3_test$scale))
  }
  scale.var.real[i] <- exp( ln.scale.resid[i] + ln.scale.means[grid_regs_df$REGION[i]] ) * (1 + var.ln.sc/2)
}

scale.var.resid <- NULL
for (i in 1:length(ln.scale.resid)){
  # value = residual + mean of region of station (index i)
  reg <- grid_regs_df$REGION[i]
  if(reg==1){
    var.ln.sc = var(region1_df$ln.scale.resid.sim1)
  }else if(reg==2){
    var.ln.sc = var(region2_df$ln.scale.resid.sim1)
  }else if(reg==3){
    var.ln.sc = var(region3_df$ln.scale.resid.sim1)
  }
  scale.var.resid[i] <- exp( ln.scale.resid[i] + ln.scale.means[grid_regs_df$REGION[i]] ) * (1 + var.ln.sc/2)
}

compare <- cbind(scale, scale.var.real, scale.var.resid)

ind_reg1 <- which(grid_regs_df$REGION==1)
ind_reg2 <- which(grid_regs_df$REGION==2)
ind_reg3 <- which(grid_regs_df$REGION==3)
apply(compare[ind_reg1, ], 2, mean)
apply(compare[ind_reg2, ], 2, mean)
apply(compare[ind_reg3, ], 2, mean)



### Holding the parameters constant within regions
constant_scale <- constant_shape <- NULL
round.scale.means <- round(scale.means, 2)
round.shape.means <- round(shape.means, 4)
for (i in 1:length(shape.resid)){
  # value = residual + mean of region of station (index i)
  constant_scale[i] <- round.scale.means[grid_regs_df$REGION[i]]
  constant_shape[i] <- round.shape.means[grid_regs_df$REGION[i]]
}

compare2 <- cbind(grid_regs_df, constant_scale, constant_shape)

# Step 3 ------------------------------------------------------------------
# Step 3: use the extRemes::revd fn to simulate daily rainfall data using the parameters given (be sure to transform scale back using exp())
# NOTE: No input for rate here. Will outputs be all above the threshold? Yes.
#       Will need to correct for this, by only generating as many observations as rate entails

# pred is a SpatialPointsDataFrame holding the parameter estimates for each station (ln.scale and shape)
# be sure to transorm scale back using exp() fn. This is done above.

# Save data as a matrix- with station numbers as the column titles and the data listed under each for the same number of days
# 14610 is the number of days in 40 years


# We decided from our other results that the rate parameter did not change by much, so to simplify modeling
# we hold the rate parameter constant at the mean, 0.0544
rate <- stations_sub@data$rate
rate <- round(mean(rate), 4)
# mean(c(.057, .0514, .0545, .0559, .0554, .0557, .0565, .0601, .0557))
# mean(ws_reg_avg[3,])

# Helpful note: it doesn't matter the order/placement of the data in the days, other than assuring they are independent/declustered


# only generate as many as I need instad of randomly subsetting!
# just use a constant rate,  0.0544
ngen <- round(14610 * rate)
sim_data <- matrix(nrow = 14610, ncol = length(scale))
for(i in 1:length(scale)){
  sub_dat <- extRemes::revd(n = ngen, type = "GP", scale = scale[i], shape = shape[i], threshold = 253) # only generate as much data as rate says
  new_dat <- c(sub_dat, rep(0, 14610-ngen)) # appending on zeros for the rest of the data
  sim_data[, i] <- new_dat
}
sim_data <- as.data.frame(sim_data)
# colnames(sim_data) <- stations_sub@data$STAT_NO


test_sim_data <- sim_data

## saving
# setwd("~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data")
# saveRDS(sim_data, file = "ex_sim_data.rds") # first example
# saveRDS(sim_data, file = "reg_sim_data.rds") # example after simulating using the regression & residuals idea
# saveRDS(test_sim_data, file = "test_sim_data.rds") # example using 3-mile grid and constant rate



## Updated to use constant shape and scale by region!
simulate_data <- function(nsim){  # , scale, shape, rate
  # only generate as many as I need instad of randomly subsetting!
  # just use a constant rate,  0.0544
  sim_dat_list <- list()
  ngen <- round(14610 * rate)
  for (j in 1:nsim){
    sim_data <- matrix(nrow = 14610, ncol = length(constant_scale))
    for(i in 1:length(constant_scale)){
      sub_dat <- extRemes::revd(n = ngen, type = "GP", scale = constant_scale[i], shape = constant_shape[i], threshold = 253) # only generate as much data as rate says
      new_dat <- c(sub_dat, rep(0, 14610-ngen)) # appending on zeros for the rest of the data
      sim_data[, i] <- new_dat
    }
    sim_data <- as.data.frame(sim_data)
    sim_dat_list[[j]] <- sim_data
  }
  
  if(nsim==1){
    return(sim_data)
  }else{
    return(sim_dat_list)
  }
}

ptm <- proc.time()
test_sims_constant <- simulate_data(nsim = 50)
proc.time() - ptm