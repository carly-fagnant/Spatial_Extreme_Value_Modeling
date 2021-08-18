#######################
# Simulation of new data
#
# by Carly Fagnant
#######################

# Using the data we have, capture the mean and covariance structure.
# The mean structure will come in at the region level - we will assume the same mean value for all points that fall within a region
# The covariance structure comes in at the point level... at least that is all we can get from our raw data.

# maybe gstat::vgmArea could provide covariance at the area level?


# Updated Steps 1 & 2: regress our data parameter values against constant means for the region, and save the residuals
# Then do variogram modeling and spatial simulations on the residuals. Add back in the mean afterwards

# Step 1: use multivariate normal (separately for each region) to simulate (shape and scale) parameter values for stations within each region
# But what covariance structure do we use here? Same for all regions? (yes, I believe so)

# Are the parameter values different for each station? Yes. But only "simulate" them once.

# Step 2: fit a cross-variogram in gstat (to this simulated data?) and use the predict fn with nsim = 1 to simulate parameter values at all of the stations at once
# This brings in more of the covariance structure?

# Note: maybe if I can access covariance matrix from gstat cross-variogram, I can combine Steps 1 and 2 into just a multivariate normal. 
#   However, I run into the issue of that I am doing it separately for each region

# Step ?: use the extRemes::revd fn to simulate daily rainfall data using the parameters given (be sure to transorm scale back using exp())
# NOTE: No input for rate here. Will outputs be all above the threshold? Yes.
#       Will need to correct for this, by randomly removing observations at the rate estimated



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


# Testing Functions -------------------------------------------------------

### Discoveries:
# - revd function to generate GPD realizations does NOT create any that are less than the threshold. (Makes sense, since we only fit to those above the threshold.)
#   This means that if we want our simulated data to look more like real rainfall data, we will need to do one of the following:
#     - randomly select a proportion of data to remove (or set to zero) or replace with uniformly distributed values < threshold
#     - subtract a set amount (quantile value?) from each observation --> NO, this would change the parameter values and return levels
#  Recall that we assumed the observations were independent (after declustering) before fitting, so do not let multiple days in a row have values above the threshold?
#  We will decluster our simulated data, but want to make sure we have the proper proportion of data left to work with

# NOTE: the avg rate of exceedance of the threshold is only 0.0543096, so about 95% of the data simulated should be under the threshold...

avg_sc <- mean(ws_reg_avg[1,])
avg_sh <- mean(ws_reg_avg[2,])
avg_rt <- mean(ws_reg_avg[3,])

simtry <- extRemes::revd(n = 10000, type = "GP", scale = avg_sc, shape = avg_sh, threshold = 253) # use avg param values
any(simtry < 253)



# Step 1 ------------------------------------------------------------------
# Step 1: use multivariate normal (separately for each region) to simulate (shape and scale) parameter values for stations within each region
# But what covariance structure do we use here? Same for all regions? (yes, I believe so)
# Note this covariance structure is strictly on the 2 parameter values, and has no spatial component

Varlogsc <- var(log(stations_sub_df$scale))
Varsh <- var(stations_sub_df$shape)
Covar <- cov(log(stations_sub_df$scale), stations_sub_df$shape)
cov_mat <- matrix(c(Varlogsc, Covar , Covar, Varsh), nrow = 2, ncol = 2, byrow = TRUE,
                                   dimnames = list(c("ln.scale", "shape"),
                                                   c("ln.scale", "shape")))

reg1mvn <- MASS::mvrnorm(n=3, mu = c(log(ws_reg_avg[1,1]), ws_reg_avg[2,1]), Sigma = cov_mat)

# Number of stations per region, from summary(region_intersect$REGION)
stats_per_reg <- c(50, 52, 64)

sim_from_mvn <- function(region){
  i <- region
  if(i %in% c(1,2,3)){
    # simulates the 2 params for the number of stations in that region, using that region's avg values
    sim_params <- MASS::mvrnorm(n = stats_per_reg[i], mu = c(log(ws_reg_avg[1,i]), ws_reg_avg[2,i]), Sigma = cov_mat)
    return(sim_params)
  }else{
    print("Invalid entry of region, please input 1, 2, or 3")
  }
}

reg1mvn <- sim_from_mvn(region = 1)
reg2mvn <- sim_from_mvn(region = 2)
reg3mvn <- sim_from_mvn(region = 3)

# Now assign these values to station locations?

# load data
station_info <- read.csv("~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/station_info.csv")

# add new location (long/lat) columns that we will transform to coordinates
station_info <- station_info %>%
  dplyr::mutate(long = LONGITUDE, lat = LATITUDE)

# subset to those 166 stations that fall within the 3 regions
# ... but also remove stations with NAs if any
# ends up being 149 stations
stations_sub_sim <- station_info %>%
  dplyr::filter(STAT_NO %in% inds) 
# %>%
#   dplyr::filter(!is.na(scale))
stations_sub_sim <- stations %>%
  dplyr::filter(STAT_NO %in% inds) 


### Function to project data given a data frame with data and coordinates
project_data <- function(df){
  coordinates(df)=~long+lat
  proj4string(df) <- "+proj=longlat +datum=WGS84"
  df <- spTransform(df, CRS("+init=epsg:2278")) # projecting data
  is.projected(df)
  return(df)
}

station_info <- project_data(station_info)

# set regions and points to have the same projection
proj4string(ws_regs) <- proj4string(station_info)





### Finding which stations fall within the 3 regions 
region_intersect <- sp::over(station_info, ws_regs)
# region_intersect lists the 601 stations in order and the integer of what region each station falls within (NA if not in any region)

region_intersect$REGION <- as.factor(region_intersect$REGION) # need these to be factors in order for the model matrix and summary to work
summary(region_intersect$REGION) # number of stations in each region (NAs are those not in any)

h <- model.matrix( ~ REGION - 1, data=region_intersect)

# this gives the numbers of the stations that fall within the 3 regions
inds <- as.integer(rownames(h))

# stations by region
reg1stat_index <- which(h[,1]==1) # index of station in stations_sub_df
reg1stat_no  <- as.integer(names(reg1stat_index)) # actual station number
reg2stat_index <- which(h[,2]==1)
reg2stat_no  <- as.integer(names(reg2stat_index))
reg3stat_index <- which(h[,3]==1)
reg3stat_no  <- as.integer(names(reg3stat_index))

stations_sub_mvn <- stations_sub_sim %>% mutate(mvln.scale = NA, mvshape = NA)

stations_sub_mvn$mvln.scale[reg1stat_index] <- reg1mvn[,1]
stations_sub_mvn$mvshape[reg1stat_index] <- reg1mvn[,2]

stations_sub_mvn$mvln.scale[reg2stat_index] <- reg2mvn[,1]
stations_sub_mvn$mvshape[reg2stat_index] <- reg2mvn[,2]

stations_sub_mvn$mvln.scale[reg3stat_index] <- reg3mvn[,1]
stations_sub_mvn$mvshape[reg3stat_index] <- reg3mvn[,2]

# getting rid of those stations which were NA in actual data. Use same stations as in stations_sub
stations_sub_mvn <- stations_sub_mvn %>% filter(X %in% stations_sub_df$X)

# Project data but also save data frame version in case it's needed
stations_sub_mvn_df <- stations_sub_mvn
stations_sub_mvn <- project_data(stations_sub_mvn)

# Have now created data set using mvrnorm data --> next run variograms and simulation


# Step 2 ------------------------------------------------------------------
# Step 2: fit a cross-variogram in gstat (to this simulated data?) and use the predict fn with nsim = 1 to simulate parameter values at all of the stations at once
# This brings in more of the covariance structure? - Yes, the spatial side of it

# Note: maybe if I can access covariance matrix from gstat cross-variogram, I can combine Steps 1 and 2 into just a multivariate normal. 
#   However, I run into the issue of that I am doing it separately for each region


# To real data:
g <- gstat(id = "ln.scale", formula = log(scale)~1, data = stations_sub)
g <- gstat(g, id = "shape", formula = shape~1, data = stations_sub)
# examine variograms and cross variogram:
vg <- variogram(g)
plot(vg)
## "Sph", "Gau", "Exp", "Mat"
cv.fit1 = fit.lmc(vg, g, vgm("Sph"), correct.diagonal = 1.01)
cv.fit = fit.lmc(vg, g, vgm("Mat"))
plot(vg, model = cv.fit)
cv.fit
plot(vg, model = cv.fit1)
cv.fit1

g.test <- gstat(g, model = vgm("Sph", range=37500), fill.all = TRUE)
test.fit1 <- fit.lmc(vg, g.test) # Does not actually set it to my range... but is a constant range of 44941.16
test.fit1
test.fit2 <- fit.lmc(vg, g.test, fit.ranges = TRUE) # singular model, fitting different ranges
test.fit2

test.fit = fit.lmc(vg, g, vgm("Mat"), fit.ranges = TRUE) # no convergence
test.fit = fit.lmc(vg, g, vgm("Sph"), fit.ranges = TRUE) # singular model - "a possible solution MIGHT be to scale semivariances and/or distances"

cv.fit$set=list(nocheck=1)
cv.fit1$set=list(nocheck=1)
pred <- predict(cv.fit1, newdata = stations_sub, nsim = 1)
pred <- predict(cv.fit, newdata = stations_sub, nsim = 1)
pred_orig <- pred
plot(ws_regs)
plot(pred, pch = "•", add = T)
spplot(pred) 


# To parameters from simulated mvnorm:
g1 <- gstat(id = "ln.scale", formula = mvln.scale~1, data = stations_sub_mvn)
g1 <- gstat(g1, id = "shape", formula = mvshape~1, data = stations_sub_mvn)
# examine variograms and cross variogram:
vg1 <- variogram(g1)
## "Sph", "Gau", "Exp", "Mat"
cv.fit1s = fit.lmc(vg1, g1, vgm("Sph"), correct.diagonal = 1.01)
cv.fit1m = fit.lmc(vg1, g1, vgm("Mat"))
plot(vg1, model = cv.fit1s)
cv.fit1s
plot(vg1, model = cv.fit1m)
cv.fit1m

cv.fit1s$set=list(nocheck=1)
cv.fit1m$set=list(nocheck=1)
pred1s <- predict(cv.fit1s, newdata = stations_sub, nsim = 1)
pred1m <- predict(cv.fit1m, newdata = stations_sub, nsim = 1)
plot(ws_regs)
plot(pred, pch = "•", add = T)
spplot(pred1s)
spplot(pred1m)

### Reasons for some errors- 
# Covariance matrix singular at location... from having multiple observations at the same location (distance=0 between pairs)
sp::zerodist(stations_sub)
sp::zerodist(stations_sub_mvn)
# (444, 491)  (471, 492)  (481, 493, 502)  (489, 495)  (508, 496)
# The stations with the same locations (were removed by removing NA values for regular stations_sub)
# the first one listed in each pair has only ~1 yr of data so can be thrown out or absorbed into the other station

# No Intrinsic Correlation or Linear Model of Coregionalization found
# Reason: ranges differ
# The strange thing is our ranges do not differ in the fit object... so add `set = list(nocheck = 1)' to ignore
# See https://stat.ethz.ch/pipermail/r-sig-geo/2008-November/004562.html



# Updated Steps 1 & 2 -----------------------------------------------------
# Updated Steps 1 & 2: regress our data parameter values against constant means for the region, and save the residuals
# Then do variogram modeling and spatial simulations on the residuals. Add back in the means afterwards

# Finding which region the subset stations fall within
region_int <- sp::over(stations_sub, ws_regs)
stations_regs <- stations_sub_df %>% mutate(REGION = region_int$REGION) # adding region information to the spatially referenced data

region_int$REGION <- as.factor(region_int$REGION)
summary(region_int$REGION)
h0 <- model.matrix( ~ REGION - 1, data=region_int)
stations_regs_cat <- stations_sub_df %>% mutate(REGION = region_int$REGION) # adding region information to the spatially referenced data
class(stations_regs_cat$REGION)

# shape
model <- lm(shape ~ REGION, data = stations_regs_cat)
summary(model)
summary(model)$coef[1]
summary(model)$coef[1] + summary(model)$coef[2]
summary(model)$coef[1] + summary(model)$coef[3]
shape.means <- c(summary(model)$coef[1],
                 summary(model)$coef[1] + summary(model)$coef[2],
                 summary(model)$coef[1] + summary(model)$coef[3])

plot(model$residuals)
hist(model$residuals)
hist(log(model$residuals))
qqnorm(model$residuals)
qqnorm(log(model$residuals))

# log(scale)
model.sc <- lm(log(scale) ~ REGION, data = stations_regs_cat)
summary(model.sc)
summary(model.sc)$coef[1]
summary(model.sc)$coef[1] + summary(model.sc)$coef[2]
summary(model.sc)$coef[1] + summary(model.sc)$coef[3]
ln.scale.means <- c(summary(model.sc)$coef[1],
                    summary(model.sc)$coef[1] + summary(model.sc)$coef[2],
                    summary(model.sc)$coef[1] + summary(model.sc)$coef[3])
model.sc$residuals
hist(model.sc$residuals)
# hist(log(model.sc$residuals))
qqnorm(model.sc$residuals)
# qqnorm(log(model.sc$residuals))

# stations_resids <- stations_regs_cat %>% mutate(ln.scale.resid = log(model.sc$residuals), shape.resid = model$residuals) # taking log after regressing for scale created NaNs
stations_resids <- stations_regs_cat %>% mutate(ln.scale.resid = model.sc$residuals, shape.resid = model$residuals)
stations_resids <- project_data(stations_resids)


# To residuals:
rg <- gstat(id = "ln.scale.resid", formula = ln.scale.resid~1, data = stations_resids)
rg <- gstat(rg, id = "shape.resid", formula = shape.resid~1, data = stations_resids)
# examine variograms and cross variogram:
rvg <- variogram(rg)
plot(rvg)
## "Sph", "Gau", "Exp", "Mat"
r.cv.fit1 = fit.lmc(rvg, rg, vgm("Sph")) #, correct.diagonal = 1.01)
r.cv.fit = fit.lmc(rvg, rg, vgm("Mat"))
plot(rvg, model = r.cv.fit)
r.cv.fit
plot(rvg, model = r.cv.fit1)
r.cv.fit1

r.cv.fit$set=list(nocheck=1)
r.cv.fit1$set=list(nocheck=1)
pred1 <- predict(r.cv.fit1, newdata = stations_sub, nsim = 1)
pred <- predict(r.cv.fit, newdata = stations_sub, nsim = 1)


station_locations <- stations_sub_df %>% dplyr::select(STAT_NO, long, lat)
station_locations <- project_data(station_locations)

pred_test <- predict(r.cv.fit, newdata = station_locations, nsim = 1)
pred@data$ln.scale.resid.sim1 == pred_test@data$ln.scale.resid.sim1
pred@data$ln.scale.resid.sim1 - pred_test@data$ln.scale.resid.sim1

plot(ws_regs)
plot(pred, pch = "•", add = T)
spplot(pred) 
spplot(pred[1]) 
spplot(pred[2]) 

spplot(pred_test)

# Make sure to add back in the means!!

ln.scale.resid <- pred@data$ln.scale.resid.sim1
shape.resid <- pred@data$shape.resid.sim1
# scale <- exp(pred@data$ln.scale.resid.sim1) # transforming it back from the log scale

scale <- shape <- NULL
for (i in 1:length(shape.resid)){
  # value = residual + mean of region of station (index i)
  scale[i] <- exp( ln.scale.resid[i] + ln.scale.means[stations_regs$REGION[i]] )  # transforming it back from the log scale
  shape[i] <- shape.resid[i] + shape.means[stations_regs$REGION[i]]
}

# Strange... the shape and scale parameters are the same as in stations_sub...
# Maybe in predict function we cannot let newdata = stations_sub? No, doesn't matter if we remove
# They are just going to be REALLY close (almost zero difference)
# I think this is because of the kriging nature and that we are going to the points that we worked off of.
# They are going to be very precise unless we move to new points/grid/area

shape - shape_orig
scale - scale_orig
shape - stations_sub@data$shape
scale - stations_sub@data$scale

# Step 3 ------------------------------------------------------------------
# Step 3: use the extRemes::revd fn to simulate daily rainfall data using the parameters given (be sure to transorm scale back using exp())
# NOTE: No input for rate here. Will outputs be all above the threshold? Yes.
#       Will need to correct for this, by randomly removing observations at the rate estimated

# pred is a SpatialPointsDataFrame holding the parameter estimates for each station (ln.scale and shape)
# be sure to transorm scale back using exp() fn

# Save data as a matrix? with station numbers as the column titles and the data listed under each for the same number of days
# 14610 is the number of days in 40 years


# # # Can ignore if did residuals version above
# scale_orig <- exp(pred_orig@data$ln.scale.sim1) # transforming it back from the log scale
# shape_orig <- pred_orig@data$shape.sim1

# Haven't simulated the rate parameters yet, so let's just use the rates from stations_sub
rate <- stations_sub@data$rate

# idea - use the rates to randomly select the data to remain above the threshold
# Helpful note: it doesn't matter the order/placement of the data in the days, other than assuring they are independent/declustered

sim_data <- matrix(nrow = 14610, ncol = length(scale))
for(i in 1:length(scale)){
  # stat_no <- stations_sub@data$STAT_NO[i]
  dat <- extRemes::revd(n = 14610, type = "GP", scale = scale[i], shape = shape[i], threshold = 253)
  n = round(rate[i]*length(dat)) # number of values to keep above threshold
  sub_dat <- sample(dat, size = n) # taking a random sample of the generated values to keep according to rate
  new_dat <- c(sub_dat, rep(0, 14610-n)) # appending on zeros for the rest of the data
  sim_data[, i] <- new_dat
}
sim_data <- as.data.frame(sim_data)
colnames(sim_data) <- stations_sub@data$STAT_NO

## saving
# setwd("~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data")
# saveRDS(sim_data, file = "ex_sim_data.rds") # first example
# saveRDS(sim_data, file = "reg_sim_data.rds") # example after simulating using the regression & residuals idea


# Model 2 - Kriging and Aggregation ---------------------------------------

# Note that one may take a shortcut in coding by setting the newdata argument in gstat::predict to be the regions of interest
# this uses spsample to draw unifrom random locations across the polygon and average them
# Not quite our gridded kriging, but very close

sim_data <- readRDS("~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/ex_sim_data.rds")

thresh <- 253

# Testing out model 2 on example simulated data
fits <- NULL
sim_fit_df <- matrix(nrow = length(scale), ncol = 6)

for (i in 1:length(scale)){
  fit <- extRemes::fevd(sim_data[, i], threshold = thresh, type = "GP", method="MLE")
  fits[[i]] <- fit # saving just in case
  # Create new data frame with station and fit information
  sim_fit_df[i, 2] <- fit$results$par[1] # scale
  sim_fit_df[i, 3] <- fit$results$par[2] # shape
  sim_fit_df[i, 4] <- fit$rate # rate
}

sim_fit_df[, 1] <- stations_sub@data$STAT_NO # 1st column is station numbers
sim_fit_df[, 5] <- stations_sub_df$long # long lat
sim_fit_df[, 6] <- stations_sub_df$lat # long lat
sim_fit <- as_data_frame(sim_fit_df)
colnames(sim_fit) <- c("STAT_NO", "scale", "shape", "rate", "long", "lat")

sim_fit <- project_data(sim_fit)

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
spplot(skrigws[3]) 

# Need to add in rate parameter somehow