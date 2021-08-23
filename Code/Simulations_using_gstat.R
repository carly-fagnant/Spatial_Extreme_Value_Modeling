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


# only generate as many as I need instad of randomly subsetting
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

ex_sim_data <- readRDS("~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/ex_sim_data.rds")
# sim_data <- readRDS("~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/reg_sim_data.rds")


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
spplot(skrigws)
spplot(skrigws[3]) # ln.scale.pred
spplot(skrigws[5]) # shape.pred

# Need to add in rate parameter somehow

# skrigws_orig <- skrigws
# skrig_orig <- skrig



# Model 3 - Combined Data Series ------------------------------------------

# Issue - we don't have dates for our simulated rainfall values, so cannot really take the max of each day
# Idea - use rate from real data version of the consolidated series, and use random selection of maximum from that?
# i.e. find number of exceedances and randomly select rows to take the maximum for 

# find the rate of exceedance in 40-yr period (1981-2020) in consolidated data
# load consolidated
consol_data <- read.csv("max_precip_regions.csv")

consol_data$Date <- lubridate::ymd(as.character(consol_data$Date)) #put Date column into date format
thresh <- 253


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

fitreg1_alt <- extRemes::fevd(sub$region1_alt, threshold=thresh, type="GP", method="MLE")
fitreg2_alt <- extRemes::fevd(sub$region2_alt, threshold=thresh, type="GP", method="MLE")
fitreg3_alt <- extRemes::fevd(sub$region3_alt, threshold=thresh, type="GP", method="MLE")
fitreg1_alt$rate; fitreg2_alt$rate; fitreg3_alt$rate
mean(c(fitreg1_alt$rate, fitreg2_alt$rate, fitreg3_alt$rate))

mean_rate <- mean(c(fitreg1$rate, fitreg2$rate, fitreg3$rate, fitreg1_alt$rate, fitreg2_alt$rate, fitreg3_alt$rate))
length(sub$Date)
mean_exceed <- round(mean_rate*length(sub$Date)) # 1898

# number of exceedances per station
num_excess <- NULL
for (i in 1:length(sim_data)){
  num_excess[i] <- length(which(sim_data[,i] > 0))
}

max(num_excess)
quantile(num_excess, 0.95) # 2200
quantile(num_excess, 0.99) # 3990
quantile(num_excess, 1)
head(sort(num_excess, decreasing = T), 20)
# stations 169, 128, 201, ...

which(num_excess==6096)
colnames(sim_data)[31]
length(which(num_excess > 2200)) / length(num_excess) # only ~5% of stations have more than 2200 exceedances, we are disregarding their "later" values

# Looking at the first 2200 (1898?) rows of data in sim_data, randomly sample 1898 of them and take the max value for those days
sub_sim_data <- sim_data[1:2200, ]
x <- 1:2200
rows_to_use <- sample(x, size = mean_exceed) # sample 1898 row numbers randomly out of 1:2200

# For those rows, take the maximum value across all the stations (THAT FALL WITHIN EACH REGION)
# region_int from earlier gives the regions that each station falls within, by index

reg1_inds <- which(region_int$REGION==1)
reg2_inds <- which(region_int$REGION==2)
reg3_inds <- which(region_int$REGION==3)

#create function that will take in a row from each region subset and find the max value
get_max_value <- function(values_vector) {
  if (all(is.na(values_vector))) {
    return(NA)
  }
  else {
    return(max(values_vector, na.rm=TRUE))
  }
}

reg1_max <- reg2_max <- reg3_max <- NULL
i = 1
for (r in rows_to_use){
  reg1_max[i] <- get_max_value(sub_sim_data[r, reg1_inds]) # taking max value over a row (of only those stations in region 1)
  reg2_max[i] <- get_max_value(sub_sim_data[r, reg2_inds])
  reg3_max[i] <- get_max_value(sub_sim_data[r, reg3_inds])
  i = i+1
}
# subsetting to region columns each time - could improve code run time by subsetting beforehand

# max_precip_regions <- data.frame(matrix(ncol=3, nrow=length()))
regions_max <- as.data.frame(cbind(reg1_max, reg2_max, reg3_max))
# appending on zeroes
consol_test_data <- rbind(cbind(reg1_max, reg2_max, reg3_max), matrix(data=0, nrow = 14610 - mean_exceed, ncol = 3))
consol_test_data <- as.data.frame(consol_test_data)
colnames(consol_test_data) <- c("region1", "region2", "region3")

fit.scale <- fit.shape <- NULL
fit.reg <- list()
# rate = 1898/14610 = 0.129911
for (reg in 1:3){
  fit <- fit.reg[[reg]] <- extRemes::fevd(consol_test_data[, reg], threshold=thresh, type="GP", method="MLE")
  fit.scale[reg] <- fit$results$par[1]
  fit.shape[reg] <- fit$results$par[2]
}

fit.scale; fit.shape
# getting some outrageous values... may have to remove unrealistic values over a certain amount.
# Maybe look back at all daily raw data I had and set the max value to be the limit

range(stations_sub_df$scale)
range(stations_sub_df$shape)




# Model 1 - Change of Support and CAR -----------------------------------------------------------------

library(spatialreg)

source('~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Test/functionV2.R') # read in updated functions
# this houses the hausMat function to get a matrix of extended Hausdorff distances

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
  }
  else {
    stop(paste("Method of normalization does not match any of the options."))
  }
}

# getting weight matrix
distMat <- hausMat(ws_regs, 0.5)
distMat

#take the inverse first, and then scalar/row normalize
(W <- normalizeMatrix(distMat, "inverse"))
isSymmetric.matrix(W) #symmetric
(W_scalar <- normalizeMatrix(W, "scalar"))
isSymmetric.matrix(W_scalar) #symmetric

hausW <- mat2listw(W_scalar)


car.ml <- spatialreg::spautolm(scale_avg ~1, data=data.frame(t(ws_reg_avg)), family="CAR",
                               listw=mat2listw(W_scalar))
summary(car.ml)$Coef
summary(car.ml)
spdep::moran.test(residuals(car.ml), hausW, zero.policy=T)
spdep::geary.test(residuals(car.ml), hausW, zero.policy=T)


# Test - Can weight matrices have non-zero values on the diagonal? Yes. Not advised, but no error is given
test_mat <- matrix(data = c(1, .9, .2,
                            .9, 1, .75,
                            .2, .75, 1), byrow = T, nrow = 3, ncol = 3)
test_mat <- matrix(data = c(1, 1, .2,
                            1, 1, .75,
                            .2, .75, 1), byrow = T, nrow = 3, ncol = 3)

solve(test_mat)

# If matrix not symmetric?  It will run, but will give warning message
test_mat <- matrix(data = c(1, .7, .2,
                            .9, 1, .75,
                            .2, .8, 1), byrow = T, nrow = 3, ncol = 3)
# if doing KNN approach, some functions (knn2nb, knnHaus) have an argument sym = T to force a matrix to be symmetric

# What if the matrix is not square? - It must be square
test_mat <- matrix(data = c(1, .7, .2,
                            .9, 1, .75), byrow = T, nrow = 2, ncol = 3)
# Error in mat2listw(test_mat) : x must be a square matrix

car.ml <- spatialreg::spautolm(scale_avg ~1, data=data.frame(t(ws_reg_avg)), family="CAR",
                               listw=mat2listw(test_mat))
car.ml <- spatialreg::spautolm(scale_avg ~1, data=data.frame(t(ws_reg_avg)), family="CAR",
                               listw=mat2listw(test_mat), )

# NOTE:
# we want to get updated estimates of the parameters for the 3 regions using the CAR model. 
# Issue- we don't know how to get that from just modeling the 3 values against an intercept.
# Idea- instead find weight matrix for point-to-area (of one point to each area) and then creating a large block diagonal matrix of distances/weights

# Then do CAR modeling where the covariates are indicator variables of the 3 regions (i.e. is a point in that region)
# then the coefficients should make more sense and give the values for the 3 regions


# to do from point to area, need to use extHaus fn instead of hausMat I believe
# hausMat only takes in one type of spatial object at a time, so I don't think it can mix points and polygons
# extHaus can handle any combination!
# but this means I'll have to create my own distance matrix


stations_sub[1, ] # to subset Spatial*DataFrame just add [index, ] to the end
extHaus(stations_sub[1, ], ws_regs[1,], f1=0.5)
class(ws_regs[1,])

extHaus(stations_sub, ws_regs[1,], f1=0.5)
extHaus(stations_sub, ws_regs[2,], f1=0.5)
extHaus(stations_sub, ws_regs[3,], f1=0.5)


extHaus(stations_sub[30,], ws_regs[1,], f1=0.5)
extHaus(stations_sub[30,], ws_regs[2,], f1=0.5)
extHaus(stations_sub[30,], ws_regs[3,], f1=0.5)

######################################
# Practicing possible solution
D <- matrix(data = c(0,   0,  1.6, 1.6,  2.1, 2.1,
                     0,   0,  1.8, 1.8,  2.2, 2.2,
                    1.3, 1.3,  0,   0,   1.8, 1.8,
                    1.8, 1.8,  0,   0,   1.2, 1.2,
                    2.3, 2.3, 1.6, 1.6,   0,   0,
                    2.1, 2.1, 1.5, 1.5,   0,   0), byrow = T, nrow = 6, ncol = 6)

D_alt_0 <- matrix(data = c(.001, .001,  1.6, 1.6,  2.1, 2.1,
                         .001, .001,  1.8, 1.8,  2.2, 2.2,
                         1.3,  1.3, .001, .001, 1.8, 1.8,
                         1.8,  1.8, .001, .001, 1.2, 1.2,
                         2.3,  2.3,  1.6, 1.6, .001, .001,
                         2.1,  2.1,  1.5, 1.5, .001, .001), byrow = T, nrow = 6, ncol = 6)

D_alt <- matrix(data = c(1, 1,  1.6, 1.6,  2.1, 2.1,
                         1, 1,  1.8, 1.8,  2.2, 2.2,
                          1.3,  1.3, 1, 1, 1.8, 1.8,
                          1.8,  1.8, 1, 1, 1.2, 1.2,
                          2.3,  2.3,  1.6, 1.6, 1, 1,
                          2.1,  2.1,  1.5, 1.5, 1, 1), byrow = T, nrow = 6, ncol = 6)
                  
# hMat <- hausMat(ws_regs, f1=0.5)
hMat_miles <- hMat/5280 # convert to miles

D <- matrix(c(rep(0,2), rep(hMat_miles[1,2], 2), rep(hMat_miles[1,3], 2),
              rep(0,2), rep(hMat_miles[1,2], 2), rep(hMat_miles[1,3], 2),
              rep(hMat_miles[1,2], 2), rep(0,2), rep(hMat_miles[2,3], 2),
              rep(hMat_miles[1,2], 2), rep(0,2), rep(hMat_miles[2,3], 2),
              rep(hMat_miles[1,3], 2), rep(hMat_miles[2,3], 2), rep(0,2),
              rep(hMat_miles[1,3], 2), rep(hMat_miles[2,3], 2), rep(0,2)), byrow = T, nrow = 6, ncol = 6)

D_alt <- matrix(c(rep(1,2), rep(hMat_miles[1,2], 2), rep(hMat_miles[1,3], 2),
                  rep(1,2), rep(hMat_miles[1,2], 2), rep(hMat_miles[1,3], 2),
                  rep(hMat_miles[1,2], 2), rep(1,2), rep(hMat_miles[2,3], 2),
                  rep(hMat_miles[1,2], 2), rep(1,2), rep(hMat_miles[2,3], 2),
                  rep(hMat_miles[1,3], 2), rep(hMat_miles[2,3], 2), rep(1,2),
                  rep(hMat_miles[1,3], 2), rep(hMat_miles[2,3], 2), rep(1,2)), byrow = T, nrow = 6, ncol = 6)

D_alt_jitter  <- apply(D_alt, 2, jitter)
D_alt_jitter2 <- apply(D_alt, 1, jitter)

solve(D_alt_jitter)
invD_alt <- 1/D_alt_jitter  # inverse distance
# W_alt <- invD/max(invD)  # scalar normalize

# setting 1's back to 1
fix_D_alt <- D_alt_jitter
fix_D_alt[1,1] <- fix_D_alt[1,2] <- fix_D_alt[2,1] <- fix_D_alt[2,2] <- 
  fix_D_alt[3,3] <- fix_D_alt[3,4] <- fix_D_alt[4,3] <- fix_D_alt[4,4] <- 
  fix_D_alt[5,5] <- fix_D_alt[5,6] <- fix_D_alt[6,5] <- fix_D_alt[6,6] <- 1
solve(fix_D_alt)
inv_fix_D_alt <- 1/fix_D_alt  # inverse distance
# When run on CAR model, it produces NaNs, including log likelihood. Also Lambda is outrageous 

extHaus(stations_sub[29,], ws_regs[1,], f1=0.5)
extHaus(stations_sub[30,], ws_regs[2,], f1=0.5)
extHaus(stations_sub[30,], ws_regs[3,], f1=0.5)

point_to_area <- matrix(nrow=6, ncol=3)
index <- c(4,5, 29,30, 1,2)
for(i in 1:6){
  point_to_area[i, 1] <- extHaus(stations_sub[index[i], ], ws_regs[1, ], f1=0.5)
  point_to_area[i, 2] <- extHaus(stations_sub[index[i], ], ws_regs[2, ], f1=0.5)
  point_to_area[i, 3] <- extHaus(stations_sub[index[i], ], ws_regs[3, ], f1=0.5)
}
point_to_area <- as.data.frame(point_to_area)
point_to_area <- point_to_area/5280 # convert to miles


# D_pa <- matrix(nrow=6, ncol=6)
# for(r in 1:6){
#   for(c in 1:6){
#     D_pa[r, c] <- point_to_area[r, c]
#   }
# }

D_pa <- matrix(c(rep(1,2), rep(point_to_area[1,2], 2), rep(point_to_area[1,3], 2),
                  rep(1,2), rep(point_to_area[2,2], 2), rep(point_to_area[2,3], 2),
                  rep(point_to_area[3,1], 2), rep(1,2), rep(point_to_area[3,3], 2),
                  rep(point_to_area[4,1], 2), rep(1,2), rep(point_to_area[4,3], 2),
                  rep(point_to_area[5,1], 2), rep(point_to_area[5,2], 2), rep(1,2),
                  rep(point_to_area[6,1], 2), rep(point_to_area[6,2], 2), rep(1,2)), byrow = T, nrow = 6, ncol = 6)

D_jitter <- as.data.frame(matrix(lapply(D_pa, jitter), nrow = 6, ncol = 6))
D_jitter2 <- as.data.frame(matrix(apply(D_pa, 2, jitter), nrow = 6, ncol = 6))
D_jitter3 <- apply(D_pa, 2, jitter)
D_jitter4 <- apply(D_pa, 1, jitter)

ex <- extent(ws_regs[1,])
ex
(ex[2]-ex[1] + ex[4]-ex[3])/2  # [1] 235897.4
235897.4/5280  # [1] 44.67754

(ex[2]-ex[1] + ex[4]-ex[3])/4 / 5280

sqrt(((ex[2]-ex[1])/2)^2 + ((ex[4]-ex[3])/2)^2) / 5280 /2

invD <- 1/D_alt  # inverse distance
W_alt <- invD/max(invD)  # scalar normalize

invJ <- 1/D_jitter3  # inverse distance
W_alt <- invD/max(invD)  # scalar normalize

solve(D)
solve(D_alt)
solve(D_jitter3)
solve(W_alt)
det(D_alt)

# make new data
test_dat <- cbind(stations_resids@data$shape[c(4,5, 29,30, 1,2)], c(1,1,0,0,0,0), c(0,0,0,0,1,1))
test_dat <- as.data.frame(test_dat)
colnames(test_dat) <- c("shape", "Reg1", "Reg3")


car_test <- spatialreg::spautolm(shape ~ Reg1 + Reg3, data = test_dat, family="CAR",
                               listw=mat2listw(invJ))
car_test_area <- spatialreg::spautolm(shape ~ Reg1 + Reg3, data = test_dat, family="CAR",
                                 listw=mat2listw(inv_fix_D_alt))
invD_alt; inv_fix_D_alt

listver <- mat2listw(W_alt)

# Jittered version of our distance matrices works, even though not symmetric
# For weight matrix, I just did inverse distance. Could go another step further to normalize?
# I said a distance of 1 for the stations within regions. Try different values
# Should also try for area-to-area version to make it symmetric - Ensor said try Cholesky option to make symmetric


######################################


car.ml.lnscale <- spatialreg::spautolm(log(scale_avg) ~1, data=data.frame(t(ws_reg_avg)), family="CAR",
                               listw=hausW)
summary(car.ml.lnscale)
spdep::moran.test(residuals(car.ml.lnscale), hausW, randomisation = F)
spdep::geary.test(residuals(car.ml.lnscale), hausW)

car.ml.sh <- spatialreg::spautolm(shape_avg ~1, data=data.frame(t(ws_reg_avg)), family="CAR",
                     listw=hausW)
summary(car.ml.sh)
spdep::moran.test(residuals(car.ml.sh), hausW)
spdep::geary.test(residuals(car.ml.sh), hausW)

car.ml.rt <- spatialreg::spautolm(rate_avg ~1, data=data.frame(t(ws_reg_avg)), family="CAR",
                     listw=hausW)
summary(car.ml.rt)
spdep::moran.test(residuals(car.ml.rt), hausW)
spdep::geary.test(residuals(car.ml.rt), hausW)


car.ml <- spatialreg::spautolm(logHV ~ logHHI + TTW + MNR + sqBDH + sqESR,
                               data=stations_sub_df, # need to have data frame at area level...
                               family="CAR", listw=kn4nb_W, # need to set W matrix my own way, using Hausdorff dist
                               zero.policy=T)

library(CARBayes)
S.CARleroux(t(ws_reg_avg[1,])~1, family="gaussian", W=W_scalar, burnin=20, n.sample=50)

S.CARleroux(t(ws_reg_avg[1,])~1, family="gaussian", W=test_mat, burnin=20, n.sample=50)



# Example from spatial short course
car.ml <- spatialreg::spautolm(logHV ~ logHHI + TTW + MNR + sqBDH + sqESR,
                               data=data.frame(dt.trans),
                               family="CAR", listw=kn4nb_W,
                               zero.policy=T)


### Need to think this through...
# we have observations at the station level, and fit rainfall to GPD distribution to get parameter estimates 
# These parameter estimates now become our Z(s)
# these estimates follow the mean distribution of the values for the 3 regions (which one it falls within) + a spatial error and measurement error
