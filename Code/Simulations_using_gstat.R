#######################
# Simulation of new data
#
# by Carly Fagnant
#######################

# Using the data we have, capture the mean and covariance structure.
# The mean structure will come in at the region level - we will assume the same mean value for all points that fall within a region
# The covariance structure comes in at the point level... at least that is all we can get from our raw data.

# maybe gstat::vgmArea could provide covariance at the area level?

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

# Step 2 ------------------------------------------------------------------
# Step 2: fit a cross-variogram in gstat (to this simulated data?) and use the predict fn with nsim = 1 to simulate parameter values at all of the stations at once
# This brings in more of the covariance structure? - Yes, the spatial side of it

# Note: maybe if I can access covariance matrix from gstat cross-variogram, I can combine Steps 1 and 2 into just a multivariate normal. 
#   However, I run into the issue of that I am doing it separately for each region


g <- gstat(id = "ln.scale", formula = log(scale)~1, data = stations_sub)
g <- gstat(g, id = "shape", formula = shape~1, data = stations_sub)
# examine variograms and cross variogram:
vg <- variogram(g)
## "Sph", "Gau", "Exp", "Mat"
cv.fit1 = fit.lmc(vg, g, vgm("Sph"), correct.diagonal = 1.01)
cv.fit = fit.lmc(vg, g, vgm("Mat"))
plot(vg, model = cv.fit)
cv.fit
plot(vg, model = cv.fit1)
cv.fit1

cv.fit$set=list(nocheck=1)
cv.fit1$set=list(nocheck=1)
pred <- predict(cv.fit1, newdata = stations_sub, nsim = 1)
pred <- predict(cv.fit, newdata = stations_sub, nsim = 1)
plot(ws_regs)
plot(pred, pch = "•", add = T)
spplot(pred) 





# Model 2 - Kriging and Aggregation ---------------------------------------

# Note that one may take a shortcut in coding by setting the newdata argument in gstat::predict to be the regions of interest
# this uses spsample 

