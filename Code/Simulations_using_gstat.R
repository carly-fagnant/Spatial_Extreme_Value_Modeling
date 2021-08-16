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

# Step ?: use the extRemes::revd fn to simulate daily rainfall data using the parameters given
# NOTE: No input for rate here. Will outputs be all above the threshold? I think so...
#       Will need to correct for this if so, by randomly removing observations at the rate estimated



# setwd("~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Test")
# wd <- getwd()



# Testing Functions -------------------------------------------------------

### Discoveries:
# - revd function to generate GPD realizations does NOT create any that are less than the threshold. (Makes sense, since we only fit to those above the threshold.)
#   This means that if we want our simulated data to look more like real rainfall data, we will need to do one of the following:
#     - randomly select a proportion of data to remove (or set to zero) or replace with uniformly distributed values < threshold
#     - subtract a set amount (quantile value?) from each observation --> NO, this would change the parameter values and return levels
#  Recall that we assumed the observations were independent (after declustering) before fitting, so do not let multiple days in a row have values above the threshold?
#  We will decluster our simulated data, but want to make sure we have the proper proportion of data left to work with

# NOTE: the avg rate of exceedance of the threshold is only 0.0543096, so about 95% of the data simulated should be under the threshold...

ws_reg_avg <- readRDS("~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/ws_reg_avg.rds")
avg_sc <- mean(ws_reg_avg[1,])
avg_sh <- mean(ws_reg_avg[2,])
avg_rt <- mean(ws_reg_avg[3,])

simtry <- extRemes::revd(n = 10000, type = "GP", scale = avg_sc, shape = avg_sh, threshold = 253) # use avg param values
any(simtry < 253)



# Step 1 ------------------------------------------------------------------
# Step 1: use multivariate normal (separately for each region) to simulate (shape and scale) parameter values for stations within each region
# But what covariance structure do we use here? Same for all regions? (yes, I believe so)




# Model 2 - Kriging and Aggregation ---------------------------------------

# Note that one may take a shortcut in coding by setting the newdata argument in gstat::predict to be the regions of interest
# this uses spsample 

