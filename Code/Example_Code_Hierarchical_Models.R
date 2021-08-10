#######################
# Hierarchical Models
#
# by Carly Fagnant
#######################

library(rgdal)
library(sp)
library(gstat)
library(dplyr)

####### 3 Hydroregions #########
setwd("~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Test")
wd <- getwd()

ws_regs <- rgdal::readOGR(dsn = paste0(wd, "/watershed_region_updated"), layer = "watershed_region_updated")
plot(ws_regs)
# plot(ws_regs[1,]) # to access the regions separately
# plot(ws_regs[2,])
# plot(ws_regs[3,])

# Alternate version of 3 regions
ws_reg_alt <- rgdal::readOGR(dsn = paste0(wd, "/watershed_region_alt"), layer = "watershed_region_alt")
plot(ws_reg_alt)

# repeat whatever you do for the first version for this version as well, labeling all results with "_alt" at the end


####### Point-referenced data #########
setwd("~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data")
# load data
station_info <- read.csv("station_info.csv")

# add new location (long/lat) columns that we will transform to coordinates
station_info <- station_info %>%
  dplyr::mutate(long = LONGITUDE, lat = LATITUDE)

# long becomes x coordinate and lat is y coordinate
coordinates(station_info)=~long+lat
proj4string(station_info) <- "+proj=longlat +datum=WGS84"
station_info <- spTransform(station_info, CRS("+init=epsg:2278")) # projecting data
is.projected(station_info)
# it should now be in the proper coordinates and in the same projection as the watersheds

# set regions and points to have the same projection
proj4string(ws_regs) <- proj4string(station_info)


# Finding which stations fall within the 3 regions ------------------------

region_intersect <- sp::over(station_info, ws_regs)
# region_intersect lists the 601 stations in order and the integer of what region each station falls within (NA if not in any region)

region_intersect$REGION <- as.factor(region_intersect$REGION) # need these to be factors in order for the model matrix and summary to work
summary(region_intersect$REGION) # number of stations in each region (NAs are those not in any)

h <- model.matrix( ~ REGION - 1, data=region_intersect)

# rownames(h)  # to get station indices (might need to do as.integer or as.numeric)
inds <- as.integer(rownames(h))
# this gives the numbers of the stations that fall within the 3 regions



# Uploading parameter data and subsetting to stations within the 3 regions --------

stations <- read.csv("station_info.csv")

# setwd("~/Documents/HCFCD Watershed Data") # wherever you have the following file saved 
# (access "window_1day_dclust_updated.rds" from zipped file in Box - should be 1.67 GB unzipped)
window <- readRDS("window_1day_dclust_updated.rds")  # list of lists (82 x 601) of GPD fits

# Accessing Model Fits and saving the parameter values (scale, shape, rate) for each station
ws_stat <- stations[,1]
ws_scale <- NULL
ws_shape <- NULL
ws_rate  <- NULL
j = 1
for(i in ws_stat){
  fit <- window[[82]][[i]] # just the last 40-year window
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
stations <- stations %>%
  dplyr::mutate(scale = ws_scale, shape = ws_shape, rate = ws_rate, long = LONGITUDE, lat = LATITUDE)

# subset to those 166 stations that fall within the 3 regions
# ... but also remove stations with NAs if any
stations_sub <- stations %>%
  dplyr::filter(STAT_NO %in% inds) %>%
  dplyr::filter(!is.na(scale))

# long becomes x coordinate and lat is y coordinate
coordinates(stations_sub)=~long+lat
proj4string(stations_sub) <- "+proj=longlat +datum=WGS84"
stations_sub <- spTransform(stations_sub, CRS("+init=epsg:2278")) # projecting data
is.projected(stations_sub)
# it should now be in the proper coordinates and in the same projection as the watersheds


# Plotting regions and the stations on top --------------------------------

### Plot to make pretty... color the regions diff colors and label them Hydrologic Regions 1, 2, & 3
plot(ws_regs)
plot(stations_sub, add=T, pch=20, cex=0.5)



# Trying variogram modeling, kriging, and simulations using gstat package  --------

### Try prediction/simulation

# Example from gstat package with meuse data
data(meuse)
data(meuse.grid)
coordinates(meuse) <- ~x+y
proj4string(meuse) <- CRS("+init=epsg:28992")
coordinates(meuse.grid) = ~x+y
proj4string(meuse.grid) <- CRS("+init=epsg:28992")
gridded(meuse.grid) = TRUE
spplot(meuse.grid)

# creating cross-variogram
meuse.g <- gstat(id = "log-zn", formula = log(zinc)~1, data = meuse) #, locations = ~x+y
meuse.g <- gstat(meuse.g, id = "log-cu", log(copper)~1, meuse) #, ~x+y

meuse.g <- gstat(meuse.g, model=vgm(1, "Sph", 900, nugget = 1), fill.all=T)
x <- variogram(meuse.g, cutoff=1000)
meuse.fit = fit.lmc(x, meuse.g)
plot(x, model = meuse.fit)

out = gstat.cv(meuse.fit, nmax = 40, nfold = 5)
summary(out)
out = gstat.cv(meuse.fit, nmax = 40, nfold = c(rep(1,100), rep(2,55)))
summary(out)
# mean error, ideally 0:
mean(out$residual)
# MSPE, ideally small
mean(out$residual^2)
# Mean square normalized error, ideally close to 1
mean(out$zscore^2)
# correlation observed and predicted, ideally 1
cor(out$observed, out$observed - out$residual)
# correlation predicted and residual, ideally 0
cor(out$observed - out$residual, out$residual)

meuse.cok <- predict(meuse.fit, newdata = meuse.grid) # cokriging?
plot(meuse.cok)

sims <- predict(meuse.fit, newdata = meuse.grid, nsim=1) ### this is the line that would not run for me...
# sims <- krige(meuse.fit, newdata = meuse.grid, nsim=5)
spplot(sims)


### Now try on our data
# I ran the following on the old data through 2017... skip commented code to continue with new data

# ## Loading in OLD Station Data for Variograms
# setwd("~/Documents/HCFCD Watershed Data")
# wd <- getwd()
# stations <- read.csv("station_info.csv")
# 
# window <- readRDS("window_1day_dclust.rds")
# ws_stat <- stations[,1]
# ws_scale <- NULL
# ws_shape <- NULL
# ws_rate  <- NULL
# j = 1
# for(i in ws_stat){
#   fit <- window[[79]][[i]]
#   if(is.na(fit)){
#     ws_scale[j] <- ws_shape[j] <- ws_rate[j]  <- NA
#   }else{
#     ws_scale[j] <- fit$results$par[1]
#     ws_shape[j] <- fit$results$par[2]
#     ws_rate[j]  <- fit$rate
#   }
#   j = j+1
# }
# 
# # add new location (long/lat) columns that we will transform to coordinates
# stations <- stations %>%
#   dplyr::mutate(scale = ws_scale, shape = ws_shape, rate = ws_rate, long = LONGITUDE, lat = LATITUDE)
# 
# # subset to those 166 stations that fall within the 3 regions
# # ... but also remove NAs # length(which(is.na(stations_sub$scale)))
# # 149 stations
# stations_sub <- stations %>%
#   dplyr::filter(STAT_NO %in% inds) %>%
#   dplyr::filter(!is.na(scale))
# 
# # long becomes x coordinate and lat is y coordinate
# coordinates(stations_sub)=~long+lat
# proj4string(stations_sub) <- "+proj=longlat +datum=WGS84"
# stations_sub <- spTransform(stations_sub, CRS("+init=epsg:2278")) # projecting data
# is.projected(stations_sub)
# # it should now be in the proper coordinates and in the same projection as the watersheds

v = variogram(log(scale)~1, data = stations_sub)
v.fit = fit.variogram(v, vgm(c("Exp", "Mat", "Sph", "Gau")))
v.fit
plot(v)
plot(v, model=v.fit, main="Log(scale)")
# plot(v.fit)

v2 = variogram(shape~1, data = stations_sub)
v2.fit = fit.variogram(v2, vgm(c("Exp", "Mat", "Sph", "Gau")))
v2.fit
plot(v2, model=v2.fit, main="shape")

v2 = variogram(log(shape)~1, data = stations_sub)
v2.fit = fit.variogram(v2, vgm(c("Exp", "Mat", "Sph", "Gau")))
# v2.fit = fit.variogram(v2, vgm(c("Exp", "Mat", "Sph", "Gau"), nugget=.01))
# v2.fit = fit.variogram(v2, vgm(nugget=.005))
v2.fit
plot(v2, model=v2.fit, main="Log(shape)")

v3 = variogram(log(rate)~1, data = stations_sub)
v3.fit = fit.variogram(v3, vgm(c("Exp", "Mat", "Sph", "Gau")))
v3.fit
plot(v3)
plot(v3, model=v3.fit, main="Log(rate)")


### Log(scale)
# Gau
# nugget: 0
# sill: 0.241171
# range: 6944.668 ft ~ 1.3 miles

### shape
# Gau
# nugget: 0
# sill: 0.2100937
# range: 5755.497 ft ~ 1.09 miles

### Log(rate)
# Sph
# nugget: 0.1923552
# sill: 0.1107377 + 0.1923552 = 0.3030929
# range: 13361.08 ft ~ 2.5 miles

### Looking for correlation between parameters - do they need to be modeled jointly?
plot(stations_sub$scale, stations_sub$shape)
cor(stations_sub$scale, stations_sub$shape)

plot(log(stations_sub$scale), stations_sub$shape)
cor(log(stations_sub$scale), stations_sub$shape)

plot(log(stations_sub$scale), stations_sub$rate)
plot(stations_sub$shape, stations_sub$rate)
cor(log(stations_sub$scale), stations_sub$rate)
cor(stations_sub$shape, stations_sub$rate)

plot(log(stations_sub$scale), log(stations_sub$rate))
plot(stations_sub$shape, log(stations_sub$rate))
cor(log(stations_sub$scale), log(stations_sub$rate))
cor(stations_sub$shape, log(stations_sub$rate))

# outlier - #31, station 169 (< 5 yrs of data)
which(stations_sub$rate > 0.4)
which(stations_sub$shape > 3)
cor(log(stations_sub$scale[-31]), stations_sub$rate[-31])  # These are the ones to go off of
cor(stations_sub$shape[-31], stations_sub$rate[-31])
cor(log(stations_sub$scale[-31]), stations_sub$shape[-31])

cor(log(stations_sub$scale[-31]), log(stations_sub$rate[-31]))
cor(stations_sub$shape[-31], log(stations_sub$rate[-31]))

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
# cv.fit = fit.lmc(vg, g, vgm("Gau"))
# cv.fit = fit.lmc(vg, g, vgm("Exp"))
# cv.fit = fit.lmc(vg, g, vgm("Mat"))
plot(vg, model = cv.fit)





###### Model 1 - Change of Support ######
# already created COS matrix H (it is h model matrix above)
# need to...
#     - figure out covariance - both measurement error and spatial 
#     - CAR modeling using median Hausdorff distance matrix

# Sigma_Z (measurement error) could be pulled from the standard errors of each parameter estimate
# Sigma_Y (spatial error) needs to model covariance using variogram?

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




