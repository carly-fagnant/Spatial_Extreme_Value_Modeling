#######################
# Testing the models on real data
#
# by Carly Fagnant
#######################
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


### Function to project data given a data frame with data and coordinates
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
# setwd("~/Documents/HCFCD Watershed Data") # wherever you have the following file saved 
## (access "window_1day_dclust_updated.rds" from zipped file in Box - should be 1.67 GB unzipped)
# window <- readRDS("window_1day_dclust_updated.rds")  # list of lists (82 x 601) of GPD fits


# Model 1 - Random Effects Model using CAR and median Hausdorff distance --------

## this was already run in "Simulations_using_gstat.R", we will bring in the saved estimates here
# car_fit <- readRDS(file = "~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/car_fit.rds")

# Finding return levels given the parameters:
car_fit

# extRemes::return.level - calculates return levels given an fevd object
# extRemes::rlevd() - calculates return levels given a distribution and parameter estimates

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


# for(reg in 1:3){
#   scale <- 
#   shape <- 
#   rate <- 
# }



# Model 2 - Kriging and Aggregation ---------------------------------------

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

pars_to_rl(krig_fit, 25)/254
pars_to_rl(krig_fit, 100)/254
pars_to_rl(krig_fit, 500)/254



# Model 3 - Consolidated Series -------------------------------------------

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

###############

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



# Extra Kriging/Variogram Code --------------------------------------------

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
