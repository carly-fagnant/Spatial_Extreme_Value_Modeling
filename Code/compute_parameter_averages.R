library(rgdal)
library(sp)
library(gstat)
library(dplyr)

memory.limit(24000)

# Loading watersheds ----------------------------------------------------------
wd <- getwd()
ws_regs <- rgdal::readOGR(
  dsn = paste0(wd, "/watershed_region_updated"),
  layer = "watershed_region_updated"
)
plot(ws_regs)

ws_reg_alt <- rgdal::readOGR(
  dsn = paste0(wd, "/watershed_region_alt"),
  layer = "watershed_region_alt"
)
plot(ws_reg_alt)

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


# Finding which stations fall within the 3 regions ----------------------------

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

# window <- readRDS(file = "window_1day_dclust_updated.rds")
# Accessing Model Fits and saving the parameter values (scale, shape, rate) for each station
# ws_stat <- stations[,1]
# ws_scale <- NULL
# ws_shape <- NULL
# ws_rate  <- NULL
# j <- 1
# for(i in ws_stat){
  #fit <- window[[82]][[i]] # just the last 40-year window
  #if(is.na(fit)){
    #ws_scale[j] <- ws_shape[j] <- ws_rate[j]  <- NA
  #} else {
    #ws_scale[j] <- fit$results$par[1]
    #ws_shape[j] <- fit$results$par[2]
    #ws_rate[j]  <- fit$rate
  #}
  #j <- j + 1
#}
ws_mat <- rbind(ws_stat, ws_scale, ws_shape, ws_rate)
saveRDS(ws_mat, file = 'ws_mat.rds') 
# store data into CRS to avoid having to load the massive 
# window_1day_dclust_updated file into RStudio.

ws_mat <- readRDS("ws_mat.rds")
ws_stat <- ws_mat[1, ]
ws_scale <- ws_mat[2, ]
ws_shape <- ws_mat[3, ]
ws_rate <- ws_mat[4, ]

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

# A function for computing the average scale, shape and rate parameter
# of each region
mean_scale_shape_rate <- function(stations_subset, regions) {
  num_regions <- length(regions@polygons)
  num_stats <- length(stations_sub)
  scale_avg <- rep(0, num_regions)
  shape_avg <- rep(0, num_regions)
  rate_avg <- rep(0, num_regions)
  # For each region, look for the stations in the regions and note their
  # parameter values to subsequently compute the region average.
  for (i in 1:num_regions) {
    print(paste0("Outer iteration: ", i))
    ws_reg <- regions[i, ]
    counter <- 0
    for (stat in 1:num_stats) {
      stat <- stations_subset[stat, ]
      if (is.null(rgeos::gDifference(stat, ws_reg))) {
        counter <- counter + 1
        scale_avg[i] <- scale_avg[i] + stat@data$scale
        shape_avg[i] <- shape_avg[i] + stat@data$shape
        rate_avg[i] <- rate_avg[i] + stat@data$rate
      }
    }
    if (counter != 0) {
      scale_avg[i] <- scale_avg[i] / counter
      shape_avg[i] <- shape_avg[i] / counter
      rate_avg[i] <- rate_avg[i] / counter
    }
  }
  scale_shape_rate_mat <- rbind(scale_avg, shape_avg, rate_avg)
  return(scale_shape_rate_mat)
}
proj4string(ws_regs) <- proj4string(stations_sub)
proj4string(ws_regs_alt) <- proj4string(stations_sub)

# Compute averages and save the results
# ws_reg_avg <- mean_scale_shape_rate(stations_sub, ws_regs)
# saveRDS(ws_reg_avg, file = 'ws_reg_avg.rds')
# ws_reg_alt_avg <- mean_scale_shape_rate(stations_sub, ws_reg_alt)
# saveRDS(ws_reg_alt_avg, file = 'ws_reg_alt_avg.rds')
# 
ws_reg_avg <- readRDS("ws_reg_avg.rds")
ws_reg_alt_avg <- readRDS("ws_reg_alt_avg.rds")

beta0 <- ws_reg_avg[, 2]
beta1 <- ws_reg_avg[, 1] - beta0
beta2 <- ws_reg_avg[, 3] - beta0

# let's do some manual fitting of two direct variograms and a cross variogram
proj4string(stations_sub) <- proj4string(ws_regs)
g <- gstat(id = "ln.scale", formula = log(scale)~1, data = stations_sub)
g <- gstat(g, id = "shape", formula = shape~1, data = stations_sub)
# examine variograms and cross variogram:
vg <- variogram(g)
plot(vg)
## "Sph", "Gau", "Exp", "Mat"
cv.fit <- fit.lmc(vg, g, vgm("Sph"))
# cv.fit = fit.lmc(vg, g, vgm("Gau"))
# cv.fit = fit.lmc(vg, g, vgm("Exp"))
# cv.fit = fit.lmc(vg, g, vgm("Mat"))
plot(vg, model = cv.fit)
cok <- predict(g, newdata = stations_sub) # cokriging?
plot(cok)
cv.fit$set=list(nocheck=1)
ws_reg_grid <- predict(cv.fit, newdata = ws_regs, nsim = 1)
# If newdata is of class SpatialPolygons or SpatialPolygonsDataFrame, 
# calculated, then the block average for each of the polygons or polygon sets 
# is using spsample to discretize the polygon(s). Argument sps.args controls 
# the  parameters used for spsample. The "location" with respect to which 
# neighbourhood selection is done is for each polygon the SpatialPolygons 
# polygon label point;

# Compare the results:
rbind(ws_reg_avg, log(ws_reg_avg[1, ]))
ws_reg_grid@data

# Ensure both stations_sub and ws_reg_alt have the same coordinate system
proj4string(ws_reg_alt) <- proj4string(stations_sub)
proj4string(stations_sub) <- proj4string(ws_reg_alt)
ws_reg_alt_grid <- predict(cv.fit, newdata = ws_reg_alt, nsim = 1)

# Compare the results:
rbind(ws_reg_alt_avg, log(ws_reg_alt_avg[1, ]))
ws_reg_alt_grid@data
