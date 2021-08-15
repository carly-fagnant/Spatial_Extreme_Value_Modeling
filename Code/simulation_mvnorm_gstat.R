library(dplyr)
library(gstat)
library(mvtnorm)
library(raster)
library(rgdal)
library(sp)

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


# Simulations at the stations -------------------------------------------------

ws_mat <- readRDS("ws_mat.rds")
ws_stat <- ws_mat[1, ]
ws_scale <- ws_mat[2, ]
ws_shape <- ws_mat[3, ]
ws_rate <- ws_mat[4, ]

var_scale <- var(ws_scale, na.rm = TRUE)
var_shape <- var(ws_shape, na.rm = TRUE)
cov <- cov(ws_scale, ws_shape, use = "pairwise.complete.obs")
sigma_mat <- matrix(c(var_scale, cov, cov, var_shape), nrow = 2, byrow = T)

stations <- read.csv("station_info.csv")

stations <- stations %>%
  dplyr::mutate(scale = ws_scale, shape = ws_shape, rate = ws_rate, long = LONGITUDE, lat = LATITUDE)

# subset to those 166 stations that fall within the 3 regions
# ... but also remove stations with NAs if any
stations_sub <- stations %>%
  dplyr::filter(STAT_NO %in% inds) %>%
  dplyr::filter(!is.na(scale))

coordinates(stations_sub)=~long+lat
proj4string(stations_sub) <- "+proj=longlat +datum=WGS84"
stations_sub <- spTransform(stations_sub, CRS("+init=epsg:2278")) # projecting data
is.projected(stations_sub)

ws_reg_avg <- readRDS("ws_reg_avg.rds")
ws_reg_alt_avg <- readRDS("ws_reg_alt_avg.rds")

proj4string(stations_sub) <- proj4string(ws_regs)
g <- gstat(id = "ln.scale", formula = log(scale)~1, data = stations_sub)
g <- gstat(g, id = "shape", formula = shape~1, data = stations_sub)
# examine variograms and cross variogram:
vg <- variogram(g)
plot(vg)
## "Sph", "Gau", "Exp", "Mat"
cv.fit <- fit.lmc(vg, g, vgm("Sph"))
cv.fit1 = fit.lmc(vg, g, vgm("Sph"), correct.diagonal = 1.01)

cv.fit$set=list(nocheck=1)
cv.fit1$set=list(nocheck=1)

n <- 1000

prediction <- predict(cv.fit, newdata = stations_sub, nsim = n)
preditction1 <- predict(cv.fit1, newdata = stations_sub, nsim = n)

n <- 1000
num_stations <- 1
prediction_region1 <- rmvnorm(
  n * num_stations,
  mean = c(ws_reg_avg[1, 1], ws_reg_avg[2, 1]),
  sigma = sigma_mat
)
prediction_region1[1, ] <- log(prediction_region1[1, ])
num_stations <- 1
prediction_region2 <- rmvnorm(
  n * num_stations,
  mean = c(ws_reg_avg[1, 2], ws_reg_avg[2, 2]),
  sigma = sigma_mat
)
prediction_region2[1, ] <- log(prediction_region2[1, ])
num_stations <- 1
prediction_region2 <- rmvnorm(
  n * num_stations,
  mean = c(ws_reg_avg[1, 3], ws_reg_avg[2, 3]),
  sigma = sigma_mat
)
prediction_region3[1, ] <- log(prediction_region3[1, ])
