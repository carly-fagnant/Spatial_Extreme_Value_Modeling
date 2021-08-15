library(dplyr)
library(geoR)
library(gstat)
library(rgdal)
library(raster)
library(sp)

# Notes -----------------------------------------------------------------------
# likfitBGCM: 
# this is a new function and still in draft format and pretty much untested

# Example code ----------------------------------------------------------------
## Testing code: Inference in the
## Bivariate Gaussian Common Component Model 
##
## 1. Simulating Y1 and Y2 with some (but not all) locations in common
##
## 1.1. Coordinates
set.seed(30)
all.coords <- round(cbind(runif(200), runif(200)), dig=3)
set.seed(111); ind1 <- sample(1:200, 120) 
coords1 <- all.coords[ind1,]
set.seed(222); ind2 <- sample(1:200, 80) 
coords2 <- all.coords[ind2,]

##
## Model parameters:
## m1 = 10, mu2 = 50
## sigma01 = 2, sigma1 = 1.5, sigma02 = 8, sigma2 = 6
## phi0 = 0.2, phi1=0.15, phi2=0.25
mu1 <- 10; mu2 <- 50
sigma01 <- 2; sigma1 <- 1.5; sigma02 <- 8; sigma2 <- 6
phi0 <- 0.2; phi1 <- 0.15; phi2 <- 0.25
##
## Model parameters (reparametrised)
## sigma = 2, nu1 = 0.75, eta = 4, sigma2 = 3
sigma <- sigma01
nu1 <- sigma1/sigma
eta <- sigma02/sigma; nu2 <- sigma2/sigma
##
## Simulating model components
##
S0 <- grf(grid=all.coords, cov.pars=c(sigma, phi0))$data
S1 <- grf(grid=coords1, cov.pars=c(sigma1, phi1))$data
S2 <- grf(grid=coords2, cov.pars=c(sigma2, phi2))$data
##
## Y1 and Y2 data
##
y1 <- as.geodata(cbind(coords1,10+S0[ind1]+S1))
y2 <- as.geodata(cbind(coords2,50+S0[ind2]+S2))
##
## maximum likelihood estimation
##
fit12 <- likfitBGCCM(y1, y2, ini.s=c(2,1.5,8,6), ini.phi=c(.2,.15,.25), control=list(trace=T))
fit12
fit12a <- likfitBGCCM(y1, y2, ini.s=c(2,1.5,8,6), ini.phi=c(.2,.15,.25), method="BFGS")
fit12a
fit12b <- likfitBGCCM(y1, y2, ini.s=c(2,1.5,8,6), ini.phi=c(.2,.15,.25), method="L-BFGS-B", control=list(trace=T))
fit12b
fit12c <- likfitBGCCM(y1, y2, ini.s=c(2,1.5,8,6), ini.phi=c(.2,.15,.25), fc.min="nlminb")
fit12c
##
## Prediction
##
locs <- cbind(c(0.2, 0.5, 0.7, 0.2), c(0.3, 0.5, 0.2, 0.8))
pred1 <- predict(fit12, loc=locs)
pred1[1:2]
pred2 <- predict(fit12, loc=locs, var=2)
pred2[1:2]

# Read in watershed regions ---------------------------------------------------
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

# Read station data -----------------------------------------------------------

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

# read parameter values
ws_mat <- readRDS("ws_mat.rds")
ws_stat <- ws_mat[1, ]
ws_scale <- ws_mat[2, ]
ws_shape <- ws_mat[3, ]
ws_rate <- ws_mat[4, ]

stations <- read.csv("station_info.csv")

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

# Fit model -------------------------------------------------------------------
proj4string(stations_sub) <- proj4string(ws_regs)

stations_sub_df <- readRDS(file = "stations_sub_df.rds")

submore <- stations_sub_df %>%
  dplyr::mutate(ln.scale = log(scale)) %>%
  dplyr::select(STAT_NO, ln.scale, shape, rate, long, lat)

extremesub <- submore %>%
  dplyr::select(ln.scale, shape, long, lat)

sub_ln.scale <- submore %>%
  dplyr::select(ln.scale, long, lat)

sub_shape <- submore %>%
  dplyr::select(shape, long, lat)

### Function to project data given a data frame with data and coordinates
project_data <- function(df){
  coordinates(df)=~long+lat
  proj4string(df) <- "+proj=longlat +datum=WGS84"
  df <- spTransform(df, CRS("+init=epsg:2278")) # projecting data
  is.projected(df)
  return(df)
}


# get scale and shape data from df and convert it into a geodata object
sub_stat <- project_data(submore)
extremesub <- project_data(extremesub)
sub_shape <- project_data(sub_shape)
obj1 <- as.geodata(extremesub)
obj2 <- as.geodata(sub_shape)

fit <- likfitBGCCM(obj1, obj2)


# Predict ---------------------------------------------------------------------
grid <- sp::makegrid(ws_regs, cellsize = 2000) # cellsize in map units
grid <- SpatialPoints(grid, proj4string = CRS(proj4string(ws_regs)))
grid <- grid[ws_regs, ] # grid points only within the polygon
plot(ws_regs)
plot(grid, pch = ".", add = T)

prediction1 <- predict(fit, grid@coords, variable.to.predict = 1)
# list of 3: 
# predict: predicted ln.scale at every coordinate
# krige.var: prediction variance at every coordinate

prediction2 <- predict(fit, grid@coords, variable.to.predict = 2)
# list of 3: 
# predict: predicted shape at every coordinate
# krige.var: prediction variance at every coordinate
