#######################
# RandomFields package vs gstat for bivariate modeling and simulation
#
# by Carly Fagnant
#######################

# RandomFields seems to have more difficult syntax, 
# and might be best only for simulation... having difficulty getting a (cross-)variogram.
# Though it might be that my graphics are messed up and not plotting correctly

### Note: I include a lot of "dev.new()" to view the plots because they were not showing up in my Rstudio plot frame
# You may or may not have the same problem

library(RandomFields)
library(gstat)
library(rgdal)
library(sp)
library(dplyr)


RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
## RFoptions(seed=NA) to make them all random again


### Examples for RandomFields

# Simulation of a bivariate linear model of coregionalization
M1 <- c(0.9, 0.6)
M2 <- c(sqrt(0.19), 0.8)
model <- RandomFields::RMmatrix(M = M1, RMwhittle(nu = 0.3)) + RMmatrix(M = M2, RMwhittle(nu = 2))
x <- y <- seq(-10, 10, 0.2)
simu <- RandomFields::RFsimulate(model, x, y)
dev.new()
plot(simu)
# dev.off()

## first example: bivariate Linear Model of Coregionalisation
x <- y <- seq(0, 10, 0.2)
model1 <- RMmatrix(M = c(0.9, 0.43), RMwhittle(nu = 0.3)) +
  RMmatrix(M = c(0.6, 0.8), RMwhittle(nu = 2))
plot(model1)
simu1 <- RFsimulate(RPdirect(model1), x, y)
dev.new()
plot(simu1)

# Simulation of a bivariate Whittle-Matern model
model <- RandomFields::RMbiwm(nudiag = c(1, 2), nured = 1, rhored = 1, cdiag = c(1, 5), s = c(1, 1, 2))
simu <- RandomFields::RFsimulate(model, x, y)

dev.new()
plot(simu)
# dev.off()

model <- RMbiwm(nudiag=c(0.3, 2), nured=1, rhored=1, cdiag=c(1, 1.5),
                s=c(1, 1, 2))
dev.new()
plot(model)
dev.new()
plot(RFsimulate(model, x, y))

model <- RMexp(var=1.6, scale=0.5) + RMnugget(var=0) #exponential + nugget
dev.new()
plot(model)

dev.off()

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


# RandomFields ------------------------------------------------------------

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


# Trying different data to see what works for variogram
sub_stat <- project_data(submore)
extremesub <- project_data(extremesub)
sub_shape <- project_data(sub_shape) # only got single variable to work to be able to see plot

RFstations <- sp2RF(sub_stat, param=list(n=1, vdim=4))
RFstations <- sp2RF(extremesub, param=list(n=1, vdim=2))
RFstations <- sp2RF(sub_shape, param=list(n=1, vdim=1))

vario <- RandomFields::RFvariogram(data = RFstations)
dev.new()
plot(vario, model=RMmatern())
plot(vario)
# not working for me :(



# gstat -------------------------------------------------------------------

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

g <- gstat(id = "ln.scale", formula = log(scale)~1, data = stations_sub)
g <- gstat(g, id = "shape", formula = shape~1, data = stations_sub)
# examine variograms and cross variogram:
plot(variogram(g))
vg <- variogram(g)
## "Sph", "Gau", "Exp", "Mat"
# cv.fit = fit.lmc(vg, g, vgm("Sph"))
cv.fit1 = fit.lmc(vg, g, vgm("Sph"), correct.diagonal = 1.01)
# cv.fit = fit.lmc(vg, g, vgm("Gau"))
# cv.fit = fit.lmc(vg, g, vgm("Exp"))
cv.fit = fit.lmc(vg, g, vgm("Mat"))
plot(vg, model = cv.fit)
cv.fit

cv.fit$set=list(nocheck=1) # added
try.cok <- predict(cv.fit, newdata = ws_regs) # cokriging?
plot(try.cok)

plot(vg, model = cv.fit1)
cv.fit1
cv.fit1$set=list(nocheck=1) # added
try.cok <- predict(cv.fit1, newdata = ws_regs) # cokriging?
plot(try.cok)

par(mfrow=c(1,1))

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

# How to plot the results of cokriging
plot(try.cok) # this only plots the first variable listed...
spplot(try.cok) # this plots all variables... then can subset to those wanted
spplot(try.cok[c(1,3)])
spplot(try.cok[1]) # coarser color legend, but easier to see differences
spplot(try.cok[3])
plot(try.cok[1]) # smoother color transitions, but almost looks blurred
plot(try.cok[3])
spplot(try.cok[5]) # covariance between variables. Does not look great, very concentrated around the station locations

gridded(grid) = TRUE # turn into spatial pixels data frame

