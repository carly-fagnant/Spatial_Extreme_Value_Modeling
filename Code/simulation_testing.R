# trying to get simulations working

library(rgdal)
library(sp)
library(dplyr)
library(gstat)

#setwd("~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data")


# Reading in data ---------------------------------------------------------
stations_sub <- readRDS("~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/stations_sub.rds")
coordinates(stations_sub)=~long+lat
proj4string(stations_sub) <- "+proj=longlat +datum=WGS84"
stations_sub <- spTransform(stations_sub, CRS("+init=epsg:2278")) # projecting data
is.projected(stations_sub)

regions <- readOGR('~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Test/watershed_region_updated')
proj4string(regions) <- CRS("+init=epsg:2278")
regions <- spTransform(regions, CRS("+init=epsg:2278"))

ws_reg_avg <- readRDS("~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/ws_reg_avg.rds")


# 1.  using station locations ---------------------------------------------
# taking code from the example of simulation
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

proj4string(ws_regs) <- proj4string(stations_sub)
cv.fit$set=list(nocheck=1) # added
try.cok <- predict(cv.fit, newdata = ws_regs) # cokriging?
plot(try.cok)
ws_reg_grid <- predict(cv.fit, newdata = ws_regs, nsim = 1)

plot(vg, model = cv.fit1)
cv.fit1
cv.fit1$set=list(nocheck=1) # added
try.cok <- predict(cv.fit1, newdata = ws_regs) # cokriging?
plot(try.cok)

##
cv.fit$set=list(nocheck=1)
test <- predict(cv.fit, newdata = stations_sub, nsim = 30)
plot(ws_regs)
plot(test, pch = "•", add = T)


# 2.  simulating data first -----------------------------------------------
## 3 univariate
stations_sub$region <- sp::over(stations_sub, regions)$REGION
mu_r1 <- mean(subset(stations_sub, region==1)$scale) #check if these match with the scale_avg in ws_reg_avg
sd_r1 <- sd(subset(stations_sub, region==1)$scale)

mu_r2 <- mean(subset(stations_sub, region==2)$scale)
sd_r2 <- sd(subset(stations_sub, region==2)$scale)

mu_r3 <- mean(subset(stations_sub, region==3)$scale)
sd_r3 <- sd(subset(stations_sub, region==3)$scale)

#not randomized
# scale_sim_results <- ifelse(stations_sub$region==1, yes=rnorm(1, mu_r1, sd_r1), no=
#                               ifelse(stations_sub$region==2, yes=rnorm(1, mu_r2, sd_r2), no=
#                                        ifelse(stations_sub$region==3, yes=rnorm(1, mu_r3, sd_r3), no=NA)))

#randomized, negative values are possible...
scale_sim_results <- rep(NA, nrow(stations_sub))
for (i in 1:nrow(stations_sub)) {
  current_region <- stations_sub$region[i]
  current_result <- ifelse(current_region==1, yes=rnorm(1, mu_r1, sd_r1), no=
                             ifelse(current_region==2, yes=rnorm(1, mu_r2, sd_r2), no=
                                      ifelse(current_region==3, yes=rnorm(1, mu_r3, sd_r3), no=NA)))
  scale_sim_results[i] <- current_result
}

## use multivariate normal, negative values included
?mvrnorm
full_scale_sim_results <- mvrnorm(n=nrow(stations_sub), mu=ws_reg_avg[1,], Sigma=cov(ws_reg_avg))
scale_station_sim_results <- rep(NA, nrow(stations_sub))
for (i in 1:nrow(stations_sub)) {
  #each station will have its own simulated result row, then select the number for the region it belongs to
  scale_station_sim_results[i] <- full_scale_sim_results[i, stations_sub$region[i]]
}

stations_sub$sim_scale <- scale_station_sim_results

# taking code from the example of simulation
v = variogram(log(sim_scale)~1, data = stations_sub)
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
#g <- gstat(g, id = "shape", formula = shape~1, data = stations_sub)
# examine variograms and cross variogram:
plot(variogram(g))
vg <- variogram(g)
## "Sph", "Gau", "Exp", "Mat"
# cv.fit = fit.lmc(vg, g, vgm("Sph"))
cv.fit1 = fit.lmc(vg, g, vgm("Sph"), correct.diagonal = 1.01)
# cv.fit = fit.lmc(vg, g, vgm("Gau"))
# cv.fit = fit.lmc(vg, g, vgm("Exp"))
cv.fit = fit.lmc(vg, g, vgm("Mat"))

cv.fit$set=list(nocheck=1)
pred <- predict(cv.fit1, newdata = stations_sub, nsim = 30)
plot(ws_regs)
plot(pred, pch = "•", add = T)
spplot(pred) 

##

grid <- sp::makegrid(ws_regs, cellsize = 2000) # cellsize in map units
grid <- SpatialPoints(grid, proj4string = CRS(proj4string(ws_regs)))
grid <- grid[ws_regs, ] # grid points only within the polygon
plot(ws_regs)
plot(grid, pch = ".", add = T)
pred <- predict(cv.fit1, newdata = grid) # cokriging?
spplot(pred[1]) 


#### repeat for shape
full_shape_sim_results <- mvrnorm(n=nrow(stations_sub), mu=ws_reg_avg[2,], Sigma=cov(ws_reg_avg)) #values are too large...
shape_station_sim_results <- rep(NA, nrow(stations_sub))
for (i in 1:nrow(stations_sub)) {
  #each station will have its own simulated result row, then select the number for the region it belongs to
  shape_station_sim_results[i] <- full_shape_sim_results[i, stations_sub$region[i]]
}

stations_sub$sim_shape <- shape_station_sim_results

v = variogram(sim_shape~1, data = stations_sub)
v.fit = fit.variogram(v, vgm(c("Exp", "Mat", "Sph", "Gau")))
v.fit
plot(v)
plot(v, model=v.fit, main="shape")

