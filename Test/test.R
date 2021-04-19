source('./Test/function.R')

shape <- readOGR('./Test/watershed_all.shp')
#plot(shape)
### https://epsg.org/crs_2278/NAD83-Texas-South-Central-ftUS.html?sessionkey=m3t6ea1o8k
ws <- spTransform(shape, CRS("+init=epsg:2278"))
plot(ws)

#################### Matrix ######################
hausMatrix <- hausMat(ws, .5, ncores = (detectCores() -1), do.parallel = T)


#################### ToPoint ######################
points <- read.csv('./Data/watershed_stations.csv')
point <- SpatialPoints(points[1,c('lat', 'lon')], proj4string=CRS("+init=epsg:2278"))
pointHaus(point, ws[1,], 0.5, tol = NULL)
