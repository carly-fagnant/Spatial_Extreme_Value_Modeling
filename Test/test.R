source('./Test/function.R')

shape <- readOGR('./watershed.shp')
#plot(shape)
ws <- spTransform(shape, CRS("+init=epsg:2278"))
plot(ws[1,])

#################### Matrix ######################
hausMatrix <- hausMat(ws[1:3,], .5, ncores = (detectCores() -1), do.parallel = T)


#################### ToPoint ######################
points <- read.csv('./Data/watershed_stations.csv')
point <- SpatialPoints(points[1,c('lat', 'lon')], proj4string=CRS("+init=epsg:2278"))
pointHaus(point, ws[1,], 0.5, tol = NULL)
