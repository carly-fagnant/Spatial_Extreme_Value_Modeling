source('./Test/function.R')
library(dplyr)

shape <- readOGR('./Test/watershed.shp')
#plot(shape)
### https://epsg.org/crs_2278/NAD83-Texas-South-Central-ftUS.html?sessionkey=m3t6ea1o8k
ws <- spTransform(shape, CRS("+init=epsg:2278"))
plot(ws)

#################### Matrix ######################
hausMatrix <- hausMat(ws, .5, ncores = (detectCores() -1), do.parallel = T)


#################### ToPoint ######################
outline <- readOGR('./Test/watershed_all.shp')
outline <- spTransform(outline, CRS("+init=epsg:2278"))
station_info <- read.csv('./Data/station_info.csv')
station_info <- station_info %>%
  dplyr::mutate(long = LONGITUDE, lat = LATITUDE)
#station_info <- station_info[,c('long', 'lat')]

pointDist <- NULL
for (i in 1:nrow(points)){
  station <- station_info[i,]
  # transforming - long becomes x coordinate and lat is y coordinate
  coordinates(station)=~long+lat
  proj4string(station) <- "+proj=longlat +datum=WGS84"
  station <- spTransform(station, CRS("+init=epsg:2278")) #  projecting data
  pointDist <- rbind(pointDist, pointHaus(station, outline, 0.5, tol = NULL))
}

write.csv(pointDist, file = './Data/station_dist.csv')


#################### Stations in Area ######################
watershed <- read.csv('./Data/watershed.csv')
station_area <- data.frame(X = station_info$X, watershed = NA)
for (i in 1:nrow(points)){
  station <- station_info[i,]
  # transforming - long becomes x coordinate and lat is y coordinate
  coordinates(station)=~long+lat
  proj4string(station) <- "+proj=longlat +datum=WGS84"
  station <- spTransform(station, CRS("+init=epsg:2278")) #  projecting data
  
  if (st_intersects(st_as_sf(station), st_as_sf(outline), sparse = F)){
    area <- st_as_sf(ws)
    for (j in 1:nrow(area)){
      a <- area[j,]
      if (st_intersects(st_as_sf(station), st_as_sf(a), sparse = F)){
        name <- a$WTSHNAME
        station_area[i, 'watershed'] <- as.character(watershed[watershed$WTSHNAME == name,2])
        break
      }
    }
  }
}

station_area$watershed <- as.factor(station_area$watershed)
write.csv(station_area, './Data/station_area.csv', row.names = FALSE)


