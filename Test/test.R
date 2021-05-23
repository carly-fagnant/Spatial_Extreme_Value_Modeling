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
watershed <- read.csv('./Data/watershed_code.csv')
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
write.csv(station_area, './Data/station_area.csv', na = "", row.names = FALSE)


#################### AreaToPoint ######################
region <- readOGR('./Test/watershed_region.shp')
region <- spTransform(region, CRS("+init=epsg:2278"))
station_info <- read.csv('./Data/station_info.csv')
station_info <- station_info %>%
  dplyr::mutate(long = LONGITUDE, lat = LATITUDE)
#station_info <- station_info[,c('long', 'lat')]

pointAreaDist <- NULL
for (i in 1:nrow(station_info)){
  station <- station_info[i,]
  # transforming - long becomes x coordinate and lat is y coordinate
  coordinates(station)=~long+lat
  proj4string(station) <- "+proj=longlat +datum=WGS84"
  station <- spTransform(station, CRS("+init=epsg:2278")) #  projecting data
  pivot <- NULL
  for (j in 1:3){
    r <- region[j,]
    pivot <- cbind(pivot, pointHaus(station, r, 0.5, tol = NULL))
  }
  pointAreaDist <- rbind(pointAreaDist, pivot)
}
colnames(pointAreaDist) <- c(1,2,3)
write.csv(pointAreaDist, file = './Data/stationToArea_dist.csv')


#################### Zipcode ######################
shape <- readOGR('./Data/Zip_Codes/Zip_Codes.shp')
zip <- spTransform(shape, CRS("+init=epsg:2278"))
plot(zip)
plot(ws, add=T, col = 'red')
zip <- st_as_sf(zip)
### https://cohgis-mycity.opendata.arcgis.com/datasets/f392021d9d2344938b0958909d690cc7_0?geometry=-98.402%2C28.930%2C-92.409%2C30.599

site_info <- read.csv('./Data/houston_sites.csv')
site_info <- site_info %>%
  dplyr::mutate(long = longitude, lat = latitude)

site_zip <- data.frame(site_code = site_info$site_code, 
                        Site_Name = site_info$Site_Name, 
                        ZIP = NA)
for (i in 1:nrow(site_zip)){
  site <- site_info[i,]
  # transforming - long becomes x coordinate and lat is y coordinate
  coordinates(site)=~long+lat
  proj4string(site) <- "+proj=longlat +datum=WGS84"
  site <- spTransform(site, CRS("+init=epsg:2278")) #  projecting data
  
  for (j in 1:nrow(zip)){
    z <- zip[j,]
    if (st_intersects(st_as_sf(site), st_as_sf(z), sparse = F)){
      site_zip[i, 'ZIP'] <- z$ZIP_CODE
      break
    }
  }
}

#site_zip$ZIP <- as.factor(site_zip$ZIP)
write.csv(site_zip, './Data/site_zip.csv', na = "", row.names = FALSE)

