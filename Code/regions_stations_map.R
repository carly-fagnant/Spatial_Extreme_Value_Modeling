## Code to plot stations that are within one of the three hydrologic regions on a map including the regions
library(rgdal)
library(sp)
library(dplyr)
library(leaflet)

# Prepare data for plotting -----------------------------------------------
regions <- rgdal::readOGR('~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Test/watershed_region_updated')
proj4string(regions) <- sp::CRS("+init=epsg:2278")
regions <- spTransform(regions, CRS("+proj=longlat +datum=WGS84")) #Leaflet uses this projection

regions_alt <- readOGR('~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Test/watershed_region_alt')
regions_alt <- spTransform(regions_alt, CRS("+proj=longlat +datum=WGS84"))

## Read in stations' information
station_info <- read.csv("~/Documents/GitHub/scrape-NOAA-data/Reference_files/station_info.csv")
# add new location (long/lat) columns that we will transform to coordinates
station_info <- station_info %>%
  dplyr::mutate(long = LONGITUDE, lat = LATITUDE)
# long becomes x coordinate and lat is y coordinate
coordinates(station_info)=~long+lat
proj4string(station_info) <- "+proj=longlat +datum=WGS84"
station_info <- spTransform(station_info, CRS(proj4string(station_info))) # projecting data

# Find the stations that are within one of the three regions
# the same station points data could be used, but this is in case it is decided to color points by region
station_info$region <- sp::over(station_info, regions)$REGION
station_info$region_alt <- sp::over(station_info, regions_alt)$REGION
stations <- subset(station_info, !is.na(region))
stations_alt <- subset(station_info, !is.na(region_alt))


# Plot the data on a map --------------------------------------------------
## create numerical labels to overlay directly upon the regions
region_labels <- data.frame(x = c(-95.8, -95.3, -95),
                            y = c(29.9, 29.85, 29.8),
                            num = c(1, 2, 3))
region_labels <- region_labels %>%
  dplyr::mutate(long = x, lat = y)
coordinates(region_labels) <- ~long+lat
proj4string(region_labels) <- "+proj=longlat +datum=WGS84"
region_labels <- spTransform(region_labels, CRS(proj4string(region_labels)))

## original hydroregions
leaflet(options = leafletOptions(zoomControl = FALSE, attributionControl=FALSE)) %>%
  addTiles() %>% 
  setView(lng = -95.45, lat = 29.9, zoom = 9) %>%
  addPolygons(data = regions, opacity = .8, weight = 1.5, fillOpacity = 0.5, 
              fillColor = c("blue", "red", "green"), color = "black") %>%
  addCircleMarkers(data = stations, lat = ~ LATITUDE, lng = ~ LONGITUDE, opacity = 1,
                   radius = 2.5, weight = 1, color = "black", fillColor = "gold", fillOpacity = 0.8) %>%
  addLegend(colors = c("blue", "red", "green"), labels = c("1", "2", "3"), title = "REGION",
            opacity = 0.6, position = "bottomleft") #%>%
  # addLabelOnlyMarkers(data = region_labels, lng = ~x, lat = ~y, label =  ~as.character(num), labelOptions = 
  #                       labelOptions(noHide = T, direction = 'top', textOnly = T,textsize = "20px", opacity = 1,
  #                                    style = list('color' = 'white')))

## alternative regions setup
leaflet(options = leafletOptions(zoomControl = FALSE, attributionControl=FALSE)) %>%
  addTiles() %>% 
  setView(lng = -95.45, lat = 29.9, zoom = 9) %>%
  addPolygons(data = regions_alt, opacity = .8, weight = 1.5, fillOpacity = 0.5, 
              fillColor = c("blue", "red", "green"), color = "black") %>%
  addCircleMarkers(data = stations_alt, lat = ~ LATITUDE, lng = ~ LONGITUDE, opacity = 1,
                   radius = 2.5, weight = 1, color = "black", fillColor = "gold", fillOpacity = 0.8) %>%
  addLegend(colors = c("blue", "red", "green"), labels = c("1", "2", "3"), title = "REGION",
            opacity = 0.6, position = "bottomleft")
