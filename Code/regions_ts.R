# Referencing Example_Code_consolidated_series.R

# Setup: Libraries and Data -----------------------------------------------
library(dplyr)
library(rgdal)
library(sp)
library(xts)

## Read in precipitation data
rain_data <- read.csv("~/Documents/GitHub/scrape-NOAA-data/Data/all_stations_precip_UDP_updated_2021.csv")

## Read in stations' information
station_info <- read.csv("~/Documents/GitHub/scrape-NOAA-data/Reference_files/station_info.csv")
# add new location (long/lat) columns that we will transform to coordinates
station_info <- station_info %>%
  dplyr::mutate(long = LONGITUDE, lat = LATITUDE)
# long becomes x coordinate and lat is y coordinate
coordinates(station_info)=~long+lat
proj4string(station_info) <- "+proj=longlat +datum=WGS84"
station_info <- spTransform(station_info, CRS("+init=epsg:2278")) # projecting data
is.projected(station_info)

## Read in region shapefiles, set them to the coordinates and in the same projection as the watersheds
regions <- readOGR('~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Test/watershed_region_updated')
#proj4string(regions) <- CRS("+init=epsg:2278")
#regions <- spTransform(regions, CRS("+init=epsg:2278"))
proj4string(regions) <- proj4string(station_info)
# it should now be in the proper coordinates and in the same projection as the watersheds
plot(regions)

regions_alt <- readOGR('~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Test/watershed_region_alt')
#regions_alt <- spTransform(regions_alt, CRS("+init=epsg:2278"))
proj4string(regions_alt) <- proj4string(station_info)
plot(regions_alt)


# Extreme Value Series by Regions -----------------------------------------
## Finding intersection of points and regions
# sp::over will let us get what region each station falls within (or NA if none) in order of the 601 stations
region_intersect <- sp::over(station_info, regions)
class(region_intersect$REGION)
region_intersect_alt <- sp::over(station_info, regions_alt)

## Aggregate station info into data by each region
#create function that will retrieve relevant stations' columns for subsetting purposes, given a region
num_stations <- 601
get_columns_index <- function(intersect_vector, region) {
  stations <- seq(1, num_stations, 1)[intersect_vector == region] #get station number for those in "region", there will be NA's
  col_nums <- stations[!is.na(stations)] + 1 #add 1 to column indices since column 1 is Date
  col_nums <- c(1, col_nums) #include Date column
  return(col_nums)
}
#make 6 subsets of the data to easily collect a max value from each date for each region subset (region1_alt technically a duplicate)
rain_region1 <- rain_data[, get_columns_index(region_intersect, 1)]
rain_region2 <- rain_data[, get_columns_index(region_intersect, 2)]
rain_region3 <- rain_data[, get_columns_index(region_intersect, 3)]
rain_region1_alt <- rain_data[, get_columns_index(region_intersect_alt, 1)]
rain_region2_alt <- rain_data[, get_columns_index(region_intersect_alt, 2)]
rain_region3_alt <- rain_data[, get_columns_index(region_intersect_alt, 3)]

## Populate consolidated time series
#create final data frame to populate
max_precip_regions <- data.frame(matrix(ncol=7, nrow=nrow(rain_data)))
colnames(max_precip_regions) <- c("Date", "region1", "region2", "region3", "region1_alt", "region2_alt", "region3_alt")
max_precip_regions$Date <- rain_data$Date

#create function that will take in a row from each region subset and find the max value
get_max_value <- function(values_vector) {
  if (all(is.na(values_vector))) {
    return(NA)
  }
  else {
    return(max(values_vector, na.rm=TRUE))
  }
}
#loop through every date in the data set and get the max value for each region subset for that day
for (i in 1:nrow(rain_data)) {
  max_precip_regions$region1[i] <- get_max_value(rain_region1[i, -1]) 
  max_precip_regions$region2[i] <- get_max_value(rain_region2[i, -1]) 
  max_precip_regions$region3[i] <- get_max_value(rain_region3[i, -1]) 
  max_precip_regions$region1_alt[i] <- get_max_value(rain_region1_alt[i, -1]) 
  max_precip_regions$region2_alt[i] <- get_max_value(rain_region2_alt[i, -1]) 
  max_precip_regions$region3_alt[i] <- get_max_value(rain_region3_alt[i, -1]) 
}

write.csv(max_precip_regions, "Data/max_precip_regions.csv")


# Plotting the Series -----------------------------------------------------
regions_max_ts <- read.csv("Data/max_precip_regions.csv")
regions_max_ts <- regions_max_ts[, -1] #remove column of row numbers

regions_max_ts$Date <- as.Date(as.character(regions_max_ts$Date), "%Y%m%d")
par(mfrow=c(3,2))
#order such that the 3 regions from the 2 setups have their own column
for (i in c(2,5,3,6,4,7)) {
  current_xts <- xts::xts(regions_max_ts[, i], order.by = regions_max_ts[, 1])
  print(plot(current_xts, main=colnames(regions_max_ts)[i]))
}
