###############################
# Example Code - Creating Consolidated Series
#
# by Carly Fagnant
###############################

# Recall PRCP = Precipitation is measured as "tenths of mm"   (254 = 1 inch)
# To get inches: x/254
# To get mm: x/10

### Goal: create one overall precipitation data series for each hydrologic region (hydroregion)
#         by taking the maximum across all daily values for stations falling within each region

#  Do this for both versions of the 3 hydrogregions in Test folder (watershed_region_updated & watershed_region_alt)

library(rgdal)
library(sp)
library(dplyr)
library(lubridate)
library(imputeTS)
library(xts)


# Step 1 - See which stations fall within each region ---------------------
#   I used function sp::over to find intersect of stations and regions, see example below

setwd("~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Test")
wd <- getwd()

ws_regs <- rgdal::readOGR(dsn = paste0(wd, "/watershed_region_updated"), layer = "watershed_region_updated")
class(ws_regs)
plot(ws_regs)
plot(ws_regs[1,])
plot(ws_regs[2,])
plot(ws_regs[3,])

# load station location data
setwd("~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data")
station_info <- read.csv("station_info.csv")

# add new location (long/lat) columns that we will transform to coordinates
station_info <- station_info %>%
  dplyr::mutate(long = LONGITUDE, lat = LATITUDE)

# long becomes x coordinate and lat is y coordinate
coordinates(station_info)=~long+lat
proj4string(station_info) <- "+proj=longlat +datum=WGS84"
station_info <- spTransform(station_info, CRS("+init=epsg:2278")) # projecting data
is.projected(station_info)
proj4string(ws_regs) <- proj4string(station_info)
# it should now be in the proper coordinates and in the same projection as the watersheds

### Finding intersection of points and regions ###
region_intersect <- sp::over(station_info, ws_regs)
class(region_intersect$REGION)
# region_intersect lists the 601 stations in order and the integer of what region each station falls within (NA if not in any region)


# Step 2 - Create a consolidated series -----------------------------------
# The example below is using different data... basically I had selected a few nearby stations to combine
# You will do something similar, but will use all stations within a region, so will have more series to combine at once

# Going into our other GitHub repository where we had the rainfall data...
# Mine uses the old data, but you will use the updated one you created
setwd("~/Documents/GitHub/scrape-NOAA-data/Data")
raw <- read.csv("all_stations_precip_UDP.csv", header = TRUE)
raw <- raw[,-1]   # delete first column since is repeat of row names (numbers)

# Function to subset to specifc station
raw_data <- function(raw, station_no){
  Galv_dat <- raw[,c(1, station_no+1)] # offset by 1 because 1st column is date
  Galv_dat[,1] <- lubridate::ymd(as.character(Galv_dat[,1])) #put Date column into date format
  return(Galv_dat)
}

# Function to plot time series
plot_ts <- function(data, station_no){
  ts <- xts::xts(data[,2], order.by = data[,1])
  plot(ts, main = station_no)
}


### Now create series of max Galveston values between all stations
galv_stats <- c(131, 138, 150, 156, 531, 589, 591) # these are the station numbers I want to combine

Galv131 <- raw_data(raw, 131)
Galv138 <- raw_data(raw, 138)
Galv150 <- raw_data(raw, 150)
Galv156 <- raw_data(raw, 156)
Galv531 <- raw_data(raw, 531)
Galv589 <- raw_data(raw, 589)
Galv591 <- raw_data(raw, 591)

### Combining stations data into one series using maximum (or NA if all NA)
# Manually make value column since the join put NAs when there were not
combo_val <- NULL
for(i in 1:length(Galv131[,1])){
  vals <- c(Galv131[i,2], Galv138[i,2], Galv150[i,2], Galv156[i,2], Galv531[i,2], Galv589[i,2], Galv591[i,2])
  if(all(is.na(vals))){
    combo_val[i] <- NA
  }else{
    combo_val[i] <- max(imputeTS::na_remove(vals))
  }
}

Galv_max <- as.data.frame(cbind(Galv131[,1], combo_val)) 
colnames(Galv_max) <- c("Date", "Value")
Galv_max$Date <- as.Date(Galv_max$Date)



### Just plotting the consolidated time series 
Galv_max_ts <- xts::xts(Galv_max[,2], order.by = Galv_max[,1])
plot(Galv_max_ts)

