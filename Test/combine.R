library(sf)
library(dplyr)
#setwd("~/Spatial_Extreme_Value_Modeling/Data/All_Watershed_Shapefiles")
wd <- "./Data/All_Watershed_Shapefiles"
folders <- list.files(wd)

catchment <- NULL
watershed <- NULL

for (dir in folders[-c(4,19)]){
  catch <- st_read(paste0(wd, "/", dir, "/", dir, "_Catchment.shp"), stringsAsFactors = FALSE)
  catch <- select(catch, WTSHNAME, geometry)
  catchment <- rbind(catchment, catch)
  water <- st_read(paste0(wd, "/", dir, "/", dir, "_Watershed.shp"), stringsAsFactors = FALSE)
  water <- select(water, WTSHNAME, geometry)
  watershed <- rbind(watershed, water)
}

# white oak bayou (WTSHNAME missing)
catch <- st_read(paste0(wd, "/", folders[4], "/", folders[4], "_Catchment.shp"), stringsAsFactors = FALSE)
names(catch)[5] <- "WTSHNAME"
catch$WTSHNAME <- rep("WHITE OAK BAYOU", length(catch$WTSHNAME))
catch <- select(catch, WTSHNAME, geometry)
catchment <- rbind(catchment, catch)
water <- st_read(paste0(wd, "/", folders[4], "/", folders[4], "_Watershed.shp"), stringsAsFactors = FALSE)
water <- select(water, WTSHNAME, geometry)
watershed <- rbind(watershed, water)

# addicks reservoir (filename does not match)
catch <- st_read(paste0(wd, "/", folders[19], "/Addicks_Modified_CAP.shp"), stringsAsFactors = FALSE)
names(catch)[8] <- "WTSHNAME"
catch$WTSHNAME <- rep("ADDICKS RESERVOIR", length(catch$WTSHNAME))
catch <- select(catch, WTSHNAME, geometry)
# catchment <- rbind(catchment, catch)
water <- st_read(paste0(wd, "/", folders[19], "/S_Subbasins.shp"), stringsAsFactors = FALSE)
water <- st_union(water)
water <- st_as_sf(data.frame("WTSHNAME" = "ADDICKS RESERVOIR", water))
st_crs(water) <- st_crs(watershed)
watershed <- rbind(watershed, water)

plot(catchment)
plot(watershed)

st_write(watershed, "~/Spatial_Extreme_Value_Modeling/watershed.shp", delete_layer = T)


