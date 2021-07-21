# Same code as combine_new_regions.R, but with Sims Bayou in region 2 instead of 3

#######################
# Combining shapefiles to fix/update the 3 hydrologic regions
#
# by Carly Fagnant
#######################


library(sf)
library(dplyr)

setwd("~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/All_Watershed_Shapefiles")
# setwd("~/Documents/HCFCD Watershed Data/All_Watershed_Shapefiles")
wd <- getwd()
folders <- list.files(wd)



# Code from combine.R - by Yue Zhuo ---------------------------------------

watershed <- NULL

# Hydrologic Region Code
one <- c("SPRING CREEK", "WILLOW CREEK", "BARKER RESERVOIR", "CYPRESS CREEK",
         "ADDICKS RESERVOIR")
two <- c("BRAYS BAYOU", "SAN JACINTO RIVER", "HUNTING BAYOU", "GREENS BAYOU",
         "LUCE BAYOU", "BUFFALO BAYOU", "WHITE OAK BAYOU", "SIMS BAYOU")
three <- c("CLEAR CREEK", "ARMAND BAYOU", "SAN JACINTO & GALVESTON BAY",
           "VINCE BAYOU", "CARPENTERS BAYOU", "SPRING GULLY & GOOSE CREEK", "CEDAR BAYOU",
           "JACKSON BAYOU")

for (dir in folders[-c(4,10, 19)]){
  water <- sf::st_read(paste0(wd, "/", dir, "/", dir, "_Watershed.shp"), stringsAsFactors = FALSE)
  water <- select(water, WTSHNAME, geometry)
  
  if (water$WTSHNAME[1] %in% one){
    code <- 1
  } else if (water$WTSHNAME[1] %in% two){
    code <- 2
  } else if (water$WTSHNAME[1] %in% three){
    code <- 3
  }
  
  water$REGION <- code
  watershed <- rbind(watershed, water)
}

# white oak bayou (WTSHNAME missing)
water <- st_read(paste0(wd, "/", folders[4], "/", folders[4], "_Watershed.shp"), stringsAsFactors = FALSE)
water <- select(water, WTSHNAME, geometry)
if (water$WTSHNAME[1] %in% one){
  code <- 1
} else if (water$WTSHNAME[1] %in% two){
  code <- 2
} else if (water$WTSHNAME[1] %in% three){
  code <- 3
}
water$REGION <- code
watershed <- rbind(watershed, water)

# cypress creek
cyp <- st_read(paste0(wd, "/", folders[10], "/", folders[10], "_Watershed.shp"), stringsAsFactors = FALSE)
cyp_sub <- select(cyp, WTSHNAME, geometry)
cyp <- st_union(cyp_sub)
cyp <- st_as_sf(data.frame("WTSHNAME" = "CYPRESS CREEK", cyp))

# buffalo bayou
baf <- st_read(paste0(wd, "/", folders[20], "/", folders[20], "_Watershed.shp"), stringsAsFactors = FALSE)
baf <- select(baf, WTSHNAME, geometry)

# barker reservoir
bar <- st_read(paste0(wd, "/", folders[18], "/", folders[18], "_Watershed.shp"), stringsAsFactors = FALSE)
bar <- select(bar, WTSHNAME, geometry)

# addicks reservoir (filename does not match)
add <- st_read(paste0(wd, "/", folders[19], "/S_Subbasins.shp"), stringsAsFactors = FALSE)
add <- st_union(add)
add <- st_as_sf(data.frame("WTSHNAME" = "ADDICKS RESERVOIR", add))
st_crs(add) <- st_crs(watershed)

pivot <- rbind(cyp, add)
pivot <- st_union(pivot)
for (i in length(pivot[[1]]):2){
  pivot[[1]][[i]] <- NULL
}
add <- st_difference(pivot, cyp)
add <- st_as_sf(data.frame("WTSHNAME" = "ADDICKS RESERVOIR", add))

pivot <- rbind(baf, add)
pivot <- st_union(pivot)
for (i in length(pivot[[1]]):2){
  pivot[[1]][[i]] <- NULL
}
add <- st_difference(pivot, baf)
add <- st_as_sf(data.frame("WTSHNAME" = "ADDICKS RESERVOIR", add))

pivot <- rbind(bar, add)
pivot <- st_union(pivot)
for (i in length(pivot[[1]]):2){
  pivot[[1]][[i]] <- NULL
}
add <- st_difference(pivot, bar)
add <- st_as_sf(data.frame("WTSHNAME" = "ADDICKS RESERVOIR", add))

pivot <- rbind(select(water, WTSHNAME, geometry), add)
pivot <- st_union(pivot)
for (i in length(pivot[[1]]):2){
  pivot[[1]][[i]] <- NULL
}
pivot[[1]][[1]] <- pivot[[1]][[1]][1]
add <- st_difference(pivot, select(water, WTSHNAME, geometry))
add <- st_as_sf(data.frame("WTSHNAME" = "ADDICKS RESERVOIR", add))

if (cyp$WTSHNAME[1] %in% one){
  code <- 1
} else if (cyp$WTSHNAME[1] %in% two){
  code <- 2
} else if (cyp$WTSHNAME[1] %in% three){
  code <- 3
}
cyp_sub$REGION <- code
watershed <- rbind(watershed, cyp_sub)

if (add$WTSHNAME[1] %in% one){
  code <- 1
} else if (add$WTSHNAME[1] %in% two){
  code <- 2
} else if (add$WTSHNAME[1] %in% three){
  code <- 3
}
add$REGION <- code
watershed <- rbind(watershed, add)

plot(watershed)




# Code to load in split San Jacinto River watershed (G) and update/correct the 3 regions - by Carly --------
setwd("~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/G_San_Jacinto_River_split")
# setwd("~/Documents/HCFCD Watershed Data/G_San_Jacinto_River_split")
wd <- getwd()

water  <- st_read(paste0(wd, "/G_San_Jacinto_River_Watershed_split_north.shp"), stringsAsFactors = FALSE)
water2 <- st_read(paste0(wd, "/G_San_Jacinto_River_Watershed_split_south.shp"), stringsAsFactors = FALSE)
water <- select(water, WTSHNAME, geometry)
water2 <- select(water2, WTSHNAME, geometry)

water$REGION <- 2
water2$REGION <- 3
watershed <- rbind(watershed, water, water2)

# reorder watersheds, removing San Jacinto River (G) and adding the north and south split of it
watershed_reg <- rbind(watershed[1:4,], watershed[19,], 
                       watershed[5,], watershed[7:9,], watershed[20:21,],
                       watershed[10:17,], watershed[22,],
                       watershed[18,], watershed[23:24,])

regionall <- NULL
for (i in 1:3){
  region <- watershed_reg %>% filter(REGION==i) %>% select(REGION, geometry)
  region <- st_as_sf(data.frame("REGION" = i, st_union(region)))
  regionall <- rbind(regionall, region)
}


### Exploring the makeup of regionall - Need to delete leftover lines from the union ###
plot(regionall[[1]]) # int 1 2 3
plot(regionall[[2]]) # list of 3
# the 3 regions   class(regionall[[2]])   "sfc_POLYGON" "sfc"  (list)

class(regionall[[2]][[1]]) # list of 2
# [1] "XY"      "POLYGON"   "sfg"    
class(regionall[[2]][[2]]) # list of 2
# [1] "XY"      "POLYGON"   "sfg"         
class(regionall[[2]][[3]]) # list of 2
# [1] "XY"      "POLYGON"   "sfg" 

class(regionall[[2]][[1]][[1]]) # matrix 14829 x 2
class(regionall[[2]][[1]][[2]]) # matrix 13 x 2

length(regionall[[2]][[2]]) # 197
class(regionall[[2]][[2]][[1]]) # matrix 16521 x 2
class(regionall[[2]][[2]][[2]]) # matrix 7 x 2
class(regionall[[2]][[2]][[3]]) # matrix 5 x 2
# ...

plot(regionall)
plot(regionall[[2]])
plot(regionall[[2]][[1]])
plot(regionall[[2]][[1]][[1]])
plot(regionall[[2]][[1]][[2]]) # delete this

plot(regionall[[2]][[2]]) # list of length 197
length(regionall[[2]][[2]])
plot(regionall[[2]][[2]][[1]])
plot(regionall[[2]][[2]][[197]]) # delete all except regionall[[2]][[2]][[1]]

plot(regionall[[2]][[3]]) # list of length 695
length(regionall[[2]][[3]])
plot(regionall[[2]][[3]][[1]])
plot(regionall[[2]][[3]][[2]]) # delete all except regionall[[2]][[3]][[1]]
plot(regionall[[2]][[3]][[695]])


### Deleting the extra lines ###
regionall[[2]][[1]][[2]] <- NULL
for (i in length(regionall[[2]][[2]]):2){ # have to index backwards in order to actually remove objects from the list
  regionall[[2]][[2]][[i]] <- NULL
}
for (i in length(regionall[[2]][[3]]):2){
  regionall[[2]][[3]][[i]] <- NULL
}
plot(regionall, axes=T)

## A few unwanted lines leftover in main boundaries - two in region 2 and one in region 3
dim(regionall[[2]][[2]][[1]]) # region 2 coordinates -  matrix 16521 x 2
dim(regionall[[2]][[3]][[1]]) # region 3 coordinates -  matrix 4406 x 2

getwd()
### SAVE THE ALTERNATE REGION ORGANIZATION
# Setting new place to save new region data
setwd("~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Test/watershed_region_alternate_temp")
#setwd("~/Documents/HCFCD Watershed Data")
st_write(regionall, "watershed_region_alternate_temp.shp", delete_layer = T)


# For overall boundary, see combine.R - this leaves unwanted lines
# st_write(st_union(regionall), "watershed_all_temp.shp", delete_layer = T)

ws_alternate_temp <- rgdal::readOGR("watershed_region_alternate_temp.shp")
plot(ws_alternate_temp)

# Clean up in ArcGIS
setwd("~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Test/watershed_region_alt")
ws_alt_final <- rgdal::readOGR("watershed_region_alt.shp")
plot(ws_alt_final)
