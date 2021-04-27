library(sf)
library(dplyr)
#setwd("~/Spatial_Extreme_Value_Modeling/Data/All_Watershed_Shapefiles")
wd <- "./Data/All_Watershed_Shapefiles"
folders <- list.files(wd)

catchment <- NULL
watershed <- NULL

# Hydrologic Region Code
one <- c("SPRING CREEK", "WILLOW CREEK", "BARKER RESERVOIR", "CYPRESS CREEK",
         "ADDICKS RESERVOIR")
two <- c("BRAYS BAYOU", "SAN JACINTO RIVER", "HUNTING BAYOU", "GREENS BAYOU",
         "LUCE BAYOU", "BUFFALO BAYOU", "WHITE OAK BAYOU")
three <- c("CLEAR CREEK", "ARMAND BAYOU", "SIMS BAYOU", "SAN JACINTO & GALVESTON BAY",
           "VINCE BAYOU", "CARPENTERS BAYOU", "SPRING GULLY & GOOSE CREEK", "CEDAR BAYOU",
           "JACKSON BAYOU")

for (dir in folders[-c(4,10, 19)]){
  # catch <- st_read(paste0(wd, "/", dir, "/", dir, "_Catchment.shp"), stringsAsFactors = FALSE)
  # catch <- select(catch, WTSHNAME, geometry)
  # catchment <- rbind(catchment, catch)
  water <- st_read(paste0(wd, "/", dir, "/", dir, "_Watershed.shp"), stringsAsFactors = FALSE)
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
# catch <- st_read(paste0(wd, "/", folders[4], "/", folders[4], "_Catchment.shp"), stringsAsFactors = FALSE)
# names(catch)[5] <- "WTSHNAME"
# catch$WTSHNAME <- rep("WHITE OAK BAYOU", length(catch$WTSHNAME))
# catch <- select(catch, WTSHNAME, geometry)
# catchment <- rbind(catchment, catch)
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

# crypress creek
cyp <- st_read(paste0(wd, "/", folders[10], "/", folders[10], "_Watershed.shp"), stringsAsFactors = FALSE)
cyp <- select(cyp, WTSHNAME, geometry)
cyp <- st_union(cyp)
cyp <- st_as_sf(data.frame("WTSHNAME" = "CYPRESS CREEK", cyp))

# buffalo bayou
baf <- st_read(paste0(wd, "/", folders[20], "/", folders[20], "_Watershed.shp"), stringsAsFactors = FALSE)
baf <- select(baf, WTSHNAME, geometry)

# barker reservoir
bar <- st_read(paste0(wd, "/", folders[18], "/", folders[18], "_Watershed.shp"), stringsAsFactors = FALSE)
bar <- select(bar, WTSHNAME, geometry)

# addicks reservoir (filename does not match)
# catch <- st_read(paste0(wd, "/", folders[19], "/Addicks_Modified_CAP.shp"), stringsAsFactors = FALSE)
# names(catch)[8] <- "WTSHNAME"
# catch$WTSHNAME <- rep("ADDICKS RESERVOIR", length(catch$WTSHNAME))
# catch <- select(catch, WTSHNAME, geometry)
# catchment <- rbind(catchment, catch)
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
cyp$REGION <- code
watershed <- rbind(watershed, cyp)

if (add$WTSHNAME[1] %in% one){
  code <- 1
} else if (add$WTSHNAME[1] %in% two){
  code <- 2
} else if (add$WTSHNAME[1] %in% three){
  code <- 3
}
add$REGION <- code
watershed <- rbind(watershed, add)


# plot(catchment)
plot(watershed)

st_write(watershed, "~/Spatial_Extreme_Value_Modeling/Test/watershed.shp", delete_layer = T)
st_write(select(watershed, WTSHNAME, geometry), "~/Spatial_Extreme_Value_Modeling/Test/watershed_name.shp", delete_layer = T)

regionall <- NULL
for (i in 1:3){
  region <- watershed %>% filter(REGION==i) %>% select(REGION, geometry)
  region <- st_as_sf(data.frame("REGION" = i, st_union(region)))
  regionall <- rbind(regionall, region)
}

regionall[[2]][[1]][[2]] <- NULL
regionall[[2]][[2]][[2]][[2]] <- NULL
for (i in length(regionall[[2]][[3]][[2]]):2){
  regionall[[2]][[3]][[2]][[i]] <- NULL
}

st_write(regionall, "~/Spatial_Extreme_Value_Modeling/Test/watershed_region.shp", delete_layer = T)

st_write(st_union(regionall), "~/Spatial_Extreme_Value_Modeling/Test/watershed_all.shp", delete_layer = T)

