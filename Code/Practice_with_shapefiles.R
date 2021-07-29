#######################
# Practice with shapefiles
#
#
# NOTES on Findings:
# - A (Clear Creek) and B (Armand Bayou) are combined
#     can be separated by subsetting [1, ] and [2, ]
#     catchments can be separated by subsetting [1:83, ] and [84:131, ]

# - K (Cypress Creek) and L (Little Cypress Creek) are combined
#     can be separated by subsetting [1, ] and [2, ]
#     catchments can be separated by subsetting [1:68, ] and [69:80, ]

# - F (San Jacinto & Galv Bay) is missing some catchments 
#      only has 9 and does not make up the same area as the watershed

# - G (San Jacinto River) is missing some catchments 
#      only has 69 and does not make up the same area as the watershed

# - S (Luce Bayou) is missing some catchments 
#      only has 3 and does not make up the same area as the watershed

# - U (Addicks Reservoir) is missing watershed shapefile
#      but does have subwatersheds that make up the total watershed area
#      Can we combine/union to get overall?
#######################




# set working directory
setwd("~/Documents/HCFCD Watershed Data/00_Example_Watershed/HEC-HMS/D_D100-00-00/Support_Data/Spatial_Data")
wd <- getwd()
# wd <- "Users/cf21/Documents/HCFCD Watershed Data/00_Example_Watershed/HEC-HMS/D_D100-00-00/Support_Data/Spatial_Data"

install.packages("sf")
library(sf)
shape <- sf::read_sf(dsn = wd, layer = "D_Brays_Bayou_Watershed")

install.packages("rgdal")
library(rgdal)
shape <- rgdal::readOGR(dsn = wd, layer = "D_Brays_Bayou_Watershed")
class(shape)
shape2 <- rgdal::readOGR(dsn = wd, layer = "D_Brays_Bayou_Catchment")
class(shape2)

plot(shape) # watershed
plot(shape2) # catchments (subwatersheds?) of watershed

# Subwatersheds?... Actually all watersheds together!!
setwd("~/Documents/HCFCD Watershed Data/HC_info")
wd2 <- getwd()
sub_ws <- rgdal::readOGR(dsn = wd2, layer = "hc_subwshds")
class(sub_ws)
plot(sub_ws) 
plot(sub_ws[2, ]) # 22 watersheds subset-able
# 1J, 2G, ..., 7P, 8R?

par(mfrow = c(1,1))
plot(sub_ws)
plot(shape)
plot(shape2)

# HEC-RAS
setwd("~/Documents/HCFCD Watershed Data/00_Example_Watershed/HEC-RAS/D133-00-00/Support_Data/Spatial_Data")
wd3 <- getwd()
ras <- rgdal::readOGR(dsn = wd3, layer = "S_XS")
class(ras)
plot(ras)



# # Read SHAPEFILE.shp from the current working directory (".")
# 
# require(rgdal)
# shape <- readOGR(dsn = ".", layer = "SHAPEFILE")
# 
# require(sf)
# shape <- read_sf(dsn = ".", layer = "SHAPEFILE")


# Removing Duplicate Files ------------------------------------------------

setwd("~/Documents/HCFCD Watershed Data/All_Watershed_Shapefiles")
wd <- getwd()

## A and B ####
# Watersheds
testA <- rgdal::readOGR(dsn = paste0(wd, "/A_Clear_Creek_B"), layer = "A_Clear_Creek_B_Armand_Bayou_Watershed")
testB <- rgdal::readOGR(dsn = paste0(wd, "/B_Armand_Bayou_A"), layer = "A_Clear_Creek_B_Armand_Bayou_Watershed")
identical(testA, testB)
plot(testA)
plot(testB)
plot(testB[1,]) # A: Clear Creek 
plot(testB[2,]) # B: Armand Bayou

# Catchments
testA <- rgdal::readOGR(dsn = paste0(wd, "/A_Clear_Creek_B"), layer = "A_Clear_Creek_B_Armand_Bayou_Catchment")
# testB <- rgdal::readOGR(dsn = paste0(wd, "/B_Armand_Bayou_A"), layer = "A_Clear_Creek_B_Armand_Bayou_Catchment")
identical(testA, testB)
plot(testA) # 131 features
plot(testA[1:83,]) # A: Clear Creek Catchments
plot(testA[84:131,]) # B: Armand Bayou Catchments

## K and L ####
# Watersheds
testK <- rgdal::readOGR(dsn = paste0(wd, "/K_Cypress_Creek_L"),        layer = "K_Cypress_Creek_L_Little_Cypress_Creek_Watershed")
testL <- rgdal::readOGR(dsn = paste0(wd, "/L_Little_Cypress_Creek_K"), layer = "K_Cypress_Creek_L_Little_Cypress_Creek_Watershed")
identical(testK, testL)
plot(testK)
plot(testK[1,]) # K: Cypress Creek
plot(testK[2,]) # L: Little Cypress Creek

# Catchments
testKc <- rgdal::readOGR(dsn = paste0(wd, "/K_Cypress_Creek_L"),        layer = "K_Cypress_Creek_L_Little_Cypress_Creek_Catchment")
testLc <- rgdal::readOGR(dsn = paste0(wd, "/L_Little_Cypress_Creek_K"), layer = "K_Cypress_Creek_L_Little_Cypress_Creek_Catchment")
identical(testKc, testLc)
plot(testKc) # 80 features
plot(testKc[1:68,]) # K: Cypress Creek Catchments
plot(testKc[69:80,]) # L: Little Cypress Creek Catchments


# Exploring Somewhat Different Files --------------------------------------

setwd("~/Documents/HCFCD Watershed Data/All_Watershed_Shapefiles")
wd <- getwd()

## F: San Jac & Galv Bay ####
# 2 versions identical, but only has 9 catchments and does not make up the same area as the watershed
testF <- rgdal::readOGR(dsn = paste0(wd, "/F_San_Jacinto_Galveston_Bay"), layer = "F_San_Jacinto_Galveston_Bay_Watershed")
testFc <- rgdal::readOGR(dsn = paste0(wd, "/F_San_Jacinto_Galveston_Bay"), layer = "F_San_Jacinto_Galveston_Bay_Catchment")
plot(testF)
plot(testFc)
testF2 <- rgdal::readOGR(dsn = "/Users/cf21/Documents/HCFCD Watershed Data/All_files_messy/F_GalvBay_FEMA_Effective/HEC-HMS/F216_00_00/Support_Data/Spatial_Data", layer = "F_San_Jacinto_Galveston_Bay_Watershed")
testF2c <- rgdal::readOGR(dsn = "/Users/cf21/Documents/HCFCD Watershed Data/All_files_messy/F_GalvBay_FEMA_Effective/HEC-HMS/F216_00_00/Support_Data/Spatial_Data", layer = "F_San_Jacinto_Galveston_Bay_Catchment")
identical(testF, testF2)

## G: San Jacinto River ####
# 19 versions in All_files_messy... have not tested them all
# only has 69 catchments and does not make up the same area as the watershed
testG <- rgdal::readOGR(dsn = paste0(wd, "/G_San_Jacinto_River"), layer = "G_San_Jacinto_River_Watershed")
testGc <- rgdal::readOGR(dsn = paste0(wd, "/G_San_Jacinto_River"), layer = "G_San_Jacinto_River_Catchment")
plot(testG)
plot(testGc)
plot(sub_ws[2, ]) # compare to old watershed boundaries we had

testG112 <- rgdal::readOGR(dsn = "/Users/cf21/Documents/HCFCD Watershed Data/All_files_messy/G_SanJac_FEMA_Effective/HEC-HMS/G112_00_00/Support_Data/Spatial_Data", layer = "G_San_Jacinto_River_Watershed")
testG112c <- rgdal::readOGR(dsn = "/Users/cf21/Documents/HCFCD Watershed Data/All_files_messy/G_SanJac_FEMA_Effective/HEC-HMS/G112_00_00/Support_Data/Spatial_Data", layer = "G_San_Jacinto_River_Catchment")
plot(testG112)
plot(testG112c)


## Q: Cedar Bayou ####
# extra files in LOMR for drainage areas - not needed
testQ <- rgdal::readOGR(dsn = paste0(wd, "/Q_Cedar_Bayou"), layer = "Q_Cedar_Bayou_Watershed")
testQc <- rgdal::readOGR(dsn = paste0(wd, "/Q_Cedar_Bayou"), layer = "Q_Cedar_Bayou_Catchment")
plot(testQ)
plot(testQc)
testQcor <- rgdal::readOGR(dsn = paste0(wd, "/Q_Cedar_Bayou/LOMR 16-06-0437P"), layer = "Corrected_Effective_DA") # same as McGee
testQmcgee <- rgdal::readOGR(dsn = paste0(wd, "/Q_Cedar_Bayou/LOMR 16-06-0437P"), layer = "DA_McGEE_GULLY")
testQprop <- rgdal::readOGR(dsn = paste0(wd, "/Q_Cedar_Bayou/LOMR 16-06-0437P"), layer = "Proposed_Drainage_Areas")
testQflow <- rgdal::readOGR(dsn = paste0(wd, "/Q_Cedar_Bayou/LOMR 16-06-0437P"), layer = "Flowpaths")
plot(testQcor)
plot(testQmcgee)
plot(testQprop)
plot(testQflow)


## S: Luce Bayou ####
# Only 3 subwatersheds, don't make up whole shape of watershed
# 2 versions identical
wd4 <- "/Users/cf21/Documents/HCFCD Watershed Data/All_files_messy/S_Luce_FEMA_Effective/HEC-HMS"
testS <- rgdal::readOGR(dsn = paste0(wd4, "/S110_00_00/Support_Data/Spatial_Data"), layer = "S_Luce_Bayou_Watershed")
testSc <- rgdal::readOGR(dsn = paste0(wd4, "/S110_00_00/Support_Data/Spatial_Data"), layer = "S_Luce_Bayou_Catchment")
testS2 <- rgdal::readOGR(dsn = paste0(wd4, "/S114_00_00/Support_Data/Spatial_Data"), layer = "S_Luce_Bayou_Watershed")
testS2c <- rgdal::readOGR(dsn = paste0(wd4, "/S114_00_00/Support_Data/Spatial_Data"), layer = "S_Luce_Bayou_Catchment")
plot(testS)
plot(testS2)
plot(testSc)
plot(testS2c)
identical(testS, testS2)


## U: Addicks Reservoir ####
# MISSING SHAPEFILE for overall watershed polygon
# has subwatersheds (subbasins) and stream? lines (CAP)
testUlines <- rgdal::readOGR(dsn = paste0(wd, "/U_Addicks_Reservoir"), layer = "Addicks_Modified_CAP")
testUc <- rgdal::readOGR(dsn = paste0(wd, "/U_Addicks_Reservoir"), layer = "S_Subbasins")
plot(testUlines)
plot(testUc)

# 2 versions are identical
testU1L <- rgdal::readOGR(dsn = "/Users/cf21/Documents/HCFCD Watershed Data/All_files_messy/U_Addicks_FEMA_Effective/HEC-HMS/U_U100-00-00/maps", layer = "Addicks_Modified_CAP")
testU1c <- rgdal::readOGR(dsn = "/Users/cf21/Documents/HCFCD Watershed Data/All_files_messy/U_Addicks_FEMA_Effective/HEC-HMS/U_U100-00-00/maps", layer = "S_Subbasins")
testU2L <- rgdal::readOGR(dsn = "/Users/cf21/Documents/HCFCD Watershed Data/All_files_messy/U_Addicks_FEMA_Effective/HEC-HMS/U_U100-00-00/Support_Data/Spatial_Data", layer = "Addicks_Modified_CAP")
testU2c <- rgdal::readOGR(dsn = "/Users/cf21/Documents/HCFCD Watershed Data/All_files_messy/U_Addicks_FEMA_Effective/HEC-HMS/U_U100-00-00/Support_Data/Spatial_Data", layer = "S_Subbasins")
identical(testU1L, testU2L)
identical(testU1c, testU2c)
plot(testU1c)
plot(testU2c)




# Hydroregions ------------------------------------------------------------

library(rgdal)
setwd("~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data")
wd <- getwd()

G_north <- rgdal::readOGR(dsn = paste0(wd,"/G_San_Jacinto_River_split"), layer = "G_San_Jacinto_River_Watershed_split_north")
G_south <- rgdal::readOGR(dsn = paste0(wd,"/G_San_Jacinto_River_split"), layer = "G_San_Jacinto_River_Watershed_split_south")

plot(G_north)
plot(G_south, axes=T)
plot(G_north, add=T)
plot(testG, axes=T)
plot(G_south, add=T)
