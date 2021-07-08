#Simulating directHaus function calculations to identify if there are 
#any significant differences in results due to spsample() randomness

#Read in the function(s)
source('./Test/function.R')
source('./Test/functionV2.R') #read in updated functions
#All simulations will be done under the old function unless otherwise indicated (f1=0, f1=1 cases)

#The following "simulate" function will be used to store simulation results.
#The parameters aside from num are those of the directHaus function:
#   A,B: region calculate the extended Hausdorff distance- SpatialPolygons or SpatialPolygonsDataFrame
#   f1: the percentage of area in B you want captured as a decimal (eg 10\% = .1)
#   f2: the percentage of area in A you want captured as a decimal (eg 10\% = .1); defaults to the value of f1
#   tol: tolerance for selecting the epsilon buffer to yield desired f1. Default is 1/10000th the sampled directional distances.
#The directHaus function returns the directional extended hausdorff distance from A to B
simulate <- function(num=1000, A, B, f1, f2=f1, tol=NULL) {
  results <- rep(NA, num)
  for (i in 1:num) {
    #results[i] <- directHaus(A, B, f1, f2, tol)
    results[i] <- directHaus(A, B, f1, tol)
  }
  return(results)
}

#Used to easily access summary statistics to paste into Excel
results_summ <- function(results) {
  summ <- rep(NA, 8)
  summ[1] <- length(results)
  summ[2] <- sd(results)
  summ[3:8] <- summary(results)
  cat(summ, sep="\n")
}


# Test Cases --------------------------------------------------------------
#### Read in Watershed data ####
shape <- readOGR('./Test/watershed.shp')
ws <- spTransform(shape, CRS("+init=epsg:2278"))
plot(ws)
ws$WTSHNAME
#plot(ws[ws$WTSHNAME=="CLEAR CREEK",])
#plot(subset(ws, WTSHNAME=="CLEAR CREEK"))

#1: regions next to each other, border shared
plot(subset(ws, WTSHNAME=="CLEAR CREEK" | WTSHNAME=="ARMAND BAYOU"))
results1 <- simulate(1000, subset(ws, WTSHNAME=="CLEAR CREEK"), subset(ws, WTSHNAME=="ARMAND BAYOU"), .90)
results1 <- unlist(results1)
summary(results1)
sd(results1)
cat(results1, sep="\n") #use to get results to paste into Excel column

#2: same as 1 but in the other direction
results2 <- unlist(simulate(1000, A=subset(ws, WTSHNAME=="ARMAND BAYOU"), B=subset(ws, WTSHNAME=="CLEAR CREEK"), f1=.90))
summary(results2)
sd(results2)
cat(results2, sep="\n") #use to get results to paste into Excel column

#3: same as 1 but with different area percentage
results3 <- unlist(simulate(A=subset(ws, WTSHNAME=="CLEAR CREEK"),
                            B=subset(ws, WTSHNAME=="ARMAND BAYOU"), f1=.50))
summary(results3)
sd(results3)
cat(results3, sep="\n") #use to get results to paste into Excel column

#4: some distance between regions, regions have a elongated shape
#this one takes a while longer to run
plot(subset(ws, WTSHNAME=="SPRING CREEK" | WTSHNAME=="BUFFALO BAYOU"))
#results4 <- unlist(simulate(100, A=subset(ws, WTSHNAME=="SPRING CREEK"),
#                            B=subset(ws, WTSHNAME=="BUFFALO BAYOU"), f1=.80))
results4 <- unlist(simulate(A=subset(ws, WTSHNAME=="SPRING CREEK"),
                            B=subset(ws, WTSHNAME=="BUFFALO BAYOU"), f1=.80))
summary(results4)
sd(results4)
cat(results4, sep="\n") #use to get results to paste into Excel column

#5: region shared, also takes a while longer
plot(subset(ws, WTSHNAME=="SAN JACINTO RIVER" | WTSHNAME=="SPRING GULLY & GOOSE CREEK"))
results5 <- unlist(simulate(1000, A=subset(ws, WTSHNAME=="SAN JACINTO RIVER"),
                            B=subset(ws, WTSHNAME=="SPRING GULLY & GOOSE CREEK"), f1=.80))
summary(results5)
sd(results5)
cat(results5, sep="\n") #use to get results to paste into Excel column

#6: detached/island type (San Jacinto & Galveston Bay)
results6 <- unlist(simulate(1000, A=subset(ws, WTSHNAME=="SAN JACINTO & GALVESTON BAY"),
                           B=subset(ws, WTSHNAME=="VINCE BAYOU"), f1=.80))
summary(results6)
sd(results6)
cat(results6, sep="\n") #use to get results to paste into Excel column

#7: two rather irregularly shaped regions
plot(subset(ws, WTSHNAME=="CYPRESS CREEK" | WTSHNAME=="BUFFALO BAYOU"))
results7 <- unlist(simulate(A=subset(ws, WTSHNAME=="CYPRESS CREEK"),
                            B=subset(ws, WTSHNAME=="BUFFALO BAYOU"), f1=.90))
summary(results7)
sd(results7)
cat(results7, sep="\n") #use to get results to paste into Excel column

#8: two neighboring more regularly shaped regions
plot(subset(ws, WTSHNAME=="ARMAND BAYOU" | WTSHNAME=="VINCE BAYOU"))
results8 <- unlist(simulate(A=subset(ws, WTSHNAME=="ARMAND BAYOU"),
                            B=subset(ws, WTSHNAME=="VINCE BAYOU"), f1=.90))
summary(results8)
sd(results8)
cat(results8, sep="\n") #use to get results to paste into Excel column

#9: same region that is split in two for A and B - USE FUNCTION V2
results9 <- unlist(simulate(A=subset(ws, WTSHNAME=="SPRING GULLY & GOOSE CREEK"),
                            B=subset(ws, WTSHNAME=="SPRING GULLY & GOOSE CREEK"), f1=.90))
results9 <- unlist(simulate(A=subset(ws, WTSHNAME=="LUCE BAYOU"),
                            B=subset(ws, WTSHNAME=="LUCE BAYOU"), f1=.60))
results_summ(results9)
cat(results9, sep="\n")


#10-20: different f1/f2 combinations for the same region
plot(subset(ws, WTSHNAME=="GREENS BAYOU" | WTSHNAME=="JACKSON BAYOU"))
results10 <- unlist(simulate(A=subset(ws, WTSHNAME=="GREENS BAYOU"),
                             B=subset(ws, WTSHNAME=="JACKSON BAYOU"), f1=.25))
results_summ(results10)
cat(results10, sep="\n")


results11 <- unlist(simulate(A=subset(ws, WTSHNAME=="GREENS BAYOU"),
                             B=subset(ws, WTSHNAME=="JACKSON BAYOU"), f1=.50))
results_summ(results11)
cat(results11, sep="\n")


results12 <- unlist(simulate(A=subset(ws, WTSHNAME=="GREENS BAYOU"),
                             B=subset(ws, WTSHNAME=="JACKSON BAYOU"), f1=.75))
results_summ(results12)
cat(results12, sep="\n")

#run simulation 13 with updated function
results13 <- unlist(simulate(A=subset(ws, WTSHNAME=="GREENS BAYOU"),
                             B=subset(ws, WTSHNAME=="JACKSON BAYOU"), f1=1))
results_summ(results13)
cat(results13, sep="\n")


#DOES NOT WORK WITH f=0 in old function: Error in slot(overlap.region@polygons[[1]], "area") : 
#  trying to get slot "polygons" from an object of a basic class ("NULL") with no slots
#results14 <- unlist(simulate(A=subset(ws, WTSHNAME=="GREENS BAYOU"),
#                             B=subset(ws, WTSHNAME=="JACKSON BAYOU"), f1=0))
results14 <- unlist(simulate(A=subset(ws, WTSHNAME=="GREENS BAYOU"),
                             B=subset(ws, WTSHNAME=="JACKSON BAYOU"), f1=.99))
results_summ(results14)
cat(results14, sep="\n")


results15 <- unlist(simulate(A=subset(ws, WTSHNAME=="GREENS BAYOU"),
                             B=subset(ws, WTSHNAME=="JACKSON BAYOU"), f1=.01))
results_summ(results15)
cat(results15, sep="\n")


results16 <- unlist(simulate(A=subset(ws, WTSHNAME=="GREENS BAYOU"),
                             B=subset(ws, WTSHNAME=="JACKSON BAYOU"), f1=.3, f2=.7))
results_summ(results16)
cat(results16, sep="\n")


results17 <- unlist(simulate(A=subset(ws, WTSHNAME=="GREENS BAYOU"),
                             B=subset(ws, WTSHNAME=="JACKSON BAYOU"), f1=.4, f2=.6))
results_summ(results17)
cat(results17, sep="\n")


results18 <- unlist(simulate(A=subset(ws, WTSHNAME=="GREENS BAYOU"),
                             B=subset(ws, WTSHNAME=="JACKSON BAYOU"), f1=.7, f2=.3))
results_summ(results18)
cat(results18, sep="\n")


results19 <- unlist(simulate(A=subset(ws, WTSHNAME=="GREENS BAYOU"),
                             B=subset(ws, WTSHNAME=="JACKSON BAYOU"), f1=.6, f2=.4))
results_summ(results19)
cat(results19, sep="\n")


results20 <- unlist(simulate(A=subset(ws, WTSHNAME=="GREENS BAYOU"),
                             B=subset(ws, WTSHNAME=="JACKSON BAYOU"), f1=.9, f2=.2))
results_summ(results20)
cat(results20, sep="\n")

#21-24: changing tolerance
plot(subset(ws, WTSHNAME=="ADDICKS RESERVOIR" | WTSHNAME=="BRAYS BAYOU"))
results21 <- unlist(simulate(A=subset(ws, WTSHNAME=="ADDICKS RESERVOIR"),
                             B=subset(ws, WTSHNAME=="BRAYS BAYOU"), f1=.80, tol=1/1000))
results_summ(results21)
cat(results21, sep="\n")


results22 <- unlist(simulate(A=subset(ws, WTSHNAME=="ADDICKS RESERVOIR"),
                             B=subset(ws, WTSHNAME=="BRAYS BAYOU"), f1=.80, tol=1/100))
results_summ(results22)
cat(results22, sep="\n")

results23 <- unlist(simulate(A=subset(ws, WTSHNAME=="ADDICKS RESERVOIR"),
                             B=subset(ws, WTSHNAME=="BRAYS BAYOU"), f1=.80, tol=1/10))
results_summ(results23)
cat(results23, sep="\n")

results24 <- unlist(simulate(A=subset(ws, WTSHNAME=="ADDICKS RESERVOIR"),
                             B=subset(ws, WTSHNAME=="BRAYS BAYOU"), f1=.50, f2=.75, tol=1/100))
results_summ(results24)
cat(results24, sep="\n")

#25-30: maybe the furthest regions
plot(subset(ws, WTSHNAME=="SPRING CREEK" | WTSHNAME=="CLEAR CREEK"))
results25 <- unlist(simulate(A=subset(ws, WTSHNAME=="SPRING CREEK"),
                             B=subset(ws, WTSHNAME=="CLEAR CREEK"), f1=.50))
results_summ(results25)
cat(results25, sep="\n")

#use updated function since f1=1
results26 <- unlist(simulate(A=subset(ws, WTSHNAME=="SPRING CREEK"),
                             B=subset(ws, WTSHNAME=="CLEAR CREEK"), f1=1))
results_summ(results26)
cat(results26, sep="\n")


plot(subset(ws, WTSHNAME=="BARKER RESERVOIR" | WTSHNAME=="LUCE BAYOU"))
results27 <- unlist(simulate(A=subset(ws, WTSHNAME=="BARKER RESERVOIR"),
                             B=subset(ws, WTSHNAME=="LUCE BAYOU"), f1=.50))
results_summ(results27)
cat(results27, sep="\n")


results28 <- unlist(simulate(A=subset(ws, WTSHNAME=="BARKER RESERVOIR"),
                             B=subset(ws, WTSHNAME=="LUCE BAYOU"), f1=.30))
results_summ(results28)
cat(results28, sep="\n")


plot(subset(ws, WTSHNAME=="BARKER RESERVOIR" | WTSHNAME=="CEDAR BAYOU"))
results29 <- unlist(simulate(A=subset(ws, WTSHNAME=="BARKER RESERVOIR"),
                             B=subset(ws, WTSHNAME=="CEDAR BAYOU"), f1=.50))
results_summ(results29)
cat(results29, sep="\n")


#use updated function since f1=1
plot(subset(ws, WTSHNAME=="BARKER RESERVOIR" | WTSHNAME=="CEDAR BAYOU"))
#Barker is on the left, Cedar on the right
plot(subset(ws, WTSHNAME=="BARKER RESERVOIR"))
plot(subset(ws, WTSHNAME=="CEDAR BAYOU"))
results30 <- unlist(simulate(A=subset(ws, WTSHNAME=="BARKER RESERVOIR"),
                             B=subset(ws, WTSHNAME=="CEDAR BAYOU"), f1=1))
results_summ(results30)
cat(results30, sep="\n")


#31-35: A and/or B are a few watershed regions instead of just 1
plot(subset(ws, WTSHNAME=="SPRING CREEK" | WTSHNAME=="GREENS BAYOU" | WTSHNAME=="BARKER RESERVOIR"))
results31 <- unlist(simulate(A=subset(ws, WTSHNAME=="SPRING CREEK" | WTSHNAME=="GREENS BAYOU"),
                             B=subset(ws, WTSHNAME=="BARKER RESERVOIR"), f1=.5))
results_summ(results31)
cat(results31, sep="\n")

plot(subset(ws, WTSHNAME=="WHITE OAK BAYOU" | WTSHNAME=="CARPENTERS BAYOU" | WTSHNAME=="GREENS BAYOU"))
results32 <- unlist(simulate(A=subset(ws, WTSHNAME=="WHITE OAK BAYOU" | WTSHNAME=="CARPENTERS BAYOU"),
                             B=subset(ws, WTSHNAME=="GREENS BAYOU"), f1=.5))
results_summ(results32)
cat(results32, sep="\n")

plot(subset(ws, WTSHNAME=="WILLOW CREEK" | WTSHNAME=="JACKSON BAYOU" |
              WTSHNAME=="SIMS BAYOU" | WTSHNAME=="ARMAND BAYOU"))
#two regions of both watersheds are not touching
results33 <- unlist(simulate(A=subset(ws, WTSHNAME=="WILLOW CREEK" | WTSHNAME=="JACKSON BAYOU"),
                             B=subset(ws, WTSHNAME=="SIMS BAYOU" | WTSHNAME=="ARMAND BAYOU"), f1=.5))
results_summ(results33)
cat(results33, sep="\n")

plot(subset(ws, WTSHNAME=="SPRING CREEK" | WTSHNAME=="SAN JACINTO RIVER" |
              WTSHNAME=="SIMS BAYOU" | WTSHNAME=="CLEAR CREEK"))
results34 <- unlist(simulate(A=subset(ws, WTSHNAME=="SPRING CREEK" | WTSHNAME=="SAN JACINTO RIVER"),
                             B=subset(ws, WTSHNAME=="SIMS BAYOU" | WTSHNAME=="CLEAR CREEK"), f1=.85))
results_summ(results34)
cat(results34, sep="\n")

plot(subset(ws, WTSHNAME=="CLEAR CREEK" |
              WTSHNAME=="CYPRESS CREEK" | WTSHNAME=="WILLOW CREEK" | WTSHNAME=="SPRING CREEK"))
results35 <- unlist(simulate(A=subset(ws, WTSHNAME=="CLEAR CREEK"),
                             B=subset(ws, WTSHNAME=="CYPRESS CREEK" | WTSHNAME=="WILLOW CREEK" | WTSHNAME=="SPRING CREEK"), f1=.5))
results_summ(results35)
cat(results35, sep="\n")


## Run this & record results for a number of times, change f parameter as needed
#use of updated function
areas <- sample(1:length(ws$WTSHNAME), 2, replace=FALSE)
(a <- ws$WTSHNAME[areas[1]])
(b <- ws$WTSHNAME[areas[2]])
ws_results1 <- simulate(A=subset(ws, WTSHNAME==a),
                        B=subset(ws, WTSHNAME==b), f1=0.5)

results_summ(ws_results1)
cat(ws_results1, sep="\n")

ws_results2 <- simulate(A=subset(ws, WTSHNAME==a),
                        B=subset(ws, WTSHNAME==b), f1=1.0)
results_summ(ws_results2)
cat(ws_results2, sep="\n")

ws_results3 <- simulate(A=subset(ws, WTSHNAME==a),
                        B=subset(ws, WTSHNAME==b), f1=0.9)
results_summ(ws_results3)
cat(ws_results3, sep="\n")



plot(subset(ws, WTSHNAME=="SAN JACINTO & GALVESTON BAY" |
              WTSHNAME=="CYPRESS CREEK" | WTSHNAME=="BARKER RESERVOIR"))
results54 <- simulate(A=subset(ws, WTSHNAME=="SAN JACINTO & GALVESTON BAY"),
                      B=subset(ws, WTSHNAME=="CYPRESS CREEK" | WTSHNAME=="BARKER RESERVOIR"), f1=.5)
results_summ(results54)
cat(results54, sep="\n")

plot(subset(ws, WTSHNAME=="CYPRESS CREEK" |
              WTSHNAME=="VINCE BAYOU" | WTSHNAME=="ARMAND BAYOU" | WTSHNAME=="CLEAR CREEK"))
results55 <- simulate(A=subset(ws, WTSHNAME=="CYPRESS CREEK"),
                      B=subset(ws, WTSHNAME=="VINCE BAYOU" | WTSHNAME=="ARMAND BAYOU" | WTSHNAME=="CLEAR CREEK"), f1=.5)
results_summ(results55)
cat(results55, sep="\n")

#### Watershed Regions ####
region <- readOGR('./Test/watershed_region.shp')
region <- spTransform(region, CRS("+init=epsg:2278"))
#examine each region
plot(region)
plot(subset(region, region$REGION==1))
plot(subset(region, region$REGION==2))
plot(subset(region, region$REGION==3))

region_results1 <- unlist(simulate(A=subset(region, region$REGION==1),
                                   B=subset(region, region$REGION==2), f1=0.5))
results_summ(region_results1)
cat(region_results1, sep="\n")

region_results2 <- unlist(simulate(A=subset(region, region$REGION==1),
                                   B=subset(region, region$REGION==3), f1=0.5))

region_results3 <- unlist(simulate(A=subset(region, region$REGION==2),
                                   B=subset(region, region$REGION==3), f1=0.5))

region_results4 <- unlist(simulate(A=subset(region, region$REGION==3),
                                   B=subset(region, region$REGION==1), f1=0.5))

#### Zip Codes ####
shape <- readOGR('./Data/Zip_Codes/Zip_Codes.shp')
zip <- spTransform(shape, CRS("+init=epsg:2278"))
plot(zip)

#Run this & record results for a number of times, change f parameter as needed
zipcodes <- sample(1:length(zip$ZIP_CODE), 2, replace=FALSE)
(a <- zip$ZIP_CODE[zipcodes[1]])
(b <- zip$ZIP_CODE[zipcodes[2]])
zip_results <- unlist(simulate(A=subset(zip, ZIP_CODE==a),
                               B=subset(zip, ZIP_CODE==b), f1=0.5))

results_summ(zip_results)
cat(zip_results, sep="\n")
