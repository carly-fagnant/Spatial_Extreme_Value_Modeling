#Simulating directHaus function calculations to identify if there are 
#any significant differences in results due to spsample() randomness

#Read in the function(s)
source('./Test/function.R')

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
    results[i] <- directHaus(A, B, f1, f2, tol)
  }
  return(results)
}


#### Test Cases ####
#Read in Watershed data
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

#9: same region that is split in two for A and B - DOES NOT WORK
results9 <- unlist(simulate(A=subset(ws, WTSHNAME=="SPRING GULLY & GOOSE CREEK"),
                            B=subset(ws, WTSHNAME=="SPRING GULLY & GOOSE CREEK"), f1=.90))

#10: 
