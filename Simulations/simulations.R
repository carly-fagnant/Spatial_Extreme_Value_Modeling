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

results1 <- simulate(1000, subset(ws, WTSHNAME=="CLEAR CREEK"), subset(ws, WTSHNAME=="ARMAND BAYOU"), .90)
results1 <- unlist(results1)
summary(results1)
sd(results1)
#cat(results1, sep="\n") #use to get results to paste into Excel column

results2 <- unlist(simulate(1000, A=subset(ws, WTSHNAME=="ARMAND BAYOU"), B=subset(ws, WTSHNAME=="CLEAR CREEK"), f1=.90))
summary(results2)
sd(results2)

