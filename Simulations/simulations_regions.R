# Continuation of simulations.R

#### Functions ####
source('./Test/functionV2.R') #updated directHaus function + necessary libraries

# the following 2 functions are both from simulations.R
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


#### Watershed Regions ####
#old shapefiles
region_orig <- readOGR('./Test/watershed_region.shp')
region_orig <- spTransform(region_orig, CRS("+init=epsg:2278"))
plot(region_orig)
region_temp <- readOGR('./Test/watershed_region_temp.shp')
region_temp <- spTransform(region_temp, CRS("+init=epsg:2278"))
plot(region_temp)

#USE THIS UPDATED ONE!
region <- readOGR('./Test/watershed_region_updated/watershed_region_updated.shp')
proj4string(region) <- CRS("+init=epsg:2278")
region <- spTransform(region, CRS("+init=epsg:2278"))
plot(region)

#examine each region
plot(subset(region, region$REGION==1))
plot(subset(region, region$REGION==2))
plot(subset(region, region$REGION==3))


#### Starting Point: Region 1 ####
#end: region 2
region_results12a <- simulate(A=subset(region, region$REGION==1),
                              B=subset(region, region$REGION==2), f1=0.5)
results_summ(region_results12a)
cat(region_results12a, sep="\n")

region_results12b <- simulate(A=subset(region, region$REGION==1),
                              B=subset(region, region$REGION==2), f1=1.0)
results_summ(region_results12b)
cat(region_results12b, sep="\n")

#end: region 3
region_results13a <- simulate(A=subset(region, region$REGION==1),
                              B=subset(region, region$REGION==3), f1=0.5)
results_summ(region_results13a)
cat(region_results13a, sep="\n")

region_results13b <- simulate(A=subset(region, region$REGION==1),
                              B=subset(region, region$REGION==3), f1=1.0)
results_summ(region_results13b)
cat(region_results13b, sep="\n")


#### Starting Point: Region 2 ####
#end: region 1
region_results21a <- simulate(A=subset(region, region$REGION==2),
                              B=subset(region, region$REGION==1), f1=0.5)
results_summ(region_results21a)
cat(region_results21a, sep="\n")

#end: region 3
region_results23a <- simulate(A=subset(region, region$REGION==2),
                              B=subset(region, region$REGION==3), f1=0.5)
results_summ(region_results23a)
cat(region_results23a, sep="\n")


#### Starting Point: Region 3 ####
#end: region 1
region_results31a <- simulate(A=subset(region, region$REGION==3),
                              B=subset(region, region$REGION==1), f1=0.5)
results_summ(region_results31a)
cat(region_results31a, sep="\n")

#end: region 2
region_results32a <- simulate(A=subset(region, region$REGION==3),
                              B=subset(region, region$REGION==2), f1=0.5)
results_summ(region_results32a)
cat(region_results32a, sep="\n")

