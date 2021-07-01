# Libraries ---------------------------------------------------------------
library(doParallel)
library(ggsn)
library(raster)
library(rgdal)
library(rgeos)
library(spatialreg)
library(spdep)

# Notes on key changes ----------------------------------------------------

#hausMat
#To do: ensure that it can correctly handle the asymmetric case where f1 is not equal to f2

# extHaus
# Added checks to determine the type of A and B in order to decide how to perform the extended hausdorff distance calculations
# Extended the function's capabilities to handle {point, point}, {point, line}, {point, area}, {line, line}, {line, area}, and {area, area} calculations
# In the cases where f1 == 1 or f2 == 1, replaced the calls to gDistance(A, B, hausdorff=T) with calls to directHaus(A, B, 1) and directHaus(B, A, 1) (respectively)
  # gDistance(A, B, hausdorff=T) returns the value of H(A, B) whereas in the cases above we are actually interested in h(A, B) (when f1 == 1) and h(B, A) (when f2 == 1)
  # Since H(A, B) = max(h(A, B), h(B, A)), it can be incorrect to use gDistance(A, B, hausdorff=T) in either of these cases.

# directHaus
# Removed the mostly unused parameter f2
# Allowed directHaus to perform computations when f1=1 instead of throwing an error
# Extendeded the function's capabilities to handle line-to-area, line-to-line, line-to-point, area-to-area, area-to-line, and area-to-point calculations
# The output was originally a list, of one number. It now directly outputs a numeric value so that it doesn't have to be accessed by $direct.haus[1] in the other functions.
# Added check if the regions are the same. At the moment 0 is returned but should 0 be returned?
# Added a case for f1=0 (simply return the minimum distance between A and B instead of performing the epsilon-buffer procedure)

# pointHaus
# Replaced gDistance(point, B, hausdorff=T, byid=T) with gDistance(point, B)
  # gDistance(point, B, hausdorff=T, byid=T) returns H(point, B) and not
  # h(point, B) (as desired). Since point is a point, h(point, B) is equal to
  # the minimum distance between point and B which can be found using
  # gDistance(point, B)
# Changed variable names to make them consistent with the input parameter names

# hausMat -----------------------------------------------------------------
#' Creating a matrix of (extended) Hausdorff distances
#' 
#' This function takes a SpatialPolygons object and a decimal in [0, 1].
#'  
#'@param shp a SpatialPolygons object.
#'@param f1 The percentage (as a decimal) of region i to retain when calculating the directional Hausdorff distance from region i to region j.
#'@param f2 The percentage (as a decimal) of region j to retain when calculating the directional Hausdorff distance from region j to i. Defaults to the value of f1. Note that specifying a different value will result in a non-symmetric matrix.
#'@param fileout Should the resulting weight matrix be written to file? Defaults to FALSE 
#'@param filename If "fileout" is TRUE, the name for the file to be outputted.
#'@param ncores If "do.parallel" is true, the number of cores to be used for parallel computation. Default is 1. 
#'@param timer If T, records "Sys.time()" when the function is called and outputs the elapsed time along with the matrix. Default is F.
#'@param do.parallel If T, uses the number of cores specified using "ncores" to set up a cluster with the "foreach" package. 
#'  
#'@return an nxn matrix of requested distances.
hausMat <- function(shp, f1, f2=f1, fileout=FALSE, filename=NULL, ncores=1, timer=F, do.parallel=T) {
  if (timer) {
    start <- Sys.time()
  }
  n <- nrow(shp@data)
  combs <- combn(1:n, 2) # n choose 2 combinations
  n.combs <- ncol(combs);
  haus.dists <- matrix(0, nrow=n, ncol=n)
  if (do.parallel) {
    #out <- matrix(-1, nrow = 20, ncol = 20)
    start <- Sys.time()
    print("Setting up parallelization")
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    print("Computing...")
    # compute the extended hausdorff distance for every combination of regions (in parallel)
    out <- foreach (i = 1:n.combs, .packages=c("rgeos", "sp", "raster"), .combine=rbind, .export=c("directHaus", "extHaus")) %dopar% {
      extHaus(shp[combs[1,i],], shp[combs[2,i],], f1=f1)
    }
    print("Completion time:")
    print(Sys.time() - start)
    # save(out, file = "hausdorff_columbus")
    stopCluster(cl)
    # closeAllConnections()  
    haus.dists[lower.tri(haus.dists, diag=F)] <- out
    # use the fact that the Hausdorff distance matrix will be symmetrical to compute the upper triangular entries
    haus.dists <- haus.dists + t(haus.dists)
  } else {
    # compute the extended hausdorff distance for every combination of regions (sequentially)
    for (i in 1:n.combs) {
      haus.dists[combs[1,i], combs[2,i]] <- extHaus(shp[combs[1,i],], shp[combs[2,i],], f1=f1)
    }
    # use the fact that the Hausdorff distance matrix will be symmetrical to compute the upper triangular entries
    haus.dists <- haus.dists + t(haus.dists)
  }
  ## add 0's to diagonal
  #if(fileout){save(haus.dists, file=filename)}
  return(haus.dists)
}

# extHaus -----------------------------------------------------------------
#' Extended Hausdorff Distance
#'
#' Calculates the extended Hausdorff distance from A to B for given quantiles and precision
#' 
#'@author Julia Schedler
#'
#'@param A region/line/point used to calculate the extended Hausdorff distance: SpatialPolygons or SpatialLines or SpatialPoints
#'@param B region/line/point used to calculate the extended Hausdorff distance: SpatialPolygons or SpatialLines or SpatialPoints
#'@param f1 the percentage of area/length in B you want captured as a decimal (eg 10% = .1)
#'@param f2 the percentage of area/length in A you want captured as a decimal (eg 10% = .1). Defaults to the same value as f1.
#'@param tol value to be passed to "tol" in "directHaus".
#'@return the extended hausdorff distance (max of directional from A to B and B to A)
#'
#'Last Edited: July 01 2021
extHaus <- function(A, B, f1, f2=f1, tol=NULL) {
  if (!is.projected(A) | !is.projected(B)) {
    stop(paste("Spatial* object (inputs ", quote(A), ", ", quote(B), ") must be projected. Try running ?spTransform().", sep=""))
  }
  if (is.null(gDifference(A, B))) {
    return(0)
  } # check if two regions are the same, if so return a hausdorff distance of zero.
  
  # get the class of A and B to perform calculations accordingly
  type.of.A <- class(A)[1]
  type.of.B <- class(B)[1]
  
  if (type.of.A == "SpatialPoints" && type.of.B == "SpatialPoints") {
    # if both A and B are points, return the cartesian minimum distance between them
    return (gDistance(A, B))
  } else if (type.of.A == "SpatialPoints") {
    # if one of A or B are points, delegate to pointHaus
    return(pointHaus(A, B, f2, tol=tol))
  } else if (type.of.B == "SpatialPoints") {
    return(pointHaus(B, A, f1, tol=tol))
  } else if (type.of.A == "SpatialLines" || type.of.A == "SpatialPolygons" || type.of.A == "SpatialPolygonsDataFrame") {
    if (f1 == 0) {
      # if f1 = 0 use the cartesian minimum distance between A and B
      A_to_B <- gDistance(A, B)
    } else {
      # use the extended directed hausdorff distance from A to B
      A_to_B <- directHaus(A, B, f1, tol=tol)
    }
  } else {
    stop("A is not SpatialPoints, SpatialLines, SpatialPolygons, or SpatialPolygonsDataFrame")
  }
  
  if (type.of.B == "SpatialLines" || type.of.B == "SpatialPolygons" || type.of.B == "SpatialPolygonsDataFrame") {
    if (f2 == 0) {
      # if f2 = 0 use the cartesian minimum distance between A and B
      B_to_A <- gDistance(A, B)
    } else {
      # use the directed extended hausdorff distance from B to A
      B_to_A <- directHaus(B, A, f2, tol=tol)
    }
  } else {
    stop("B is not SpatialPoints, SpatialLines, SpatialPolygons, or SpatialPolygonsDataFrame")
  }
  haus.dist <- max(A_to_B, B_to_A)
  return(haus.dist)
}

# directHaus --------------------------------------------------------------
#' Calculate the directional extended Hausdorff Distance
#' 
#' This function takes a SpatialPolygons/SpatialLines object, a SpatialPolygons/SpatialLines/SpatialPoints object and a decimal in [0, 1]
#'  
#'@author Julia Schedler
#'
#'@param A region/line used to calculate the directed extended Hausdorff distance: SpatialPolygons or SpatialLines
#'@param B region/line/point used to calculate the directed extended Hausdorff distance: SpatialPolygons or SpatialLines or SpatialPoints
#'@param f1 the percentage of area in A you want captured as a decimal (eg 10\% = .1)
#'@param tol tolerance for selecting the epsilon buffer to yield desired f1. Default is 1/10000th the sampled directional distances.
#'
#'@return The directional extended hausdorff distance from A to B
#'
#'Last Edited: July 01 2021
directHaus <- function(A, B, f1, tol=NULL) {
  if (!is.projected(A) | !is.projected(B)) {
    stop(paste("Spatial* object (inputs ", quote(A),", ", quote(B), ") must be projected. Try running ?spTransform().", sep = ""))
  }
  if(is.null(gDifference(A,B))) {
    return(0) ## check if two regions are the same, if so return HD of zero
  }
  if (f1 == 0) {
    return(gDistance(A, B))
  }
  # generate points
  n <- 10000
  a.coords <- spsample(A, n=n, type="regular") # sample points within A
  # sp::spsample should work with both polygons and lines: https://cran.r-project.org/web/packages/sp/sp.pdf
  # compute minimum distance of points a.coords to a point in B
  dists <- gDistance(a.coords, B, byid=T)
  
  if(is.null(tol)){
    tol <- sd(dists[1,]) / 100
  }
  # find desired quantile of distances
  epsilon <- as.numeric(quantile(dists[1,], f1)) #buffer width
  buff <- buffer(B, width=epsilon, dissolve=T) 
  # raster::buffer should be point/polgon/line friendly
  overlap.region <- gIntersection(buff, A) # overlap region of buffer+B with A
  # rgeos:gIntersection should be polygon/line friendly https://cran.rstudio.com/web/packages/rgeos/rgeos.pdf
  
  type.of.A <- class(A)[1]
  # calculate fraction of A in the overlap region
  if (type.of.A == "SpatialLines") {
    #use length if A is a line
    overlap <- rgeos::gLength(overlap.region) / rgeos::gLength(A) 
  } else if (type.of.A == "SpatialPolygons" | type.of.A == "SpatialPolygonsDataFrame") {
    # use area if A is a polygon
    overlap <- slot(overlap.region@polygons[[1]], "area") / slot(A@polygons[[1]],"area")
  }
  eps_diff <- abs(f1 - overlap)
  i <- 2
  k <- 1
  # expand the buffer width until |f1 - fraction of A in the overlap| is within the tolerance
  while (eps_diff > tol) {
    if (abs(f1 - overlap) < 10^(-i) || k > 10) {
      i <- i + 1 
      k <- 1
    }
    epsilon <- epsilon * (1 + sign(f1 - overlap) * 10^(-i))
    
    buff <- buffer(B, width=epsilon, dissolve=T)
    overlap.region <- gIntersection(buff, A);
    if (type.of.A == "SpatialLines") {
      #use length if A is a line
      overlap <- rgeos::gLength(overlap.region) / rgeos::gLength(A)  
    } else if (type.of.A == "SpatialPolygons" | type.of.A == "SpatialPolygonsDataFrame") {
      # use area if A is a polygon
      overlap <- slot(overlap.region@polygons[[1]], "area") / slot(A@polygons[[1]],"area")
    }
    eps_diff <- abs(f1 - overlap);
    k <- k + 1
  }
  
  #Get buffer coordinates
  if (type.of.A == "SpatialLines") {
    buff.coords <- SpatialPoints(overlap.region@lines[[1]]@Lines[[1]]@coords, proj4string=CRS(proj4string(A))) 
  } else if (type.of.A == "SpatialPolygons" | type.of.A == "SpatialPolygonsDataFrame") {
    buff.coords <- SpatialPoints(slot(overlap.region@polygons[[1]]@Polygons[[1]], "coords"), proj4string=CRS(proj4string(A)))
  }
  epsilon <- max(gDistance(buff.coords, B, byid=T))
  ## visualize how much area was captured?
  return(epsilon)
}

# pointHaus ---------------------------------------------------------------
#' Calculate the extended Hausdorff Distance from a point to an area
#' 
#' This function takes a SpatialPoints object, a SpatialLines/SpatialPolygons objec, and a decimal in [0, 1].
#'  
#'@author ???
#'
#'@param point a SpatialPolints object representing a point
#'@param B region calculate the extended Hausdorff distance: SpatialPolygons or SpatialLines
#'@param f2 the percentage of area in B you want captured as a decimal (eg 10% = .1)
#'@param tol tolerance for selecting the epsilon buffer to yield desired f2. Default is NULL.
#'
#'@return the extended hausdorff distance between point and B
#'
#'Last Edited: July 01 2021
pointHaus <- function(point, B, f2, tol=NULL) {
    # calculate the directed Hausdorff distance between point and B (i.e. the min distance between point and B)
    point_to_B <- gDistance(point, B)
    if (f2 == 0) {
        # if f2=0, h(B, "point") equals the minimum distance between "point"
        # and B which equals the value above.
        return(point_to_B)
    } else {
        # calculate the extended directed Hausdorff distance from B to "point"
        B_to_point <- directHaus(B, point, f2, tol=tol)
        return(max(point_to_B, B_to_point))
    }
}
