library(rgdal)
library(doParallel)
library(rgeos)
library(spdep)
library(raster)
library(spatialreg)
library(ggsn)

#' Creating a matrix of (extended) Hausdorff distances
#' 
#' This function takes a SpatialPolygons object 
#'  
#'@param shp a SpatialPolygons object.
#'@param f1 The percentage (as a decimal) of region i to retain when calculating the directional Hausdorff distance from region i to region j.
#'@param f2 The percentage (as a decimal) of region j to retain when calculating the directional Hausdorff distance from region j to i. Defaults to the value of f1. Note that specifying a different value will result in a non-symmetric matrix.
#'@param fileout Should the resulting weight matrix be written to file? Defaults to FALSE 
#'@param filename If \code{fileout} is TRUE, the name for the file to be outputted.
#'@param ncores If \code{do.parallel} is true, the number of cores to be used for parallel computation. Default is 1. 
#'@param timer If T, records \code{Sys.time()} when the function is called and outputs the elapsed time along with the matrix. Default is F.
#'@param do.parallel If T, uses the number of cores specified using \code{ncores} to set up a cluster with the \code{foreach} package. 
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

#' Extended Hausdorff Distance
#'
#' Calculates the extended Hausdorff distance from A to B for given quantiles and precision
#' 
#'@author Julia Schedler
#'
#'@param A, B region calculate the extended Hausdorff distance- SpatialPolygons or SpatialPolygonsDataFrame
#'@param f1 the percentage of area in B you want captured as a decimal (eg 10\% = .1)
#'@param f2 the percentage of area in A you want captured as a decimal (eg 10\% = .1)
#'@param tol value to be passed to \code{tol} in \code{\link{directHaus}}.
#'
#'@return haus.dist: the extended hausdorff distance (max of directional from A to B and B to A)
extHaus <- function(A, B, f1, f2=f1, tol=NULL) {
    if (!is.projected(A) | !is.projected(B)) {
        stop(paste("Spatial* object (inputs ", quote(A), ", ", quote(B), ") must be projected. Try running ?spTransform().", sep=""))
    }
    if (is.null(gDifference(A, B))) {
        return(0)
    } ## check if two regions are the same, if so return HD of zero.
  
    if (f1 == 1) {
        # if f1 = 1 use the hausdorff distance between A and B (needs correcting?)
        A_to_B <- gDistance(A, B, hausdorff=T);
    } else if (f1 == 0) {
        # if f1 = 0 use the cartesian minimum distance between A and B
        A_to_B <- gDistance(A, B)
    } else {
        # use the extended directed hausdorff distance from A to B
        A_to_B <- directHaus(A, B, f1, tol=tol)$direct.haus[1]
    }
    if (f2 == 1) {
        # if f2 = 1 use the hausdorff distance between A and B (needs correcting?)
        B_to_A <- gDistance(A, B, hausdorff=T)
    } else if (f2 == 0) {
        # if f2 = 0 use the cartesian minimum distance between A and B
        B_to_A <- gDistance(A, B)
    } else {
        # use the directed extended hausdorff distance from B to A
        B_to_A <- directHaus(B, A, f2, tol=tol)$direct.haus[1]
    }
    haus.dist <- max(A_to_B, B_to_A)
    return(haus.dist)
}

#' Calculate the directional extended Hausdorff Distance
#' 
#' This function takes a SpatialPolygons object 
#'  
#'@author Julia Schedler
#'
#'@param A, B regions used to calculate the directed extended Hausdorff distance: SpatialPolygons or SpatialPolygonsDataFrame
#'@param f1 the percentage of area in B you want captured as a decimal (eg 10\% = .1)
#'@param f2 the percentage of area in A you want captured as a decimal (eg 10\% = .1); defaults to the value of f1
#'@param tol tolerance for selecting the epsilon buffer to yield desired f1. Default is 1/10000th the sampled directional distances.
#'
#'@return The directional extended hausdorff distance from A to B
directHaus<- function(A, B, f1, f2=f1, tol=NULL) {
    if (!is.projected(A) | !is.projected(B)) {
        stop(paste("Spatial* object (inputs ", quote(A),", ", quote(B), ") must be projected. Try running ?spTransform().", sep = ""))
    }
    if (f1 == 1 || f2 == 1) {
        stop("f = 1 equivalent to Hausdorff distance; Use gDistance with hausdorff = T.")
    }
    #A.area = slot(A@polygons[[1]], "area")
    # generate points
    n <- 10000
    a.coords <- spsample(A, n=n, type="regular") # sample points within A
    # compute minimum distance of points a.coords to a point in B
    dists <- gDistance(a.coords, B, byid=T)
  
    if(is.null(tol)){
        tol <- sd(dists[1,]) / 100
    }
    ## find desired quantile of distances
    epsilon <- as.numeric(quantile(dists[1,], f1)) #buffer width
    buff <- buffer(B, width=epsilon, dissolve=T)
    overlap.region <- gIntersection(buff, A) # overlap region of buffer+B with A
    overlap <- slot(overlap.region@polygons[[1]], "area") / slot(A@polygons[[1]],"area") # calculate fraction of A in the overlap region
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
        slot(overlap.region@polygons[[1]], "area") / slot(A@polygons[[1]], "area")
        overlap <- slot(overlap.region@polygons[[1]], "area") / slot(A@polygons[[1]], "area");
        eps_diff <- abs(f1 - overlap);
        k <- k + 1
    }
    buff.coords <- SpatialPoints(slot(overlap.region@polygons[[1]]@Polygons[[1]], "coords"), proj4string=CRS(proj4string(A)))
    epsilon <- max(gDistance(buff.coords, B, byid=T))
    ## visualize how much area was captured?
  
    out <- list(epsilon)
    names(out) <- c("direct.haus")
    return(out)
}

#' Calculate the extended Hausdorff Distance from a point to an area
#' 
#' This function takes a SpatialPolygons object 
#'  
#'@author ???
#'
#'@param point a SpatialPolygons object representing a point
#'@param B region calculate the extended Hausdorff distance- SpatialPolygons or SpatialPolygonsDataFrame
#'@param f2 the percentage of area in B you want captured as a decimal (eg 10\% = .1)
#'@param tol tolerance for selecting the epsilon buffer to yield desired f2. Default is NULL.
#'
#'@return The directional extended hausdorff distance from point to B
pointHaus <- function(point, B, f2, tol=NULL) {
    # calculate the Hausdorff distance between A and B (needs updating?)
    A_to_B <- gDistance(point, B, hausdorff=T, byid=T) # replace with gDistance(point, B, byid=T)?
    if (f2 == 0) {
        return(A_to_B)
    } else {
        # calculate the extended directed Hausdorff distance from B to point using f2
        B_to_A <- directHaus(B, point, f2, tol=tol)$direct.haus[1]
        return(max(A_to_B, B_to_A))
    }
}


