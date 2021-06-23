library(rgdal)
library(doParallel)
library(rgeos)
library(spdep)
library(raster)
library(spatialreg)
library(ggsn)


hausMat <- function(shp, f1, f2=f1, fileout=FALSE, filename=NULL, ncores=1, timer=F, do.parallel=T) {
    if (timer) {
      start <- Sys.time()
    }
    n <- nrow(shp@data)
    combs <- combn(1:n,2)
    n.combs <- ncol(combs);
    haus.dists <- matrix(0, nrow=n, ncol=n)
    if (do.parallel) {
    
        #out <- matrix(-1, nrow = 20, ncol = 20)
        start <- Sys.time()
        print("Setting up parallelization")
        cl <- makeCluster(ncores)
        registerDoParallel(cl)
        print("Computing...")
        out <- foreach (i = 1:n.combs, .packages=c("rgeos", "sp", "raster"), .combine=rbind, .export=c("directHaus", "extHaus")) %dopar% {
            extHaus(shp[combs[1,i],], shp[combs[2,i],], f1 = f1)
        }
        print("Completion time:")
        print(Sys.time()-start)
        #save(out, file = "hausdorff_columbus")
        stopCluster(cl)
        #closeAllConnections()  
        haus.dists[lower.tri(haus.dists, diag=F)] <-  out
        haus.dists <- haus.dists + t(haus.dists)
    } else {
        for (i in 1:n.combs) {
            haus.dists[combs[1,i], combs[2,i]] <- extHaus(shp[combs[1,i],], shp[combs[2,i],], f1=f1)
        }
        haus.dists <- haus.dists + t(haus.dists)
    }
    ## add 0's to diagonal
    #if(fileout){save(haus.dists, file =filename)}
    return(haus.dists)
}


extHaus <- function(A, B, f1, f2=f1, tol = NULL) {
    if (!is.projected(A) | !is.projected(B)) {
        stop(paste("Spatial* object (inputs ", quote(A),", ", quote(B), ") must be projected. Try running ?spTransform().", sep = ""))
    }
    if (is.null(gDifference(A,B))) {
      return(0)
    } ## check if two regions are the same, if so return HD of zero.
  
    if (f1 == 1) {
      A_to_B <- gDistance(A,B, hausdorff=T);
    } else if (f1 == 0) {
      A_to_B <- gDistance(A,B)
    } else {
      A_to_B <- directHaus(A,B, f1, tol=tol)$direct.haus[1]
    }
    if (f2 == 1) {
      B_to_A <- gDistance(A,B, hausdorff=T)
    } else if (f2 == 0) {
      B_to_A <- gDistance(A,B)
    } else {
      B_to_A <- directHaus(B,A, f2, tol=tol)$direct.haus[1]
    }
    haus.dist <- max(A_to_B, B_to_A)
    return(haus.dist)
}


directHaus<- function(A, B, f1, f2=f1, tol=NULL) {
    if (!is.projected(A) | !is.projected(B)) {
        stop(paste("Spatial* object (inputs ", quote(A),", ", quote(B), ") must be projected. Try running ?spTransform().", sep = ""))
    }
    if (f1 == 1 || f2 == 1) {
        stop("f = 1 equivalent to Hausdorff distancce; Use gDistance with hausdorff = T.")
    }
    #A.area = slot(A@polygons[[1]], "area")
    # generate points
    n <- 10000
    a.coords <- spsample(A, n=n,type="regular")
    ## points from A to B
    dists <- gDistance(a.coords, B, byid=T)
  
    if(is.null(tol)){
        tol = sd(dists[1,])/100
    }
    ## find desired quantile of distances
    eps <- as.numeric(quantile(dists[1,], f1))
    buff<- buffer(B, width=eps, dissolve=T)
    overlap.region = gIntersection(buff, A);
    overlap = slot(overlap.region@polygons[[1]], "area")/slot(A@polygons[[1]],"area");
    eps_diff = abs(f1 - overlap)
    i <- 2
    k <- 1
    while (eps_diff > tol){
        if (abs(f1 - overlap) < 10^(-i) || k > 10) {
            i <- i + 1 
            k <- 1
        }
        eps <- eps * (1 + sign(f1 - overlap) * 10^(-i))
    
        buff <- buffer(B, width=eps, dissolve=T)
        overlap.region <- gIntersection(buff, A);
        slot(overlap.region@polygons[[1]], "area")/slot(A@polygons[[1]], "area")
        overlap <- slot(overlap.region@polygons[[1]], "area")/slot(A@polygons[[1]],"area");
        eps_diff <- abs(f1-overlap);
        k <- k + 1
    }
    buff.coords <- SpatialPoints(slot(overlap.region@polygons[[1]]@Polygons[[1]], "coords"), proj4string=CRS(proj4string(A)))
    eps <- max(gDistance(buff.coords, B, byid=T))
    ## visualize how much area was captured?
  
    out <- list(eps);names(out) <- c("direct.haus")
    return(out)
}


pointHaus <- function(point, B, f2, tol=NULL) {
    A_to_B <- gDistance(point, B, hausdorff=T, byid=T)
    if (f2 == 0) {
        return(A_to_B)
    } else {
        B_to_A <- directHaus(B, point, f2, tol=tol)$direct.haus[1]
        return(max(A_to_B, B_to_A))
    }
}
  
  
  
  
