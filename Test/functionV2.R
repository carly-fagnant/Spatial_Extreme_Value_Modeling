# Libraries ---------------------------------------------------------------
library(doParallel)
library(foreach)
library(ggsn)
library(raster)
library(rgdal)
library(rgeos)
library(spatialreg)
library(spdep)

# Notes on key changes ----------------------------------------------------

#hausMat
# Testing: ensure that it can correctly handle the asymmetric case where f1 is not equal to f2
# Testing: ensure that file output feature works correctly
# Improved modularity and code readability by creating a sequential helper function and a parallel helper function
# Doing each operation in parallel is probably not worth it as the cost of creating a new thread is likely higher than the benefits of the additional parallelism.
# Its probably only worth creating as many tasks as there are cores (or rows): use a forallchunked equivalent.
# From the foreach package vignette: "But for the kinds of quick running operations that we’ve been doing, there wouldn’t be much point to executing them in parallel. 
# Running many tiny tasks in parallel will usually take more time to execute than running them sequentially, and if it already runs fast, there’s no motivation to 
# make it run faster anyway. But if the operation that we’re executing in parallel takes a minute or longer, there starts to be some motivation."
# Added file output capabilities.

# extHaus
# Added checks to determine the type of A and B in order to decide how to perform the extended hausdorff distance calculations
# Extended the function's capabilities to handle {point, point}, {point, line}, {point, area}, {line, line}, {line, area}, and {area, area} calculations
# In the cases where f1 == 1 xor f2 == 1, replaced the calls to gDistance(A, B, hausdorff=T) with calls to directHaus(A, B, 1) and directHaus(B, A, 1) (respectively)
  # gDistance(A, B, hausdorff=T) returns the value of H(A, B) whereas in the cases above we are actually interested in h(A, B) (when f1 == 1) and h(B, A) (when f2 == 1)
  # Since H(A, B) = max(h(A, B), h(B, A)), it can be incorrect to use gDistance(A, B, hausdorff=T) in either of these cases.
# In the case of f1 = f2 = 1, replaced the two calls to directHaus with a single call to gDistance(A, B, hausdorff=T) as it is faster and more accurate

# directHaus
# Removed the mostly unused parameter f2
# Allowed directHaus to perform computations when f1=1 instead of throwing an error
# Extended the function's capabilities to handle line-to-area, line-to-line, line-to-point, area-to-area, area-to-line, and area-to-point calculations
# The output was originally a list, of one number. It now directly outputs a numeric value so that it doesn't have to be accessed by $direct.haus[1] in the other functions.
# Added a check for the case where the two input regions are the same.
# Added a case for f1=0 (simply return the minimum distance between A and B instead of performing the epsilon-buffer procedure)
# Ensured that the function does not fail for small values of f1 (e.g. f1 < 0.001) by adding a is.null check on the buffer + A overlap region.
# Question: inside the main while loop, should we increment k twice if we have to increase the buffer size to avoid zero overlap between the buffer and A? 
  # In other words, do we need the 'k <- k + 1' in line 396
# Question: should the default tolerance be sd(dists) / 100?

# pointHaus
# Replaced gDistance(point, B, hausdorff=T, byid=T) with gDistance(point, B)
  # gDistance(point, B, hausdorff=T, byid=T) returns H(point, B) and not
  # h(point, B) (as desired). Since point is a point, h(point, B) is equal to
  # the minimum distance between point and B which can be found using
  # gDistance(point, B)
# Changed variable names to make them consistent with the input parameter names

# Next steps:
# Accuracy problems with extHaus
  # High % error when working on line (rivers) data (particularly at extreme values of max(f1, f2))
  # Is a different testing approach necessary or does the method used in directHaus need change?
  # Increasing number of points sampled did not fix the problem
  # Decreasing the tolerance (by 1/2) did not fix the problem
  # Could try changing the way the length is measured?
# Address directHaus questions and make modifications accordingly
  # Test directHaus using a greater number of sampled points
  # Test directHaus using a different default tolerance
# Test hausMat
  # Finish testing hausMat correctness
  # Finish testing file I/O capabilities
# Re-run entire test suite (functionV2Testing.R)
# Edit code style: refer to http://adv-r.had.co.nz/Style.html

# hausMat -----------------------------------------------------------------
#' Creates a matrix of (extended) Hausdorff distances
#' 
#' This function computes a matrix of extended Hausdorff distances using the
#' "shp" and the input decimal "f1".
#'
#'@param shp an nxn matrix of Spatial* objects or a Spatial*DataFrame object
#'   with n rows.
#'@param f1 The percentage (as a decimal) of region i to retain when
#'   calculating the directional Hausdorff distance from region i to region j.
#'@param f2 The percentage (as a decimal) of region j to retain when
#'   calculating the directional Hausdorff distance from region j to i. 
#'   Defaults to the value of f1. Note that specifying a different  value will
#'   result in a non-symmetric matrix.
#'@param fileout Should the resulting weight matrix be written to file? 
#'   Defaults to FALSE 
#'@param filename If "fileout" is TRUE, the name for the file to be outputted.
#'@param ncores If "do.parallel" is true, ncores is the number of cores to be
#'   used for  parallel computation. Defaults to 1. 
#'@param timer If T, records "Sys.time()" when the function is called and 
#'   outputs the elapsed time along with the matrix. Default is F.
#'@param do.parallel If TRUE, uses the number of cores specified using "ncores"
#'   to set up a cluster with the "foreach" package.
#'@param tol a tolerance value to be passed onto extHaus. Defaults to NULL.  
#'  
#'@return an nxn matrix of extended Hausdorff distances
#'
#' Last edited: July 8 2020
hausMat <- function(shp, f1, f2 = f1, fileout = FALSE, filename = NULL,
                    ncores = 1, timer = F, do.parallel = T, tol = NULL) {
  if (timer) {
    start <- Sys.time()
  }
  if (do.parallel) {
    # compute the extended Hausdorff distance in parallel
    dists <- par_haus_mat(shp, f1 = f1, f2 = f2, ncores = ncores, tol = tol)
  } else {
    # compute the extended Hausdorff distance sequentially
    dists <- seq_haus_mat(shp, f1 = f1, f2 = f2, tol = tol)
  }
  if (timer) {
    print("Completion time:")
    print(Sys.time() - start)
  }
  if (fileout) {
    # write nxn matrix to filename
    write(t(dists), file = filename)
  }
  return(dists)
}

seq_haus_mat <- function(shp, f1, f2 = f1, tol = NULL) {
  n <- nrow(shp@data)
  dists <- matrix(0, nrow = n, ncol = n)
  if (f1 == f2) {
    combs <- combn(1:n, 2) # n choose 2 combinations
    n_combs <- ncol(combs)
    for (i in 1:n_combs) {
      # compute the extended Hausdorff distance for each possible combination
      # of indices
      index1 <- combs[1, i]
      index2 <- combs[2, i] 
      dists[index1, index2] <- extHaus(
        shp[index1, ],
        shp[index2, ],
        f1 = f1,
        tol = tol
      )
    }
    # since f1 = f2 the matrix is symmetric and so we can compute the missing
    # entries of dist by adding its transpose to itself.
    dists <- dists + t(dists)
  } else {
    # since f1 does not equal f2 the matrix may not be symmetric and so we
    # must compute the extended Hausdorff distance for each entry
    for (i in 1:n) {
      for (j in 1:n) {
        dists[i, j] <- extHaus(shp[i, ], shp[j, ], f1 = f1, f2 = f2, tol = tol)
      }
    }
  }
  return(dists)
}

par_haus_mat_slow <- function(shp, f1, f2 = f1, ncores = 1, tol = NULL) {
  #start <- Sys.time() used only for time testing
  n <- nrow(shp@data)
  dists <- matrix(0, nrow = n, ncol = n)
  #print("Setting up parallelization")
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  #print("Computing...")
  # compute the extended Hausdorff distance for every combination of regions
  if (f1 == f2) {
    combs <- combn(1:n, 2) # n choose 2 combinations
    n_combs <- ncol(combs)
    out <- foreach(
      i = 1:n_combs,
      .packages = c("rgeos", "sp", "raster"),
      .combine = rbind,
      .export = c("directHaus", "extHaus")
    ) %dopar% {
      extHaus(shp[combs[1, i], ], shp[combs[2, i], ], f1 = f1, tol = tol)
    } 
  } else {
    out <- foreach (
      i = 1:n,
      .packages = c("rgeos", "sp", "raster", "doParallel", "foreach"),
      .combine = rbind,
      .export = c("directHaus", "extHaus")
    ) %dopar% {
      foreach (
        j = 1:n,
        .packages = c("rgeos", "sp", "raster"),
        .combine = rbind,
        .export=c("directHaus", "extHaus")
      ) %dopar% {
        extHaus(shp[i, ], shp[j, ], f1 = f1, f2 = f2, tol = tol) 
      }
    }
  }
  stopCluster(cl)
  print("Parallel computation complete")
  if (f1 == f2) {
    dists[lower.tri(dists, diag = F)] <- out
    # use the fact that the Hausdorff distance matrix will be symmetric to 
    # compute the upper triangular entries
    dists <- dists + t(dists)
  } else {
    dists <- matrix(out, nrow = n, byrow = T)
  }
  return(dists)
  #return(Sys.time() - start) used only for time testing
}

# Modified version of parHausMat_slow that attempts to minimize the number of 
# parallel tasks created by only creating  "ncores" tasks and dividing 
# the work evenly across tasks
par_haus_mat <- function(shp, f1, f2 = f1, ncores = 1, tol = NULL) {
  #start <- Sys.time() used only for time testing
  n <- nrow(shp@data)
  dists <- matrix(0, nrow = n, ncol = n)
  print("Setting up parallelization")
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  print("Computing...")
  # compute the extended Hausdorff distance for every combination of regions
  if (f1 == f2) {
    combs <- combn(1:n, 2) # n choose 2 combinations
    n_combs <- ncol(combs)
    # calculate how many sequential iterations we need to do per core to
    # avoid creating more parallel tasks than cores
    iterPerCore <- ceiling(n_combs / ncores)
    out <- foreach(
      k = 1:ncores,
      .packages = c("rgeos", "sp", "raster"),
      .combine = "c",
      .export = c("directHaus", "extHaus")
    ) %dopar% {
      val <- rep(-1, iterPerCore)
      for (iter in 1:iterPerCore) {
        remaining <- n_combs - ((k - 1) * iterPerCore) - (iter - 1)
        if (remaining > 0) {
          i <- ((k - 1) * iterPerCore) + iter
          val[iter] <- extHaus(
            shp[combs[1, i], ],
            shp[combs[2, i], ],
            f1 = f1,
            tol = tol
          ) 
        }
      }
      val
    }
  } else {
    # calculate how many sequential iterations we need to do per core to 
    # avoid creating more parallel tasks than cores
    iterPerCore <- ceiling(n^2 / ncores)
    out <- foreach(
      k = 1:ncores,
      .packages = c("rgeos", "sp", "raster"),
      .combine = "c",
      .export = c("directHaus", "extHaus")
    ) %dopar% {
      val <- rep(-1, iterPerCore)
      for (iter in 1:iterPerCore) {
        remaining <- n^2 - ((k - 1) * iterPerCore) - (iter - 1)
        if (remaining > 0) {
          index <- ((k - 1) * iterPerCore) + iter
          i <- ceiling(index / n)
          j <- index %% n
          if (index %% n == 0) {
            j <- n
          }
          val[iter] <- extHaus(shp[i, ], shp[j, ], f1 = f1, f2 = f2, tol = tol) 
        }
      }
      val
    }
  }
  stopCluster(cl)
  print("Parallel computation complete")
  if (f1 == f2) {
    dists[lower.tri(dists, diag = F)] <- out[1:n_combs]
    # use the fact that the Hausdorff distance matrix will be symmetric to
    # compute the upper triangular entries
    dists <- dists + t(dists)
  } else {
    dists <- matrix(out[1:n^2], nrow = n, byrow = T)
  }
  return(dists)
  #return(Sys.time() - start) used only for time testing
}

# extHaus -----------------------------------------------------------------
#' Calculates the extended Hausdorff distance between two spatial objects
#'
#' Calculates the extended Hausdorff distance from "A" to "B" based on the
#' input decimals "f1" and "f2"
#' 
#'@author Julia Schedler
#'
#'@param A region/line/point used to calculate the extended Hausdorff distance.
#'   Must be Spatial* or Spatial*DataFrame
#'@param B region/line/point used to calculate the extended Hausdorff distance. 
#'   Must be Spatial* or Spatial*DataFrame
#'@param f1 the percentage of area or length in A you want captured as a 
#'   decimal. E.g. 10% = 0.1
#'@param f2 the percentage of area or length in B you want captured as a 
#'   decimal. E.g. 10% = 0.1. Defaults to the same value as f1.
#'@param tol value to be passed to the parameter "tol" in "directHaus".
#'
#'@return the extended Hausdorff distance (i.e. the max of the directional 
#'   Hausdorff distance from A to B and from B to A)
#'
#' Last Edited: July 14 2021
extHaus <- function(A, B, f1, f2 = f1, tol = NULL) {
  if (!sp::is.projected(A) | !sp::is.projected(B)) {
    stop(paste("Spatial* object (inputs ", quote(A), ", ", quote(B),
               ") must be projected. Try running ?spTransform().", sep = ""))
  }
  if (is.null(gDifference(A, B))) {
    # if two regions are the same return a Hausdorff distance of zero.
    return(0) 
  } 
  
  # get the class of A and B to perform calculations accordingly
  type_of_A <- class(A)[1]
  type_of_B <- class(B)[1]
  A_is_points <- (type_of_A == "SpatialPoints") ||
                 (type_of_A == "SpatialPointsDataFrame")
  B_is_points <- (type_of_B == "SpatialPoints") ||
                 (type_of_B == "SpatialPointsDataFrame")
  A_is_lines <- (type_of_A == "SpatialLines") ||
                (type_of_A == "SpatialLinesDataFrame")
  B_is_lines <- (type_of_B == "SpatialLines") ||
                (type_of_B == "SpatialLinesDataFrame")
  A_is_polygons <- (type_of_A == "SpatialPolygons") ||
                   (type_of_A == "SpatialPolygonsDataFrame")
  B_is_polygons <- (type_of_B == "SpatialPolygons") ||
                   (type_of_B == "SpatialPolygonsDataFrame")
  
  if (f1 == 1 && f2 == 1) {
    # if f1 = f2 = 1 then we are simply interested in the Hausdorff distance
    # between A and B. Delegate to gDistance as it is faster than directHaus
    return(rgeos::gDistance(A, B, hausdorff = T))
  }
  
  if (A_is_points || B_is_points) {
    if (A_is_points && B_is_points) {
      # if both A and B are points, the extended Hausdorff distance between
      # them is the minimum Cartesian distance between them
      return(rgeos::gDistance(A, B)) 
    } else if (A_is_points) {
      # if one of A or B are points, delegate to pointHaus
      return(pointHaus(A, B, f2, tol = tol))
    } else {
      # if one of A or B are points, delegate to pointHaus
      return(pointHaus(B, A, f1, tol = tol))
    }
  } 
  
  if (A_is_lines || A_is_polygons) {
    if (f1 == 0) {
      # if f1 = 0 use the Cartesian minimum distance between A and B
      A_to_B <- rgeos::gDistance(A, B)
    } else {
      # use the extended directed Hausdorff distance from A to B
      A_to_B <- directHaus(A, B, f1, tol = tol)
    }
  } else {
    stop("A is not SpatialPoints, SpatialLines, SpatialPolygons, SpatialPointsDataFrame, SpatialLinesDataFrame or SpatialPolygonsDataFrame")
  }
  
  if (B_is_lines || B_is_polygons) {
    if (f2 == 0) {
      # if f2 = 0 use the Cartesian minimum distance between A and B
      B_to_A <- rgeos::gDistance(A, B)
    } else {
      # use the directed extended Hausdorff distance from B to A
      B_to_A <- directHaus(B, A, f2, tol = tol)
    }
  } else {
    stop("B is not SpatialPoints, SpatialLines, SpatialPolygons, SpatialPointsDataFrame, SpatialLinesDataFrame or SpatialPolygonsDataFrame")
  }
  dist <- max(A_to_B, B_to_A)
  return(dist)
}

# directHaus --------------------------------------------------------------
#' Calculates the directional extended Hausdorff distance
#' 
#' This function calculates the directional extended Hausdorff distance from 
#' "A" to "B" using the input decimal f1
#' 
#'@author Julia Schedler
#'
#'@param A region/line used to calculate the directed extended Hausdorff 
#'   distance: SpatialPolygons, SpatialLines or SpatialLinesDataFrame or
#'   SpatialPolygonsDataFrame
#'@param B region/line/point used to calculate the directed extended Hausdorff 
#'   distance: Spatial* or Spatial*DataFrame
#'@param f1 the percentage of area or length in A you want captured as a 
#'   decimal. E.g. 10% = 0.1
#'@param tol tolerance for selecting the epsilon buffer to yield desired f1. 
#'   Default is 1/10000th the sampled directional distances.
#'
#'@return The directional extended Hausdorff distance from A to B
#'
#' Last Edited: July 14 2021
directHaus <- function(A, B, f1, tol = NULL) {
  if (!sp::is.projected(A) || !sp::is.projected(B)) {
    stop(paste("Spatial* object (inputs ", quote(A),", ", quote(B),
               ") must be projected. Try running ?spTransform().", sep = ""))
  }
  
  if (is.null(rgeos::gDifference(A,B))) {
    # check if two regions are the same, if so return distance of zero
    return(0) 
  }

  if (f1 == 0) {
    # if f1 = 0, return the minimum euclidean distance between A and B
    return(rgeos::gDistance(A, B))
  }
  
  # generate points
  n <- 10000
  a_coords <- sp::spsample(A, n = n, type = "regular") # sample points within A
  # compute minimum distance of points a.coords to a point in B
  dists <- rgeos::gDistance(a_coords, B, byid = T)
  
  if (is.null(tol)) {
    # choose a default tolerance if none is specified
    tol <- sd(dists[1, ]) / 100
  }
  
  # find desired quantile of distances
  epsilon <- as.numeric(quantile(dists[1, ], f1)) #buffer width
  buff <- raster::buffer(B, width = epsilon, dissolve = T)
  overlap_region <- rgeos::gIntersection(buff, A) 
  # gIntersection returns the overlap region of buffer+B with A
  # returns null if buff and A do not intersect
  
  # if the buffer isn't wide enough to create a region that overlaps with A
  # increment epsilon until buff and A overlap
  while (is.null(overlap_region)) {
    epsilon <- epsilon * 1.01
    buff <- buffer(B, width = epsilon, dissolve = T)
    overlap_region <- rgeos::gIntersection(buff, A)
  }
  
  type_of_A <- class(A)[1]
  A_is_lines <- (type_of_A == "SpatialLines") || 
                (type_of_A == "SpatialLinesDataFrame")
  A_is_polygons <- (type_of_A == "SpatialPolygons") || 
                   (type_of_A == "SpatialPolygonsDataFrame")
  
  # calculate fraction of A in the overlap region
  if (A_is_lines) {
    #use length if A is a line
    overlap <- rgeos::gLength(overlap_region) / rgeos::gLength(A) 
  } else if (A_is_polygons) {
    # use area if A is a polygon
    overlap <- slot(overlap_region@polygons[[1]], "area") / 
      slot(A@polygons[[1]],"area")
  } else {
    stop("A is neither SpatialLines, SpatialLinesDataFrame, SpatialPolygons or SpatialPolygonsDataFrame")
  }
  eps_diff <- abs(f1 - overlap)
  
  # expand/contract the buffer width until |f1 - fraction of A in the overlap| 
  # is within the specified tolerance
  i <- 2
  k <- 1
  while (eps_diff > tol) {
    # shrink the amount by which the buffer width is modified as the number of
    # iterations increase and we approach tol
    if (abs(f1 - overlap) < 10^(-i) || k > 10) {
      i <- i + 1 
      k <- 1
    }
    epsilon <- epsilon * (1 + sign(f1 - overlap) * 10^(-i))
    buff <- raster::buffer(B, width = epsilon, dissolve = T)
    overlap_region <- rgeos::gIntersection(buff, A)
    
    # if the buffer isn't wide enough to create a region that overlaps with A
    # increment epsilon until buff and A overlap
    while (is.null(overlap_region)) {
      epsilon <- epsilon * (1 + 10^(-i))
      buff <- buffer(B, width = epsilon, dissolve = T)
      overlap_region <- rgeos::gIntersection(buff, A)
      k <- k + 1
    }
    
    # update eps_diff
    if (A_is_lines) {
      #use length if A is a line
      overlap <- rgeos::gLength(overlap_region) / rgeos::gLength(A)  
    } else if (A_is_polygons) {
      # use area if A is a polygon
      overlap <- slot(overlap_region@polygons[[1]], "area") / 
                 slot(A@polygons[[1]],"area")
    }
    eps_diff <- abs(f1 - overlap);
    k <- k + 1
  }
  
  #Get buffer coordinates
  if (A_is_lines) {
    buff_coords <- SpatialPoints(
      overlap_region@lines[[1]]@Lines[[1]]@coords, 
      proj4string = CRS(proj4string(A))
    ) 
  } else if (A_is_polygons) {
    buff_coords <- sp::SpatialPoints(
      slot(overlap_region@polygons[[1]]@Polygons[[1]], "coords"),
      proj4string = CRS(proj4string(A))
    )
  }
  epsilon <- max(rgeos::gDistance(buff_coords, B, byid = T))
  return(epsilon)
}

# pointHaus ---------------------------------------------------------------
#' Calculate the extended Hausdorff distance from a point to an area or line
#'  
#' Calculate the extended Hausdorff distance from "point" to "B" using "f2"  
#'  
#'@author ???
#'
#'@param point a SpatialPolints object representing a point
#'@param B region used to calculate the extended Hausdorff distance. Must be a
#'   SpatialPolygons, SpatialLines, SpatialLinesDataFrame or 
#'   SpatialPolygonsDataFrame
#'@param f2 the percentage of area or length in B you want captured as a 
#'   decimal. E.g. 10% = 0.1
#'@param tol tolerance for selecting the epsilon buffer to yield desired f2.
#'   Defaults to NULL.
#'
#'@return the extended Hausdorff distance between point and B
#'
#' Last Edited: July 5 2021
pointHaus <- function(point, B, f2, tol = NULL) {
  # The directed Hausdorff distance between "point" and B is the minimum
  # cartesian distance between "point" and B
  point_to_B <- rgeos::gDistance(point, B)
  if (f2 == 0) {
    # if f2 = 0, h(B, "point") equals the minimum distance between "point"
    # and B which equals the value above.
    return(point_to_B)
  } else {
    # calculate the extended directed Hausdorff distance from B to "point"
    B_to_point <- directHaus(B, point, f2, tol = tol)
    return(max(point_to_B, B_to_point))
  }
}
