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
# Improved modularity and code readability by creating a sequential helper function and a parallel helper function
# Ensured that hausMat can correctly handle the asymmetric case where f1 is not equal to f2
# Added file output capabilities
# Improved speed by changing parallel implementation to minimize the number of parallel tasks created
  # Doing each operation in parallel is probably not worth it as the cost of creating a new thread is likely higher than the benefits of the additional parallelism.
  # Its probably only worth creating as many tasks as there are cores (or rows): use a forallchunked equivalent.
  # From the foreach package vignette: "But for the kinds of quick running operations that we’ve been doing, there wouldn’t be much point to executing them in parallel. 
  # Running many tiny tasks in parallel will usually take more time to execute than running them sequentially, and if it already runs fast, there’s no motivation to 
  # make it run faster anyway. But if the operation that we’re executing in parallel takes a minute or longer, there starts to be some motivation."

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

# pointHaus
# Replaced gDistance(point, B, hausdorff=T, byid=T) with gDistance(point, B)
  # gDistance(point, B, hausdorff=T, byid=T) returns H(point, B) and not
  # h(point, B) (as desired). Since point is a point, h(point, B) is equal to
  # the minimum distance between point and B which can be found using
  # gDistance(point, B)
# Changed variable names to make them consistent with the input parameter names

# Next steps:
# Edit code style: refer to http://adv-r.had.co.nz/Style.html

# hausMat -----------------------------------------------------------------
#' Computes a matrix of extended Hausdorff distances using the
#' input Spatial objects "shp" and the input decimals "f1" and "f2".
#'
#'@param shp a Spatial*DataFrame object with n rows.
#'@param f1 The percentage (as a decimal) of region i to retain when
#'   calculating the directional Hausdorff distance from region i to region j.
#'   E.g. 10% = 0.1.
#'@param f2 The percentage (as a decimal) of region j to retain when
#'   calculating the directional Hausdorff distance from region j to i. 
#'   E.g. 10% = 0.1. Defaults to the value of f1. Note that specifying a 
#'   value of f2 that does equal f1 will result in a non-symmetric matrix.
#'@param fileout Should the resulting weight matrix be written to file? 
#'   Defaults to FALSE.
#'@param filename If "fileout" is TRUE, the name for the file to be outputted.
#'@param ncores If "do.parallel" is true, "ncores" is the number of cores to
#'   be used for  parallel computation. Defaults to NULL, although this value
#'   is then updated to detectCores() - 1 if "do.parallel" is true.
#'@param timer If TRUE, records "Sys.time()" when the function is called and 
#'   outputs the elapsed time along with the matrix. Default is FALSE.
#'@param do.parallel If TRUE, uses the number of cores specified by "ncores"
#'   to set up a cluster with the "foreach" package. Defaults to TRUE.
#'@param tol a tolerance value to be passed onto extHaus. Defaults to NULL.  
#'  
#'@return an nxn matrix of extended Hausdorff distances.
#'
#' Last edited: July 24 2020
hausMat <- function(shp, f1, f2 = f1, fileout = FALSE, filename = NULL,
                    ncores = NULL, timer = FALSE, do.parallel = TRUE, tol = NULL) {
  if (timer) {
    start <- Sys.time()
  }
  if (do.parallel) {
    # compute the extended Hausdorff distance in parallel
    if (is.null(ncores)) {
      ncores <- detectCores() - 1 # calculate default number of cores to use
    }
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
    # write n x n matrix of Hausdorff distances to the file specified
    write(t(dists), file = filename)
  }
  return(dists)
}

# seq_haus_mat
#' Computes a matrix of extended Hausdorff distances using the
#' input Spatial objects "shp" and the input decimals "f1" and "f2".
#'
#'@param shp a Spatial*DataFrame object with n rows.
#'@param f1 The percentage (as a decimal) of region i to retain when
#'   calculating the directional Hausdorff distance from region i to region j.
#'   E.g. 10% = 0.1.
#'@param f2 The percentage (as a decimal) of region j to retain when
#'   calculating the directional Hausdorff distance from region j to i. 
#'   E.g. 10% = 0.1. Defaults to the value of f1. Note that specifying a 
#'   value of f2 that does equal f1 will result in a non-symmetric matrix.
#'@param tol a tolerance value to be passed onto extHaus. Defaults to NULL.  
#'  
#'@return an nxn matrix of extended Hausdorff distances.
#'
#' Last edited: July 24 2020
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
    # must compute the extended Hausdorff distance for each entry in the
    # output matrix.
    for (i in 1:n) {
      for (j in 1:n) {
        dists[i, j] <- extHaus(shp[i, ], shp[j, ], f1 = f1, f2 = f2, tol = tol)
      }
    }
  }
  return(dists)
}

# Deprecated
par_haus_mat_slow <- function(shp, f1, f2 = f1, ncores, tol = NULL) {
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
}

# Modified version of par_haus_mat_slow that attempts to minimize the number 
# of parallel tasks created by only creating  "ncores" tasks and dividing 
# the work evenly across tasks
# par_haus_mat
#' Computes a matrix of extended Hausdorff distances using the
#' input Spatial objects "shp" and the input decimals "f1" and "f2".
#'
#'@param shp a Spatial*DataFrame object with n rows.
#'@param f1 The percentage (as a decimal) of region i to retain when
#'   calculating the directional Hausdorff distance from region i to region j.
#'   E.g. 10% = 0.1.
#'@param f2 The percentage (as a decimal) of region j to retain when
#'   calculating the directional Hausdorff distance from region j to i. 
#'   E.g. 10% = 0.1. Defaults to the value of f1. Note that specifying a 
#'   value of f2 that does equal f1 will result in a non-symmetric matrix.
#'@param ncores The number of cores to be used for  parallel computation.
#'@param tol a tolerance value to be passed onto extHaus. Defaults to NULL.  
#'  
#'@return an nxn matrix of extended Hausdorff distances.
#'
#' Last edited: July 24 2020
par_haus_mat <- function(shp, f1, f2 = f1, ncores, tol = NULL) {
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
}

# extHaus -----------------------------------------------------------------
#' Calculates the extended Hausdorff distance between two spatial objects
#' 
#'@author Julia Schedler
#'
#'@param A region/line/point used to calculate the extended Hausdorff distance.
#'   Must be Spatial* or Spatial*DataFrame.
#'@param B region/line/point used to calculate the extended Hausdorff distance. 
#'   Must be Spatial* or Spatial*DataFrame.
#'@param f1 the percentage of area or length in A you want captured as a 
#'   decimal. E.g. 10% = 0.1.
#'@param f2 the percentage of area or length in B you want captured as a 
#'   decimal. E.g. 10% = 0.1. Defaults to the same value as f1.
#'@param tol value to be passed to the parameter "tol" in "directHaus".
#'   Defaults to NULL.
#'   
#'@return the extended Hausdorff distance between "A" and "B" using the
#'   specified decimal values for "f1" and "f2".
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
#' This function calculates the directional extended Hausdorff distance from 
#' "A" to "B" using the input decimal "f1".
#' 
#'@author Julia Schedler
#'
#'@param A region/line used to calculate the directed extended Hausdorff 
#'   distance: SpatialPolygons, SpatialLines or SpatialLinesDataFrame or
#'   SpatialPolygonsDataFrame.
#'@param B region/line/point used to calculate the directed extended Hausdorff 
#'   distance: Spatial* or Spatial*DataFrame.
#'@param f1 the percentage of area or length in A you want captured as a 
#'   decimal. E.g. 10% = 0.1.
#'@param tol tolerance for selecting the epsilon buffer to yield desired f1. 
#'   Default is NULL but this value is then updated to be 1/10000th of the 
#'   sampled directional distances.
#'
#'@return The directional extended Hausdorff distance from A to B.
#'
#' Last Edited: July 22 2021
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
  
  type_of_A <- class(A)[1]
  A_is_lines <- (type_of_A == "SpatialLines") || 
    (type_of_A == "SpatialLinesDataFrame")
  A_is_polygons <- (type_of_A == "SpatialPolygons") || 
    (type_of_A == "SpatialPolygonsDataFrame")
  
  if (A_is_lines) {
    size_A <- rgeos::gLength(A)
  } else if (A_is_polygons) {
    size_A <- slot(A@polygons[[1]],"area")
  } else {
    stop("A is not SpatialLines, SpatialLinesDataFrame, SpatialPolygons, SpatialPolygonsDataFrame")
  }
  
  # sample points from A
  n_pts <- 10000
  a_coords <- sp::spsample(A, n = n_pts, type = "regular")
  # compute minimum distance of each point in a_coords to a point in B
  dists <- rgeos::gDistance(a_coords, B, byid = T)
  
  if (is.null(tol)) {
    # choose a default tolerance if none is specified
    tol <- sd(dists[1, ]) / 200
  }
  
  # find desired quantile of distances
  epsilon <- as.numeric(quantile(dists[1, ], f1)) #buffer width
  buff <- raster::buffer(B, width = epsilon, dissolve = T)
  overlap_region <- rgeos::gIntersection(buff, A) 
  # gIntersection returns the overlap region of buffer + B with A
  # returns null if buff and A do not intersect
  
  
  # if the buffer isn't wide enough to create a region that overlaps with A
  # increment epsilon until buff and A overlap
  while (is.null(overlap_region)) {
    epsilon <- epsilon * 1.01
    buff <- buffer(B, width = epsilon, dissolve = T)
    overlap_region <- rgeos::gIntersection(buff, A)
  }
  
  # calculate fraction of A in the overlap region
  if (A_is_lines) {
    #use length if A is a line
    overlap <- rgeos::gLength(overlap_region) / size_A
  } else if (A_is_polygons) {
    # use area if A is a polygon
    overlap <- slot(overlap_region@polygons[[1]], "area") / size_A
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
    # iterations increase and we approach "tol"
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
      overlap <- rgeos::gLength(overlap_region) / size_A
    } else if (A_is_polygons) {
      # use area if A is a polygon
      overlap <- slot(overlap_region@polygons[[1]], "area") / size_A
    } else {
      stop("A is neither SpatialLines, SpatialLinesDataFrame, SpatialPolygons or SpatialPolygonsDataFrame")
    }
    eps_diff <- abs(f1 - overlap);
    k <- k + 1
  }
  
  if (A_is_lines) {
    overlap_coords <- sp::spsample(overlap_region, n = n_pts, type = "regular")
  } else if (A_is_polygons) {
    overlap_coords <- sp::SpatialPoints(
      slot(overlap_region@polygons[[1]]@Polygons[[1]], "coords"),
      proj4string = CRS(proj4string(overlap_region))
    )
  }
  dists <- rgeos::gDistance(overlap_coords, B, byid = T)
  epsilon <- max(dists)
  
  return(epsilon)
}

# pointHaus ---------------------------------------------------------------
#' Calculate the extended Hausdorff distance from a point to an area or line
#'  
#'@author ???
#'
#'@param point a SpatialPolints object representing a point.
#'@param B region used to calculate the extended Hausdorff distance. Must be a
#'   SpatialPolygons, SpatialLines, SpatialLinesDataFrame or 
#'   SpatialPolygonsDataFrame.
#'@param f2 the percentage of area or length in B you want captured as a 
#'   decimal. E.g. 10% = 0.1.
#'@param tol tolerance for selecting the epsilon buffer to yield desired f2.
#'   Defaults to NULL.
#'
#'@return the extended Hausdorff distance between point and B.
#'
#' Last Edited: July 5 2021
pointHaus <- function(point, B, f2, tol = NULL) {
  # h("point", B) is the minimum cartesian distance between "point" and B
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
