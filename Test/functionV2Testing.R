# A file for testing the functions in functionV2.R

# Creating Spatial Objects ------------------------------------------------

# Create spatial objects
pt1_coords <- matrix(c(0.5, 0.5), ncol = 2, byrow = T)
pt2_coords <- matrix(c(2, 2), ncol = 2, byrow = T)
l1_coords <- matrix(c(0, 10, 1, 1), ncol = 2, byrow = T)
l2_coords <- matrix(c(3, 2, 4, 10), ncol = 2, byrow = T)
p1_coords <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1, 0, 0), ncol = 2, byrow = T)
p2_coords <- matrix(c(2, 0, 3, 1, 4, 0, 2, 0), ncol = 2, byrow = T)
crdref <- sp::CRS("+init=epsg:2278")

pt1 <- sp::SpatialPoints(pt1_coords, proj4string = crdref)
pt2 <- sp::SpatialPoints(pt2_coords, proj4string = crdref)
l1 <- raster::spLines(l1_coords, crs = crdref)
l2 <- raster::spLines(l2_coords, crs = crdref)
p1 <- raster::spPolygons(p1_coords, crs = crdref)
p2 <- raster::spPolygons(p2_coords, crs = crdref)

inputs <- c(pt1, pt2, l1, l2, p1, p2)

# Project the spatial objects with spTransform
for (i in 1:6) {
  print(i)
  print(sp::is.projected(inputs[[i]]))
  sp::proj4string(inputs[[i]]) <- crdref
  inputs[[i]] <- sp::spTransform(inputs[[i]], proj4string(inputs[[i]]))
  print(sp::is.projected(inputs[[i]]))
  print("------------------------------------------------------")
}

# Testing: Points, Lines, and Polygons ------------------------------------

cases1 <- c(pt1, pt1, pt1, l1, l1, p1)
cases2 <- c(pt2, l1, p1, l2, p1, p2)

test_extHaus <- function(A, B, n_samp, f1, f2 = f1) {
  A_is_points <- (class(A)[1] == "SpatialPoints")
  B_is_points <- (class(B)[1] == "SpatialPoints")
  
  flag <- 0
  val <- extHaus(A, B, f1, f2)
  
  if (is.null(gDifference(A, B))) {
    ans <- 0
  } else {
    if (f1 == 0 || A_is_points) {
      A_to_B <- rgeos::gDistance(A, B)
    } else {
      # sample points from A and get the desired distance quantile. 
      # Use thus to estimate the directed Hausdorff distance from A to B.
      a_coords <- sp::spsample(A, n = n_samp, type = "regular")
      distsf1 <- rgeos::gDistance(a_coords, B, byid = T)
      A_to_B <- quantile(distsf1, probs = f1)
    }
    if (f2 == 0 || B_is_points) {
      B_to_A <- rgeos::gDistance(A, B)
    } else {
      # sample points from B and get the desired distance quantile. 
      # Use thus to estimate the directed Hausdorff distance from B to A.
      b_coords <- sp::spsample(B, n = n_samp, type = "regular")
      distsf2 <- rgeos::gDistance(b_coords, A, byid = T) 
      B_to_A <- quantile(distsf2, probs = f2)
    }
    ans <- max(A_to_B, B_to_A)  
  }
  
  error <- abs(val - ans)
  pct_error <- 100 * error / ans
  if (error > 0.01 * ans) {
    flag <- 1
    print("  Error")
    print(paste0("  Expected ", ans, " but computed ", val))
    print(paste0("  Percentage error: ", pct_error, "%"))
    print(paste0("  f1: ", f1))
    print(paste0("  f2: ", f2))
  }
  return(flag)
}

for (i in 1:6) {
  flag_main <- 0
  print(paste0("Testing input set ", i))
  print(paste0("Class of A: ",  class(cases1[[i]])[1]))
  print(paste0("Class of B: ",  class(cases2[[i]])[1]))
  
  ## Case 1: f1 = 0
  print("Test case 1: f1 = 0")
  flag <- test_extHaus(cases1[[i]], cases2[[i]], 1000000, f1 = 0, f2 = 0.5)
  flag_main <- max(flag_main, flag)

  ## Case 2: f2 = 0
  print("Test case 2: f2 = 0")
  flag <- test_extHaus(cases1[[i]], cases2[[i]], 1000000, f1 = 0.5, f2 = 0)
  flag_main <- max(flag_main, flag)
  
  ## Case 3: f1 = f2 = 0
  print("Test case 3: f1 = f2 = 0")
  flag <- test_extHaus(cases1[[i]], cases2[[i]], 1000000, f1 = 0)
  flag_main <- max(flag_main, flag)
  
  ## Case 4: 0 < f1, f2 < 1
  print("Test case 4: 0 < f1, f2 < 1")
  # Test case 4.1
  flag <- test_extHaus(cases1[[i]], cases2[[i]], 1000000, f1 = 0.25, f2 = 0.75)
  flag_main <- max(flag_main, flag)
  # Test case 4.2
  flag <- test_extHaus(cases1[[i]], cases2[[i]], 1000000, f1 = 0.75, f2 = 0.25)
  flag_main <- max(flag_main, flag)
  
  ## Case 5: 0 < (f1 = f2) < 1
  print("Test case 5: 0 < (f1 = f2) < 1")
  # Test case 5.1
  flag <- test_extHaus(cases1[[i]], cases2[[i]], 1000000, f1 = 0.25)
  flag_main <- max(flag_main, flag)
  # Test case 5.2
  flag <- test_extHaus(cases1[[i]], cases2[[i]], 1000000, f1 = 0.5)
  flag_main <- max(flag_main, flag)
  # Test case 5.3
  flag <- test_extHaus(cases1[[i]], cases2[[i]], 1000000, f1 = 0.75)
  flag_main <- max(flag_main, flag)
  
  ## Case 6: f1 = 1
  print("Test case 6: f1 = 1")
  flag <- test_extHaus(cases1[[i]], cases2[[i]], 1000000, f1 = 1, f2 = 0.5)
  flag_main <- max(flag_main, flag)
  
  ## Case 7: f2 = 1
  print("Test case 7: f2 = 1")
  flag <- test_extHaus(cases1[[i]], cases2[[i]], 1000000, f1 = 0.5, f2 = 1)
  flag_main <- max(flag_main, flag)
  
  ## Case 8: f1 = f2 = 1
  print("Test case 8: f1 = f2 = 1")
  flag <- test_extHaus(cases1[[i]], cases2[[i]], 1000000, f1 = 1)
  flag_main <- max(flag_main, flag)
  
  if (flag_main == 0) {
    print("All test cases passed!")
  }
  print("--------------------------------------------------------------------")
}
# Test observations:
# All input sets passed all test cases!

# Testing: SpatialDataFrames ----------------------------------------------

setwd("Data")

# Read shp file
tracts_houston <- rgdal::readOGR("Census_2010_Tracts.shp")

# Specifying projection information for SpatialPolygonsDataFrame object
tracts_houston <- sp::spTransform(
  tracts_houston,
  sp::CRS(sp::proj4string(tracts_houston))
)

# Extract tracts of Harris County
tracts_harris <- tracts_houston[grepl(c("201"), tracts_houston@data$COUNTY), ]
tracts_harris@data <- tracts_harris@data[, c("OBJECTID", "STATE",
                                             "COUNTY","TRACT", "SUM_TotPop")]
tracts_harris <- sp::spTransform(tracts_harris, CRS("+init=epsg:2278"))

# Extract river data
rivers <- rgdal::readOGR("Major_Rivers.gdb")
rivers <- sp::spTransform(rivers, CRS("+init=epsg:2278"))

# Generate random values of f1 and f2
set.seed(123)
f1_vec <- rbeta(length(tracts_harris)^2, shape1 = 0.5, shape2 = 0.5)
f2_vec <- rbeta(length(tracts_harris)^2, shape1 = 0.5, shape2 = 0.5)
# using the specified shape parameters we should be able to sample
# values of f1 and f2 close to 0 and 1 to verify the function's correctness
# near the edge cases.

# Generate a SpatialPointsDataFrame using random coordinates
n_points <- 100
coords <- matrix(c(1, 1), nrow = 1, ncol = 2, byrow = T)
data <- as.data.frame(matrix(rep(0, n_points * 2), nrow = n_points, ncol = 2))
for (i in 2:n_points) {
  coords <- rbind(coords, matrix(runif(2), ncol = 2, nrow = 1, byrow = T))
}
coords <- as.data.frame(coords)
crdref <- sp::CRS("+init=epsg:2278")
spdf <- SpatialPointsDataFrame(
  coords = coords,
  data = data,
  proj4string = crdref
)

test_extHausv2 <- function(A, B, a_coords, n_samp, f1, f2) {
  flag <- 0
  
  val <- extHaus(A, B, f1, f2)
  if (is.null(gDifference(A, B))) {
    ans <- 0
  } else {
    distsf1 <- rgeos::gDistance(a_coords, B, byid = T)
    b_coords <- sp::spsample(B, n = n_samp, type = "regular")
    distsf2 <- rgeos::gDistance(b_coords, A, byid = T)
    ans <- max(quantile(distsf1, probs = f1),
               quantile(distsf2, probs = f2)) 
  }
  error <- abs(val - ans)
  pct_error <- 100 * error / ans
  if (error > 0.01 * ans) {
    flag <- 1
    print("Error")
    print(paste0("Expected ", ans, " but computed ", val))
    print(paste0("Percentage error: ", pct_error, "%"))
    print(paste0("f1_vec: ", f1))
    print(paste0("f2_vec: ", f2))
    print(paste0("Type of A: ", class(A)[1]))
    print(paste0("Type of B: ", class(B)[1]))
  }
  return(flag)
}

# Test polygon-to-(polygon / line / point) calculations
n_poly <- 10
n_other <- 10
n_samp <- 1000000
for (i in 1:n_poly) {
  print(paste0("Outer iteration: ", i))
  flag <- 0
  a_coords <- sp::spsample(tracts_harris[i, ], n = n_samp, type = "regular")
  for (j in 1:n_other) {
    print(paste0("Iteration (", i, ", ", j, ")"))
    f1val <- f1_vec[i * j]
    f2val <- f2_vec[i * j]
    
    flag1 <- test_extHausv2(
      A = tracts_harris[i, ],
      B = tracts_harris[j, ],
      a_coords = a_coords,
      n_samp = n_samp,
      f1 = f1val,
      f2 = f2val
    )
    print("Finished flag 1---------------")
    
    flag2 <- test_extHausv2(
      A = tracts_harris[i, ],
      B = rivers[j, ],
      a_coords = a_coords,
      n_samp = n_samp,
      f1 = f1val,
      f2 = f2val
    )
    print("Finished flag 2---------------")
    
    flag3 <- 0
    val <- extHaus(tracts_harris[i, ], spdf[j, ], f1val, f2val)
    dists <- rgeos::gDistance(a_coords, spdf[j, ], byid = T)
    min_dist <- rgeos::gDistance(spdf[j, ], tracts_harris[i, ])
    ans <- max(quantile(dists, f1val), min_dist)
    error <- abs(val - ans)
    pct_error <- 100 * error / ans
    if (error > 0.01 * ans) {
      flag3 <- 1
      print(paste0("Error"))
      print(paste0("Expected ", ans, " but computed ", val))
      print(paste0("Percentage error: ", pct_error, "%"))
      print(paste0("f1_vec: ", f1val))
      print(paste0("f2_vec: ", f2val))
    }
    print("Finished flag 3---------------")
    
    flag_mid <- max(flag1, flag2, flag3)
    if (flag_mid == 1) {
      flag <- flag_mid
    }
    print("------------------------------------------------------------------")
  }
  if (flag == 0) {
    print("Completed outer iteration without errors.")
  }
  print("********************************************************************")
}
# Old version: 
  # 13 / 100 iterations with more than 2% error
  # 9 / 13 of those had more than 30% error
  # Increasing the number of points sampled (by 10x) did not fix this
  # Halving the tolerance also did not affect the error 
# The percentage error did not even decrease a little bit in either case
# New version:
  # No iteration with more than 2% error
  # Iterations with more than 1% error: (3, 2), (3, 6), (4, 10), (5, 4), (8, 2), (8, 10), (10, 4)

# Test line-to-(polygon / line / point) calculations
n_line <- 10
n_other <- 10
n_samp <- 1000000
for (i in 1:n_line) {
  print(paste0("Outer iteration: ", i))
  flag <- 0
  a_coords <- sp::spsample(rivers[i, ], n = n_samp, type = "regular")
  for (j in 1:n_other) {
    print(paste0("Iteration (", i, ", ", j, ")"))
    f1val <- f1_vec[i * j]
    f2val <- f2_vec[i * j]
    
    flag1 <- test_extHausv2(
      A = rivers[i, ],
      B = tracts_harris[j, ],
      a_coords = a_coords,
      n_samp = n_samp,
      f1 = f1val,
      f2 = f2val
    )
    print("Finished flag 1---------------")
    
    flag2 <- test_extHausv2(
      A = rivers[i, ],
      B = rivers[j, ],
      a_coords = a_coords,
      n_samp = n_samp,
      f1 = f1val,
      f2 = f2val
    )
    print("Finished flag 2---------------")
    
    flag3 <- 0
    val <- extHaus(rivers[i, ], spdf[j, ], f1val, f2val)
    dists <- rgeos::gDistance(a_coords, spdf[j, ], byid = T)
    min_dist <- rgeos::gDistance(spdf[j, ], rivers[i, ])
    ans <- max(quantile(dists, f1val), min_dist)
    error <- abs(val - ans)
    pct_error <- 100 * error / ans
    if (error > 0.01 * ans) {
      flag3 <- 1
      print(paste0("Error"))
      print(paste0("Expected ", ans, " but computed ", val))
      print(paste0("Percentage error: ", pct_error, "%"))
      print(paste0("f1_vec: ", f1val))
      print(paste0("f2_vec: ", f2val))
    }
    print("Finished flag 3---------------")
    
    flag_mid <- max(flag1, flag2, flag3)
    if (flag_mid == 1) {
      flag <- flag_mid
    }
    print("------------------------------------------------------------------")
  }
  if (flag == 0) {
    print("Completed outer iteration without errors.")
  }
  print("********************************************************************")
}
# Old version: a lot of iterations with big percentage errors (> 25%)
# New version: no iterations with more than 2% error


# Test point to point calculations
n_point <- 10
n_other <- 10
for (i in 1:n_point) {
  print(paste0("Outer iteration: ", i))
  flag <- 0
  for (j in 1:n_other) {
    print(paste0("Iteration (", i, ", ", j, ")"))
    f1val <- f1_vec[i * j]
    f2val <- f2_vec[i * j]
    val <- extHaus(spdf[i, ], spdf[j, ], f1val, f2val)
    ans <- rgeos::gDistance(spdf[i, ], spdf[j, ])
    error <- abs(val - ans)
    pct_error <- 100 * error / ans
    if (error > 0.01 * ans) {
      flag <- 1
      print(paste0("Error"))
      print(paste0("Expected ", ans, " but computed ", val))
      print(paste0("Percentage error: ", pct_error, "%"))
      print(paste0("f1_vec: ", f1val))
      print(paste0("f2_vec: ", f2val))
    }
    print("------------------------------------------------------------------")
  }
  if (flag == 0) {
    print("Completed outer iteration without errors.")
  }
  print("********************************************************************")
}
# Test observations: no errors!

# Testing: hausMat --------------------------------------------------------

test_haus <- function(mat, mat_size, type, f1, f2 = f1) {
  if (type == "line") {
    inputs <- rivers[1:mat_size, ]
  } else if (type == "area") {
    inputs <- tracts_harris[1:mat_size, ]
  } else if (type == "point") {
    inputs <- spdf[1:mat_size, ]
  } else {
    stop("type must be either 'line', 'area' or 'point'.")
  }
  
  for (i in 1:mat_size) {
    print(paste0("Outer iteration: ", i))
    flag <- 0
    for (j in 1: mat_size) {
      # compute real value and absolute error
      val1 <- extHaus(inputs[i, ], inputs[j, ], f1 = f1, f2 = f2)
      error1 <- abs(mat[i, j] - val1)
      
      # compare the value in the matrix to the value computed by extHaus to
      # check for correctness
      if (error1 > 0.01 * val1) {
        flag <- 1
        print(paste0("Error in mat1 iteration (", i, ", ", j, ")"))
        print(paste0("Expected ", val1, " but computed ", error1))
        print(paste0("Percentage error: ", 100 * error1 / val1, "%"))
        print("----------------------------------------------------------------")
      }
    }
    if (flag == 0) {
      print("Completed outer iteration without errors.")
    }
    print("********************************************************************")
  }
}

# Test area-area calculations
mat_size <- 5
mat1 <- hausMat(
  tracts_harris[1:mat_size, ],
  f1 = 0.5,
  ncores = detectCores() - 1
)
mat2 <- hausMat(
  tracts_harris[1:mat_size, ],
  f1 = 0.75,
  f2 = 0.25,
  ncores = detectCores() - 1
)
test_haus(mat1, mat_size, "area", f1 = 0.5)
test_haus(mat2, mat_size, "area", f1 = 0.75, f2 = 0.25)
# Test observations: no errors

# Test line-to-line calculations
mat_size <- 5
mat1 <- hausMat(
  rivers[1:mat_size, ],
  f1 = 0.5,
  ncores = detectCores() - 1
)
mat2 <- hausMat(
  rivers[1:mat_size, ],
  f1 = 0.75,
  f2 = 0.25,
  ncores = detectCores() - 1
)
test_haus(mat1, mat_size, "line", f1 = 0.5)
test_haus(mat2, mat_size, "line", f1 = 0.75, f2 = 0.25)
# Test observations: no errors

# Test point-to-point calculations
mat_size <- 5
mat1 <- hausMat(
  spdf[1:mat_size, ],
  f1 = 0.5,
  ncores = detectCores() - 1
)
mat2 <- hausMat(
  spdf[1:mat_size, ],
  f1 = 0.75,
  f2 = 0.25,
  ncores = detectCores() - 1
)
test_haus(mat1, mat_size, "point", f1 = 0.5)
test_haus(mat2, mat_size, "point", f1 = 0.75, f2 = 0.25)
# Test observations: no errors

# Testing: hausMat file I/O -----------------------------------------------

# Case 1: f1 = f2
mat <- hausMat(
  spdf[1:5, ],
  f1 = 0.5,
  fileout = TRUE,
  filename = "results.csv",
  ncores = detectCores() - 1
)
file_data <- read.csv("results.csv")
print(file_data)
print(mat)

# Case 2: f1 does not equal f2
mat <- hausMat(
  spdf[1:5, ],
  f1 = 0.25,
  f2 = 0.75,
  fileout = TRUE,
  filename =  "results.csv",
  ncores = detectCores() - 1
)
file_data <- read.csv("results.csv")
print(file_data)
print(mat)
# Testing: par_haus_mat_slow vs par_haus_mat ------------------------------

# If par_haus_mat is better than par_haus_mat_slow then we would expect for
# the time difference between par_haus_mat_slow and par_haus_mat to become
# increasingly positive as the size of the input matrix grows. This is because
# the number of parallel tasks created by par_haus_mat_slow would grow whereas the
# number of parallel tasks created by par_haus_mat would remain constant
# and equal to ncores.

n <- 100
coords <- matrix(c(1, 1), nrow = 1, ncol = 2, byrow = T)
data <- as.data.frame(matrix(rep(0, n * 2), nrow = n, ncol = 2))
for (i in 2:n) {
  coords <- rbind(coords, matrix(runif(2), ncol = 2, nrow = 1, byrow = T))
}
coords <- as.data.frame(coords)
crdref <- sp::CRS("+init=epsg:2278")
spdf <- SpatialPointsDataFrame(
  coords = coords,
  data = data,
  proj4string = crdref
)

hm_run_time <- rep(0, n - 1)
fast_run_time <- rep(0, n - 1)
diff <- rep(0, n - 1)
f1vec <- seq(from = 0, to = 1, length.out = n) 

#changing input size f1 does not equal f2
maxN <- 100
num_obs <- 10
for (size in 77:maxN) {
  print(size)
  for (i in 1:num_obs) {
    print(i)
    hm_run_time[size - 1] <- hm_run_time[size - 1] + par_haus_mat_slow(
      spdf[1:size, ], 
      f1 = 0.25,
      f2 = 0.75,
      ncores = detectCores() - 1,
      tol = 0.01
    )
    fast_run_time[size - 1] <- fast_run_time[size - 1] + par_haus_mat(
      spdf[1:size, ],
      f1 = 0.25,
      f2 = 0.75,
      ncores = detectCores() - 1,
      tol = 0.01
    )
  }
  hm_run_time[size - 1] <- hm_run_time[size - 1] / num_obs
  fast_run_time[size - 1] <- fast_run_time[size - 1] / num_obs
  diff[size - 1] <- hm_run_time[size - 1] - fast_run_time[size - 1]
  print("--------------------------------------------------------------------")
}
plot(hm_run_time[1:maxN - 1])
plot(fast_run_time[1:maxN - 1])
plot(diff[1:maxN - 1])
min(diff)
mean(diff)
max(diff)
hist(hm_run_time[1:maxN - 1])
hist(fast_run_time[1:maxN - 1])
hist(diff[1:maxN - 1])
x <- 2:maxN
summary(lm(diff[1:maxN-1] ~ x))
# is there a statistically significant slope?

# using x = 2:81 and the first 80 entries of diff we get:
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept) -0.171584   0.275523  -0.623 0.535262    
#x            0.021226   0.005801   3.659 0.000459 ***
# slope is statistically significant

#changing input size f1 = f2
maxN <- 100
num_obs <- 10
for (size in 2:maxN) {
  print(size)
  for (i in 1:num_obs) {
    print(i)
    hm_run_time[size - 1] <- hm_run_time[size - 1] + par_haus_mat_slow(
      spdf[1:size, ],
      f1 = 0.5,
      ncores = detectCores() - 1,
      tol = 0.01
    )
    fast_run_time[size - 1] <- fast_run_time[size - 1] + par_haus_mat(
      spdf[1:size, ],
      f1 = 0.5,
      ncores = detectCores() - 1,
      tol = 0.01
    )
  }
  hm_run_time[size - 1] <- hm_run_time[size - 1] / num_obs
  fast_run_time[size - 1] <- fast_run_time[size - 1] / num_obs
  diff[size - 1] <- hm_run_time[size - 1] - fast_run_time[size - 1]
  print("--------------------------------------------------------------------")
}
plot(hm_run_time[1:maxN - 1])
plot(fast_run_time[1:maxN - 1])
plot(diff[1:maxN - 1])
min(diff) # -2.72
mean(diff) # 1.02
max(diff) # 5.81
hist(hm_run_time[1:maxN - 1])
hist(fast_run_time[1:maxN - 1])
hist(diff[1:maxN - 1])
x <- 2:maxN
summary(lm(diff[1:maxN-1] ~ x))

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept) -0.614333   0.174692  -3.517 0.000666 ***
#x            0.032040   0.002988  10.722  < 2e-16 ***

num_obs <- 40
hm_run_time <- rep(0, num_obs)
fast_run_time <- rep(0, num_obs)
diff <- rep(0, num_obs)
for (i in 1:num_obs) {
  print(i)
  hm_run_time[i] <- par_haus_mat_slow(
    spdf[1:99, ], 
    f1 = 0.5,
    ncores = detectCores() - 1,
    tol = 0.01
  )
  fast_run_time[i] <- par_haus_mat(
    spdf[1:99, ],
    f1 = 0.5,
    ncores = detectCores() - 1,
    tol = 0.01
  )
  diff[i] = hm_run_time[i] - fast_run_time[i]
}
# histograms look skewed
hist(hm_run_time)
hist(fast_run_time)
mean(hm_run_time) # 27.25494
mean(fast_run_time) #25.15338
mean(diff) #2.101565
min(diff) # -11.13924
max(diff) #16.06571

# log transform looks more normally distributed
hist(log(hm_run_time))
hist(log(fast_run_time))

# Hypothesis testing
# H0: mean(log_hm_run_time) <= mean(log_fast_run_time)
# H1: mean(log_hm_run_time) > mean(log_fast_run_time)
log_hm_run_time <- log(hm_run_time)
log_fast_run_time <- log(fast_run_time)
log_hm_run_time <- log_hm_run_time[is.finite(log_hm_run_time)]
log_fast_run_time <- log_fast_run_time[is.finite(log_fast_run_time)]
s1 <- sd(log_hm_run_time)
s2 <- sd(log_fast_run_time)
s1_adj <- s1^2 / length(log_hm_run_time)
s2_adj <- s2^2 / length(log_fast_run_time)
spool <- sqrt(s1_adj + s2_adj)
test_statistic <- (mean(log_hm_run_time) - mean(log_fast_run_time)) / spool
pval <- 1 - pnorm(test_statistic)
#pval: 0.002393715
# Reject null hypothesis

#changing ncores

# We would expect the performance gap between par_haus_mat and 
# par_haus_mat_slow to decrease as the number of cores increases. This is 
# because as more cores are used, the number of parallel tasks per core created
# by par_haus_mat_slow will decrease.

n <- 100
maxcores <- detectCores() - 1
hm_run_time <- matrix(rep(0, maxcores * n), nrow = maxcores)
fast_run_time <- matrix(rep(0, maxcores * n), nrow = maxcores)
diff <- matrix(rep(0, maxcores * n), nrow = maxcores)
for (core in 1:maxcores) {
  for (i in 1:n) {
    hm_run_time[core, i] <- par_haus_mat_slow(
      spdf[1:50, ],
      f1 = 1,
      ncores = core,
      tol = 0.01
    )
    fast_run_time[core, i] <- par_haus_mat(
      spdf[1:50, ],
      f1 = 1,
      ncores = core,
      tol = 0.01
    )
    diff[core, i] <- hm_run_Time[core, i] - fast_run_time[core, i]
  }
}
hm_means <- apply(hm_run_time, 1, mean)
fast_means <- apply(fast_run_time, 1, mean)
diff_means <- apply(diff, 1, mean)
plot(diff_means)

# Testing: gDistance(hausdorff=T) vs. directHaus --------------------------

# Initialize test parameters
n <- length(tracts_harris) + length(rivers) # number of iterations
test_tol <- NULL # tolerance
npts <- 100000 
# number of points that will be sampled from tracts/rivers to compute the 
# "real" distance
a_coords <- spsample(tracts_harris[1, ], n = npts, type = "regular")

# Initialize data collection vectors
time_gDist <- rep(0, n)
time_dirHaus <- rep(0, n)
time_diff <- rep(0, n)

for (i in 1:n) {
  flag <- 0 # value is changed to 1 if an error or strange event occurs
  if (i <= length(tracts_harris)) {
    # Compute distance using gDistance
    start1 <- Sys.time()
    val1 <- rgeos::gDistance(
      tracts_harris[1, ],
      tracts_harris[i, ],
      hausdorff = T
    )
    time_gDist[i] <- Sys.time() - start1
    
    # Compute distance using directHaus
    start2 <- Sys.time()
    A_to_B <- directHaus(
      tracts_harris[1, ],
      tracts_harris[i, ],
      f1 = 1,
      tol = tol
    )
    B_to_A <- directHaus(
      tracts_harris[i, ],
      tracts_harris[1, ],
      f1 = 1,
      tol = tol
    )
    val2 <- max(A_to_B, B_to_A)
    time_dirHaus[i] <- Sys.time() - start2
    time_diff[i] <- time_dirHaus[i] - time_gDist[i]
    
    # Compute distance using sampling
    distsf1 <- rgeos::gDistance(a_coords, tracts_harris[i, ], byid = T)
    b_coords <- sp::spsample(tracts_harris[i, ], n = npts, type = "regular")
    distsf2 <- rgeos::gDistance(b_coords, tracts_harris[1, ], byid = T)
    val3 <- max(max(distsf1), max(distsf2))
  } else {
    # Compute distance using gDistance
    start1 <- Sys.time()
    val1 <- rgeos::gDistance(
      tracts_harris[1, ],
      rivers[i - length(tracts_harris), ],
      hausdorff = T
    )
    time_gDist[i] <- Sys.time() - start1
    
    # Compute distance using directHaus
    start2 <- Sys.time()
    A_to_B <- directHaus(
      tracts_harris[1, ],
      rivers[i - length(tracts_harris), ],
      f1 = 1,
      tol = tol
    )
    B_to_A <- directHaus(
      rivers[i - length(tracts_harris), ],
      tracts_harris[1, ],
      f1 = 1,
      tol = tol
    )
    val2 <- max(A_to_B, B_to_A)
    time_dirHaus[i] <- Sys.time() - start2
    time_diff[i] <- time_dirHaus[i] - time_gDist[i]
    
    # Compute distance using sampling
    distsf1 <- rgeos::gDistance(
      a.coords,
      rivers[i - length(tracts_harris), ],
      byid = T
    )
    b_coords <- sp::spsample(
      rivers[i - length(tracts_harris), ],
      n = npts,
      type = "regular"
    )
    distsf2 <- rgeos::gDistance(b_coords, tracts_harris[1, ], byid = T)
    val3 <- max(max(distsf1), max(distsf2))
  }
  # Compute error using the distance computed using the sampling method as
  # the "true" value
  abs_error_gDistance <- abs(val3 - val1)
  pct_error_gDistance <- 100 * abs_error_gDistance / val3
  abs_error_direct <- abs(val3 - val2)
  pct_error_direct <- 100 * abs_error_direct / val3
  
  # Check for significant errors or stange events
  if (abs_error_direct > 0 && pct_error_direct > 2) {
    print(paste0("Error at iteration ", i))
    print(paste0("sampling method: ", val3))
    print(paste0("Absolute gDistance error: ", abs_error_gDistance))
    print(paste0("Percentage gDistance error: ", pct_error_gDistance, "%"))
    print(paste0("gDistance: ", val1))
    print(paste0("Absolute directHaus error: ", abs_error_direct))
    print(paste0("Percentage directHaus error: ", pct_error_direct, "%"))
    print(paste0("directHaus: ", val2))
    flag <- 1
  }
  if (val2 > val1) {
    print(paste0("Iteration: ", i))
    print("directHaus distance > gDistance distance")
    print(paste0("directHaus: ", val2))
    print(paste0("gDistance: ", val1))
    flag <- 1
  }
  if (val2 > val3) {
    print(paste0("Iteration: ", i))
    print("directHaus distance > sampling distance")
    print(paste0("sampling method: ", val3))
    print(paste0("directHaus: ", val2))
    flag <- 1
  }
  if (val3 > val1) {
    print(paste0("Iteration: ", i))
    print("sampling distance > gDistance distance")
    print(paste0("sampling method: ", val3))
    print(paste0("gDistance: ", val1))
    flag <- 1
  }
  if (abs_error_gDistance > abs_error_direct) {
    print(paste0("Iteration: ", i))
    print("gDistance has greater error than directHaus")
    print(paste0("sampling method: ", val3))
    print(paste0("Absolute gDistance error: ", abs_error_gDistance))
    print(paste0("Percentage gDistance error: ", pct_error_gDistance, "%"))
    print(paste0("gDistance: ", val1))
    print(paste0("Absolute directHaus error: ", abs_error_direct))
    print(paste0("Percentage directHaus error: ", pct_error_direct, "%"))
    print(paste0("directHaus: ", val2))
    flag <- 1
  }
  if (flag == 1) {
    print("------------------------------------------------------------------")
  }
}
min(time_diff) # 0.00509
mean(time_diff) # 0.7610596
max(time_diff) # 10.00405
hist(time_diff)

# histograms look skewed
hist(time_dirHaus)
hist(time_gDist)

# log transform looks more normally distributed
hist(log(time_dirHaus))
hist(log(time_gDist))

log_time_dirHaus <- log(time_dirHaus)
log_time_gDist <- log(time_gDist)
log_time_dirHaus <- log_time_dirHaus[is.finite(log_time_dirHaus)]
log_time_gDist <- log_time_gDist[is.finite(log_time_gDist)]

# Do hypothesis testing
# H0: mean(log_time_dirHaus) <= mean(log_time_gDist)
# H1: mean(log_time_dirHaus) > mean(log_time_gDist)
s1 <- sd(log_time_dirHaus)
s2 <- sd(log_time_gDist)
s1_adj <- s1^2 / length(log_time_dirHaus)
s2_adj <- s2^2 / length(log_time_gDist)
spool <- sqrt(s1_adj + s2_adj)
test_statistic <- (mean(log_time_dirHaus) - mean(log_time_gDist)) / spool
pval <- 1 - pnorm(test_statistic)
# pval: 0 
# statistically significant: reject null hypothesis

# Test observations using NULL tolerance:
# Distance computed by directHaus is smaller than distance computed by gDistance
# Most of the calculations between census tracts have less than 2% error
# Only 3 cases (iter: 199, 201, 748) have greater than 2% error
# These cases are still below 5% error
# A lot of calculations between census tract 1 and rivers have more than 2% error
# Distance computed by sampling method is generally (but not always) smaller than the distance computed by gDistance
# Distance computed by sampling method is generally (but not always) larger than the distance computed by directHaus
# The distance computed by the sampling method is generally (but not always) closer to the distance computed by gDistance than it is to the distance computed by directHaus
# directHaus performed poorly on the tracts.harris[1,]-to-rivers[i - length(tracts.harris),] trials
# ~20 iterations with more than 2% error
# ~15 iterations with more than 20% error
# gDistance did not perform poorly on the tracts-to-rivers trials (no iteration with more than 2% error)
# should maybe increase the number of points sampled in directHaus to get better results
