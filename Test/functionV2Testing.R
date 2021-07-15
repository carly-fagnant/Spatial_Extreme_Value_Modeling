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

# Testing Hausdorff: Points, Lines, and Polygons --------------------------

cases1 <- c(pt1, pt1, pt1, l1, l1, p1)
cases2 <- c(pt2, l1, p1, l2, p1, p2)

# Compute preliminary results: sample points from each of the input areas and
# get the distance quantiles. Use these to estimate the directed 
# Hausdorff distance between the input areas
f1_mat <- matrix(rep(-1, 6 * 4), nrow = 6)
f2_mat <- matrix(rep(-1, 6 * 4), nrow = 6)
for (input_case in 1:6) {
  input1 <- cases1[[input_case]]
  input2 <- cases2[[input_case]]
  n <- 100000
  a_coords <- sp::spsample(input1, n = n, type = "regular")
  distsf1 <- rgeos::gDistance(a_coords, input2, byid = T)
  b_coords <- sp::spsample(input2, n = n, type = "regular")
  distsf2 <- rgeos::gDistance(b_coords, input1, byid = T)
  f1_mat[input_case,] <- quantile(distsf1, probs = c(0.25, 0.5, 0.75, 1))
  f2_mat[input_case,] <- quantile(distsf2, probs = c(0.25, 0.5, 0.75, 1))
}

# Compute final answers by getting results from f1_mat and f2_mat
answers <- matrix(rep(-1, 6 * 11), nrow = 6)
for (input_case in 1:6) {
  input1 <- cases1[[input_case]]
  input2 <- cases2[[input_case]]
  f1.25 <- f1_mat[input_case, 1]
  f1.50 <- f1_mat[input_case, 2]
  f1.75 <- f1_mat[input_case, 3]
  f1.100 <- f1_mat[input_case, 4]
  f2.25 <- f2_mat[input_case, 1]
  f2.50 <- f2_mat[input_case, 2]
  f2.75 <- f2_mat[input_case, 3]
  f2.100 <- f2_mat[input_case, 4]
  for (test_case in 1:11) {
    if (test_case == 1) {
      val <- max(rgeos::gDistance(input1, input2), f2.50)
    } else if (test_case == 2) {
      val <- max(f1.50, rgeos::gDistance(input1, input2))
    } else if (test_case == 3) {
      val <- gDistance(input1, input2)
    } else if (test_case == 4) {
      val <- max(f1.25, f2.75)
    } else if (test_case == 5) {
      val <- max(f1.75, f2.25)
    } else if (test_case == 6) {
      val <- max(f1.25, f2.25)
    } else if (test_case == 7) {
      val <- max(f1.50, f2.50)
    } else if (test_case == 8) {
      val <- max(f1.75, f2.75)
    } else if (test_case == 9) {
      val <- max(f1.100, f2.50)
    } else if (test_case == 10) {
      val <- max(f1.50, f2.100)
    } else {
      val <- rgeos::gDistance(input1, input2, hausdorff = T)
    }
    answers[input_case, test_case] <- val
  }
}
answers[3, ] <- rep(0, 11)

# Iterate through the input cases and see if the computed values match the
# estimated final answers
for (i in 1:6) {
  flag <- 0
  print(paste0("Testing input set ", i))
  print(paste0("Class of A: ",  class(cases1[[i]])[1]))
  print(paste0("Class of B: ",  class(cases2[[i]])[1]))
  
  ## Case 1: f1 = 0
  test1.1 <- extHaus(cases1[[i]], cases2[[i]], f1 = 0, f2 = 0.5)
  if (abs(test1.1 - answers[i, 1]) > 0.01) {
    flag <- 1
    print(paste0("Test case 1 failed on input pair set: ", i))
    print(paste0("Expected ", answers[i, 1], " but computed ", test1.1))
  }
  
  ## Case 2: f2 = 0
  test2.1 <- extHaus(cases1[[i]], cases2[[i]], f1 = 0.5, f2 = 0)
  if (abs(test2.1 - answers[i, 2]) > 0.01) {
    flag <- 1
    print(paste0("Test case 2 failed on input pair set: ", i))
    print(paste0("Expected ", answers[i, 2], " but computed ", test2.1))
  }
  
  ## Case 3: f1 = f2 = 0
  test3.1 <- extHaus(cases1[[i]], cases2[[i]], f1 = 0)
  if (abs(test3.1 - answers[i, 3]) > 0.01) {
    flag <- 1
    print(paste0("Test case 3 failed on input pair set: ", i))
    print(paste0("Expected ", answers[i, 3], " but computed ", test3.1))
  }
  
  ## Case 4: 0 < f1, f2 < 1
  test4.1 <- extHaus(cases1[[i]], cases2[[i]], f1 = 0.25, f2 = 0.75)
  test4.2 <- extHaus(cases1[[i]], cases2[[i]], f1 = 0.75, f2 = 0.25)
  if (abs(test4.1 - answers[i, 4]) > 0.01) {
    flag <- 1
    print(paste0("Test case 4 failed on input pair set: ", i))
    print(paste0("Expected ", answers[i, 4], " but computed ", test4.1))
  }
  if (abs(test4.2 - answers[i, 5]) > 0.01) {
    flag <- 1
    print(paste0("Test case 4 failed on input pair set: ", i))
    print(paste0("Expected ", answers[i, 5], " but computed ", test4.2))
  }
  
  ## Case 5: 0 < (f1 = f2) < 1
  test5.1 <- extHaus(cases1[[i]], cases2[[i]], f1 = 0.25)
  test5.2 <- extHaus(cases1[[i]], cases2[[i]], f1 = 0.5)
  test5.3 <- extHaus(cases1[[i]], cases2[[i]], f1 = 0.75)
  if (abs(test5.1 - answers[i, 6]) > 0.01) {
    flag <- 1
    print(paste0("Test case 5 failed on input pair set: ", i))
    print(paste0("Expected ", answers[i, 6], " but computed ", test5.1))
  }
  if (abs(test5.2 - answers[i, 7]) > 0.01) {
    flag <- 1
    print(paste0("Test case 5 failed on input pair set: ", i))
    print(paste0("Expected ", answers[i, 7], " but computed ", test5.2))
  }
  if (abs(test5.3 - answers[i, 8]) > 0.01) {
    flag <- 1
    print(paste0("Test case 5 failed on input pair set: ", i))
    print(paste0("Expected ", answers[i, 8], " but computed ", test5.3))
  }
  
  ## Case 6: f1 = 1
  test6.1 <- extHaus(cases1[[i]], cases2[[i]], f1 = 1, f2 = 0.5)
  if (abs(test6.1 - answers[i, 9]) > 0.01) {
    flag <- 1
    print(paste0("Test case 6 failed on input pair set: ", i))
    print(paste0("Expected ", answers[i, 9], " but computed ", test6.1))
  }
  
  ## Case 7: f2 = 1
  test7.1 <- extHaus(cases1[[i]], cases2[[i]], f1 = 0.5, f2 = 1)
  if (abs(test7.1 - answers[i, 10]) > 0.01) {
    flag <- 1
    print(paste0("Test case 6 failed on input pair set: ", i))
    print(paste0("Expected ", answers[i, 10], " but computed ", test7.1))
  }
  
  ## Case 8: f1 = f2 = 1
  test8.1 <- extHaus(cases1[[i]], cases2[[i]], f1 = 1)
  if (abs(test8.1 - answers[i, 11]) > 0.01) {
    flag <- 1
    print(paste0("Test case 6 failed on input pair set: ", i))
    print(paste0("Expected ", answers[i, 11], " but computed ", test8.1))
  }
  
  if (flag == 0) {
    print(paste0("All tests for input set ", i, " passed!"))
  }
  
  print("--------------------------------------------------------------------")
}

# Testing Hausdorff: SpatialDataFrames -----------------------------
setwd("Data")
## Read shp file
tracts_houston <- rgdal::readOGR("Census_2010_Tracts.shp")
## specifying projection information for SpatialPolygonsDataFrame object
tracts_houston <- sp::spTransform(
  tracts_houston,
  sp::CRS(sp::proj4string(tracts_houston))
)
## extract tracts of Harris County
tracts_harris <- tracts_houston[grepl(c("201"), tracts_houston@data$COUNTY), ]
tracts_harris@data <- tracts_harris@data[, c("OBJECTID", "STATE",
                                              "COUNTY","TRACT", "SUM_TotPop")]
tracts_harris <- sp::spTransform(tracts_harris, CRS("+init=epsg:2278"))

rivers <- rgdal::readOGR("Major_Rivers.gdb")
rivers <- sp::spTransform(rivers, CRS("+init=epsg:2278"))

set.seed(123)
f1_vec <- rbeta(length(tracts_harris)^2, shape1 = 0.5, shape2 = 0.5)
f2_vec <- rbeta(length(tracts_harris)^2, shape1 = 0.5, shape2 = 0.5)
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

# Polygon to Polygon / Line / Point
nPoly <- 10
nOther <- 10
for (i in 1:nPoly) {
  for (j in 1:nOther) {
    f1val <- f1_vec[i * j]
    f2val <- f2_vec[i * j]
    print(paste0("i: ", i))
    print(paste0("j: ", j))
    print(paste0("f1_vec: ", f1val))
    print(paste0("f2_vec: ", f2val))
    # compute distances using extHaus
    val1 <- extHaus(tracts_harris[i, ], tracts_harris[j, ], f1val, f2val)
    val2 <- extHaus(tracts_harris[i, ], rivers[j, ], f1val, f2val)
    val3 <- extHaus(tracts_harris[i, ], spdf[j, ], f1val, f2val)
    
    # compute "true" values by sampling points from the input areas and using
    # the corresponding distance quantiles as an estimate for the directed
    # Hausdorff distances between the two objects
    n <- 100000
    a_coords <- sp::spsample(tracts_harris[i, ], n = n, type = "regular")
    distsf1.1 <- rgeos::gDistance(a_coords, tracts_harris[j, ], byid = T)
    b_coords <- sp::spsample(tracts_harris[j, ], n = n, type = "regular")
    distsf1.2 <- rgeos::gDistance(b_coords, tracts_harris[i, ], byid = T)
    # Extended Hausdorff distance is the max of the two directed Hausdorff 
    # distance balues compute above
    ans1 <- max(quantile(distsf1.1, probs = f1val),
                quantile(distsf1.2, probs = f2val))
    
    if (abs(val1 - ans1) > 0.01) {
      print(paste0("Error at iteration (", i, ", ", j, ")"))
      print("Area to area calculation error")
      print(paste0("Expected ", ans1, " but computed ", val1))
    }
    
    distsf2.1 <- rgeos::gDistance(a_coords, rivers[j, ], byid = T)
    b_coords <- sp::spsample(rivers[j, ], n = n, type = "regular")
    distsf2.2 <- rgeos::gDistance(b_coords, tracts_harris[i, ], byid = T)
    ans2 <- max(quantile(distsf2.1, probs = f1val),
                quantile(distsf2.2, probs = f2val))
    
    if (abs(val2 - ans2) > 0.01) {
      print(paste0("Error at iteration (", i, ", ", j, ")"))
      print("Area to line calculation error")
      print(paste0("Expected ", ans2, " but computed ", val2))
    }
    
    distsf3.1 <- rgeos::gDistance(a_coords, spdf[j, ], byid = T)
    distsf3.2 <- rgeos::gDistance(spdf[j, ], tracts_harris[i, ])
    ans3 <- max(quantile(distsf3.1, probs = f1val), distsf3.2)
    
    if (abs(val3 - ans3) > 0.01) {
      print(paste0("Error at iteration (", i, ", ", j, ")"))
      print("Area to point calculation error")
      print(paste0("Expected ", ans3, " but computed ", val3))
    }
    print("------------------------------------------------------------------")
  }
  print("********************************************************************")
}

# Line to Polygon / Line / Point
nLine <- 10
nOther <- 10
for (i in 1:nLine) {
  for (j in 1:nOther) {
    f1val <- f1_vec[i * j]
    f2val <- f2_vec[i * j]
    print(paste0("i: ", i))
    print(paste0("j: ", j))
    print(paste0("f1_vec: ", f1val))
    print(paste0("f2_vec: ", f2val))
    val1 <- extHaus(rivers[i, ], tracts_harris[j, ], f1val, f2val)
    val2 <- extHaus(rivers[i, ], rivers[j, ], f1val, f2val)
    val3 <- extHaus(rivers[i, ], spdf[j, ], f1val, f2val)
    
    n <- 100000
    a_coords <- sp::spsample(rivers[i, ], n = n, type = "regular")
    
    distsf1.1 <- rgeos::gDistance(a_coords, tracts_harris[j, ], byid = T)
    b_coords <- sp::spsample(tracts_harris[j, ], n = n, type = "regular")
    distsf1.2 <- rgeos::gDistance(b_coords, tracts_harris[i, ], byid = T)
    ans1 <- max(quantile(distsf1.1, probs = f1val),
                quantile(distsf1.2, probs = f2val))
    
    if (abs(val1 - ans1) > 0.01) {
      print(paste0("Error at iteration (", i, ", ", j, ")"))
      print("Area to area calculation error")
      print(paste0("Expected ", ans1, " but computed ", val1))
    }
    
    distsf2.1 <- rgeos::gDistance(a_coords, rivers[j, ], byid = T)
    b_coords <- sp::spsample(rivers[j, ], n = n, type = "regular")
    distsf2.2 <- rgeos::gDistance(b_coords, tracts_harris[i, ], byid = T)
    ans2 <- max(quantile(distsf2.1, probs = f1val),
                quantile(distsf2.2, probs = f2val))
    
    if (abs(val2 - ans2) > 0.01) {
      print(paste0("Error at iteration (", i, ", ", j, ")"))
      print("Area to line calculation error")
      print(paste0("Expected ", ans2, " but computed ", val2))
    }
    
    distsf3.1 <- rgeos::gDistance(a_coords, spdf[j, ], byid = T)
    distsf3.2 <- rgeos::gDistance(spdf[j, ], tracts_harris[i, ])
    ans3 <- max(quantile(distsf3.1, probs = f1val), distsf3.2)
    
    if (abs(val3 - ans3) > 0.01) {
      print(paste0("Error at iteration (", i, ", ", j, ")"))
      print("Area to point calculation error")
      print(paste0("Expected ", ans3, " but computed ", val3))
    }
    print("------------------------------------------------------------------")
  }
  print("********************************************************************")
}

# Point to Point
nPoint <- 10
nOther <- 10
for (i in 1:nPoint) {
  for (j in 1:nOther) {
    f1val <- f1_vec[i * j]
    f2val <- f2_vec[i * j]
    print(paste0("i: ", i))
    print(paste0("j: ", j))
    print(paste0("f1_vec: ", f1val))
    print(paste0("f2_vec: ", f2val))
    val <- extHaus(spdf[i, ], spdf[j, ], f1val, f2val)
    error <- abs(val - rgeos::gDistance(spdf[i, ], spdf[j, ]))
    if ( error > 0.01 ) {
      print(paste0("Error in iteration (", i, ", ", j, ")"))
    }
    print("------------------------------------------------------------------")
  }
  print("********************************************************************")
}

# Testing Hausdorff: hausMat -----------------------------------------------

# line-line testing
mat1 <- hausMat(
  rivers[1:5, ],
  f1 = 0.5,
  ncores = detectCores() - 1,
  tol = 0.001
)
mat2 <- hausMat(
  rivers[1:5, ],
  f1 = 0.75,
  f2 = 0.25,
  ncores = detectCores() - 1,
  tol = 0.001
)
for (i in 1:5) {
  for (j in 1:5) {
    # test case 1: f1 = f2
    val1 <- extHaus(rivers[i, ], rivers[j, ], 0.5, tol = 0.001)
    quant1 <- abs(mat1[i, j] - val1)
    # compare the value in the matrix to the value computed by extHaus to
    # check for correctness
    if (quant1 > 0.001 * val1) {
      print(paste0("Error in mat1 iteration (", i, ", ", j, ")"))
      print(paste0("Absolute error: ", quant1))
      print(paste0("Fractional error: ", quant1 / val1))
    }
    
    # test case 2: f1 does not equal f2
    val2 <- extHaus(rivers[i, ], rivers[j, ], 0.75, 0.25, tol = 0.001)
    quant2 <- abs(mat2[i, j] - val2)
    if (quant2 > 0.001 * val2) {
      print(paste0("Error in mat2 iteration (", i, ", ", j, ")"))
      print(paste0("Absolute error: ", quant2))
      print(paste0("Fractional error: ", quant2 / val2))
    }
    print("---------------------------------------------------")
  }
}

# area-area testing
mat1 <- hausMat(
  tracts_harris[1:10, ],
  f1 = 0.5,
  ncores = detectCores() - 1,
  tol = 0.001
)
mat2 <- hausMat(
  tracts_harris[1:10, ],
  f1 = 0.75,
  f2=0.25,
  ncores = detectCores() - 1,
  tol = 0.001
)
for (i in 1:10) {
  for (j in 1:10) {
    
    # test case 1: f1 = f2
    val1 <- extHaus(tracts_harris[i, ], tracts_harris[j, ], 0.5, tol = 0.001)
    quant1 <- abs(mat1[i, j] - val1)
    if (quant1 > 0.001 * val1) {
      print(paste0("Error in mat1 iteration (", i, ", ", j, ")"))
      print(paste0("Absolute error: ", quant1))
      print(paste0("Fractional error: ", quant1 / val1))
    }
    
    # test case 2: f1 does not equal f2
    val2 <- extHaus(
      tracts_harris[i, ],
      tracts_harris[j, ],
      0.75,
      0.25,
      tol = 0.001
    )
    quant2 <- abs(mat2[i, j] - val2)
    if (quant2 > 0.001 * val2) {
      print(paste0("Error in mat2 iteration (", i, ", ", j, ")"))
      print(paste0("Absolute error: ", quant2))
      print(paste0("Fractional error: ", quant2 / val2))
    }
    print("---------------------------------------------------")
  }
}

# Testing Hausdorff: par_haus_mat vs par_haus_mat_fast ----------------------------

# If parHausMatFastBoi is better than parHausMat then we would expect for
# the time difference between parHausMat and parHausMatFastBOi to become
# increasingly positive as the size of the input matrix grows. This is because
# the number of parallel tasks created by parHausMat would grow whereas the
# number of parallel tasks created by parHausMatFastBoi would remain constant
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
for (size in 2:maxN) {
  for (i in 1:num_obs) {
    hm_run_time[size - 1] <- hm_run_time[size - 1] + par_haus_mat(
      spdf[1:size,], 
      f1 = 0.25,
      f2 = 0.75,
      ncores = detectCores() - 1,
      tol = 0.01
    )
    fast_run_time[size - 1] <- fast_run_time[size - 1] + par_haus_mat_fast(
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

#changing input size f1 = f2
maxN <- 100
num_obs <- 10
for (size in 2:maxN) {
  for (i in 1:num_obs) {
    hm_run_time[size - 1] <- hm_run_time[size - 1] + par_haus_mat(
      spdf[start:end, ],
      f1 = 0.5,
      ncores = detectCores() - 1,
      tol = 0.01
    )
    fast_run_time[size - 1] <- fast_run_time[size - 1] + par_haus_mat_fast(
      spdf[start:end, ],
      f1 = 0.5,
      ncores = detectCores() - 1,
      tol = 0.01
    )
  }
  hm_run_time[size - 1] <- hm_run_time[size - 1] / num_obs
  fast_run_time[size - 1] <- fast_run_time[size - 1] / num_obs
  diff[size - 1] <- hm_run_time[size - 1] - fast_run_time[size - 1]
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

# using x = 2:61 and the first 60 entries of diff we get:
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
#(Intercept) -0.112486   0.142879  -0.787    0.434    
# x            0.020592   0.003975   5.181 2.91e-06 ***
# The slope is statistically significant and positive

#changing ncores

# We would expect the performance gap between parHausMatFastBoi and parHausMat
# to decrease as the number of cores increases. This is because as more cores
# are used, the number of parallel tasks per core created by parHausMat will
# decrease.

n <- 100
maxcores <- detectCores() - 1
hm_run_time <- matrix(rep(0, maxcores * n), nrow = maxcores)
fast_run_time <- matrix(rep(0, maxcores * n), nrow = maxcores)
diff <- matrix(rep(0, maxcores * n), nrow = maxcores)
for (core in 1:maxcores) {
  for (i in 1:n) {
    hm_run_time[core, i] <- par_haus_mat(
      spdf[1:50, ],
      f1 = 1,
      ncores = core,
      tol = 0.01
    )
    fast_run_time[core, i] <- par_haus_mat_fast(
      spdf[1:50, ],
      f1 = 1,
      ncores = core,
      tol = 0.01
    )
    diff[core,i] <- hm_run_Time[core, i] - fast_run_time[core, i]
  }
}
hm_means <- apply(hm_run_time, 1, mean)
fast_means <- apply(fast_run_time, 1, mean)
diff_means <- apply(diff, 1, mean)
plot(diff_means)

# Testing gDistance(hausdorff=T) vs. directHaus ---------------------------

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
     tol = tol)
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

# do hypothesis testing to check if mean(log_time_dirHaus) > mean(log_time_gDist)
# H0: mean(log_time_dirHaus) <= mean(log_time_gDist)
# H1: mean(log_time_dirHaus) > mean(log_time_gDist)
s1 <- sd(log_time_dirHaus)
s2 <- sd(log_time_gDist)
spool <- sqrt((s1^2/length(log_time_dirHaus)) + (s2^2/length(log_time_gDist)))
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
