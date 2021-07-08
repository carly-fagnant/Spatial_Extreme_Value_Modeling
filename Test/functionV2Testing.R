# A file for testing the functions in functionV2.R

# Creating Spatial Objects ------------------------------------------------
pt1_coords <- matrix(c(0.5, 0.5), ncol=2, byrow=T)
pt2_coords <- matrix(c(2, 2), ncol=2, byrow=T)
l1_coords <- matrix(c(0, 10, 1, 1), ncol=2, byrow=T)
l2_coords <- matrix(c(3, 2, 4, 10), ncol=2, byrow=T)
p1_coords <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1, 0, 0), ncol=2, byrow=T)
p2_coords <- matrix(c(2, 0, 3, 1, 4, 0, 2, 0), ncol=2, byrow=T)
crdref <- sp::CRS("+init=epsg:2278")

pt1 <- sp::SpatialPoints(pt1_coords, proj4string=crdref)
pt2 <- sp::SpatialPoints(pt2_coords, proj4string=crdref)
l1 <- raster::spLines(l1_coords, crs=crdref)
l2 <- raster::spLines(l2_coords, crs=crdref)
p1 <- raster::spPolygons(p1_coords, crs=crdref)
p2 <- raster::spPolygons(p2_coords, crs=crdref)

inputs <- c(pt1, pt2, l1, l2, p1, p2)
#specify CRS using the EPSG code but without using 
# proj4strings like +init=epsg:????, +proj=longlat

for (i in 1:6) {
  print(i)
  print(sp::is.projected(inputs[[i]]))
  sp::proj4string(inputs[[i]]) <- crdref
  inputs[[i]] <- sp::spTransform(inputs[[i]], proj4string(inputs[[i]]))
  print(sp::is.projected(inputs[[i]]))
  print("------------------------------------------------------")
}

# Testing Hausdorff -------------------------------------------------------

cases1 <- c(pt1, pt1, pt1, l1, l1, p1)
cases2 <- c(pt2, l1, p1, l2, p1, p2)

# Compute preliminary results
f1_mat <- matrix(rep(-1, 6*4), nrow=6)
f2_mat <- matrix(rep(-1, 6*4), nrow=6)
for (input_case in 1:6) {
  input1 <- cases1[[input_case]]
  input2 <- cases2[[input_case]]
  n <- 100000
  a.coords <- spsample(input1, n=n, type="regular") # sample points within A
  distsf1 <- gDistance(a.coords, input2, byid=T)
  b.coords <- spsample(input2, n=n, type="regular")
  distsf2 <- gDistance(b.coords, input1, byid=T)
  f1_mat[input_case,] <- quantile(distsf1, probs=c(0.25, 0.5, 0.75, 1))
  f2_mat[input_case,] <- quantile(distsf2, probs=c(0.25, 0.5, 0.75, 1))
}

# Compute final answers
answers <- matrix(rep(-1, 6*11), nrow=6)
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
      val <- max(gDistance(input1, input2), f2.50)
    } else if (test_case == 2) {
      val <- max(f1.50, gDistance(input1, input2))
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
      val <- gDistance(input1, input2, hausdorff=T)
    }
    answers[input_case, test_case] <- val
  }
}
answers[3, ] <- rep(0, 11)

for (i in 1:6) {
  flag <- 0
  print(paste0("Testing input set ", i))
  print(paste0("Class of A: ",  class(cases1[[i]])[1]))
  print(paste0("Class of B: ",  class(cases2[[i]])[1]))
  
  ## Case 1: f1 = 0
  test1.1 <- extHaus(cases1[[i]], cases2[[i]], f1=0, f2=0.5)
  if (abs(test1.1 - answers[i, 1]) > 0.01) {
    flag <- 1
    print(paste0("Test case 1 failed on input pair set: ", i))
    print(paste0("Expected ", answers[i, 1], " but computed ", test1.1))
  }
  
  ## Case 2: f2 = 0
  test2.1 <- extHaus(cases1[[i]], cases2[[i]], f1=0.5, f2=0)
  if (abs(test2.1 - answers[i, 2]) > 0.01) {
    flag <- 1
    print(paste0("Test case 2 failed on input pair set: ", i))
    print(paste0("Expected ", answers[i, 2], " but computed ", test2.1))
  }
  
  ## Case 3: f1 = f2 = 0
  test3.1 <- extHaus(cases1[[i]], cases2[[i]], f1=0)
  if (abs(test3.1 - answers[i, 3]) > 0.01) {
    flag <- 1
    print(paste0("Test case 3 failed on input pair set: ", i))
    print(paste0("Expected ", answers[i, 3], " but computed ", test3.1))
  }
  
  ## Case 4: 0 < f1, f2 < 1
  test4.1 <- extHaus(cases1[[i]], cases2[[i]], f1=0.25, f2=0.75)
  test4.2 <- extHaus(cases1[[i]], cases2[[i]], f1=0.75, f2=0.25)
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
  test5.1 <- extHaus(cases1[[i]], cases2[[i]], f1=0.25)
  test5.2 <- extHaus(cases1[[i]], cases2[[i]], f1=0.5)
  test5.3 <- extHaus(cases1[[i]], cases2[[i]], f1=0.75)
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
  test6.1 <- extHaus(cases1[[i]], cases2[[i]], f1=1, f2=0.5)
  if (abs(test6.1 - answers[i, 9]) > 0.01) {
    flag <- 1
    print(paste0("Test case 6 failed on input pair set: ", i))
    print(paste0("Expected ", answers[i, 9], " but computed ", test6.1))
  }
  
  ## Case 7: f2 = 1
  test7.1 <- extHaus(cases1[[i]], cases2[[i]], f1=0.5, f2=1)
  if (abs(test7.1 - answers[i, 10]) > 0.01) {
    flag <- 1
    print(paste0("Test case 6 failed on input pair set: ", i))
    print(paste0("Expected ", answers[i, 10], " but computed ", test7.1))
  }
  
  ## Case 8: f1 = f2 = 1
  test8.1 <- extHaus(cases1[[i]], cases2[[i]], f1=1)
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

# Testing Hausdorff: Census Tracts and Rivers -----------------------------
setwd("Data")
## Read shp file
tracts.houston <- rgdal::readOGR("Census_2010_Tracts.shp")
## specifying projection information for SpatialPolygonsDataFrame object
tracts.houston <- sp::spTransform(tracts.houston, sp::CRS(sp::proj4string(tracts.houston)))
## extract tracts of Harris County
tracts.harris <- tracts.houston[grepl(c("201"),
                                      tracts.houston@data$COUNTY), ]
tracts.harris@data <- tracts.harris@data[ ,c("OBJECTID", "STATE",
                                             "COUNTY","TRACT", "SUM_TotPop")]
tracts.harris <- sp::spTransform(tracts.harris, CRS("+init=epsg:2278"))

rivers <- rgdal::readOGR("Major_Rivers.gdb")
rivers <- sp::spTransform(rivers, CRS("+init=epsg:2278"))
set.seed(123)
f1_vec <- rbeta(length(tracts.harris)^2, shape1=0.5, shape2=0.5)
f2_vec <- rbeta(length(tracts.harris)^2, shape1=0.5, shape2=0.5)

for (i in 1:10) {
  #for (j in 1:length(tracts.harris)) {
    #print(paste0("i: ", i))
    #print(paste0("j: ", j))
    #print(paste0("f1_vec: ", f1_vec[i*j]))
    #print(paste0("f2_vec: ", f2_vec[i*j]))
    #extHaus(tracts.harris[i,], tracts.harris[j,], f1_vec[i*j], f2_vec[i*j])
  #}
  
  for (k in 1:10) {
    print(paste0("i: ", i))
    print(paste0("k: ", k))
    print(paste0("f1_vec: ", f1_vec[i*k]))
    print(paste0("f2_vec: ", f2_vec[i*k]))
    extHaus(tracts.harris[i,], rivers[k,], f1_vec[i*k], f2_vec[i*k])
  }
}

# extHaus(tracts.harris[5,], rivers[5,], 10^-18) = 249125.8
# extHaus(tracts.harris[5,], rivers[5,], 10^(-18), tol=0.001) = 247208.6 
# gDistance(tracts.harris[5,], rivers[5,]) = 246964.7

# Testing Hausdorff: hausMat -----------------------------------------------

# line-line
mat1 <- hausMat(rivers[1:5,], f1=0.5, ncores=detectCores()-1, tol=0.001)
mat2 <- hausMat(rivers[1:5,], f1=0.75, f2=0.25, ncores=detectCores()-1, tol=0.001)
for (i in 1:5) {
  for (j in 1:5) {
    val1 <- extHaus(rivers[i,], rivers[j,], 0.5, tol=0.001)
    quant1 <- abs(mat1[i,j] - val1)
    if (quant1 > 0.001 * val1) {
      print(paste0("Error in mat1 iteration (", i, ", ", j, ")"))
      print(paste0("Absolute error: ", quant1))
      print(paste0("Fractional error: ", quant1 / val1))
    }
    val2 <- extHaus(rivers[i,], rivers[j,], 0.75, 0.25, tol=0.001)
    quant2 <- abs(mat2[i,j] - val2)
    if (quant2 > 0.001 * val2) {
      print(paste0("Error in mat2 iteration (", i, ", ", j, ")"))
      print(paste0("Absolute error: ", quant2))
      print(paste0("Fractional error: ", quant2 / val2))
    }
    print("---------------------------------------------------")
  }
}

# area-area
mat1 <- hausMat(tracts.harris[1:10,], f1=0.5, ncores=detectCores()-1, tol=0.001)
mat2 <- hausMat(tracts.harris[1:10,], f1=0.75, f2=0.25, ncores=detectCores()-1, tol=0.001)
for (i in 1:10) {
  for (j in 1:10) {
    val1 <- extHaus(tracts.harris[i,], tracts.harris[j,], 0.5, tol=0.001)
    quant1 <- abs(mat1[i,j] - val1)
    if (quant1 > 0.001 * val1) {
      print(paste0("Error in mat1 iteration (", i, ", ", j, ")"))
      print(paste0("Absolute error: ", quant1))
      print(paste0("Fractional error: ", quant1 / val1))
    }
    val2 <- extHaus(tracts.harris[i,], tracts.harris[j,], 0.75, 0.25, tol=0.001)
    quant2 <- abs(mat2[i,j] - val2)
    if (quant2 > 0.001 * val2) {
      print(paste0("Error in mat2 iteration (", i, ", ", j, ")"))
      print(paste0("Absolute error: ", quant2))
      print(paste0("Fractional error: ", quant2 / val2))
    }
    print("---------------------------------------------------")
  }
}
# Testing Hausdorff: hausMat vs hausMatFastBoi ----------------------------

n <- 100
coords <- matrix(c(1, 1), nrow=n, ncol=n, byrow=T)
data <- as.data.frame(matrix(rep(0, n*n), nrow=n, ncol=n))
for (i in 1:n) {
  coords <- rbind(coords, matrix(runif(2), ncol=2, byrow=T))
}
coords <- as.data.frame(coords)
crdref <- sp::CRS("+init=epsg:2278")
spdf <- SpatialPointsDataFrame(coords=coords, data=data, proj4string=crdref)

hmRunTime <- rep(0, n-1)
fastBoiRunTime <- rep(0, n-1)
diff <- rep(0, n-1)
f1vec <- seq(from=0, to=1, length.out=n) 

#changing input size f1 does not equal f2
maxN <- 100
for (size in 2:maxN) {
  start <- 1
  end <- start + size
  hmRunTime[size - 1] <- parHausMat(spdf[start:end,], f1 = f1vec[size - 1],
                                    f2 = 1 - f1vec[size - 1],
                                    ncores = detectCores() - 1, tol = 0.01)
  fastBoiRunTime[size - 1] <- parHausMatFastBoi(spdf[start:end,],
                                              f1 = f1vec[size - 1],
                                              f2 = 1 - f1vec[size - 1],
                                              ncores=detectCores() - 1,
                                              tol = 0.01)
  diff[size - 1] <- hmRunTime[size - 1] - fastBoiRunTime[size - 1]
}
plot(hmRunTime[1:maxN - 1])
plot(fastBoiRunTime[1:maxN - 1])
plot(diff[1:maxN - 1])
min(diff)
mean(diff)
max(diff)
boxplot(hmRunTime[1:maxN - 1])
boxplot(fastBoiRunTime[1:maxN - 1])
boxplot(diff[1:maxN - 1])
hist(diff[1:maxN - 1])

#changing input size f1 = f2
maxN <- 100
for (size in 2:maxN) {
  start <- 1
  end <- start + size
  hmRunTime[size - 1] <- parHausMat(spdf[start:end,], f1 = f1vec[size - 1],
                                    ncores = detectCores() - 1, tol = 0.01)
  fastBoiRunTime[size - 1] <- parHausMatFastBoi(spdf[start:end,],
                                                f1 = f1vec[size - 1],
                                                ncores = detectCores() - 1,
                                                tol = 0.01)
  diff[size - 1] <- hmRunTime[size - 1] - fastBoiRunTime[size - 1]
}
plot(hmRunTime[1:maxN - 1])
plot(fastBoiRunTime[1:maxN - 1])
plot(diff[1:maxN - 1])
min(diff)
mean(diff)
max(diff)
boxplot(hmRunTime[1:maxN - 1])
boxplot(fastBoiRunTime[1:maxN - 1])
boxplot(diff[1:maxN - 1])
hist(diff[1:maxN - 1])
