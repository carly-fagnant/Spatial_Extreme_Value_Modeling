###############################
# Example Code - Fitting Data to GPD and Moving Windows
#
# by Carly Fagnant
###############################

# Recall PRCP = Precipitation is measured as "tenths of mm"   (254 = 1 inch)
# To get inches: x/254
# To get mm: x/10

library(lubridate)
library(extRemes)
library(stats)
library(dplyr)
library(xts)
library(gnFit)
# library(ismev)
# library(tseries)
# library(trend)
# library(astsa)
# library(ggmap)


setwd("~/Documents/GitHub/scrape-NOAA-data/Data")
# precip <- read.csv("all_stations_precip_UDP_updated_2021.csv")  # new data you will use

precip <- read.csv("all_stations_precip_UDP.csv", header = TRUE)
precip <- precip[,-1]   # delete first column since is repeat of row names (numbers)

### Declustering
precip_dclust = precip
for ( i in seq(from=2, to=length(precip[1,]))){
  station = precip[,i]
  non_na = which(station>=0)
  station[is.na(station)] <- 0
  
  dec = extRemes::decluster(station, threshold=0, clusterfun = "max", replace.with = 0)
  stat_clus = as.numeric(dec)
  
  station = stat_clus
  precip_dclust[non_na,i] = station[non_na]
}

### Save as data and format
data <- precip_dclust
data[,1] <- lubridate::ymd(as.character(data[,1])) #put Date column into date format
data <- as.data.frame(data) #make it a data frame

### Set threshold
thresh <- 253   # for 1-day


# Distribution Fits -------------------------------------------------------

numstat <- numrow <- 601 #stations
numcol <- 79 #moving windows (this number will change with the new data since we have more years)

# GPD Fit function
# Given a matrix of data, fits a GPD to each column of data and saves the fit objects in a list
fitgpdR <- function(data, thresh){
  gpdfits <- list() #list of gpd fits
  
  for(i in 1:numstat){
    column <- 1+i # add 1 since first column is the date
    station_i <- na.omit(data[, column])
    if(class(try(extRemes::fevd(station_i, threshold=thresh, type="GP", method="MLE")))=="try-error"){ 
      gpdfits[[i]] <- NA
    }else{
      gpdfits[[i]] <- extRemes::fevd(station_i, threshold=thresh, type="GP", method="MLE")
    }
  }  
  return(gpdfits)
}

# Fit to entire station length --------------------------------------------

# Fitting GPD to entire length of data for each station
fullfits <- fitgpdR(data, thresh)
# which(is.na(fullfits))  # only a few are NA





# Goodness of Fit (gof) functions ------------------------------------------------------------

# gof is the gnFit::gnfit function altered to not plot QQ and hist, and to save p-value as 0 if significantly small (not just give warning)
gof <- function (dat, dist, df = NULL, pr = NULL, threshold = NULL) 
{
  dat <- as.numeric(dat)
  x <- NULL
  z <- list()
  op <- par(mfrow = c(1, 2))
  if (is.null(pr)) {
    if (dist == "gev" | dist == "gpd") 
      stop("Enter Parameters!")
    else if (dist == "t") {
      if (df > 2) {
        loc <- mean(dat)
        sc <- sqrt((df - 2) * var(dat)/df)
        xdat <- (dat - loc)/sc
        prob <- pt(xdat, df)
        qqplot(qt(ppoints(500), df), xdat, main = paste0("Q-Q plot for ", 
                                                         dist, ". distribution-", "DF:", df), xlab = "", 
               ylab = "")
        qqline(xdat, distribution = function(p) qt(p, 
                                                   df), probs = c(0.1, 0.6), col = "blue")
        hist(xdat, breaks = 20, prob = TRUE, xlab = "x-normalized value", 
             main = paste0(dist, ". distribution curve over histogram"))
        curve(dt(x, df), col = "blue", add = TRUE, yaxt = "n")
        points(xdat, rep(0, length(xdat)))
      }
      else stop("DF must be > 2")
    }
    else if (dist == "gum") {
      sc <- sqrt(6 * var(dat)/pi^2)
      loc <- mean(dat) - 0.577 * sc
      pr <- c(loc, sc, 0)
      prob <- gevf(pr, dat)
      gev.qq(pr, dat)
      gev.his(pr, dat)
    }
    else {
      loc <- mean(dat)
      ifelse(dist == "norm", sc <- sd(dat), NA)
      ifelse(dist == "laplace", sc <- sqrt(var(dat)/2), 
             NA)
      ifelse(dist == "logis", sc <- sqrt(3 * var(dat)/pi^2), 
             NA)
      prob <- get(paste0("p", dist))(dat, loc, sc)
      qqplot(get(paste0("q", dist))(ppoints(500), loc,
                                    sc), dat, main = paste0("Q-Q plot for ", dist,
                                                            ". distribution"), xlab = "", ylab = "")
      qqline(dat, distribution = function(p) get(paste0("q",
                                                        dist))(p, loc, sc), probs = c(0.1, 0.6), col = "blue")
      hist(dat, breaks = 15, prob = TRUE, xlab = "x-normalized value",
           main = paste0(dist, ". distribution curve over histogram"))
      curve(get(paste0("d", dist))(x, loc, sc), col = "blue",
            add = TRUE, yaxt = "n")
      points(dat, rep(0, length(dat)))
    }
  }
  else {
    if (dist == "gev") {
      prob <- gevf(pr, dat)
      gev.qq(pr, dat)
      gev.his(pr, dat)
    }
    else if (dist == "gum") {
      pr[3] <- 0
      prob <- gevf(pr, dat)
      gev.qq(pr, dat)
      gev.his(pr, dat)
    }
    else if (dist == "gpd") 
      if (!is.null(threshold)) {  # if pr given in evir format
        # pr <- c(pr[[2]], pr[[1]]) ####this is the line added to convert form parameters given in
        u <- threshold
        dat <- dat[dat > u]
        prob <- gpdf(pr, u, dat)
        # prob[which(prob==1)] <- 0.99999999  # added b/c prob=1 gives Inf
        par(mfrow=c(1,2),oma=c(0,0,2,0)) 
        # gpd.qq(pr, u, dat)
        # gpd.his(pr, u, dat)
        # title(stationTitle, outer=TRUE)
        
      }
    else stop("threshold is missing!")
  }
  n <- length(dat)
  k <- seq(1:n)
  qnor <- qnorm(sort(prob))
  pnor <- pnorm((qnor - mean(qnor))/sd(qnor))
  w <- round((sum((pnor - (2 * k - 1)/(2 * n))^2) + 1/(12 * 
                                                         n)) * (1 + 0.5/n), 4)
  if (w < 0.0275) {
    pval <- 1 - exp(-13.953 + 775.5 * w - 12542.61 * w^2)
  }
  else if (w < 0.051) {
    pval <- 1 - exp(-5.903 + 179.546 * w - 1515.29 * w^2)
  }
  else if (w < 0.092) {
    pval <- exp(0.886 - 31.62 * w + 10.897 * w^2)
  }
  else if (w < 1.1) {
    pval <- exp(1.111 - 34.242 * w + 12.832 * w^2)
  }
  else {
    warning("p-value is smaller than 7.37e-10")
    pval <- 0 ## Added in so can calculate percentages of stations with good fit
  }
  z$Wpval <- pval
  A <- (-n - sum((2 * k - 1) * log(pnor) + (2 * n + 1 - 2 * 
                                              k) * log(1 - pnor))/n) * (1 + 0.75/n + 2.25/n^2)
  A <- round((1 + 0.75/n + 2.25/n^2) * A, 4)
  if (A < 0.2) {
    pval <- 1 - exp(-13.436 + 101.14 * A - 223.73 * A^2)
  }
  else if (A < 0.34) {
    pval <- 1 - exp(-8.318 + 42.796 * A - 59.938 * A^2)
  }
  else if (A < 0.6) {
    pval <- exp(0.9177 - 4.279 * A - 1.38 * A^2)
  }
  else if (A < 10) {
    pval <- exp(1.2937 - 5.709 * A + 0.0186 * A^2)
  }
  else {
    warning("p-value is smaller than 7.37e-10")
    pval <- 0 ## Added in so can calculate percentages of stations with good fit
  }
  z$Apval <- pval
  z$Cram <- w
  z$Ander <- A
  message("Test of Hypothesis for ", dist, " distribution")
  message("Cramer-von Misses Statistics:  ", z$Cram, "   P-Value:  ", 
          round(z$Wpval, 5))
  message("Anderson-Darling Statistics:   ", z$Ander, "   P-Value:  ", 
          round(z$Apval, 5))
  class(z) <- "gnfit"
  invisible(z)
}


gof_fix_error <- function (dat, dist, df = NULL, pr = NULL, threshold = NULL) 
{
  dat <- as.numeric(dat)
  x <- NULL
  z <- list()
  op <- par(mfrow = c(1, 2))
  if (is.null(pr)) {
    if (dist == "gev" | dist == "gpd") 
      stop("Enter Parameters!")
    else if (dist == "t") {
      if (df > 2) {
        loc <- mean(dat)
        sc <- sqrt((df - 2) * var(dat)/df)
        xdat <- (dat - loc)/sc
        prob <- pt(xdat, df)
        qqplot(qt(ppoints(500), df), xdat, main = paste0("Q-Q plot for ", 
                                                         dist, ". distribution-", "DF:", df), xlab = "", 
               ylab = "")
        qqline(xdat, distribution = function(p) qt(p, 
                                                   df), probs = c(0.1, 0.6), col = "blue")
        hist(xdat, breaks = 20, prob = TRUE, xlab = "x-normalized value", 
             main = paste0(dist, ". distribution curve over histogram"))
        curve(dt(x, df), col = "blue", add = TRUE, yaxt = "n")
        points(xdat, rep(0, length(xdat)))
      }
      else stop("DF must be > 2")
    }
    else if (dist == "gum") {
      sc <- sqrt(6 * var(dat)/pi^2)
      loc <- mean(dat) - 0.577 * sc
      pr <- c(loc, sc, 0)
      prob <- gevf(pr, dat)
      gev.qq(pr, dat)
      gev.his(pr, dat)
    }
    else {
      loc <- mean(dat)
      ifelse(dist == "norm", sc <- sd(dat), NA)
      ifelse(dist == "laplace", sc <- sqrt(var(dat)/2), 
             NA)
      ifelse(dist == "logis", sc <- sqrt(3 * var(dat)/pi^2), 
             NA)
      prob <- get(paste0("p", dist))(dat, loc, sc)
      qqplot(get(paste0("q", dist))(ppoints(500), loc,
                                    sc), dat, main = paste0("Q-Q plot for ", dist,
                                                            ". distribution"), xlab = "", ylab = "")
      qqline(dat, distribution = function(p) get(paste0("q",
                                                        dist))(p, loc, sc), probs = c(0.1, 0.6), col = "blue")
      hist(dat, breaks = 15, prob = TRUE, xlab = "x-normalized value",
           main = paste0(dist, ". distribution curve over histogram"))
      curve(get(paste0("d", dist))(x, loc, sc), col = "blue",
            add = TRUE, yaxt = "n")
      points(dat, rep(0, length(dat)))
    }
  }
  else {
    if (dist == "gev") {
      prob <- gevf(pr, dat)
      gev.qq(pr, dat)
      gev.his(pr, dat)
    }
    else if (dist == "gum") {
      pr[3] <- 0
      prob <- gevf(pr, dat)
      gev.qq(pr, dat)
      gev.his(pr, dat)
    }
    else if (dist == "gpd") 
      if (!is.null(threshold)) {  # if pr given in evir format
        # pr <- c(pr[[2]], pr[[1]]) ####this is the line added to convert from parameters given in evir format?
        u <- threshold
        dat <- dat[dat > u]
        prob <- gpdf(pr, u, dat)
        if(any(prob==1) | any(is.nan(prob))){
          prob <- gpdf(round(pr, 12), u, dat)
          if(any(prob==1) | any(is.nan(prob))){
            prob <- gpdf(round(pr, 11), u, dat)
            if(any(prob==1) | any(is.nan(prob))){
              prob <- gpdf(round(pr, 10), u, dat)
              if(any(prob==1) | any(is.nan(prob))){
                prob <- gpdf(round(pr, 9), u, dat)
                if(any(prob==1) | any(is.nan(prob))){
                  prob <- gpdf(round(pr, 8), u, dat)
                  if(any(prob==1) | any(is.nan(prob))){
                    prob <- gpdf(round(pr, 7), u, dat)
                    if(any(prob==1) | any(is.nan(prob))){
                      prob <- gpdf(round(pr, 6), u, dat)
                    }
                  }
                }
              }
            }
          }
        }
        # prob[which(prob==1)] <- 0.99999999  # added b/c prob=1 gives Inf
        par(mfrow=c(1,2),oma=c(0,0,2,0)) 
        # gpd.qq(pr, u, dat)
        # gpd.his(pr, u, dat)
        # title(stationTitle, outer=TRUE)
        
      }
    else stop("threshold is missing!")
  }
  n <- length(dat)
  k <- seq(1:n)
  qnor <- qnorm(sort(prob))
  pnor <- pnorm((qnor - mean(qnor))/sd(qnor))
  w <- round((sum((pnor - (2 * k - 1)/(2 * n))^2) + 1/(12 * 
                                                         n)) * (1 + 0.5/n), 4)
  if (w < 0.0275) {
    pval <- 1 - exp(-13.953 + 775.5 * w - 12542.61 * w^2)
  }
  else if (w < 0.051) {
    pval <- 1 - exp(-5.903 + 179.546 * w - 1515.29 * w^2)
  }
  else if (w < 0.092) {
    pval <- exp(0.886 - 31.62 * w + 10.897 * w^2)
  }
  else if (w < 1.1) {
    pval <- exp(1.111 - 34.242 * w + 12.832 * w^2)
  }
  else {
    warning("p-value is smaller than 7.37e-10")
    pval <- 0 ## Added in so can calculate percentages of stations with good fit
  }
  z$Wpval <- pval
  A <- (-n - sum((2 * k - 1) * log(pnor) + (2 * n + 1 - 2 * 
                                              k) * log(1 - pnor))/n) * (1 + 0.75/n + 2.25/n^2)
  A <- round((1 + 0.75/n + 2.25/n^2) * A, 4)
  if (A < 0.2) {
    pval <- 1 - exp(-13.436 + 101.14 * A - 223.73 * A^2)
  }
  else if (A < 0.34) {
    pval <- 1 - exp(-8.318 + 42.796 * A - 59.938 * A^2)
  }
  else if (A < 0.6) {
    pval <- exp(0.9177 - 4.279 * A - 1.38 * A^2)
  }
  else if (A < 10) {
    pval <- exp(1.2937 - 5.709 * A + 0.0186 * A^2)
  }
  else {
    warning("p-value is smaller than 7.37e-10")
    pval <- 0 ## Added in so can calculate percentages of stations with good fit
  }
  z$Apval <- pval
  z$Cram <- w
  z$Ander <- A
  message("Test of Hypothesis for ", dist, " distribution")
  message("Cramer-von Misses Statistics:  ", z$Cram, "   P-Value:  ", 
          round(z$Wpval, 5))
  message("Anderson-Darling Statistics:   ", z$Ander, "   P-Value:  ", 
          round(z$Apval, 5))
  class(z) <- "gnfit"
  invisible(z)
}



# Now going to do GOF for ALL stations, for their fits to the entire station length

data_no_na <- function(data, thresh){
  dat <- list() #list of data for each
  
  for(i in 1:numstat){
    column <- 1+i # add 1 since first column is the date
    station_i <- na.omit(data[, column])
    if(class(try(extRemes::fevd(station_i, threshold=thresh, type="GP", method="MLE")))=="try-error"){ 
      dat[[i]] <- NA
    }else{
      dat[[i]] <- station_i
    }
  }  
  return(dat)
}

fulldata_no_na <- data_no_na(data, thresh)  # saving data corresponding to fullfits


CVMp <- ADp <- NULL
for(h in c(1:numstat)){ # for some reason 481 was not working with gof fn, so use gof_fix_error  #but it has <1 year in length so will be eliminated anyway
  if(!is.na(fullfits[[h]])){
    gof_h <- gof_fix_error(fulldata_no_na[[h]], dist="gpd", pr=fullfits[[h]]$results$par, threshold=thresh)
    CVMp[h] <- gof_h$Wpval
    ADp[h]  <- gof_h$Apval
  }else{
    CVMp[h] <- ADp[h]  <- NA
  }
}



# All Stations - 40-Year Windows ------------------------------------------

# Eventually we will want these fits for all of the 40-year windows, but for now we can just fit to the most recent 40 years.


### Uncomment if you want to run from scratch- Note: will take about an hour
### You can instead load the .rds file below to load a saved copy (see readRDS line below)

# ### 40-YEAR MOVING WINDOWS ///////////////
numcol <- 79 #windows of 40
# startday <- "-01-01"
# endday <- "-12-31"
# window <- list()
# labelyr <- NULL   # end year of window
# 
# # to access window j, gpdfit for station i, use window[[j]][[i]]
# 
# startyr <- 1900
# for(j in 1:numcol){
#   start <- paste0(startyr, startday)
#   endyr <- startyr + 39
#   end <- paste0(endyr, endday)
#   labelyr[j] <- lubridate::year(as.Date(end)) 
#   
#   start <- which(data$Date==start)  #finding indexes corresponding to start & end dates
#   end   <- which(data$Date==end)
#   sub <- data[start:end, ]  #subset data to those 40 years
#   window[[j]] <- fitgpdR(sub, thresh)   # GPD fit with extRemes package
#   
#   startyr <- startyr + 1
# }
# 
# # saveRDS(window, file="window_1day_dclust.rds")

window <- readRDS(file="window_1day_dclust.rds") # get this file from Box, should be 1.67 GB .rds file when unzipped


# Saving Return Levels for Easy Plotting ----------------------------------

# in mm
# this is set up differently for easier plotting
# saving vector of RLs for each station across 79 windows

trend_RL <- list()
for(i in 1:numstat){
  RL <- NULL
  for(j in 1:numcol){
    fit <- window[[j]][[i]]
    if(!is.na(fit) == TRUE){
      RL[j] <- extRemes::return.level(fit, return.period=100)/10  # divide by 10 to get mm
    }else{
      RL[j] <- NA
    }
  }
  trend_RL[[i]] <- RL
}

# trend_RL[[588]] # Hobby

trend_RL_25 <- list()
trend_RL_100 <- list()
trend_RL_500 <- list()
for(i in 1:numstat){
  RL_25 <- RL_100 <- RL_500 <- NULL
  for(j in 1:numcol){
    fit <- window[[j]][[i]]
    if(!is.na(fit) == TRUE){
      RL_25[j] <- extRemes::return.level(fit, return.period=25)/10  # divide by 10 to get mm
      RL_100[j] <- extRemes::return.level(fit, return.period=100)/10  # divide by 10 to get mm
      RL_500[j] <- extRemes::return.level(fit, return.period=500)/10  # divide by 10 to get mm
    }else{
      RL_25[j] <- RL_100[j] <- RL_500[j] <- NA
    }
  }
  trend_RL_25[[i]]  <- RL_25
  trend_RL_100[[i]] <- RL_100
  trend_RL_500[[i]] <- RL_500
}


labelyr <- 1939:2017