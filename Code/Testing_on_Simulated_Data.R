


thresh = 253
# stations <- read.csv("~/Documents/GitHub/Spatial_Extreme_Value_Modeling/Data/station_info.csv")
# stat_nos <- stations[,1]

# Function to format data to stations within the regions, given...
#   window_fit: the window list subset to the moving window of choice (e.g. input window[[82]] for the final 40-yr window)
#   reg_polygons: the SpatialPolygonsDataFrame for the regions of interest
#   sim_data: 
format_data_from_sim <- function(sim_data, reg_polygons, grid_regs_df){
  # Accessing Model Fits and saving the parameter values (scale, shape, rate) for each station
  ws_scale <- NULL
  ws_shape <- NULL
  ws_rate  <- NULL
  # j = 1
  for(i in 1:ncol(sim_data)){
    fit <- extRemes::fevd(sim_data[,i], threshold=thresh, type="GP", method="MLE")
    if(is.na(fit)){
      ws_scale[i] <- ws_shape[i] <- ws_rate[i]  <- NA
    }else{
      ws_scale[i] <- fit$results$par[1]
      ws_shape[i] <- fit$results$par[2]
      ws_rate[i]  <- fit$rate
    }
    # j = j+1
  }
  
  grid_regs_df <- grid_regs_df %>%
    dplyr::mutate(scale = ws_scale, shape = ws_shape, rate = ws_rate, long = LONGITUDE, lat = LATITUDE)
  
  # add new location (long/lat) columns that we will transform to coordinates, as well as parameter values
  stations_df <- stations %>%
    dplyr::mutate(scale = ws_scale, shape = ws_shape, rate = ws_rate, long = LONGITUDE, lat = LATITUDE)
  
  stations <- project_data(stations_df)
  region_int <- sp::over(stations, reg_polygons)
  region_int$REGION <- as.factor(region_int$REGION)
  h <- model.matrix( ~ REGION - 1, data=region_int)
  inds <- as.integer(rownames(h))
  
  stations_sub_df <- stations_df %>%
    dplyr::filter(STAT_NO %in% inds) %>%
    dplyr::filter(!is.na(scale))
  
  stations_sub <- project_data(stations_sub_df)
  
  ## Have to redo intersect, because now NAs are removed. 
  ## If remove them before the intersect, the station indicies get messed up. So have to do intersect both before and after
  region_int <- sp::over(stations_sub, reg_polygons)
  region_int$REGION <- as.factor(region_int$REGION)
  h0 <- model.matrix( ~ REGION - 1, data=region_int)
  
  # Using stations_sub_df and h0 is indicator of regions
  stations_sub_by_reg <- stations_sub_df %>% dplyr::mutate(Reg_fac = region_int$REGION, Reg1 = h0[,1], Reg2 = h0[,2], Reg3 = h0[,3])
  sort_dat <- stations_sub_by_reg %>% dplyr::arrange(Reg3, Reg2) # put stations in order of which regions they fall within
  
  return(sort_dat)
}
