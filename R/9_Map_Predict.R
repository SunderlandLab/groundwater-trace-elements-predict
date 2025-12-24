###############################################################################
#File name: Map_Predict.R
#Author: Jennifer Sun
#Date: Nov 2025
#Purpose: Predict concentrations on a continuous 0.1x0.1 degree grid across the US
##############################################################################


#### Attach Libraries and Set Working Directory ####
packages <- c('plyr','dplyr','pROC','rpart','stringr','tidyverse','tidymodels','xgboost','rgdal','sp')
lapply(packages, library, character.only=TRUE)

setwd(here::here("data"))

#### Define Model Version ####
metal.codes <- c("As", "Cd", "Li", "Mn", "Sr")
metals <- c("Arsenic", "Cadmium", "Lithium", "Manganese", "Strontium")
names(metals) <- metal.codes

if (metal.code=='Cd') {
  load(paste0(metal.code, "_ModelPackage_2step.RData"))
} else {
  load(paste0(metal.code, "_ModelPackage.RData"))
}
load(filename)

## 1. Create 0.1 x 0.1 degree grid with predictors across the US ----------------------------------------------------------------------------------------------
# Load in Map of US
mapUSm<- readOGR(dsn="CoVar/US48", layer="US_48states") # projected map of the U.S.(in meters)

create_grid <- function(metal.code) {
  # Create Grid over CONUS
  grid = makegrid(mapUSm, cellsize = 0.1) # cell size is in map units, which here is lat/long degrees
  grid = SpatialPoints(grid, proj4string = crs(mapUSm)) # 775 x 1254 points, grid dimensions
  
  # Cookie cut out extra points
  inside.US<-!is.na(over(grid,as(mapUSm,'SpatialPolygons')))
  grid$insideUS<-inside.US
  grid<-grid[grid$insideUS==1,]
  
  # Extract predictor variables (use 2_Predictors_Extract)
  grid <- readRDS('R_Output/MapGrid/mapgrid_Model_Ready.rds') # saved grid with predictors
  
  ## 2. Prepare separate grid for each metal to be compatible with xgboost model inputs ----------------------------------------------------------------------------
  # Rename TRI and SEMS columns you want to use 
  grid@data = grid@data %>% dplyr::rename('TRI.total.impact' = paste0('TRI.total.impact.',metal.code),
                                          'TRI.water.impact' = paste0('TRI.water.impact.',metal.code),
                                          'SEMS' = paste0('SEMS.',metal.code))
  
  # Add back in non-predictor columns so the model workflow can run; these columns are discarded before prediction
  newcols = setdiff(colnames(df), colnames(grid@data))
  for (col in newcols) {
    grid@data = grid@data %>% dplyr::mutate(!!col := 0)
  }
  
  # Remove unnecessary new columns
  grid@data = grid@data %>% dplyr::select(-gridIndex)
  
  # Manually impute factor values that don't exist in training data
  if (metal.code %in% c('As','Mn','Cd','Li','Sr')) {
    grid@data['aquifer'][grid@data['aquifer'] == 402] <- 999 # no information
    grid@data['aquifer'][grid@data['aquifer'] == 406] <- 999 # no information
    grid@data['aquifer'][grid@data['aquifer'] == 601] <- 999 # no information
    grid@data['str_int'][grid@data['str_int'] == 5] <- 4 # Next closest category
  } 
  
  if (metal.code %in% c('As','Mn','Cd','Li')) {
    grid@data['landcover_500m'][grid@data['landcover_500m'] == 12] <- 11 # there are not many of these values; perennial ice/snow counted as open water
    grid@data['landcover_500m'][grid@data['landcover_500m'] == 0] <- 11 # there are not many of these values; no data counted as open water 
  }
  
  # add in dummy columns to match MICE model input
  grid$lon <- grid$long
  grid$.imp <- NA
  grid$.id <- NA
  
  # match data types with model input
  grid$location.id = as.character(grid$location.id) 
  grid$censored = as.logical(grid$location.id) 
  grid$data.source = as.factor(grid$location.id) 
  # grid$HLR = as.factor(grid$HLR) 
  grid$date = as.POSIXct(grid$date)
  grid$landcover_500m = as.factor(grid$landcover_500m)
  grid$aquifer = as.factor(grid$aquifer)
  grid$aq_rocktype = as.factor(grid$aq_rocktype)
  grid$drainage_class_int = as.factor(grid$drainage_class_int)
  grid$hydgrp_int = as.factor(grid$hydgrp_int)
  grid$weg_int = as.factor(grid$weg_int)
  grid$str_int = as.factor(grid$str_int)
  grid$soilorder_int = as.factor(grid$soilorder_int)
  grid$lith = as.factor(grid$lith)
  
  # save grid
  df.grid = grid@data
  df.grid$longitude = grid@coords[,1]
  df.grid$latitude = grid@coords[,2]
  
  saveRDS(df.grid, paste0(metal.code, '_MapGrid.rds'))
  
}

## Predict on map grid --------------------------------------------------------------------------------------------------------------------------------------------------------------
predict_on_grid <- function(metal.code, grid) {
  grid = readRDS(paste0(metal.code, '_MapGrid.rds'))
  
  ## 3. Predict detections on grid (Cd only) ----------------------------------------------------------------------------------------------------------------------------------------
  
  if (metal.code == 'Cd') {
    load('Cd_DetectionModel.RData')
    load(paste0(metal.code, "_ModelPackage_2step.RData"))
    
    grid.detectprobs = detection_model %>% 
      map( ~ predict(.x, new_data=grid, type='prob') %>% pull(.pred_TRUE)) %>%  # model prediction
      as_tibble(.name_repair = ~ paste0("p_model", seq_along(.))) %>%  # save results in a tibble
      mutate(p_mean = rowMeans(across(everything()))) %>% # take the average of the predictions 
      mutate(detect_class = as.integer(p_mean >=0.29)) # create filter for detections
    
    # add predictions back to the grid
    grid$pred_detectpmean = grid.detectprobs$p_mean
    grid$pred_detectclass = grid.detectprobs$detect_class
  }


  ## 4. Predict concentrations on grid ------------------------------------------------------------------------------------------------------------------------------------------------------

  # A. Predict --------------------------------------------------------------------------------------------------------------------------------------------------------------
  grid.predictions = final_full_model %>% 
    map( ~ predict(.x, new_data=grid)) # note that the .pred column is a log concentration
  
  # B. EDM adjustment --------------------------------------------------------------------------------------------------------------------------------------------------------
  
  train_predict <- function(model, df.train) {
    # prepare training predictions
    df.train.filtered = subset(df.train, is.imputed==0) # test predictions similarly done only on un-imputed values in script #5
    train.pred.filtered = predict(model, df.train.filtered) %>% bind_cols(df.train.filtered)
    
    return(train.pred.filtered) # also named train.pred.adj
  }
  
  EDM_transform <- function(test.pred, train.pred.filtered) { # grid predictions, train.pred.filtered; 
    # NOTE: changed inputs from train.pred.adj and test.pred.adj to current versions to reflect the contents of the INPUT rather than the intended contents of the OUTPUT 
    
    # order train values
    train.pred.filtered$logconc_ordered = sort(train.pred.filtered$logconc) # order concentrations
    train.pred.filtered$.pred_ordered = sort(train.pred.filtered$.pred) # order train predictions
    tbl_ordered = train.pred.filtered[,c('logconc_ordered','.pred_ordered')]
    
    # transformation: adjust test predictions based on the difference between train predictions & obs skew
    test_adj = data.frame(approx(tbl_ordered$.pred_ordered, tbl_ordered$logconc_ordered, xout=test.pred$.pred)) 
    
    return(test_adj$y)
    
  }
  
  # perform EDM adjustment
  train.pred.filtered = map2(final_model, df_train, train_predict) # predict on train datasets using models trained on only the same train datasets
  
  grid.pred.adj = map2(grid.predictions, train.pred.filtered, EDM_transform) # EDM transformation
  df.grid.pred.adj = imap_dfc(grid.pred.adj, ~tibble(!!paste0('logpred_', .y) := as.numeric(.x))) %>%  # create tibble with named columns
    mutate(pred_logconcadj_mean = rowMeans(across(everything()), na.rm = TRUE)) %>% # take average
    mutate(pred_concadj_mean = 10^pred_logconcadj_mean) # back transform to geometric mean in regular concentration space
  
  # add adjusted concentrations back to the grid
  grid = bind_cols(grid,df.grid.pred.adj) 
  
  # C. Save results
  write.csv(grid, paste0(metal.code, '_GridPredicts.csv'))

}

## 5. Run grid predictions ------------------------------------------------------------------------------------------------------------------------------------------
map(metal.codes, create_grid)
map(metal.codes, predict_on_grid)
