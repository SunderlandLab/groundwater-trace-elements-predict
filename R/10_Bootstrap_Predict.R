######################################
#File name: Bootstrap_Predict.R
#Author: Jennifer Sun
#Date: Nov 2025
#Purpose: Bootstrap map predictions
#####################################

#### Attach Libraries and Set Working Directory ####
packages <- c('plyr','dplyr','pROC','rpart','stringr','tidyverse','tidymodels','xgboost')
lapply(packages, library, character.only=TRUE)

setwd(here::here("data"))

#### Define Model Version ####
metal.codes <- c("As", "Cd", "Li", "Mn", "Sr")
metals <- c("Arsenic", "Cadmium", "Lithium", "Manganese", "Strontium")
names(metals) <- metal.codes

### 1. Read in and format grid data -----------------------------------------------------------------------------------------------------------

grid = readRDS(paste0(metal.code, "_MapGrid.rds"))

# match data types with model input
grid$location.id = as.character(grid$location.id) 
grid$censored = as.logical(grid$location.id) 
grid$data.source = as.factor(grid$location.id) 
grid$HLR = as.factor(grid$HLR) 
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

### 2. Set up regression model -----------------------------------------------------------------------------------------------------

bootstrap_regression_model <- function(metal.code) {
  ## A. format input data -------------------------------------------------------------------------------------------------------------
  df <- readRDS(paste0("R_Output/", metal.code, "_imputed_data.rds")) %>%
    mutate(logconc = log10(censored.conc),
           is.imputed = as.integer(censored)) 
  
  # remove any rows with missing landcover data (for Cd)
  if (metal.code=='Cd') {
    df = subset(df, !is.na(VALUE_0_pct))
  }
  
  ## B. create bootstraps ---------------------------------------------------------------------------------------------------------------
  bootstraps <- df %>% group_split(.imp, .keep = TRUE) %>% # split by imputation group, map over each of these
    imap_dfr(~ {
      imp_val <- unique(.x$.imp) # extract imputation version to save alongside the bootstrap split
      rsample::bootstraps(.x, times = 20, strata = is.imputed) %>% # run bootstraps x20 for each .imp value
        mutate(.imp = imp_val)
    })
  
  ## C. select model parameters ------------------------------------------------------------------------------------------------------------
  # set hyperparameters (extract from best selected model)
  if (metal.code=='Cd') {
    filename = paste0(metal.code, "_ModelPackage_2step.RData") 
  } else {
    filename = paste0(metal.code, "_ModelPackage.RData") 
  }
  load(filename)
  
  xgb.engine <- extract_fit_engine(final_full_model[[1]])
  xgb.params = xgb.engine$params
  df_hp = list('min_n'= xgb.params$min_child_weight,
               'mtry'= xgb.params$colsample_bytree,
               'trees'= xgb.engine$niter,
               'tree_depth'=xgb.params$max_depth,
               'learn_rate'=xgb.params$eta,
               'loss_reduction'=xgb.params$gamma,
               'sample_size'=xgb.params$subsample)
  
  # model recipe 
  preprocess_recipe <- function(data) {
    recipes::recipe(logconc ~ ., data = data) %>%
      recipes::step_rm(
        date,
        conc,
        DL.missing,
        detect.limit,
        censored,
        is.imputed,
        data.source,
        lon,
        lat,
        location.id,
        well.depth,
        .imp,
        .id
      ) %>%
      recipes::step_zv(all_predictors()) %>%
      recipes::step_normalize(all_numeric_predictors()) %>%
      recipes::step_dummy(all_nominal_predictors())
  }
  
  tree_model <- boost_tree(mode='regression',stop_iter=50) %>%
    set_args(min_n=df_hp[['min_n']], 
             trees=df_hp[['trees']],
             tree_depth=df_hp[['tree_depth']],
             learn_rate=df_hp[['learn_rate']], 
             mtry=df_hp[['mtry']], 
             sample_size=df_hp[['sample_size']],
             loss_reduction=df_hp[['loss_reduction']]) %>%   
    set_engine('xgboost', seed=123, counts=FALSE) %>% 
    set_mode('regression')

  ## D. Run bootstrap regression model -------------------------------------------------------------------------------------------------------------
  
  xgboost_model_predict <- function(split, impID, bootstrapID) {
    ### a. split data from bootstraps -----------------------------------------------------------------------------------------------------------
    
    # Extract train and test sets for each imputation
    df_train <- training(split)
    df_test <- testing(split)
    
    ### b. set up model workflow -----------------------------------------------------------------------------------------------------
    
    # Apply preprocessing to each imputed dataset
    train_recipe <- preprocess_recipe(df_train)
    
    # Compile workflow for each imputed dataset
    model_workflow <- workflows::workflow() %>%
      workflows::add_recipe(train_recipe) %>% 
      workflows::add_model(tree_model)
    
    ### c. train  models -------------------------------------------------------------------------------------------------------------
    set.seed(123)
    # for Cadmium, filter out censored data
    if (metal.code == 'Cd') {
      
      # train model
      xgb_fit = fit(model_workflow, data=subset(df_train, is.imputed==0)) # train the model on filtered data 
      
      # print("For Cd, limit to detected samples")
    } else{
      xgb_fit = fit(model_workflow, data=train_data) # train the model on ALL train data
      
    }
    
    ### d. predict -----------------------------------------------------------------------------------------------------------------------
    # predict on map grid  
    predictions <- predict(xgb_fit, new_data = grid) %>%
      bind_cols(grid) %>%     
      rename(prediction = .pred) %>% 
      select(c(longitude, latitude, prediction))
    
    # predict on the test data for metrics
    test_predictions <- predict(xgb_fit, new_data = df_test) %>% # can subset by truly detected or not later 
      bind_cols(df_test) %>%
      rename(prediction = .pred) %>%
      select(c(location.id, prediction))
    
    # predict on the train data for metrics
    train_predictions <- predict(xgb_fit, new_data = df_train) %>% # in case you want this for EDM adjustment 
      bind_cols(df_train) %>%
      rename(prediction = .pred) %>%
      select(c(location.id, prediction))
    
    # create unique bootstrap ID 
    bootID = paste0(impID, '_', bootstrapID)
    
    # Print progress
    message(bootID)
    
    # Return predictions & metrics 
    return(list(impID = impID, 
                bootID = bootID, 
                grid_predictions = predictions, 
                test_predictions = test_predictions, 
                train_predictions = train_predictions))

  }
  
  # run model
  results <- pmap(list(split=bootstraps$splits, impID=bootstraps$.imp, bootstrapID = bootstraps$id), 
                  xgboost_model_predict)
  
  # save results 
  save(results, file = paste0(metal.code,"_BootModels.RData"))
  print(paste("Bootstraps for", metal.code, "saved"))
}

  
### 3. Set up detection model for cadmium ------------------------------------------------------------------------------------------------------

bootstrap_detection_model <- function(metal.code) {
  ## A. format input data -------------------------------------------------------------------------------------------------------------
  df <- readRDS(paste0("R_Output/", metal.code, "_imputed_data.rds")) %>%
    mutate(logconc = log10(censored.conc),
           is.imputed = as.integer(censored)) 
  
  # remove any rows with missing landcover data (for Cd)
  if (metal.code=='Cd') {
    df = subset(df, !is.na(VALUE_0_pct))
  }
  
  ## B. create bootstraps ---------------------------------------------------------------------------------------------------------------
  bootstraps <- df %>% group_split(.imp, .keep = TRUE) %>% # split by imputation group, map over each of these
    imap_dfr(~ {
      imp_val <- unique(.x$.imp) # extract imputation version to save alongside the bootstrap split
      rsample::bootstraps(.x, times = 20, strata = is.imputed) %>% # run bootstraps x20 for each .imp value
        mutate(.imp = imp_val)
    })
  
  ## C. select model parameters ------------------------------------------------------------------------------------------------------------
  preprocess_detect_recipe <- function(data) {
    recipes::recipe(detect ~ ., data = data) %>%
      recipes::step_rm(
        date, conc, DL.missing, detect.limit, censored, logconc, 
        data.source, lon, lat, location.id, well.depth, .imp, .id, is.imputed,
      ) %>%
      recipes::step_zv(all_predictors()) %>%
      recipes::step_normalize(all_numeric_predictors()) %>%
      recipes::step_dummy(all_nominal_predictors())
  }
  
  clf_model <- boost_tree(mode = 'classification', stop_iter = 50) %>%
    set_args(
      trees = 200,
      tree_depth = 6,
      learn_rate = 0.03,
      sample_size = 0.8,
      loss_reduction = 1
    ) %>%
    set_engine("xgboost")
  
  ## D. Run bootstrap detection model -----------------------------------------------------------------------------------------------------------------
  xgboost_detection_model <- function(split, impID, bootstrapID) {
    
    ## a. split data from bootstraps -----------------------------------------------------------------------------------------------------------
    
    # Extract train and test sets for each imputation
    df_train <- training(split)
    df_test <- testing(split)
    
    ## b. set up model-----------------------------------------------------------------------------------------------------
    
    # Apply preprocessing to each imputed dataset
    train_recipe <- preprocess_detect_recipe(df_train)
    
    model_workflow <- workflows::workflow() %>%
      workflows::add_recipe(train_recipe) %>%
      workflows::add_model(clf_model)
    
    ## c. fit detection model----------------------------------------------------------------------------------------------------------
    xgb_detect_fit = fit(model_workflow, data=df_train)
    
    ## d. predict model ------------------------------------------------------------------------------------------------------------------------
    
    # predict on map grid  
    predictions <- predict(xgb_detect_fit, new_data = grid, type='prob') %>%
      bind_cols(grid) %>%     
      rename(pred_true = .pred_TRUE) %>% 
      select(c(longitude, latitude, pred_true))
    
    # predict on the test data for metrics
    test_predictions <- predict(xgb_detect_fit, new_data = df_test, type='prob') %>% # can subset by truly detected or not later 
      bind_cols(df_test) %>%
      rename(pred_true = .pred_TRUE) %>%
      select(c(location.id, pred_true))
    
    # predict on the train data for metrics
    train_predictions <- predict(xgb_detect_fit, new_data = df_train, type='prob') %>% # in case you want this for EDM adjustment 
      bind_cols(df_train) %>%
      rename(pred_true = .pred_TRUE) %>%
      select(c(location.id, pred_true))
    
    # create unique bootstrap ID 
    bootID = paste0(impID, '_', bootstrapID)
    
    # Print progress
    message(bootID)
    
    # Return predictions & metrics 
    return(list(impID = impID, 
                bootID = bootID, 
                grid_predictions = predictions, 
                test_predictions = test_predictions, 
                train_predictions = train_predictions)) 
  }
  
  # run model
  results <- pmap(list(split=bootstraps$splits, impID=bootstraps$.imp, bootstrapID = bootstraps$id), 
                  xgboost_detection_model)
  
  # save results 
  save(results, file = paste0(metal.code,"_detects_BootModels.RData"))
  print(paste("Bootstraps for", metal.code, "saved"))
  
}


## 4. Compile and summarize regression predictions ----------------------------------------------------------------------------------------------------------------

summarize_regression_bootstraps <- function(metal.code) {
  ## a. Read in model results ------------------------------------------------------------------------------------------------------------------------------------
  load(paste0(metal.code,"_BootModels.RData"))
  
  # original model with train predictions for EDM adjustment
  if (metal.code=='Cd') {
    load(paste0(metal.code, "_ModelPackage_2step.RData")) 
  } else {
    load(paste0(metal.code, "_ModelPackage.RData"))
  }
  
  ## b. Extract and compile bootstrap predictions ------------------------------------------------------------------------------------------------------------------
  grid_list = map(results, 3) # results
  bootID_list = map_chr(results, 2) # column names
  impID_list = map_int(results, 1) # identify adjusted train predictions to use for EDM adjustment
  
  # Rename the prediction column in each dataframe
  rename_predictions <- function(df, name) {  
    df %>% rename(!!name := prediction)
  }
  grid_list2 <- map2(grid_list, bootID_list, rename_predictions)
  
  grid_predicts <- reduce(grid_list2, left_join, by = c("longitude", "latitude"))
  
  ## c. EDM adjustment --------------------------------------------------------------------------------------------------------------------------------------------------
  # calculate train predictions for 5 imputed train datasets
  train_predict <- function(model, df.train) {
    # prepare training predictions
    df.train.filtered = subset(df.train, is.imputed==0) # test predictions similarly done only on un-imputed values in script #5
    train.pred.filtered = predict(model, df.train.filtered) %>% bind_cols(df.train.filtered)
    
    return(train.pred.filtered) # also named train.pred.adj
  }
  
  # EDM transformation
  EDM_transform <- function(test.pred, train.pred.filtered) { # grid predictions (as single vector), train.pred.filtered; 
    # NOTE: changed inputs from train.pred.adj and test.pred.adj to current versions to reflect the contents of the INPUT rather than the intended contents of the OUTPUT 
    
    # order train values
    train.pred.filtered$logconc_ordered = sort(train.pred.filtered$logconc) # order concentrations
    train.pred.filtered$.pred_ordered = sort(train.pred.filtered$.pred) # order train predictions
    tbl_ordered = train.pred.filtered[,c('logconc_ordered','.pred_ordered')]
    
    # transformation: adjust test predictions based on the difference between train predictions & obs skew
    test_adj = data.frame(approx(tbl_ordered$.pred_ordered, tbl_ordered$logconc_ordered, xout=test.pred)) 
    
    return(test_adj$y)
  }
  
  # create train dataset predictions
  train.pred.filtered = map2(final_model, df_train, train_predict) # predict on train datasets using models trained on only the same train datasets
  
  # apply EDM transformation to grid (or test) datasets predicted using model trained on corresponding train dataset 
  grid.pred.adj = map2(bootID_list, impID_list, ~EDM_transform(grid_predicts[[.x]], train.pred.filtered[[.y]])) %>% # EDM transformation (five separate grid transformations)
    set_names(bootID_list) %>%
    as_tibble() 
  
  # combine adjusted results with lat/longs
  grid.pred.adj <- bind_cols(select(grid_predicts, c('longitude','latitude')), grid.pred.adj)
  
  ## d. Add detection filter results for Cd model 
  if (metal.code == 'Cd') {
    load('Cd_DetectionModel.RData')
    
    grid = grid %>% mutate(HLR = as.numeric(HLR))
    grid.detectprobs = detection_model %>% 
      map( ~ predict(.x, new_data=grid, type='prob') %>% pull(.pred_TRUE)) %>%  # model prediction
      as_tibble(.name_repair = ~ paste0("p_model", seq_along(.))) %>%  # save results in a tibble
      mutate(p_mean = rowMeans(across(everything()))) %>% # take the average of the predictions 
      mutate(detect_class = as.integer(p_mean >=0.29)) # create filter for detections
    
    # add detection filter to adjusted prediction results
    grid.results = bind_cols(grid.pred.adj, grid.detectprobs)
  } else {
    grid.results = grid.pred.adj
  }
  
  ## e. Calculate summary values and compile into a single dataframe 
  
  grid.summary = grid.results %>%  # create tibble with named columns
    rowwise() %>%
    mutate(
      pred_log_mean = mean(c_across(contains('Bootstrap')), na.rm = TRUE), 
      pred_log_p025 = as.numeric(quantile(c_across(contains('Bootstrap')), 0.025, na.rm = TRUE)),
      pred_log_p50 = as.numeric(median( c_across(contains('Bootstrap')), na.rm = TRUE)),
      pred_log_p975 = as.numeric(quantile(c_across(contains('Bootstrap')), 0.975, na.rm = TRUE)),
      pred_log_sd = apply(as.matrix(across(contains('Bootstrap'))), 1, sd, na.rm = TRUE),
      pred_log_lower = pred_log_mean - 1.96 * pred_log_sd,
      pred_log_upper = pred_log_mean + 1.96 * pred_log_sd,
    ) %>% 
    mutate(
      pred_mean = 10^pred_log_mean,
      pred_p025 = 10^pred_log_p025,
      pred_p50  = 10^pred_log_p50,
      pred_p975 = 10^pred_log_p975,
      pred_lower = 10^pred_log_lower,
      pred_upper = 10^pred_log_upper
    ) %>%
    ungroup() %>% 
    mutate(across(where(is.numeric), ~ replace_na(., -1))) %>% # repalce all NA values with -1 so they can be mapped 
    select(-contains("Bootstrap"), -starts_with("p_model")) # remove individual model results
  
  ## D. Save results and summary grids
  write.csv(grid.results, paste0('Cd_BootGridPredicts.csv'))
  write.csv(grid.summary, paste0('Cd_BootGridPredictSummary.csv'))
 
}

### 5. Compile and summarize detection predictions -----------------------------------------------------------------------------------------------------------

summarize_detection_bootstraps <- function(metal.code) {
  results = load(paste0(metal.code, "_detects_BootModels.RData"))
  
  grid_list = map(results, 3) # results
  bootID_list = map_chr(results, 2) # column names
  impID_list = map_int(results, 1) # identify adjusted train predictions to use for EDM adjustment
  
  # Rename the prediction column in each dataframe
  rename_predictions <- function(df, name) {  
    df %>% rename(!!name := pred_true)
  }
  grid_list2 <- map2(grid_list, bootID_list, rename_predictions)
  
  grid_predicts <- reduce(grid_list2, left_join, by = c("longitude", "latitude"))
  
  grid.detects.summary = grid_predicts %>%  # create tibble with named columns
    rowwise() %>%
    mutate(
      pred_mean = mean(c_across(contains('Bootstrap')), na.rm = TRUE), 
      pred_p025 = as.numeric(quantile(c_across(contains('Bootstrap')), 0.025, na.rm = TRUE)),
      pred_p50 = as.numeric(median( c_across(contains('Bootstrap')), na.rm = TRUE)),
      pred_p975 = as.numeric(quantile(c_across(contains('Bootstrap')), 0.975, na.rm = TRUE)),
      pred_sd = apply(as.matrix(across(contains('Bootstrap'))), 1, sd, na.rm = TRUE),
      pred_lower = pred_mean - 1.96 * pred_sd,
      pred_upper = pred_mean + 1.96 * pred_sd,
    ) %>% 
    ungroup() %>% 
    mutate(across(where(is.numeric), ~ replace_na(., -1))) # repalce all NA values with -1 so they can be mapped 
  
  write.csv(grid.detects.summary, paste0(metal.code, '_BootGridDetects.csv')) 
  
}



### 6. Formally run and save models --------------------------------------------------------------------------------------------------------------

# run bootstrap regression models 
map(metal.codes, bootstrap_regression_model)
map(metal.codes, summarize_regression_bootstraps)

# run detection model for Cd 
map(metal.codes, bootstrap_regression_model)
map(metal.codes, summarize_detection_bootstraps)
