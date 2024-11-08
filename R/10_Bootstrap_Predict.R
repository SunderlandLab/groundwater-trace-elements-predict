######################################
#File name:Bootstrap_Grid.R
#Author: Jennifer Sun
#Date: July 2024
#Purpose: Bootstrap map predictions
#####################################

#### Attach Libraries and Set Working Directory ####
packages <- c('plyr','dplyr','pROC','rpart','stringr','tidyverse','tidymodels','xgboost')
lapply(packages, library, character.only=TRUE)

setwd(here::here("data"))

metal = 'Arsenic'
metal.code = "As"
version.number = '7a'
predictor.version = 'v1'
grid.version = 'tune8_std'
set.seed(123)

### 1. Read in and format data -----------------------------------------------------------------------------------------------------------

## A. training data 
df = readRDS(paste0(metal.code,"_df_PredictorsSelected_",version.number,"_",predictor.version,".rds")) 

# remove any rows with missing landcover data (for Cd)
if (metal=='Cadmium') {
  df = subset(df, !is.na(VALUE_0_pct))
}

df$logconc = log10(df$censored.conc) 
df$is.imputed = as.integer(df$censored)


## B. grid data 
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

### 2. set up model-----------------------------------------------------------------------------------------------------

## select hyperparameters (best selected model)
df_hp <- data.frame(metals=c('As','Mn','Sr','Li','Cd'),
                    trees=c(500, 500, 500, 100, 100),
                    tree_depth=c(12, 12, 9, 12, 6),
                    learn_rate=c(0.03, 0.01, 0.03, 0.05, 0.05),
                    mtry=c(0.6, 0.6, 0.6, 0.6, 0.6),
                    sample_size=c(0.8, 0.8, 0.6, 1.0, 0.6),
                    loss_red=c(1.0, 1.0, 0, 0.5, 5))
df_hp <- df_hp %>% column_to_rownames(., var='metals')

## regression model
model.recipe = recipes::recipe(logconc~., data=df) %>%
  recipes::step_rm(date, ros.conc, censored.conc, conc, DL.missing, censored,is.imputed, data.source, long, lat, location.id) %>%
  recipes::step_rm(well.depth) %>%
  recipes::step_rm(detect.limit) %>%
  recipes::step_zv(all_predictors()) %>% # added 8/4/23
  recipes::step_normalize(all_numeric_predictors()) %>% # includes numeric and integer types; this is actually standardization (mean=0, std=1) 
  recipes::step_dummy(all_nominal_predictors()) %>%
  recipes::prep()

## set model parameters
tree.model <- boost_tree(mode='regression',stop_iter=50, min_n=10, trees=df_hp[metal.code, 'trees'],tree_depth=df_hp[metal.code, 'tree_depth'],
                         learn_rate=df_hp[metal.code, 'learn_rate'], mtry=df_hp[metal.code, 'mtry'], sample_size=df_hp[metal.code, 'sample_size'],
                         loss_reduction=df_hp[metal.code, 'loss_red']) %>% #  
  set_engine('xgboost', counts=FALSE) %>% # objective='reg:tweedie'
  set_mode('regression')

## compile model workflow
model.workflow <- workflows::workflow() %>%
  workflows::add_recipe(model.recipe) %>%
  workflows::add_model(tree.model) 

### 3. Run bootstrap model  -------------------------------------------------------------------------------------------------------------

## bootstrap resample
df.folds = bootstraps(data = df, times = 100, strata = is.imputed)

train_predict_metrics <- function(split, counter) {
  train_data <- analysis(split) 
  test_data <- assessment(split)
  
  # Train model
  xgb_fit = fit(model.workflow, data=subset(train_data, is.imputed==0))
  
  # Predict on the map grid
  predictions <- predict(xgb_fit, new_data = grid) %>%
    bind_cols(grid) %>%     
    rename(prediction = .pred) %>% 
    select(c(longitude, latitude, prediction))
  
  # Predict on the test data for metrics
  test_predictions <- predict(xgb_fit, new_data = test_data) %>%
    bind_cols(test_data) %>%
    rename(prediction = .pred)
  
  # Calculate performance metrics
  metrics <- metrics(test_predictions, truth = logconc, estimate = prediction) # what is target? 
  
  # Print progress
  print(counter)
  
  # Return predictions & metrics 
  return(list(predictions = predictions, metrics = metrics))
}

# Run the models and gather predictions and metrics
results <- map2(df.folds$splits, seq_along(df.folds$splits), train_predict_metrics)


## 4. Compile predictions ----------------------------------------------------------------------------------------------------------------
rename_predictions <- function(df, name) {
  df %>% rename(!!name := prediction)
}

df_list = map(results, 1) # extract the list of predictions from results
df_list <- map(df_list, as.data.frame) # extract the dataframes from the tibbles

# Rename the prediction column in each dataframe
df_list2 <- imap(df_list, ~ rename_predictions(.x, paste0("model", .y)))

# Merge all dataframes using reduce and left_join
merged_df <- reduce(df_list2, left_join, by = c("longitude", "latitude"))


### 5. EDM transformations ----------------------------------------------------------------------------------------------------------------
EDM_transform <- function(test.pred.adj, train.pred.adj) {  
  # order train values
  train.pred.adj$logconc_ordered = sort(train.pred.adj$logconc) # order concentrations
  train.pred.adj$.pred_ordered = sort(train.pred.adj$.pred) # order predictions
  tbl_ordered = train.pred.adj[,c('logconc_ordered','.pred_ordered')]
  
  # transformation: adjust test predictions based on the difference between train predictions & obs skew
  test_adj = data.frame(approx(tbl_ordered$.pred_ordered, tbl_ordered$logconc_ordered, xout=test.pred.adj)) 

  return(test_adj$y)
}


# EDM adjustment
adj_df = merged_df[,-c(1,2)] %>%
  map_dfc(~ EDM_transform(.x, train.predictions.uncens))

grid.predictions = bind_cols(merged_df[,c(1,2)], adj_df)

# Save all prediction values
write.csv(grid.predictions, paste0(folder.name, bootstraps.name, metal.code,"_BootGridPredictsAll_",version.number,"_",predictor.version,'.csv'))

### 6. Calculate summary values ----------------------------------------------------------------------------------------------------------------
# Calculate the average
pred_means <- rowMeans(adj_df, na.rm = TRUE) # in a limited number of cases, values may be outside model bounds and results in NA predictions
pred_sd <- apply(adj_df, 1, sd, na.rm=TRUE)
lower_bound <- pred_means - 1.96 * pred_sd
upper_bound <- pred_means + 1.96 * pred_sd
pred_median <- apply(adj_df, 1, median, na.rm=TRUE)
na_count <- rowSums(is.na(adj_df))

summary_df = data.frame(
  mean = 10**pred_means,
  sd = 10**pred_sd,
  median = 10**pred_median,
  lower = 10**lower_bound,
  upper = 10**upper_bound, 
  NAs = na_count
)

# Replace all NA values with -1 so they can be mapped, but they are still recorded separately from predicted 0s
summary_df[is.na(summary_df)] = -1

# Compile a dataframe
grid.summary = bind_cols(merged_df[,c(1,2)], summary_df)

# Save bootstrap summary values
write.csv(grid.summary, paste0(folder.name, bootstraps.name, metal.code,"_BootGridPredictSummary",version.number,"_",predictor.version,'.csv'))
