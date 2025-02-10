###########################################################################
#File name:Regression_Model.R
#Author: Jennifer Sun, Cindy Hu
#Date: Jan 2025
#Purpose: Conduct cv/grid search and select/save best regression model
###########################################################################

source(here::here('R/0_helper_fct.R'))

setwd(here::here("data"))

#### Define Model Version ####
metal.codes <- c("As", "Cd", "Li", "Mn", "Sr")
metals <- c("Arsenic", "Cadmium", "Lithium", "Manganese", "Strontium")
names(metals) <- metal.codes

tune_xgboost_models <- function(metal.code) {
  print(paste("Start hyperparameter tuning for", metal.code))
  ### 1. set up and split data randomlly -----------------------------------------------------------------------------------------------------------
  df <- readRDS(paste0("R_Output/", metal.code, "_df_PredictorsSelected.rds")) #%>%
  # group_by(.imp) %>%
  # sample_n(size = 100, replace = FALSE) %>%
  # ungroup()
  df$logconc = log10(df$conc)
  df$is.imputed = as.integer(df$censored)
  
  set.seed(123)
  # Split data separately for each imputation
  df_splits <- df %>%
    group_split(.imp) %>%  # Split into separate imputed datasets
    map( ~ rsample::initial_split(.x, prop = 4 / 5, strata = is.imputed))  # Stratify by is.imputed
  
  # Extract train and test sets for each imputation
  df_train <- map(df_splits, rsample::training)
  df_test <- map(df_splits, rsample::testing)
  
  ### 2. set up model-----------------------------------------------------------------------------------------------------
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
  
  # Apply preprocessing to each imputed dataset
  train_recipes <- map(df_train, preprocess_recipe)
  
  # with grid search parameters
  tree_model <- boost_tree(mode = 'regression', stop_iter = 50) %>%
    set_args(
      trees = tune(),
      tree_depth = tune(),
      learn_rate = tune(),
      sample_size = tune(),
      loss_reduction = tune()
    ) %>%
    set_engine('xgboost') %>%
    parsnip::set_mode('regression')
  
  
  # Compile workflow for each imputed dataset
  model_workflows <- map(
    train_recipes,
    ~ workflows::workflow() %>%
      workflows::add_recipe(.x) %>%
      workflows::add_model(tree_model)
  )
  
  #### 3. Conduct parameter tuning with cross-validation ------------------------------------------------------------------------------------
  # Set up cross-validation for each imputed dataset
  cv_folds <- map(df_train, ~ vfold_cv(.x, v = 10, repeats = 1))
  
  # # Define parameter grid
  model_grid <- expand.grid(
    trees = c(100, 200, 500),
    tree_depth = c(6, 9, 12),
    learn_rate = c(0.01, 0.03, 0.05),
    sample_size = c(0.6, 0.8, 1.0),
    loss_reduction = c(0, 0.5, 1, 5)
  )
  
  # Perform hyperparameter tuning for each imputation
  system.time({
    tune_results <- map2(
      model_workflows,
      cv_folds,
      ~ tune_grid(
        object = .x,
        resamples = .y,
        grid = model_grid,
        metrics = metric_set(rsq, rmse, mae)
      )
    )
  })
  
  # Collect tuning results for all imputations
  all_tune_results <- map_dfr(tune_results, collect_metrics, .id = "imputation") %>%
    # average across five imputations to find the best hyperparameters
    group_by(.metric, .config) %>%
    summarise(
      mean_val = mean(mean),
      # pooled std_err, same n, sum of squares divided by the number of configurations, and then square root
      mean_p_sd = sqrt(sum(std_err ^ 2) / n()),
      # same configuration has the same hyperparameters
      across(
        c(trees, tree_depth, learn_rate, sample_size, loss_reduction, n),
        ~ first(.x)
      ),
      .groups = 'keep'
    )
  
  tune_summary <- all_tune_results %>%
    group_by(.metric) %>%
    summarise(min_val = min(mean_val),
              pooled_sd = sqrt(sum(mean_p_sd ^ 2) / n())) %>%
    mutate(min_plus_sd = min_val + pooled_sd)
  
  ## 6. select best and simplest models from grid search -------------------------------------------------------------------------------------------------------------
  ## select all models within 1 SE of the best RMSE
  good_models <- all_tune_results %>%
    filter(.metric == 'rmse') %>%
    left_join(tune_summary %>% select(.metric, min_val, min_plus_sd),
              by = c('.metric')) %>%
    filter(mean_val < min_plus_sd, mean_val > min_val)
  
  if (nrow(good_models) > 1) {
    good_models <- good_models %>%
      arrange(trees, tree_depth, sample_size, desc(learn_rate)) %>%
      filter(row_number() == 1)
  }
  
  
  write_csv(good_models,
            paste0("R_Output/", metal.code, "_df_hyperparameter.csv"))
  print(paste0('Variable selection for ', metal.code, ' complete'))
  
}

purrr::map(metal.codes, tune_xgboost_models)

# xgboost_without_tuning <- function(i){
#   # Define the XGBoost model with fixed hyperparameters
#   xgb_spec <- boost_tree(
#     trees = 500,           # Fixed number of trees
#     tree_depth = 9,        # Fixed tree depth
#     learn_rate = 0.03,     # Fixed learning rate
#     loss_reduction = 0, # Minimum loss reduction
#     sample_size = 0.6,     # Row sampling
#     mtry = floor(0.6 * ncol(df_train[[i]])-1),             # Feature sampling
#     mode = "regression"    # Change to "classification" if needed
#   ) %>%
#     set_engine("xgboost")
#
#   # Create a workflow
#   xgb_wf <- workflow() %>%
#     add_recipe(train_recipes[[i]]) %>%
#     add_model(xgb_spec)
#
#   # Fit the model to the training data
#   xgb_fit <- fit(xgb_wf, data = df_train[[i]])
#
#   # Make predictions on test data
#   predictions <- predict(xgb_fit, df_test[[i]]) %>%
#     bind_cols(df_test[[i]])
#
#   # Evaluate model performance
#   metrics <- predictions %>%
#     metrics(truth = logconc, estimate = .pred)
#
#   # Print evaluation metrics
#   print(metrics)
#
#   # Plot variable importance
#   xgb_fit %>%
#     extract_fit_parsnip() %>%
#     vip::vip(num_features = 10) +
#     ggtitle(paste0("Variable Importance for ", metal.code, " Prediction, Imputation [", i, ']'))
#
# }
#
# purrr::map(1:5, xgboost_without_tuning)
# ## hyperparameters for 'best selected' model from good_models (best judgement from previous step)
# df_hp <- data.frame(metals=c('As','Mn','Sr','Li','Cd'),
#                     trees=c(500, 500, 500, 100, 100),
#                     tree_depth=c(12, 12, 9, 12, 6),
#                     learn_rate=c(0.03, 0.01, 0.03, 0.05, 0.05),
#                     mtry=c(0.6, 0.6, 0.6, 0.6, 0.6),
#                     sample_size=c(0.8, 0.8, 0.6, 1.0, 0.6),
#                     loss_reduction=c(1.0, 1.0, 0, 0.5, 5))
# df_hp <- df_hp %>% column_to_rownames(., var='metals')
#
# #### 5. update model workflow -------------------------------------------------------------------------------------
# # Update model workflow with best selected hyperparameters
# best_model_workflow <- model_workflows[[1]] %>%
#   update_model(
#     boost_tree(
#       trees = df_hp[metal.code, 'trees'],
#       tree_depth = df_hp[metal.code, 'tree_depth'],
#       learn_rate = df_hp[metal.code, 'learn_rate'],
#       mtry = df_hp[metal.code, 'mtry'],
#       sample_size = df_hp[metal.code, 'sample_size'],
#       loss_reduction = df_hp[metal.code, 'loss_reduction'],
#       engine = 'xgboost'
#     )
#   )

#
# ### 6. train final models -------------------------------------------------------------------------------------------------------------
# set.seed(123)
# final.model = fit(best.model.workflow, df.train[df.train$is.imputed == 0, ]) # for collecting test model metrics
# final.full.model <- fit(best.model.workflow, df[df$is.imputed==0,]) # for feature analyses and prediction
#
# ## local
# filename = paste0("R_Output/",modelfolder.name,metal.code,"_ModelPackage_",version.number,"_",predictor.version,"_",model.run,'.RData')
# save(final.model, final.full.model, df.test, df.train, model.recipe, model.workflow, file = filename) # without tuning
#
#
#
