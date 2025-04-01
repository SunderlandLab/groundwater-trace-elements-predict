###########################################################################
#File name:Regression_Model.R
#Author: Jennifer Sun, Cindy Hu
#Date: Feb 26, 2025
#Purpose: Conduct cv/grid search and select/save best regression model
###########################################################################

source(here::here('R/0_helper_fct.R'))
plan(multisession, workers = future::availableCores())
print(paste("The number of workers are", nbrOfWorkers()))
options(future.rng.onMisuse = "ignore")

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
    future_map( ~ rsample::initial_split(.x, prop = 4 / 5, strata = is.imputed))  # Stratify by is.imputed
  
  # Extract train and test sets for each imputation
  df_train <- future_map(df_splits, rsample::training)
  df_test <- future_map(df_splits, rsample::testing)
  
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
  train_recipes <- future_map(df_train, preprocess_recipe)
  
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
  cv_folds <- future_map(df_train, ~ vfold_cv(.x, v = 10, repeats = 1))
  
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
    tune_results <- future_map2(
      model_workflows,
      cv_folds,
      ~ tune_grid(
        object = .x,
        resamples = .y,
        grid = model_grid,
        metrics = metric_set(rsq, rmse, mae),
        control = control_grid(parallel_over = "everything")
      ),
      seed = NULL
    )
  })
  
  # Collect tuning results for all imputations
  all_tune_results <- future_map_dfr(tune_results, collect_metrics, .id = "imputation") %>%
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
  
  ## 4. select best and simplest models from grid search -------------------------------------------------------------------------------------------------------------
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

#keep this as purrr because furrr likes to error out, won't matter if we submit one metal at a time
purrr::map(metal.codes, tune_xgboost_models)

# The hyperparameter tuning was run on Odyssey HPC from February 28 to March 24, 2025; each trace element took about one week.

# Keep the lowest RMSE model for each metal
get_best_model <- function(metal.code) {
  read_csv(paste0("R_Output/", metal.code, "_df_hyperparameter.csv")) %>%
    filter(mean_val == min(mean_val)) %>%
    mutate(metal = metal.code)
}

df_hp <- purrr::map(metal.codes, get_best_model) %>%
  bind_rows()

#### 5. update model workflow -------------------------------------------------------------------------------------
update_xgboost_models <- function(metal.code) {
  print(paste("Update model workflow for", metal.code))
  ### 1. set up and split data randomlly -----------------------------------------------------------------------------------------------------------
  df <- readRDS(paste0("R_Output/", metal.code, "_df_PredictorsSelected.rds"))
  df$logconc = log10(df$conc)
  df$is.imputed = as.integer(df$censored)
  
  set.seed(123)
  # Split data separately for each imputation
  df_splits <- df %>%
    group_split(.imp) %>%  # Split into separate imputed datasets
    future_map( ~ rsample::initial_split(.x, prop = 4 / 5, strata = is.imputed))  # Stratify by is.imputed
  
  # Extract train and test sets for each imputation
  df_train <- future_map(df_splits, rsample::training)
  df_test <- future_map(df_splits, rsample::testing)
  
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
  train_recipes <- future_map(df_train, preprocess_recipe)
  
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
  
  # Get the names of valid arguments for boost_tree()
  valid_args <- names(formals(boost_tree))
  
  # Subset only valid columns
  params <- df_hp %>% filter(metal == metal.code) %>%
    select(intersect(names(df_hp), valid_args)) %>%
    as.list()
  # Update model workflow with best selected hyperparameters
  best_model_workflow <-  purrr::map(
    model_workflows,
    ~ update_model(.x, do.call(boost_tree, c(params, list(mode = "regression", engine = "xgboost"))))
  )


  ### 6. train final models -------------------------------------------------------------------------------------------------------------
  set.seed(123)
  final_model <- map2(best_model_workflow, df_train, fit) # for collecting test model metrics
  final_full_model <- map2(best_model_workflow, df %>%group_split(.imp), fit) # for feature analyses and prediction

  filename <- paste0("R_Output/", metal.code,"_ModelPackage.RData")
  save(final_model, final_full_model, df_test, df_train, train_recipes, model_workflows, file = filename) 
  print(paste("Model package for", metal.code, "saved"))
}

# iterate over five trace elements and save model packages
purrr::map(metal.codes, update_xgboost_models)