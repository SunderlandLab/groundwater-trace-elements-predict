###########################################################################
#File name:Regression_Model.R
#Author: Jennifer Sun
#Date: Feb 2023
#Purpose: Conduct cv/grid search and select/save best regression model
###########################################################################

source(here::here('R/0_helper_fct.R'))

setwd(here::here("data"))

#### Define Model Version ####
metal.codes <- c("As", "Cd", "Li", "Mn", "Sr")
metals <- c("Arsenic", "Cadmium", "Lithium", "Manganese", "Strontium")
names(metals) <- metal.codes

metal.code <- metal.codes[1]

### 1. set up and split data randomlly -----------------------------------------------------------------------------------------------------------
df <- readRDS(paste0("R_Output/", metal.code, "_df_PredictorsSelected.rds")) # local
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
  set_args(trees = tune(),
           tree_depth = tune(),
           learn_rate = tune()) %>%
  set_engine('xgboost') %>%
  set_mode('regression')


# Compile workflow for each imputed dataset
model_workflows <- map(
  train_recipes,
  ~ workflows::workflow() %>%
    workflows::add_recipe(.x) %>%
    workflows::add_model(tree_model)
)

### 3. Conduct parameter tuning with cross-validation ------------------------------------------------------------------------------------
# Set up cross-validation for each imputed dataset
cv_folds <- map(df_train, ~ vfold_cv(.x, v = 10, repeats = 1))

# # Define parameter grid
# model_grid <- expand.grid(trees = c(100, 200, 500), 
#                           tree_depth = c(6, 9, 12), 
#                           learn_rate = c(0.01, 0.03, 0.05))

# # Perform hyperparameter tuning for each imputation
# tune_results <- map2(model_workflows, cv_folds, ~ tune_grid(
#   object = .x,
#   resamples = .y,
#   grid = model_grid,
#   metrics = metric_set(rsq, rmse, mae)
# ))
# 
# ## print results
# tune.results = model.tune.results %>% collect_metrics(summmarize=TRUE) %>% print()
# tune.results.all = model.tune.results %>% collect_metrics(summarize=FALSE) %>% print()
# 
# tune.results %>% subset(.metric=='rmse') %>% print(n=(dim(tune.results)[1]/3))
# tune.results %>% subset(.metric=='rsq') %>% print(n=(dim(tune.results)[1]/3))
# 
# ## 6. select best and simplest models from grid search -------------------------------------------------------------------------------------------------------------
# autoplot(model.tune.results) 
# tune.results = model.tune.results %>% collect_metrics(summmarize=FALSE) 
# 
# # select 'best-performing' model
# model.tune.results %>% select_best(metric='rmse')
# model.tune.results %>% select_best(metric='rsq')
# 
# ## select all models within 1 SE of the best RMSE
# tbl_test = as.data.frame(tune.results %>% subset(.metric=='rmse' & n==10))
# best_model = tbl_test %>% filter(mean == min(mean)) # why did I not use "select_best(metrics='rmse' here? because of ties?)
# best_model
# good_models = tbl_test %>% filter(mean < (best_model$mean + best_model$std_err) & mean > (best_model$mean - best_model$std_err))
# good_models %>% arrange(trees, tree_depth, mtry, desc(learn_rate), desc(loss_reduction))

## hyperparameters for 'best selected' model from good_models (best judgement from previous step)
df_hp <- data.frame(metals=c('As','Mn','Sr','Li','Cd'),
                    trees=c(500, 500, 500, 100, 100),
                    tree_depth=c(12, 12, 9, 12, 6),
                    learn_rate=c(0.03, 0.01, 0.03, 0.05, 0.05),
                    mtry=c(0.6, 0.6, 0.6, 0.6, 0.6),
                    sample_size=c(0.8, 0.8, 0.6, 1.0, 0.6),
                    loss_red=c(1.0, 1.0, 0, 0.5, 5))
df_hp <- df_hp %>% column_to_rownames(., var='metals')

#### 5. update model workflow -------------------------------------------------------------------------------------
# Update model workflow with best selected hyperparameters
best.model.workflow <- model.workflow %>% 
  update_model(
    boost_tree(
      trees = df_hp[metal.code, 'trees'],
      tree_depth = df_hp[metal.code, 'tree_depth'],
      learn_rate = df_hp[metal.code, 'learn_rate'],
      mtry = df_hp[metal.code, 'mtry'],
      sample_size = df_hp[metal.code, 'sample_size'],
      loss_red = df_hp[metal.code, 'loss_red'],
      engine = 'xgboost'
    )
  )


### 6. train final models -------------------------------------------------------------------------------------------------------------
set.seed(123)
final.model = fit(best.model.workflow, df.train[df.train$is.imputed == 0, ]) # for collecting test model metrics
final.full.model <- fit(best.model.workflow, df[df$is.imputed==0,]) # for feature analyses and prediction

## local
filename = paste0("R_Output/",modelfolder.name,metal.code,"_ModelPackage_",version.number,"_",predictor.version,"_",model.run,'.RData')
save(final.model, final.full.model, df.test, df.train, model.recipe, model.workflow, file = filename) # without tuning



