
source(here::here('R/0_helper_fct.R'))
plan(multisession, workers = future::availableCores())
print(paste("The number of workers are", nbrOfWorkers()))
options(future.rng.onMisuse = "ignore")

setwd(here::here("data"))

#### Define Model Version ####
metal.codes <- c("As", "Cd", "Li", "Mn", "Sr")
metals <- c("Arsenic", "Cadmium", "Lithium", "Manganese", "Strontium")
names(metals) <- metal.codes
metal.code <- 'Cd'
### 1. set up and split data randomlly -----------------------------------------------------------------------------------------------------------
df <- readRDS(paste0("Data_Files/", metal.code, "_df_PredictorsSelected.rds")) 
# remove any rows with missing landcover data (for Cd)
if (metal.code=='Cd') {
  df = subset(df, !is.na(VALUE_0_pct))
}

df$logconc = log10(df$censored.conc) 
df$is.imputed = as.integer(df$censored)

set.seed(123)
df.split = rsample::initial_split(df, prop=4/5, strata=is.imputed) # stratify by whether or not the data are imputed
df.train = rsample::training(df.split)
df.test = rsample::testing(df.split)

### 2. set up model-----------------------------------------------------------------------------------------------------
model.recipe = recipes::recipe(logconc~., data=df) %>%
  recipes::step_rm(date, ros.conc, censored.conc, conc, DL.missing, censored,is.imputed, data.source, long, lat, location.id) %>%
  recipes::step_rm(well.depth) %>%
  recipes::step_zv(all_predictors()) %>% 
  recipes::step_normalize(all_numeric_predictors()) %>% # includes numeric and integer types; this is actually standardization (mean=0, std=1) 
  recipes::step_dummy(all_nominal_predictors()) 

# with grid search parameters
tree.model <- boost_tree(mode='regression',stop_iter=50) %>%
  set_args(trees=tune(), tree_depth=tune(),learn_rate=tune()) %>%
  set_engine('xgboost') %>%
  set_mode('regression')


## compile model workflow
model.workflow <- workflows::workflow() %>%
  workflows::add_recipe(model.recipe) %>%
  workflows::add_model(tree.model)

## hyperparameters for 'best selected' model from good_models (best judgement from previous step)
df_hp <- data.frame(metals=c('As','Mn','Sr','Li','Cd'),
                    trees=c(500, 500, 500, 100, 100),
                    tree_depth=c(12, 12, 9, 12, 6),
                    learn_rate=c(0.03, 0.01, 0.03, 0.05, 0.05),
                    mtry=c(0.6, 0.6, 0.6, 0.6, 0.6),
                    sample_size=c(0.8, 0.8, 0.6, 1.0, 0.6),
                    loss_reduction=c(1.0, 1.0, 0, 0.5, 5))
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
      loss_reduction = df_hp[metal.code, 'loss_reduction']
    )%>%
      set_engine('xgboost', counts = FALSE) %>%
      set_mode('regression')
  )


### 6. train final models -------------------------------------------------------------------------------------------------------------
set.seed(123)
final.model = fit(best.model.workflow, df.train[df.train$is.imputed == 0, ]) # for collecting test model metrics
final.full.model <- fit(best.model.workflow, df[df$is.imputed==0,]) # for feature analyses and prediction
# # Split data separately for each imputation
# df_splits <- df %>%
#   group_split(.imp) %>%  # Split into separate imputed datasets
#   future_map( ~ rsample::initial_split(.x, prop = 4 / 5))  # Stratify by is.imputed
# 
# # Extract train and test sets for each imputation
# df_train <- future_map(df_splits, rsample::training)
# df_test <- future_map(df_splits, rsample::testing)
# 
# ### 2. set up model-----------------------------------------------------------------------------------------------------
# preprocess_recipe <- function(data) {
#   recipes::recipe(logconc ~ ., data = data) %>%
#     recipes::step_rm(
#       date,
#       conc,
#       DL.missing,
#       censored,
#       censored.conc,
#       ros.conc,
#       data.source,
#       long,
#       lat,
#       location.id,
#       well.depth,
#       .imp
#     ) %>%
#     recipes::step_zv(all_predictors()) %>%
#     recipes::step_normalize(all_numeric_predictors()) %>%
#     recipes::step_dummy(all_nominal_predictors())
# }
# 
# # Apply preprocessing to each imputed dataset
# train_recipes <- future_map(df_train, preprocess_recipe)
# 
# # with grid search parameters
# tree_model <- boost_tree(mode = 'regression', stop_iter = 50) %>%
#   set_args(
#     trees = tune(),
#     tree_depth = tune(),
#     learn_rate = tune(),
#     sample_size = tune(),
#     loss_reduction = tune()
#   ) %>%
#   set_engine('xgboost', counts = FALSE) %>%
#   parsnip::set_mode('regression')
# 
# 
# # Compile workflow for each imputed dataset
# model_workflows <- map(
#   train_recipes,
#   ~ workflows::workflow() %>%
#     workflows::add_recipe(.x) %>%
#     workflows::add_model(tree_model)
# )
# 
# # Get the names of valid arguments for boost_tree()
# valid_args <- names(formals(boost_tree))
# 
# # Subset only valid columns
# params <- df_hp %>% filter(metal == metal.code) %>%
#   select(intersect(names(df_hp), valid_args)) %>%
#   as.list()
# # Update model workflow with best selected hyperparameters
# best_model_workflow <-  purrr::map(
#   model_workflows,
#   ~ update_model(.x, do.call(boost_tree, c(params, list(mode = "regression", engine = "xgboost"))))
# )
# 
# 
# ### 6. train final models -------------------------------------------------------------------------------------------------------------
# set.seed(123)
# final_model <- map2(best_model_workflow, df_train, fit) # for collecting test model metrics
# final_full_model <- map2(best_model_workflow, df %>%group_split(.imp), fit) # for feature analyses and prediction
# 
#load model project
#filename <- paste0('R_Output/Cd_ModelPackage_2a_v1_xgboost_tune2_std_simple.RData')
#load(filename)
final_model <- final.model
df_train <- df.train
df_test <- df.test
train_recipes <- model.recipe
### 7. begine SHAP analysis
#### 1. Prepare data --------------------------------------------------------------------------------------------------------
## extract xgboost object from model trained on training data, for diagnostics
final_xg_list <- extract_fit_engine(final_model)
# Initialize output lists
df_pred_list <- list()
factorcols_list <- list()

# Loop over imputed datasets

  df <- rbind(df_train, df_test)

  # Process predictors using recipe
  df_pred <- bake(
    prep(train_recipes), # bake with the reference recipe to ensure consistency
    has_role('predictor'),
    new_data = df,
    composition = 'matrix'
  )

  # Extract factor names
  factors <- df %>%
    dplyr::select(-c('location.id', 'data.source')) %>%
    dplyr::select_if(is.factor) %>%
    colnames()

  # Extract factor-related columns
  factorcols <- list()
  for (factor in factors) {
    cols <- as.data.frame(df_pred) %>%
      dplyr::select(starts_with(factor)) %>%
      colnames()
    factorcols[[factor]] <- cols
  }

  # Align df_pred with model features
  remove_features <- c(setdiff(colnames(df_pred), final_xg_list$feature_names))
  print(paste0("features to remove from df_pred are ", remove_features))
  remove_index <- which(colnames(df_pred) %in% remove_features)
  if (length(remove_index) > 0) {
    df_pred <- df_pred[, -remove_index]
  }

  # Save results
  df_pred_list <- df_pred
  factorcols_list <- factorcols

#### 2a. Calculate xgboost shapley values --------------------------------------------------------------------------------------------------------

  
    df <- rbind(df_train, df_test)
    df_pred <- df_pred_list
    factorcols <- factorcols_list

    shap_list <- shapviz(
      object = final_xg_list,
      X_pred = df_pred,
      X = df,
      collapse = factorcols
    )


# 1. Extract SHAP matrices
shap_matrices <- shap_list$S
# 2. Pool SHAP values by averaging across imputations
# Find common column names
common_cols <- Reduce(intersect, lapply(shap_matrices, colnames))

# # Subset each matrix to only those common columns
# shap_matrices_common <- lapply(shap_matrices, function(mat) mat[, common_cols, drop = FALSE])
# 
# # Average them
# pooled_S <- Reduce(`+`, shap_matrices_common) / length(shap_matrices_common)
# 3. Pool baseline values
baseline_values <- shap_list$baseline
pooled_baseline <- mean(baseline_values)
# 4. Use the feature matrix from the first object
X_ref <- shap_list$X
# 5. Reconstruct a new shapviz object
shap_pooled <- shapviz(
  shap_matrices,
  X = X_ref,
  baseline_value = pooled_baseline,
  model_class = "xgboost"
)
varImp.plot <- sv_importance(shap_pooled, kind = "both", show_numbers = TRUE) +
  scale_y_discrete(limits=rev(calc_meanabs_shap(shap_pooled))) +
  theme_classic() +
  theme(axis.text.x = element_text(size=12), axis.text.y=element_text(size=12))
varImp.plot
ggsave(paste0('R_Output/',metal.code,'_SHAP_importance_plot_jennifer_rep.png'), width = 8, height = 6)
