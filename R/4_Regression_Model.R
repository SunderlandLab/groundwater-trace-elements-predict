###########################################################################
#File name:Regression_Model.R
#Author: Jennifer Sun
#Date: Feb 2023
#Purpose: Conduct cv/grid search and select/save best regression model
###########################################################################

#### Attach Libraries and Set Working Directory ####
packages <- c('beepr','plyr','dplyr','pROC','rpart','stringr','viridis','tidyverse','tidymodels','fastDummies','doParallel')
lapply(packages, library, character.only=TRUE)

setwd('/Users/jennifer/Documents/Harvard/Drinking Water/Metals_Modeling.nosync/DWheavymetal')

#### Define Global Variables ####
metal = 'Manganese'
MCL = 300
metal.code = "Mn"
version.number = '7a'
predictor.version = 'v1'
folder.name = paste0(metal.code,'_',version.number,'/')


### 1a. set up and split data randomlly -----------------------------------------------------------------------------------------------------------
df = readRDS(paste0("R_Output/",folder.name,metal.code,"_df_PredictorsSelected_",version.number,"_",predictor.version,".rds")) # local

# remove any rows with missing landcover data (for Cd)
if (metal=='Cadmium') {
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
  recipes::step_dummy(all_nominal_predictors()) %>%
  recipes::prep()


# with grid search parameters
tree.model <- boost_tree(mode='regression',stop_iter=50) %>%
  set_args(trees=tune(), tree_depth=tune(),learn_rate=tune()) %>%
  set_engine('xgboost') %>%
  set_mode('regression')


## compile model workflow
model.workflow <- workflows::workflow() %>%
  workflows::add_recipe(model.recipe) %>%
  workflows::add_model(tree.model)

### 3. conduct parameter tuning grid search with CV----------- ---------------------------------------------------------------------------
model.cv = group_vfold_cv(df.train, folds, v=10, repeats=1)

model.grid = expand.grid(trees=c(100, 150, 200), tree_depth=c(6,12,18), learn_rate=c(0.1, 0.05, 0.01)) # make sure this is even working

model.tune.results <- model.workflow %>%
  tune::tune_grid(resamples=model.cv, grid=model.grid, metrics=yardstick::metric_set(rsq, rmse, mae))

## print results
tune.results = model.tune.results %>% collect_metrics(summmarize=TRUE) %>% print()
tune.results.all = model.tune.results %>% collect_metrics(summarize=FALSE) %>% print()

tune.results %>% subset(.metric=='rmse') %>% print(n=(dim(tune.results)[1]/3))
tune.results %>% subset(.metric=='rsq') %>% print(n=(dim(tune.results)[1]/3))

## 6. select best and simplest models from grid search -------------------------------------------------------------------------------------------------------------
autoplot(model.tune.results) 
tune.results = model.tune.results %>% collect_metrics(summmarize=FALSE) 

# select 'best-performing' model
model.tune.results %>% select_best(metric='rmse')
model.tune.results %>% select_best(metric='rsq')

## select all models within 1 SE of the best RMSE
tbl_test = as.data.frame(tune.results %>% subset(.metric=='rmse' & n==10))
best_model = tbl_test %>% filter(mean == min(mean)) # why did I not use "select_best(metrics='rmse' here? because of ties?)
best_model
good_models = tbl_test %>% filter(mean < (best_model$mean + best_model$std_err) & mean > (best_model$mean - best_model$std_err))
good_models %>% arrange(trees, tree_depth, mtry, desc(learn_rate), desc(loss_reduction))

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



