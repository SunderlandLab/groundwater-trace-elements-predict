####################################################################################################
#File name: Model_Eval_Reg.R
#Author: Jennifer Sun
#Date: Feb 2023
#Purpose: Model evaluation, including EDM correction and both regression and classification metrics
####################################################################################################

#### Attach Libraries and Set Working Directory ####
packages <- c('dplyr','ggplot2','plyr','pROC','ranger','RColorBrewer','stringr', 'VIM', 'pdp',
              'viridis','formattable','tidyverse','tidymodels','RFpredInterval')
lapply(packages, library, character.only=TRUE)

setwd('/Users/jennifer/Documents/Harvard/Drinking Water/Metals_Modeling.nosync/DWheavymetal')

### Read in full data package
metal = 'Arsenic'
metal.code = 'As'
MCL = 10
version.number = '7a' # 5a
grid.version = 'tune8_std'
predictor.version = 'v1'
model.run = 'xgboost_tune8_std_simple' # 'xgboost_tune8_std_simple
cluster.folder = 'tune8_std/' # 'scv_50km_5folds/scv_50km_5folds_tune2/'
folder.name = paste0(metal.code,'_',version.number,'/')
modelfolder.name = paste0(metal.code,'_',version.number,'/',grid.version,'/')

filename = paste0('R_Output/',modelfolder.name,metal.code,"_ModelPackage_",version.number,"_",predictor.version,"_",model.run,'.RData' )
load(filename)


#### Helper functions  -------------------------------------------------------------------------------------------------------

## print result metrics
calc_reg_metrics <- function(df_results) {
  logregMetrics = df_results %>% yardstick::metrics(logconc, .pred) %>% as.data.frame() %>%
                  dplyr::select(-c(.estimator)) %>% dplyr::rename('metric'='.metric','estimate'='.estimate')
  return(logregMetrics)
}

calc_class_metrics <- function(df_results, threshold) {
  multi_met = metric_set(accuracy, spec, sens)
  df_results$pred_exceeds = as.numeric(df_results$antilog_pred >= threshold)
  df_results$conc_exceeds = as.numeric(df_results$conc >= threshold)
  logclassResults = df_results %>% multi_met(truth=as.factor(conc_exceeds), estimate=as.factor(pred_exceeds), event_level='second') %>%  # first level = 0, second level = 1
                    as.data.frame() %>% dplyr::select(-c(.estimator)) %>% dplyr::rename('metric'='.metric','estimate'='.estimate')
  return(logclassResults)
}

## EDM transformation
EDM_transform <- function(test.pred.adj, train.pred.adj) {
  # order train values
  train.pred.adj$logconc_ordered = sort(train.pred.adj$logconc) # order concentrations
  train.pred.adj$.pred_ordered = sort(train.pred.adj$.pred) # order predictions
  tbl_ordered = train.pred.adj[,c('logconc_ordered','.pred_ordered')]
  
  # transformation: adjust test predictions based on the difference between train predictions & obs skew
  test_adj = data.frame(approx(tbl_ordered$.pred_ordered, tbl_ordered$logconc_ordered, xout=test.pred.adj$.pred)) 
  # warning: in regularize.value() collapsing to unique 'x' values -- may be when test data is outside the bounds of train data and has to find nearest value
  test.pred.adj$.pred_adj = test_adj$y # append transformed data to the original dataset
  
  # check that transformation was completed as expected
  test.pred.adj$.pred_check = test_adj$x # this should match .pred
  # head(test.pred.adj[c('.pred','.pred_adj','.pred_check')]) # double check that .pred matches .pred_check
  
  # order test dataset for plotting
  test.pred.adj$logconc_ordered = sort(test.pred.adj$logconc)
  test.pred.adj$.pred_ordered = sort(test.pred.adj$.pred)
  
  # rename columns so the adjusted predictions are '.pred' and can be read into other functions
  test.pred.adj2 = test.pred.adj %>% dplyr::rename(.pred_orig=.pred,.pred=.pred_adj,antilog_pred_orig=antilog_pred)
  
  # calculate predictions in regular concentration units
  test.pred.adj2$antilog_pred = 10**test.pred.adj2$.pred
  
  return(test.pred.adj2)
}

## print figures
plot_trans_ordered <- function(df.test, df.train) { # changed both of these from '.adj' to not adjusted
  # plot ordered predicted and observed values (compare distributions)
  df.train$logconc_ordered = sort(df.train$logconc) # order concentrations
  df.train$.pred_ordered = sort(df.train$.pred) # order predictions
  ggplot(df.train, aes(.pred_ordered, logconc_ordered)) + geom_point(aes(color='train')) +
    geom_abline(slope=1, intercept=0, color='red') +
    geom_point(data=df.test, aes(logconc_ordered, .pred_ordered, color='test')) +
    scale_color_manual(name='',values=c('train'='red','test'='blue'), 
                       breaks=c('train','test')) + 
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
    ylab('Modeled values (log ug/L), ordered') + xlab('Observed values (log ug/L), ordered') + 
    ggtitle(paste0(metal))
}

plot_trans_cdf <- function(df.test.adj, df.train) { # changed df.train to not include adjusted
  # cdf of test pred/obs/adjusted values (view how much adjustment was done)
  ggplot(df.test.adj, aes(logconc)) + stat_ecdf(geom='line', aes(color='Observed')) +
    stat_ecdf(data=df.train, aes(.pred, color='Train Predictions'), linewidth=1.1) + 
    stat_ecdf(data=df.test.adj, aes(.pred_orig, color='Predicted')) +  # this should match .pred_check
    stat_ecdf(data=df.test.adj, aes(.pred, color='Adjusted')) + 
    scale_color_manual(name='',values=c('Train Predictions'='chartreuse3','Observed'='black','Predicted'='red','Adjusted'='blue'), 
                       breaks=c('Observed','Train Predictions', 'Predicted','Adjusted')) + 
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
    ylab('CDF') + xlab('log ug/L') + 
    ggtitle(paste0(metal))
}


## calculate errors
calculate_reg_errors <- function(df) {
  df = df %>% mutate(residuals = ifelse(is.imputed == 0, (conc - antilog_pred), NA)) # residuals, calculated for uncensored data only
  df = df %>% mutate(factor_error = ifelse(is.imputed == 0, (antilog_pred / conc), NA)) # factor error, calculated for uncensored data only
  df = df %>% mutate(pct_error = ifelse(is.imputed == 0, ((antilog_pred - conc) * 100 / antilog_pred), NA)) # percent error, calculated for uncensored data only
  df = df %>% dplyr::mutate(factor_error_plot = ifelse(factor_error < 1, (-1/factor_error), factor_error))
  
  return(df)
}

calculate_class_errors <- function(df, MCL) {
  df$pred_exceeds = as.numeric(df$antilog_pred >= MCL)
  df$conc_exceeds = as.numeric(df$conc >= MCL)
  df = df %>% mutate(class_error = ifelse(pred_exceeds == conc_exceeds, 0, 1)) # indicates a misclassification 
  df = df %>% mutate(class_error_type = ifelse((pred_exceeds==1 & conc_exceeds==0), 'FP',
                                               ifelse((pred_exceeds==0 & conc_exceeds==1), 'FN',
                                               ifelse((pred_exceeds==1 & conc_exceeds==1), 'TP', 'TN')))) # indicates misclassification type
  return(df)
}



#### Step 1. Make model predictions -----------------------------------------------------------------------------------------------
# test data (full dataset for classification; uncensored data for regression)
test.predictions = predict(final.model, df.test) %>% bind_cols(df.test)
test.predictions$logconc = log10(test.predictions$censored.conc)
test.predictions$antilog_pred = 10**test.predictions$.pred

test.predictions.uncens = subset(test.predictions, is.imputed==0)

# train data 
train.predictions = predict(final.model, df.train) %>% bind_cols(df.train)
train.predictions$logconc = log10(train.predictions$censored.conc)
train.predictions$antilog_pred = 10**train.predictions$.pred

# uncensored train data subset (easy access for EDM adjustment)
train.predictions.uncens = subset(train.predictions, is.imputed==0)


#### Step 2. Evaluate original and adjusted model results -----------------------------------------------------------------------------

# EDM transformation

# version using full train predictions
test.predictions.adj = EDM_transform(test.predictions, train.predictions)

# version using uncensored train predictions only 
test.predictions.adj.uncens = EDM_transform(test.predictions, train.predictions.uncens) 

test.predictions.adj.plot = EDM_transform(subset(test.predictions, is.imputed==0), train.predictions) # showing uncensored only
plot_trans_ordered(test.predictions.adj.plot, train.predictions.uncens) 
plot_trans_cdf(test.predictions.adj.plot, train.predictions.uncens)

# A. ORIGINAL test data
test_reg = calc_reg_metrics(subset(test.predictions, is.imputed==0)) %>% mutate(group='test') # uncensored data only
test_class = calc_class_metrics(test.predictions, MCL) %>% mutate(group='test') # whole dataset

# B. ADJUSTED test data 
test_regAdj = calc_reg_metrics(subset(test.predictions.adj.uncens, is.imputed==0)) %>% mutate(group='test adj') # uncensored data only # some of the predictions are NA 
test_classAdj = calc_class_metrics(test.predictions.adj, MCL) %>% mutate(group='test adj') # whole dataset

test_regAdj = calc_reg_metrics(subset(test.predictions.adj.uncens, (conc = 10 | conc > 10))) %>% mutate(group='test adj') 

# C. TRAIN data
train_reg = calc_reg_metrics(train.predictions.uncens) %>% mutate(group='train')
train_class = calc_class_metrics(train.predictions, MCL) %>% mutate(group='train') # how well does the model classify censored data < MCL?

# D. Compile all results into a table
df_metrics = rbind(test_reg, test_class, test_regAdj, test_classAdj, train_reg, train_class) %>% 
  spread(group, estimate) %>% 
  arrange(factor(metric, levels=c('sens','spec','accuracy','rsq','rmse','mae'))) %>%
  dplyr::mutate_if(is.numeric, round,3) %>% 
  dplyr::mutate(metric = recode(metric,'sens'='sensitivity','spec'='specificity','rsq'='R2','rmse'='RMSE','mae'='MAE'))

formattable(df_metrics,align = c("l", rep("r", ncol(df_metrics) - 1)))


#### Step 3. Calculate residuals/errors and write out results ------------------------------------------------------------------------------------------
# original results 
test.predictions = calculate_reg_errors(test.predictions)
test.predictions = calculate_class_errors(test.predictions, MCL)

# adjusted results
test.predictions.adj = calculate_reg_errors(test.predictions.adj.uncens) # added uncens Jan 2024
test.predictions.adj = calculate_class_errors(test.predictions.adj, MCL)

# train data
train.predictions = calculate_reg_errors(train.predictions.uncens) # added uncens Jan 2024
train.predictions = calculate_class_errors(train.predictions, MCL)

# write out data 
writename = paste0('R_Output/',modelfolder.name,metal.code,"_Predicts.RData")
save(test.predictions, test.predictions.adj, train.predictions, file=writename)

