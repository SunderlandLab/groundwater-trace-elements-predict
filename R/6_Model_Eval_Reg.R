####################################################################################################
#File name: Model_Eval_Reg.R
#Author: Jennifer Sun, Cindy Hu
#Date: March 2025
#Purpose: Model evaluation, including EDM correction and both regression and classification metrics
####################################################################################################

#### Attach Libraries and Set Working Directory ####

source(here::here('R/0_helper_fct.R'))
plan(multisession, workers = future::availableCores())
print(paste("The number of workers are", nbrOfWorkers()))
options(future.rng.onMisuse = "ignore")

setwd(here::here("data"))

#### Define Model Version ####
metal.codes <- c("As", "Cd", "Li", "Mn", "Sr")
metals <- c("Arsenic", "Cadmium", "Lithium", "Manganese", "Strontium")
names(metals) <- metal.codes
MCLs <- c(10, 5, 60, 50, 4000)
names(MCLs) <- metal.codes

#### Helper functions  -------------------------------------------------------------------------------------------------------

## print result metrics
calc_reg_metrics <- function(list_of_dfs, group_label = "test") {
  list_of_dfs %>%
    purrr::map(~ {
      yardstick::metrics(.x, logconc, .pred) %>%
        as.data.frame() %>%
        dplyr::select(-.estimator) %>%
        dplyr::rename(metric = .metric, estimate = .estimate)
    }) %>%
    bind_rows(.id = "imputation") %>%
    dplyr::group_by(metric) %>%
    dplyr::summarise(pooled_estimate = mean(estimate), .groups = "drop") %>%
    dplyr::mutate(group = group_label)
}

calc_class_metrics <- function(list_df_results, threshold,  group_label = "test") {
  multi_met <- metric_set(accuracy, spec, sens)
  
  # function for a single df
  calc_metrics <- function(df_results) {
    df_results$pred_exceeds <- as.factor(df_results$antilog_pred >= threshold)
    df_results$conc_exceeds <- as.factor(df_results$conc >= threshold)
    
    metrics_df <- df_results %>%
      multi_met(
        truth = conc_exceeds,
        estimate = pred_exceeds,
        event_level = "second"
      ) %>%
      as.data.frame() %>%
      dplyr::select(-c(.estimator)) %>%
      dplyr::rename(metric = .metric, estimate = .estimate)
    
    return(metrics_df)
  }
  
  # Apply to each imputed dataset
  all_metrics <- lapply(list_df_results, calc_metrics)
  combined <- bind_rows(all_metrics, .id = "imputation_id")
  
  # Pool by averaging across imputations
  pooled <- combined %>%
    group_by(metric) %>%
    summarise(pooled_estimate = mean(estimate), .groups = "drop")%>%
    dplyr::mutate(group = group_label)
  
  return(pooled)
}


## EDM transformation
EDM_transform <- function(test.pred.adj, train.pred.adj) {
  # order train values
  tbl_ordered <- train.pred.adj %>%
    transmute(logconc_ordered = sort(logconc),
              .pred_ordered = sort(.pred))
  
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
    ggtitle(paste0(metal.code))
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
    ggtitle(paste0(metal.code))
}


## calculate errors
calculate_reg_errors <- function(df) {
  df <- df %>%
    mutate(
      residuals = if_else(is.imputed == 0, conc - antilog_pred, NA_real_),
      factor_error = if_else(is.imputed == 0, antilog_pred / conc, NA_real_),
      pct_error = if_else(is.imputed == 0, (antilog_pred - conc) * 100 / antilog_pred, NA_real_),
      factor_error_plot = if_else(factor_error < 1, -1 / factor_error, factor_error)
    )
  return(df)
}

calculate_class_errors <- function(df, MCL) {
  df <- df %>%
    mutate(
      pred_exceeds = as.numeric(antilog_pred >= MCL),
      conc_exceeds = as.numeric(conc >= MCL),
      class_error = if_else(pred_exceeds == conc_exceeds, 0, 1),
      class_error_type = case_when(
        pred_exceeds == 1 & conc_exceeds == 0 ~ "FP",
        pred_exceeds == 0 & conc_exceeds == 1 ~ "FN",
        pred_exceeds == 1 & conc_exceeds == 1 ~ "TP",
        pred_exceeds == 0 & conc_exceeds == 0 ~ "TN"
      )
    )
  return(df)
}

evaluate_and_predict <- function(metal.code){
  print(paste("Begin evaluation and prediction for", metal.code))
  MCL <- MCLs[metal.code]
  filename <- paste0('R_Output/', metal.code,"_ModelPackage.RData")
  load(filename)
  
  #### Step 1. Make model predictions -----------------------------------------------------------------------------------------------
  # Create predictions for test data
  test_predictions <- pmap(
    list(final_model, df_test),
    ~ predict(..1, ..2) %>%
      bind_cols(..2) %>%
      mutate(
        logconc = log10(conc),
        antilog_pred = 10^.pred
      )
  )
  
  # Subset uncensored test data
  test_predictions_uncens <- map(test_predictions, ~ filter(.x, is.imputed == 0))
  
  #  Create predictions for train data
  train_predictions <- pmap(
    list(final_model, df_train),
    ~ predict(..1, ..2) %>%
      bind_cols(..2) %>%
      mutate(
        logconc = log10(conc),
        antilog_pred = 10^.pred
      )
  )
  
  # Subset uncensored train data
  train_predictions_uncens <- map(train_predictions, ~ filter(.x, is.imputed == 0))
  
  #### Step 2. Evaluate original and adjusted model results -----------------------------------------------------------------------------
  
  # EDM transformation
  test_predictions_adj<- map2(test_predictions, train_predictions, EDM_transform)
  # Uncensored version
  test_predictions_adj_uncens <- map2(test_predictions, train_predictions_uncens, EDM_transform)
  test_predictions_adj_plot <- map2(test_predictions_uncens, train_predictions, EDM_transform)
  # diagnostic plots
  # map2(test_predictions_adj_plot, train_predictions_uncens, plot_trans_ordered) 
  map2(test_predictions_adj_plot, train_predictions_uncens, plot_trans_cdf) %>%
    imap( ~ {
      # Get original title if it exists
      original_title <- .x$labels$title %||% ""  # use `%||%` to avoid NULL
      
      # Create new title by appending
      new_title <- paste(original_title, "-", paste("Imputation", .y))
      
      # Add updated title
      .x + ggtitle(new_title)
    }) %>%
    patchwork::wrap_plots(plot_list_labeled, ncol = 2)
  # save ggplot object
  ggsave(paste0('R_Output/',metal.code,'_EDM_transformation_plot.png'), width = 8, height = 6)
  # A. ORIGINAL test data
  # pool across 5 imputations
  test_reg <- calc_reg_metrics(test_predictions_uncens, group_label = "test")
  test_class <- calc_class_metrics(test_predictions, MCL, group_label = "test")
  # B. ADJUSTED test data 
  test_regAdj  <- calc_reg_metrics(test_predictions_adj_uncens, group_label = "test adj")
  test_classAdj <- calc_class_metrics(test_predictions_adj, MCL, group_label = "test adj")
  # C. TRAIN data
  train_reg <- calc_reg_metrics(train_predictions_uncens, group_label = "train")
  train_class  <- calc_class_metrics(train_predictions, MCL, group_label ='train')
  # D. Compile all results into a table
  df_metrics <- rbind(test_reg, test_class, test_regAdj, test_classAdj, train_reg, train_class) %>% 
    spread(group, pooled_estimate) %>% 
    arrange(factor(metric, levels=c('sens','spec','accuracy','rsq','rmse','mae'))) %>%
    dplyr::mutate_if(is.numeric, round,3) %>% 
    dplyr::mutate(metric = recode(metric,'sens'='sensitivity','spec'='specificity','rsq'='R2','rmse'='RMSE','mae'='MAE'))
  # write out results 
  write_csv(df_metrics, paste0('R_Output/',metal.code,'_Model_Eval_Metrics.csv'))
  
  #### Step 3. Calculate residuals/errors and write out results ------------------------------------------------------------------------------------------
  # original results 
  test_predictions <- map(test_predictions, calculate_reg_errors)
  test_predictions <- map(test_predictions, ~calculate_class_errors(.x, MCL))
  
  # adjusted results
  test_predictions_adj <- map(test_predictions_adj, calculate_reg_errors) 
  test_predictions_adj <- map(test_predictions_adj, ~calculate_class_errors(.x, MCL))
  
  # train data
  train_predictions <- map(train_predictions, calculate_reg_errors) 
  train_predictions <- map(train_predictions, ~calculate_class_errors(.x, MCL))
  
  # write out data 
  writename <- paste0('R_Output/', metal.code,"_Predicts.RData")
  save(test_predictions, test_predictions_adj, train_predictions, file=writename)
  print(paste("evaluation and prediction for", metal.code, "is finished."))
}

# iterate over five trace elements and save model packages
purrr::map(metal.codes, evaluate_and_predict)
