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
