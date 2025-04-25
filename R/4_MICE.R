##################################################################################
#File name: MICE.R
#Author: Cindy Hu
#Date: Jan 2025
#Purpose: MICE imputation for censored concentrations
# R v4.4.2
##################################################################################

source(here::here('R/0_helper_fct.R'))

setwd(here::here("data"))

#### Define Model Version ####
metal.codes <- c("As", "Cd", "Li", "Mn", "Sr")
metals <- c("Arsenic", "Cadmium", "Lithium", "Manganese", "Strontium")
names(metals) <- metal.codes
MCLs <- c(10, 5, 60, 300, 4000)
names(MCLs) <- metal.codes

impute_missing_values <- function(metal.code) {
  # Load Data
  master <- readRDS(paste0("Data_Files/", metal.code, "_df_PredictorsSelected.rds")) %>%
    select(-c(censored.conc, ros.conc)) %>%
    # make aquifer a factor variable
    mutate(aquifer = as.factor(aquifer))
  if(!"detect.limit"%in%colnames(master)){
    # bring back limit of detection
    sup_master <- readRDS(paste0("R_Output/", metal.code, "_RawPredictors.rds")) %>%
      st_drop_geometry() %>%
      dplyr::select(location.id, detect.limit)
    master <- master %>%
      left_join(sup_master, by = "location.id")
  }
  
  percent_missing <- data.frame(
    column = names(master),
    percent_missing = sapply(master, function(x)
      mean(is.na(x)) * 100)
  ) %>%
    arrange(desc(percent_missing))
  
  ### Step 1, missingness has meaning, set missing drainage to zero
  master$drainage[is.na(master$drainage)] = 0
  
  ### Step 2, factor variables, create a missing category
  factor_column_names <- names(master)[sapply(master, is.factor)]
  character_column_names <- names(master)[sapply(master, is.character)]
  factor_column_to_fill <- percent_missing %>%
    filter(column %in% factor_column_names |
             column %in% character_column_names) %>%
    filter(percent_missing > 0) %>%
    pull(column)
  
  master <- master %>%
    mutate(across(
      factor_column_to_fill,
      ~ forcats::fct_na_value_to_level(.x, "missing")
    ))
  
  ### Step 3, numeric variables, impute missing values with KNN
  numeric_column_to_fill <- percent_missing %>%
    filter(!(column %in% factor_column_names) &
             !(column %in% character_column_names)) %>%
    filter(percent_missing > 0) %>%
    pull(column) %>%
    # do not impute conc, well.depth
    setdiff(c('conc', 'well.depth'))
  
  if(length(numeric_column_to_fill)>0){
    master_imp <- master %>%
    rename(lon = long) %>%
    VIM::kNN(
      variable = numeric_column_to_fill,
      k = 5,
      dist_var = c("lon", "lat"),
      impNA = FALSE
    ) %>%
    dplyr::select(-ends_with("_imp"))
  
  percent_missing_imp <- data.frame(
    column = names(master_imp),
    percent_missing = sapply(master_imp, function(x)
      mean(is.na(x)) * 100)
  ) %>%
    arrange(desc(percent_missing))
  } else{
    master_imp <- master
  }
  
  ### Step 4, impute censored conc with MICE
  
  # if censored samples don't have detect.limit, drop the observation
  master_imp <- master_imp %>%
    filter(!(censored & is.na(detect.limit)))
  
  if(nrow(master_imp)<nrow(master)){
    print("pause to check sample size in imputation input.")
  }
  
  print(paste("percent censored is ",
              round(sum(master_imp$censored) / nrow(master_imp) * 100, 2)))
  # We run the mice code with 0 iterations
  imp <- mice(master_imp, maxit = 0)
  # extractr predictorMatrix and methods of imputation
  predM <- as_tibble(imp$predictorMatrix)
  # specify variables not to impute
  cols_to_zero <- percent_missing_imp %>%
    # columns with too much missingness
    filter(percent_missing > 10) %>%
    pull(column) %>%
    # these need to be imputed, so exclude them from cols_to_zero
    setdiff(c("conc"))
  # Use mutate(across()) to set them to 0
  predM <- predM %>%
    mutate(across(all_of(cols_to_zero), ~ 0))
  # If you need to convert back to a matrix
  predM <- as.matrix(predM)
  rownames(predM) <- colnames(master_imp)
  colnames(predM) <- colnames(master_imp)
  meth <- imp$method
  meth[cols_to_zero] <- ""
  meth["conc"] <- "conc_below_limit"  # Assign custom method
  
  mice.impute.conc_below_limit <- function(y, ry, x, ...) {
    # Ensure "detect.limit" exists in x
    if (!"detect.limit" %in% colnames(x)) {
      stop("Error: 'detect.limit' column is missing from predictor matrix.")
    }
    detect_limit <- x[!ry, "detect.limit"]  # Extract detection limits for missing rows
    
    # Default PMM imputation
    imputed_values <- mice.impute.pmm(y, ry, x, ...)
    
    # Generate random values below the detection limit
    imputed_values <- EnvStats::rlnormTrunc(
      length(imputed_values),
      min = 0.001,
      # Lower bound (assumes non-negative conc)
      max = detect_limit,
      # Upper bound (detection limit)
      meanlog = log(imputed_values),
      sd = sd(log(y), na.rm = TRUE)
    )  # Use observed SD
    
    return(imputed_values)
  }
  
  # MICE imputation step
  imp2 <- tryCatch({
    message("Attempting MICE imputation (first try)...")
    
    # First attempt
    mice(
      master_imp,
      maxit = 5,
      predictorMatrix = predM,
      method = meth,
      print = TRUE,
      seed = 123
    )
  }, error = function(e) {
    message("Error encountered: ", e$message)
    message("Handling singularity: Removing near-zero variance variables and retrying MICE...")
    
    
    # if failed due to singularity, try again after removing the vars with near zero variance
    master_imp <- master_imp %>%
      dplyr::select(-(caret::nearZeroVar(master_imp, freqCut = 999 / 1)))
    # We run the mice code with 0 iterations
    imp <- mice(master_imp, maxit = 0)
    # extractr predictorMatrix and methods of imputation
    predM <- as_tibble(imp$predictorMatrix)
    # specify variables not to impute
    cols_to_zero <- percent_missing_imp %>%
      # columns with too much missingness
      filter(percent_missing > 10) %>%
      pull(column) %>%
      # these need to be imputed, so exclude them from cols_to_zero
      setdiff(c("conc", "detect.limit"))
    # Use mutate(across()) to set them to 0
    predM <- predM %>%
      mutate(across(all_of(cols_to_zero), ~ 0))
    # If you need to convert back to a matrix
    predM <- as.matrix(predM)
    rownames(predM) <- colnames(master_imp)
    colnames(predM) <- colnames(master_imp)
    meth <- imp$method
    meth[cols_to_zero] <- ""
    meth["conc"] <- "conc_below_limit"  # Assign custom method
    
    message("Retrying MICE imputation (after adjustments)...")
    
    # second attempt
    mice(
      master_imp,
      maxit = 5,
      predictorMatrix = predM,
      method = meth,
      print =  TRUE,
      seed = 123
    )
  })
  
  percent_missing_imp2 <- complete(imp2, "long") %>%
    group_by(.imp) %>%
    summarise(across(everything(), ~ mean(is.na(.x)) * 100)) %>%
    arrange(desc(conc))
  
  # inspect quality of imputations
  bind_rows(master_imp, complete(imp2, "long")) %>%
    mutate(.imp = ifelse(is.na(.imp), "original", .imp)) %>%
    #visualize density plot, by .imp and censored
    ggplot(aes(x = .imp, y = conc)) +
    geom_jitter(aes(color = censored), width = 0.25, alpha = 0.5) +
    # log transform y axis
    scale_y_sqrt() +
    # rename x-axis to imputation, rename y-axis to concentration
    labs(x = "Imputation", y = "Concentration") +
    # add the title metal.code +
    ggtitle(paste0("Imputation of ", metals[metal.code])) +
    # add a horizontal line at the detection limit 5
    geom_hline(
      yintercept = quantile(master_imp$detect.limit, 0.95, na.rm = TRUE),
      linetype = "dashed"
    ) +
    theme_minimal(base_size = 9)
  
  ggsave(
    paste0(
      "R_Output/",
      metal.code,
      "_imputation_qualitycheck_plot.png"
    ),
    width = 6,
    height = 4,
    dpi = 300
  )
  # save imputed data
  saveRDS(complete(imp2, "long"),
          paste0("R_Output/", metal.code, "_imputed_data.rds"))
}

# vectorize across all five metals
purrr::map(metal.codes, impute_missing_values)

# Create Table 1
# read in imputated data, compute summary statistcs, such as number of wells, percent not censorsored, concentration 
# min 25th percentile, median, 75th percentile, and max
# % above MCL
imputed_data <- purrr::map_dfr(metal.codes, ~ readRDS(paste0("R_Output/", .x, "_imputed_data.rds")), .id = "metal") %>%
  mutate(metal = metal.codes[as.integer(metal)]) %>%
  group_by(metal) %>%
  summarise(
    n_wells = n_distinct(location.id),
    n_samples = n(),
    percent_not_censored = sum(!censored) / n() * 100,
    conc_min = min(conc, na.rm = TRUE),
    conc_25th = quantile(conc, 0.25, na.rm = TRUE),
    conc_median = median(conc, na.rm = TRUE),
    conc_75th = quantile(conc, 0.75, na.rm = TRUE),
    conc_max = max(conc, na.rm = TRUE),
    percent_above_MCL = sum(conc > MCLs[metal]) / n() * 100
  )
# save table 1
write.csv(imputed_data, "R_Output/Table1_imputed_data_summary.csv", row.names = FALSE)
