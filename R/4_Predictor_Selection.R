############################################################################################
#File name: Predictor_Selection.R
#Author: Jennifer Sun, Cindy Hu
#Date: Jan 2025
#Purpose: Final dataframe cleaning steps, including manual and automatic predictor selection
############################################################################################

source(here::here('R/0_helper_fct.R'))

setwd(here::here("data"))

#### Define Model Version ####
metal.codes <- c("As", "Cd", "Li", "Mn", "Sr")
metals <- c("Arsenic", "Cadmium", "Lithium", "Manganese", "Strontium")
names(metals) <- metal.codes

select_predictors <- function(metal.code) {
  print(paste0('Begin variable selection for ', metal.code))
  # Read the refined data set
  df <- readRDS(paste0("R_Output/", metal.code, "_imputed_data.rds"))
  percent_missing <- data.frame(column = names(df),
                                percent_missing = sapply(df, function(x)
                                  mean(is.na(x)) * 100)) %>%
    arrange(desc(percent_missing))
  
  ## 1. Remove unwanted predictor variables -----------------------------------------------------------------------------------------------------------------
  variables_missing_gt10pct = percent_missing %>%
    filter(percent_missing > 10) %>%
    pull(column) %>%
    setdiff("well.depth") # keep well depth
  
  df <- df %>%
    dplyr::select(-all_of(variables_missing_gt10pct)) # variables with > 10% missingness
  
  ## 2. Remove variables that are replicated within raw datasets -------------------------------------------------------------------------------------
  # duplicates from statsgo soil properties; use calculated mean instead
  prefixes <- str_replace(colnames(df), "(L|H|mean)$", "") %>%
    # find the ones that appear three times
    tibble(value = .) %>%
    count(value, name = 'count') %>%
    filter(count == 3) %>%
    pull(value)
  df <- df %>%
    dplyr::select(-c(
      paste0(prefixes, 'L'),
      paste0(prefixes, 'H'),
      'PFLATLOW',
      'PFLATUP'
    ))
  
  # df = subset(df, select=-c(cec_025, cec_05, cec_050, clay_025, clay_05, clay_2550, clay_3060,
  #                           ec_025, ec_05, ksat_05, ph_025, ph_05, ph_2550, ph_3060, sand_025, sand_05,
  #                           sand_2550, sand_3060, silt_025, silt_05, silt_2550, silt_3060, paws_025, paws_050,
  #                           min_ksat, max_ksat)) # duplicates from #16. soil property rasters
  
  
  ## 3. Remove highly correlated variables -----------------------------------------------------------------------------------------------------------------------
  # remove irrelevant variables before running correlation analysis
  df_corData <- df %>%
    dplyr::select(where(is.numeric)) %>%
    dplyr::select(-any_of(c(
      'lon', 'lat', 'conc', 'well.depth', 'DL.missing'
    ))) %>% # keep key vars
    dplyr::select(-any_of(c(
      'TRI.total.impact', 'TRI.water.impact', 'SEMS'
    ))) # anthropogenic input variables that may be correlated but we will keep
  
  # identify highly correlated variables
  df_corMat = cor(df_corData, method = 'pearson') # manually inspect
  # automatic removal of correlated variables - did not use (instead, manually select which of the correlated variables to remove)
  corVars <- caret::findCorrelation(df_corMat, cutoff = 0.9)
  corVars <- colnames(df_corData)[corVars]
  
  # manually selected parameters to remove based on correlation coefficient (view correlation matrix)
  if (metal.code == 'Sr') {
    df <- df %>%
      dplyr::select(-any_of(
        c(
          'C_Ni',
          'C_Tot_Flds',
          'AVG_POR',
          'Hydrate_Y',
          'AVG_NO10',
          'AVG_SAND',
          'AVG_CLAY',
          'AVG_SILT',
          'sand',
          'silt',
          'no3_pub',
          'AVG_KV'
        )
      ))
  } else if (metal.code == 'Li') {
    df <- df %>%
      dplyr::select(-any_of(
        c(
          'C_Ni',
          'C_Tot_Flds',
          'AVG_POR',
          'Hydrate_Y',
          'AVG_NO10',
          'AVG_SAND',
          'AVG_CLAY',
          'AVG_SILT',
          'sand',
          'silt',
          'AVG_NO200',
          'no3_pub'
        )
      ))
  } else if (metal.code == 'As') {
    df <- df %>%
      dplyr::select(-any_of(
        c(
          'C_Ni',
          'C_Tot_Flds',
          'AVG_POR',
          'Hydrate_Y',
          'AVG_NO10',
          'AVG_SAND',
          'AVG_CLAY',
          'AVG_SILT',
          'sand',
          'silt',
          'aq_rocktype'
        )
      ))
  } else if (metal.code %in% c('Cd', 'Mn')) {
    df <- df %>%
      dplyr::select(-any_of(
        c(
          'C_Ni',
          'C_Tot_Flds',
          'AVG_POR',
          'Hydrate_Y',
          'AVG_NO10',
          'AVG_SAND',
          'AVG_CLAY',
          'AVG_SILT',
          'sand',
          'silt'
        )
      ))
  }
  
  # Save final dataframe
  saveRDS(df,
          paste0("R_Output/", metal.code, "_df_PredictorsSelected.rds"))
  print(paste0('Variable selection for ', metal.code, ' complete'))
}


# vectorize across all five metals
purrr::map(metal.codes, select_predictors)
