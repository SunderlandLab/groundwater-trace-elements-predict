###########################################################################
#File name: Shapley_Analysis.R
#Author: Jennifer Sun, Cindy Hu
#Date: April 2023, March 2025
#Purpose: Analyze shapley values
###########################################################################

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

#### read in auxilary geospatial data for regional shapley value 
mapUSm <- read_sf(dsn="CoVar/US48", layer="US_48states") # load projected map of the US (in m)
US.df <- fortify(mapUSm)
aquifershp <- read_sf(dsn="CoVar/aquifrp025_nt00003", layer="aquifrp025") %>% st_set_crs(4269)  
# CRS info comes from meta data https://catalog.data.gov/dataset/aquifers1
# https://catalog.data.gov/harvest/object/a92b0a2a-2a5e-4709-879d-e78666d4c471
pennsylvanian_aquifers <- aquifershp %>% filter(AQ_NAME == "Pennsylvanian aquifers") %>%
  dplyr::select(geometry) 
mississippi_river_valley_aquifers <- aquifershp %>% filter(str_detect(AQ_NAME, "Mississippi River Valley"))%>%
  dplyr::select(geometry) 


#### rename features 
relabels <- c('ph_2550' = 'pH',
             'SLOPE'='slope',
             'KB'='bedrock',
             'rech'='recharge',
             'caco3_kg_sq_m'='Calcite',
             'C_Hornbl' = 'Hornblende, C horizon',
             'no3_pub' = 'public well NO3',
             'no3_dom' = 'domestic well NO3',
             'C_Sr' = 'Strontium, C horizon',
             'C_Hornbl' = 'Hornblende, C horizon',
             'Hydrate_N' = '% non-hydric soils',
             'WTDEPmean' = 'depth to water table',
             'VALUE_52_pct' = '% shrubland',
             'ppt' = 'precipitation',
             'C_Tot_Plag' = 'total plaggen',
             'C_Sb' = 'Antimony, C horizon',
             'AQPERMNEW'='aquifer permeability',
             'C_Pb' = 'Lead, C horizon',
             'stm_den'='stream density',
             'quickflow'='quick overland flow',
             'C_Calcite' = 'Calcite, C horizon',
             'VALUE_21_pct' = '% Developed, Open Space',
             'sar' = 'sodium adsorption ratio',
             'VALUE_31_pct' = '% Barren Land',
             'SEMS' = 'distance to Superfund site',
             'C_Cd' = 'Cadmium, C horizon',
             'C_Ti' = 'Titanium, C horizon',
             'HLR'= 'Hydrologic Landscape Region',
             'mean_ksat' = 'hydraulic conductivity',
             'OMmean'='Organic matter', 
             'tmean' = 'temperature',
             'bfi' = 'base flow index',
             'soiloc' = 'soil organic carbon',
             'C_As' = 'Arsenic, C horizon',
             'unsatWC' = 'unsat zone water content',
             'C_Li' = 'Lithium, C horizon',
             'surfgeo' = 'Surficial geology',
             'C_Na' = 'sodium, C horizon',
             'VALUE_11_pct' = '% open water',
             'well.depth' = 'well depth',
             'trans' = 'transmissivity',
             'SLOPEmean' = 'slope',
             'C_K' = 'Potassium, C horizon',
             'VALUE_42_pct' = '% Evergreen Forest')

run_shapley_analysis <- function(metal.code, region_shapefile = NULL){
  print(paste("Begin shapley analysis for", metal.code))
  filename <- paste0('R_Output/', metal.code,"_ModelPackage_2step.RData")
  load(filename)
  
  #### 1. Prepare data --------------------------------------------------------------------------------------------------------
  ## extract xgboost object from model trained on training data, for diagnostics
  final_xg_list <- map(final_model, extract_fit_engine) 
  # Initialize output lists
  df_pred_list <- list()
  factorcols_list <- list()
  
  # Loop over imputed datasets
  for (i in seq_along(df_train)) {
    df <- rbind(df_train[[i]], df_test[[i]])  
    
    # Process predictors using recipe
    df_pred <- bake(
      prep(train_recipes[[i]]), # bake with the reference recipe to ensure consistency
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
    remove_features <- c(setdiff(colnames(df_pred), final_xg_list[[i]]$feature_names))
    print(paste0("features to remove from df_pred are ", remove_features))
    remove_index <- which(colnames(df_pred) %in% remove_features)
    if (length(remove_index) > 0) {
      df_pred <- df_pred[, -remove_index]
    }
    
    # Save results
    df_pred_list[[i]] <- df_pred
    factorcols_list[[i]] <- factorcols
  }
  
  #### 2a. Calculate xgboost shapley values --------------------------------------------------------------------------------------------------------
  if(is.null(region_shapefile)){
    # national analysis
    shap_list <- map(seq_along(df_train), function(i) {
      print(i)
      df <- rbind(df_train[[i]], df_test[[i]]) 
      df_pred <- df_pred_list[[i]] 
      factorcols <- factorcols_list[[i]]

      shapviz(
        object = final_xg_list[[i]],
        X_pred = df_pred,
        X = df,
        collapse = factorcols
      )
    })

  }else{
    # regional analysis, filter by aquifer shapefile
      shap_list <- map(seq_along(df_train), function(i) {
        df_pred <- df_pred_list[[i]]
        df_region <- rbind(df_train[[i]], df_test[[i]]) %>%
          mutate(row_id = row_number()) %>%
          st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
          st_transform(st_crs(region_shapefile)) %>%
          st_join(region_shapefile, join = st_within, left = FALSE)  %>%# inner join
          st_drop_geometry() 
        
        matched_idx <- df_region$row_id #keep the index to filter the df_pred_list
        df_region <- df_region %>% dplyr::select(-c('row_id')) #df_region must have the same number of columns as df
        factorcols <- factorcols_list[[i]]
        
        shapviz(
          object = final_xg_list[[i]],
          X_pred = df_pred[matched_idx, ],
          X = df_region,
          collapse = factorcols
        )
    })
  }
  saveRDS(shap_list, paste0('R_Output/', metal.code, "_", deparse(substitute(region_shapefile)), '_shap_list_2step.rds'))
}

# 5 national shapley value analyses
purrr::map(metal.codes[2], ~run_shapley_analysis(.x, region_shapefile = NULL))
# 2 regional shapley value analyses
# run_shapley_analysis(metal.code = "Mn", region_shapefile = pennsylvanian_aquifers)
# run_shapley_analysis(metal.code = "Mn", region_shapefile = mississippi_river_valley_aquifers)

visualize_shapley_analysis <- function(shap_list_name){
  shap_list <- readRDS(paste0('R_Output/', shap_list_name, '_shap_list_2step.rds'))
  # 1. Extract SHAP matrices
  shap_matrices <- lapply(shap_list, function(sv) sv$S)
  # 2. Pool SHAP values by averaging across imputations
  # Find common column names
  common_cols <- Reduce(intersect, lapply(shap_matrices, colnames))
  
  # Subset each matrix to only those common columns
  shap_matrices_common <- lapply(shap_matrices, function(mat) mat[, common_cols, drop = FALSE])
  
  # Average them
  pooled_S <- Reduce(`+`, shap_matrices_common) / length(shap_matrices_common)
  # 3. Pool baseline values
  baseline_values <- sapply(shap_list, function(sv) sv$baseline)
  pooled_baseline <- mean(baseline_values)
  # 4. Use the feature matrix from the first object
  X_ref <- shap_list[[1]]$X
  # 5. Reconstruct a new shapviz object
  shap_pooled <- shapviz(
    pooled_S,
    X = X_ref,
    baseline_value = pooled_baseline,
    model_class = "xgboost"
  )
  varImp.plot <- sv_importance(shap_pooled, kind = "both", show_numbers = TRUE) + 
    scale_y_discrete(limits=rev(calc_meanabs_shap(shap_pooled)), labels=relabels) +  
    theme_classic() + 
    theme(axis.text.x = element_text(size=12), axis.text.y=element_text(size=12))  
  ggsave(varImp.plot, filename = paste0('R_Output/', shap_list_name, '_SHAP_varImp_plot_2step.png'), width = 6, height = 5)
}

# visualize 5 national shapley values
purrr::map(paste0(metal.codes[2], "_NULL"), visualize_shapley_analysis)
# visualize 2 regional shapley values
# visualize_shapley_analysis("Mn_pennsylvanian_aquifers")
# visualize_shapley_analysis("Mn_mississippi_river_valley_aquifers")

