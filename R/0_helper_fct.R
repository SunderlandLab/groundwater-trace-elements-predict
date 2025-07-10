packages <- c(
  "akima",
  "beepr",
  "chron",
  "caret",
  "EnvStats",
  "foreign",
  "furrr",
  "future",
  "geosphere",
  "ggplot2",
  "ggpmisc",
  "gstat",
  "lattice",
  "NADA",
  "plyr",
  "parallel",
  "raster",
  "readxl",
  "recipes",
  "randomForest",
  "sf",
  "shapviz",
  "SHAPforxgboost",
  "stringr",
  "terra",
  "tiff",
  "VIM",
  "dplyr",
  "assertthat",
  "data.table",
  "tidyr",
  "tidyverse",
  "tidymodels",
  "mice",
  "purrr",
  "sp",
  "xgboost",
  "yardstick"
)

lapply(packages, library, character.only = TRUE)
#### Helper functions

#' Extract Raster Values at Given Point Locations
#'
#' This function extracts raster values at specified point locations and assigns them to a new column
#' in a spatial points dataset. It can perform either direct extraction or buffered extraction using
#' a specified interpolation method.
#'
#' @param raster_path Character. The file path to the raster dataset.
#' @param points_sf sf object. A spatial points dataset (simple feature) where raster values will be extracted.
#' @param column_name Character. The name of the new column in `points_sf` to store the extracted raster values.
#' @param method Character. Interpolation method for extracting raster values. Default is `'bilinear'`.
#'   Other options include `'simple'` (nearest neighbor).
#' @param buffer_on Logical. If `TRUE`, extracts values within a buffer (500m radius) using the mean function.
#'   Default is `FALSE`.
#'
#' @return An `sf` object with the extracted raster values added as a new column.
#'
#' @details
#' - The function ensures that `points_sf` is an `sf` object.
#' - If the raster lacks a CRS, it assigns EPSG:5070 (Albers Equal Area).
#' - Points are reprojected to match the raster CRS before extraction.
#' - If `buffer_on = TRUE`, values are averaged within a 500m buffer around each point.
#' - Calls `gc()` for garbage collection to free memory after raster extraction.
#'
#' @examples
#' \dontrun{
#' # Load required libraries
#' library(sf)
#' library(raster)
#'
#' # Define file path
#' raster_path <- "path/to/raster.tif"
#'
#' # Load sample points (ensure it's an sf object)
#' points_sf <- st_read("path/to/points.geojson")
#'
#' # Extract raster values
#' points_sf <- extract_raster_values(raster_path, points_sf, column_name = "elevation")
#'
#' # Extract raster values with buffer
#' points_sf <- extract_raster_values(raster_path, points_sf, column_name = "avg_elevation", buffer_on = TRUE)
#' }
#'
#' @export
extract_raster_values <- function(raster_path,
                                  points_sf,
                                  column_name,
                                  method = 'bilinear',
                                  buffer_on = FALSE) {
  # Load the raster
  raster_layer <- raster(raster_path)
  # Ensure points_sf is an sf object
  if (!inherits(points_sf, "sf"))
    stop("points_sf must be an sf object.")
  # if raster is missing CRS, assign the CRS to it
  if (is.na(raster::crs(raster_layer))) {
    raster::crs(raster_layer) <- CRS('+init=epsg:5070')
  }
  
  # Reproject points to match the raster CRS
  points_sf <- st_transform(points_sf, crs = st_crs(raster_layer))
  
  if (!buffer_on) {
    # Extract raster values to the points
    points_sf[[column_name]] <- raster::extract(raster_layer, points_sf, method = method)
  } else {
    # Extract raster values to the points with a buffer
    points_sf[[column_name]] <- raster::extract(raster_layer,
                                                points_sf,
                                                buffer = 500,
                                                fun = mean)
  }
  
  # Clean up
  rm(raster_layer)
  gc()  # Optional: Call garbage collection
  
  # Notify completion
  message(paste(column_name, "completed"))
  
  return(points_sf)
}

#' Extract Attribute Values from a Vector Layer to Points
#'
#' This function performs a spatial join to extract attribute values from a specified vector dataset
#' (e.g., shapefile, GeoPackage) and assigns them to a column in a spatial points dataset (`sf` object).
#'
#' @param points_sf sf object. A spatial points dataset to which the attribute values will be added.
#' @param vector_dsn Character. The data source name (file path) of the vector dataset.
#' @param vector_layer Character. The specific layer name within the vector data source.
#' @param join_column Character. The column in the vector dataset to be joined with `points_sf`.
#'
#' @return An `sf` object with the extracted attribute values added as a new column.
#'
#' @details
#' - The function reads the vector data as an `sf` object.
#' - If the vector dataset lacks a CRS, it assigns the CRS of `mapUSm` (assumed to be defined elsewhere).
#' - Ensures vector geometries are valid using `st_make_valid()`.
#' - Reprojects `points_sf` to match the CRS of the vector data before performing the join.
#' - Performs a spatial join (`st_join()`) to merge the selected attribute from the vector dataset.
#' - Removes temporary variables after execution.
#'
#' @examples
#' \dontrun{
#' # Load required libraries
#' library(sf)
#' library(dplyr)
#'
#' # Define file paths and layers
#' vector_dsn <- "path/to/vector_data.gpkg"
#' vector_layer <- "land_use_layer"
#' join_column <- "land_cover_type"
#'
#' # Load sample points (ensure it's an sf object)
#' points_sf <- st_read("path/to/points.geojson")
#'
#' # Extract land cover type from vector data
#' points_sf <- extract_vector_values(points_sf, vector_dsn, vector_layer, join_column)
#' }
#'
#' @export
extract_vector_values <- function(points_sf,
                                  vector_dsn,
                                  vector_layer,
                                  join_column) {
  # Read the vector data as sf
  vector_data <- st_read(dsn = vector_dsn, layer = vector_layer)
  # if vector is missing CRS, assign the CRS of mapUSm to it
  if (is.na(st_crs(vector_data))) {
    st_crs(vector_data) <- st_crs(mapUSm)
  }
  
  vector_data <- st_make_valid(vector_data)
  
  # Transform input CRS to match the vector data CRS
  points_sf <- st_transform(points_sf, crs = st_crs(vector_data))
  
  # Perform a spatial join to add the specified column from vector data
  points_sf <- points_sf %>%
    st_join(vector_data %>% dplyr::select(geometry, all_of(join_column)))
  
  # Clean up temporary variables
  rm(vector_data)
  
  # Print completion message
  print(paste(vector_layer, 'Vector data extraction completed.'))
  
  return(points_sf)
}

#' Prepare Well Data for Anthropogenic Impact Calculations
#'
#' This function processes well location data by attaching hydrologic unit codes (HUC) and elevation
#' data, preparing it for further analysis of anthropogenic impacts.
#'
#' @param M_stn sf object. A spatial dataset of well locations.
#'
#' @return An `sf` object with the following additional attributes:
#' - `huc12`: Hydrologic unit code (HUC12) assigned via spatial join.
#' - `elevation`: Extracted elevation values from a raster dataset.
#'
#' @details
#' - The function first reprojects `M_stn` to match the CRS of the `huc12` dataset.
#' - Performs a spatial join to attach `huc12` codes to well locations.
#' - Calls `extract_raster_values()` to assign elevation data from `"Covar/Elevation_US.tif"`.
#' - Prints a message upon successful preprocessing.
#'
#' @examples
#' \dontrun{
#' # Load required libraries
#' library(sf)
#'
#' # Load sample well locations (ensure it's an sf object)
#' wells_sf <- st_read("path/to/well_locations.geojson")
#'
#' # Prepare well data for impact assessment
#' wells_sf <- prepare_well_impact(wells_sf)
#' }
#'
#' @export
prepare_well_impact <- function(M_stn) {
  ## Prepare well data for anthropogenic impact calculations
  # Attach HUC data to sites
  M_stn <- st_transform(M_stn, st_crs(huc12))
  M_stn <- st_join(M_stn, huc12["huc12"])
  
  # Attach elevation data to sites
  M_stn <- extract_raster_values(raster_path = "Covar/Elevation_US.tif",
                                 points_sf = M_stn,
                                 column_name = "elevation")
  print('position data preprocessed for impact calculation')
  return(M_stn)
}

# Function to calculate impact for a given well
calc.TRI.impact <- function(xrow, df, SEMS = FALSE) {
  # for each well, select sources in the same HUC12
  xlist = xrow %>% st_drop_geometry()
  ptsrc <- df %>%
    filter(elevation > xlist['elevation'], huc12 == as.character(xlist['huc12']))
  impact = 0
  if (nrow(ptsrc) > 0) {
    for (j in 1:nrow(ptsrc)) {
      # calculate the distance between the well and the contamination site
      dist = distm(st_coordinates(xrow), st_coordinates(ptsrc[j, ]), fun =
                     distHaversine) / 1000
      # turn matrix to a number
      dist = as.numeric(dist)
      if (SEMS) {
        # SEMS site release is a constant number of 1
        release = 1
      } else {
        # for TRI site, find/calculate release volume for that well and year
        release = ptsrc[j, ] %>%
          st_drop_geometry() %>%
          dplyr::select(contains(as.character(xlist['year']))) %>%
          as.numeric()# dplyr; if this is slow data.table is faster
      }
      # calculate the impact at the well for the given contamination site
      w = release / exp(dist)
      # sum impacts across contamination sites associated with a well
      impact = impact + w
    } # else, there is no impact so leave it at 0
  }
  xrow <- xrow %>% mutate(TRI.impact = ifelse(length(impact) == 0, 0, impact))
  return(xrow)
}

# Function to extract TRI data and add it to the dataframe
extract_TRI_impact <- function(M_stn, metal.code) {
  ## 19) TRI onsite releases
  # A) prepare site data
  # read in dataframe of cumulative TRI impacts (see TRI_US_Sums.R)
  if (!file.exists(paste0('CoVar/TRI/TRI.', metal.code, '.onsite.cumulativeSums.csv'))) {
    return(M_stn)
  } else {
    TRI.onsite <- read.csv(paste0(
      'CoVar/TRI/TRI.',
      metal.code,
      '.onsite.cumulativeSums.csv'
    )) %>%
      st_as_sf(coords = c("LONGITUDE", "LATITUDE"),
               crs = 4326) %>%
      prepare_well_impact()
    
    overlapping_huc12 <- unique(intersect(M_stn$huc12, TRI.onsite$huc12))
    #simplify TRI data to the ones that are in the same huc12 as the well data
    TRI.onsite <- subset(TRI.onsite, huc12 %in% overlapping_huc12)
    
    # B) calculate impact
    # create list of wells to iterate over
    wells.list <- M_stn %>%
      filter(huc12 %in% overlapping_huc12)
    
    if (nrow(wells.list) > 0) {
      wells.list <- wells.list %>%
        mutate(year = as.integer(format(date, "%Y"))) %>%
        dplyr::select(location.id, huc12, elevation, year) %>% # Keep relevant columns
        rowwise() %>%
        group_split() %>% # Split rows into a list
        # apply calc.TRI.impact to each element in wells.list
        lapply(function(x)
          calc.TRI.impact(x, TRI.onsite)) %>%
        # convert back to a data frame
        bind_rows()
      
      # Add impacts as a column to M_stn
      M_stn <- M_stn %>%
        left_join(
          wells.list %>% st_drop_geometry() %>% dplyr::select("location.id", TRI.impact),
          by = "location.id"
        ) %>%
        dplyr::rename(TRI.total.impact = TRI.impact) %>%
        # fill na values in TRI.total.impact with 0
        mutate(TRI.total.impact = ifelse(is.na(TRI.total.impact), 0, TRI.total.impact))
    } else{
      M_stn <- M_stn %>%
        mutate(TRI.total.impact = 0)
    }
    print('TRI onsite impacts completed')
    
    # 20) TRI water releases
    # A) prepare site data
    # read in dataframe
    TRI.water <- read.csv(paste0(
      'CoVar/TRI/TRI.',
      metal.code,
      '.water.cumulativeSums.csv'
    )) %>%
      st_as_sf(coords = c("LONGITUDE", "LATITUDE"),
               crs = 4326) %>%
      prepare_well_impact()
    
    overlapping_huc12 <- unique(intersect(M_stn$huc12, TRI.water$huc12))
    #simplify TRI data to the ones that are in the same huc12 as the well data
    TRI.water <- subset(TRI.water, huc12 %in% overlapping_huc12)
    
    # B) calculate impact
    # create list of wells to iterate over
    wells.list <- M_stn %>%
      filter(huc12 %in% overlapping_huc12)
    
    if (nrow(wells.list) > 0) {
      wells.list <- wells.list %>%
        mutate(year = as.integer(format(date, "%Y"))) %>%
        dplyr::select(location.id, huc12, elevation, year) %>% # Keep relevant columns
        rowwise() %>%
        group_split() %>% # Split rows into a list
        # apply calc.TRI.impact to each element in wells.list
        lapply(function(x)
          calc.TRI.impact(x, TRI.water)) %>%
        # convert back to a data frame
        bind_rows()
      
      # Add impacts as a column to M_stn
      M_stn <- M_stn %>%
        left_join(
          wells.list %>% st_drop_geometry() %>% dplyr::select("location.id", TRI.impact),
          by = "location.id"
        ) %>%
        dplyr::rename(TRI.water.impact = TRI.impact) %>%
        # fill na values in TRI.total.impact with 0
        mutate(TRI.water.impact = ifelse(is.na(TRI.water.impact), 0, TRI.water.impact))
    } else{
      M_stn <- M_stn %>%
        mutate(TRI.water.impact = 0)
    }
    
    print('TRI water impacts completed')
    
    return(M_stn)
  }
}

# Function to extract SEMS data and add it to the dataframe
extract_SEMS_impact <- function(M_stn, metal.code) {
  # 21) SEMS Impact
  # A) Attach HUC and elevation data to sources
  # read in and clean data
  SEMS <- read.csv(paste0('CoVar/SEMS/SEMS_', metal.code, '_cleaned.csv')) %>%
    mutate(NPL.STATUS = as.character(NPL.STATUS)) %>%
    mutate(NPL.STATUS = if_else(NPL.STATUS == 'FALSE', 'F', NPL.STATUS)) %>%
    filter(NPL.STATUS == 'F') %>%
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
    prepare_well_impact()
  
  overlapping_huc12 <- unique(intersect(M_stn$huc12, SEMS$huc12))
  #simplify SEMS data to the ones that are in the same huc12 as the well data
  SEMS <- subset(SEMS, huc12 %in% overlapping_huc12)
  
  # B) Calculate impact for each well
  # create list of wells to iterate over
  wells.list <- M_stn %>%
    filter(huc12 %in% overlapping_huc12)
  if (nrow(wells.list) > 0) {
    wells.list <- wells.list %>%
      mutate(year = as.integer(format(date, "%Y"))) %>%
      dplyr::select(location.id, huc12, elevation, year) %>% # Keep relevant columns
      rowwise() %>%
      group_split() %>% # Split rows into a list
      # apply calc.TRI.impact to each element in wells.list
      lapply(function(x)
        calc.TRI.impact(x, SEMS, SEMS = TRUE)) %>%
      # convert back to a data frame
      bind_rows()
    
    # Add impacts as a column to M_stn
    M_stn <- M_stn %>%
      left_join(
        wells.list %>% st_drop_geometry() %>% dplyr::select("location.id", TRI.impact),
        by = "location.id"
      ) %>%
      dplyr::rename(SEMS = TRI.impact) %>%
      # fill na values in TRI.total.impact with 0
      mutate(SEMS = ifelse(is.na(SEMS), 0, SEMS))
  } else{
    M_stn <- M_stn %>%
      mutate(SEMS = 0)
  }
  
  return(M_stn)
  
}

# Helper functions for model performance evaluation
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

#' Calculate Top Predictors by Mean Absolute SHAP Value
#'
#' This function takes a shapviz object and computes the mean absolute SHAP value 
#' for each predictor. It returns the names of the top 10 predictors 
#' with the highest mean absolute SHAP values.
#'
#' @param shap A `shapviz` object generated by the `shapviz()` function.
#' @return A character vector of the top 10 predictor names ranked by importance.
#' @import dplyr
#' @export
calc_meanabs_shap <- function(shap) {
  library(dplyr)
  
  shap_absmean <- shap$S %>%
    as.data.frame() %>%
    summarise(across(everything(), ~ mean(abs(.)))) %>%
    pivot_longer(everything(), names_to = "predictor", values_to = "SHAP") %>%
    arrange(desc(SHAP)) %>%
    slice_head(n = 10) %>%
    pull(predictor)
  
  return(shap_absmean)
}

get_important_predictors <- function(metal.code){
  important_predictors <- case_when(
    metal.code == 'As' ~ c('ph', 'ppt', 'C_Sr', 'rech', 'C_Hornbl', 'C_As', 
                           'C_Sb', 'C_Na', 'caco3_kg_sq_m', 'AQPERMNEW', 'detect.limit'),
    metal.code == 'Cd' ~ c("SEMS", "TRI.water.impact", "C_Cd", "C_Pb", "HLR", "no3_pub", 
                           "C_Calcite", "C_Ca", "om_kg_sq_m", 'VALUE_21_pct', "detect.limit"),
    metal.code == 'Li' ~ c('ph', 'rech', 'ppt', 'no3_dom', 'caco3_kg_sq_m', 'C_Kaolinit', 'cec', 
                           'C_Tot_Plag', 'C_Cd', 'bfi', 'detect.limit'),
    metal.code == 'Mn' ~ c('no3_dom', 'unsatWC', 'bfi', 'tmean', 'no3_pub', 'soiloc', 
                           'C_Hornbl', 'AQPERMNEW', 'Hydrate_N', 'C_Na', 'detect.limit'),
    metal.code == 'Sr' ~ c('ph', 'caco3_kg_sq_m', 'bfi', 'rech', 'tmean', 'C_Cd', 'ppt', 
                           'C_Ca', 'ec', 'C_Tot_Plag', 'detect.limit'),
    TRUE ~ NA_character_
  )
  return(important_predictors)
}

print('helper functions loaded.')