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
  "gstat",
  "lattice",
  "NADA",
  "plyr",
  "parallel",
  "raster",
  "readxl",
  "recipes",
  "sf",
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
  "xgboost"
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

print('helper functions loaded.')