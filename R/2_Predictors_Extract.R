##################################################################################
#File name: Predictors_Extract.R
#Author: Jennifer Sun, Jonas LaPier, Cindy Hu
#Date: Jan 2025
#Purpose: Attach explanatory variable data (i.e. predictors) to well locations
# R v4.4.2
##################################################################################

packages <- c('akima','beepr','chron','foreign','geosphere','ggplot2','gstat', 'lattice','NADA','plyr', 
              'parallel','raster', 'readxl','sf','stringr','terra','tiff','VIM','dplyr')
lapply(packages, library, character.only=TRUE)

setwd(here::here("data"))
mapUSm <- st_read("CoVar/US48/US_48states.shp") # load projected map of the US (in m)

#### Define Model Version ####
metal.codes <- c("As", "Cd", "Li", "Mn", "Sr")
metals <- c("Arsenic", "Cadmium",  "Lithium", "Manganese", "Strontium")
names(metals) <- metal.codes

#### Helper functions
extract_raster_values <- function(raster_path, points_sf, column_name, method = 'bilinear', buffer_on = FALSE) {
  # Load the raster
  raster_layer <- raster(raster_path)
  # Ensure points_sf is an sf object
  if (!inherits(points_sf, "sf")) stop("points_sf must be an sf object.")
  # if raster is missing CRS, assign the CRS to it
  if(is.na(raster::crs(raster_layer))){
    raster::crs(raster_layer) <- CRS('+init=epsg:5070')
  }
  
  # Reproject points to match the raster CRS
  points_sf <- st_transform(points_sf, crs = st_crs(raster_layer))
  
  if(!buffer_on){
    # Extract raster values to the points
    points_sf[[column_name]] <- raster::extract(raster_layer, points_sf, method = method)
  } else {
    # Extract raster values to the points with a buffer
    points_sf[[column_name]] <- raster::extract(raster_layer, points_sf, buffer = 500, fun = mean)
  }
  
  # Clean up
  rm(raster_layer)
  gc()  # Optional: Call garbage collection
  
  # Notify completion
  message(paste(column_name, "completed"))
  
  return(points_sf)
}

extract_vector_values <- function(points_sf, vector_dsn, vector_layer, join_column) {
  # Read the vector data as sf
  vector_data <- st_read(dsn = vector_dsn, layer = vector_layer)
  # if vector is missing CRS, assign the CRS of mapUSm to it
  if(is.na(st_crs(vector_data))){
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

#### Main function
predictor_extract_function <- function(metal.code){
  ##  Step 1. Load in WQP well location data and well measurements ------------------------------------------------------------------------------- 
  M_stn = readRDS(paste0("Data_Files/", metal.code,"_Combined_WQP.rds")) # read in preprocessed data
  M_stn <- st_as_sf(M_stn, coords = c("longitude", "latitude"), crs = st_crs(mapUSm))# setting initial crs

  ##  Step 2. Environmental predictors: extract values to wells--------------------------------------------------------------------------
  # 1) Precipitation (mm)
  M_stn <- extract_raster_values(
    raster_path = 'CoVar/PRISM_ppt_30yr_normal_800mM2_annual_asc/PRISM_ppt_30yr_normal_800mM2_annual_asc.asc',
    points_sf = M_stn,
    column_name = "ppt"
  )

  # 2) Annual Mean Temp (F)
  M_stn <- extract_raster_values(
    raster_path = 'CoVar/PRISM_tmean_30yr_normal_800mM2_annual_asc/PRISM_tmean_30yr_normal_800mM2_annual_asc.asc',
    points_sf = M_stn,
    column_name = "tmean"
  )

  # 3) Groundwater flow parameters from Reitz et al. 2017 (800 m resolution)
  M_stn <- extract_raster_values(
    raster_path = 'CoVar/reitz_hydrology/QuickFlow_0013/RO_0013.tif',
    points_sf = M_stn,
    column_name = "quickflow"
  )
  
  # 4) Recharge
  M_stn <- extract_raster_values(
    raster_path = 'CoVar/reitz_hydrology/EffRecharge_0013_v2/0013/RC_eff_0013.tif',
    points_sf = M_stn,
    column_name = "effrech"
  )
  
  # 5) ET
  M_stn <- extract_raster_values(
    raster_path = 'CoVar/reitz_hydrology/ET_0013/ET_0013.tif',
    points_sf = M_stn,
    column_name = "ET"
  )

  # 3b) Original recharge from Wolock et al. for comparison
  M_stn <- extract_raster_values(
    raster_path = 'CoVar/rech48grd/rech48grd/w001001x.adf',
    points_sf = M_stn,
    column_name = "rech"
  )

  # 4) Soil Organic Content, (tons/Hectare)
  M_stn <- extract_raster_values(
    raster_path = 'CoVar/HWSDa_OC_Dens_30SEC_US.asc',
    points_sf = M_stn,
    column_name = "soiloc"
  )

  # 5) Soil Geochemistry, C-horizon
  print("started geochem")
  # We iterate across a list of soil geochemical properties
  chemNameList<-c('C_As', 'C_C_Inorg','C_C_Org', 'C_Hornbl','C_C_Tot','C_Na','C_Tot_Clay',
                  'C_P','C_U','C_Mo','C_Ca','C_Sr','C_Fe','C_Mn','C_V','C_Ti','C_K','C_Be',
                  'C_Tot_Plag','C_Tot_10A','C_Tot_14A','C_Tot_K_fs','C_Ni','C_Cr','C_Tot_Flds',
                  'C_Hg','C_Kaolinit','C_Pb','C_Cd','C_Li','C_Sb','C_Mg','C_Rb','C_Gypsum','C_Aragon','C_Calcite')
  
  geochem<-read.table('CoVar/soil_geochemistry/Appendix_4b_Chorizon_18Sept2013.txt', header = T, # first row is unit
                      sep = '\t',stringsAsFactors = F) 
  # Drop the first observation (row)
  geochem <- geochem[-1, ]
  # only stations with viable Latitude and Longitude are useful to us, so we exclude any station that is missing the coordinate.
  # convert to sf object
  geochem_sf <- geochem |>
    mutate(across(all_of(chemNameList), as.numeric)) |>  # Convert all columns to numeric
    mutate(Longitude = as.numeric(as.character(Longitude)),
           Latitude = as.numeric(as.character(Latitude))) %>%
    filter(!is.na(Longitude) & !is.na(Latitude)) |>
    st_as_sf(coords = c("Longitude", "Latitude"),  # Specify the columns for coordinates
    crs = 4326)  # Set CRS to WGS84 (EPSG:4326))
  M_stn <- st_transform(M_stn, st_crs(geochem_sf)) # extract crs and set well projection
  
  fx = function(chemName) {
    system(paste("echo '",chemName,"'"))
    geochemOut <- as(geochem_sf, "Spatial")
    # Create an empty raster
    r <- raster(extent(geochemOut), res = 0.5, crs = proj4string(geochemOut))
    # Rasterize using raster
    r <- rasterize(geochemOut, r, field = chemName, fun = mean)
    res <- terra::extract(r, M_stn, method='bilinear')
    return(res)  # extract value for each variable
    }
  list = mclapply(chemNameList, fx, mc.cores = 6) # runs parallel processing
  names(list) = chemNameList # double check that no columns were skipped or returned errors
  M_stn = bind_cols(M_stn, data.frame(list))
  rm(geochem)
  print('geochem completed')
  

  # 6) Hydrologic Landscape Regions of the US
  hlrus<-st_read(dsn="CoVar/hlrshape/hlrus.shp")
  hlrus <- st_transform(hlrus, st_crs(mapUSm))
  M_stn <- st_transform(M_stn, st_crs(mapUSm)) # extract crs and set well projection
  M_stn<-st_join(M_stn, hlrus[c("SLOPE","SAND","RELIEF","PFLATTOT","PFLATLOW","PFLATUP","HLR","AQPERMNEW")])
  rm(hlrus)
  print('hlrus completed')

  # 7) Base Flow Index (integer)
  M_stn <- extract_raster_values(raster_path ='CoVar/bfi48grd/bfi48grd/w001001x.adf',    
                                 points_sf = M_stn,
                                 column_name = "bfi")
  # 8) Soil Properties from STATSGO
  statsgo<-st_read(dsn="CoVar/statsgo/muid_poly_with_data.shp")
  M_stn <- st_transform(M_stn, st_crs(statsgo)) # extract crs and set well projection
  proplist = c('PERML','PERMH','AWCL','AWCH','BDL','BDH','OML','OMH','SLOPEL','SLOPEH','WTDEPL',
             'WTDEPH','ROCKDEPL','ROCKDEPH','KFACT','TFACT','WEG')
  M_stn<-st_join(M_stn, statsgo[proplist])
  M_stn <- M_stn %>%
    mutate(
      PERMmean = (PERML + PERMH) / 2,
      AWCmean = (AWCL + AWCH) / 2,
      BDmean = (BDL + BDH) / 2,
      OMmean = (OML + OMH) / 2,
      SLOPEmean = (SLOPEL + SLOPEH) / 2,
      WTDEPmean = (WTDEPL + WTDEPH) / 2,
      ROCKDEPmean = (ROCKDEPL + ROCKDEPH) / 2
    )
  rm(statsgo)
  print('statsgo completed')


  # 9) Bedrock Geology: King and Beikman 
  M_stn <- extract_vector_values(M_stn, "CoVar/kbgeology", "kbge", "UNIT") |>
    rename(KB = UNIT)
  print('bedrock geology completed')


  # 10) Surficial Geology
  M_stn <- extract_vector_values(M_stn, "CoVar/Soller_surfgeo/USGS_DS_425_SHAPES", "Surficial_materials", "UNIT_NAME")|>
    rename(surfgeo = UNIT_NAME)
  print('surface geology completed')

  # 10) Generalized Lithology
  rocktype <- data.table::fread("CoVar/Anning_Schweitzer_lithology/geol_poly_rocktypes.csv") %>%
    distinct(UNIT_LINK, rocktype) %>%     # Retain unique UNIT_LINK and rocktype pairs
    group_by(UNIT_LINK) %>%              # Group by UNIT_LINK
    slice_min(order_by = rocktype, n = 1) %>% # Select the first rocktype per group
    ungroup() 
  lith<-st_read("CoVar/Anning_Schweitzer_lithology/geol_poly/geol_poly.shp") |>
    left_join(rocktype, by=c("UNIT_LINK"="UNIT_LINK")) 
  M_stn <- st_transform(M_stn, st_crs(lith)) |>
    st_join(lith['rocktype'])
  rm(lith)
  rm(rocktype)
  print ('generalized lithology completed')

  # 11) Land Cover  
  nlcd <- read.csv('CoVar/nlcd2011/AllUniqueWells_20230306_nlcd2011.csv')
  nlcd_cat <- read.csv('CoVar/nlcd2011/AllUniqueWells_20230306_nlcd2011_bycat.csv')
  
  # join majority landcover type
  nlcd <- nlcd %>% 
    dplyr::select(location_id, MAJORITY) %>% 
    dplyr::rename('landcover_500m' = 'MAJORITY','location.id'='location_id')
  M_stn <- dplyr::left_join(M_stn, nlcd, by='location.id')

  # join % of each landcover type
  nlcd_cat <- nlcd_cat %>% 
    dplyr::select(-c(OBJECTID)) %>% 
    dplyr::mutate(sum = rowSums(across(where(is.numeric)), na.rm=TRUE)) 
  nlcd_pct <- nlcd_cat %>% 
    dplyr::mutate_at(vars(VALUE_0:sum), list(pct=~./sum)) %>% 
    dplyr::select(LOCATION_ID, contains('pct'), -sum_pct) %>%
    dplyr::rename('location.id'='LOCATION_ID')
  M_stn <- dplyr::left_join(M_stn, nlcd_pct, by='location.id')
  rm(nlcd, nlcd_cat, nlcd_pct)

  # # 12) County level drainage data, divide best guess drainage acre by total area, percent
  # drainage <- raster('CoVar/usgs_tiledrainage/SubsurfaceDrainExtentUS_90s.tif')
  # # original resolution is 30 meters, too fine, aggregate by 20 times
  # drainage <- aggregate(drainage, fact = 20, fun = mean)
  # writeRaster(drainage, "CoVar/usgs_tiledrainage/SubsurfaceDrainExtentUS_90s_agg.tif", overwrite = TRUE)

  M_stn <- extract_raster_values('CoVar/usgs_tiledrainage/SubsurfaceDrainExtentUS_90s_agg.tif',
                                 M_stn,
                                 'drainage',
                                 method = 'bilinear')

  # 13) stream density
  M_stn <- extract_raster_values('CoVar/stmdenhuc12/w001001.adf', 
                                 M_stn,
                                 'stn_den',
                                 method = 'bilinear')

  # 14) Surficial groundwater flow parameters (MODFLOW 6) (250 m resolution)
  M_stn <- extract_raster_values('CoVar/modflow6_surfgw/Output_CONUS_trans_dtw/conus_MF6_SS_Unconfined_250_dtw.tif', 
                                 M_stn,
                                 'dtw',
                                 buffer_on = TRUE)
  M_stn <- extract_raster_values('CoVar/modflow6_surfgw/Output_CONUS_trans_dtw/conus_MF6_SS_Unconfined_250_trans.tif', 
                                 M_stn,
                                 'trans',
                                 buffer_on = TRUE)
  M_stn <- extract_raster_values('CoVar/modflow6_surfgw/Output_CONUS_unsat_traveltime/conus_MF6_SS_Unconfined_250_tt_total.tif', 
                                 M_stn,
                                 'unsatTT',
                                 buffer_on = TRUE)
  M_stn <- extract_raster_values('CoVar/modflow6_surfgw/Output_CONUS_unsat_watercontent/conus_MF6_SS_Unconfined_250_wc_avg.tif', 
                                 M_stn,
                                 'unsatWC',
                                 buffer_on = TRUE)
  print('groundwater flow parameters completed')


  # 15) Aquifer
  M_stn <- extract_vector_values(M_stn, 
                        "CoVar/aquifrp025_nt00003", 
                        "aquifrp025",
                        c("AQ_CODE", "ROCK_TYPE"))

  # 16) Soil Properties Rasters (UCDAVIS - combined SSURGO and STATSGO 800m rasters)
  # http://soilmap2-1.lawr.ucdavis.edu/soil-properties/download.php
  for (filename in list.files(path="CoVar/soil_properties")) {
    if (grepl("int",filename)){
      M_stn <- extract_raster_values(paste0("CoVar/soil_properties/",filename), 
                                     M_stn, 
                                     tools::file_path_sans_ext(filename),
                                     method = 'simple')
    } else {
      M_stn <- extract_raster_values(paste0("CoVar/soil_properties/",filename), 
                                     M_stn, 
                                     tools::file_path_sans_ext(filename))
    }
  }


  # 17) Soil Property Rasters Part II (USGS SSURGO area- and depth-weighted 90m rasters)
  # read in parameter values for each well as calculated in ArcGIS & merge based on well name / location ID 
  for (filename in list.files(path="CoVar/SSURGO90m/well_data")) {
    tbl = read_excel(paste0('CoVar/SSURGO90m/well_data/',filename))
    name = str_replace(filename, 'AllUniqueWells_mean_', '') %>% str_remove('.xlsx')
    tbl = tbl %>% dplyr::select('location_id','MEAN') %>% dplyr::rename(!!name := 'MEAN','location.id'='location_id')
    M_stn <- dplyr::left_join(M_stn, tbl, by='location.id')
    print(paste(name,'completed'))
  } 
  
  rm(tbl)


# 18) Predicted Groundwater NO3 (Ransom et al. 2021, 1 km resolution)
# groundwater no3 predicted at domestic well depths
  M_stn <- extract_raster_values(
    raster_path = 'CoVar/USGS_NO3/no3_doms.asc',
    points_sf = M_stn,
    column_name = "no3_dom"
  )
  # groundwater no3 predicted at public well depths
  M_stn <- extract_raster_values(
    raster_path = 'CoVar/USGS_NO3/no3_pubs.asc',
    points_sf = M_stn,
    column_name = "no3_pub"
  )

  ###  Step 3. Anthropogenic impacts: extract values to wells --------------------------------------------------------------------------------------
  ## Read in functions for each calculation
  # Function to calculate impact for a given well
  calc.TRI.impact <- function(xrow, df) { 
    # for each well, select sources in the same HUC12
    xlist = unlist(xrow) # longitude, latitude, huc12, elevation, year
    ptsrc <- subset(df, (huc12==xlist[3]) & (elevation > xlist[4])) 
    impact = 0
    if (nrow(ptsrc)>0) {
      #print('there are impacts')
      for (j in 1:nrow(ptsrc)) {
        # calculate the distance between the well and the contamination site
        dist = distm(c(xlist[1],xlist[2]), coordinates(ptsrc[j,]), fun=distHaversine)/1000
        # find/calculate release volume for that well and year
        release = ptsrc[j,]@data %>% dplyr::select(contains(as.character(xlist[5]))) # dplyr; if this is slow data.table is faster
        # calculate the impact at the well for the given contamination site
        w = release/exp(dist)
        # sum impacts across contamination sites associated with a well
        impact = impact + w
        #print(impact)
      } # else, there is no impact so leave it at 0
    }
    return(impact)
  }
  
  # Function to calculate impact for a given well
  calc.SEMS.impact <- function(xrow, df) {
    # for each well, select sources in the same HUC12
    xlist = unlist(xrow) # longitude, latitude, huc12, elevation
    ptsrc <- subset(df, (huc12==xlist[3]) & (elevation > xlist[4]))
    impact = 0
    if (nrow(ptsrc)>0) {
      for (j in 1:nrow(ptsrc)) {
        dist = distm(c(xlist[1],xlist[2]), coordinates(ptsrc[j,]), fun=distHaversine)/1000
        w = 1/exp(dist)
        impact = impact + w
      } # else, there is no impact so leave it at 0
    }
    return(impact)
  }

  # Function to prepare well data for anthropogenic impact calculations 
  gdb_path <- "CoVar/WBD_National_GDB/WBD_National_GDB.gdb"
  huc12 <- st_read(dsn = gdb_path, layer = "WBDHU12") |>
    st_make_valid()
  elevation <- raster("Covar/Elevation_US.tif")

  prepare_well_impact <- function(M_stn) {
    ## Prepare well data for anthropogenic impact calculations
    # Attach HUC data to sites
    M_stn <- st_transform(M_stn, st_crs(huc12))
    M_stn <- st_join(M_stn, huc12["huc12"]) 
    
    # Attach elevation data to sites
    M_stn <- extract_raster_values(
      raster_path = "Covar/Elevation_US.tif",
      points_sf = M_stn,
      column_name = "elevation"
    )
    print('well position data preprocessed')
    return(M_stn)
  }

# Function to extract TRI data and add it to the dataframe
extract_TRI_impact <- function(M_stn, metal.code) {
  ## 19) TRI onsite releases
  # A) prepare site data
  # read in dataframe of cumulative TRI impacts (see TRI_US_Sums.R)
  print(metal.code)
  TRI.onsite <- read.csv(paste0('TRI_US/TRI.', metal.code,'.onsite.cumulativeSums.csv'))
  TRI.onsite <- st_as_sf(TRI.onsite,
                         coords = c("LONGITUDE", "LATITUDE"),
                         crs = 4326)
  # attach huc12 to TRI sites
  TRI.onsite <- st_transform(TRI.onsite, st_crs(huc12))
  TRI.onsite <- st_join(TRI.onsite, huc12['huc12'])
  # attach elevations to TRI sites
  M_stn <- extract_raster_values(
    raster_path = "Covar/Elevation_US.tif",
    points_sf = TRI.onsite,
    column_name = "elevation"
  )
  
  # B) calculate impact
  # create list of wells to iterate over
  wells.list = split(c(coordinates(M_stn), as.numeric(M_stn$huc12), M_stn$elevation, as.integer(format(M_stn$date, format="%Y"))), seq(nrow(M_stn)))
  # Iterate over each well
  impacts = lapply(wells.list, function(x) calc.TRI.impact(x, TRI.onsite))
  # Add impacts as a column to M_stn
  M_stn@data$TRI.total.impact <- unlist(impacts) 
  print('TRI onsite impacts completed')
  
  #ggplot(M_stn@data) + geom_density(aes(log10(M_stn$TRI.onsite)))
  
  # 20) TRI water releases
  # A) prepare site data
  # read in dataframe 
  TRI.water <- read.csv(paste0('TRI_US/TRI.', metal.code,'.water.cumulativeSums.csv'))
  coordinates(TRI.water) = ~LONGITUDE + LATITUDE
  # attach huc12 to TRI sites
  proj4string(TRI.water) <- crs(huc12)
  TRI.water$huc12 <- unlist(over(TRI.water, huc12[,'huc12']))
  # attach elevations to TRI sites
  TRI.water <- spTransform(TRI.water, crs(elevation)) # extract crs and set well projection
  TRI.water$elevation <- raster::extract(elevation, TRI.water, method='bilinear')
  print('TRI water site position data preprocessed')
  
  # B) calculate impact
  # create list of wells to iterate over
  wells.list = split(c(coordinates(M_stn), as.numeric(M_stn$huc12), M_stn$elevation, as.integer(format(M_stn$date, format="%Y"))), seq(nrow(M_stn)))
  # Iterate over each well
  impacts = lapply(wells.list, function(x) calc.TRI.impact(x, TRI.water))
  # Add impacts as a column to M_stn
  M_stn@data$TRI.water.impact <- unlist(impacts) 
  print('TRI water impacts completed')
  
  return(M_stn)

}

# Function to extract SEMS data and add it to the dataframe
extract_SEMS_impact <- function(M_stn, metal.code) {
  # 21) SEMS Impact 
  # A) Attach HUC and elevation data to sources
  # read in and clean data
  SEMS <- read.csv(paste0('CoVar/NPL/SEMS_',metal.code,'_cleaned.csv'))
  SEMS$NPL.STATUS[SEMS$NPL.STATUS == 'FALSE'] = 'F' # when is this needed? 
  SEMS <- subset(SEMS, NPL.STATUS=='F') # F means "final"
  coordinates(SEMS) = ~Longitude + Latitude 
  SEMS.coords <- coordinates(SEMS)
  # attach huc12
  proj4string(SEMS) <- crs(huc12)
  SEMS$huc12 <- unlist(over(SEMS, huc12[,'huc12']))
  # attach elevation
  SEMS <- spTransform(SEMS, crs(elevation)) # extract crs and set well projection
  SEMS$elevation <- raster::extract(elevation, SEMS, method='bilinear')
  print('SEMS site position data preprocessed')

  # B) Calculate impact for each well
  # create list of wells to iterate over
  wells.list = split(c(coordinates(M_stn), as.numeric(M_stn$huc12), M_stn$elevation), seq(nrow(M_stn)))
  # Iterate over each well
  impacts = lapply(wells.list, function(x) calc.SEMS.impact(x, SEMS))
  # Add impacts as a column to M_stn
  M_stn@data$SEMS <- as.numeric(impacts) # count: for Cd, 906 wells have some type of SEMS impact 
  print('SEMS completed')
  
  return(M_stn)
  
}

## Calculate impact
M_stn = prepare_well_impact(M_stn)
M_stn = extract_TRI_impact(M_stn, metal.code)
M_stn = extract_SEMS_impact(M_stn, metal.code)

# Drop huc12 and elevation as columns
M_stn <-  M_stn |>
  dplyr::select(-c(huc12,elevation)) 


## Step 4. Post-Process Extracted Data --------------------------------------------------------------------------------------------------------------------
# Convert Categorical Variables to Factors
# For some reason, this only works with the intermediate variable "fac"
fac = as.factor(M_stn$landcover_500m)
M_stn$landcover_500m = fac
fac = as.factor(M_stn$aquifer)
M_stn$aquifer = fac
fac = as.factor(M_stn$aq_rocktype)
M_stn$aq_rocktype = fac
fac = as.factor(M_stn$lith)
M_stn$lith = fac
fac = as.factor(M_stn$soilorder_int)
M_stn$soilorder_int = fac
fac = as.factor(M_stn$str_int)
M_stn$str_int = fac
fac = as.factor(M_stn$weg_int)
M_stn$weg_int = fac
fac = as.factor(M_stn$hydgrp_int)
M_stn$hydgrp_int = fac
fac = as.factor(M_stn$drainage_class_int)
M_stn$drainage_class_int = fac
fac = as.factor(M_stn$dHLR)
M_stn$HLR = fac

# save the data file out
saveRDS(M_stn, paste0("R_Output/",folder.name,metal.code,"_RawPredictors.rds"))
print("Covariables Extracted")
}

## Step 5. Refine the Data Frame with Data Imputation -----------------------------------------------------------------------------------------------------
# This section fills in missing values with imputed values using a K-nearest neighbors imputation

# Load Data
master = readRDS(paste0("R_Output/",folder.name, metal.code,"_RawPredictors.rds"))

# Calculate how many missing values there are per predictor
is.missing <- function(df) {
  missingness = matrix(ncol=3, nrow = length(names(df))) 
  for (i in 1:length(names(df))) {
    Var = names(df)[i]
    missing <- sum(is.na(df[,Var]))
    missingpct <- round(missing/ length(df[,Var]),5)
    missingness[i,] = c(Var, missing, missingpct)
  }
  print(missingness)
}

missingness = as.data.frame(is.missing(master@data))
colnames(missingness) = c(c('Variable','MissingRecords','MissingPct'))
missingness[order(missingness$MissingPct, decreasing=TRUE),]

## Remove Extra Columns
# master<-subset(master,select=-c(conc,location.id,insideUS,TRI.As.Ct,TRI.Pb.Ct,TRI.Cd.Ct)) 
master<-subset(master,select=-c(insideUS)) 

# Set missing impacts and drainage to zero
master$TRI.total.impact[is.na(master$TRI.total.impact)] = 0 # no missing TRI data on the map grid, so just skip this
master$TRI.water.impact[is.na(master$TRI.water.impact)] = 0
master$drainage[is.na(master$drainage)] = 0

## Imputation 
# Create a data frame with location information as columns
master_locs = as.data.frame(master@data)
master_locs$long = coordinates(master)[,1] # master$longitude [map grid coordinate is not named longitude]
master_locs$lat = coordinates(master)[,2] # master$latitude [map grid coordinate is not named longitude]
  

# set integer values to numeric
master_locs = master_locs %>% mutate_if(is.integer, as.numeric)

# fill in missing values manually for factor variables - just create a missing category
master_locs$KB[is.na(master_locs$KB)] = 'ms'

lithlevels = levels(master_locs$lith)
lithlevels = c(lithlevels, 'NONE')
levels(master_locs$lith) = lithlevels
master_locs$lith[is.na(master_locs$lith)] = 'NONE'

int_predictors = c('str_int','hydgrp_int','drainage_class_int','soilorder_int','weg_int')
for (x in int_predictors) {
  pred_levels = levels(master_locs[[x]])
  pred_levels = c(pred_levels, 0)
  levels(master_locs[[x]]) = pred_levels
  master_locs[[x]][is.na(master_locs[[x]])] = 0
}

# calculate missingness in remaining predictors
missingness = as.data.frame(is.missing(master_locs))
colnames(missingness) = c(c('Variable','MissingRecords','MissingPct'))
missingness[order(missingness$MissingPct, decreasing=TRUE),]

# use KNN imputation to impute missing values for numeric columns
if (metal.code %in% c('As','Mn','Li','Sr')) { 
  vars = names(dplyr::select(ungroup(master@data), -c('location.id','date','conc','censored.conc','censored','well.depth','data.source','ros.conc','DL.missing','TRI.total.impact','TRI.water.impact','SEMS')))
  extra.remove = list('drainage','C_Aragon','WTDEPAMJ','C_Gypsum','WTDEP_MIN') # do not impute for parameters like > 10% missing
  vars = setdiff(vars, extra.remove)
} else if (metal.code == 'Cd') {
  vars = names(dplyr::select(ungroup(master@data), -c('location.id','date','conc','censored.conc','censored','well.depth','data.source','ros.conc','DL.missing','TRI.total.impact','TRI.water.impact','SEMS')))
  landcover_predictors = M_stn@data %>% select(contains('VALUE'), landcover_500m) %>% colnames()
  vars = setdiff(vars,landcover_predictors) # do not impute landcover values for the 56 wells that are missing this information
  extra.remove = list('drainage','C_Aragon','WTDEPAMJ','C_Gypsum','WTDEP_MIN') # do not impute for parameters like > 10% missing
  vars = setdiff(vars, extra.remove)
} else {
  print ('error in metal code selection')
}

master_locs = VIM::kNN(master_locs, variable=vars, k=5, dist_var=c('long','lat'), impNA=FALSE) 
master@data = master_locs 

# remove the location information from the data
master = master[,-grep('_imp',names(master))] # removes columns with the suffix "imp" which indicates whether the row was imputed or not for that column
# master <- subset(master, select=-c(long,lat))

# Save the refined data set
saveRDS(master, paste0("R_Output/",folder.name,metal.code,"_Model_Ready.rds"))
