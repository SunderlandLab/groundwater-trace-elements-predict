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
source(here::here('R/0_helper_fct.R'))

setwd(here::here("data"))
mapUSm <- st_read("CoVar/US48/US_48states.shp") # load projected map of the US (in m)
huc12 <- st_read(dsn = "CoVar/WBD_National_GDB/WBD_National_GDB.gdb", layer = "WBDHU12") |> st_make_valid()
elevation <- raster("Covar/Elevation_US.tif")

#### Define Model Version ####
metal.codes <- c("As", "Cd", "Li", "Mn", "Sr")
metals <- c("Arsenic", "Cadmium",  "Lithium", "Manganese", "Strontium")
names(metals) <- metal.codes

#### Main function
predictor_extract_function <- function(metal.code){
  print(paste('Start extracting predictor values for', metal.code))
  
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

  ## 19) Anthropogenic impacts
  
  M_stn <- M_stn %>%
    prepare_well_impact() 
  
  M_stn <- M_stn %>%
    extract_TRI_impact(metal.code) %>%
    extract_SEMS_impact(metal.code)

  # Drop huc12 and elevation as columns
  M_stn <-  M_stn |>
    dplyr::select(-c(huc12,elevation)) 


  ## Step 4. Post-Process Extracted Data --------------------------------------------------------------------------------------------------------------------
  # Convert Categorical Variables to Factors
  M_stn <- M_stn %>%
    mutate(across(c(landcover_500m, AQ_CODE, ROCK_TYPE, rocktype, 
                    soilorder_int, str_int, weg_int, hydgrp_int, 
                    drainage_class_int, HLR), as.factor))
  
  # save the data file out
  saveRDS(M_stn, paste0("R_Output/",metal.code,"_RawPredictors.rds"))
  print("Covariables Extracted")
}

# vectorize across all five metals
purrr::map(metal.codes[-c(1,2)], predictor_extract_function)
