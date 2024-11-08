###########################################################################
#File name: All_Data_Prep.R
#Author: Jennifer Sun, Jonas LaPier, Cindy Hu
#Date: March 2023
#Purpose: All steps for WQP well data preparation 
########################################################################### 

#Please note that rgdal will be retired during 2023,
#plan transition to sf/stars/terra functions using GDAL and PROJ

#### Attach Libraries and Set Working Directory ####
packages <- c('akima','assertthat', 'chron','data.table','tidyr','tidyverse','plyr','geosphere','ggplot2', 
              'lattice','NADA','raster','rgdal','rgeos','sp', 'sf', 'stringr','terra','data.table','dplyr')
lapply(packages, library, character.only=TRUE)


setwd(here::here("data"))
mapUSm <- readOGR(dsn="CoVar/US48", layer="US_48states") # load projected map of the US (in m)

#### Define Model Version ####
metal = "Arsenic"
metal.code = "As"
version.number = '7a'
MCL = 10
folder.name = paste0(metal.code,'_',version.number,'/')

## MCL Values ----
# Arsenic -> 10 ug/L
# Cadmium -> 5 ug/L
# Lead -> 15 ug/L (action level for drinking water)
# Manganese -> 50 ug/L (secondary standard)
# Uranium -> 30 ug/L (MCL)
# Strontium -> 4000 ug/L (HBSL)
# Lithium -> 60 ug/L (drinking water HBSL) 

## 1. Load WQP Files -----------------------------------------------------------------------------------------------------------------
#  Read in the station and sample result csv files obtained from the WQP. Locations outside the contiguous US are excluded.

## Load Data
metal.stn = read.csv(paste0("WQP/", metal.code, "-station.csv"))
metal.smpl = read.csv(paste0("WQP/", metal.code, "-result.csv"))

# Merge files based on location identifier
metal.large = metal.smpl |>
  left_join(metal.stn, by=c('MonitoringLocationIdentifier', 'ProviderName'))


## 2. Format WQP data ------------------------------------------------------------------------------------------------------------------
# Rename columns and match up data types and units

metal.df = metal.large |>
  # Keep only important columns (all well depth units are in feet)
  dplyr::select("MonitoringLocationIdentifier", "MonitoringLocationTypeName", "LatitudeMeasure",
                              "LongitudeMeasure","AquiferName", "AquiferTypeName", "WellDepthMeasure.MeasureValue",
                              "ActivityIdentifier","ActivityStartDate","ActivityTypeCode", "ResultMeasureValue","ResultMeasure.MeasureUnitCode",
                              "ResultDetectionConditionText",
                              "DetectionQuantitationLimitMeasure.MeasureValue","DetectionQuantitationLimitMeasure.MeasureUnitCode",
                              "ActivityMediaSubdivisionName",'ResultSampleFractionText','ResultStatusIdentifier','ResultValueTypeName',
                              'ResultLaboratoryCommentText','DetectionQuantitationLimitTypeName','ProviderName') |>
  # Rename columns
  dplyr::rename("location.id" = "MonitoringLocationIdentifier", 
         "location.type" = "MonitoringLocationTypeName",
         "latitude" = "LatitudeMeasure",
         "longitude" = "LongitudeMeasure",
         "aquifer.name" = "AquiferName",
         "aquifer.type" = "AquiferTypeName",
         "well.depth" = "WellDepthMeasure.MeasureValue",
         "activity" = "ActivityIdentifier",
         "date" = "ActivityStartDate",
         "activity.type" = "ActivityTypeCode",
         "conc" = "ResultMeasureValue",
         "conc.unit" = "ResultMeasure.MeasureUnitCode",
         "result.condition" = "ResultDetectionConditionText",
         "detect.limit" = "DetectionQuantitationLimitMeasure.MeasureValue",
         "detect.unit" = "DetectionQuantitationLimitMeasure.MeasureUnitCode",
         "sample.media" = "ActivityMediaSubdivisionName",
         "sample.fraction" = "ResultSampleFractionText",
         "result.status" = "ResultStatusIdentifier",
         "result.type" = "ResultValueTypeName",
         "result.comments" = "ResultLaboratoryCommentText",
         "DL.type" = "DetectionQuantitationLimitTypeName",
         "data.source" = "ProviderName") |>
  # convert conc, well.depth, detect.limit from character to numeric
  mutate(conc = as.numeric(conc),
         well.depth = as.numeric(as.character(well.depth)),
         detect.limit = as.numeric(as.character(detect.limit))) |>
  # add a year column, and convert date to date type
  mutate(date = lubridate::ymd(date),
         year = lubridate::year(date)) |>
  # if conc.unit is mg/l or ppm, convert conc to 1000*conc, then reset conc.unit to ug/L
  mutate(conc = ifelse(conc.unit %in% c("mg/l", "ppm"), 1000*conc, conc),
         conc.unit = ifelse(conc.unit %in% c("mg/l", "ppm"), "ug/l", conc.unit)) |>
  # if detect.unit is mg/l or ppm, convert detect.limit to 1000*detect.limit, then reset detect.unit to ug/L
  mutate(detect.limit = ifelse(detect.unit %in% c("mg/l", "ppm"), 1000*detect.limit, detect.limit),
         detect.unit = ifelse(detect.unit %in% c("mg/l", "ppm"), "ug/l", detect.unit))
  


#### 3. Data Cleaning ---------------------------------------------------------------------------------------------------------------------
metal.df1 <- metal.df |>
  # subset to samples collected in 1990 and later
  filter(year >= 1990) |>
  # Remove rows with other units (i.e. mg/kg, pCi/L)
  filter(conc.unit == "ug/l" | (conc.unit=="" & detect.unit == "ug/l"))  |>
  # Replace NA in conc variable with zero
  mutate(conc = ifelse(is.na(conc), 0, conc)) |>
  # Remove negative conc
  filter(conc >= 0) |>
  # keep rows that match string "Well" in location.type
  filter(grepl("Well", location.type)) |>
  # keep rows where activity.type is Sample or Sample-Routine
  filter(activity.type %in% c("Sample", "Sample-Routine")) |>
  # Remove records that are not labeled as groundwater in sample.media column
  filter(sample.media %in% c("Ground Water", "Groundwater")) |>
  # Remove records that are the wrong sample fraction type
  filter(!sample.fraction %in% c("Suspended", "", "Bed Sediment")) |>
  # remove records in which the concentration AND the detection limit are both missing
  filter(!(conc==0 & (is.na(detect.limit)|detect.limit==0))) |>
  # remove records with result.condition in '*Present','*Present >QL','Detected Not Quantified','Present Above Quantification Limit','Systematic Contamination'
  #???# Original comment says "Remove records with missing detection values (present > QL or detected but no reported value)", not sure if the code matches
  filter(!(result.condition %in% c('*Present','*Present >QL','Detected Not Quantified','Present Above Quantification Limit','Systematic Contamination')))

# Cut to values within CONUS
coordinates(metal.df1) = c("longitude", "latitude") # convert to spatial points dataframe
proj4string(metal.df1) <- crs(mapUSm) # setting initial crs
inside.US = !is.na(over(metal.df1, as(mapUSm,'SpatialPolygons')))
metal.df1$insideUS = inside.US
metal.df2 = metal.df1[metal.df1$insideUS==TRUE,]
#drop geometry
metal.df2 <- st_drop_geometry(metal.df2) |>
  as.data.frame()

## 4. Detection Limit Handling -------------------------------------------------------------------------------------------------------------------

# remove records with problematic detection limits 
if (metal.code %in% c('Cd')) {  
  metal.df2 <- metal.df2 %>%
    # remove records with no reported detection limit
    # remove records where detection limit is 0 or above the MCL
    filter(!is.na(detect.limit), detect.limit > 0, detect.limit < MCL) %>%
    mutate(DL.missing = as.numeric(detect.limit==0|is.na(detect.limit)))
} else if (metal.code %in% c('U','Sr','Li','Mn','As')) {    
  # Include As here for regression
  # uniform detection limit determination for arsenic hurdle model
  metal.df2 <- metal.df2 %>%
    filter(detect.limit >= 0 | is.na(detect.limit))%>%
    mutate(DL.missing = as.numeric(detect.limit==0|is.na(detect.limit)))
}

## 5. ROS imputation + data cleaning for regression model ---------------------------------------------------------------------------------
print(paste0('Unique wells: ',length(unique(metal.df2$location.id))))
print(paste0('Total records: ',length(metal.df2$conc)))

if (metal.code %in% c('Mn','Li','As','Cd','Sr')) {
  AS <- metal.df2
  
  metal.df2 <- metal.df2 |>
    mutate(censored = conc==0) |>
    # for censored concentrations, use detection limit values
    # for not censored concentrations, use the actual values
    mutate(ros.conc = ifelse(censored, detect.limit, conc)) |>
    sample_frac(1) |>
    arrange(ros.conc)
   
  # perform ROS imputation
  ros.object = ros(metal.df2$ros.conc, metal.df2$censored)
  
  # Add the ROS imputation results back into the dataframe
  metal.df2 <- metal.df2 |>
    mutate(censored.conc = ros.object$modeled,
           ros.censored = ros.object$censored)
  
  # double check that the imputation  was done correctly 
  # For all uncensored values, censored.conc should equal the original conc
  uncensoredDF_errors = subset(metal.df2, censored==FALSE & censored.conc!=conc)
  assert_that(nrow(uncensoredDF_errors) == 0)
  
  # for all censored values, the imputed value should be different from the DL 
  censoredDF_errors = subset(metal.df2, censored==TRUE & censored.conc==ros.conc)
  assert_that(nrow(censoredDF_errors) == 0)
  
  print(paste0('imputation alignment errors: ', sum(metal.df2$ros.censored - metal.df2$censored))) # should be zero b/c they should all match 
  print(paste0('errors: ', sum(length(uncensoredDF_errors$conc), length(censoredDF_errors$conc))))
  
  # clean up the dataframe to simplify after the above steps 
  metal.df2 <- metal.df2 |>
    dplyr::select(-c(ros.censored))

} 


## 5. Select well value (merge multiple measurements at the same well location) -----------------------------------------------------------------------------
# select most recent value (and use median value when there are multiple measurements on the same day)
metal.df3 = metal.df2 |> 
  group_by(location.id) |>
  dplyr::filter(date==max(date)) |> # select data from the most recent date for each well 
  dplyr::mutate(median.conc = median(censored.conc)) |> # calculate the median value for each group (if there are > 1 value reported per date) 
  dplyr::slice(which.min(abs(censored.conc-median.conc))) |> # choose the row whose value is closest to the summary conc calculated above
  dplyr::select(-c(median.conc)) |> # drop the helper column
  as.data.frame() |>
  dplyr::select(c(location.id, date, conc, censored.conc, censored, well.depth, data.source, ros.conc, DL.missing, detect.limit))

print("WQP Processing Complete")


## 6. Save cleaned data ---------------------------------------------------------------------------------------------------------------------------
saveRDS(metal.df3, paste0("Data_Files/", metal.code,"_Combined_WQP.rds"))
print(paste("WQP Data Saved for", metal))
