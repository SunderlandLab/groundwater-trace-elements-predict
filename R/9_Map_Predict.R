###############################################################################
#File name: Map_Predict.R
#Author: Jennifer Sun
#Date: July 2024
#Purpose: Predict concentrations on a continuous 0.1x0.1 degree grid across the US
##############################################################################


#### Attach Libraries and Set Working Directory ####
packages <- c('plyr','dplyr','pROC','rpart','stringr','tidyverse','tidymodels','xgboost','rgdal','sp')
lapply(packages, library, character.only=TRUE)

setwd(here::here("data"))

metal = 'Arsenic'
metal.code = "As"
version.number = '7a'
grid.version = 'tune8_std'
predictor.version = 'v1'
model.run = 'xgboost_tune8_std_simple' 
cluster.folder = 'tune8/'
folder.name = paste0(metal.code,'_',version.number,'/')
modelfolder.name = paste0(metal.code,'_',version.number,'/',grid.version,'/')

filename = paste0("R_Output/",modelfolder.name,metal.code,"_ModelPackage_",version.number,"_",predictor.version,"_",model.run,'.RData')
load(filename)


## 1. Create 0.1 x 0.1 degree grid with predictors across the US ----------------------------------------------------------------------------------------------
# Load in Map of US
mapUSm<- readOGR(dsn="CoVar/US48", layer="US_48states") # projected map of the U.S.(in meters)

# Create Grid over CONUS
grid = makegrid(mapUSm, cellsize = 0.1) # cell size is in map units, which here is lat/long degrees
grid = SpatialPoints(grid, proj4string = crs(mapUSm)) # 775 x 1254 points, grid dimensions

# Cookie cut out extra points
inside.US<-!is.na(over(grid,as(mapUSm,'SpatialPolygons')))
grid$insideUS<-inside.US
grid<-grid[grid$insideUS==1,]

# Extract predictor variables (use 2_Predictors_Extract)
grid <- readRDS('R_Output/MapGrid/mapgrid_Model_Ready.rds') # saved grid with predictors

## 2. Prepare separate grid for each metal to be compatible with xgboost model inputs ----------------------------------------------------------------------------
# Rename TRI and SEMS columns you want to use 
grid@data = grid@data %>% dplyr::rename('TRI.total.impact' = paste0('TRI.total.impact.',metal.code),
                                        'TRI.water.impact' = paste0('TRI.water.impact.',metal.code),
                                        'SEMS' = paste0('SEMS.',metal.code))

# Add back in non-predictor columns so the model workflow can run; these columns are discarded before prediction
newcols = setdiff(colnames(df), colnames(grid@data))
for (col in newcols) {
  grid@data = grid@data %>% dplyr::mutate(!!col := 0)
}

# Remove unnecessary new columns
grid@data = grid@data %>% dplyr::select(-gridIndex)

# Manually impute factor values that don't exist in training data
if (metal.code %in% c('As','Mn','Cd','Li','Sr')) {
  grid@data['aquifer'][grid@data['aquifer'] == 402] <- 999 # no information
  grid@data['aquifer'][grid@data['aquifer'] == 406] <- 999 # no information
  grid@data['aquifer'][grid@data['aquifer'] == 601] <- 999 # no information
  grid@data['str_int'][grid@data['str_int'] == 5] <- 4 # Next closest category
} 

if (metal.code %in% c('As','Mn','Cd','Li')) {
  grid@data['landcover_500m'][grid@data['landcover_500m'] == 12] <- 11 # there are not many of these values; perennial ice/snow counted as open water
  grid@data['landcover_500m'][grid@data['landcover_500m'] == 0] <- 11 # there are not many of these values; no data counted as open water 
}

df.grid = grid@data
df.grid$longitude = grid@coords[,1]
df.grid$latitude = grid@coords[,2]

saveRDS(df.grid, paste0(metal.code, '_MapGrid.rds'))

## 3. Predict on grid ------------------------------------------------------------------------------------------------------------------------------------------------------
# predict
grid.predictions = data.frame(predict(final.full.model, new_data = grid@data)) # for regression # changed to final full model 9/9/24

# add predictions back to grid
grid$pred_logconc = grid.predictions$.pred
grid$pred_conc = 10^grid$pred_logconc


## 4. EDM adjustment --------------------------------------------------------------------------------------------------------------------------------------------------------
# prepare training predictions 
df.train.uncens = subset(df.train, is.imputed==0)
train.predictions.uncens = predict(final.full.model, df.train.uncens) %>% bind_cols(df.train.uncens) # changed to final.full.model 9/9/24

train.pred.adj = train.predictions.uncens
grid.predictions.adj = grid.predictions

train.pred.adj$logconc_ordered = sort(train.pred.adj$logconc)
train.pred.adj$.pred_ordered = sort(train.pred.adj$.pred)

# transformation: adjust test predictions based on train predictions vs. obs skew
tbl_ordered = train.pred.adj[,c('logconc_ordered','.pred_ordered')]
grid_adj = data.frame(approx(tbl_ordered$.pred_ordered, tbl_ordered$logconc_ordered, xout=grid.predictions.adj$.pred))

grid.predictions.adj$.pred_adj = grid_adj$y # append transformed data to the original dataset
grid.predictions.adj$.pred_check = grid_adj$x # this should match .pred


# add adjusted predictions back to the grid
grid$pred_logconc_adj = grid.predictions.adj$.pred_adj
grid$pred_conc_adj = 10^grid$pred_logconc_adj

## 5. Save as CSV for plotting in ArcGIS -------------------------------------------------------------------------------------------------------------------------------------

# extract lat/long
grid$longitude = grid@coords[,1]
grid$latitude = grid@coords[,2]

# R_Output
write.csv(grid@data, paste0('R_Output/',modelfolder.name, metal.code, '_GridPredicts_',version.number,"_",predictor.version,"_",model.run,'_finalfullmodel.csv'))

