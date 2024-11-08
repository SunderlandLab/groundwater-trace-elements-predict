############################################################################################
#File name: Predictor_Selection.R
#Author: Jennifer Sun
#Date: June 2022
#Purpose: Final dataframe cleaning steps, including manual and automatic predictor selection
############################################################################################

#### Attach Libraries and Set Working Directory ####
packages <- c('akima','chron','corrplot','data.table','dismo','dplyr','geosphere','ggplot2','gstat',
              'lattice','maptools','NADA','ncdf4','parallel','plyr','pROC','ranger','raster','RColorBrewer',
              'rgdal','rgeos','stringr','tiff','VIM','pdp','caret','viridis','formattable','tidyverse','tidymodels',
              'car','fastDummies','doParallel','xlsx')
lapply(packages, library, character.only=TRUE)

setwd('/Users/jennifer/Documents/Harvard/Drinking Water/Metals_Modeling.nosync/DWheavymetal')

#### Define Global Variables ####
metal = 'Manganese'
metal.code = "Mn"
version.number = '7a'
MCL = 300
model.run = ''
set.seed(123)
folder.name = paste0(metal.code,'_',version.number,'/')

# Read the refined data set
data = readRDS(paste0("R_Output/",folder.name,metal.code,"_Model_Ready",version.number,".rds"))
df = data@data

## 1. Remove unwanted predictor variables -----------------------------------------------------------------------------------------------------------------
df = subset(df, select=-c(C_Gypsum, C_Aragon, WTDEP_MIN, WTDEPAMJ, Hydrate_NA)) # variables with >> 10% missingness
df = subset(df, select=-c(survey_type_int)) # this designates statsgo or ssurgo in the soil property rasters; not a soil property


## 2. Remove variables that are replicated within raw datasets ------------------------------------------------------------------------------------- 
df = subset(df, select=-c(PERML,PERMH,AWCL,AWCH,BDL,BDH,OML,OMH,SLOPEL,SLOPEH, 
                          WTDEPL,WTDEPH,ROCKDEPL,ROCKDEPH, PFLATLOW, PFLATTOT, PFLATUP)) # duplicates from statsgo soil properties; use calculated mean instead
df = subset(df, select=-c(cec_025, cec_05, cec_050, clay_025, clay_05, clay_2550, clay_3060, 
                          ec_025, ec_05, ksat_05, ph_025, ph_05, ph_2550, ph_3060, sand_025, sand_05,
                          sand_2550, sand_3060, silt_025, silt_05, silt_2550, silt_3060, paws_025, paws_050,
                          min_ksat, max_ksat)) # duplicates from #16. soil property rasters 


## 3. Remove highly correlated variables -----------------------------------------------------------------------------------------------------------------------

# remove irrelevant variables before running correlation analysis
df_corData = df %>% dplyr::select(where(is.numeric), -c(long, lat, conc, censored.conc, well.depth, ros.conc, DL.missing)) # categorical variables
df_corData = df_corData %>% dplyr::select(-c(TRI.total.impact, TRI.water.impact, SEMS)) # anthropogenic input variables that may be correlated but we will keep
df_corData = subset(df_corData, !is.na(VALUE_0_pct)) # remove any rows where there are missing landcover values still (i.e. cadmium dataset)

# identify highly correlated variables
df_corMat = cor(df_corData, method='pearson') # manually inspect 
# automatic removal of correlated variables - did not use (instead, manually select which of the correlated variables to remove)
# corVars <- caret::findCorrelation(df_corMat, cutoff=0.9) 
# corVars <- colnames(df_corData)[corVars]

# manually selected parameters to remove based on correlation coefficient (view correlation matrix)
if (metal.code == 'Sr') {
  df = subset(df, select = -c(C_Ni,C_Tot_Flds,AVG_POR,Hydrate_Y,AVG_NO10,AVG_SAND,AVG_CLAY,AVG_SILT,sand,silt, no3_pub, AVG_KV)) 
} else if (metal.code == 'Li') {
  df = subset(df, select = -c(C_Ni,C_Tot_Flds,AVG_POR,Hydrate_Y,AVG_NO10,AVG_SAND,AVG_CLAY,AVG_SILT,sand,silt, AVG_NO200, no3_pub)) 
} else if (metal.code == 'As') {
  df = subset(df, select = -c(C_Ni,C_Tot_Flds,AVG_POR,Hydrate_Y,AVG_NO10,AVG_SAND,AVG_CLAY,AVG_SILT,sand,silt,aq_rocktype)) 
} else if (metal.code %in% list('Cd', 'Mn')) {
  df = subset(df, select = -c(C_Ni,C_Tot_Flds,AVG_POR,Hydrate_Y,AVG_NO10,AVG_SAND,AVG_CLAY,AVG_SILT,sand,silt)) 
}

# Save final dataframe
saveRDS(df,paste0("R_Output/",folder.name,metal.code,"_df_PredictorsSelected.rds")) 

