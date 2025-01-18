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
