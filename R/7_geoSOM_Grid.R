##############################################################################
# geoSOM_Grid.R
# Author: Jennifer Sun and Cindy Hu 
# Date: March 2024
# Purpose: Produce geographic SOM grids to identify clusters of metal mixtures 
##############################################################################

packages <- c('tidyverse','cowplot','NatParksPalettes','kohonen','Rcpp','sf')
lapply(packages, library, character.only=TRUE)

setwd(here::here("data"))

set.seed(11670)

## Functions for executing geoSOM clustering -----------------------------------------------------------------------------------

# scale data
ScaleMax <- function(x){
  c_dat <- scale(x, center = T, scale = F)
  s_dat <- c_dat / max(abs(c_dat), na.rm = T)
  return(s_dat)
} 

# geoSOM
"geosom" <- function(data, coords, k, norm, ...) {
  # subset spatial data
  geodata = data.matrix(data[,(names(data) %in% coords)])
  # subset non-spatial data
  data = data.matrix(data[,!(names(data) %in% coords)])
  spatial.weight = c(1, k)
  # if data to be normalised, we do it here
  # xyz cannot be scaled disproportionately
  if (isTRUE(norm)) {
    # data = scale(data, center = T, scale = T)
    data = lapply(data, ScaleMax) %>% data.frame()
  }
  # disregard geography if spatial weight is zero
  if (k==0) {
    supersom(list(data), normalizeDataLayers=FALSE, ...)
  } else {
    supersom(list(data, geodata), user.weights=spatial.weight, normalizeDataLayers=FALSE, ...)
  }
}


# plot function 3 x 2
"show_results" = function(kohobj) {
  par(mfrow = c(2, 3))
  codes = getCodes(kohobj)
  # in the case of >1 attribute, only the first is plotted in "property" plot
  if (is.list(codes)) {
    plot(kohobj, type = "codes", whatmap = 1, shape="straight")
    plot(kohobj, type = "property", property = codes[[1]][, 1], shape="straight")
  } else {
    plot(kohobj, type = "codes", shape="straight")
    plot(kohobj, type = "property", property = codes[, 1], shape="straight")
  }
  plot_types = c("mapping", "counts", "dist.neighbours", "changes")
  for (plot_type in plot_types) {
    plot(kohobj, type = plot_type, shape="straight")
  }
}


# plot clusters
volcanoes_pal <- colorRampPalette(colors = NatParksPalettes$Volcanoes[[1]])

plotSOM <- function(som.model, clusters, k, t.ype = 'mapping'){
  plot(som.model,
       type = t.ype,
       keepMargins = F,
       bgcol = volcanoes_pal(k)[clusters],
       # bgcol = viridis::viridis(n=k)[clusters],
       pchs = NA,
       # col = NA
  )
  add.cluster.boundaries(som.model, clusters, col = "white")
}

# replace values with percentiles
Replace_percentile <- function(vec1, vec2) {
  percentiles <- data.frame(Value = quantile(vec1, probs = seq(0, 1, by = 1e-3), 
                                             names = T, na.rm = TRUE)) %>% unique()
  percentiles$Percentile <- row.names(percentiles) 
  percentiles$Percentile <- as.numeric(gsub("%", "", percentiles$Percentile))
  percentiles$Percentile[is.na(percentiles$Percentile)] <- 0
  
  ranks  <- cut(vec2, breaks = percentiles$Value, labels = percentiles$Percentile[-1],
                include.lowest = F, right = T, na.pass = T)
  
  ranks <- ranks %>% as.character() 
  ranks[which(is.na(ranks) == T)] <- 0
  return(ranks)
}


## 1. Read in contaminant mixture data -----------------------------------------------------------------------------------------------------------
data = read.csv('Combined_GridPredicts_xgboost_std_simple_predlogconcadj.csv', row.names=NULL)
data = subset(data, select = -c(X))

## 2. Create SOM grid -----------------------------------------------------------------------------------------------------------------------------

data_train <- data
data_train[c('As', 'Mn','Sr','Li','Cd')] <-
  lapply(data %>% select(c(
    'As', 'Mn','Sr','Li','Cd'
  )), ScaleMax) %>% data.frame()
summary(data_train)

# How many samples in data?
# sample.size = nrow(som_train)

# Choose grid size
grid.size <- 16 # assuming square grid
k.val <- 15

## Create grid
som.grid  <- somgrid(xdim = grid.size, ydim = grid.size, 
                     topo = 'hexagonal', 
                     neighbourhood.fct = "gaussian", # neighbourhood.fct = "bubble",
                     toroidal = T)

# Create and save model.
geosom.model <- geosom(data   = data_train,
                       coords = c("longitude", "latitude"),
                       norm   = F,
                       k      =  k.val, # spatial weight
                       rlen   = 1000,
                       #maxNA.fraction = 0.5,
                       grid   = som.grid)

# Plot model results.
show_results(geosom.model)

# Save geosom.model
saveRDS(geosom.model,
        paste0('SOM_grids/geoSOM_grid', '_k', k.val, '.rds'))

## 3. Visualize results -----------------------------------------------------------------------------------------------------------------------------

## Plotting training progress
plot(geosom.model, type = "changes")

## Plot number of observations in ea node in map
plot(geosom.model, type = "counts",
     palette.name = viridis::viridis_pal()) # red = least # observations, yellow = most observations, gray = no observations

## Mean distance; smaller the distance, better the object is represented by codebook vectors
plot(geosom.model, type = "quality",
     palette.name = viridis::viridis_pal())

## Sum of distance to all immediate neighbors (U-matrix); low values indicates clusters, hi values indicates separation b/w clusters 
plot(geosom.model, type = "dist.neighbours",
     palette.name = viridis::viridis_pal())

plot(geosom.model, type='codes')

