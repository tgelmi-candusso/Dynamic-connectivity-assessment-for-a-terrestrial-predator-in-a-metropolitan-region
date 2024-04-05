####### Step 5: create COST SURFACE MAPS with RSS coefficients #######

#the estimate of cost surface value or resistance values, is based on the selection strength of coyote across the range of values of each landscape variable.
#Each cell in the map will have a combination of the landscape variables, and each will contribute a level of resistance to the map cell. 
# This script is for the habitat model that has no interaction factors

#cost coefficient when there is no interaction factor (day = 0, i.e. night) = 1-(b0+(by*Y))
#cost coefficient when there is an interaction factor (day = 1, i.e. day)= 1-(b0+by*Y+bx+byx*Y)
#where b0 = intercept of logRSS, by= coefficient of the urban variable, bx = coefficient of the second variable and byx = coefficient of the interaction between 1 and second variable in this case day and the urban variable.

## libraries 
require(sf)
library(maptools)  
library(ggplot2)
library(raster)
library(lubridate)
library(amt)
require(dplyr)
library(parallel)
library(glmmTMB)
library(tidyverse)
library(buildmer)
##setwd()

##base map to make sure all maps align
shape <- read_sf(dsn = ".", layer = "studyarea_fixed") #object with a polygon of the study area, mainly to clean little water infiltrations from surrounding open water in the study area that were creating noise in the connectivity. An option is also to mask the areas with water with a high resistance value within the range. 
shape1 <- st_transform(shape, "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")
plot(shape1)
####raster###

#alternative to read directly from .tif, use the following code, recalled from step 2
# NDVI <- raster("TO_NDVI_clean.tif") #tighter boundary, surrounding water excess deleted hoping to delete artifact
# names(NDVI) <- "NDVI"

land_var <- c('NDVI', 'POP', 'BUP', 'DT_LT', 'DT_MT', 'DT_HT', 'DT_NT', 'DT_PS')
i=8
#read stats_rss table from step 4 if it's not present in the environment.
stats_table_all <- readRDS ("stats_table_all.rds")
resistance_layers <- stack()

for (i in 1:length(land_var)){
  #get values needed to estimate resistance
  stats <- stats_table_all[,grepl(land_var[i], colnames(stats_table_all))]
  b0 <- stats["intercept"]
  by <- stats["coefficient"]
  
  #fix map so they all assemble nicely 
  landscape_layer <- readRDS(paste0(land_var[i],"_clean.rds"))

  #alternative to read directly from .tif, use the following code, recalled from step 2
  # landscape_layer <- raster("*.tif") #tighter boundary, surrounding water excess deleted hoping to delete artifact
  # names(landscape_layer) <- paste0(land_var[i]) #be aware of the crs of the raster vs reference map (shape1)
  
  r2_1 <- crop(landscape_layer, extent(shape1))
  if (i==1){
    r2 <- r2_1 
  }
  resample(r2_1, r2) -> r3
  Y <- mask(r3, shape1)
  
  #estimate resitance surface across the map for the landscape variable
  resistance_layer  <- exp(1-(b0+(by*Y)))
  names(resistance_layer) <- paste0(land_var[i])
  resistance_layers <- stack(resistance_layers, resistance_layer )
  
  print(paste(land_var[i], "- done"))
  }

Cost_surface_layer <- calc(resistance_layers,sum)
Cost_surface_layer <- Cost_surface_layer - minValue(Cost_surface_layer) +1 #bring all to zero and add 1 so it can run in omniscape
#plot(Cost_surface_layer)

writeRaster(Cost_surface_layer, "Landscape_only_model_Cost_surface_layer.tif", overwrite=TRUE) #differentiate for each model

