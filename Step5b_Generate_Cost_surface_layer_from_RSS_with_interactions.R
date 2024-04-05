####### Step 5b: create COST SURFACE MAPS with RSS coefficients from models with interaction factors#######

#this is for models WITH interaction factors, note changes for THREE LEVEL factors included as comments throughout the scripts
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
  bx <- stats["coefficient2"]
  byx <- stats["interaction"]

  ##models with THREE LEVEL interaction factors:
  # b0 <- stats["intercept"]
  # by <- stats["coefficient"] #pup-rearing
  # bx <- stats["coefficient2.1"]
  # bx1 <- stats["coefficient2.1"]
  # byx <- stats["interaction1"] #breeding
  # byx1 <- stats["interaction2"] #pup-rearing (keep track of these >2-level-model levels)
  

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
  #-- TWO level interaction models
  resistance_layer_level0 <- exp(1-(b0+(by*Y))) #day = 0 (i.e. night)
  resistance_layer_level1 <-  exp(1-(b0+by*Y+bx+byx*Y)) #day = 1 (i.e. day)
  
## for THREE level interaction models do this instead:
# resistance_layer_level0 <- exp(1-(b0+(by*Y))) #pup-rearing
# resistance_layer_level1 <- exp(1-(b0+by*Y+bx+byx*Y)) #breeding
# resistance_layer_level2 <- exp(1-(b0+by*Y+bx1+byx1*Y)) #dispersal
  
  names(resistance_layer_level0) <- paste0(land_var[i])
  resistance_layers_level0 <- stack(resistance_layers_level0, resistance_layer_level0 )

  names(resistance_layer_level1) <- paste0(land_var[i])
  resistance_layers_level1 <- stack(resistance_layers_level1, resistance_layer_level1 )

# run this two for THREE level models
# names(resistance_layer_level2) <- paste0(land_var[i])
# resistance_layers_level2 <- stack(resistance_layers_level2, resistance_layer_level2 )
  
  print(paste(land_var[i], "- done"))
}

### now sum up the values across landscape variables -- we did an equal weight arithmetic sum, since the contribution of each landscape variable to the overall selection is already included in the math, we dont need to weight them further.
Cost_surface_layer_0 <- calc(resistance_layers_level0,sum)
Cost_surface_layer_0 <- Cost_surface_layer_0 - minValue(Cost_surface_layer_0) +1 #bring all to zero and add 1 so it can run in omniscape
writeRaster(Cost_surface_layer_0, "Resistance_layer_level0.tif", overwrite=TRUE) #differentiate for the interaction factor level in each model # eg. this would be for daymodel day = 0 i.e. night, resistance at night
                                                                                #differentiate for the interaction factor level in each model eg. this would be for breedingmodel this would be breeding = 0 and dispersal = 0, i.e. pup-rearing, resistance during pup-rearing season
Cost_surface_layer_1 <- calc(resistance_layers_level1,sum)
Cost_surface_layer_1 <- Cost_surface_layer_1 - minValue(Cost_surface_layer_1) +1 #bring all to zero and add 1 so it can run in omniscape
writeRaster(Cost_surface_layer_1, "Resistance_layer_level0.tif", overwrite=TRUE) #differentiate for the interaction factor level in each model # eg. this would be for daymodel day = 1 i.e. day, resistance during daytime
                                                                         #differentiate for the interaction factor level in each model eg. this would be for breedingmodel this would be breeding = 1 and dispersal = 0, i.e. breeding, resistance during breeding season
#  for THREE level interaction models add this:
# Cost_surface_layer_2 <- calc(resistance_layers_level2,sum)
# Cost_surface_layer_2 <- Cost_surface_layer_2 - minValue(Cost_surface_layer_2) +1 #bring all to zero and add 1 so it can run in omniscape
# writeRaster(Cost_surface_layer_2, "Resistance_layer_level0.tif", overwrite=TRUE) #differentiate for the interaction factor level in each model eg. this would be for breedingmodel this would be breeding = 0 and dispersal = 1, i.e. dispersal, resistance during dispersal season

### IMPORTANT NOTE ####

#### NOTE ON HOW TO DEAL WITH LANDSCAPE-VARIABLE SPECIFIC NON-SIGNIFICANT INTERACTIONS WITHIN A BEST-FIT MODEL

## In our analysis, WE only added the interaction component for the landscape variables that had a significant interaction with the demographic or temporal factor. 
## This is not included in the code above, because it is a manual selection following the model results and really it is at the scientist's discretion on whether to consider this or not in their analysis
## We decided to do this because our math relies on the model's coefficients, but this coefficient has a confidence interval and with a non-significan p-value, we observed this confident interval was too large, so modifying adding a context-specific weight based on that interaction coefficient value alone wouldn't have been representative of the system. Either way we tried both ways and the final results doesnt change much, see how it looks in your system. 

## add this code to the loop above if exluding non-significant interaction weights, replacing code where we generate resistance_layers in lines 73-81. 

  # define non-significant interaction coefficients
  non_significant_interactions_withLevel1 <- c("POP","DT_MT","DT_LT") #add here all landscape variables that didnt have a significant interaction with a spatial or demographic variable. This object will change with each model of course.

  # estimate resistance for each landscape variable based on the above.  
  resistance_layer_level0 <- exp(1-(b0+(by*Y))) #day = 0 (i.e. night)
  if(land_var[i] %in% non_significant_interactions_withLevel1){ #day = 1 (i.e. day) 
    resistance_layer_level1 <- exp(1-(b0+by*Y))#+bx+byx*Y)) #interaction component excluded
  }else{
    resistance_layer_level1 <- exp(1-(b0+by*Y+bx+byx*Y))
  } 

# for THREE LEVEL INTERACTION MODELS also include:
# non_significant_interactions_withLevel2 <- c("POP", "DT_LT") #add here all landscape variables that didnt have a significant interaction with a spatial or demographic variable. This object will change with each model of course.

# if(land_var[i] %in% non_significant_interactions_withLevel2){ #dispersal = 1 and breeding = 0, i.e. dispersal season
#  resistance_layer_level2 <- exp(1-(b0+by*Y))#+bx1+byx1*Y)) #interaction component excluded
# } else{
#  resistance_layer_level2 <- exp(1-(b0+by*Y+bx1+byx1*Y))
# }

