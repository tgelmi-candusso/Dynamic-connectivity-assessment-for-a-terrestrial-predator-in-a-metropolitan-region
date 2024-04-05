#### Step 6b: Habitat suitability index from models with interactions components ##########

# we used the landscape-only model to generate habitat patches for the connectivity analysis using Omniscape
# these patches were used for all the other models too. 
# However, this is the code to generate habitat suitability maps from models with interaction components


# Estimate habitat suitability across all cells in the map

daymodel37 <- readRDS("daymodel37.rds") # call best fit landscape-only model output from step 3
m <- fixef(daymodel37)

land_var <- c('NDVI', 'POP', 'BUP', 'DT_LT', 'DT_MT', 'DT_HT', 'DT_NT', 'DT_PS')
#study area outline
shape <- read_sf(dsn = ".", layer = "studyarea_fixed") #object with a polygon of the study area, mainly to clean little water infiltrations from surrounding open water in the study area that were creating noise in the connectivity. An option is also to mask the areas with water with a high resistance value within the range. 
shape1 <- st_transform(shape, "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")

suit_stack <- stack()
for (i in 1:length(land_var)){
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
  
  suit_level0<- m$cond[grepl(land_var[i], names(m$cond))][1]*Y
  suit_stack_level0<- stack(suit_stack_level0, suit_level0)
  
  suit_level1<- (m$cond[grepl(land_var[i], names(m$cond))][1]+
                   m$cond[grepl(land_var[i], names(m$cond))][2])*Y
  suit_stack_level1<- stack(suit_stack_level1, suit_level1)
  
  # for THREE LEVEL INTERACTION FACTORS add:
  # suit_level2<- (m$cond[grepl(land_var[i], names(m$cond))][1]+
  #                  m$cond[grepl(land_var[i], names(m$cond))][3])*Y
  # suit_stack_level2<- stack(suit_stack_level2, suit_level2)
  
  print(paste(land_var[i], "- done"))
  
}

# estiamte habitat suitability index with all landscape layers

habsuit_model_level0 <- exp(calc(suit_stack_level0,sum))
writeRaster(habsuit_model_level0, "Habitat_suitability_map_level0.tif") #change map name to include the specific model it is based on and the level description. e.g. here daymodel_night

habsuit_model_level1 <- exp(calc(suit_stack_level1,sum))
writeRaster(habsuit_model_level1, "Habitat_suitability_map_level1.tif") # e.g. add to map name daymodel_day

# for THREE level interaction factors add:
# habsuit_model_level2 <- exp(calc(suit_stack_level2,sum))
# writeRaster(habsuit_model_level2, "Habitat_suitability_map_level2.tif") # e.g. add to map name daymodel_day


plot(habsuit_habonly,
     main="Habitat Suitability\nEnvironmental attributes only",
     xlab = "UTM Westing Coordinate (m)", 
     ylab = "UTM Northing Coordinate (m)")

#if extracting potential habitat patches from the habitat suitability map here is the code, if you want to base your habitat patches on context specific situations. However, We were more interest on the effect of movement across the urban landscape rather than the movement between specific patches so for comparison sake we kept patches constant across context
# select only top 1% as habitat patches (0.99 quantile) change quantile number if needed.
HS01 <- habsuit_model_level0
HS01[HS01[] < quantile(HS01, probs = c( 0.99))] = NA
writeRaster(HS01, "habitat_patches_level0.tif")

HS01 <- habsuit_model_level1
HS01[HS01[] < quantile(HS01, probs = c( 0.99))] = NA
writeRaster(HS01, "habitat_patches_level1.tif")

### IMPORTANT NOTE ####

#### NOTE ON HOW TO DEAL WITH LANDSCAPE-VARIABLE SPECIFIC NON-SIGNIFICANT INTERACTIONS WITHIN A BEST-FIT MODEL

## In our analysis, WE only added the interaction component for the landscape variables that had a significant interaction with the demographic or temporal factor. 
## This is not included in the code above, because it is a manual selection following the model results and really it is at the scientist's discretion on whether to consider this or not in their analysis
## We decided to do this because our math relies on the model's coefficients, but this coefficient has a confidence interval and with a non-significan p-value, we observed this confident interval was too large, so modifying adding a context-specific weight based on that interaction coefficient value alone wouldn't have been representative of the system. Either way we tried both ways and the final results doesnt change much, see how it looks in your system. 

## add this code to the loop above if exluding non-significant interaction weights, replacing code where we generate resistance_layers in lines 32-43.

# define non-significant interaction coefficients
non_significant_interactions_withLevel1 <- c("POP","DT_MT","DT_LT") #add here all landscape variables that didnt have a significant interaction with a spatial or demographic variable. This object will change with each model of course.

suit_level0<- m$cond[grepl(land_var[i], names(m$cond))][1]*Y
suit_stack_level0<- stack(suit_stack_level0, suit_level0)

if(land_var[i] %in% non_significant_interactions_withLevel1){ #day = 1 (i.e. day) 
  suit_level1<- (m$cond[grepl(land_var[i], names(m$cond))][1])*Y
}else{
  suit_level1<- (m$cond[grepl(land_var[i], names(m$cond))][1]+
                   m$cond[grepl(land_var[i], names(m$cond))][2])*Y
}
suit_stack_level1<- stack(suit_stack_level1, suit_level1)

# for THREE LEVEL INTERACTION MODELS also include:
# non_significant_interactions_withLevel2 <- c("POP", "DT_LT") #add here all landscape variables that didnt have a significant interaction with a spatial or demographic variable. This object will change with each model of course.
# if(land_var[i] %in% non_significant_interactions_withLevel2){ #e.g. in breeding model dispersal = 1 and breeding = 0, i.e. dispersal season
#   suit_level2<- (m$cond[grepl(land_var[i], names(m$cond))][1])*Y
# } else{
#   suit_level2<- (m$cond[grepl(land_var[i], names(m$cond))][1]+
#                  m$cond[grepl(land_var[i], names(m$cond))][3])*Y
# }
# suit_stack_level2<- stack(suit_stack_level2, suit_level2)


