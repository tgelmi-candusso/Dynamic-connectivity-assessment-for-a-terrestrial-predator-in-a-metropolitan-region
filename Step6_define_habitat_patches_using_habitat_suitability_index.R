#### Step 6: Generating map with habitat patches from habitat suitability index ##########

# we used the landscape-only model to generate habitat patches for the connectivity analysis using Omniscape
# these patches were used for all the other models too. 
# The top 1% of the habitat suitability index values were used as habitat.

# Estimate habitat suitability across all cells in the map

habmodel22 <- readRDS("habmodel22.rds") # call best fit landscape-only model output from step 3
m <- fixef(habmodel22)

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
  
  suit<- m$cond[grepl(land_var[i], names(m$cond))]*Y
  suit_stack<- stack(suit_stack, suit)
  
  print(paste(land_var[i], "- done"))
  
}

# estiamte habitat suitability index with all landscape layers

habsuit_habonly <- exp(calc(suit_stack,sum))

writeRaster(habsuit_habonly, "Habitat_suitability_map.tif")

plot(habsuit_habonly,
     main="Habitat Suitability\nEnvironmental attributes only",
     xlab = "UTM Westing Coordinate (m)", 
     ylab = "UTM Northing Coordinate (m)")

# select only top 1% as habitat patches (0.99 quantile)
HS01 <- habsuit_habonly
HS01[HS01[] < quantile(habsuit_habonly, probs = c( 0.99))] = NA
writeRaster(HS01, "habitat_patches.tif")

