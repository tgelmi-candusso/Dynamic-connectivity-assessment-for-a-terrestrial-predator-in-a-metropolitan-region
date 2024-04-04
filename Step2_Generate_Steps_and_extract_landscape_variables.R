##### Step 2: Extract landscape information at each step #####


## add libraries 
library(amt)
library(sf)
library(raster)
library(tidyverse)
library(ggspatial)
library(fasterize)
library(tmap)
library(maptools)
library(lubridate)

## set working directory
setwd("C:/Users/tizge/Documents/Toronto project/Coyote tracking by MNRF/Resource Selection Function/ResourceSelectionFunction")


#create object for NAD83 crs because I am lazy and it is long
NAD83 <- "+proj=utm +zone=17 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#Read in coyote GPS locations 
coyote <- readRDS("coyotes_sf.rds") 

#Extract coordinates and convert to tibble for making amt track object
coyote <- coyote %>% 
  mutate(x = st_coordinates(coyote)[row_number(), 1], 
         y = st_coordinates(coyote)[row_number(), 2]) %>% 
  st_set_geometry(NULL)


#clean column names with spaces
coyote$`Temp [C]` -> coyote$temp

#Convert coyotes into an amt track object ## CREATING STEPS
coyote_amt <-  mk_track(coyote, .x = x, .y = y, .t = LOCDateTime, 
                        #crs = CRS(NAD83),
                        Id = Id, 
                        sex = Sex.x, age = Age.x, status = life_stage, 
                        day = daytime, 
                        altitude = Altitude,
                        temp = temp,
                        breeding = breeding, pup_rearing = pup_rearing, dispersal = dispersal,
                        winter=veg_bare, summer=veg_covered
)

saveRDS(coyote_amt, "coyote_amt.rds")



###need raster here#####

#vegetation
NDVI <- raster("TO_NDVI_clean.tif") #tighter boundary, surrounding water excess deleted hoping to delete artifact
names(NDVI) <- "NDVI"
NDVI <- scale(NDVI) #scaling directly here the layer, scaling so the mean is the middle point.
plot(NDVI)
saveRDS(NDVI, "NDVI_clean.rds") # saving for using when generating the resistance map

#population density
POP <- raster("TO_Pop_clean.tif") #this is number of people per cell which was originally 250m but I resampled to 30m 
names(POP) <- "POP"
POP <- scale(POP)
plot(POP)
saveRDS(POP, "POP_clean.rds")


#built infrastructure
BUP <- raster("TO_BUP_clean.tif") #this is percentage of built surface area within a 250m cell, and has been resampled to 30m
names(BUP) <- "BUP"
BUP <- scale(BUP)
plot(BUP)
saveRDS(BUP, "BUP_clean.rds")


#Distance to Linear Features
DT_LT <- raster("TO_LT_LFT.tif") 
names(DT_LT) <- "DT_LT"
DT_LT <- scale(DT_LT) #when scaling it doesnt matter whether it is in meters or kilometers.
plot(DT_LT)
saveRDS(DT_LT, "DT_LT_clean.rds")


DT_MT <- raster("TO_MT_LFT.tif") 
names(DT_MT) <- "DT_MT"
DT_MT <- scale(DT_MT)
plot(DT_MT)
saveRDS(DT_MT, "DT_MT_clean.rds")


DT_HT <- raster("TO_HT_LFT.tif") 
names(DT_HT) <- "DT_HT"
DT_HT <- scale(DT_HT)
plot(DT_HT)
saveRDS(DT_HT, "DT_HT_clean.rds")


DT_NT <- raster("TO_NT_LFT.tif") 
names(DT_NT) <- "DT_NT"
DT_NT <- scale(DT_NT)
plot(DT_NT)
saveRDS(DT_NT, "DT_NT_clean.rds")


DT_PS <- raster("TO_PS_LFT.tif") 
names(DT_PS) <- "DT_PS"
DT_PS <- scale(DT_PS)
plot(DT_PS)
saveRDS(DT_PS, "DT_PS_clean.rds")


#Start loop for selecting random steps for availability in SFF
#Unique id
ids <- unique(coyote_amt$Id)
#ids1 <- ids[c(1:27,29:34)] ###take urban020 out and do it separately, was giving trouble.
ids1 <- ids[c(1:2, 5:34)] ###take urban020 out and do it separately, was giving trouble,
## take also MISS003 nd MISS004 they were resampling at 5 instead of 3, so ettter using their median of 2.5
i=5
for (i in 1:length(ids)) {
  
  #Create sub dataframe
  subdf <- filter(coyote_amt, Id == ids[i])
  #print(paste(i,ids[i],summarize_sampling_rate(subdf)$min))
  
  #resample track so it is regular
  subdf <- track_resample(subdf, 
                          rate = minutes(as.integer(summarize_sampling_rate(subdf)$median * 60)), 
                          tolerance = minutes(15), start = 1)
  
  #set seed so it is reproducible
  set.seed(2143)
  #create 9 random available steps for each used step
  sub_track <- steps_by_burst(subdf, keep_cols = "end", lonlat = FALSE) %>% 
    random_steps(n = 9)
  
  #View(subdf)
  #extract landscape information for each step (used and available)
  sub_track <- sub_track %>% 
    extract_covariates(NDVI, where = "both")%>% 
    extract_covariates(POP, where = "both")%>% 
    extract_covariates(BUP, where = "both")%>% 
    extract_covariates(DT_LT, where = "both")%>% 
    extract_covariates(DT_MT, where = "both")%>% 
    extract_covariates(DT_HT, where = "both")%>% 
    extract_covariates(DT_NT, where = "both")%>% 
    extract_covariates(DT_PS, where = "both")
  
  
  #View(sub_track)
  #Create a function for converting points into a line
  make_line <- function(x1_, y1_, x2_, y2_) {
    st_linestring(matrix(c(x1_, y1_, 
                           x2_, y2_), 2, 2, byrow = TRUE))
  }
  
  #Convert points from track into lines for each step and attach relevant attributes
  sub_step <- sub_track %>%
    dplyr::select(x1_, y1_, x2_, y2_) %>% 
    pmap(make_line) %>% 
    st_as_sfc(crs = st_crs(NAD83)) %>% 
    {tibble(animal = sub_track$Id, 
            sex = sub_track$sex, 
            age = sub_track$age,
            status = sub_track$status,
            burst = sub_track$burst_, 
            step = sub_track$step_id_,
            used = sub_track$case_, 
            interval = sub_track$dt_, 
            step_length = sub_track$sl_/1000,
            turn_angle = sub_track$ta_,
            NDVI_start = sub_track$NDVI_start,
            NDVI_end = sub_track$NDVI_end,
            POP_start = sub_track$POP_start ,
            POP_end = sub_track$POP_end ,
            BUP_start =  sub_track$BUP_start ,
            BUP_end =  sub_track$BUP_end ,
            DT_LT_start =  sub_track$DT_LT_start ,
            DT_LT_end =  sub_track$DT_LT_end ,
            DT_MT_start =  sub_track$DT_MT_start ,
            DT_MT_end =  sub_track$DT_MT_end ,
            DT_HT_start =  sub_track$DT_HT_start ,
            DT_HT_end =  sub_track$DT_HT_end ,
            DT_PS_start =  sub_track$DT_PS_start ,
            DT_PS_end =  sub_track$DT_PS_end ,
            DT_NT_start =  sub_track$DT_NT_start ,
            DT_NT_end =  sub_track$DT_NT_end,
            day = sub_track$day,  
            breeding = sub_track$breeding, 
            pup_rearing = sub_track$pup_rearing, 
            dispersal = sub_track$dispersal, 
            winter = sub_track$winter, 
            summer = sub_track$summer, 
            geometry = .)} %>% 
    st_sf() 
  
  #if else loop to save output
  if (i == 1) {
    step <- sub_step 
  } else {
    step <- rbind(step, sub_step)
  }
  
  #print animal ID currently being worked on as process indicator
  cat("Processed animal", ids[i], "\n")
}


## check time interval between steps
step1 <- step

int <- step1 %>% 
  st_set_geometry(NULL) %>% 
  group_by(animal) %>% 
  summarise(interval = mean(interval))

#create a unique id for every step (both used and available)
step1[,"UID"] <- 1:nrow(step1) 


#Create  version of all used steps to be included in available pool and set UID to -999
temp_step <- filter(step1, used == 1)
temp_step[, "used"] <- 0
temp_step[, "UID"] <- -999
#Combine to overall data set
step_final <- rbind(step1, temp_step)

#Convert Habitat column into dummy variables

#Convert factor variables to dummy binary variables
coyote_steps_9 <- step_final %>% 
  mutate(adult = ifelse(age == "A", 1, 0), 
         male = ifelse(sex == "M", 1, 0),
         day = ifelse(day == "day", 1, 0),
         transient = ifelse(status == "Transient", 1,0))


####filter coyotes outside of region of interest

coyote_steps_9 <- coyote_steps_9 %>% 
  filter(!is.na(NDVI_end)) %>% 
  filter(!is.na(BUP_end)) %>%
  filter(!is.na(POP_end)) %>%  
  filter(!is.na(DT_LT_end))%>%  
  filter(!is.na(DT_MT_end))%>%  
  filter(!is.na(DT_HT_end))%>%  
  filter(!is.na(DT_NT_end))%>%  
  filter(!is.na(DT_PS_end))


####filter animals outside of region

##animals under30days, or completely outside of study site included in the raster files 
ids2 <- unique(coyote_steps_9$animal)

coyote_steps_9 <- coyote_steps_9 %>% filter(!(animal %in% c(
  "MISS021",
  "URBAN006",
  "URBAN008",
  "URBAN010",
  "URBAN015",
  "URBAN021",
  "URBAN024")))

coyote_steps_9 -> coyote_steps #This file is in the repository.

saveRDS(coyote_steps, "coyote_steps_clean_2.rds") 

#save file for repository
coyotes_steps <- read_rds("coyote_steps_clean_2.rds")




