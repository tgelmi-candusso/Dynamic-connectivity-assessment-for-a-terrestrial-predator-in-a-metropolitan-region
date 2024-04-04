## IDENTIFYING SOCIAL STATUS ##

## libraries ##

library(sf)
require(tidyverse)
library(lubridate)
library(gridExtra)
library(ggspatial)
library(tmap)
library(amt)
require(dplyr)
library(maptools)
library(ggmap)

#create directory for location plots
dir.create("location_plots_june_trial_3")


#create location plot for each collar_animal to look for outlier locations
for (i in 1:length(ids)) {
  subdf <- coyotes1 %>% 
    filter(Id == ids[i])
  subplot <- ggplot(subdf, aes(Longitude, Latitude)) +
    geom_point() +
    coord_equal() +
    theme_bw() +
    labs(title = ids[i])
  cat("Processing animal and collar", as.character(ids[i]), "\n")
  ggsave(paste("location_plots_june_trial_3\\", as.character(ids[i]), 
               ".png", sep = ""),
         subplot,
         width = 12,
         height = 6)
}

#warnings()
####Identify Transients####

#define season variables for use in net displacement plot axis
season.start <- seq.POSIXt(from = as.POSIXct("2012-04-01 00:00:00", 
                                             format = "%F %T"),
                           to = as.POSIXct("2021-05-15 00:00:00",  
                                           format = "%F %T"), 
                           by = "2 months")


#calculate net displacement in km for each canid ID
coyotes_sf <- coyotes_sf %>% 
  group_by(Id) %>% 
  mutate(net.disp = (st_distance(geometry, geometry[1]))) %>% 
  ungroup()
#View(coyotes_sf)

#Create track for each canid
coyote_track <- coyotes_sf %>% 
  group_by(Id) %>% 
  summarize(do_union = FALSE) %>% 
  sf::st_cast("LINESTRING") %>% 
  ungroup()
#class(coyotes_sf)
#create folder for plots
dir.create("net_displacement_june_trial_3")

#Start plot loop #########BUG: Units, will try loading ggforce

#install.packages("ggforce")
library(ggforce) ##fixed but now I get error unexpected symbol on net displacement.

##add background map



for (i in 1:length(ids)) {
  
  subdf_loc <- coyotes_sf %>% 
    filter(Id == ids[i])
  
  subdf_track <- coyote_track %>% 
    filter(Id == ids[i])
  

  map <- ggplot(subdf_loc, aes(color = as.numeric(difftime(LOCDateTime,
                                                           LOCDateTime[1])))) +
    geom_sf(data = subdf_track, inherit.aes = FALSE) +
    # ggmap()+
    annotation_scale(location = "br") +
    geom_sf() +
    scale_color_gradient2(low = "#91cf60", mid = "#ffffbf", 
                          high = "#de2d26", name = "Time",
                          labels = NULL) +
    theme_bw()
  
  
  
  #Build net displacement plot
  plot <- ggplot(subdf_loc, aes(LOCDateTime, net.disp, 
                                color = as.numeric(difftime(LOCDateTime, 
                                                            LOCDateTime[1])))) +
    geom_line(size = 1) +
    scale_color_gradient2(low = "#91cf60", mid = "#ffffbf", 
                          high = "#de2d26", name = "Time",
                          labels = NULL) +
    scale_x_datetime(breaks = season.start, 
                     date_minor_breaks = "1 month", 
                     date_labels = "%b-%Y") +
    labs(x = "Date", y = "Net_Displacement", title = as.character(ids[i])) +
    theme_linedraw() +
    theme(panel.grid.minor.x = element_line(color = alpha("black", 0.7)), 
          panel.grid.major.y =  element_line(color = alpha("black", 0)),
          panel.grid.minor.y = element_line(color = alpha("black", 0)))
  
  subplot <- grid.arrange(map, plot, nrow = 2)
  
  cat("Processing animal", as.character(ids[i]), "\n")
  ggsave(paste("net_displacement_june_trial_3//", as.character(ids[i]), 
               ".png", sep = ""),
         subplot,
         width = 12,
         height = 10)
  
  rm(subdf_loc, subdf_track, subplot, map, plot)
}


#set interactive view of maps
tmap::tmap_mode("view")

#Look through every animal to determine if/when there is a status change
test <- filter(coyotes_sf, Id == "URBAN011")

status_test <- test

#filter by datetime. start with break identified in net_displacement plots for looks more in-depth of specific period where the status might've changed
#status_test <- filter(test,GMDateTime >= as.POSIXct("2021-01-09 09:00:00", tz = "UTC"))

#Create track for each canid
status_track <- status_test %>% 
  group_by(Id) %>% 
  summarize(do_union = FALSE) %>% 
  st_cast("LINESTRING")



#Check each animal with a status change. 
#Adjust "status test" until animal status consistent 
tm_shape(status_track) +
  tm_lines(alpha = 0.3) +
  tm_shape(status_test) +
  tm_basemap(server = c("OpenStreetMap.Mapnik")) +
  tm_symbols(col = "GMDateTime", size = 0.01, palette = "YlOrRd", 
             legend.col.show = FALSE)



#Make sure locations are not excluded by start and end dates
for (i in 1:length(ids)) {
  
  status_key <- filter(coyote_status, Id == ids[i])
  sub_df <- filter(coyotes, Id == ids[i])
  sub_status <- filter(sub_df, GMDateTime < min(status_key$startGMT) |
                         GMDateTime > max(status_key$endGMT)) %>% 
    mutate(id = ids[i])
  
  if ( i == 1) {
    status_test <- sub_status
  } else {
    status_test <- rbind(status_test, sub_status)
  }
  
}

#Visually check data summaries and clean workspace
View(select(status_test, Id, GMDateTime))
unique(status_test$Id)
rm(status_test, status_track, status_key, sub_df, sub_status, subdf)

#Create index vector for status id
status_id <- unique(coyote_status$status_id)


coyote_status <- read_csv("2-GTA_coyote_status.csv")


#Identify status by start and end dates
for (i in 1:length(status_id)) {
  
  status_key <- filter(coyote_status, status_id == status_id[i])
  sub_df <- filter(coyotes_sf, Id == status_key$Id)
  sub_status <- filter(sub_df, GMDateTime >= status_key$startGMT &
                         GMDateTime <= status_key$endGMT) %>% 
    mutate(status_id = status_id[i], 
           status = status_key$status)
  
  if ( i == 1) {
    status_change <- sub_status
  } else {
    status_change <- rbind(status_change, sub_status)
  }
  
}

#Final check to make sure no locations are missing
setdiff(paste(coyotes_sf$Id, coyotes_sf$GMDateTime), 
        paste(status_change$Id, status_change$GMDateTime))

#Make sure no locations were included more than once
View(status_change[duplicated(paste(status_change$Id, status_change$GMDateTime)),])

status_change <- distinct(status_change, GMDateTime, Id, .keep_all = TRUE)

rm(status_key, sub_status, sub_df)


####FINAL CHECK###

#Create animal Summmary and remove animals with too few days
animals <- status_change %>% 
  group_by(status_id) %>% 
  summarize(loc = n(), 
            stat = unique(status),
            start = first(GMDateTime), 
            end = last(GMDateTime), 
            int = mean(difftime(GMDateTime, lag(GMDateTime), units = "hours"), 
                       na.rm = TRUE), 
            days = length(unique(as.Date(GMDateTime)))) %>% 
  filter(stat %in% c("Resident", "Transient"), 
         days > 30, 
         status_id  != "URBAN008a")
View(animals)

#Remove animals with too few days from master file
status_change <- status_change %>% 
  filter(status_id %in% animals$status_id)

####Save Files####

#Save output
write_csv(coyote_status, "2-GTA_coyote_status.csv")
saveRDS(status_change, "2-GTA_coyote_locations_cleaned.rds")
readRDS("2-GTA_coyote_locations_cleaned.rds") -> status_change
