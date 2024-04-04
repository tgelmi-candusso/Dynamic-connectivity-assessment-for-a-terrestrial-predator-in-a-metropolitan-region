#### Step 1: Clean GPS collar data ###


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

####Clean Data####

#Read in data
#remove loacations that are far away
#make datetime a datetime object

coyotes <- read_csv("Locations_Updated_01May21_n_fixed.csv", na=c("NA","N/A", "<NA>"), 
                    col_type = cols(
                      Id = col_character(),
                      "Device ID" = col_double(),
                      DOP = col_double(),
                      "Temp [C]" = col_double(),
                      x = col_double(),
                      y = col_double(),
                      GMDateTime = col_datetime(format = "%Y-%m-%d %H:%M:%S"),
                      LOCDateTime = col_datetime(format = "%Y-%m-%d %H:%M:%S"),
                      "Fix Status" = col_character(),
                      "Main [V]" = col_character(),
                      "Back [V]" = col_character(),
                      LOCDateTime_c = col_character(),
                      GMDateTime_c = col_character(),
                      GMTime = col_time(format = "%H:%M:%S"),
                      LOCTime = col_time(format = "%H:%M:%S"),
                      GMDate = col_date(format = "%Y-%m-%d"),
                      LOCDate = col_date(format = "%Y-%m-%d")))



#Create a vector of unique ids
ids <- unique(coyotes$Id)


#Summary of how many days in each collar
coyotes_summary <- coyotes %>% 
  group_by(Id) %>% 
  filter(!is.na(LOCDateTime))%>%
  dplyr::summarize(loc = n(), 
                   Id = unique(Id),
                   start = first(LOCDateTime), 
                   end = last(LOCDateTime),
                   days = length(unique(as.Date(LOCDateTime))))
View(coyotes_summary)


#Check Data
head(coyotes)
str(coyotes)

##clean location errors beyond obvious movement range seen in the map.
coyotes1 <- coyotes %>% 
  filter(!is.na(Latitude))%>% 
  filter(!is.na(Longitude))%>%
  filter(Latitude>40)%>% 
  filter(!(Id == "MISS004" & (Latitude >43.68 | Longitude>-79.4))) %>% 
  filter(!(Id == "MISS011" & Latitude <43.52)) %>% 
  filter(!(Id == "URBAN001" & Latitude >43.8)) %>% 
  filter(!(Id == "URBAN005" & Latitude >43.6)) %>% 
  filter(!(Id == "URBAN009" & Longitude >-79.5)) %>% 
  filter(!(Id == "URBAN012" & Longitude > -79.25)) %>% 
  filter(!(Id == "URBAN014" & Longitude >-79.5))%>%   
  filter(!(Id == "URBAN020" & Longitude >-79.25)) %>% 
  filter(!(Id == "URBAN023" & Longitude>-79.23))

##clean capture-recaptures from information provided by the MNRF of when these were captured and recaptured e.g. for mange:
coyotes1 <- coyotes1 %>% 
  filter(!(Id == "URBAN016" & between(LOCDate, as.Date("2020-07-16"), as.Date("2020-08-01"))))%>%
  filter(!(Id == "URBAN016" & between(LOCDate, as.Date("2020-10-22"), as.Date("2020-12-15"))))%>%
  filter(!(Id == "URBAN014" & between(LOCDate, as.Date("2019-12-14"), as.Date("2020-01-22"))))%>%
  filter(!(Id == "URBAN014" & between(LOCDate, as.Date("2020-12-18"), as.Date("2021-02-22"))))%>%
  filter(!(Id == "URBAN011" & between(LOCDate, as.Date("2019-10-01"), as.Date("2019-11-30"))))%>%
  filter(!(Id == "URBAN011" & between(LOCDate, as.Date("2020-10-21"), as.Date("2020-12-22"))))


#### Add demographics and seasonality to the steps
## 
# 
##DEFINE social status  (aka life stage). These were observed and identified manually using a spatiotemporal graph of their steps. See social_status_analysis.R

coyotes1$life_stage <- NA
#set as date time to be able to filter the life stage ranges
coyotes1$GMDateTime_c  <- as_datetime(coyotes1$GMDateTime) #transient resident table are in GMT
#check  date time is not messed up
coyotes1 %>% filter(!(GMDateTime == GMDateTime_c))

#manually defined from side analysis on their social status.
coyotes1 <- coyotes1 %>% mutate(life_stage = case_when(
  (Id=="MISS001" & GMDateTime_c >=as_datetime("2012-05-05T15:02:00Z")  & GMDateTime_c <=as_datetime("2012-07-18T07:34:00Z")) ~ "Transient",
  (Id=="MISS001" & GMDateTime_c >=as_datetime("2012-07-18T07:35:00Z")  & GMDateTime_c <=as_datetime("2013-04-24T15:00:00Z")) ~ "Resident",
  (Id=="MISS002" & GMDateTime_c >=as_datetime("2012-07-07T06:03:00Z")  & GMDateTime_c <=as_datetime("2013-02-10T17:30:00Z")) ~ "Resident",
  (Id=="MISS003" & GMDateTime_c >=as_datetime("2012-11-14T15:01:00Z")  & GMDateTime_c <=as_datetime("2013-02-25T00:00:00Z")) ~ "Resident",
  (Id=="MISS003" & GMDateTime_c >=as_datetime("2013-02-25T00:00:00Z")  & GMDateTime_c <=as_datetime("2013-07-21T02:30:00Z")) ~ "Resident",
  (Id=="MISS003" & GMDateTime_c >=as_datetime("2013-07-21T02:30:00Z")  & GMDateTime_c <=as_datetime("2013-11-14T15:00:00Z")) ~ "Resident",
  (Id=="MISS004" & GMDateTime_c >=as_datetime("2012-11-15T15:01:00Z")  & GMDateTime_c <=as_datetime("2014-02-05T03:01:00Z")) ~ "Resident",
  (Id=="MISS007" & GMDateTime_c >=as_datetime("2013-07-03T12:01:00Z")  & GMDateTime_c <=as_datetime("2013-12-25T04:31:00Z")) ~ "Resident",
  (Id=="MISS011" & GMDateTime_c >=as_datetime("2013-10-24T12:03:00Z")  & GMDateTime_c <=as_datetime("2014-02-20T04:29:00Z")) ~ "Resident",
  (Id=="MISS011" & GMDateTime_c >=as_datetime("2014-02-20T04:30:00Z")  & GMDateTime_c <=as_datetime("2014-04-30T10:33:00Z")) ~ "Transient",
  (Id=="MISS011" & GMDateTime_c >=as_datetime("2014-04-30T10:33:00Z")  & GMDateTime_c <=as_datetime("2014-06-25T04:30:00Z")) ~ "Resident",
  (Id=="MISS013" & GMDateTime_c >=as_datetime("2013-07-24T11:07:00Z")  & GMDateTime_c <=as_datetime("2014-07-18T04:31:00Z")) ~ "Resident",
  (Id=="MISS014" & GMDateTime_c >=as_datetime("2013-07-29T00:08:00Z")  & GMDateTime_c <=as_datetime("2013-09-09T06:01:00Z")) ~ "Resident",
  (Id=="MISS014" & GMDateTime_c >=as_datetime("2013-09-09T06:01:00Z")  & GMDateTime_c <=as_datetime("2013-11-23T21:01:00Z")) ~ "Transient",
  (Id=="MISS014" & GMDateTime_c >=as_datetime("2013-11-23T21:01:00Z")  & GMDateTime_c <=as_datetime("2014-04-28T22:59:00Z")) ~ "Resident",
  (Id=="MISS018" & GMDateTime_c >=as_datetime("2013-10-10T15:01:00Z")  & GMDateTime_c <=as_datetime("2014-07-15T01:30:00Z")) ~ "Transient",
  (Id=="MISS020" & GMDateTime_c >=as_datetime("2013-11-08T18:00:00Z")  & GMDateTime_c <=as_datetime("2014-08-28T06:00:00Z")) ~ "Resident",
  (Id=="MISS021" & GMDateTime_c >=as_datetime("2014-02-05T15:05:00Z")  & GMDateTime_c <=as_datetime("2014-02-22T07:30:00Z")) ~ "Transient",
  (Id=="MISS022" & GMDateTime_c >=as_datetime("2014-03-21T15:01:00Z")  & GMDateTime_c <=as_datetime("2015-03-25T12:02:00Z")) ~ "Transient",
  (Id=="URBAN001" & GMDateTime_c >=as_datetime("2014-06-06T13:33:00Z")  & GMDateTime_c <=as_datetime("2015-06-18T12:00:00Z")) ~ "Resident",
  (Id=="URBAN002" & GMDateTime_c >=as_datetime("2014-07-18T13:31:00Z")  & GMDateTime_c <=as_datetime("2015-07-08T04:31:00Z")) ~ "Resident",
  (Id=="URBAN005" & GMDateTime_c >=as_datetime("2016-04-04T22:31:00Z")  & GMDateTime_c <=as_datetime("2016-10-14T16:30:00Z")) ~ "Transient",
  (Id=="URBAN006" & GMDateTime_c >=as_datetime("2016-11-11T00:03:00Z")  & GMDateTime_c <=as_datetime("2016-11-16T13:31:00Z")) ~ "Resident",
  (Id=="URBAN007" & GMDateTime_c >=as_datetime("2017-05-10T00:07:00Z")  & GMDateTime_c <=as_datetime("2017-06-04T15:00:00Z")) ~ "Transient",
  (Id=="URBAN008" & GMDateTime_c >=as_datetime("2017-06-01T03:02:00Z")  & GMDateTime_c <=as_datetime("2017-06-27T21:01:00Z")) ~ "Transient",
  (Id=="URBAN008" & GMDateTime_c >=as_datetime("2017-06-27T21:01:00Z")  & GMDateTime_c <=as_datetime("2019-12-12T08:33:00Z")) ~ "Resident",
  (Id=="TOR01" & GMDateTime_c >=as_datetime("2011-05-31T17:01:00Z")  & GMDateTime_c <=as_datetime("2011-10-27T14:01:00Z")) ~ "Resident",
  (Id=="URBAN009" & GMDateTime_c >=as_datetime("2017-11-16T18:03:00Z")  & GMDateTime_c <=as_datetime("2017-11-27T16:30:00Z")) ~ "Transient",
  (Id=="URBAN009" & GMDateTime_c >=as_datetime("2017-11-27T16:30:00Z")  & GMDateTime_c <=as_datetime("2019-12-14T04:31:00Z")) ~ "Resident",
  (Id=="URBAN010" & GMDateTime_c >=as_datetime("2018-01-16T00:01:00Z")  & GMDateTime_c <=as_datetime("2018-01-17T12:00:00Z")) ~ "Transient",
  (Id=="URBAN010" & GMDateTime_c >=as_datetime("2018-01-17T12:00:00Z")  & GMDateTime_c <=as_datetime("2018-08-31T15:01:30Z")) ~ "Resident",
  (Id=="URBAN011" & GMDateTime_c >=as_datetime("2020-12-22T23:56:00Z")  & GMDateTime_c <=as_datetime("2021-05-01T23:46:00Z")) ~ "Resident",
  (Id=="URBAN012" & GMDateTime_c >=as_datetime("2019-12-19T21:05:00Z")  & GMDateTime_c <=as_datetime("2020-03-17T09:03:00Z")) ~ "Resident",
  (Id=="URBAN012" & GMDateTime_c >=as_datetime("2020-03-17T09:03:00Z")  & GMDateTime_c <=as_datetime("2020-03-18T09:01:00Z")) ~ "Transient",
  (Id=="URBAN013" & GMDateTime_c >=as_datetime("2019-12-19T21:04:00Z")  & GMDateTime_c <=as_datetime("2020-06-09T18:01:00Z")) ~ "Transient",
  (Id=="URBAN014" & GMDateTime_c >=as_datetime("2020-01-23T06:02:18Z")  & GMDateTime_c <=as_datetime("2021-05-01T23:45:00Z")) ~ "Resident",
  (Id=="URBAN015" & GMDateTime_c >=as_datetime("2020-02-12T00:03:00Z")  & GMDateTime_c <=as_datetime("2020-08-18T03:00:00Z")) ~ "Transient",
  (Id=="URBAN015" & GMDateTime_c >=as_datetime("2020-08-18T03:00:00Z")  & GMDateTime_c <=as_datetime("2021-01-09T18:01:00Z")) ~ "Resident",
  (Id=="URBAN016" & GMDateTime_c >=as_datetime("2020-08-02T06:00:29Z")  & GMDateTime_c <=as_datetime("2021-05-01T23:45:00Z")) ~ "Transient",
  (Id=="URBAN017" & GMDateTime_c >=as_datetime("2020-09-29T21:01:00Z")  & GMDateTime_c <=as_datetime("2021-02-18T15:01:00Z")) ~ "Resident",
  (Id=="URBAN018" & GMDateTime_c >=as_datetime("2020-09-29T20:59:00Z")  & GMDateTime_c <=as_datetime("2021-05-01T23:45:13Z")) ~ "Transient",
  (Id=="URBAN019" & GMDateTime_c >=as_datetime("2020-11-11T20:58:24Z")  & GMDateTime_c <=as_datetime("2021-01-09T09:00:00Z")) ~ "Transient",
  (Id=="URBAN019" & GMDateTime_c >=as_datetime("2021-01-09T09:00:00Z")  & GMDateTime_c <=as_datetime("2021-05-01T23:45:29Z")) ~ "Resident",
  (Id=="URBAN020" & GMDateTime_c >=as_datetime("2020-11-11T20:58:00Z")  & GMDateTime_c <=as_datetime("2021-05-01T23:45:29Z")) ~ "Transient",
  (Id=="URBAN021" & GMDateTime_c >=as_datetime("2020-11-19T21:01:00Z")  & GMDateTime_c <=as_datetime("2020-12-31T18:01:00Z")) ~ "Transient",
  (Id=="URBAN022" & GMDateTime_c >=as_datetime("2020-12-22T20:58:30Z")  & GMDateTime_c <=as_datetime("2021-02-19T09:00:18Z")) ~ "Transient",
  (Id=="URBAN023" & GMDateTime_c >=as_datetime("2021-02-04T00:03:00Z")  & GMDateTime_c <=as_datetime("2021-05-01T23:45:16Z")) ~ "Resident",
  (Id=="URBAN024" & GMDateTime_c >=as_datetime("2021-02-24T01:30:38Z")  & GMDateTime_c <=as_datetime("2021-05-01T23:45:05Z")) ~ "Transient",
  (Id=="URBAN025" & GMDateTime_c >=as_datetime("2021-04-10T00:01:00Z")  & GMDateTime_c <=as_datetime("2021-05-11T22:30:00Z")) ~ "Transient",
  TRUE~ "1"))


### DEFINE SEASONS ###

coyotes1$LOCDateTime_c <- with_tz(coyotes1$GMDateTime, tz = "America/Toronto")

## function for defining season based on seaosn start and end dates
define.season_LOC <- function(date, start.day, end.day) {
  Year <- as.character(year(date))
  start <- as.POSIXct(paste(start.day, Year, sep = "-"),
                      format = "%d-%m-%Y", tz = "America/Toronto")
  end <- as.POSIXct(paste(paste(end.day, Year, sep = "-"), "23:59:59", sep = " "),
                    format = "%d-%m-%Y %H:%M:%S", tz = "America/Toronto")
  ifelse(date >= start &
           date <= end,
         select.time <- 1,
         select.time <- 0)
  
}

## ADD variables created so far as binary information to each step.
## Each level of a variable is a column and it is either 1 or 0 if that step falls in that level. 

coyotes1 <- coyotes1 %>% 
  mutate(
    breeding = define.season_LOC(LOCDateTime_c, "01-01", "30-04"), 
    pup_rearing = define.season_LOC(LOCDateTime_c, "01-05", "31-08"),
    dispersal = define.season_LOC(LOCDateTime_c, "01-09", "31-12"),
    resident = case_when(life_stage == "Resident" ~ 1, TRUE  ~0),
    transient = case_when(life_stage == "Transient" ~ 1, TRUE  ~0),
    veg_covered = define.season_LOC(LOCDateTime_c, "01-04", "30-09")) #summer


#bypassing problems: opposite time that is not covered is bare...crossing over years was being an issue somehow. 
coyotes1 <- coyotes1 %>% 
  mutate(veg_bare = case_when(veg_covered == 1 ~ 0, #winter
                              TRUE~1))

##DEFINE day/night 
# and add as binary variable to the previous table


coyotes1 <- coyotes1 %>% 
  mutate(day = ifelse(sunriset(matrix(c(-79.39276124616116, 43.65259573314474), nrow = 1), 
                               LOCDateTime_c, direction = "sunrise", 
                               POSIXct.out = TRUE)$time <= LOCDateTime_c &
                        sunriset(matrix(c(-79.39276124616116, 43.65259573314474), nrow = 1), 
                                 LOCDateTime_c, direction = "sunset", 
                                 POSIXct.out = TRUE)$time >= LOCDateTime_c, 1, 0))

###make categorical variables with the previous information 

coyotes1 <- coyotes1 %>% 
  mutate(behavioral_season = ifelse(breeding==1,"breeding", 
                                    ifelse(dispersal==1,"dispersal",
                                           ifelse(pup_rearing==1,"pup_rearing",NA))))

coyotes1 <- coyotes1 %>% 
  mutate(vegetation_season = ifelse(veg_bare==1,"veg_bare", 
                                    ifelse(veg_covered==1,"veg_covered",NA)))

coyotes1 <- coyotes1 %>% 
  mutate(daytime = ifelse(day==1,"day", 
                          ifelse(day==0,"night",NA)))

### ADD Sex and Age ## information obtained from MNRF notes when collaring

coy_lh <- read_csv("coyotes_lifehistory_mayupdate.csv")

coyotes_df <- left_join(coyotes1, coy_lh, by="Id", copy = TRUE)
#View(coyotes_df)


coyotes_df <- coyotes_df%>% select(-Sex.y, -Age.y, -Sex.x.x, -Age.x.x, -Age.y.y, -Sex.y.y)
write.csv(coyotes_df, file="coyotes_update_fix_variables.csv", na="NA")




#Convert coyotes into an sf object
coyotes_sf <- st_as_sf(coyotes_df, coords = c("Longitude", "Latitude"), 
                       crs = st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>% 
  st_transform(crs = 26917)



write_rds(coyotes_sf, file="coyotes_sf.rds") # this object is used in the next step




