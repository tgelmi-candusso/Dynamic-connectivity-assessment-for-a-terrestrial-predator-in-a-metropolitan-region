### Step 4: Relative Selection Strength estimates #### 

# This script uses habmodel22, but the same was repeated for all models 
# This script is specific for models without interactions as the RSS function changes
# The estimate follows Avgar et al., 2017 function.

## libraries ##

library(ggplot2)
library(raster)
library(lubridate)
library(amt)
require(dplyr)
library(parallel)
library(glmmTMB)
library(sf)
library(tidyverse)
library(buildmer)

##set working directory
getwd()
setwd("C:/Users/tizge/Documents/Toronto project/Coyote tracking by MNRF/Resource Selection Function/ResourceSelectionFunction")

####raster###
##call rasters for landscape variables
NDVI <- readRDS("NDVI_clean.rds")
POP <- readRDS("POP_clean.rds")
BUP <- readRDS("BUP_clean.rds")
DT_LT <- readRDS("DT_LT_clean.rds")
DT_MT <- readRDS("DT_MT_clean.rds")
DT_HT <- readRDS("DT_HT_clean.rds")
DT_NT <- readRDS("DT_NT_clean.rds")
DT_PS <- readRDS("DT_PS_clean.rds")


############RSS############ 

###call model rds object 
#habmodel22 <- readRDS("habmodel22_clean_NDVI.BUP.PS.NT.POP.MT.LT.HT_nocovid.rds")

habmodel22 <- readRDS("habmodel22.rds")

###dataset####
coyote_steps <- readRDS("coyote_steps_clean.rds")


####Relative Selection Strength estimates###

# To estimate relative selection strength, for each variable, 
# we generate a dummy table recreating the value range of the variable while all others remain constant
#and use the coefficient of the model to estimate the relative selection 
#strength across each value
# not the most efficient code

###define "m" - with the model being used for the RSS #######
m <- fixef(habmodel22)

##making this inefficient code a loop ###
land_var <- c('NDVI', 'POP', 'BUP', 'DT_LT', 'DT_MT', 'DT_HT', 'DT_NT', 'DT_PS')
i=2
for (i in 1:length(land_var)){
  ### create dummy tables for estimate rss
  s2 <-data.frame(matrix(, nrow=200, ncol=0))
  s2 <- s2 %>% mutate(
    NDVI_nd =  mean(coyote_steps$NDVI_end),
    POP_end =  mean(coyote_steps$POP_end),
    BUP_end = mean(coyote_steps$BUP_end),
    DT_LT_n = mean(coyote_steps$DT_LT_end),
    DT_MT_n = mean(coyote_steps$DT_MT_end),
    DT_HT_n = mean(coyote_steps$DT_HT_end),
    DT_NT_n = mean(coyote_steps$DT_NT_end),
    DT_PS_n = mean(coyote_steps$DT_PS_end))

    s1 <- s2
  
    # generate the columns with values within the range of each landscape variable (scaled - as in model)
  if (land_var[i] == 'NDVI'){
    s1$NDVI_nd = seq(from = -2.8, to =2.7, length.out = 200) 
    s1$NDVI_x2 <- s2$NDVI_nd
  } else if(land_var[i] == 'POP'){
    s1$POP_end = seq(from = -0.7, to =15.8, length.out = 200)
    s1$POP_x2 <- s2$POP_end
  } else if(land_var[i] == 'BUP'){
    s1$BUP_end = seq(from = -1.8, to =3.1, length.out = 200)
    s1$BUP_x2 <- s2$BUP_end
  } else if(land_var[i] == 'DT_LT'){
    s1$DT_LT_n = seq(from = -0.6, to =11.1, length.out = 200)
    s1$DT_LT_n_x2 <- s2$DT_LT_n
  } else if(land_var[i] == 'DT_MT'){
    s1$DT_MT_n = seq(from = -0.8, to = 12, length.out = 200)
    s1$DT_MT_n_x2 <- s2$DT_MT_n
  } else if(land_var[i] == 'DT_HT'){
    s1$DT_HT_n = seq(from = -0.8, to = 11.2, length.out = 200)
    s1$DT_HT_n_x2 <- s2$DT_HT_n
  }else if(land_var[i] == 'DT_NT'){
    s1$DT_NT_n = seq(from = -0.9, to = 10.4, length.out = 200)
    s1$DT_NT_n_x2 <- s2$DT_NT_n
  }else if(land_var[i] == 'DT_PS'){
    s1$DT_PS_n = seq(from = -1, to = 5, length.out = 200)
    s1$DT_PS_n_x2 <- s2$DT_PS_n
  }

  s1$log_rss_recalc<-0 
  s1 <- remove_rownames(s1) 
  
  #summarize the table to what we need for each variable
  table_juice <- s1[grepl(land_var[i], names(s1))]
  names(table_juice)<-c("var", "cos")
  
  #Estimate RSS
  for (j in 1:length(s1$log_rss_recalc)){
    s1$log_rss_recalc[j] <- m$cond[grepl(land_var[i], names(m$cond))]*(table_juice$var[j]-table_juice$cos[j])
  }
  
  #collect the regression coefficients of the log RSS curve
  #simplify table to run the regression
  table_juice$log_rss_recalc <- s1$log_rss_recalc
  #linear regression
  lm_rss<-lm(log_rss_recalc ~ var,  data = table_juice) # linear regression of log RSS to find invert and find resistance values
  #collect stats for generating resistance map in the next step
  summary(lm_rss)$coefficients -> Stats_lm_rss #coefficient of the linear regression
  
  intercept  <- Stats_lm_rss[1,1]
  coefficient  <- Stats_lm_rss[2,1]
  stats_table <- rbind(intercept, coefficient)
  colnames(stats_table)<-paste0(land_var[i])
  
  if(i==1){
    stats_table_all <- stats_table
  }else {
    stats_table_all <- cbind(stats_table_all, stats_table)
  }
  
  #RDS object, can also be a csv
  write_rds(stats_table_all, "stats_table_all.rds")
  
  ##write table with rss of all landscape variables for  plotting
  table_juice$variable <- paste(land_var[i])
  if(i==1){
    all_rss_tables <- table_juice
  }else{
    all_rss_tables <- rbind(all_rss_tables, table_juice)
  }
  
  ##RDS OBJECT can also be a csv
  write_rds(all_rss_tables, "all_rss_tables_for_plots.rds")
  
  ### plots - single variable, just for checking everything is working
  
  ggplot(table_juice, aes(x = var, y = log_rss_recalc )) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
    xlab("Distance to public service linear features (scaled)") +
    ylab("logRSS") +
    theme_bw()+
    scale_y_continuous()+
    ggtitle(land_var[i])
  #+ geom_ribbon(aes(ymin = lwr, ymax = upr),fill = "gray80", alpha = 0.1)
  
  
  ggplot(table_juice, aes(x = var, y = exp(log_rss_recalc) )) +
    geom_line(size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
    xlab("Distance to public service linear features (scaled)") +
    ylab("RSS") +
    theme_bw()+
    scale_y_continuous()+
    ggtitle(land_var[i])
}

#note rds objects will be used in next step



