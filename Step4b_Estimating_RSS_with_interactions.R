### Step 4b: Relative Selection Strength estimates from HSF with interactions #### 

# This script uses daymodel22, but the same was repeated for all interaction models.
# Note: the breeding model is composed of THREE levels,therefore there is a commented section within the code to accomodate this. 
# This script is specific for models with interactions as the RSS function changes
# The estimate follows Avgar et al., 2017 mathematical formula


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


############RSS############ 

###call model rds object 
#habmodel22 <- readRDS("habmodel22_clean_NDVI.BUP.PS.NT.POP.MT.LT.HT_nocovid.rds")

daymodel37 <- readRDS("daymodel22.rds")
#daymodel37 <- readRDS("daymodel37_clean.rds")

###dataset####
coyote_steps <- readRDS("coyote_steps_clean.rds")


####Relative Selection Strength estimates###

# To estimate relative selection strength, for each variable, 
# we generate a dummy table recreating the value range of the variable while all others remain constant
#and use the coefficient of the model to estimate the relative selection 
#strength across each value

###define "m" - with the model being used for the RSS #######
m <- fixef(daymodel37)

##making this inefficient code a loop ###
land_var <- c('NDVI', 'POP', 'BUP', 'DT_LT', 'DT_MT', 'DT_HT', 'DT_NT', 'DT_PS')
i=2
for (i in 1:length(land_var)){
  ### create dummy tables for estimate rss
  s2 <-data.frame(matrix(nrow=200, ncol=0))
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
    s1$NDVI_nd = seq(from = -2.8, to =2.7, length.out = 200) #min max obtained from min max of landscape variable scaled raster layers 
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
  
  # for models with an interaction component we have to repeated the table with range values, adding a column describing the level of the landscape variable they represent,
  # and we just cbind these two tables.
  
  # add factor levels (2 level factors)
  s1$day = factor("0", levels = levels(as.factor(coyote_steps$day))) #0 or 1  and has t be a factor #column name must be the landscape_variable name as per the model
  s1b <- s1
  s1b$day = factor("1", levels = levels(as.factor(coyote_steps$day))) #0 or 1  and has t be a factor #column name must be the landscape_variable name as per the model
  s1<- rbind(s1, s1b)
  
  #see commented out sections for models with THREE level interaction factors, with more levels the logic remains, just add whatever was added for the extra level
  # note the interaction levels are always in binary form, so it is the combination of the columns that reflects the level of the interaction factor
  ### for models with THREE level interaction factors (eg. behavioral season model) do this instead -- we need to duplicate the table one more time and add a second column
  # s1$breeding = factor("0", levels = levels(as.factor(coyote_steps$breeding)))
  # s1$dispersal = factor("0", levels = levels(as.factor(coyote_steps$dispersal)))
  # s1b<- s1
  # s1$breeding = factor("1", levels = levels(as.factor(coyote_steps$breeding)))
  # s1$dispersal = factor("0", levels = levels(as.factor(coyote_steps$dispersal)))
  # s1c<-s1
  # s1$breeding = factor("0", levels = levels(as.factor(coyote_steps$breeding)))
  # s1$dispersal = factor("1", levels = levels(as.factor(coyote_steps$dispersal)))
  # s1<- rbind(s1, s1b, s1c)
  # 
   
  s1$log_rss_recalc<-0 
  s1 <- remove_rownames(s1) 
  
  #summarize the table to what we need for each variable - TWO LEVEL INTERCTION FACTORS
  table_juice <- s1%>% select(contains(land_var[i]), "day") #make sure you change the column name to the interaction factor name of whatever model you are running, e.g. "transient" in the transientmodel, the loop can be improved to include a nested loop for each model where this changes to model[i]
  names(table_juice)<-c("var", "cos", "level")
  
  #for THREE level factors, add the extra two columns
  #table_juice <- s1%>% select(contains(land_var[i]), "breeding", "dispersal")
  #names(table_juice)<-c("var", "cos", "level1", "level2")
  
  
  #Estimate RSS
  #no interaction from script Step4a:
  # for (j in 1:length(s1$log_rss_recalc)){
  #   s1$log_rss_recalc[j] <- m$cond[grepl(land_var[i], names(m$cond))]*(table_juice$var[j]-table_juice$cos[j])
  # }
  
  #for models with interaction factor -- TWO level interaction factors:
  for (j in 1:length(table_juice$log_rss_recalc)){
    if (table_juice$level[j] == "0") {
      table_juice$log_rss_recalc[j] <- 
        (m$cond[grepl(land_var[i], names(m$cond))][[1]]*(table_juice$var[j]-table_juice$cos[j]))
    } else if (table_juice$level[j] == "1") {
      table_juice$log_rss_recalc[j] <- 
        ((table_juice$var[j]-table_juice$cos[j])*(m$cond[grepl(land_var[i], names(m$cond))][[1]] #predictor coeff
                                                +m$cond[grepl(land_var[i], names(m$cond))][[2]] #interaction coefficient
        ))
    }
    
    ## FOR THREE LEVEL INTERACTION FACTORS - estimate RSS this way instead #note level 3 is implied in the level1&level2 = 0
    for (j in 1:length(table_juice$log_rss_recalc)){ #s1 and table_juice have the same number of rows, they just vary in columns, since table_juice is the clean version of s1, so this looping works because of that
      if (table_juice$level1[j] == "0" && table_juice$level2[j] == "0") {
        table_juice$log_rss_recalc[j] <- 
          (m$cond[grepl(land_var[i], names(m$cond))][[1]]*(table_juice$var[j]-table_juice$cos[j]))
      } else if (table_juice$level1[j] == "1" ) { #breeding =1
        table_juice$log_rss_recalc[j] <- 
          ((table_juice$var[j]-table_juice$cos[j])*(m$cond[grepl(land_var[i], names(m$cond))][[1]]+
                                                      m$cond[grepl(land_var[i], names(m$cond))][[2]]))
      } else if (table_juice$level2[j] == "1" ) { #dispersal = 1
        table_juice$log_rss_recalc[j] <- 
          ((table_juice$var[j]-table_juice$cos[j])*(m$cond[grepl(land_var[i], names(m$cond))][[1]]+
                                                      m$cond[grepl(land_var[i], names(m$cond))][[3]]))
      }
    }
    
  }
  

  #linear regression -- TWO LEVEL INTERACTION FACTORS
  lm_rss<-lm(log_rss_recalc ~ var*level,  data = table_juice) # linear regression of log RSS to find invert and find resistance values
  
  ## FOR THREE LEVEL FACTORS:
  # lm_rss<-lm(log_rss_recalc ~ var*level1 + var*level2,  data = table_juice) # linear regression of log RSS to find invert and find resistance values
  
  #collect stats for generating resistance map in the next step
  summary(lm_rss)$coefficients -> Stats_lm_rss #coefficient of the linear regression
  
  #for TWO level interaction factors:
  intercept  <- Stats_lm_rss[1,1]
  coefficient  <- Stats_lm_rss[2,1]
  coefficient2  <- Stats_lm_rss[3,1]
  interaction  <- Stats_lm_rss[4,1]
  stats_table <- rbind(intercept, coefficient, coefficient2, interaction)
  colnames(stats_table)<-paste0(land_var[i])
  
  ## FOR THREE LEVEL INTERACTION FACTOR do this instead:
  # intercept  <- Stats_lm_rss[1,1]
  # coefficient  <- Stats_lm_rss[2,1] #pup-rearing #you can change these object names to reflect this if numbered coefficient gets confusing, you'll need this info for interpretation
  # coefficient2.1  <- Stats_lm_rss[3,1] #breeding
  # coefficient2.2  <- Stats_lm_rss[4,1] #dispersal
  # interaction1  <- Stats_lm_rss[5,1]
  # interaction2  <- Stats_lm_rss[6,1]
  # stats_table <- rbind(intercept, coefficient, coefficient2.1, coefficient2.2, interaction1, interaction2)
  # colnames(stats_table)<-paste0(land_var[i])
  
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
  
  ## Plotting models with THREE level interactin factors
  # ##to plot THREE LEVEL factors you can add a column with interpretation of what the binary columns mean"
  # table_juice <- table_juice %>% mutate(
  #   level = ifelse(level1=="1" & level2 =="0", "breeding",
  #                  ifelse(level1=="0" & level2 =="1", "dispersal",
  #                         "pup-rearing")))
  
  #you can also do the previous part for models with TWO level interaction factors if you dont want to worry about changing the legend, for example if looping through models

  
  ggplot(table_juice, aes(x = var, y = log_rss_recalc, color=level )) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
    xlab("Landscape variable range (scaled)") +
    ylab("logRSS") +
    theme_bw()+
    scale_y_continuous()+
    ggtitle(land_var[i])+
    scale_color_discrete(name="Daytime", labels = c("night", "day"))

  print(ggplot(table_juice, aes(x = var, y = exp(log_rss_recalc), color =level )) +
    geom_line(size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
    xlab("Landscape variable range (scaled)") +
    ylab("RSS") +
    theme_bw()+
    scale_y_continuous()+
    ggtitle(land_var[i])+
    scale_color_discrete(name="Daytime", labels = c("night", "day"))) #comment this out for three level factors, so you see the names we created above
}

#note rds objects will be used in next step



