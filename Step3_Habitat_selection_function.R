### Step 3: Habitat selection model  ###

#This script contains the top models only, selected through the specification process in "Model_specification_process.R"

## libraries ##

library(sf)
library(tidyverse)
library(glmmTMB)
library(buildmer)
library(dplyr)
library(regclass)
library(car)
library(bbmle)
library(AICcmodavg)
library(sjPlot)

## load data from Step 2, included in repository

readRDS("coyote_steps_clean_2.rds")  -> coyote_steps

##### LANDSCAPE-ONLY MODEL####

habmodel22 <- glmmTMB(used ~ -1  + NDVI_end + BUP_end + DT_PS_end  + DT_NT_end + POP_end + DT_MT_end + DT_LT_end + DT_HT_end +
                                         (1|animal) + 
                                         (0+NDVI_end|animal) + 
                                         (0+BUP_end|animal) +
                                         (0+DT_PS_end|animal) +
                                         (0+DT_NT_end|animal) +  
                                         (0+POP_end|animal) +
                                         (0+DT_MT_end|animal) + 
                                         (0+DT_LT_end|animal) + 
                                         (0+DT_HT_end|animal) + 
                                         (1|step), data = coyote_steps, 
                                       map = list(theta = factor(c(1:9, NA))),
                                       start = list(theta = c(rep(0, 9), log(1e3))),
                                       family = poisson)


saveRDS(NDVI.BUP.PS.NT.POP.MT.LT.HT, "habmodel22.rds")
capture.output(summary(NDVI.BUP.PS.NT.POP.MT.LT.HT),file="habmodel22.doc")

## DIEL CYCLE MODEL ####

daymodel37<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                      DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + 
                      DT_LT_end:day + NDVI_end:day + DT_NT_end:day + BUP_end:day + DT_PS_end:day + DT_HT_end:day +  POP_end:day  +  DT_MT_end:day +
                      (1|animal) + 
                      (0+NDVI_end|animal) + 
                      (0+POP_end|animal) +
                      (0+BUP_end|animal) + 
                      (0+DT_LT_end|animal) + 
                      (0+DT_MT_end|animal) + 
                      (0+DT_HT_end|animal) + 
                      (0+DT_NT_end|animal) + 
                      (0+DT_PS_end|animal) +
                      (0+DT_LT_end:day|animal) + (0+NDVI_end:day |animal) + (0+DT_NT_end:day |animal) + 
                      (0+BUP_end:day|animal) + (0+DT_PS_end:day |animal) + (0+DT_HT_end:day |animal) + (0+POP_end:day |animal) + (0+DT_MT_end:day |animal)  + (1|step), data = coyote_steps,  
                    map = list(theta = factor(c(1:17, NA))),
                    start = list(theta = c(rep(0, 17), log(1e3))),
                    family = poisson)
saveRDS(daymodel37, "daymodel37.rds") 
capture.output(summary(daymodel37),file="daymodel37.doc")

### SOCIAL STATUS MODEL ####

transientmodelac <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                              DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                              BUP_end:transient + NDVI_end:transient +
                              (1|animal) + 
                              (0+NDVI_end|animal) + 
                              (0+POP_end|animal) +
                              (0+BUP_end|animal) + 
                              (0+DT_LT_end|animal) + 
                              (0+DT_MT_end|animal) + 
                              (0+DT_HT_end|animal) + 
                              (0+DT_NT_end|animal) + 
                              (0+DT_PS_end|animal) + 
                              (0+BUP_end:transient|animal) + 
                              (0+NDVI_end:transient|animal) + 
                              (1|step), data = coyote_steps, 
                            map = list(theta = factor(c(1:11, NA))),
                            start = list(theta = c(rep(0, 11), log(1e3))),
                            family = poisson)

saveRDS(transientmodelac, "transientmodelac.rds")
capture.output(summary(transientmodelac),file="transientmodelac.doc")

### BEHAVIORAL SEASON MODEL ####

breedingmodel25<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                           DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                           DT_NT_end:breeding + DT_NT_end:dispersal+
                           BUP_end:breeding + BUP_end:dispersal+
                           NDVI_end:breeding +NDVI_end:dispersal+  
                           DT_PS_end:breeding +DT_PS_end:dispersal+     
                           DT_HT_end:breeding +DT_HT_end:dispersal+ 
                           (1|animal) + 
                           (0+NDVI_end|animal) + 
                           (0+POP_end|animal) +
                           (0+BUP_end|animal) + 
                           (0+DT_LT_end|animal) + 
                           (0+DT_MT_end|animal) + 
                           (0+DT_HT_end|animal) + 
                           (0+DT_NT_end|animal) + 
                           (0+DT_PS_end|animal) + 
                           (0+DT_NT_end:dispersal|animal) +  (0+DT_NT_end:breeding|animal) + 
                           (0+BUP_end:dispersal|animal) + (0+BUP_end:breeding|animal) + 
                           (0+NDVI_end:dispersal|animal) + (0+NDVI_end:breeding|animal) + 
                           (0+DT_PS_end:dispersal|animal) + (0+DT_PS_end:breeding|animal) + 
                           (0+DT_HT_end:dispersal|animal) + (0+DT_HT_end:breeding|animal) +  
                           (1|step), data = coyote_steps, 
                         map = list(theta = factor(c(1:19, NA))),
                         start = list(theta = c(rep(0, 19), log(1e3))),
                         family = poisson)
saveRDS(breedingmodel25, "breedingmodel25.rds")
capture.output(summary(breedingmodel25),file="breedingmodel25.doc")


### CLIMATE SEASON MODEL ####

wintermodel15<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                         DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + 
                         DT_NT_end:winter+
                         NDVI_end:winter+
                         DT_MT_end:winter +
                         (1|animal) + 
                         (0+NDVI_end|animal) + 
                         (0+POP_end|animal) +
                         (0+BUP_end|animal) + 
                         (0+DT_LT_end|animal) + 
                         (0+DT_MT_end|animal) + 
                         (0+DT_HT_end|animal) + 
                         (0+DT_NT_end|animal) + 
                         (0+DT_PS_end|animal) +(0+DT_NT_end:winter|animal) +  (0+NDVI_end:winter|animal) +  (0+DT_MT_end:winter|animal) + (1|step), data = coyote_steps, 
                       map = list(theta = factor(c(1:12, NA))),
                       start = list(theta = c(rep(0, 12), log(1e3))),
                       family = poisson)
saveRDS(wintermodel15, "wintermodel15.rds")
capture.output(summary(wintermodel15),file="wintermodel15.doc")


### SEX MODEL ####


malemodelb <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                        DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                        POP_end:male +
                        (1|animal) + 
                        (0+NDVI_end|animal) + 
                        (0+POP_end|animal) +
                        (0+BUP_end|animal) + 
                        (0+DT_LT_end|animal) + 
                        (0+DT_MT_end|animal) + 
                        (0+DT_HT_end|animal) + 
                        (0+DT_NT_end|animal) + 
                        (0+DT_PS_end|animal) + 
                        (0+POP_end:male|animal) + 
                        (1|step), data = coyote_steps, 
                      map = list(theta = factor(c(1:10, NA))),
                      start = list(theta = c(rep(0, 10), log(1e3))),
                      family = poisson)
saveRDS(malemodelb, "malemodelb.rds")
capture.output(summary(malemodelb),file="malemodelb.doc")

### AGE MODEL ####

adultmodelb <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                         DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                         POP_end:adult +
                         (1|animal) + 
                         (0+NDVI_end|animal) + 
                         (0+POP_end|animal) +
                         (0+BUP_end|animal) + 
                         (0+DT_LT_end|animal) + 
                         (0+DT_MT_end|animal) + 
                         (0+DT_HT_end|animal) + 
                         (0+DT_NT_end|animal) + 
                         (0+DT_PS_end|animal) + 
                         (0+POP_end:adult|animal) + 
                         (1|step), data = coyote_steps, 
                       map = list(theta = factor(c(1:10, NA))),
                       start = list(theta = c(rep(0, 10), log(1e3))),
                       family = poisson)
saveRDS(adultmodelb, "adultmodelb.rds")
capture.output(summary(adultmodelb),file="adultmodelb.doc")
