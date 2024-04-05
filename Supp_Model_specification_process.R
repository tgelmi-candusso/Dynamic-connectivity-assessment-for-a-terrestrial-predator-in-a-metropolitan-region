### Model specification process ###
# This code runs the habitat selection function repeatedly hierarchically adding variables, following our model specification protocol

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

#set working directory
setwd("C:/Users/tizge/Documents/Toronto project/Coyote tracking by MNRF/Resource Selection Function/ResourceSelectionFunction")

readRDS("coyote_steps_2.rds")  -> coyote_steps



#####VIF#### #check for correlation between variables
vf1 <- lm(coyote_steps$NDVI_end~coyote_steps$BUP_end)
vf2 <- lm(coyote_steps$NDVI_end~coyote_steps$POP_end)
vf3 <- lm(coyote_steps$POP_end~coyote_steps$BUP_end)
vf4 <- lm(coyote_steps$NDVI_end~coyote_steps$POP_end+coyote_steps$BUP_end)
vf5 <- lm(coyote_steps$POP_end~coyote_steps$NDVI_end+coyote_steps$BUP_end)
vf6 <- lm(coyote_steps$BUP_end~coyote_steps$NDVI_end+coyote_steps$POP_end)


summary(vf1) 
vf1_r2 <- 0.1291
vf1_res <- 1/(1-vf1_r2)

summary(vf2)
vf2_r2 <- 0.02185
vf2_res <- 1/(1-vf2_r2)

summary(vf3) 
vf3_r2 <- 0.1443
vf3_res <- 1/(1-vf3_r2)

summary(vf4) 
vf4_r2 <- 0.1292
vf4_res <- 1/(1-vf4_r2)

summary(vf5) 
vf5_r2 <-0.1444
vf5_res <- 1/(1-vf5_r2)

summary(vf6) 
vf6_r2 <- 0.2382
vf6_res <- 1/(1-vf6_r2)

VIFS <- c(vf1_res, vf2_res, vf3_res, vf4_res, vf5_res, vf6_res)
##1 < VIF < 5 = Moderately correlated


####Null Model####
#Null Model
nullmodel <- glmmTMB(used ~ -1 +  (1|animal) + 
                       (1 | step), data = coyote_steps, 
                     map = list(theta = factor(c(1, NA))),
                     start = list(theta = c(0, log(1e3))),
                     family = poisson)

#summary(nullmodel)
#capture.output(#summary(nullmodel),file="nullmodel.doc")

##### Environmental Models #####  
#just habitat, so in my case just urban variables

NDVI_m <- glmmTMB(used ~ -1 + NDVI_end + (1|animal) + 
                    (0+NDVI_end|animal)  + 
                    (1 | step), data = coyote_steps, 
                  map = list(theta = factor(c(1:2, NA))),
                  start = list(theta = c(rep(0, 2), log(1e3))),
                  family = poisson) 

#summary(NDVI_m)
#capture.output(#summary(NDVI_m),file="NDVI.doc")



POP_m <- glmmTMB(used ~ -1 + POP_end + (1|animal)  + 
                   (0+POP_end|animal)+ 
                   (1 | step), data = coyote_steps, 
                 map = list(theta = factor(c(1:2, NA))),
                 start = list(theta = c(rep(0, 2), log(1e3))),
                 family = poisson) 
#summary(POP_m)
#capture.output(#summary(POP_m),file="POP.doc")


BUP_m <- glmmTMB(used ~ -1 + BUP_end + (1|animal)  +
                   (0+BUP_end|animal)+ 
                   (1 | step), data = coyote_steps, 
                 map = list(theta = factor(c(1:2, NA))),
                 start = list(theta = c(rep(0, 2), log(1e3))),
                 family = poisson)

#summary(BUP_m)
#capture.output(#summary(BUP_m),file="BUP.doc")


DTLT_m <- glmmTMB(used ~ -1 + DT_LT_end + (1|animal)  + 
                    (0+ DT_LT_end|animal)+ 
                    (1 | step), data = coyote_steps, 
                  map = list(theta = factor(c(1:2, NA))),
                  start = list(theta = c(rep(0, 2), log(1e3))),
                  family = poisson)
#summary(DTLT_m)
#capture.output(#summary(DTLT_m),file="DTLT.doc")


DTMT_m <- glmmTMB(used ~ -1 + DT_MT_end + (1|animal)  + 
                    (0+ DT_MT_end|animal)+ 
                    (1 | step), data = coyote_steps, 
                  map = list(theta = factor(c(1:2, NA))),
                  start = list(theta = c(rep(0, 2), log(1e3))),
                  family = poisson)
#summary(DTMT_m)
#capture.output(#summary(DTMT_m),file="DTMT.doc")

DTHT_m <- glmmTMB(used ~ -1 + DT_HT_end + (1|animal)  + 
                    (0+ DT_HT_end|animal)+ 
                    (1 | step), data = coyote_steps, 
                  map = list(theta = factor(c(1:2, NA))),
                  start = list(theta = c(rep(0, 2), log(1e3))),
                  family = poisson)

#summary(DTHT_m)
#capture.output(#summary(DTHT_m),file="DTHT.doc")



DTNT_m <- glmmTMB(used ~ -1 + DT_NT_end + (1|animal)  + 
                    (0+ DT_NT_end|animal)+ 
                    (1 | step), data = coyote_steps, 
                  map = list(theta = factor(c(1:2, NA))),
                  start = list(theta = c(rep(0, 2), log(1e3))),
                  family = poisson)

#summary(DTNT_m)
#capture.output(#summary(DTNT_m),file="DTNT.doc")


DTPS_m <- glmmTMB(used ~ -1 + DT_PS_end + (1|animal)  + 
                    (0+ DT_PS_end|animal)+ 
                    (1 | step), data = coyote_steps, 
                  map = list(theta = factor(c(1:2, NA))),
                  start = list(theta = c(rep(0, 2), log(1e3))),
                  family = poisson)


#summary(DTPS_m)
#capture.output(#summary(DTPS_m),file="DTPS.doc")

AIC1 <- bbmle::AICtab(nullmodel, NDVI_m, POP_m, BUP_m, DTLT_m, DTMT_m, DTHT_m, DTNT_m, DTPS_m)
capture.output(AIC1,file="AIC1_singleHabitatComponents.doc")
# capture.output(AIC1,file="AIC1_singleHabitatComponents_nocovid.doc")


####three main urban variables


NDVI.POP <- glmmTMB(used ~ -1 +  NDVI_end +POP_end + 
                      (1|animal) + 
                      (0+NDVI_end|animal) + 
                      (0+POP_end|animal) +
                      (1|step), data = coyote_steps, 
                    map = list(theta = factor(c(1:3, NA))),
                    start = list(theta = c(rep(0, 3), log(1e3))),
                    family = poisson)
#capture.output(#summary(habmodel4),file="habmodel4.doc")


NDVI.BUP <- glmmTMB(used ~ -1 + NDVI_end + BUP_end + 
                      (1|animal) + 
                      (0+NDVI_end|animal) + 
                      (0+BUP_end|animal) + 
                      (1|step), data = coyote_steps, 
                    map = list(theta = factor(c(1:3, NA))),
                    start = list(theta = c(rep(0, 3), log(1e3))),
                    family = poisson)
#capture.output(#summary(habmodel5),file="habmodel5.doc")

NDVI.LT <- glmmTMB(used ~ -1 +  NDVI_end +DT_LT_end + 
                     (1|animal) + 
                     (0+NDVI_end|animal) + 
                     (0+DT_LT_end|animal) +
                     (1|step), data = coyote_steps, 
                   map = list(theta = factor(c(1:3, NA))),
                   start = list(theta = c(rep(0, 3), log(1e3))),
                   family = poisson)
#capture.output(#summary(habmodel4),file="habmodel4.doc")


NDVI.MT <- glmmTMB(used ~ -1 + NDVI_end + DT_MT_end + 
                     (1|animal) + 
                     (0+NDVI_end|animal) + 
                     (0+DT_MT_end|animal) + 
                     (1|step), data = coyote_steps, 
                   map = list(theta = factor(c(1:3, NA))),
                   start = list(theta = c(rep(0, 3), log(1e3))),
                   family = poisson)
#capture.output(#summary(habmodel5),file="habmodel5.doc")

NDVI.HT <- glmmTMB(used ~ -1 + NDVI_end + DT_HT_end + 
                     (1|animal) + 
                     (0+NDVI_end|animal) + 
                     (0+DT_HT_end|animal) + 
                     (1|step), data = coyote_steps, 
                   map = list(theta = factor(c(1:3, NA))),
                   start = list(theta = c(rep(0, 3), log(1e3))),
                   family = poisson)
#capture.output(#summary(habmodel5),file="habmodel5.doc")


NDVI.NT <- glmmTMB(used ~ -1 + NDVI_end + DT_NT_end + 
                     (1|animal) + 
                     (0+NDVI_end|animal) + 
                     (0+DT_NT_end|animal) + 
                     (1|step), data = coyote_steps, 
                   map = list(theta = factor(c(1:3, NA))),
                   start = list(theta = c(rep(0, 3), log(1e3))),
                   family = poisson)
#capture.output(#summary(habmodel5),file="habmodel5.doc")

NDDVI.PS <- glmmTMB(used ~ -1 + NDVI_end + DT_PS_end + 
                      (1|animal) + 
                      (0+NDVI_end|animal) + 
                      (0+DT_PS_end|animal) + 
                      (1|step), data = coyote_steps, 
                    map = list(theta = factor(c(1:3, NA))),
                    start = list(theta = c(rep(0, 3), log(1e3))),
                    family = poisson)
AIV_env2 <- bbmle::AICtab(NDVI_m, NDVI.POP, NDVI.BUP, NDVI.LT, NDVI.MT, NDVI.HT, NDVI.NT, NDVI.PS)
capture.output(AIV_env2,file="habmodel_stage2.doc")



NDVI.BUP.POP <- glmmTMB(used ~ -1 + NDVI_end + BUP_end + POP_end  +
                          (1|animal) + 
                          (0+NDVI_end|animal) +
                          (0+BUP_end|animal) + 
                          (0+POP_end|animal) +
                          (1|step), data = coyote_steps, 
                        map = list(theta = factor(c(1:4, NA))),
                        start = list(theta = c(rep(0, 4), log(1e3))),
                        family = poisson)
#capture.output(#summary(habmodel1),file="habmodel1.doc")

NDVI.BUP.NT <- glmmTMB(used ~ -1 + NDVI_end + BUP_end + DT_NT_end  +
                         (1|animal) + 
                         (0+NDVI_end|animal) +
                         (0+BUP_end|animal) + 
                         (0+DT_NT_end|animal) +
                         (1|step), data = coyote_steps, 
                       map = list(theta = factor(c(1:4, NA))),
                       start = list(theta = c(rep(0, 4), log(1e3))),
                       family = poisson)
NDVI.BUP.LT <- glmmTMB(used ~ -1 + NDVI_end + BUP_end + DT_LT_end  +
                         (1|animal) + 
                         (0+NDVI_end|animal) +
                         (0+BUP_end|animal) + 
                         (0+DT_LT_end|animal) +
                         (1|step), data = coyote_steps, 
                       map = list(theta = factor(c(1:4, NA))),
                       start = list(theta = c(rep(0, 4), log(1e3))),
                       family = poisson)
NDVI.BUP.MT <- glmmTMB(used ~ -1 + NDVI_end + BUP_end + DT_MT_end  +
                         (1|animal) + 
                         (0+NDVI_end|animal) +
                         (0+BUP_end|animal) + 
                         (0+DT_MT_end|animal) +
                         (1|step), data = coyote_steps, 
                       map = list(theta = factor(c(1:4, NA))),
                       start = list(theta = c(rep(0, 4), log(1e3))),
                       family = poisson)
NDVI.BUP.HT <- glmmTMB(used ~ -1 + NDVI_end + BUP_end + DT_HT_end  +
                         (1|animal) + 
                         (0+NDVI_end|animal) +
                         (0+BUP_end|animal) + 
                         (0+DT_HT_end|animal) +
                         (1|step), data = coyote_steps, 
                       map = list(theta = factor(c(1:4, NA))),
                       start = list(theta = c(rep(0, 4), log(1e3))),
                       family = poisson)
NDVI.BUP.PS <- glmmTMB(used ~ -1 + NDVI_end + BUP_end + DT_PS_end  +
                         (1|animal) + 
                         (0+NDVI_end|animal) +
                         (0+BUP_end|animal) + 
                         (0+DT_PS_end|animal) +
                         (1|step), data = coyote_steps, 
                       map = list(theta = factor(c(1:4, NA))),
                       start = list(theta = c(rep(0, 4), log(1e3))),
                       family = poisson)

NDVI.BUP<-habmodel5
AIV_env3 <- bbmle::AICtab(NDVI.BUP, NDVI.BUP.POP, NDVI.BUP.LT, NDVI.BUP.MT, NDVI.BUP.HT, NDVI.BUP.NT, NDVI.BUP.PS)
capture.output(AIV_env3,file="habmodel_stage3.doc")

NDVI.BUP.PS.NT <- glmmTMB(used ~ -1 + NDVI_end + BUP_end + DT_PS_end  + DT_NT_end +
                            (1|animal) + 
                            (0+NDVI_end|animal) +
                            (0+BUP_end|animal) + 
                            (0+DT_PS_end|animal) +
                            (0+DT_NT_end|animal) + 
                            (1|step), data = coyote_steps, 
                          map = list(theta = factor(c(1:5, NA))),
                          start = list(theta = c(rep(0, 5), log(1e3))),
                          family = poisson)

NDVI.BUP.PS.POP <- glmmTMB(used ~ -1 + NDVI_end + BUP_end + DT_PS_end  + POP_end +
                             (1|animal) + 
                             (0+NDVI_end|animal) +
                             (0+BUP_end|animal) + 
                             (0+DT_PS_end|animal) +
                             (0+POP_end|animal) + 
                             (1|step), data = coyote_steps, 
                           map = list(theta = factor(c(1:5, NA))),
                           start = list(theta = c(rep(0, 5), log(1e3))),
                           family = poisson)
NDVI.BUP.PS.MT <- glmmTMB(used ~ -1 + NDVI_end + BUP_end + DT_PS_end  + DT_MT_end +
                            (1|animal) + 
                            (0+NDVI_end|animal) +
                            (0+BUP_end|animal) + 
                            (0+DT_PS_end|animal) +
                            (0+DT_MT_end|animal) + 
                            (1|step), data = coyote_steps, 
                          map = list(theta = factor(c(1:5, NA))),
                          start = list(theta = c(rep(0, 5), log(1e3))),
                          family = poisson)
NDVI.BUP.PS.LT <- glmmTMB(used ~ -1 + NDVI_end + BUP_end + DT_PS_end  + DT_LT_end +
                            (1|animal) + 
                            (0+NDVI_end|animal) +
                            (0+BUP_end|animal) + 
                            (0+DT_PS_end|animal) +
                            (0+DT_LT_end|animal) + 
                            (1|step), data = coyote_steps, 
                          map = list(theta = factor(c(1:5, NA))),
                          start = list(theta = c(rep(0, 5), log(1e3))),
                          family = poisson)
NDVI.BUP.PS.HT <- glmmTMB(used ~ -1 + NDVI_end + BUP_end + DT_PS_end  + DT_HT_end +
                            (1|animal) + 
                            (0+NDVI_end|animal) +
                            (0+BUP_end|animal) + 
                            (0+DT_PS_end|animal) +
                            (0+DT_HT_end|animal) + 
                            (1|step), data = coyote_steps, 
                          map = list(theta = factor(c(1:5, NA))),
                          start = list(theta = c(rep(0, 5), log(1e3))),
                          family = poisson)

AIV_env4 <- bbmle::AICtab(NDVI.BUP.PS, NDVI.BUP.PS.POP, NDVI.BUP.PS.LT, NDVI.BUP.PS.MT, NDVI.BUP.PS.HT, NDVI.BUP.PS.NT)
capture.output(AIV_env4,file="habmodel_stage4.doc")

NDVI.BUP.PS.NT.POP <- glmmTMB(used ~ -1  + NDVI_end + BUP_end + DT_PS_end  + DT_NT_end + POP_end + 
                                (1|animal) + 
                                (0+NDVI_end|animal) + 
                                (0+BUP_end|animal) +
                                (0+DT_PS_end|animal) +
                                (0+DT_NT_end|animal) +  
                                (0+POP_end|animal) +  
                                (1|step), data = coyote_steps, 
                              map = list(theta = factor(c(1:6, NA))),
                              start = list(theta = c(rep(0, 6), log(1e3))),
                              family = poisson)

NDVI.BUP.PS.NT.MT <- glmmTMB(used ~ -1  + NDVI_end + BUP_end + DT_PS_end  + DT_NT_end + DT_MT_end + 
                               (1|animal) + 
                               (0+NDVI_end|animal) + 
                               (0+BUP_end|animal) +
                               (0+DT_PS_end|animal) +
                               (0+DT_NT_end|animal) +  
                               (0+DT_MT_end|animal) +  
                               (1|step), data = coyote_steps, 
                             map = list(theta = factor(c(1:6, NA))),
                             start = list(theta = c(rep(0, 6), log(1e3))),
                             family = poisson)

NDVI.BUP.PS.NT.LT <- glmmTMB(used ~ -1  + NDVI_end + BUP_end + DT_PS_end  + DT_NT_end + DT_LT_end + 
                               (1|animal) + 
                               (0+NDVI_end|animal) + 
                               (0+BUP_end|animal) +
                               (0+DT_PS_end|animal) +
                               (0+DT_NT_end|animal) +  
                               (0+DT_LT_end|animal) +  
                               (1|step), data = coyote_steps, 
                             map = list(theta = factor(c(1:6, NA))),
                             start = list(theta = c(rep(0, 6), log(1e3))),
                             family = poisson)


NDVI.BUP.PS.NT.HT <- glmmTMB(used ~ -1  + NDVI_end + BUP_end + DT_PS_end  + DT_NT_end + DT_HT_end + 
                               (1|animal) + 
                               (0+NDVI_end|animal) + 
                               (0+BUP_end|animal) +
                               (0+DT_PS_end|animal) +
                               (0+DT_NT_end|animal) +  
                               (0+DT_HT_end|animal) +  
                               (1|step), data = coyote_steps, 
                             map = list(theta = factor(c(1:6, NA))),
                             start = list(theta = c(rep(0, 6), log(1e3))),
                             family = poisson)

AIV_env5 <- bbmle::AICtab(NDVI.BUP.PS.NT, NDVI.BUP.PS.NT.POP, NDVI.BUP.PS.NT.LT, NDVI.BUP.PS.NT.MT, NDVI.BUP.PS.NT.HT)
capture.output(AIV_env5,file="habmodel_stage5.doc")

NDVI.BUP.PS.NT.POP.LT <- glmmTMB(used ~ -1 + NDVI_end + BUP_end + DT_PS_end  + DT_NT_end + POP_end + DT_LT_end +
                                   (1|animal) + 
                                   (0+NDVI_end|animal) + 
                                   (0+BUP_end|animal) +
                                   (0+DT_PS_end|animal) +
                                   (0+DT_NT_end|animal) +  
                                   (0+POP_end|animal) +
                                   (0+DT_LT_end|animal) + 
                                   (1|step), data = coyote_steps, 
                                 map = list(theta = factor(c(1:7, NA))),
                                 start = list(theta = c(rep(0, 7), log(1e3))),
                                 family = poisson)

NDVI.BUP.PS.NT.POP.MT <- glmmTMB(used ~ -1 + NDVI_end + BUP_end + DT_PS_end  + DT_NT_end + POP_end + DT_MT_end +
                                   (1|animal) + 
                                   (0+NDVI_end|animal) + 
                                   (0+BUP_end|animal) +
                                   (0+DT_PS_end|animal) +
                                   (0+DT_NT_end|animal) +  
                                   (0+POP_end|animal) +
                                   (0+DT_MT_end|animal) + 
                                   (1|step), data = coyote_steps, 
                                 map = list(theta = factor(c(1:7, NA))),
                                 start = list(theta = c(rep(0, 7), log(1e3))),
                                 family = poisson)
NDVI.BUP.PS.NT.POP.HT <- glmmTMB(used ~ -1 + NDVI_end + BUP_end + DT_PS_end  + DT_NT_end + POP_end + DT_HT_end +
                                   (1|animal) + 
                                   (0+NDVI_end|animal) + 
                                   (0+BUP_end|animal) +
                                   (0+DT_PS_end|animal) +
                                   (0+DT_NT_end|animal) +  
                                   (0+POP_end|animal) +
                                   (0+DT_HT_end|animal) + 
                                   (1|step), data = coyote_steps, 
                                 map = list(theta = factor(c(1:7, NA))),
                                 start = list(theta = c(rep(0, 7), log(1e3))),
                                 family = poisson)

AIV_env6 <- bbmle::AICtab(NDVI.BUP.PS.NT.POP, NDVI.BUP.PS.NT.POP.LT, NDVI.BUP.PS.NT.POP.MT, NDVI.BUP.PS.NT.POP.HT)
capture.output(AIV_env6,file="habmodel_stage6.doc")

NDVI.BUP.PS.NT.POP.MT.LT <- glmmTMB(used ~ -1 + NDVI_end + BUP_end + DT_PS_end  + DT_NT_end + POP_end + DT_MT_end + DT_LT_end +
                                      (1|animal) + 
                                      (0+NDVI_end|animal) + 
                                      (0+BUP_end|animal) +
                                      (0+DT_PS_end|animal) +
                                      (0+DT_NT_end|animal) +  
                                      (0+POP_end|animal) +
                                      (0+DT_MT_end|animal) + 
                                      (0+DT_LT_end|animal) + 
                                      (1|step), data = coyote_steps, 
                                    map = list(theta = factor(c(1:8, NA))),
                                    start = list(theta = c(rep(0, 8), log(1e3))),
                                    family = poisson)

NDVI.BUP.PS.NT.POP.MT.HT <- glmmTMB(used ~ -1 + NDVI_end + BUP_end + DT_PS_end  + DT_NT_end + POP_end + DT_MT_end + DT_HT_end +
                                      (1|animal) + 
                                      (0+NDVI_end|animal) + 
                                      (0+BUP_end|animal) +
                                      (0+DT_PS_end|animal) +
                                      (0+DT_NT_end|animal) +  
                                      (0+POP_end|animal) +
                                      (0+DT_MT_end|animal) + 
                                      (0+DT_HT_end|animal) + 
                                      (1|step), data = coyote_steps, 
                                    map = list(theta = factor(c(1:8, NA))),
                                    start = list(theta = c(rep(0, 8), log(1e3))),
                                    family = poisson)

AIV_env7 <- bbmle::AICtab(NDVI.BUP.PS.NT.POP.MT, NDVI.BUP.PS.NT.POP.MT.LT, NDVI.BUP.PS.NT.POP.MT.HT)
capture.output(AIV_env7,file="habmodel_stage7.doc")

#####BEST FIT ENVIRONMENTAL HSF MODEL####
NDVI.BUP.PS.NT.POP.MT.LT.HT <- glmmTMB(used ~ -1  + NDVI_end + BUP_end + DT_PS_end  + DT_NT_end + POP_end + DT_MT_end + DT_LT_end + DT_HT_end +
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

AIV_env8 <- bbmle::AICtab(NDVI.BUP.PS.NT.POP.MT.LT, NDVI.BUP.PS.NT.POP.MT.LT.HT)
capture.output(AIV_env8,file="habmodel_stage8.doc")
#capture.output(AIV_env8,file="habmodel_stage8_nocovid.doc")

saveRDS(NDVI.BUP.PS.NT.POP.MT.LT.HT, "habmodel22_NDVI.BUP.PS.NT.POP.MT.LT.HT.rds")
capture.output(summary(NDVI.BUP.PS.NT.POP.MT.LT.HT),file="habmodel22_NDVI.BUP.PS.NT.POP.MT.LT.HT.doc")
# saveRDS(NDVI.BUP.PS.NT.POP.MT.LT.HT, "habmodel22_NDVI.BUP.PS.NT.POP.MT.LT.HT_nocovid.rds")
# capture.output(summary(NDVI.BUP.PS.NT.POP.MT.LT.HT),file="habmodel22_NDVI.BUP.PS.NT.POP.MT.LT.HT_nocovid.doc")


####day round 1 ######

daymodela <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                       DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                       NDVI_end:day +
                       (1|animal) + 
                       (0+NDVI_end|animal) + 
                       (0+POP_end|animal) +
                       (0+BUP_end|animal) + 
                       (0+DT_LT_end|animal) + 
                       (0+DT_MT_end|animal) + 
                       (0+DT_HT_end|animal) + 
                       (0+DT_NT_end|animal) + 
                       (0+DT_PS_end|animal) + 
                       (0+NDVI_end:day|animal) + 
                       (1|step), data = coyote_steps, 
                     map = list(theta = factor(c(1:10, NA))),
                     start = list(theta = c(rep(0, 10), log(1e3))),
                     family = poisson)

saveRDS(daymodela, "daymodela.rds")
capture.output(summary(daymodela),file="daymodela.doc")

daymodelb <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                       DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                       POP_end:day +
                       (1|animal) + 
                       (0+NDVI_end|animal) + 
                       (0+POP_end|animal) +
                       (0+BUP_end|animal) + 
                       (0+DT_LT_end|animal) + 
                       (0+DT_MT_end|animal) + 
                       (0+DT_HT_end|animal) + 
                       (0+DT_NT_end|animal) + 
                       (0+DT_PS_end|animal) + 
                       (0+POP_end:day|animal) + 
                       (1|step), data = coyote_steps, 
                     map = list(theta = factor(c(1:10, NA))),
                     start = list(theta = c(rep(0, 10), log(1e3))),
                     family = poisson)

saveRDS(daymodelb, "daymodelb.rds")
capture.output(summary(daymodelb),file="daymodelb.doc")

daymodelc <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                       DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                       BUP_end:day +
                       (1|animal) + 
                       (0+NDVI_end|animal) + 
                       (0+POP_end|animal) +
                       (0+BUP_end|animal) + 
                       (0+DT_LT_end|animal) + 
                       (0+DT_MT_end|animal) + 
                       (0+DT_HT_end|animal) + 
                       (0+DT_NT_end|animal) + 
                       (0+DT_PS_end|animal) + 
                       (0+BUP_end:day|animal) + 
                       (1|step), data = coyote_steps, 
                     map = list(theta = factor(c(1:10, NA))),
                     start = list(theta = c(rep(0, 10), log(1e3))),
                     family = poisson)

saveRDS(daymodelc, "daymodelc.rds")
capture.output(summary(daymodelc),file="daymodelc.doc")

daymodeld <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                       DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                       DT_LT_end:day +
                       (1|animal) + 
                       (0+NDVI_end|animal) + 
                       (0+POP_end|animal) +
                       (0+BUP_end|animal) + 
                       (0+DT_LT_end|animal) + 
                       (0+DT_MT_end|animal) + 
                       (0+DT_HT_end|animal) + 
                       (0+DT_NT_end|animal) + 
                       (0+DT_PS_end|animal) + 
                       (0+DT_LT_end:day|animal) + 
                       (1|step), data = coyote_steps, 
                     map = list(theta = factor(c(1:10, NA))),
                     start = list(theta = c(rep(0, 10), log(1e3))),
                     family = poisson)

saveRDS(daymodeld, "daymodeld.rds")
capture.output(summary(daymodeld),file="daymodeld.doc")


daymodele <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                       DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                       DT_MT_end:day +
                       (1|animal) + 
                       (0+NDVI_end|animal) + 
                       (0+POP_end|animal) +
                       (0+BUP_end|animal) + 
                       (0+DT_LT_end|animal) + 
                       (0+DT_MT_end|animal) + 
                       (0+DT_HT_end|animal) + 
                       (0+DT_NT_end|animal) + 
                       (0+DT_PS_end|animal) + 
                       (0+DT_MT_end:day|animal) + 
                       (1|step), data = coyote_steps, 
                     map = list(theta = factor(c(1:10, NA))),
                     start = list(theta = c(rep(0, 10), log(1e3))),
                     family = poisson)
saveRDS(daymodele, "daymodele.rds")
capture.output(summary(daymodele),file="daymodele.doc")

daymodelf <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                       DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                       DT_HT_end:day +
                       (1|animal) + 
                       (0+NDVI_end|animal) + 
                       (0+POP_end|animal) +
                       (0+BUP_end|animal) + 
                       (0+DT_LT_end|animal) + 
                       (0+DT_MT_end|animal) + 
                       (0+DT_HT_end|animal) + 
                       (0+DT_NT_end|animal) + 
                       (0+DT_PS_end|animal) + 
                       (0+DT_HT_end:day|animal) + 
                       (1|step), data = coyote_steps, 
                     map = list(theta = factor(c(1:10, NA))),
                     start = list(theta = c(rep(0, 10), log(1e3))),
                     family = poisson)
saveRDS(daymodelf, "daymodelf.rds")
capture.output(summary(daymodelf),file="daymodelf.doc")

daymodelg <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                       DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                       DT_NT_end:day +
                       (1|animal) + 
                       (0+NDVI_end|animal) + 
                       (0+POP_end|animal) +
                       (0+BUP_end|animal) + 
                       (0+DT_LT_end|animal) + 
                       (0+DT_MT_end|animal) + 
                       (0+DT_HT_end|animal) + 
                       (0+DT_NT_end|animal) + 
                       (0+DT_PS_end|animal) + 
                       (0+DT_NT_end:day|animal) + 
                       (1|step), data = coyote_steps, 
                     map = list(theta = factor(c(1:10, NA))),
                     start = list(theta = c(rep(0, 10), log(1e3))),
                     family = poisson)

saveRDS(daymodelg, "daymodelg.rds")
capture.output(summary(daymodelg),file="daymodelg.doc")

daymodelh <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                       DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                       DT_PS_end:day +
                       (1|animal) + 
                       (0+NDVI_end|animal) + 
                       (0+POP_end|animal) +
                       (0+BUP_end|animal) + 
                       (0+DT_LT_end|animal) + 
                       (0+DT_MT_end|animal) + 
                       (0+DT_HT_end|animal) + 
                       (0+DT_NT_end|animal) + 
                       (0+DT_PS_end|animal) + 
                       (0+DT_PS_end:day|animal) + 
                       (1|step), data = coyote_steps, 
                     map = list(theta = factor(c(1:10, NA))),
                     start = list(theta = c(rep(0, 10), log(1e3))),
                     family = poisson)
saveRDS(daymodelh, "daymodelh.rds")
capture.output(summary(daymodelh),file="daymodelh.doc")
AICday1 <- bbmle::AICtab(habmodel22, daymodela, daymodelb, daymodelc, daymodeld, daymodele, daymodelf, daymodelg, daymodelh)
capture.output(AICday1,file="AIC_daymodel_stage1.doc")

###daymodel round 2####

daymodel10<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                      DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + DT_LT_end:day+NDVI_end:day+(1|animal) + 
                      (0+NDVI_end|animal) + 
                      (0+POP_end|animal) +
                      (0+BUP_end|animal) + 
                      (0+DT_LT_end|animal) + 
                      (0+DT_MT_end|animal) + 
                      (0+DT_HT_end|animal) + 
                      (0+DT_NT_end|animal) + 
                      (0+DT_PS_end|animal) +(0+DT_LT_end:day|animal) + (0+NDVI_end:day|animal) + (1|step), data = coyote_steps, 
                    map = list(theta = factor(c(1:11, NA))),
                    start = list(theta = c(rep(0, 11), log(1e3))),
                    family = poisson)
saveRDS(daymodel10, "daymodel10.rds")
capture.output(summary(daymodel10),file="daymodel10.doc")
daymodel11<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                      DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + DT_LT_end:day+BUP_end:day +(1|animal) + 
                      (0+NDVI_end|animal) + 
                      (0+POP_end|animal) +
                      (0+BUP_end|animal) + 
                      (0+DT_LT_end|animal) + 
                      (0+DT_MT_end|animal) + 
                      (0+DT_HT_end|animal) + 
                      (0+DT_NT_end|animal) + 
                      (0+DT_PS_end|animal) +(0+DT_LT_end:day|animal) + (0+BUP_end:day |animal) + (1|step), data = coyote_steps, 
                    map = list(theta = factor(c(1:11, NA))),
                    start = list(theta = c(rep(0, 11), log(1e3))),
                    family = poisson)
saveRDS(daymodel11, "daymodel11.rds")
capture.output(summary(daymodel11),file="daymodel11.doc")
daymodel12<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                      DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + DT_LT_end:day+POP_end:day +(1|animal) + 
                      (0+NDVI_end|animal) + 
                      (0+POP_end|animal) +
                      (0+BUP_end|animal) + 
                      (0+DT_LT_end|animal) + 
                      (0+DT_MT_end|animal) + 
                      (0+DT_HT_end|animal) + 
                      (0+DT_NT_end|animal) + 
                      (0+DT_PS_end|animal) +(0+DT_LT_end:day|animal) + (0+POP_end:day |animal) + (1|step), data = coyote_steps, 
                    map = list(theta = factor(c(1:11, NA))),
                    start = list(theta = c(rep(0, 11), log(1e3))),
                    family = poisson)
saveRDS(daymodel12, "daymodel12.rds")
capture.output(summary(daymodel12),file="daymodel12.doc")
daymodel13<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                      DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + DT_LT_end:day+DT_MT_end:day +(1|animal) + 
                      (0+NDVI_end|animal) + 
                      (0+POP_end|animal) +
                      (0+BUP_end|animal) + 
                      (0+DT_LT_end|animal) + 
                      (0+DT_MT_end|animal) + 
                      (0+DT_HT_end|animal) + 
                      (0+DT_NT_end|animal) + 
                      (0+DT_PS_end|animal) +(0+DT_LT_end:day|animal) + (0+DT_MT_end:day |animal) + (1|step), data = coyote_steps, 
                    map = list(theta = factor(c(1:11, NA))),
                    start = list(theta = c(rep(0, 11), log(1e3))),
                    family = poisson)
saveRDS(daymodel13, "daymodel13.rds")
capture.output(summary(daymodel13),file="daymodel13.doc")
daymodel14<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                      DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + DT_LT_end:day+DT_HT_end:day +(1|animal) + 
                      (0+NDVI_end|animal) + 
                      (0+POP_end|animal) +
                      (0+BUP_end|animal) + 
                      (0+DT_LT_end|animal) + 
                      (0+DT_MT_end|animal) + 
                      (0+DT_HT_end|animal) + 
                      (0+DT_NT_end|animal) + 
                      (0+DT_PS_end|animal) +(0+DT_LT_end:day|animal) + (0+DT_HT_end:day |animal) + (1|step), data = coyote_steps, 
                    map = list(theta = factor(c(1:11, NA))),
                    start = list(theta = c(rep(0, 11), log(1e3))),
                    family = poisson)
saveRDS(daymodel14, "daymodel14.rds")
capture.output(summary(daymodel14),file="daymodel14.doc")
daymodel15<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                      DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + DT_LT_end:day+DT_NT_end:day +(1|animal) + 
                      (0+NDVI_end|animal) + 
                      (0+POP_end|animal) +
                      (0+BUP_end|animal) + 
                      (0+DT_LT_end|animal) + 
                      (0+DT_MT_end|animal) + 
                      (0+DT_HT_end|animal) + 
                      (0+DT_NT_end|animal) + 
                      (0+DT_PS_end|animal) +(0+DT_LT_end:day|animal) + (0+DT_NT_end:day |animal) + (1|step), data = coyote_steps, 
                    map = list(theta = factor(c(1:11, NA))),
                    start = list(theta = c(rep(0, 11), log(1e3))),
                    family = poisson)
saveRDS(daymodel15, "daymodel15.rds")
capture.output(summary(daymodel15),file="daymodel15.doc")
daymodel16<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                      DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + DT_LT_end:day+DT_PS_end:day +(1|animal) + 
                      (0+NDVI_end|animal) + 
                      (0+POP_end|animal) +
                      (0+BUP_end|animal) + 
                      (0+DT_LT_end|animal) + 
                      (0+DT_MT_end|animal) + 
                      (0+DT_HT_end|animal) + 
                      (0+DT_NT_end|animal) + 
                      (0+DT_PS_end|animal) +(0+DT_LT_end:day|animal) + (0+DT_PS_end:day |animal) + (1|step), data = coyote_steps, 
                    map = list(theta = factor(c(1:11, NA))),
                    start = list(theta = c(rep(0, 11), log(1e3))),
                    family = poisson)
saveRDS(daymodel16, "daymodel16.rds")
capture.output(summary(daymodel16),file="daymodel16.doc")
day2 <- bbmle::AICtab(daymodeld, daymodel10, daymodel11, daymodel12, daymodel13,daymodel14, daymodel15, daymodel16)
capture.output(day2,file="AIC_daymodel_stage2.doc")

###day round 3####
daymodel17<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                      DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + DT_LT_end:day+NDVI_end:day +DT_PS_end:day+(1|animal) + 
                      (0+NDVI_end|animal) + 
                      (0+POP_end|animal) +
                      (0+BUP_end|animal) + 
                      (0+DT_LT_end|animal) + 
                      (0+DT_MT_end|animal) + 
                      (0+DT_HT_end|animal) + 
                      (0+DT_NT_end|animal) + 
                      (0+DT_PS_end|animal) +(0+DT_LT_end:day|animal) + (0+NDVI_end:day |animal) + (0+DT_PS_end:day|animal) + (1|step), data = coyote_steps, 
                    map = list(theta = factor(c(1:12, NA))),
                    start = list(theta = c(rep(0, 12), log(1e3))),
                    family = poisson)
saveRDS(daymodel17, "daymodel17.rds")
capture.output(summary(daymodel17),file="daymodel17.doc")
daymodel18<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                      DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + DT_LT_end:day+NDVI_end:day +BUP_end:day +(1|animal) + 
                      (0+NDVI_end|animal) + 
                      (0+POP_end|animal) +
                      (0+BUP_end|animal) + 
                      (0+DT_LT_end|animal) + 
                      (0+DT_MT_end|animal) + 
                      (0+DT_HT_end|animal) + 
                      (0+DT_NT_end|animal) + 
                      (0+DT_PS_end|animal) +(0+DT_LT_end:day|animal) + (0+NDVI_end:day |animal) + (0+BUP_end:day |animal) + (1|step), data = coyote_steps, 
                    map = list(theta = factor(c(1:12, NA))),
                    start = list(theta = c(rep(0, 12), log(1e3))),
                    family = poisson)
saveRDS(daymodel18, "daymodel18.rds")
capture.output(summary(daymodel18),file="daymodel18.doc")
daymodel19<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                      DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + DT_LT_end:day+NDVI_end:day +DT_MT_end:day +(1|animal) + 
                      (0+NDVI_end|animal) + 
                      (0+POP_end|animal) +
                      (0+BUP_end|animal) + 
                      (0+DT_LT_end|animal) + 
                      (0+DT_MT_end|animal) + 
                      (0+DT_HT_end|animal) + 
                      (0+DT_NT_end|animal) + 
                      (0+DT_PS_end|animal) +(0+DT_LT_end:day|animal) + (0+NDVI_end:day |animal) + (0+DT_MT_end:day |animal) + (1|step), data = coyote_steps, 
                    map = list(theta = factor(c(1:12, NA))),
                    start = list(theta = c(rep(0, 12), log(1e3))),
                    family = poisson)
saveRDS(daymodel19, "daymodel19.rds")
capture.output(summary(daymodel19),file="daymodel19.doc")
daymodel20<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                      DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + DT_LT_end:day+NDVI_end:day +DT_HT_end:day +(1|animal) + 
                      (0+NDVI_end|animal) + 
                      (0+POP_end|animal) +
                      (0+BUP_end|animal) + 
                      (0+DT_LT_end|animal) + 
                      (0+DT_MT_end|animal) + 
                      (0+DT_HT_end|animal) + 
                      (0+DT_NT_end|animal) + 
                      (0+DT_PS_end|animal) +(0+DT_LT_end:day|animal) + (0+NDVI_end:day |animal) + (0+DT_HT_end:day |animal) + (1|step), data = coyote_steps, 
                    map = list(theta = factor(c(1:12, NA))),
                    start = list(theta = c(rep(0, 12), log(1e3))),
                    family = poisson)
saveRDS(daymodel20, "daymodel20.rds")
capture.output(summary(daymodel20),file="daymodel20.doc")
daymodel21<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                      DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + DT_LT_end:day+NDVI_end:day +DT_NT_end:day +(1|animal) + 
                      (0+NDVI_end|animal) + 
                      (0+POP_end|animal) +
                      (0+BUP_end|animal) + 
                      (0+DT_LT_end|animal) + 
                      (0+DT_MT_end|animal) + 
                      (0+DT_HT_end|animal) + 
                      (0+DT_NT_end|animal) + 
                      (0+DT_PS_end|animal) +(0+DT_LT_end:day|animal) + (0+NDVI_end:day |animal) + (0+DT_NT_end:day |animal) + (1|step), data = coyote_steps, 
                    map = list(theta = factor(c(1:12, NA))),
                    start = list(theta = c(rep(0, 12), log(1e3))),
                    family = poisson)
saveRDS(daymodel21, "daymodel21.rds")
capture.output(summary(daymodel21),file="daymodel21.doc")
daymodel22<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                      DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + DT_LT_end:day+NDVI_end:day +POP_end:day +(1|animal) + 
                      (0+NDVI_end|animal) + 
                      (0+POP_end|animal) +
                      (0+BUP_end|animal) + 
                      (0+DT_LT_end|animal) + 
                      (0+DT_MT_end|animal) + 
                      (0+DT_HT_end|animal) + 
                      (0+DT_NT_end|animal) + 
                      (0+DT_PS_end|animal) +(0+DT_LT_end:day|animal) + (0+NDVI_end:day |animal) + (0+POP_end:day |animal) + (1|step), data = coyote_steps, 
                    map = list(theta = factor(c(1:12, NA))),
                    start = list(theta = c(rep(0, 12), log(1e3))),
                    family = poisson)
saveRDS(daymodel22, "daymodel22.rds")
capture.output(summary(daymodel22),file="daymodel22.doc")
day3 <- bbmle::AICtab(daymodel10, daymodel17, daymodel18, daymodel19, daymodel20,daymodel21, daymodel22)
capture.output(day3,file="AIC_daymodel_stage3.doc")


###day model round 4####
daymodel23<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                      DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + DT_LT_end:day+NDVI_end:day +DT_NT_end:day +BUP_end:day +(1|animal) + 
                      (0+NDVI_end|animal) + 
                      (0+POP_end|animal) +
                      (0+BUP_end|animal) + 
                      (0+DT_LT_end|animal) + 
                      (0+DT_MT_end|animal) + 
                      (0+DT_HT_end|animal) + 
                      (0+DT_NT_end|animal) + 
                      (0+DT_PS_end|animal) +(0+DT_LT_end:day|animal) + (0+NDVI_end:day |animal) + (0+DT_NT_end:day |animal) + (0+BUP_end:day |animal) + (1|step), data = coyote_steps, 
                    map = list(theta = factor(c(1:13, NA))),
                    start = list(theta = c(rep(0, 13), log(1e3))),
                    family = poisson)
saveRDS(daymodel23, "daymodel23.rds")
capture.output(summary(daymodel23),file="daymodel23.doc")
daymodel24<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                      DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + DT_LT_end:day+NDVI_end:day +DT_NT_end:day +DT_MT_end:day +(1|animal) + 
                      (0+NDVI_end|animal) + 
                      (0+POP_end|animal) +
                      (0+BUP_end|animal) + 
                      (0+DT_LT_end|animal) + 
                      (0+DT_MT_end|animal) + 
                      (0+DT_HT_end|animal) + 
                      (0+DT_NT_end|animal) + 
                      (0+DT_PS_end|animal) +(0+DT_LT_end:day|animal) + (0+NDVI_end:day |animal) + (0+DT_NT_end:day |animal) + (0+DT_MT_end:day |animal) + (1|step), data = coyote_steps, 
                    map = list(theta = factor(c(1:13, NA))),
                    start = list(theta = c(rep(0, 13), log(1e3))),
                    family = poisson)
saveRDS(daymodel24, "daymodel24.rds")
capture.output(summary(daymodel24),file="daymodel24.doc")
daymodel25<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                      DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + DT_LT_end:day+NDVI_end:day +DT_NT_end:day +DT_HT_end:day +(1|animal) + 
                      (0+NDVI_end|animal) + 
                      (0+POP_end|animal) +
                      (0+BUP_end|animal) + 
                      (0+DT_LT_end|animal) + 
                      (0+DT_MT_end|animal) + 
                      (0+DT_HT_end|animal) + 
                      (0+DT_NT_end|animal) + 
                      (0+DT_PS_end|animal) +(0+DT_LT_end:day|animal) + (0+NDVI_end:day |animal) + (0+DT_NT_end:day |animal) + (0+DT_HT_end:day |animal) + (1|step), data = coyote_steps, 
                    map = list(theta = factor(c(1:13, NA))),
                    start = list(theta = c(rep(0, 13), log(1e3))),
                    family = poisson)
saveRDS(daymodel25, "daymodel25.rds")
capture.output(summary(daymodel25),file="daymodel25.doc")
daymodel26<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                      DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + DT_LT_end:day+NDVI_end:day +DT_NT_end:day +DT_PS_end:day +(1|animal) + 
                      (0+NDVI_end|animal) + 
                      (0+POP_end|animal) +
                      (0+BUP_end|animal) + 
                      (0+DT_LT_end|animal) + 
                      (0+DT_MT_end|animal) + 
                      (0+DT_HT_end|animal) + 
                      (0+DT_NT_end|animal) + 
                      (0+DT_PS_end|animal) +(0+DT_LT_end:day|animal) + (0+NDVI_end:day |animal) + (0+DT_NT_end:day |animal) + (0+DT_PS_end:day |animal) + (1|step), data = coyote_steps, 
                    map = list(theta = factor(c(1:13, NA))),
                    start = list(theta = c(rep(0, 13), log(1e3))),
                    family = poisson)
saveRDS(daymodel26, "daymodel26.rds")
capture.output(summary(daymodel26),file="daymodel26.doc")

daymodel27<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                      DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + DT_LT_end:day+NDVI_end:day +DT_NT_end:day +POP_end:day+(1|animal) + 
                      (0+NDVI_end|animal) + 
                      (0+POP_end|animal) +
                      (0+BUP_end|animal) + 
                      (0+DT_LT_end|animal) + 
                      (0+DT_MT_end|animal) + 
                      (0+DT_HT_end|animal) + 
                      (0+DT_NT_end|animal) + 
                      (0+DT_PS_end|animal) +(0+DT_LT_end:day|animal) + (0+NDVI_end:day |animal) + (0+DT_NT_end:day |animal) + (0+POP_end:day|animal) + (1|step), data = coyote_steps, 
                    map = list(theta = factor(c(1:13, NA))),
                    start = list(theta = c(rep(0, 13), log(1e3))),
                    family = poisson)
saveRDS(daymodel27, "daymodel27.rds")
capture.output(summary(daymodel27),file="daymodel27.doc")


day4 <- bbmle::AICtab(daymodel21, daymodel23, daymodel24, daymodel25, daymodel26,daymodel27)
capture.output(day4,file="AIC_daymodel_stage4.doc")


###day round 5 ####
daymodel28<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                      DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + DT_LT_end:day+NDVI_end:day +DT_NT_end:day +BUP_end:day +DT_MT_end:day +(1|animal) + 
                      (0+NDVI_end|animal) + 
                      (0+POP_end|animal) +
                      (0+BUP_end|animal) + 
                      (0+DT_LT_end|animal) + 
                      (0+DT_MT_end|animal) + 
                      (0+DT_HT_end|animal) + 
                      (0+DT_NT_end|animal) + 
                      (0+DT_PS_end|animal) +(0+DT_LT_end:day|animal) + (0+NDVI_end:day |animal) + (0+DT_NT_end:day |animal)  + (0+BUP_end:day |animal) + (0+DT_MT_end:day|animal) + (1|step), data = coyote_steps, 
                    map = list(theta = factor(c(1:14, NA))),
                    start = list(theta = c(rep(0, 14), log(1e3))),
                    family = poisson)
saveRDS(daymodel28, "daymodel28.rds")
capture.output(summary(daymodel28),file="daymodel28.doc")
daymodel29<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + 
                      DT_LT_end:day+ NDVI_end:day +DT_NT_end:day +BUP_end:day + POP_end:day +(1|animal) + 
                      (0+NDVI_end|animal) + 
                      (0+POP_end|animal) +
                      (0+BUP_end|animal) + 
                      (0+DT_LT_end|animal) + 
                      (0+DT_MT_end|animal) + 
                      (0+DT_HT_end|animal) + 
                      (0+DT_NT_end|animal) + 
                      (0+DT_PS_end|animal) +
                      (0+DT_LT_end:day|animal) + (0+NDVI_end:day |animal) + (0+DT_NT_end:day |animal) + 
                      (0+BUP_end:day |animal) + (0+POP_end:day |animal) + (1|step), data = coyote_steps, 
                    map = list(theta = factor(c(1:14, NA))),
                    start = list(theta = c(rep(0, 14), log(1e3))),
                    family = poisson)
saveRDS(daymodel29, "daymodel29.rds")
capture.output(summary(daymodel29),file="daymodel29.doc")
daymodel30<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + 
                      DT_LT_end:day + NDVI_end:day +DT_NT_end:day +BUP_end:day +DT_HT_end:day +
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
                      (0+BUP_end:day|animal) + (0+DT_HT_end:day |animal) + (1|step), data = coyote_steps, 
                    map = list(theta = factor(c(1:14, NA))),
                    start = list(theta = c(rep(0, 14), log(1e3))),
                    family = poisson)
saveRDS(daymodel30, "daymodel30.rds")
capture.output(summary(daymodel30),file="daymodel30.doc")
daymodel31<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                      DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + 
                      DT_LT_end:day + NDVI_end:day + DT_NT_end:day + BUP_end:day + DT_PS_end:day +(1|animal) + 
                      (0+NDVI_end|animal) + 
                      (0+POP_end|animal) +
                      (0+BUP_end|animal) + 
                      (0+DT_LT_end|animal) + 
                      (0+DT_MT_end|animal) + 
                      (0+DT_HT_end|animal) + 
                      (0+DT_NT_end|animal) + 
                      (0+DT_PS_end|animal) +
                      (0+DT_LT_end:day|animal) + (0+NDVI_end:day |animal) + (0+DT_NT_end:day |animal) + 
                      (0+BUP_end:day|animal) + (0+DT_PS_end:day |animal) + (1|step), data = coyote_steps, 
                    map = list(theta = factor(c(1:14, NA))),
                    start = list(theta = c(rep(0, 14), log(1e3))),
                    family = poisson)
saveRDS(daymodel31, "daymodel31.rds")
capture.output(summary(daymodel31),file="daymodel31.doc")
day5 <- bbmle::AICtab(daymodel23, daymodel28, daymodel29, daymodel30, daymodel31)
capture.output(day5,file="AIC_daymodel_stage5.doc")

###day round 6 ####
daymodel32<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + 
                      DT_LT_end:day + 
                      NDVI_end:day + DT_NT_end:day + BUP_end:day + DT_PS_end:day + POP_end:day +
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
                      (0+BUP_end:day|animal) + (0+DT_PS_end:day |animal) + (0+POP_end:day |animal) + (1|step), data = coyote_steps, 
                    map = list(theta = factor(c(1:15, NA))),
                    start = list(theta = c(rep(0, 15), log(1e3))),
                    family = poisson)
saveRDS(daymodel32, "daymodel32.rds")
capture.output(summary(daymodel32),file="daymodel32.doc")
daymodel33<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                      DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + 
                      DT_LT_end:day + NDVI_end:day + DT_NT_end:day + BUP_end:day + DT_PS_end:day + DT_MT_end:day +
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
                      (0+BUP_end:day|animal) + (0+DT_PS_end:day |animal) + (0+DT_MT_end:day |animal) + (1|step), data = coyote_steps,  
                    map = list(theta = factor(c(1:15, NA))),
                    start = list(theta = c(rep(0, 15), log(1e3))),
                    family = poisson)
saveRDS(daymodel33, "daymodel33.rds")
capture.output(summary(daymodel33),file="daymodel33.doc")
daymodel34<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                      DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + 
                      DT_LT_end:day + NDVI_end:day + DT_NT_end:day + BUP_end:day + DT_PS_end:day + DT_HT_end:day +
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
                      (0+BUP_end:day|animal) + (0+DT_PS_end:day |animal) + (0+DT_HT_end:day |animal) + (1|step), data = coyote_steps,  
                    map = list(theta = factor(c(1:15, NA))),
                    start = list(theta = c(rep(0, 15), log(1e3))),
                    family = poisson)
saveRDS(daymodel34, "daymodel34.rds")
capture.output(summary(daymodel34),file="daymodel34.doc")

day6 <- bbmle::AICtab(daymodel31, daymodel32, daymodel33, daymodel34)
capture.output(day6,file="AIC_daymodel_stage6.doc")

###day round 7 ####
daymodel35<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                      DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + 
                      DT_LT_end:day + NDVI_end:day + DT_NT_end:day + BUP_end:day + DT_PS_end:day + DT_HT_end:day +  DT_MT_end:day +
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
                      (0+BUP_end:day|animal) + (0+DT_PS_end:day |animal) + (0+DT_HT_end:day |animal) + (0+DT_MT_end:day |animal) + (1|step), data = coyote_steps,  
                    map = list(theta = factor(c(1:16, NA))),
                    start = list(theta = c(rep(0, 16), log(1e3))),
                    family = poisson)
saveRDS(daymodel35, "daymodel35.rds")
capture.output(summary(daymodel35),file="daymodel35.doc")

daymodel36<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                      DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + 
                      DT_LT_end:day + NDVI_end:day + DT_NT_end:day + BUP_end:day + DT_PS_end:day + DT_HT_end:day +  POP_end:day +
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
                      (0+BUP_end:day|animal) + (0+DT_PS_end:day |animal) + (0+DT_HT_end:day |animal) + (0+POP_end:day |animal) + (1|step), data = coyote_steps,  
                    map = list(theta = factor(c(1:16, NA))),
                    start = list(theta = c(rep(0, 16), log(1e3))),
                    family = poisson)
saveRDS(daymodel36, "daymodel36.rds")
capture.output(summary(daymodel36),file="daymodel36.doc")

day7 <- bbmle::AICtab(daymodel34, daymodel35, daymodel36)
capture.output(day7,file="AIC_daymodel_stage7.doc")

####day round 8#####
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
saveRDS(daymodel37, "daymodel37.rds") #previous name daymodel37.rds
capture.output(summary(daymodel37),file="daymodel37.doc")
# capture.output(summary(daymodel37),file="daymodel37_nocovid.doc")

daymodel34 <- readRDS("daymodel34.rds") 
daymodel37 <- readRDS("daymodel37.rds")

day8 <- bbmle::AICtab(daymodel34, daymodel37)
capture.output(day8,file="AIC_daymodel_stage8.doc")
#BEST FIT: daymodel37(abcdefgh)

####transient round 1 ####


transientmodela <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                             DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                             NDVI_end:transient +
                             (1|animal) + 
                             (0+NDVI_end|animal) + 
                             (0+POP_end|animal) +
                             (0+BUP_end|animal) + 
                             (0+DT_LT_end|animal) + 
                             (0+DT_MT_end|animal) + 
                             (0+DT_HT_end|animal) + 
                             (0+DT_NT_end|animal) + 
                             (0+DT_PS_end|animal) + 
                             (0+NDVI_end:transient|animal) + 
                             (1|step), data = coyote_steps, 
                           map = list(theta = factor(c(1:10, NA))),
                           start = list(theta = c(rep(0, 10), log(1e3))),
                           family = poisson)
saveRDS(transientmodela, "transientmodela.rds")
capture.output(summary(transientmodela),file="transientmodela.doc")


transientmodelb <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                             DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                             POP_end:transient +
                             (1|animal) + 
                             (0+NDVI_end|animal) + 
                             (0+POP_end|animal) +
                             (0+BUP_end|animal) + 
                             (0+DT_LT_end|animal) + 
                             (0+DT_MT_end|animal) + 
                             (0+DT_HT_end|animal) + 
                             (0+DT_NT_end|animal) + 
                             (0+DT_PS_end|animal) + 
                             (0+POP_end:transient|animal) + 
                             (1|step), data = coyote_steps, 
                           map = list(theta = factor(c(1:10, NA))),
                           start = list(theta = c(rep(0, 10), log(1e3))),
                           family = poisson)
saveRDS(transientmodelb, "transientmodelb.rds")
capture.output(summary(transientmodelb),file="transientmodelb.doc")

transientmodelc <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                             DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                             BUP_end:transient +
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
                             (1|step), data = coyote_steps, 
                           map = list(theta = factor(c(1:10, NA))),
                           start = list(theta = c(rep(0, 10), log(1e3))),
                           family = poisson)

saveRDS(transientmodelc, "transientmodelc.rds")
capture.output(summary(transientmodelc),file="transientmodelc.doc")

transientmodeld <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                             DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                             DT_LT_end:transient +
                             (1|animal) + 
                             (0+NDVI_end|animal) + 
                             (0+POP_end|animal) +
                             (0+BUP_end|animal) + 
                             (0+DT_LT_end|animal) + 
                             (0+DT_MT_end|animal) + 
                             (0+DT_HT_end|animal) + 
                             (0+DT_NT_end|animal) + 
                             (0+DT_PS_end|animal) + 
                             (0+DT_LT_end:transient|animal) + 
                             (1|step), data = coyote_steps, 
                           map = list(theta = factor(c(1:10, NA))),
                           start = list(theta = c(rep(0, 10), log(1e3))),
                           family = poisson)
saveRDS(transientmodeld, "transientmodeld.rds")
capture.output(summary(transientmodeld),file="transientmodeld.doc")

transientmodele <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                             DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                             DT_MT_end:transient +
                             (1|animal) + 
                             (0+NDVI_end|animal) + 
                             (0+POP_end|animal) +
                             (0+BUP_end|animal) + 
                             (0+DT_LT_end|animal) + 
                             (0+DT_MT_end|animal) + 
                             (0+DT_HT_end|animal) + 
                             (0+DT_NT_end|animal) + 
                             (0+DT_PS_end|animal) + 
                             (0+DT_MT_end:transient|animal) + 
                             (1|step), data = coyote_steps, 
                           map = list(theta = factor(c(1:10, NA))),
                           start = list(theta = c(rep(0, 10), log(1e3))),
                           family = poisson)
saveRDS(transientmodele, "transientmodele.rds")
capture.output(summary(transientmodele),file="transientmodele.doc")

transientmodelf <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                             DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                             DT_HT_end:transient +
                             (1|animal) + 
                             (0+NDVI_end|animal) + 
                             (0+POP_end|animal) +
                             (0+BUP_end|animal) + 
                             (0+DT_LT_end|animal) + 
                             (0+DT_MT_end|animal) + 
                             (0+DT_HT_end|animal) + 
                             (0+DT_NT_end|animal) + 
                             (0+DT_PS_end|animal) + 
                             (0+DT_HT_end:transient|animal) + 
                             (1|step), data = coyote_steps, 
                           map = list(theta = factor(c(1:10, NA))),
                           start = list(theta = c(rep(0, 10), log(1e3))),
                           family = poisson)
saveRDS(transientmodelf, "transientmodelf.rds")
capture.output(summary(transientmodelf),file="transientmodelf.doc")

transientmodelg <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                             DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                             DT_NT_end:transient +
                             (1|animal) + 
                             (0+NDVI_end|animal) + 
                             (0+POP_end|animal) +
                             (0+BUP_end|animal) + 
                             (0+DT_LT_end|animal) + 
                             (0+DT_MT_end|animal) + 
                             (0+DT_HT_end|animal) + 
                             (0+DT_NT_end|animal) + 
                             (0+DT_PS_end|animal) + 
                             (0+DT_NT_end:transient|animal) + 
                             (1|step), data = coyote_steps, 
                           map = list(theta = factor(c(1:10, NA))),
                           start = list(theta = c(rep(0, 10), log(1e3))),
                           family = poisson)

saveRDS(transientmodelg, "transientmodelg.rds")
capture.output(summary(transientmodelg),file="transientmodelg.doc")

transientmodelh <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                             DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                             DT_PS_end:transient +
                             (1|animal) + 
                             (0+NDVI_end|animal) + 
                             (0+POP_end|animal) +
                             (0+BUP_end|animal) + 
                             (0+DT_LT_end|animal) + 
                             (0+DT_MT_end|animal) + 
                             (0+DT_HT_end|animal) + 
                             (0+DT_NT_end|animal) + 
                             (0+DT_PS_end|animal) + 
                             (0+DT_PS_end:transient|animal) + 
                             (1|step), data = coyote_steps, 
                           map = list(theta = factor(c(1:10, NA))),
                           start = list(theta = c(rep(0, 10), log(1e3))),
                           family = poisson)
saveRDS(transientmodelh, "transientmodelh.rds")
capture.output(summary(transientmodelh),file="transientmodelh.doc")

habmodel22 <- readRDS("habmodel22_NDVI.BUP.PS.NT.POP.MT.LT.HT.rds")
transientmodela <- readRDS("transientmodela.rds") 
transientmodelb <- readRDS("transientmodelb.rds") 
transientmodelc <- readRDS("transientmodelc.rds") 
transientmodeld <- readRDS("transientmodeld.rds") 
transientmodele <- readRDS("transientmodele.rds") 
transientmodelf <- readRDS("transientmodelf.rds") 
transientmodelg <- readRDS("transientmodelg.rds") 
transientmodelh <- readRDS("transientmodelh.rds")

AICtransient1 <- bbmle::AICtab(#habmodel22, 
  transientmodela, transientmodelb, transientmodelc, transientmodeld, transientmodele, transientmodelf, transientmodelg, transientmodelh)
capture.output(AICtransient1,file="AIC_transientmodel_stage1.doc")

#### transient round 2 ####
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
summary(transientmodelac)
saveRDS(transientmodelac, "transientmodelac.rds")
capture.output(summary(transientmodelac),file="transientmodelac.doc")
# capture.output(summary(transientmodelac),file="transientmodelac_nocovid.doc")

AICtransient1 <- bbmle::AICtab(habmodel22, transientmodelc, transientmodelac)
capture.output(AICtransient1,file="AIC_transientmodel_stage2.doc")

#BEST FIT: transientmodelac

####winter round 1####
wintermodela <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                          DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                          NDVI_end:winter +
                          (1|animal) + 
                          (0+NDVI_end|animal) + 
                          (0+POP_end|animal) +
                          (0+BUP_end|animal) + 
                          (0+DT_LT_end|animal) + 
                          (0+DT_MT_end|animal) + 
                          (0+DT_HT_end|animal) + 
                          (0+DT_NT_end|animal) + 
                          (0+DT_PS_end|animal) + 
                          (0+NDVI_end:winter|animal) + 
                          (1|step), data = coyote_steps, 
                        map = list(theta = factor(c(1:10, NA))),
                        start = list(theta = c(rep(0, 10), log(1e3))),
                        family = poisson)
saveRDS(wintermodela, "wintermodela.rds")
capture.output(summary(wintermodela),file="wintermodela.doc")

wintermodelb <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                          DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                          POP_end:winter +
                          (1|animal) + 
                          (0+NDVI_end|animal) + 
                          (0+POP_end|animal) +
                          (0+BUP_end|animal) + 
                          (0+DT_LT_end|animal) + 
                          (0+DT_MT_end|animal) + 
                          (0+DT_HT_end|animal) + 
                          (0+DT_NT_end|animal) + 
                          (0+DT_PS_end|animal) + 
                          (0+POP_end:winter|animal) + 
                          (1|step), data = coyote_steps, 
                        map = list(theta = factor(c(1:10, NA))),
                        start = list(theta = c(rep(0, 10), log(1e3))),
                        family = poisson)
saveRDS(wintermodelb, "wintermodelb.rds")
capture.output(summary(wintermodelb),file="wintermodelb.doc")

wintermodelc <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                          DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                          BUP_end:winter +
                          (1|animal) + 
                          (0+NDVI_end|animal) + 
                          (0+POP_end|animal) +
                          (0+BUP_end|animal) + 
                          (0+DT_LT_end|animal) + 
                          (0+DT_MT_end|animal) + 
                          (0+DT_HT_end|animal) + 
                          (0+DT_NT_end|animal) + 
                          (0+DT_PS_end|animal) + 
                          (0+BUP_end:winter|animal) + 
                          (1|step), data = coyote_steps, 
                        map = list(theta = factor(c(1:10, NA))),
                        start = list(theta = c(rep(0, 10), log(1e3))),
                        family = poisson)

saveRDS(wintermodelc, "wintermodelc.rds")
capture.output(summary(wintermodelc),file="wintermodelc.doc")

wintermodeld <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                          DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                          DT_LT_end:winter +
                          (1|animal) + 
                          (0+NDVI_end|animal) + 
                          (0+POP_end|animal) +
                          (0+BUP_end|animal) + 
                          (0+DT_LT_end|animal) + 
                          (0+DT_MT_end|animal) + 
                          (0+DT_HT_end|animal) + 
                          (0+DT_NT_end|animal) + 
                          (0+DT_PS_end|animal) + 
                          (0+DT_LT_end:winter|animal) + 
                          (1|step), data = coyote_steps, 
                        map = list(theta = factor(c(1:10, NA))),
                        start = list(theta = c(rep(0, 10), log(1e3))),
                        family = poisson)
saveRDS(wintermodeld, "wintermodeld.rds")
capture.output(summary(wintermodeld),file="wintermodeld.doc")

wintermodele <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                          DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                          DT_MT_end:winter +
                          (1|animal) + 
                          (0+NDVI_end|animal) + 
                          (0+POP_end|animal) +
                          (0+BUP_end|animal) + 
                          (0+DT_LT_end|animal) + 
                          (0+DT_MT_end|animal) + 
                          (0+DT_HT_end|animal) + 
                          (0+DT_NT_end|animal) + 
                          (0+DT_PS_end|animal) + 
                          (0+DT_MT_end:winter|animal) + 
                          (1|step), data = coyote_steps, 
                        map = list(theta = factor(c(1:10, NA))),
                        start = list(theta = c(rep(0, 10), log(1e3))),
                        family = poisson)
saveRDS(wintermodele, "wintermodele.rds")
capture.output(summary(wintermodele),file="wintermodele.doc")

wintermodelf <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                          DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                          DT_HT_end:winter +
                          (1|animal) + 
                          (0+NDVI_end|animal) + 
                          (0+POP_end|animal) +
                          (0+BUP_end|animal) + 
                          (0+DT_LT_end|animal) + 
                          (0+DT_MT_end|animal) + 
                          (0+DT_HT_end|animal) + 
                          (0+DT_NT_end|animal) + 
                          (0+DT_PS_end|animal) + 
                          (0+DT_HT_end:winter|animal) + 
                          (1|step), data = coyote_steps, 
                        map = list(theta = factor(c(1:10, NA))),
                        start = list(theta = c(rep(0, 10), log(1e3))),
                        family = poisson)
saveRDS(wintermodelf, "wintermodelf.rds")
capture.output(summary(wintermodelf),file="wintermodelf.doc")

wintermodelg <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                          DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                          DT_NT_end:winter +
                          (1|animal) + 
                          (0+NDVI_end|animal) + 
                          (0+POP_end|animal) +
                          (0+BUP_end|animal) + 
                          (0+DT_LT_end|animal) + 
                          (0+DT_MT_end|animal) + 
                          (0+DT_HT_end|animal) + 
                          (0+DT_NT_end|animal) + 
                          (0+DT_PS_end|animal) + 
                          (0+DT_NT_end:winter|animal) + 
                          (1|step), data = coyote_steps, 
                        map = list(theta = factor(c(1:10, NA))),
                        start = list(theta = c(rep(0, 10), log(1e3))),
                        family = poisson)
saveRDS(wintermodelg, "wintermodelg.rds")
capture.output(summary(wintermodelg),file="wintermodelg.doc")

wintermodelh <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                          DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                          DT_PS_end:winter +
                          (1|animal) + 
                          (0+NDVI_end|animal) + 
                          (0+POP_end|animal) +
                          (0+BUP_end|animal) + 
                          (0+DT_LT_end|animal) + 
                          (0+DT_MT_end|animal) + 
                          (0+DT_HT_end|animal) + 
                          (0+DT_NT_end|animal) + 
                          (0+DT_PS_end|animal) + 
                          (0+DT_PS_end:winter|animal) + 
                          (1|step), data = coyote_steps, 
                        map = list(theta = factor(c(1:10, NA))),
                        start = list(theta = c(rep(0, 10), log(1e3))),
                        family = poisson)
saveRDS(wintermodelh, "wintermodelh.rds")
capture.output(summary(wintermodelh),file="wintermodelh.doc")


habmodel22 <- readRDS("habmodel22_NDVI.BUP.PS.NT.POP.MT.LT.HT.rds")
wintermodela <- readRDS("wintermodela.rds") 
wintermodelb <- readRDS("wintermodelb.rds") 
wintermodelc <- readRDS("wintermodelc.rds") 
wintermodeld <- readRDS("wintermodeld.rds") 
wintermodele <- readRDS("wintermodele.rds") 
wintermodelf <- readRDS("wintermodelf.rds") 
wintermodelg <- readRDS("wintermodelg.rds") 
wintermodelh <- readRDS("wintermodelh.rds")

AICwinter1 <- bbmle::AICtab(habmodel22, wintermodela, wintermodelb, wintermodelc, wintermodeld, wintermodele, wintermodelf, wintermodelg, wintermodelh)
capture.output(AICwinter1,file="AIC_wintermodel_stage1.doc")

#####winter round 2 ############
wintermodel10<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                         DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + 
                         DT_NT_end:winter+
                         NDVI_end:winter+(1|animal) + 
                         (0+NDVI_end|animal) + 
                         (0+POP_end|animal) +
                         (0+BUP_end|animal) + 
                         (0+DT_LT_end|animal) + 
                         (0+DT_MT_end|animal) + 
                         (0+DT_HT_end|animal) + 
                         (0+DT_NT_end|animal) + 
                         (0+DT_PS_end|animal) +
                         (0+DT_NT_end:winter|animal) +  
                         (0+NDVI_end:winter|animal) + (1|step), data = coyote_steps, 
                       map = list(theta = factor(c(1:11, NA))),
                       start = list(theta = c(rep(0, 11), log(1e3))),
                       family = poisson)
saveRDS(wintermodel10, "wintermodel10.rds")
capture.output(summary(wintermodel10),file="wintermodel10.doc")


wintermodel12<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                         DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + 
                         DT_NT_end:winter+
                         DT_MT_end:winter +(1|animal) + 
                         (0+NDVI_end|animal) + 
                         (0+POP_end|animal) +
                         (0+BUP_end|animal) + 
                         (0+DT_LT_end|animal) + 
                         (0+DT_MT_end|animal) + 
                         (0+DT_HT_end|animal) + 
                         (0+DT_NT_end|animal) + 
                         (0+DT_PS_end|animal) +(0+DT_NT_end:winter|animal) +  (0+DT_MT_end:winter|animal) + (1|step), data = coyote_steps, 
                       map = list(theta = factor(c(1:11, NA))),
                       start = list(theta = c(rep(0, 11), log(1e3))),
                       family = poisson)
saveRDS(wintermodel12, "wintermodel12.rds")
capture.output(summary(wintermodel12),file="wintermodel12.doc")

wintermodel13<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                         DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + 
                         DT_NT_end:winter+
                         DT_HT_end:winter +(1|animal) + 
                         (0+NDVI_end|animal) + 
                         (0+POP_end|animal) +
                         (0+BUP_end|animal) + 
                         (0+DT_LT_end|animal) + 
                         (0+DT_MT_end|animal) + 
                         (0+DT_HT_end|animal) + 
                         (0+DT_NT_end|animal) + 
                         (0+DT_PS_end|animal) +(0+DT_NT_end:winter|animal) +  (0+DT_HT_end:winter|animal) + (1|step), data = coyote_steps, 
                       map = list(theta = factor(c(1:11, NA))),
                       start = list(theta = c(rep(0, 11), log(1e3))),
                       family = poisson)
saveRDS(wintermodel13, "wintermodel13.rds")
capture.output(summary(wintermodel13),file="wintermodel13.doc")

wintermodelg <- readRDS("wintermodelg.rds")
AICwinter2 <- bbmle::AICtab(wintermodelg, wintermodel10, wintermodel12, wintermodel13)
capture.output(AICwinter2,file="AIC_wintermodel_stage2.doc")


####winter round 3#######
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

wintermodel16<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                         DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end + 
                         DT_NT_end:winter+
                         NDVI_end:winter+
                         DT_HT_end:winter +
                         (1|animal) + 
                         (0+NDVI_end|animal) + 
                         (0+POP_end|animal) +
                         (0+BUP_end|animal) + 
                         (0+DT_LT_end|animal) + 
                         (0+DT_MT_end|animal) + 
                         (0+DT_HT_end|animal) + 
                         (0+DT_NT_end|animal) + 
                         (0+DT_PS_end|animal) +(0+DT_NT_end:winter|animal) +  (0+NDVI_end:winter|animal) +  (0+DT_HT_end:winter|animal) + (1|step), data = coyote_steps, 
                       map = list(theta = factor(c(1:12, NA))),
                       start = list(theta = c(rep(0, 12), log(1e3))),
                       family = poisson)
saveRDS(wintermodel16, "wintermodel16.rds")
capture.output(summary(wintermodel16),file="wintermodel16.doc")

AICwinter3 <- bbmle::AICtab(wintermodel10, wintermodel15, wintermodel16) #BEST FIT WINTERMODEL15
capture.output(AICwinter3,file="AIC_wintermodel_stage3.doc")



#BEST FIT: wintermodel15(aeg)

###breeding round 1####
breedingmodela <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                            DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                            NDVI_end:breeding + NDVI_end:dispersal + 
                            (1|animal) + 
                            (0+NDVI_end|animal) + 
                            (0+POP_end|animal) +
                            (0+BUP_end|animal) + 
                            (0+DT_LT_end|animal) + 
                            (0+DT_MT_end|animal) + 
                            (0+DT_HT_end|animal) + 
                            (0+DT_NT_end|animal) + 
                            (0+DT_PS_end|animal) + 
                            (0+NDVI_end:breeding|animal) + 
                            (0+NDVI_end:dispersal|animal) + 
                            (1|step), data = coyote_steps, 
                          map = list(theta = factor(c(1:11, NA))),
                          start = list(theta = c(rep(0, 11), log(1e3))),
                          family = poisson)
saveRDS(breedingmodela, "breedingmodela.rds")
capture.output(summary(breedingmodela),file="breedingmodela.doc")

breedingmodelb <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                            DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                            POP_end:breeding +POP_end:dispersal +
                            (1|animal) + 
                            (0+NDVI_end|animal) + 
                            (0+POP_end|animal) +
                            (0+BUP_end|animal) + 
                            (0+DT_LT_end|animal) + 
                            (0+DT_MT_end|animal) + 
                            (0+DT_HT_end|animal) + 
                            (0+DT_NT_end|animal) + 
                            (0+DT_PS_end|animal) + 
                            (0+POP_end:breeding|animal) + 
                            (0+POP_end:dispersal|animal) +
                            (1|step), data = coyote_steps, 
                          map = list(theta = factor(c(1:11, NA))),
                          start = list(theta = c(rep(0, 11), log(1e3))),
                          family = poisson)
saveRDS(breedingmodelb, "breedingmodelb.rds")
capture.output(summary(breedingmodelb),file="breedingmodelb.doc")

breedingmodelc <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                            DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                            BUP_end:breeding + BUP_end:dispersal +
                            (1|animal) + 
                            (0+NDVI_end|animal) + 
                            (0+POP_end|animal) +
                            (0+BUP_end|animal) + 
                            (0+DT_LT_end|animal) + 
                            (0+DT_MT_end|animal) + 
                            (0+DT_HT_end|animal) + 
                            (0+DT_NT_end|animal) + 
                            (0+DT_PS_end|animal) + 
                            (0+BUP_end:breeding|animal) + 
                            (0+BUP_end:dispersal|animal) + 
                            (1|step), data = coyote_steps, 
                          map = list(theta = factor(c(1:11, NA))),
                          start = list(theta = c(rep(0, 11), log(1e3))),
                          family = poisson)
saveRDS(breedingmodelc, "breedingmodelc.rds")
capture.output(summary(breedingmodelc),file="breedingmodelc.doc")

breedingmodeld <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                            DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                            DT_LT_end:breeding + DT_LT_end:dispersal +
                            (1|animal) + 
                            (0+NDVI_end|animal) + 
                            (0+POP_end|animal) +
                            (0+BUP_end|animal) + 
                            (0+DT_LT_end|animal) + 
                            (0+DT_MT_end|animal) + 
                            (0+DT_HT_end|animal) + 
                            (0+DT_NT_end|animal) + 
                            (0+DT_PS_end|animal) + 
                            (0+DT_LT_end:breeding|animal) + 
                            (0+DT_LT_end:dispersal|animal) + 
                            (1|step), data = coyote_steps, 
                          map = list(theta = factor(c(1:11, NA))),
                          start = list(theta = c(rep(0, 11), log(1e3))),
                          family = poisson)
saveRDS(breedingmodeld, "breedingmodeld.rds")
capture.output(summary(breedingmodeld),file="breedingmodeld.doc")

breedingmodele <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                            DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                            DT_MT_end:breeding + DT_MT_end:dispersal +
                            (1|animal) + 
                            (0+NDVI_end|animal) + 
                            (0+POP_end|animal) +
                            (0+BUP_end|animal) + 
                            (0+DT_LT_end|animal) + 
                            (0+DT_MT_end|animal) + 
                            (0+DT_HT_end|animal) + 
                            (0+DT_NT_end|animal) + 
                            (0+DT_PS_end|animal) + 
                            (0+DT_MT_end:breeding|animal) + 
                            (0+DT_MT_end:dispersal|animal) +
                            (1|step), data = coyote_steps, 
                          map = list(theta = factor(c(1:11, NA))),
                          start = list(theta = c(rep(0, 11), log(1e3))),
                          family = poisson)
saveRDS(breedingmodele, "breedingmodele.rds")
capture.output(summary(breedingmodele),file="breedingmodele.doc")

breedingmodelf <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                            DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                            DT_HT_end:breeding + DT_HT_end:dispersal +
                            (1|animal) + 
                            (0+NDVI_end|animal) + 
                            (0+POP_end|animal) +
                            (0+BUP_end|animal) + 
                            (0+DT_LT_end|animal) + 
                            (0+DT_MT_end|animal) + 
                            (0+DT_HT_end|animal) + 
                            (0+DT_NT_end|animal) + 
                            (0+DT_PS_end|animal) + 
                            (0+DT_HT_end:breeding|animal) + 
                            (0+DT_HT_end:dispersal|animal) + 
                            (1|step), data = coyote_steps, 
                          map = list(theta = factor(c(1:11, NA))),
                          start = list(theta = c(rep(0, 11), log(1e3))),
                          family = poisson)
saveRDS(breedingmodelf, "breedingmodelf.rds")
capture.output(summary(breedingmodelf),file="breedingmodelf.doc")

breedingmodelg <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                            DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                            DT_NT_end:breeding + DT_NT_end:dispersal +
                            (1|animal) + 
                            (0+NDVI_end|animal) + 
                            (0+POP_end|animal) +
                            (0+BUP_end|animal) + 
                            (0+DT_LT_end|animal) + 
                            (0+DT_MT_end|animal) + 
                            (0+DT_HT_end|animal) + 
                            (0+DT_NT_end|animal) + 
                            (0+DT_PS_end|animal) + 
                            (0+DT_NT_end:breeding|animal) + 
                            (0+DT_NT_end:dispersal|animal) + 
                            (1|step), data = coyote_steps, 
                          map = list(theta = factor(c(1:11, NA))),
                          start = list(theta = c(rep(0, 11), log(1e3))),
                          family = poisson)
saveRDS(breedingmodelg, "breedingmodelg.rds")
capture.output(summary(breedingmodelg),file="breedingmodelg.doc")

breedingmodelh <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                            DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                            DT_PS_end:breeding + DT_PS_end:dispersal +
                            (1|animal) + 
                            (0+NDVI_end|animal) + 
                            (0+POP_end|animal) +
                            (0+BUP_end|animal) + 
                            (0+DT_LT_end|animal) + 
                            (0+DT_MT_end|animal) + 
                            (0+DT_HT_end|animal) + 
                            (0+DT_NT_end|animal) + 
                            (0+DT_PS_end|animal) + 
                            (0+DT_PS_end:breeding|animal) + 
                            (0+DT_PS_end:dispersal|animal) + 
                            (1|step), data = coyote_steps, 
                          map = list(theta = factor(c(1:11, NA))),
                          start = list(theta = c(rep(0, 11), log(1e3))),
                          family = poisson)
saveRDS(breedingmodelh, "breedingmodelh.rds")
capture.output(summary(breedingmodelh),file="breedingmodelh.doc")

habmodel22 <- readRDS("habmodel22_NDVI.BUP.PS.NT.POP.MT.LT.HT.rds")
breedingmodela <- readRDS("breedingmodela.rds") 
breedingmodelb <- readRDS("breedingmodelb.rds") 
breedingmodelc <- readRDS("breedingmodelc.rds") 
breedingmodeld <- readRDS("breedingmodeld.rds") 
breedingmodele <- readRDS("breedingmodele.rds") 
breedingmodelf <- readRDS("breedingmodelf.rds") 
breedingmodelg <- readRDS("breedingmodelg.rds") 
breedingmodelh <- readRDS("breedingmodelh.rds")

AICbreeding1 <- bbmle::AICtab(habmodel22, breedingmodela, breedingmodelb, breedingmodelc, breedingmodeld, breedingmodele, breedingmodelf, breedingmodelg, breedingmodelh)
capture.output(AICbreeding1,file="AIC_breedingmodel_stage1.doc")
##PASS: g. h. c. d. f. a. 
########breeding round 2 ##########

breedingmodel10<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                           DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                           DT_NT_end:breeding + DT_NT_end:dispersal+
                           BUP_end:breeding +BUP_end:dispersal+
                           (1|animal) + 
                           (0+NDVI_end|animal) + 
                           (0+POP_end|animal) +
                           (0+BUP_end|animal) + 
                           (0+DT_LT_end|animal) + 
                           (0+DT_MT_end|animal) + 
                           (0+DT_HT_end|animal) + 
                           (0+DT_NT_end|animal) + 
                           (0+DT_PS_end|animal) + 
                           (0+DT_NT_end:breeding|animal) + (0+DT_NT_end:dispersal|animal) + 
                           (0+BUP_end:breeding|animal) +  (0+BUP_end:dispersal|animal) + 
                           (1|step), data = coyote_steps, 
                         map = list(theta = factor(c(1:13, NA))),
                         start = list(theta = c(rep(0, 13), log(1e3))),
                         family = poisson)
saveRDS(breedingmodel10, "breedingmodel10.rds")
capture.output(summary(breedingmodel10),file="breedingmodel10.doc")

breedingmodel11<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                           DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                           DT_NT_end:breeding + DT_NT_end:dispersal+
                           DT_LT_end:breeding +DT_LT_end:dispersal+
                           (1|animal) + 
                           (0+NDVI_end|animal) + 
                           (0+POP_end|animal) +
                           (0+BUP_end|animal) + 
                           (0+DT_LT_end|animal) + 
                           (0+DT_MT_end|animal) + 
                           (0+DT_HT_end|animal) + 
                           (0+DT_NT_end|animal) + 
                           (0+DT_PS_end|animal) + 
                           (0+DT_NT_end:breeding|animal) + (0+DT_NT_end:dispersal|animal) + 
                           (0+DT_LT_end:breeding|animal) +  (0+DT_LT_end:dispersal|animal) + 
                           (1|step), data = coyote_steps, 
                         map = list(theta = factor(c(1:13, NA))),
                         start = list(theta = c(rep(0, 13), log(1e3))),
                         family = poisson)
saveRDS(breedingmodel11, "breedingmodel11.rds")
capture.output(summary(breedingmodel11),file="breedingmodel11.doc")

breedingmodel12<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                           DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                           DT_NT_end:breeding + DT_NT_end:dispersal+
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
                           (0+DT_NT_end:breeding|animal) + (0+DT_NT_end:dispersal|animal) + 
                           (0+DT_HT_end:breeding|animal) +  (0+DT_HT_end:dispersal|animal) + 
                           (1|step), data = coyote_steps, 
                         map = list(theta = factor(c(1:13, NA))),
                         start = list(theta = c(rep(0, 13), log(1e3))),
                         family = poisson)
saveRDS(breedingmodel12, "breedingmodel12.rds")
capture.output(summary(breedingmodel12),file="breedingmodel12.doc")

breedingmodel13<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                           DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                           DT_NT_end:breeding + DT_NT_end:dispersal+
                           NDVI_end:breeding +NDVI_end:dispersal+
                           (1|animal) + 
                           (0+NDVI_end|animal) + 
                           (0+POP_end|animal) +
                           (0+BUP_end|animal) + 
                           (0+DT_LT_end|animal) + 
                           (0+DT_MT_end|animal) + 
                           (0+DT_HT_end|animal) + 
                           (0+DT_NT_end|animal) + 
                           (0+DT_PS_end|animal) + 
                           (0+DT_NT_end:breeding|animal) + (0+DT_NT_end:dispersal|animal) + 
                           (0+NDVI_end:breeding|animal) +  (0+NDVI_end:dispersal|animal) + 
                           (1|step), data = coyote_steps, 
                         map = list(theta = factor(c(1:13, NA))),
                         start = list(theta = c(rep(0, 13), log(1e3))),
                         family = poisson)
saveRDS(breedingmodel13, "breedingmodel13.rds")
capture.output(summary(breedingmodel13),file="breedingmodel13.doc")

breedingmodel14<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                           DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                           DT_NT_end:breeding + DT_NT_end:dispersal+
                           DT_PS_end:breeding +DT_PS_end:dispersal+
                           (1|animal) + 
                           (0+NDVI_end|animal) + 
                           (0+POP_end|animal) +
                           (0+BUP_end|animal) + 
                           (0+DT_LT_end|animal) + 
                           (0+DT_MT_end|animal) + 
                           (0+DT_HT_end|animal) + 
                           (0+DT_NT_end|animal) + 
                           (0+DT_PS_end|animal) + 
                           (0+DT_NT_end:breeding|animal) + (0+DT_NT_end:dispersal|animal) + 
                           (0+DT_PS_end:breeding|animal) +  (0+DT_PS_end:dispersal|animal) + 
                           (1|step), data = coyote_steps, 
                         map = list(theta = factor(c(1:13, NA))),
                         start = list(theta = c(rep(0, 13), log(1e3))),
                         family = poisson)
saveRDS(breedingmodel14, "breedingmodel14.rds")
capture.output(summary(breedingmodel14),file="breedingmodel14.doc")

breedingmodelg <- readRDS("breedingmodelg.rds")
AICbreeding1 <- bbmle::AICtab(breedingmodelg, breedingmodel10, breedingmodel11, breedingmodel12, breedingmodel13, breedingmodel14)
capture.output(AICbreeding1,file="AIC_breedingmode2_stage1.doc")

##PASS: all of them best: 10
####breeding round 3#####
breedingmodel17<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                           DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                           DT_NT_end:breeding + DT_NT_end:dispersal+
                           BUP_end:breeding + BUP_end:dispersal+
                           NDVI_end:breeding +NDVI_end:dispersal+ 
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
                           (1|step), data = coyote_steps, 
                         map = list(theta = factor(c(1:15, NA))),
                         start = list(theta = c(rep(0, 15), log(1e3))),
                         family = poisson)
saveRDS(breedingmodel17, "breedingmodel17.rds")
capture.output(summary(breedingmodel17),file="breedingmodel17.doc")

breedingmodel18<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                           DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                           DT_NT_end:breeding + DT_NT_end:dispersal+
                           BUP_end:breeding + BUP_end:dispersal+
                           DT_LT_end:breeding +DT_LT_end:dispersal+ 
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
                           (0+DT_LT_end:dispersal|animal) + (0+DT_LT_end:breeding|animal) +  
                           (1|step), data = coyote_steps, 
                         map = list(theta = factor(c(1:15, NA))),
                         start = list(theta = c(rep(0, 15), log(1e3))),
                         family = poisson)
saveRDS(breedingmodel18, "breedingmodel18.rds")
capture.output(summary(breedingmodel18),file="breedingmodel18.doc")

breedingmodel19<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                           DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                           DT_NT_end:breeding + DT_NT_end:dispersal+
                           BUP_end:breeding + BUP_end:dispersal+
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
                           (0+DT_HT_end:dispersal|animal) + (0+DT_HT_end:breeding|animal) +  
                           (1|step), data = coyote_steps, 
                         map = list(theta = factor(c(1:15, NA))),
                         start = list(theta = c(rep(0, 15), log(1e3))),
                         family = poisson)
saveRDS(breedingmodel19, "breedingmodel19.rds")
capture.output(summary(breedingmodel19),file="breedingmodel19.doc")

breedingmodel20<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                           DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                           DT_NT_end:breeding + DT_NT_end:dispersal+
                           BUP_end:breeding + BUP_end:dispersal+
                           DT_PS_end:breeding +DT_PS_end:dispersal+ 
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
                           (0+DT_PS_end:dispersal|animal) + (0+DT_PS_end:breeding|animal) +  
                           (1|step), data = coyote_steps, 
                         map = list(theta = factor(c(1:15, NA))),
                         start = list(theta = c(rep(0, 15), log(1e3))),
                         family = poisson)
saveRDS(breedingmodel20, "breedingmodel20.rds")
capture.output(summary(breedingmodel20),file="breedingmodel20.doc")

AICbreeding1 <- bbmle::AICtab(breedingmodel10, breedingmodel17, breedingmodel18, breedingmodel19, breedingmodel20)
capture.output(AICbreeding1,file="AIC_breedingmode2_stage3.doc")



####breedingmodel round 4####

breedingmodel21<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                           DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                           DT_NT_end:breeding + DT_NT_end:dispersal+
                           BUP_end:breeding + BUP_end:dispersal+
                           NDVI_end:breeding +NDVI_end:dispersal+  
                           DT_PS_end:breeding +DT_PS_end:dispersal+ 
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
                           (1|step), data = coyote_steps, 
                         map = list(theta = factor(c(1:17, NA))),
                         start = list(theta = c(rep(0, 17), log(1e3))),
                         family = poisson)
saveRDS(breedingmodel21, "breedingmodel21.rds")
capture.output(summary(breedingmodel21),file="breedingmodel21.doc")


breedingmodel22<-glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                           DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                           DT_NT_end:breeding + DT_NT_end:dispersal+
                           BUP_end:breeding + BUP_end:dispersal+
                           NDVI_end:breeding +NDVI_end:dispersal+  
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
                           (0+DT_HT_end:dispersal|animal) + (0+DT_HT_end:breeding|animal) +  
                           (1|step), data = coyote_steps, 
                         map = list(theta = factor(c(1:17, NA))),
                         start = list(theta = c(rep(0, 17), log(1e3))),
                         family = poisson)
saveRDS(breedingmodel22, "breedingmodel22.rds")
capture.output(summary(breedingmodel22),file="breedingmodel22.doc")

AICbreeding1 <- bbmle::AICtab(breedingmodel17, breedingmodel21, breedingmodel22)
capture.output(AICbreeding1,file="AIC_breedingmode2_stage4.doc")



###breedingmodel round 5####
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

AICbreeding1 <- bbmle::AICtab(breedingmodel21, breedingmodel25)
capture.output(AICbreeding1,file="AIC_breedingmode2_stage5.doc")

##best breeding model breedingmodel25

####male round 1####
malemodela <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                        DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                        NDVI_end:male +
                        (1|animal) + 
                        (0+NDVI_end|animal) + 
                        (0+POP_end|animal) +
                        (0+BUP_end|animal) + 
                        (0+DT_LT_end|animal) + 
                        (0+DT_MT_end|animal) + 
                        (0+DT_HT_end|animal) + 
                        (0+DT_NT_end|animal) + 
                        (0+DT_PS_end|animal) + 
                        (0+NDVI_end:male|animal) + 
                        (1|step), data = coyote_steps, 
                      map = list(theta = factor(c(1:10, NA))),
                      start = list(theta = c(rep(0, 10), log(1e3))),
                      family = poisson)
saveRDS(malemodela, "malemodela.rds")
capture.output(summary(malemodela),file="malemodela.doc")

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

malemodelc <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                        DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                        BUP_end:male +
                        (1|animal) + 
                        (0+NDVI_end|animal) + 
                        (0+POP_end|animal) +
                        (0+BUP_end|animal) + 
                        (0+DT_LT_end|animal) + 
                        (0+DT_MT_end|animal) + 
                        (0+DT_HT_end|animal) + 
                        (0+DT_NT_end|animal) + 
                        (0+DT_PS_end|animal) + 
                        (0+BUP_end:male|animal) + 
                        (1|step), data = coyote_steps, 
                      map = list(theta = factor(c(1:10, NA))),
                      start = list(theta = c(rep(0, 10), log(1e3))),
                      family = poisson)

saveRDS(malemodelc, "malemodelc.rds")
capture.output(summary(malemodelc),file="malemodelc.doc")

malemodeld <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                        DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                        DT_LT_end:male +
                        (1|animal) + 
                        (0+NDVI_end|animal) + 
                        (0+POP_end|animal) +
                        (0+BUP_end|animal) + 
                        (0+DT_LT_end|animal) + 
                        (0+DT_MT_end|animal) + 
                        (0+DT_HT_end|animal) + 
                        (0+DT_NT_end|animal) + 
                        (0+DT_PS_end|animal) + 
                        (0+DT_LT_end:male|animal) + 
                        (1|step), data = coyote_steps, 
                      map = list(theta = factor(c(1:10, NA))),
                      start = list(theta = c(rep(0, 10), log(1e3))),
                      family = poisson)
saveRDS(malemodeld, "malemodeld.rds")
capture.output(summary(malemodeld),file="malemodeld.doc")

malemodele <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                        DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                        DT_MT_end:male +
                        (1|animal) + 
                        (0+NDVI_end|animal) + 
                        (0+POP_end|animal) +
                        (0+BUP_end|animal) + 
                        (0+DT_LT_end|animal) + 
                        (0+DT_MT_end|animal) + 
                        (0+DT_HT_end|animal) + 
                        (0+DT_NT_end|animal) + 
                        (0+DT_PS_end|animal) + 
                        (0+DT_MT_end:male|animal) + 
                        (1|step), data = coyote_steps, 
                      map = list(theta = factor(c(1:10, NA))),
                      start = list(theta = c(rep(0, 10), log(1e3))),
                      family = poisson)
saveRDS(malemodele, "malemodele.rds")
capture.output(summary(malemodele),file="malemodele.doc")

malemodelf <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                        DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                        DT_HT_end:male +
                        (1|animal) + 
                        (0+NDVI_end|animal) + 
                        (0+POP_end|animal) +
                        (0+BUP_end|animal) + 
                        (0+DT_LT_end|animal) + 
                        (0+DT_MT_end|animal) + 
                        (0+DT_HT_end|animal) + 
                        (0+DT_NT_end|animal) + 
                        (0+DT_PS_end|animal) + 
                        (0+DT_HT_end:male|animal) + 
                        (1|step), data = coyote_steps, 
                      map = list(theta = factor(c(1:10, NA))),
                      start = list(theta = c(rep(0, 10), log(1e3))),
                      family = poisson)
saveRDS(malemodelf, "malemodelf.rds")
capture.output(summary(malemodelf),file="malemodelf.doc")

malemodelg <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                        DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                        DT_NT_end:male +
                        (1|animal) + 
                        (0+NDVI_end|animal) + 
                        (0+POP_end|animal) +
                        (0+BUP_end|animal) + 
                        (0+DT_LT_end|animal) + 
                        (0+DT_MT_end|animal) + 
                        (0+DT_HT_end|animal) + 
                        (0+DT_NT_end|animal) + 
                        (0+DT_PS_end|animal) + 
                        (0+DT_NT_end:male|animal) + 
                        (1|step), data = coyote_steps, 
                      map = list(theta = factor(c(1:10, NA))),
                      start = list(theta = c(rep(0, 10), log(1e3))),
                      family = poisson)
saveRDS(malemodelg, "malemodelg.rds")
capture.output(summary(malemodelg),file="malemodelg.doc")

malemodelh <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                        DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                        DT_PS_end:male +
                        (1|animal) + 
                        (0+NDVI_end|animal) + 
                        (0+POP_end|animal) +
                        (0+BUP_end|animal) + 
                        (0+DT_LT_end|animal) + 
                        (0+DT_MT_end|animal) + 
                        (0+DT_HT_end|animal) + 
                        (0+DT_NT_end|animal) + 
                        (0+DT_PS_end|animal) + 
                        (0+DT_PS_end:male|animal) + 
                        (1|step), data = coyote_steps, 
                      map = list(theta = factor(c(1:10, NA))),
                      start = list(theta = c(rep(0, 10), log(1e3))),
                      family = poisson)
saveRDS(malemodelh, "malemodelh.rds")
capture.output(summary(malemodelh),file="malemodelh.doc")

habmodel22 <- readRDS("habmodel22_NDVI.BUP.PS.NT.POP.MT.LT.HT.rds")
malemodela <- readRDS("malemodela.rds") 
malemodelb <- readRDS("malemodelb.rds") 
malemodelc <- readRDS("malemodelc.rds") 
malemodeld <- readRDS("malemodeld.rds") 
malemodele <- readRDS("malemodele.rds") ## here is where it crashed.
malemodelf <- readRDS("malemodelf.rds") 
malemodelg <- readRDS("malemodelg.rds") 
malemodelh <- readRDS("malemodelh.rds")

AICmale1 <- bbmle::AICtab(habmodel22, malemodela, malemodelb, malemodelc, malemodeld, malemodele, malemodelf, malemodelg, malemodelh)
capture.output(AICmale1,file="AIC_malemodel_stage1.doc")
#BEST FIT: malemodelb #BEST FIT: malemodelb
###adult round 1########
adultmodela <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                         DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                         NDVI_end:adult + 
                         (1|animal) + 
                         (0+NDVI_end|animal) + 
                         (0+POP_end|animal) +
                         (0+BUP_end|animal) + 
                         (0+DT_LT_end|animal) + 
                         (0+DT_MT_end|animal) + 
                         (0+DT_HT_end|animal) + 
                         (0+DT_NT_end|animal) + 
                         (0+DT_PS_end|animal) + 
                         (0+NDVI_end:adult|animal) + 
                         (1|step), data = coyote_steps, 
                       map = list(theta = factor(c(1:10, NA))),
                       start = list(theta = c(rep(0, 10), log(1e3))),
                       family = poisson)
saveRDS(adultmodela, "adultmodela.rds")
capture.output(summary(adultmodela),file="adultmodela.doc")

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

adultmodelc <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                         DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                         BUP_end:adult + 
                         (1|animal) + 
                         (0+NDVI_end|animal) + 
                         (0+POP_end|animal) +
                         (0+BUP_end|animal) + 
                         (0+DT_LT_end|animal) + 
                         (0+DT_MT_end|animal) + 
                         (0+DT_HT_end|animal) + 
                         (0+DT_NT_end|animal) + 
                         (0+DT_PS_end|animal) + 
                         (0+BUP_end:adult|animal) + 
                         (1|step), data = coyote_steps, 
                       map = list(theta = factor(c(1:10, NA))),
                       start = list(theta = c(rep(0, 10), log(1e3))),
                       family = poisson)
saveRDS(adultmodelc, "adultmodelc.rds")
capture.output(summary(adultmodelc),file="adultmodelc.doc")

adultmodeld <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                         DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                         DT_LT_end:adult + 
                         (1|animal) + 
                         (0+NDVI_end|animal) + 
                         (0+POP_end|animal) +
                         (0+BUP_end|animal) + 
                         (0+DT_LT_end|animal) + 
                         (0+DT_MT_end|animal) + 
                         (0+DT_HT_end|animal) + 
                         (0+DT_NT_end|animal) + 
                         (0+DT_PS_end|animal) + 
                         (0+DT_LT_end:adult|animal) +  
                         (1|step), data = coyote_steps, 
                       map = list(theta = factor(c(1:10, NA))),
                       start = list(theta = c(rep(0, 10), log(1e3))),
                       family = poisson)
saveRDS(adultmodeld, "adultmodeld.rds")
capture.output(summary(adultmodeld),file="adultmodeld.doc")

adultmodele <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                         DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                         DT_MT_end:adult +
                         (1|animal) + 
                         (0+NDVI_end|animal) + 
                         (0+POP_end|animal) +
                         (0+BUP_end|animal) + 
                         (0+DT_LT_end|animal) + 
                         (0+DT_MT_end|animal) + 
                         (0+DT_HT_end|animal) + 
                         (0+DT_NT_end|animal) + 
                         (0+DT_PS_end|animal) + 
                         (0+DT_MT_end:adult|animal) + 
                         (1|step), data = coyote_steps, 
                       map = list(theta = factor(c(1:10, NA))),
                       start = list(theta = c(rep(0, 10), log(1e3))),
                       family = poisson)
saveRDS(adultmodele, "adultmodele.rds")
capture.output(summary(adultmodele),file="adultmodele.doc")

adultmodelf <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                         DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                         DT_HT_end:adult + 
                         (1|animal) + 
                         (0+NDVI_end|animal) + 
                         (0+POP_end|animal) +
                         (0+BUP_end|animal) + 
                         (0+DT_LT_end|animal) + 
                         (0+DT_MT_end|animal) + 
                         (0+DT_HT_end|animal) + 
                         (0+DT_NT_end|animal) + 
                         (0+DT_PS_end|animal) + 
                         (0+DT_HT_end:adult|animal) + 
                         (1|step), data = coyote_steps, 
                       map = list(theta = factor(c(1:10, NA))),
                       start = list(theta = c(rep(0, 10), log(1e3))),
                       family = poisson)
saveRDS(adultmodelf, "adultmodelf.rds")
capture.output(summary(adultmodelf),file="adultmodelf.doc")

adultmodelg <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                         DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                         DT_NT_end:adult +
                         (1|animal) + 
                         (0+NDVI_end|animal) + 
                         (0+POP_end|animal) +
                         (0+BUP_end|animal) + 
                         (0+DT_LT_end|animal) + 
                         (0+DT_MT_end|animal) + 
                         (0+DT_HT_end|animal) + 
                         (0+DT_NT_end|animal) + 
                         (0+DT_PS_end|animal) + 
                         (0+DT_NT_end:adult|animal) + 
                         (1|step), data = coyote_steps, 
                       map = list(theta = factor(c(1:10, NA))),
                       start = list(theta = c(rep(0, 10), log(1e3))),
                       family = poisson)
saveRDS(adultmodelg, "adultmodelg.rds")
capture.output(summary(adultmodelg),file="adultmodelg.doc")

adultmodelh <- glmmTMB(used ~ -1 + NDVI_end + POP_end + BUP_end + 
                         DT_LT_end + DT_MT_end + DT_HT_end + DT_NT_end + DT_PS_end +
                         DT_PS_end:adult + 
                         (1|animal) + 
                         (0+NDVI_end|animal) + 
                         (0+POP_end|animal) +
                         (0+BUP_end|animal) + 
                         (0+DT_LT_end|animal) + 
                         (0+DT_MT_end|animal) + 
                         (0+DT_HT_end|animal) + 
                         (0+DT_NT_end|animal) + 
                         (0+DT_PS_end|animal) + 
                         (0+DT_PS_end:adult|animal) +  
                         (1|step), data = coyote_steps, 
                       map = list(theta = factor(c(1:10, NA))),
                       start = list(theta = c(rep(0, 10), log(1e3))),
                       family = poisson)
saveRDS(adultmodelh, "adultmodelh.rds")
capture.output(summary(adultmodelh),file="adultmodelh.doc")

habmodel22 <- readRDS("habmodel22_NDVI.BUP.PS.NT.POP.MT.LT.HT.rds")
adultmodela <- readRDS("adultmodela.rds") 
adultmodelb <- readRDS("adultmodelb.rds") 
adultmodelc <- readRDS("adultmodelc.rds") 
adultmodeld <- readRDS("adultmodeld.rds") 
adultmodele <- readRDS("adultmodele.rds") 
adultmodelf <- readRDS("adultmodelf.rds") 
adultmodelg <- readRDS("adultmodelg.rds") 
adultmodelh <- readRDS("adultmodelh.rds")

AICadult1 <- bbmle::AICtab(habmodel22, adultmodela, adultmodelb, adultmodelc, adultmodeld, adultmodele, adultmodelf, adultmodelg, adultmodelh)
capture.output(AICadult1,file="AIC_adultmodel_stage1.doc")

#BEST FIT: adultmodelb

###summary tables for paper####
readRDS("habmodel22_NDVI.BUP.PS.NT.POP.MT.LT.HT.rds")  -> habmodel22
readRDS("daymodel37.rds")  -> daymodel37
readRDS("transientmodelac.rds")  -> transientmodelac
readRDS("wintermodel15.rds")  -> wintermodel15
readRDS("breedingmodel25.rds")  -> breedingmodel25

habmodel22 -> base
daymodel37 -> day
transientmodelac -> soc
wintermodel15 -> cli
breedingmodel25 -> bio

pl <- c(
  `NDVI_end` = "Vegetation density",
  `BUP_end` = "Built-up density",
  `DT_PS_end` = "Distance to Public service lines",
  `DT_NT_end` = "Distance to Hiking trails",
  `POP_end` = "Human Population density",
  `DT_MT_end` = "Distance to Medium-traffic roads",
  `DT_LT_end` = "Distance to Low-traffic roads",
  `DT_HT_end` = "Distance to High-traffic roads",
  `NDVI_end:day` = "Vegetation:day",
  `BUP_end:day` = "Built-up:day",
  `DT_PS_end:day` = "Public service:day",
  `DT_NT_end:day` = "Hiking trails:day",
  `POP_end:day` = "Human Population:day",
  `DT_MT_end:day` = "Medium-traffic:day",
  `DT_LT_end:day` = "Low-traffic:day",
  `DT_HT_end:day` = "High-traffic:day",
  `NDVI_end:transient` = "Vegetation:transient",
  `BUP_end:transient` = "Built-up:transient",
  `DT_PS_end:transient` = "Public service:transient",
  `DT_NT_end:transient` = "Hiking trails:transient",
  `POP_end:transient` = "Human Population:transient",
  `DT_MT_end:transient` = "Medium-traffic:transient",
  `DT_LT_end:transient` = "Low-traffic:transient",
  `DT_HT_end:transient` = "High-traffic:transient",
  `NDVI_end:winter` = "Vegetation:winter",
  `BUP_end:winter` = "Built-up:winter",
  `DT_PS_end:winter` = "Public service:winter",
  `DT_NT_end:winter` = "Hiking trails:winter",
  `POP_end:winter` = "Human Population:winter",
  `DT_MT_end:winter` = "Medium-traffic:winter",
  `DT_LT_end:winter` = "Low-traffic:winter",
  `DT_HT_end:winter` = "High-traffic:winter",
  `NDVI_end:dispersal` = "Vegetation:dispersal",
  `BUP_end:dispersal` = "Built-up:dispersal",
  `DT_PS_end:dispersal` = "Public service:dispersal",
  `DT_NT_end:dispersal` = "Hiking trails:dispersal",
  `POP_end:dispersal` = "Human Population:dispersal",
  `DT_MT_end:dispersal` = "Medium-traffic:dispersal",
  `DT_LT_end:dispersal` = "Low-traffic:dispersal",
  `DT_HT_end:dispersal` = "High-traffic:dispersal",
  `NDVI_end:breeding` = "Vegetation:breeding",
  `BUP_end:breeding` = "Built-up:breeding",
  `DT_PS_end:breeding` = "Public service:breeding",
  `DT_NT_end:breeding` = "Hiking trails:breeding",
  `POP_end:breeding` = "Human Population:breeding",
  `DT_MT_end:breeding` = "Medium-traffic:breeding",
  `DT_LT_end:breeding` = "Low-traffic:breeding",
  `DT_HT_end:breeding` = "High-traffic:breeding"
)
tab_model(base,day,soc,cli,bio, 
          transform = NULL,
          pred.labels = pl, 
          dv.labels = c("Base model", "Diel cycle", "Social status", "Climate seasons", "biological seasons"),
          p.style = "stars",
          file = "SSFplot.doc")

b<-plot_model(base,
              show.values = TRUE,
              transform=NULL,
              pred.labels = pl, 
              value.offset = 0.45, 
              value.size = 3,
              dot.size = 1.5,
              line.size = .5,
              width = .5,
              vline.color = "black",
              show.intercept = FALSE
) + theme_sjplot()

d<-plot_model(day,
              show.values = TRUE,
              transform=NULL,
              pred.labels = pl, 
              value.offset = 0.45, 
              value.size = 3,
              dot.size = 1.5,
              line.size = .5,
              width = .5,
              vline.color = "black",
              show.intercept = FALSE
) + theme_sjplot()

s<-plot_model(soc,
              show.values = TRUE,
              transform=NULL,
              pred.labels = pl, 
              value.offset = 0.45, 
              value.size = 3,
              dot.size = 1.5,
              line.size = .5,
              width = .5,
              vline.color = "black",
              show.intercept = FALSE
) + theme_sjplot()


c<-plot_model(cli,
              show.values = TRUE,
              transform=NULL,
              pred.labels = pl, 
              value.offset = 0.45, 
              value.size = 3,
              dot.size = 1.5,
              line.size = .5,
              width = .5,
              vline.color = "black",
              show.intercept = FALSE
) + theme_sjplot()


bi<-plot_model(bio,
               show.values = TRUE,
               transform=NULL,
               pred.labels = pl, 
               value.offset = 0.45, 
               value.size = 3,
               dot.size = 1.5,
               line.size = .5,
               width = .5,
               vline.color = "black",
               show.intercept = FALSE
) + theme_sjplot()

library(ggpubr)
ggarrange(b, d, s, c, bi)

