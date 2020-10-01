

####################
#
#
# plot_define_groups.R
#
# 5/16/19
# Anika Petach
# 
#
#####################

require(dplyr)
library(lme4)
require(ggplot2)
require(GGally)
library(sjPlot)
library(fields)
require(car)
require(data.table) #for fread
require(coefplot)
library(cowplot)
library(ggthemes)
require(fields)
require(raster)
require(sp)

setwd("/Users/Anika/Documents/GradSchool/FIA_CompetitionModel/")

options(scipen = 6)

source("scripts/functions.R")


# read data for models ------------------------------
FIA_plot <- readRDS("output/fia_plot_foranalysis.RDS")


# Soil moisture -----------------------------------
#soilm <- raster("/Users/Anika/Downloads/w.full.200811/w.full.200811/w.full.200811.tif") # could be mm moisture, Nov 2008

#summary(soilm)

#soilm[soilm == -999] <- NA #recode -999 as NA

plot <- fread("data/raw/plot.csv", select = c('CN', 
                                              'LAT',
                                              'LON'))
plot$CN <- as.numeric(plot$CN)

FIA_plot <- merge(FIA_plot, plot, by.x = "pcn1", by.y = "CN", all.x = TRUE, all.y = FALSE)

coords <- FIA_plot[,c("LON","LAT")]
#sm <- raster::extract(soilm, coords)

#FIA_plot <- cbind(FIA_plot, sm)


# disturbance code ---------------------------------
cond <- fread("data/raw/cond.csv", select = c('PLT_CN', 
                                              'DSTRBCD1',
                                              'FORTYPCD',
                                              'STATECD'))
cond$PLT_CN <- as.numeric(cond$PLT_CN)
cond <- with(cond, cond[order(PLT_CN),])
cond <- cond[!duplicated(cond$PLT_CN), ] #remove duplication condition information

FIA_plot <- merge(FIA_plot, cond, by.x = "pcn1", by.y = 'PLT_CN', all.x = T, all.y = F)
# 
# dstrb$dtype <- ifelse(dstrb$DSTRBCD1 == 0, 'none',
#                       ifelse(dstrb$DSTRBCD1 >= 10 & dstrb$DSTRBCD1 < 20, 'insect',
#                              ifelse(dstrb$DSTRBCD1 >=20 & dstrb$DSTRBCD1 < 30, 'disease',
#                                     ifelse(dstrb$DSTRBCD1 >=30 & dstrb$DSTRBCD1 < 40, 'fire', 
#                                            ifelse(dstrb$DSTRBCD1 >=40 & dstrb$DSTRBCD1 < 50, 'animal',
#                                                   ifelse(dstrb$DSTRBCD1 >=50 & dstrb$DSTRBCD1 < 60, 'weather',
#                                                          ifelse(dstrb$DSTRBCD1 == 60, 'vegetation',
#                                                                 ifelse(dstrb$DSTRBCD1 == 70, 'unknown',
#                                                                        ifelse(dstrb$DSTRBCD1 == 80, 'human',
#                                                                               ifelse(dstrb$DSTRBCD1 >= 90 & dstrb$DSTRBCD1 < 100, 'geology',
#                                                                                      NA))))))))))



# add groups ------------------------------------------------------
coords <- FIA_plot[,c("LON","LAT")]


###
# Define age class
FIA_plot <- FIA_plot %>%
  mutate(YOU0OLD1 = case_when(stdage < 60 ~ 0,
                              stdage >= 60 ~ 1))
FIA_plot$YOU0OLD1 <- as.factor(FIA_plot$YOU0OLD1)

###
# Define N dep class
Ndep <- raster("/Users/Anika/Documents/GradSchool/FIA_EstimateProj/output/N_dep.grd")

dep <- raster::extract(Ndep, coords)

FIA_plot <- cbind(FIA_plot, dep)

FIA_plot <- FIA_plot %>%
  mutate(NDEP = case_when(dep > 0 & dep < 1 ~ 0,
                          dep >= 1 & dep < 2 ~ 1,
                          dep >= 2 & dep < 3 ~ 2,
                          dep >= 3 & dep < 4 ~ 3,
                          dep >= 4 & dep < 5 ~ 4,
                          dep >= 5 & dep < 6 ~ 5,
                          dep >= 6 & dep < 7 ~ 6,
                          dep >= 7 & dep < 8 ~ 7,
                          dep >= 8 & dep < 9 ~ 8,
                          dep >= 9 ~ 9))
FIA_plot$NDEP <- as.factor(FIA_plot$NDEP)


###
# Define soil moisture class
us_soil <- raster("output/smap_proj.grd")

sm <- raster::extract(us_soil, coords)

FIA_plot <- cbind(FIA_plot, sm)

FIA_plot <- FIA_plot %>%
  mutate(SM = case_when(sm > 0 & sm < 0.1 ~ 0,
                        sm >= 0.1 & sm < 0.2 ~ 0.1,
                        sm >= 0.2 & sm < 0.3 ~ 0.2,
                        sm >= 0.3 & sm < 0.4 ~ 0.3,
                        sm >= 0.4 & sm < 0.5 ~ 0.4,
                        sm >= 0.5 & sm < 0.6 ~ 0.5,
                        sm >= 0.6 & sm < 0.7 ~ 0.6,
                        sm >= 0.7 ~ 0.7))
FIA_plot$SM <- as.factor(FIA_plot$SM)


###
# Define soil texture classes, soilgrids
texture <- raster("/Users/Anika/Downloads/TEXMHT_M_sl4_250m_ll.tif")
texture <- raster::extract(texture, coords)
#require(GSIF)
#data(soil.legends)

FIA_plot <- cbind(FIA_plot, texture)
FIA_plot$texture <- as.factor(FIA_plot$texture)


###
# Define forest classes

FIA_plot$fortypsm <- trunc(FIA_plot$FORTYPCD/10)*10
FIA_plot$fortypsm <- as.integer(FIA_plot$fortypsm)
FIA_plot[(FIA_plot$fortypsm %in% c(510, 520)), "fortypsm"] <- 500
FIA_plot[FIA_plot$fortypsm %in% c(720), "fortypsm"] <- 700
FIA_plot[FIA_plot$fortypsm %in% c(930), "fortypsm"] <- 920
FIA_plot$fortypsm <- as.factor(FIA_plot$fortypsm)


###
# Disturbance

FIA_plot$dtype <- ifelse(FIA_plot$DSTRBCD1 == 0, 'none',
                         ifelse(FIA_plot$DSTRBCD1 >= 10 & FIA_plot$DSTRBCD1 < 20, 'insect',
                                ifelse(FIA_plot$DSTRBCD1 >=20 & FIA_plot$DSTRBCD1 < 30, 'disease',
                                       ifelse(FIA_plot$DSTRBCD1 >=30 & FIA_plot$DSTRBCD1 < 40, 'fire', 
                                              ifelse(FIA_plot$DSTRBCD1 >=40 & FIA_plot$DSTRBCD1 < 50, 'animal',
                                                     ifelse(FIA_plot$DSTRBCD1 >=50 & FIA_plot$DSTRBCD1 < 60, 'weather',
                                                            ifelse(FIA_plot$DSTRBCD1 == 60, 'vegetation',
                                                                   ifelse(FIA_plot$DSTRBCD1 == 70, 'unknown',
                                                                          ifelse(FIA_plot$DSTRBCD1 == 80, 'human',
                                                                                 ifelse(FIA_plot$DSTRBCD1 >= 90 & FIA_plot$DSTRBCD1 < 100, 'geology',
                                                                                        NA))))))))))

FIA_plot <- FIA_plot %>%
  mutate(dtype = case_when(DSTRBCD1 == 0 ~ 'none',
                           DSTRBCD1 >= 10 & DSTRBCD1 < 20 ~ 'insect',
                           DSTRBCD1 >= 20 & DSTRBCD1 < 30 ~ 'disease',
                           DSTRBCD1 >= 30 & DSTRBCD1 < 40 ~ 'fire',
                           DSTRBCD1 >= 40 & DSTRBCD1 < 50 ~ 'animal',
                           DSTRBCD1 >= 50 & DSTRBCD1 < 60 ~ 'weather',
                           DSTRBCD1 == 60 ~ 'vegetation',
                           DSTRBCD1 == 70 ~ 'unknown',
                           DSTRBCD1 == 80 ~ 'human',
                           DSTRBCD1 >= 90 & DSTRBCD1 < 100 ~ 'geology',
                           is.na(DSTRBCD1) ~ 'NA'))
FIA_plot$dtype <- as.factor(FIA_plot$dtype)

FIA_plot <- FIA_plot %>%
  mutate(NDEP2 = case_when(dep > 0 & dep < 2 ~ "0-2",
                           dep >= 2 & dep < 4 ~ "2-4",
                           dep >= 4 & dep < 6 ~ "4-6",
                           dep >= 6 & dep < 8 ~ "6-8",
                           dep >= 8 ~ "8+"))
FIA_plot$NDEP2 <- as.factor(FIA_plot$NDEP2)

FIA_plot <- FIA_plot %>%
  mutate(Cunder = case_when(understoryC1 > 0 & understoryC1 < 0.5 ~ "0-0.5",
                            understoryC1 >= 0.5 & understoryC1 < 1 ~ "0.5-1",
                            understoryC1 >= 1 & understoryC1 < 1.5 ~ "1-1.5",
                            understoryC1 >= 1.5 & understoryC1 < 2 ~ "1.5-2",
                            understoryC1 >= 8 ~ "2+"))
FIA_plot$Cunder <- as.factor(FIA_plot$Cunder)

#get MAT and MAP (code snippet from Ben)
clim <- getData("worldclim",var="bio",res=10)
clim <- clim[[c(1,12)]]
names(clim) <- c("MAT","MAP")
#coords <- DATA[,c("longitude","latitude")]
coords <- FIA_plot[,c("LON","LAT")]
points <- SpatialPoints(coords,proj4string=clim@crs)
climvals <- extract(clim,points)
FIA_plot <- cbind(FIA_plot,climvals)
FIA_plot$MAT <- (FIA_plot$MAT/10)


saveRDS(FIA_plot, "output/FIA_plot_withgroups.RDS")
