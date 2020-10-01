

###################
#
#
# Add groupings to individual level data
# Growth rate, recruitment, survival
# Decid, N dep, soil moisture, MR, traits
#
# Author: Anika Petach
# Date: 3/25/19
#
#
#####################

require(dplyr)
library(lme4)
library(nlme)
require(ggplot2)
require(GGally)
require(dplyr)
library(DataExplorer)
library(MASS)
library(sjPlot)
library(fields)
require(car)
require(data.table) #for fread
require(coefplot)
require(raster)
require(tidyverse)

options(scipen=6)

setwd("/Users/Anika/Documents/GradSchool/FIA_CompetitionModel/")

source("scripts/functions.R")

# Growth rate ------------------------------------
gs <- readRDS("output/nci_formodels_growth.RDS") #made in growth_prep.R

spdata <- read.csv("/Users/Anika/Documents/GradSchool/FIA_EstimateProj/ACGCA_FIAv51_from_JWL/REF_SPECIES_jwl.csv")
spdata <- spdata %>% dplyr::select(SPCD, COMMON_NAME, GENUS, SPECIES, JENKINS_SPGRPCD) #get biomass parameters if needed

gs <- merge(gs, spdata, by.x="spcd", by.y="SPCD", all.x=T, all.y=F)
gs$Decid0Ever1 <- as.factor(gs$Decid0Ever1)

###
# Define evergreen vs deciduous
classten <- read.csv("data/raw/Compiled_decid0ever1.csv")
classten <- classten[, c("spnames", "AccSpeciesID", "Decid0Ever1")]

gs <- gs %>% tidyr::unite(gensp, GENUS, SPECIES, sep = " ", remove=FALSE) #get col to merge classten on

gs <- merge(gs, classten, by.x="gensp", by.y="spnames", all.x=T)


###
# Define age class
gs <- gs %>%
  mutate(YOU0OLD1 = case_when(STDAGE < 60 ~ 0,
                              STDAGE >= 60 ~ 1))
gs$YOU0OLD1 <- as.factor(gs$YOU0OLD1)

###
# Define canopy class
gs <- gs %>%
  mutate(NON0CAN1 = case_when(cclcd1 < 4 ~ 1,
                              cclcd1 >= 4 ~ 0))
gs$NON0CAN1 <- as.factor(gs$NON0CAN1)

###
# Define N dep class
Ndep <- raster("/Users/Anika/Documents/GradSchool/FIA_EstimateProj/output/N_dep.grd")

coords <- gs[,c("LON","LAT")]
dep <- raster::extract(Ndep, coords)

gs <- cbind(gs, dep)

gs <- gs %>%
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
gs$NDEP <- as.factor(gs$NDEP)


###
# Define soil moisture class
us_soil <- raster("output/smap_proj.grd")

coords <- gs[,c("LON","LAT")]
sm <- raster::extract(us_soil, coords)

gs <- cbind(gs, sm)

gs <- gs %>%
  mutate(SM = case_when(sm > 0 & sm < 0.1 ~ 0,
                        sm >= 0.1 & sm < 0.2 ~ 0.1,
                        sm >= 0.2 & sm < 0.3 ~ 0.2,
                        sm >= 0.3 & sm < 0.4 ~ 0.3,
                        sm >= 0.4 & sm < 0.5 ~ 0.4,
                        sm >= 0.5 & sm < 0.6 ~ 0.5,
                        sm >= 0.6 & sm < 0.7 ~ 0.6,
                        sm >= 0.7 ~ 0.7))
gs$SM <- as.factor(gs$SM)


###
# Define soil texture classes, soilgrids
texture <- raster("/Users/Anika/Downloads/TEXMHT_M_sl4_250m_ll.tif")
texture <- raster::extract(texture, coords)
#require(GSIF)
#data(soil.legends)

gs <- cbind(gs, texture)
gs$texture <- as.factor(gs$texture)

###
# Define MR type
phil <- read.csv("data/raw/Phillips_AM_EM_data.csv")
phil <- phil[,3:5]
colnames(phil) <- c("GENUS", "SPECIES", "Mycostatus")
phil$spnames <- paste(phil$GENUS, phil$SPECIES)

splist <- gs[,c("GENUS","SPECIES")]
#splist <- unique(gs[,"gensp"])
splist$spnames <- paste(splist$GENUS, splist$SPECIES)
splist <- splist[!duplicated(splist$spnames), ]

mrlist0 <- merge(splist, phil, by="spnames", all.x=T, all.y=F) #used to be spnames
mrlist1 <- mrlist0[!is.na(mrlist0$Mycostatus), ]
mrlist1 <- mrlist1[,c(1:3,6)] #keep spnames, species, genus, mycostatus
colnames(mrlist1) <- c("spnames", "GENUS", "SPECIES", "Mycostatus")
mrlist2 <- mrlist1 %>% 
  group_by(GENUS) %>% 
  summarize(myco_phil = first(Mycostatus)) #MR for genus level (can merge back on to mrlist)

gs <- merge(gs, mrlist2, by="GENUS", all.x=T, all.y=F)



###
# Define trait classes
usda <- readRDS("data/derived/usda_numeric.RDS")
gs <- merge(gs, usda, by.x="spcd", by.y="SPCD", all.x=T, all.y=F)

gs <- gs %>%
  mutate(Tcnratio = case_when(CNRatio == 1 ~ "high",
                              STDAGE > 1 ~ "low"),
         Tallelopathy = case_when(Allelopath == 1 ~ "yes",
                                  Allelopath == 0 ~ "no"),
         Tleafret = case_when(LeafRetention == 1 ~ "yes",
                              LeafRetention == 0 ~ "no"),
         Tleaflife = case_when(Lifespan == 1 ~ "long",
                               Lifespan > 1 ~ "short"),
         Tfertilityreq = case_when(FerilityReq == 1 ~ "high",
                                   FerilityReq > 1 ~ "low"),
         Tdroughttol = case_when(DroughtTolerance <= 2 ~ "high",
                                 DroughtTolerance > 2 ~ "low"),
         Tshadetol = case_when(ShadeTolerance == 1 ~ "high",
                               ShadeTolerance > 1 ~ "low"))

gs$Tcnratio <- as.factor(gs$Tcnratio)
gs$Tallelopathy <- as.factor(gs$Tallelopathy)
gs$Tleafret <- as.factor(gs$Tleafret)
gs$Tleaflife <- as.factor(gs$Tleaflife)
gs$Tfertilityreq <- as.factor(gs$Tfertilityreq)
gs$Tdroughttol <- as.factor(gs$Tdroughttol)
gs$Tshadetol <- as.factor(gs$Tshadetol)


###
# Define forest classes
cond <- fread("data/raw/cond.csv", select = c('PLT_CN', 
                                              'FORTYPCD'))
cond <- data.table(cond)
cond <- with(cond, cond[order(PLT_CN, FORTYPCD),])
cond <- cond[!duplicated(cond$PLT_CN), ] #remove duplication condition information
gs <- merge(gs, cond, by.x="pcn1", by.y="PLT_CN", all.x=T, all.y=F)
gs$fortypsm <- trunc(gs$FORTYPCD/10)*10
gs$fortypsm <- as.integer(gs$fortypsm)
gs[(gs$fortypsm %in% c(510, 520)), "fortypsm"] <- 500
gs[gs$fortypsm %in% c(720), "fortypsm"] <- 700
gs[gs$fortypsm %in% c(930), "fortypsm"] <- 920
gs$fortypsm <- as.factor(gs$fortypsm)


# Save gs
saveRDS(gs, "output/gs_withgroups.RDS")




# Recruitment rate -------------------------------
rs <- readRDS("output/nci_formodels_recr.RDS") #made in growth_prep.R

spdata <- read.csv("/Users/Anika/Documents/GradSchool/FIA_EstimateProj/ACGCA_FIAv51_from_JWL/REF_SPECIES_jwl.csv")
spdata <- spdata %>% dplyr::select(SPCD, COMMON_NAME, GENUS, SPECIES, JENKINS_SPGRPCD) #get biomass parameters if needed

rs <- merge(rs, spdata, by.x="spcd", by.y="SPCD", all.x=T, all.y=F)


# Add whether recruited or not
rs$recir <- 0
rs[rs$FATE==1,]$recir <- 1
rs$recir <- rs$recir/rs$t #probability of recruitment per year

###
# Define evergreen vs deciduous
classten <- read.csv("data/raw/Compiled_decid0ever1.csv")
classten <- classten[, c("spnames", "AccSpeciesID", "Decid0Ever1")]

rs <- rs %>% tidyr::unite(gensp, GENUS, SPECIES, sep = " ", remove=F) #get col to merge classten on

rs <- merge(rs, classten, by.x="gensp", by.y="spnames", all.x=T)


###
# Define age class
rs <- rs %>%
  mutate(YOU0OLD1 = case_when(STDAGE < 60 ~ 0,
                              STDAGE >= 60 ~ 1))
rs$YOU0OLD1 <- as.factor(rs$YOU0OLD1)

###
# Define canopy class
rs <- rs %>%
  mutate(NON0CAN1 = case_when(cclcd2 < 4 ~ 1,
                              cclcd2 >= 4 ~ 0))
rs$NON0CAN1 <- as.factor(rs$NON0CAN1)

###
# Define N dep class
Ndep <- raster("/Users/Anika/Documents/GradSchool/FIA_EstimateProj/output/N_dep.grd")

coords <- rs[,c("LON","LAT")]
dep <- raster::extract(Ndep, coords)

rs <- cbind(rs, dep)

rs <- rs %>%
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
rs$NDEP <- as.factor(rs$NDEP)


###
# Define soil moisture class
us_soil <- raster("output/smap_proj.grd")

coords <- rs[,c("LON","LAT")]
sm <- raster::extract(us_soil, coords)

rs <- cbind(rs, sm)

rs <- rs %>%
  mutate(SM = case_when(sm > 0 & sm < 0.1 ~ 0,
                        sm >= 0.1 & sm < 0.2 ~ 0.1,
                        sm >= 0.2 & sm < 0.3 ~ 0.2,
                        sm >= 0.3 & sm < 0.4 ~ 0.3,
                        sm >= 0.4 & sm < 0.5 ~ 0.4,
                        sm >= 0.5 & sm < 0.6 ~ 0.5,
                        sm >= 0.6 & sm < 0.7 ~ 0.6,
                        sm >= 0.7 ~ 0.7))
rs$SM <- as.factor(rs$SM)


###
# Define soil texture classes, soilgrids
texture <- raster("/Users/Anika/Downloads/TEXMHT_M_sl4_250m_ll.tif")
texture <- raster::extract(texture, coords)
#require(GSIF)
#data(soil.legends)

rs <- cbind(rs, texture)
rs$texture <- as.factor(rs$texture)

###
# Define MR type
phil <- read.csv("data/raw/Phillips_AM_EM_data.csv")
phil <- phil[,3:5]
colnames(phil) <- c("GENUS", "SPECIES", "Mycostatus")
phil$spnames <- paste(phil$GENUS, phil$SPECIES)

splist <- rs[,c("GENUS","SPECIES")]
splist$spnames <- paste(splist$GENUS, splist$SPECIES)
splist <- splist[!duplicated(splist$spnames), ]

mrlist0 <- merge(splist, phil, by="spnames", all.x=T, all.y=F)
mrlist1 <- mrlist0[!is.na(mrlist0$Mycostatus), ]
mrlist1 <- mrlist1[,c(1:3,6)] #keep spnames, species, genus, mycostatus
colnames(mrlist1) <- c("spnames", "GENUS", "SPECIES", "Mycostatus")
mrlist2 <- mrlist1 %>% 
  group_by(GENUS) %>% 
  summarize(myco_phil = first(Mycostatus)) #MR for genus level (can merge back on to mrlist)

rs <- merge(rs, mrlist2, by="GENUS", all.x=T, all.y=F)



###
# Define trait classes
usda <- readRDS("data/derived/usda_numeric.RDS")
rs <- merge(rs, usda, by.x="spcd", by.y="SPCD", all.x=T, all.y=F)

rs <- rs %>%
  mutate(Tcnratio = case_when(CNRatio == 1 ~ "high",
                              STDAGE > 1 ~ "low"),
         Tallelopathy = case_when(Allelopath == 1 ~ "yes",
                                  Allelopath == 0 ~ "no"),
         Tleafret = case_when(LeafRetention == 1 ~ "yes",
                              LeafRetention == 0 ~ "no"),
         Tleaflife = case_when(Lifespan == 1 ~ "long",
                               Lifespan > 1 ~ "short"),
         Tfertilityreq = case_when(FerilityReq == 1 ~ "high",
                                   FerilityReq > 1 ~ "low"),
         Tdroughttol = case_when(DroughtTolerance <= 2 ~ "high",
                                 DroughtTolerance > 2 ~ "low"),
         Tshadetol = case_when(ShadeTolerance == 1 ~ "high",
                               ShadeTolerance > 1 ~ "low"))

rs$Tcnratio <- as.factor(rs$Tcnratio)
rs$Tallelopathy <- as.factor(rs$Tallelopathy)
rs$Tleafret <- as.factor(rs$Tleafret)
rs$Tleaflife <- as.factor(rs$Tleaflife)
rs$Tfertilityreq <- as.factor(rs$Tfertilityreq)
rs$Tdroughttol <- as.factor(rs$Tdroughttol)
rs$Tshadetol <- as.factor(rs$Tshadetol)


###
# Define forest classes
cond <- fread("data/raw/cond.csv", select = c('PLT_CN', 
                                              'FORTYPCD'))
cond <- data.table(cond)
cond <- with(cond, cond[order(PLT_CN, FORTYPCD),])
cond <- cond[!duplicated(cond$PLT_CN), ] #remove duplication condition information
rs <- merge(rs, cond, by.x="pcn2", by.y="PLT_CN", all.x=T, all.y=F)
rs$fortypsm <- trunc(rs$FORTYPCD/10)*10
rs$fortypsm <- as.integer(rs$fortypsm)
rs[(rs$fortypsm %in% c(510, 520)), "fortypsm"] <- 500
rs[rs$fortypsm %in% c(720), "fortypsm"] <- 700
rs[rs$fortypsm %in% c(930), "fortypsm"] <- 920
rs$fortypsm <- as.factor(rs$fortypsm)

# Save rs
saveRDS(rs, "output/rs_withgroups.RDS")




# Survival rate ----------------------------------
ms <- readRDS("output/nci_formodels_mort.RDS") #made in growth_prep.R

spdata <- read.csv("/Users/Anika/Documents/GradSchool/FIA_EstimateProj/ACGCA_FIAv51_from_JWL/REF_SPECIES_jwl.csv")
spdata <- spdata %>% dplyr::select(SPCD, COMMON_NAME, GENUS, SPECIES, JENKINS_SPGRPCD) #get biomass parameters if needed

ms <- merge(ms, spdata, by.x="spcd", by.y="SPCD", all.x=T, all.y=F)

source("scripts/functions.R")

# Add whether survived or not
ms$mortr <- 0
ms[ms$FATE==3,]$mortr <- 1
ms$mortr <- ms$mortr/ms$t #probability of mortality per year
ms$survr <- ms$surv/ms$t

###
# Define evergreen vs deciduous
classten <- read.csv("data/raw/Compiled_decid0ever1.csv")
classten <- classten[, c("spnames", "AccSpeciesID", "Decid0Ever1")]

ms <- ms %>% tidyr::unite(gensp, GENUS, SPECIES, sep = " ", remove=FALSE) #get col to merge classten on

ms <- merge(ms, classten, by.x="gensp", by.y="spnames", all.x=T)


###
# Define age class
ms <- ms %>%
  mutate(YOU0OLD1 = case_when(STDAGE < 60 ~ 0,
                              STDAGE >= 60 ~ 1))
ms$YOU0OLD1 <- as.factor(ms$YOU0OLD1)

###
# Define canopy class
ms <- ms %>%
  mutate(NON0CAN1 = case_when(CCLCD < 4 ~ 1,
                              CCLCD >= 4 ~ 0))
ms$NON0CAN1 <- as.factor(ms$NON0CAN1)

###
# Define N dep class
Ndep <- raster("/Users/Anika/Documents/GradSchool/FIA_EstimateProj/output/N_dep.grd")

coords <- ms[,c("LON","LAT")]
dep <- raster::extract(Ndep, coords)

ms <- cbind(ms, dep)

ms <- ms %>%
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
ms$NDEP <- as.factor(ms$NDEP)


###
# Define soil moisture class
us_soil <- raster("output/smap_proj.grd")

coords <- ms[,c("LON","LAT")]
sm <- raster::extract(us_soil, coords)

ms <- cbind(ms, sm)

ms <- ms %>%
  mutate(SM = case_when(sm > 0 & sm < 0.1 ~ 0,
                        sm >= 0.1 & sm < 0.2 ~ 0.1,
                        sm >= 0.2 & sm < 0.3 ~ 0.2,
                        sm >= 0.3 & sm < 0.4 ~ 0.3,
                        sm >= 0.4 & sm < 0.5 ~ 0.4,
                        sm >= 0.5 & sm < 0.6 ~ 0.5,
                        sm >= 0.6 & sm < 0.7 ~ 0.6,
                        sm >= 0.7 ~ 0.7))
ms$SM <- as.factor(ms$SM)


###
# Define soil texture classes, soilgrids
texture <- raster("/Users/Anika/Downloads/TEXMHT_M_sl4_250m_ll.tif")
texture <- raster::extract(texture, coords)
#require(GSIF)
#data(soil.legends)

ms <- cbind(ms, texture)
ms$texture <- as.factor(ms$texture)

###
# Define MR type
phil <- read.csv("data/raw/Phillips_AM_EM_data.csv")
phil <- phil[,3:5]
colnames(phil) <- c("GENUS", "SPECIES", "Mycostatus")
phil$spnames <- paste(phil$GENUS, phil$SPECIES)

splist <- ms[,c("GENUS","SPECIES")]
splist$spnames <- paste(splist$GENUS, splist$SPECIES)
splist <- splist[!duplicated(splist$spnames), ]

mrlist0 <- merge(splist, phil, by="spnames", all.x=T, all.y=F)
mrlist1 <- mrlist0[!is.na(mrlist0$Mycostatus), ]
mrlist1 <- mrlist1[,c(1:3,6)] #keep spnames, species, genus, mycostatus
colnames(mrlist1) <- c("spnames", "GENUS", "SPECIES", "Mycostatus")
mrlist2 <- mrlist1 %>% 
  group_by(GENUS) %>% 
  summarize(myco_phil = first(Mycostatus)) #MR for genus level (can merge back on to mrlist)

ms <- merge(ms, mrlist2, by="GENUS", all.x=T, all.y=F)



###
# Define trait classes
usda <- readRDS("data/derived/usda_numeric.RDS")
usda <- usda[,c(1:5,7:37)] #remove repeated cols like genus and common name
ms <- merge(ms, usda, by.x="spcd", by.y="SPCD", all.x=T, all.y=F)

ms <- ms %>%
  mutate(Tcnratio = case_when(CNRatio == 1 ~ "high",
                              CNRatio > 1 ~ "low"),
         Tallelopathy = case_when(Allelopath == 1 ~ "yes",
                                  Allelopath == 0 ~ "no"),
         Tleafret = case_when(LeafRetention == 1 ~ "yes",
                              LeafRetention == 0 ~ "no"),
         Tleaflife = case_when(Lifespan == 1 ~ "long",
                               Lifespan > 1 ~ "short"),
         Tfertilityreq = case_when(FerilityReq == 1 ~ "high",
                                   FerilityReq > 1 ~ "low"),
         Tdroughttol = case_when(DroughtTolerance <= 2 ~ "high",
                                 DroughtTolerance > 2 ~ "low"),
         Tshadetol = case_when(ShadeTolerance == 1 ~ "high",
                               ShadeTolerance > 1 ~ "low"))

ms$Tcnratio <- as.factor(ms$Tcnratio)
ms$Tallelopathy <- as.factor(ms$Tallelopathy)
ms$Tleafret <- as.factor(ms$Tleafret)
ms$Tleaflife <- as.factor(ms$Tleaflife)
ms$Tfertilityreq <- as.factor(ms$Tfertilityreq)
ms$Tdroughttol <- as.factor(ms$Tdroughttol)
ms$Tshadetol <- as.factor(ms$Tshadetol)


###
# Define forest classes
cond <- fread("data/raw/cond.csv", select = c('PLT_CN', 
                                              'FORTYPCD'))
cond <- data.table(cond)
cond <- with(cond, cond[order(PLT_CN, FORTYPCD),])
cond <- cond[!duplicated(cond$PLT_CN), ] #remove duplication condition information
ms <- merge(ms, cond, by.x="pcn1", by.y="PLT_CN", all.x=T, all.y=F)
ms$fortypsm <- trunc(ms$FORTYPCD/10)*10
ms$fortypsm <- as.integer(ms$fortypsm)
ms[(ms$fortypsm %in% c(510, 520)), "fortypsm"] <- 500
ms[ms$fortypsm %in% c(720), "fortypsm"] <- 700
ms[ms$fortypsm %in% c(930), "fortypsm"] <- 920
ms$fortypsm <- as.factor(ms$fortypsm)

###
# Assign names to groups (canopy, fixer, soil texture, forest type, age, decid)
ms <- ms %>%
  mutate(canopy = case_when(NON0CAN1 == 0 ~ "noncanopy",
                              NON0CAN1 == 1 ~ "canopy"),
         fixer = case_when(FIX == 0 ~ "nonfixer",
                           FIX == 1 ~ "fixer"),
         soiltexture = case_when(texture == 1 ~ "clay",
                                 texture == 2 ~ "silty clay",
                                 texture == 3 ~ "sandy clay",
                                 texture == 4 ~ "clay loam",
                                 texture == 5 ~ "silty clay loam",
                                 texture == 6 ~ "sandy clay loam",
                                 texture == 7 ~ "loam",
                                 texture == 8 ~ "silty loam",
                                 texture == 9 ~ "sandy loam",
                                 texture == 10 ~ "silt",
                                 texture == 11 ~ "loamy sand",
                                 texture == 12 ~ "sand"),
         fortype = case_when(fortypsm == 100 ~ "white/red/jack pine",
                             fortypsm == 120 ~ "spruce/fir",
                             fortypsm == 140 ~ "longleaf/slash pine",
                             fortypsm == 150 ~ "tropical softwood",
                             fortypsm == 160 ~ "loblolly/shortleaf pine",
                             fortypsm == 170 ~ "other eastern softwoods",
                             fortypsm == 180 ~ "pinyon/juniper",
                             fortypsm == 200 ~ "douglas-fir",
                             fortypsm == 220 ~ "ponderosa pine",
                             fortypsm == 240 ~ "western white pine",
                             fortypsm == 260 ~ "fir/spruce/mountain hemlock",
                             fortypsm == 280 ~ "lodgepole pine",
                             fortypsm == 300 ~ "hemlock/sitka spruce",
                             fortypsm == 320 ~ "western larch",
                             fortypsm == 340 ~ "redwood",
                             fortypsm == 360 ~ "other western softwoods",
                             fortypsm == 370 ~ "CA mixed conifer",
                             fortypsm == 380 ~ "exotic softwood",
                             fortypsm == 390 ~ "other softwood",
                             fortypsm == 400 ~ "oak/pine",
                             fortypsm == 500 ~ "oak/hickory",
                             fortypsm == 600 ~ "oak/gum/cypress",
                             fortypsm == 700 ~ "elm/ash/cottonwood",
                             fortypsm == 800 ~ "maple/beech/birch",
                             fortypsm == 900 ~ "aspen/birch",
                             fortypsm == 910 ~ "alder/maple",
                             fortypsm == 920 ~ "western oak",
                             fortypsm == 940 ~ "tanoak/laurel"),
         age = case_when(YOU0OLD1 == 0 ~ "young",
                         YOU0OLD1 == 1 ~ "old"),
         decid = case_when(Decid0Ever1 == 0 ~ "deciduous",
                           Decid0Ever1 == 0.5 ~ "semi-deciduous",
                           Decid0Ever1 == 1 ~ "evergreen")
  )


# Save gs
saveRDS(ms, "output/ms_withgroups.RDS")


