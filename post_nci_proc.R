
###########################
### post NCI processing ###
###########################


require(data.table)
require(dplyr)


options(scipen=6)
setwd("/Users/Anika/Documents/GradSchool/FIA_CompetitionModel/")


# Calc growth, survival, and recruitment --------------------------------
# LOAD DATA
FIA <- readRDS("output/nci_fromHab.RDS") #from nci in habanero


# GROWTH RATE
FIA$lg <- log(FIA$dbh2) - log(FIA$dbh1)
FIA$gr <- (FIA$dbh2 - FIA$dbh1)/(FIA$mt2 - FIA$mt1)

# Initialzize positivized relative growth rates
FIA$lgp <- log(FIA$dbh2) - log(FIA$dbh1)

#Compute "positivized" relative growth rates following approach in Condit et al. 2006 Science (Eqn.S5)
MDL <- 0.05
hMDL <- MDL/2
FIA[FIA$lg <=0 & !is.na(FIA$lg), ]$lgp <- 
  log(FIA[FIA$lg <= 0 & !is.na(FIA$lg), ]$dbh1 + hMDL) -
  log(FIA[FIA$lg <= 0 & !is.na(FIA$lg), ]$dbh1)

FIA$LGR <- FIA$lg/(FIA$mt2 - FIA$mt1)
FIA$LGRp <- FIA$lgp/(FIA$mt2 - FIA$mt1)
FIA$t <- FIA$mt2 - FIA$mt1

# remove trees where t is NA or t is 0
FIA <- FIA[FIA$t > 0 & !is.na(FIA$t), ]

# check growth rates make sense
summary(FIA[FIA$FATE==2, ]$LGRp)
hist(FIA[FIA$FATE==2, ]$LGRp, breaks=50)

# check outliers
nrow(FIA[which(FIA$FATE==2 & FIA$LGRp > 0.5),]) #29

# RECRUITMENT (individual)
FIA <- data.table(FIA)
FIA[, NEWtrees := NULL]
FIA[, ORIGtrees := NULL]
FIA[, FINtrees := NULL]
FIA[, lambda := NULL]
FIA[, surv := NULL]
FIA[, recr := NULL]

FIA[FATE==1, reci := 1]
FIA[FATE==2 | FATE==3, reci := 0]

# SURVIVAL (individual)
FIA[FATE==1 | FATE==2, survi := 1]
FIA[FATE==3, survi := 0]

cond <- fread("data/raw/cond.csv", select = c('PLT_CN',
                                              'STDAGE'))
cond$PLT_CN <- as.numeric(cond$PLT_CN)
cond <- with(cond, cond[order(PLT_CN, STDAGE),])
cond <- cond[!duplicated(cond$PLT_CN), ] #remove duplication condition information
FIA <- merge(FIA, cond, by.x="pcn1", by.y="PLT_CN", all.x=T, all.y=F)

saveRDS(FIA, "output/indlevel_data.RDS")

# demographic rates (plot scale) ----------------------------------------
FIA <- readRDS("output/nci_fromHab.RDS") #from nci in habanero

cond <- fread("data/raw/cond.csv", select = c('PLT_CN', 
                                              'DSTRBCD1',
                                              'CARBON_LITTER',
                                              'CARBON_SOIL_ORG',
                                              'CARBON_UNDERSTORY_AG',
                                              'STDAGE'))
cond$PLT_CN <- as.numeric(cond$PLT_CN)
cond <- with(cond, cond[order(PLT_CN, STDAGE),])
cond <- cond[!duplicated(cond$PLT_CN), ] #remove duplication condition information

FIA <- merge(FIA, cond, by.x="pcn1", by.y="PLT_CN", all.x=T, all.y=F)
FIA <- rename(FIA, DSTRBCD_TYPE1_1 = DSTRBCD1,
              CARBON_SOIL_ORG1 = CARBON_SOIL_ORG,
              CARBON_LITTER1 = CARBON_LITTER,
              CARBON_UNDERSTORY_AG1 = CARBON_UNDERSTORY_AG,
              STDAGE1 = STDAGE)# rename new cols
FIA <- merge(FIA, cond, by.x="pcn2", by.y="PLT_CN", all.x=T, all.y=F)
FIA <- rename(FIA, DSTRBCD_TYPE1_2 = DSTRBCD1,
              CARBON_SOIL_ORG2 = CARBON_SOIL_ORG,
              CARBON_LITTER2 = CARBON_LITTER,
              CARBON_UNDERSTORY_AG2 = CARBON_UNDERSTORY_AG,
              STDAGE2 = STDAGE)# rename new cols

FIA$BAm2ha1 <- FIA$dbh1/10*FIA$tpha1
FIA$BAm2ha2 <- FIA$dbh2/10*FIA$tpha2

FIA$STDAGE <- rowMeans(FIA[,c("STDAGE1","STDAGE2")], na.rm=T)

saveRDS(FIA, "output/plotlevel_data.RDS")

