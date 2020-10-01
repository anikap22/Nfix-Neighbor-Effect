




######################
# 
# mort_prep.R
#
#
# prepare mortality data for models_ind.R
# produces "output/nci_formodels.RDS"
#
#
# Anika Petach
# 12/17/18
#
######################

require(data.table)
require(dplyr)

options(scipen=6)

setwd("/Users/Anika/Documents/GradSchool/FIA_CompetitionModel/")
FIA <- readRDS("output/indlevel_data.RDS") #import data table

# Remove islands -----------------------------------------------------------
plot <- fread("data/raw/plot.csv", select = c('CN', 
                                              'STATECD'))
plot$CN <- as.numeric(plot$CN)
FIA <- merge(FIA, plot, by.x="pcn1", by.y="CN", all.x=T, all.y=F)
FIA <- FIA[STATECD != 15 & STATECD < 60,] #remove islands

# NCI_prop and NCI_total ---------------------------------------------------
# NCI total
FIA[, NCI_total1 := NCI_Nfixer1 + NCI_nonfixer1]
FIA[, NCI_total2 := NCI_Nfixer2 + NCI_nonfixer2]

# NCI prop
FIA[, NCI_prop1 := NCI_Nfixer1 / NCI_total1]
FIA[, NCI_prop2 := NCI_Nfixer2 / NCI_total2]
# When NCI_total1 or NCI_total2 is 0 then NCI_prop should be 0 not NA
FIA[NCI_total1==0, NCI_prop1 := 0]
FIA[NCI_total2==0, NCI_prop2 := 0]

# remove DISTm > 7.3152 #should remove after you filter DESIGNCD
# remove all trees from plots where any tree is over 7.3152m (was different design code)
#FIA <- FIA[DISTm < 7.3152]
# temp <- FIA[DISTm > 7.3152]
# plotsout <- unique(temp$PCN1)
# FIA <- FIA[!PCN1 %in% plotsout]

#remove diacm1 = NA (no data for these trees)
#FIA <- FIA[!is.na(FIA$DIAcm1), ]

# get whichever CCLCD is not NA
FIA <- FIA %>% mutate(CCLCD = coalesce(cclcd1, cclcd2))

# mort <- dplyr::select(FIA, PCN1, PCN2, TCN1, MEASYEAR1, SPCD, DIAcm1, TPA_UNADJ, Genus,
#                       LAT, LON,
#                       FIX, ACT1vsRIZ0, NCI_weighted_Nfixer, NCI_weighted_nonfixer, ELEVm, 
#                       NCI_total1, NCI_prop1, CCLCD, STDAGE, SLOPE, ASPECT, CARBON_LITTER1, 
#                       CARBON_SOIL_ORG1, CARBON_UNDERSTORY_AG1, DSTRBCD_TYPE1_1, FATE, t,
#                       SUBP1)

#xtabs(~FATE + FIX, data=mort)

# transform and scale covariates
# ms <- mort %>%
#   mutate(NCI_props = as.vector(scale(asin(sqrt(NCI_prop1)))),
#          dbhs = as.vector(scale(log(DIAcm1))),
#          NCIs = as.vector(scale(log(NCI_total1 + 0.001))))
ms <- FIA %>%
  mutate(NCI_props = as.vector(scale(NCI_prop1, center=T, scale=T)),
         dbhs = as.vector(scale(dbh1, center=T, scale=T)),
         NCIs = as.vector(scale(NCI_total1, center=T, scale=T)))

ms <- dplyr::select(ms, -spcd1, -spcd2, -tpha1, -tpha2)


saveRDS(ms, "output/nci_formodels_mort.RDS")
