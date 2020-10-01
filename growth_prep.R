

######################
# 
# growth_prep.R
#
#
# prepare growth data for models_ind.R
# produces "output/nci_formodels.RDS"
#
#
# Anika Petach
# 12/13/18
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
#FIA <- merge(FIA, plot, by.x="PCN1", by.y="CN", all.x=T, all.y=F)
FIA <- merge(FIA, plot, by.x="pcn1", by.y="CN", all.x=T, all.y=F)
FIA <- data.table(FIA)
FIA <- FIA[STATECD != 15 & STATECD < 60,] # remove HI and islands (PR and VI)

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
#temp <- FIA[DISTm > 7.3152]
#plotsout <- unique(temp$PCN1)
#FIA <- FIA[!PCN1 %in% plotsout]

# Growth analysis (ind) -----------------------------------------------------
growths <- FIA[FATE==2,]

# Average values from two points
g <- growths %>%
  mutate(dbh = (dbh1 + dbh2)/2, 
         NCI = (NCI_total1 + NCI_total2)/2,
         NCI_prop = (NCI_prop1 + NCI_prop2)/2)
         #CARBON_LITTER = (CARBON_LITTER1 + CARBON_LITTER2)/2,
         #CARBON_SOIL_ORG = (CARBON_SOIL_ORG1 + CARBON_SOIL_ORG2)/2,
         #CARBON_UNDER_AG = (CARBON_UNDERSTORY_AG1 + CARBON_UNDERSTORY_AG2)/2)

rm(growths)

# fix Genus NA
levels(g$Genus) <- c(levels(g$Genus), "NA")
g[is.na(g$Genus), ]$Genus <- "NA"

# transform and scale covariates
# gs <- g %>%
#   mutate(NCI_props = as.vector(scale(asin(sqrt(NCI_prop)))),
#          dbhs = as.vector(scale(log(dbh))),
#          NCIs = as.vector(scale(log(NCI + 0.001))))
gs <- g %>%
  mutate(NCI_props = as.vector(scale(NCI_prop, center=T, scale=T)),
         dbhs = as.vector(scale(dbh, center=T, scale=T)),
         NCIs = as.vector(scale(NCI_total2, center=T, scale=T)))

rm(g)

counts <- gs %>% group_by(Genus) %>% summarise(n=n_distinct(tcn1))
print(counts, n=nrow(out))

# gs <- dplyr::select(gs, PCN1, PCN2, SPCD, TCN1, TCN2, MEASYEAR2, WATERCD2, MEASYEAR1, WATERCD1, SUBP1,
#        SUBP2, CCLCD1, CCLCD2, CARBON_AG2, CARBON_BG2, CARBON_AG1, CARBON_BG2, FATE, LAT,
#        LON, AZIMUTH, TPA_UNADJ, SUBP1, CARBON_LITTER, CARBON_SOIL_ORG, CARBON_UNDER_AG, 
#        DSTRBCD_TYPE1_1, DSTRBCD_TYPE1_2, SLOPE, ASPECT, STDAGE, DISTm, 
#        ELEVm, Genus, COMMON_NAME, FIX, ACT1vsRIZ0, area, area2, LGRp, t, reci, survi, 
#        BAm2ha1, BAm2ha2, STATECD, NCI_total1, NCI_total2, dbh, NCI, NCI_prop, NCI_props,
#        dbhs, NCIs)

gs <- dplyr::select(gs, -spcd1, -spcd2, -tpha1, -tpha2, -lg)

saveRDS(gs, "output/nci_formodels_growth.RDS")


# investigate alder density -------------------------
temp <- FIA[FIA$Genus == "Alnus",] #used FIA from first lines of script
temp1 <- temp %>% dplyr::group_by(pcn1) %>% summarize(density = sum(tpha2,na.rm=T))
hist(temp1$density, breaks=30)
abline(v=400, col='red') #fang 2019 value up to which doug fir were facilitated
abline(v=100, col='blue') #fang 2019 value up to which red cedar were facilitated
mean(temp1$density)