



######################
# 
# usda_prep.R
#
#
# prepare usda data for use in models_ind.R
# produces "data/derived/usda.RDS"
#
#
# Anika Petach
# 1/10/19
#
######################

require(dplyr)

setwd("/Users/Anika/Documents/GradSchool/FIA_CompetitionModel/")
usdaraw <- read.csv("data/raw/USDAplants_raw.csv", header=T, na.strings=c("","NA"))


usda1 <- usdaraw[!is.na(usdaraw$Active.Growth.Period) | !is.na(usdaraw$C.N.Ratio),]

nrow(usdaraw)
nrow(usda1)

# merge with jwl to get spcd numbers
spdata <- read.csv("/Users/Anika/Documents/GradSchool/FIA_EstimateProj/ACGCA_FIAv51_from_JWL/REF_SPECIES_jwl.csv")
spdata <- spdata %>% dplyr::select(SPCD, COMMON_NAME, GENUS, SPECIES, SPECIES_SYMBOL) #get biomass parameters if needed

usda <- merge(usda1, spdata, by.x="Accepted.Symbol", by.y="SPECIES_SYMBOL", all.x=T, all.y=F)

# rename columns
colnames(usda) <- c("Symbol","AltSymbol","ScientificName","CommonName","Area","Genus","Family",
                    "FamilySymbol","Duration","ActiveGrowth","CNRatio","FireResistance",
                    "FoliagePorosity","GrowthRate","Height20yrsFt","HeightMatureFt",
                    "Allelopath","LeafRetention","Lifespan","Resprout","AnaerobicTolerance",
                    "CaCO3Tolerance","DroughtTolerance","FerilityReq","MoistureUse",
                    "pHmin","pHmax","PlantingDensityMinPerAcre","PlantingDensityMaxPerAcre",
                    "MaxRootDepthIn","SalinityTolerance","ShadeTolerance","BloomPeriod",
                    "SeedAbundance","SeedsPerPound","VegSpreadRate","SPCD","COMMON_NAME",
                    "GENUS","SPECIES")

saveRDS(usda, "data/derived/usda.RDS")

# convert factors to numbers
levels(usda$CNRatio)[1] <- 1 #High
levels(usda$CNRatio)[2] <- 3 #Low
levels(usda$CNRatio)[3] <- 2 #Medium
usda$CNRatio <- as.numeric(as.character(usda$CNRatio))

levels(usda$FireResistance)[1] <- 0 #No
levels(usda$FireResistance)[2] <- 1 #Yes
usda$FireResistance <- as.numeric(as.character(usda$FireResistance))

levels(usda$FoliagePorosity)[1] <- 3 #Dense
levels(usda$FoliagePorosity)[2] <- 2 #Moderate
levels(usda$FoliagePorosity)[3] <- 1 #Porous
usda$FoliagePorosity <- as.numeric(as.character(usda$FoliagePorosity))

levels(usda$GrowthRate)[1] <- 2 #Moderate
levels(usda$GrowthRate)[2] <- 1 #Rapid
levels(usda$GrowthRate)[3] <- 2 #Slow
usda$GrowthRate <- as.numeric(as.character(usda$GrowthRate))

levels(usda$Allelopath)[1] <- 0 #No
levels(usda$Allelopath)[2] <- 1 #Yes
usda$Allelopath <- as.numeric(as.character(usda$Allelopath))

levels(usda$LeafRetention)[1] <- 0 #No
levels(usda$LeafRetention)[2] <- 1 #Yes
usda$LeafRetention <- as.numeric(as.character(usda$LeafRetention))

levels(usda$Lifespan)[1] <- 1 #Long
levels(usda$Lifespan)[2] <- 2 #Moderate
levels(usda$Lifespan)[3] <- 3 #Short
usda$Lifespan <- as.numeric(as.character(usda$Lifespan))

levels(usda$Resprout)[1] <- 0 #No
levels(usda$Resprout)[2] <- 1 #Yes
usda$Resprout <- as.numeric(as.character(usda$Resprout))

levels(usda$AnaerobicTolerance)[1] <- 1 #High
levels(usda$AnaerobicTolerance)[2] <- 3 #Low
levels(usda$AnaerobicTolerance)[3] <- 2 #Medium
levels(usda$AnaerobicTolerance)[4] <- 4 #None
usda$AnaerobicTolerance <- as.numeric(as.character(usda$AnaerobicTolerance))

levels(usda$CaCO3Tolerance)[1] <- 1 #High
levels(usda$CaCO3Tolerance)[2] <- 3 #Low
levels(usda$CaCO3Tolerance)[3] <- 2 #Medium
levels(usda$CaCO3Tolerance)[4] <- 4 #None
usda$CaCO3Tolerance <- as.numeric(as.character(usda$CaCO3Tolerance))

levels(usda$DroughtTolerance)[1] <- 1 #High
levels(usda$DroughtTolerance)[2] <- 3 #Low
levels(usda$DroughtTolerance)[3] <- 2 #Medium
levels(usda$DroughtTolerance)[4] <- 4 #None
usda$DroughtTolerance <- as.numeric(as.character(usda$DroughtTolerance))

levels(usda$FerilityReq)[1] <- 1 #High
levels(usda$FerilityReq)[2] <- 3 #Low
levels(usda$FerilityReq)[3] <- 2 #Medium
usda$FerilityReq <- as.numeric(as.character(usda$FerilityReq))

levels(usda$MoistureUse)[1] <- 1 #High
levels(usda$MoistureUse)[2] <- 3 #Low
levels(usda$MoistureUse)[3] <- 2 #Medium
usda$MoistureUse <- as.numeric(as.character(usda$MoistureUse))

levels(usda$SalinityTolerance)[1] <- 1 #High
levels(usda$SalinityTolerance)[2] <- 3 #Low
levels(usda$SalinityTolerance)[3] <- 2 #Medium
levels(usda$SalinityTolerance)[4] <- 4 #None
usda$SalinityTolerance <- as.numeric(as.character(usda$SalinityTolerance))

levels(usda$ShadeTolerance)[1] <- 2 #Intermediate
levels(usda$ShadeTolerance)[2] <- 3 #Intolerant
levels(usda$ShadeTolerance)[3] <- 1 #Tolerant
usda$ShadeTolerance <- as.numeric(as.character(usda$ShadeTolerance))

levels(usda$BloomPeriod)[1] <- 3 #Early spring
levels(usda$BloomPeriod)[2] <- 7 #Early summer
levels(usda$BloomPeriod)[3] <- 11 #Fall
levels(usda$BloomPeriod)[4] <- 12 #Indeterminate
levels(usda$BloomPeriod)[5] <- 6 #Late spring
levels(usda$BloomPeriod)[6] <- 10 #Late summer
levels(usda$BloomPeriod)[7] <- 2 #Late winter
levels(usda$BloomPeriod)[8] <- 5 #Mid spring
levels(usda$BloomPeriod)[9] <- 9 #Mid summer
levels(usda$BloomPeriod)[10] <- 4 #Spring
levels(usda$BloomPeriod)[11] <- 8 #Summer
levels(usda$BloomPeriod)[12] <- 1 #Winter
usda$BloomPeriod <- as.numeric(as.character(usda$BloomPeriod))

levels(usda$SeedAbundance)[1] <- 1 #High
levels(usda$SeedAbundance)[2] <- 3 #Low
levels(usda$SeedAbundance)[3] <- 2 #Medium
levels(usda$SeedAbundance)[4] <- 4 #None
usda$SeedAbundance <- as.numeric(as.character(usda$SeedAbundance))

levels(usda$VegSpreadRate)[1] <- 2 #Moderate
levels(usda$VegSpreadRate)[2] <- 4 #None
levels(usda$VegSpreadRate)[3] <- 1 #Rapid
levels(usda$VegSpreadRate)[4] <- 3 #Slow
usda$VegSpreadRate <- as.numeric(as.character(usda$VegSpreadRate))

levels(usda$ActiveGrowth)[1] <- 9 #Fall
levels(usda$ActiveGrowth)[2] <- 3 #Fall, winter, spring
levels(usda$ActiveGrowth)[3] <- 7 #Spring
levels(usda$ActiveGrowth)[4] <- 6 #Spring, fall
levels(usda$ActiveGrowth)[5] <- 4 #Spring, summer
levels(usda$ActiveGrowth)[6] <- 2 #Spring, summer, fall
levels(usda$ActiveGrowth)[7] <- 8 #Summer
levels(usda$ActiveGrowth)[8] <- 5 #Summer, fall
levels(usda$ActiveGrowth)[9] <- 1 #Year round
usda$ActiveGrowth <- as.numeric(as.character(usda$ActiveGrowth))

saveRDS(usda, "data/derived/usda_numeric.RDS")
