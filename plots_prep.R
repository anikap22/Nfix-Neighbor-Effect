
######################
# 
# growth_prep.R
#
#
# prepare data for models_plot.R
# produces "output/nci_formodels_plot.RDS"
#
#
# Anika Petach
# 1/17/19
#
######################

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

source("scripts/functions.R")

setwd("/Users/Anika/Documents/GradSchool/FIA_CompetitionModel/")

options(scipen = 6)


# Make plot time series: --------------------------------------------

FIA <- readRDS("output/plotlevel_data.RDS")

#get columns that are for fixers and nonfixers lumped
FIAp <- FIA %>%
  group_by(pcn1) %>%
  summarize(BA1 = sum(dbh1/10*tpha1, na.rm=T),
            BA2 = sum(dbh2/10*tpha2, na.rm=T),
            BA_Nfixer1 = sum(BAm2ha1[FIX==1], na.rm=T),
            BA_Nfixer2 = sum(BAm2ha2[FIX==1], na.rm=T),
            BA_nonfixer1 = sum(BAm2ha1[FIX==0], na.rm=T),
            BA_nonfixer2 = sum(BAm2ha2[FIX==0], na.rm=T),
            n = n(),
            #recp = length(LAT[FATE==1])/(length(LAT[FATE==2|FATE==3])*first(t)),
            #survp = 1-(log(length(LAT[FATE==2|FATE==3]))-log(length(LAT[FATE==1])))/first(t),
            stdage = mean(STDAGE, na.rm=T),
            #  slope = mean(SLOPE, na.rm=T),
            #  aspect = mean(ASPECT, na.rm=T),
            elev = mean(ELEVm, na.rm=T),
            dstrb1 = mean(DSTRBCD_TYPE1_1, na.rm=T),
            dstrb2 = mean(DSTRBCD_TYPE1_2, na.rm=T),
            soilC1 = mean(CARBON_SOIL_ORG1, na.rm=T),
            soilC2 = mean(CARBON_SOIL_ORG2, na.rm=T),
            understoryC1 = mean(CARBON_UNDERSTORY_AG1, na.rm=T),
            understoryC2 = mean(CARBON_UNDERSTORY_AG2, na.rm=T),
            litterC1 = mean(CARBON_LITTER1, na.rm=T),
            litterC2 = mean(CARBON_LITTER2, na.rm=T),
            ndied = sum(FATE[FATE==3])/3,
            nrec = sum(FATE[FATE==1]),
            mt1 = mean(mt1, na.rm=T),
            mt2 = mean(mt2, na.rm=T))

# Get plots with stdage:
# fia_plot <- fia[!is.na(fia$STDAGE),]
# 
# #unique_pcn1 <- unique(fia_plot$PCN1) #list of unique pcn numbers
# unique_combos <- unique(fia_plot[,c('PCN1','PCN2')]) #list of unique pcn combos
# 
# FIA_plot <- data.frame(unique_pcn1 = integer64(),
#                        BA_nonfixer1 = numeric(),
#                        BA_Nfixer1 = numeric(),
#                        BA_nonfixer2 = numeric(),
#                        BA_Nfixer2 = numeric(),
#                        mt1 = integer(),
#                        mt2 = integer(),
#                        pcn2 = integer64(),
#                        STDAGE = integer(),
#                        died = integer(),
#                        rec = integer(),
#                        num0 = integer(),
#                        fracD = numeric())
# 
# FIA_subplot <- data.frame(unique_pcn1 = integer64(),
#                           supblot = integer(),
#                           BA_nonfixer1 = numeric(),
#                           BA_Nfixer1 = numeric(),
#                           BA_nonfixer2 = numeric(),
#                           BA_Nfixer2 = numeric(),
#                           mt1 = integer(),
#                           mt2 = integer(),
#                           pcn2 = integer64(),
#                           STDAGE = integer(),
#                           died = integer(),
#                           rec = integer(),
#                           num0 = integer(),
#                           fracD = numeric())
# 
# 
# for(i in 1:nrow(unique_combos)){
#   print(i/nrow(unique_combos)*100)
#   #Plot:
#   plot_df <- fia[fia$PCN1 == unique_combos$PCN1[i] & fia$PCN2 == unique_combos$PCN2[i], ]
#   if(nrow(plot_df) > 1 & 
#      sd(plot_df$MEASYEAR1, na.rm=T) != 0 | 
#      is.na(sd(plot_df$MEASYEAR1, na.rm=T))){warning("mt")}
#   #if(nrow(plot_df)>1&sd(plot_df$prev_plt_cn1)!=0){warning("prev_plt_cn")}
#   FIA_plot[i,] <- list(unique_combos$PCN1[i],
#                        sum(plot_df[plot_df$FIX==0, "DIAcm1"]^2*pi/4, na.rm=T),
#                        sum(plot_df[plot_df$FIX==1, "DIAcm1"]^2*pi/4, na.rm=T),
#                        sum(plot_df[plot_df$FIX==0, "DIAcm2"]^2*pi/4, na.rm=T),
#                        sum(plot_df[plot_df$FIX==1, "DIAcm2"]^2*pi/4, na.rm=T),
#                        mean(plot_df$MEASYEAR1, na.rm=T),
#                        mean(plot_df$MEASYEAR2, na.rm=T),
#                        mean(plot_df$PCN2, na.rm=T),
#                        mean(plot_df$STDAGE, na.rm=T),
#                        nrow(plot_df[plot_df$FATE==3,]),
#                        nrow(plot_df[plot_df$FATE==1,]),
#                        nrow(plot_df),
#                        nrow(plot_df[plot_df$DEC0EVER1 == 0,])/nrow(plot_df))
#   
#   #Subplot:
#   subplot_list <- unique(plot_df$SUBP1)
#   for(j in 1:length(subplot_list)){
#     subplot_df <- plot_df[plot_df$SUBP1 == subplot_list[j], ]
#     FIA_subplot[j, ] <- list(unique_combos$PCN1[i],
#                              subplot_list[j],
#                              sum(subplot_df[subplot_df$FIX==0, "DIAcm1"]^2*pi/4, na.rm=T),
#                              sum(subplot_df[subplot_df$FIX==1, "DIAcm1"]^2*pi/4, na.rm=T),
#                              mean(subplot_df$MEASYEAR1, na.rm=T),
#                              mean(subplot_df$PCN2, na.rm=T),
#                              mean(subplot_df$STDAGE, na.rm=T),
#                              nrow(subplot_df[subplot_df$FATE==3,]),
#                              nrow(subplot_df[subplot_df$FATE==1,]),
#                              nrow(subplot_df),
#                              nrow(subplot_df[subplot_df$DEC0EVER1 == 0,])/nrow(subplot_df))
#     
#   }}

saveRDS(FIAp, "output/fia_plot.RDS")
#saveRDS(FIA_subplot, "output/fia_subplot.RDS")
# Calculate totals

# Read data ---------------------------------------------------------------
FIA_plot <- readRDS("output/fia_plot.RDS")
#FIA_subplot <- readRDS("output/fia_subplot.RDS")


# colnames(FIA_plot) <- c("pcn","BA_nonfixer1","BA_Nfixer1","BA_nonfixer2",
#                         "BA_Nfixer2","mt1","mt2","prev_plt_cn","stdage","died","rec",
#                         "num0","fracD")
FIA_plot$BA_total1 <- FIA_plot$BA_nonfixer1 + FIA_plot$BA_Nfixer1
FIA_plot$BA_Nfixer_prop1 <- FIA_plot$BA_Nfixer1/FIA_plot$BA_total1
FIA_plot$BA_total2 <- FIA_plot$BA_nonfixer2 + FIA_plot$BA_Nfixer2
FIA_plot$BA_Nfixer_prop2 <- FIA_plot$BA_Nfixer2/FIA_plot$BA_total2

# colnames(FIA_subplot) <- c("pcn","BA_nonfixer1","BA_Nfixer1","BA_nonfixer2",
#                            "BA_Nfixer2","mt1","mt2","prev_plt_cn","stdage","died","rec",
#                            "num0","fracD")
# FIA_subplot$BA_total1 <- FIA_subplot$BA_nonfixer1 + FIA_subplot$BA_Nfixer1
# FIA_subplot$BA_Nfixer_prop1 <- FIA_subplot$BA_Nfixer1/FIA_subplot$BA_total1
# FIA_subplot$BA_total2 <- FIA_subplot$BA_nonfixer2 + FIA_subplot$BA_Nfixer2
# FIA_subplot$BA_Nfixer_prop2 <- FIA_subplot$BA_Nfixer2/FIA_subplot$BA_total2

# BA change -----------------------------------
FIA_plot_raw <- FIA_plot #make a copy of the original data
# Plot growth
FIA_plot <- FIA_plot[!is.na(FIA_plot$mt1) & !is.na(FIA_plot$mt2),] #no plots excluded

FIA_plot$BA_change <- (FIA_plot$BA_total1 - FIA_plot$BA_total2) / 
  (FIA_plot$mt1 - FIA_plot$mt2) / FIA_plot$BA_total2 * 100
FIA_plot[FIA_plot$BA_total2 == 0, ]$BA_change <- 0

FIA_plot$BA_change_nonfixer <- (FIA_plot$BA_nonfixer1 - FIA_plot$BA_nonfixer2) /
  (FIA_plot$mt1 - FIA_plot$mt2) / FIA_plot$BA_nonfixer2 * 100
FIA_plot[FIA_plot$BA_nonfixer2 == 0, ]$BA_change_nonfixer <- 0

FIA_plot$BA_change_Nfixer <- (FIA_plot$BA_Nfixer1 - FIA_plot$BA_Nfixer2) /
  (FIA_plot$mt1 - FIA_plot$mt2) / FIA_plot$BA_Nfixer2 * 100
FIA_plot[FIA_plot$BA_Nfixer2 == 0, ]$BA_change_Nfixer <- 0

FIA_plot$BA_Nfixer_prop_total <- (FIA_plot$BA_Nfixer1 + FIA_plot$BA_Nfixer2) / 
  (FIA_plot$BA_total1 + FIA_plot$BA_total2)

# Subplot growth
# FIA_subplot <- FIA_subplot[!is.na(FIA_subplot$mt1) & !is.na(FIA_subplot$mt2),]
# 
# FIA_subplot$BA_change <- (FIA_subplot$BA_total1 - FIA_subplot$BA_total2) /
#   (FIA_subplot$mt1 - FIA_subplot$mt2) / FIA_subplot$BA_total2 * 100
# 
# FIA_subplot$BA_change_nonfixer <- (FIA_subplot$BA_nonfixer1 - FIA_subplot$BA_nonfixer2) /
#   (FIA_subplot$mt1 - FIA_subplot$mt2) / FIA_subplot$BA_nonfixer2 * 100
# 
# FIA_subplot$BA_change_Nfixer <- (FIA_subplot$BA_Nfixer1 - FIA_subplot$BA_Nfixer2) /
#   (FIA_subplot$mt1 - FIA_subplot$mt2) / FIA_subplot$BA_Nfixer2 * 100
# 
# FIA_subplot$BA_Nfixer_prop_total <- (FIA_subplot$BA_Nfixer1 + FIA_subplot$BA_Nfixer2) /
#   (FIA_subplot$BA_total1 + FIA_subplot$BA_total2)

# Convert Nan to NA
FIA_plot$BA_Nfixer_prop1[is.nan(FIA_plot$BA_Nfixer_prop1)] <- NA
FIA_plot$BA_Nfixer_prop2[is.nan(FIA_plot$BA_Nfixer_prop2)] <- NA

# Merge on state names
# fia_small <- fia[,c("STATECD", "PCN1", "STDAGE", "CARBON_LITTER", "CARBON_SOIL_ORG",
#                     "CARBON_UNDER_AG", "WATERCD1")] #add slope, aspect, elevm?
# fia_small <- fia_small[!duplicated(fia_small$PCN1),]
# FIA_plot <- merge(FIA_plot, fia_small, by.x="pcn", by.y="PCN1", all.x=T, all.y=F)

# Remove outliers from plot:
FIA_plot <- FIA_plot[FIA_plot$BA_change > quantile(FIA_plot$BA_change)[2] -
                       1.5*IQR(FIA_plot$BA_change) &
                       FIA_plot$BA_change < quantile(FIA_plot$BA_change)[4] + 
                       1.5*IQR(FIA_plot$BA_change), ]
FIA_plot <- FIA_plot[FIA_plot$BA_change_nonfixer > quantile(FIA_plot$BA_change_nonfixer, na.rm=T)[2] -
                       1.5*IQR(FIA_plot$BA_change_nonfixer, na.rm=T) & 
                       FIA_plot$BA_change_nonfixer < quantile(FIA_plot$BA_change_nonfixer, na.rm=T)[4] +
                       1.5*IQR(FIA_plot$BA_change_nonfixer, na.rm=T),]

# Recruitment and mortality rates:
#FIA_plot$mortyr <- FIA_plot$died / (FIA_plot$mt2 - FIA_plot$mt1)
#FIA_plot$recyr <- FIA_plot$rec / (FIA_plot$mt2 - FIA_plot$mt1)
FIA_plot$nt1 <- FIA_plot$n - FIA_plot$nrec
FIA_plot$nt2 <- FIA_plot$n - FIA_plot$ndied
#FIA_plot$mortyr <- (log(FIA_plot$nt) - log(FIA_plot$nt1)) / (FIA_plot$mt2 - FIA_plot$mt1)
FIA_plot$recyr <- FIA_plot$nrec / (FIA_plot$nt1 * (FIA_plot$mt2 - FIA_plot$mt1))
#FIA_plot$mortyr <- FIA_plot$died / (FIA_plot$nt1 * (FIA_plot$mt2 - FIA_plot$mt1))
FIA_plot$survyr <- (FIA_plot$nt1 - FIA_plot$ndied) / (FIA_plot$nt1 * (FIA_plot$mt2 - FIA_plot$mt1))

# Remove plots that started with 0 trees
#FIA_plot <- FIA_plot[FIA_plot$nt != 0, ]
FIA_plot <- FIA_plot[FIA_plot$nt1 != 0, ]

FIA_plot <- data.frame(FIA_plot)

# get number that survived
FIA_plot$nsurv <- FIA_plot$n - FIA_plot$ndied - FIA_plot$nrec

saveRDS(FIA_plot, "output/fia_plot_foranalysis.RDS")
