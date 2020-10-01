
############################
#
# This code is for calculating dbh growth rate of trees in FIA census.
# Used to calculate N-fixer crowding, graph, and evaluate significance
#
# adapted from Sian KG & Anika Petach
# 7/8/18
#
############################
# Note: TCN, PCNs messed up b/c expecting integer64 but can't fread an RDS

rm(list=ls())

library(dplyr)
library(data.table)
library(bit64)
library(svMisc)
library(profvis)

setwd("/Users/Anika/Documents/GradSchool/FIA_CompetitionModel/")


#----------------------------------------------------------------------------
# Load tree time series (made in surv_rec_lib.R)
FIA0 <- readRDS("data/derived/tree_ts.RDS")

FIA <- FIA0[!is.na(FIA0$DIST1) | !is.na(FIA0$DIST2),] #get rid of rows where all NA
rm(FIA0)
FIA$SPCD <- FIA$SPCD1
FIA$SPCD[!is.na(FIA$SPCD2)] = FIA$SPCD2[!is.na(FIA$SPCD2)]

# or load here
FIA <- readRDS("output/t_ts_new.RDS")

#merge fixr list 
# Load fixer list
fixer <- read.csv("data/raw/FIA_Fixer_spcd.csv", header = T) 
tree_spcd <- merge(FIA, fixer, 
                   by = "SPCD", 
                   all.x = TRUE, all.y = FALSE)
FIA <- tree_spcd
rm(tree_spcd)
FIA[is.na(FIA$FIX),"FIX"] <- 0

#New tree #:
FIA$ts3 <- seq(1, nrow(FIA), 1)

#Convert to data.table
FIA <- as.data.table(FIA)

#Constants
radius <- 7.3152 #(m)
a <- pi*radius^2

#Take only plots after 2000
#FIA <- FIA[FIA$MEASYEAR1.x > 2000, ]
FIA <- FIA[FIA$INVYR1 > 2000, ]

#Standard plot numbers
#FIA$PCN2 <- FIA$CN2.x
#FIA$PCN2[!is.na(FIA$CN2.y)] = FIA$CN2.y[!is.na(FIA$CN2.y)]
#FIA$PCN1 <- FIA$CN1.x
#FIA$PCN1[!is.na(FIA$CN1.y)] = FIA$CN1.y[!is.na(FIA$CN1.y)]

#Delete unnecessary data
# FIA <- FIA %>% 
#   dplyr::select(-MEASYEAR2.x, -WATERCD2.x, -CN0.x, -PLOT1.x, -MEASYEAR1.x,
#                              -tPLOT2, -CN1.x, -CN1.y, -CN2.x, -CN2.y, -CN0.y,
#                              -MEASYEAR2.y, -MEASYEAR1.y, -LAT1.y, -LON1.y, -ELEV1.y)

FIAa <- FIA
rm(FIA)

#subset just t2 (gets growth)
#FIA <- FIAa %>% 
#  dplyr::select(TCN2, AZIMUTH2, DIST2, DIA2, FIX, ts3, SUBP2, PCN2, Genus)

#select only trees where DIST2 and DIA2 are not na
#FIA <- FIA %>% 
#  dplyr::filter(!is.na(DIST2) & !is.na(DIA2))

#
saveRDS(FIA, "output/FIA_forNCIprop.RDS")

# NCI calc -----------------------------------------------------------------
#load data
FIA <- readRDS("output/FIA_forNCIprop.RDS")


#for PCN2, TREE2
output <- NULL
plots <- as.numeric(levels(factor(FIA$PCN2))) #PCN2 on grow, mort, surv
#run analysis on quarters of dataset
aa <- nrow(FIA)/4
#for(i in 1:aa){
#for(i in 1:length(plots)){

#profvis({
for(i in 1:10){
  #progress(i, length(plots))
  print(i/length(plots)*100)
  #plot_df <- dplyr::filter(FIA, PCN2==plots[i]) # FIA for plot
  plot_df <- FIA[FIA$PCN2 == plots[i], ]
  #plot_df <- dplyr::filter(plot_df, !is.na(DIST2))
  plot_df <- plot_df[!is.na(plot_df$DIST2), ]
  subplots <- as.numeric(levels(factor(plot_df$SUBP2))) #vector of subplots in plot
  
  if(nrow(plot_df) == 0){
    print("skip")
  }else{
    for(j in 1:length(subplots)){
      #subplot_df <- dplyr::filter(plot_df, SUBP2==subplots[j])   # FIA for subplot
      subplot_df <- plot_df[plot_df$SUBP2 == subplots[j], ]
      subplot_df$NCI_unweighted_Nfixer <- 0
      subplot_df$NCI_unweighted_nonfixer <- 0
      trees <- as.numeric(levels(factor(subplot_df$ts3)))  # vector of trees in subplot
      
      for(k in 1:length(trees)){
        #tree_wo_k <- subplot_df %>% filter(ts3 != trees[k])  # FIA for subplot excluding focaltree(k)
        tree_wo_k <- subplot_df[subplot_df$ts3 != trees[k], ]
        n <- which(subplot_df$ts3 == trees[k])
        tree_wo_k$dist_to_k <- sqrt((tree_wo_k$DIST2)^2
                                    + (subplot_df$DIST2[n])^2
                                    - 2*tree_wo_k$DIST2*
                                      subplot_df$DIST2[n]*
                                      cos(pi/180*(tree_wo_k$AZIMUTH2-
                                                    subplot_df$AZIMUTH2[n])))
        
        #tree_wo_k_neighbour <- tree_wo_k %>% filter(dist_to_k<radius & dist_to_k>0)
        tree_wo_k_neighbour <- 
          tree_wo_k[tree_wo_k$dist_to_k < radius & tree_wo_k$dist_to_k > 0, ]   # tree_wo_k within radius (m) of k
        tree_wo_k_neighbour_Nfixer <- 
          tree_wo_k_neighbour[tree_wo_k_neighbour$FIX == 1, ]
        tree_wo_k_neighbour_nonfixer <- 
          tree_wo_k_neighbour[tree_wo_k_neighbour$FIX == 0, ]
      
        subplot_df$NCI_unweighted_Nfixer[n]<-
          sum((tree_wo_k_neighbour_Nfixer$DIA2/100)^2/tree_wo_k_neighbour_Nfixer$dist_to_k^2,na.rm=T)
        subplot_df$NCI_unweighted_nonfixer[n]<-
          sum((tree_wo_k_neighbour_nonfixer$DIA2/100)^2/tree_wo_k_neighbour_nonfixer$dist_to_k^2,na.rm=T)
      } #end k loop through trees
      
      subplot_df$area <- 
        2*radius^2*acos(subplot_df$DIST2^2/(2*subplot_df$DIST2*radius))-
        0.5*sqrt((-subplot_df$DIST2+2*radius)
                 *(subplot_df$DIST2)
                 *(subplot_df$DIST2)
                 *(subplot_df$DIST2+2*radius))
      
      subplot_df$NCI_weighted_Nfixer <- 
        subplot_df$NCI_unweighted_Nfixer*subplot_df$area/(a)
      subplot_df$NCI_weighted_nonfixer <- 
        subplot_df$NCI_unweighted_nonfixer*subplot_df$area/(a)
      
     # output1 <- subplot_df[,c("PCN2", "TCN2", "SUBP2", "NCI_weighted_Nfixer",
      #                         "NCI_weighted_nonfixer")]
      ##output <- rbind(output,output1)
    #  output <- rbindlist(list(output, output1))
      
      for(k in 1:length(trees)){
        l <- which(FIA$ts3==trees[k])
        n <- which(subplot_df$ts3==trees[k])
        FIA[l,"NCI_Nfixer2"] <- subplot_df$NCI_weighted_Nfixer[n]
        FIA[l,"NCI_nonfixer2"] <- subplot_df$NCI_weighted_nonfixer[n]
      }
    } #end j loop through subplots
  } #end else
  
} #end plot loop
#}) #end profvis


#subset just t1
#select only trees where DIST1 and DIA1 are not na
#run analysis on quarters of dataset
#for PCN1, TREE1

plots <- as.numeric(levels(factor(FIA$PCN1))) #PCN2 on grow, mort, surv

output0 <- NULL
for(i in 1:length(plots)){
  #for(i in 1:5){
  print(i/length(plots)*100)
  #progress(i, length(plots))
  plot_df <- dplyr::filter(FIA, PCN1 == plots[i])  # FIA for plot
  plot_df <- dplyr::filter(plot_df, !is.na(DIST1))
  subplots <- as.numeric(levels(factor(plot_df$SUBP1))) #vector of subplots in plot
  
  if(nrow(plot_df) == 0){
    print("skip")
  }else{
    for(j in 1:length(subplots)){
      subplot_df <- dplyr::filter(plot_df, SUBP1 == subplots[j])   # FIA for subplot
      subplot_df$NCI_unweighted_Nfixer <- 0
      subplot_df$NCI_unweighted_nonfixer <- 0
      trees <- as.numeric(levels(factor(subplot_df$ts3))) # vector of trees in subplot
      
      for(k in 1:length(trees)){
        tree_wo_k <- subplot_df %>% filter(ts3 != trees[k])  # FIA for subplot excluding focaltree(k)
        n <- which(subplot_df$ts3 == trees[k])
        tree_wo_k$dist_to_k <- 
          sqrt((tree_wo_k$DIST1)^2 + 
                 (subplot_df$DIST1[n])^2 - 
                 2*tree_wo_k$DIST1*subplot_df$DIST1[n]*
                                      cos(pi/180*(tree_wo_k$AZIMUTH1-
                                                    subplot_df$AZIMUTH1[n])))
        
        tree_wo_k_neighbour <- tree_wo_k %>% 
          filter(dist_to_k < radius & dist_to_k > 0)  # tree_wo_k within radius (m) of k
        tree_wo_k_neighbour_Nfixer <- 
          tree_wo_k_neighbour[tree_wo_k_neighbour$FIX == 1, ]
        tree_wo_k_neighbour_nonfixer <- 
          tree_wo_k_neighbour[tree_wo_k_neighbour$FIX == 0, ]
        
        subplot_df$NCI_unweighted_Nfixer[n] <-
          sum((tree_wo_k_neighbour_Nfixer$DIA1/100)^2/tree_wo_k_neighbour_Nfixer$dist_to_k^2, na.rm = T)
        subplot_df$NCI_unweighted_nonfixer[n] <-
          sum((tree_wo_k_neighbour_nonfixer$DIA1/100)^2/tree_wo_k_neighbour_nonfixer$dist_to_k^2, na.rm = T)
      } #end k loop through trees
      
      subplot_df$area <- 
        2*radius^2*acos(subplot_df$DIST1^2/(2*subplot_df$DIST1*radius)) -
        0.5*sqrt((-subplot_df$DIST1+2*radius)
                 *(subplot_df$DIST1)
                 *(subplot_df$DIST1)
                 *(subplot_df$DIST1+2*radius))
      subplot_df$NCI_weighted_Nfixer <- 
        subplot_df$NCI_unweighted_Nfixer*subplot_df$area/(a)
      subplot_df$NCI_weighted_nonfixer <- 
        subplot_df$NCI_unweighted_nonfixer*subplot_df$area/(a)
      
      output2 <- subplot_df[,c("PCN1", "TREE1", "PREV_TCN2", "SUBP1", 
                               "NCI_weighted_Nfixer", "NCI_weighted_nonfixer")]
      output0 <- rbind(output0, output2)
      
      
      # for(k in 1:length(trees)){
      #   l <- which(FIA$ts3==trees[k])
      #   n <- which(subplot_df$ts3==trees[k])
      #   FIA[l,"NCI_Nfixer1"] <- subplot_df$NCI_weighted_Nfixer[n]
      #   FIA[l,"NCI_nonfixer1"] <- subplot_df$NCI_weighted_nonfixer[n]
      # }
    } #end j loop through subplots
  } #end else
  
} #end i plots loop

saveRDS(output0,"output/output0.RDS")

colnames(output0) <- c("PCN1","TREE1","PREV_TCN2","SUBP1","NCI_Nfixer1","NCI_nonfixer1")
test <- merge(FIA, output0, 
              by = "PREV_TCN2", 
              all = TRUE)

colnames(output) <- c("PCN2", "TREE2", "TCN2", "SUBP2", "NCI_Nfixer2", "NCI_nonfixer2")
test1 <- merge(test, output, 
               by = "TCN2",
               all = T)

saveRDS(test1,"output/FIA_withNCI.RDS")


# NCI by genus calc ------------------------------------------------------
genuslist <- c("Acacia", "Albizia", "Alnus", 
               "Cercocarpus", "Elaeagnus", "Olneya", 
               "Prosopis", "Robinia")

FIA <- readRDS("output/FIA_forNCIprop.RDS")

genusi <- genuslist[1]

#Constants
radius <- 7.3152 #(m)
a <- pi*radius^2

#for PCN2, TREE2
output <- NULL
plots <- as.numeric(levels(factor(FIA$PCN2))) #PCN2 on grow, mort, surv
#run analysis on quarters of dataset
aa <- nrow(FIA)/4
#for(i in 1:aa){
#for(i in 1:length(plots)){

#profvis({
for(i in 1:length(plots)){
  #progress(i, length(plots))
  print(i/length(plots)*100)
  plot_df <- FIA[FIA$PCN2 == plots[i], ]
  plot_df <- plot_df[!is.na(plot_df$DIST2), ]
  subplots <- as.numeric(levels(factor(plot_df$SUBP2))) #vector of subplots in plot
  
  if(nrow(plot_df) == 0){
    print("skip")
  }else{
    for(j in 1:length(subplots)){
      subplot_df <- plot_df[plot_df$SUBP2 == subplots[j], ]
      subplot_df$NCI_unweighted_Nfixer <- 0
      subplot_df$NCI_unweighted_nonfixer <- 0
      trees <- as.numeric(levels(factor(subplot_df$ts3)))  # vector of trees in subplot
      
      for(k in 1:length(trees)){
        tree_wo_k <- subplot_df[subplot_df$ts3 != trees[k], ]
        n <- which(subplot_df$ts3 == trees[k])
        tree_wo_k$dist_to_k <- sqrt((tree_wo_k$DIST2)^2
                                    + (subplot_df$DIST2[n])^2
                                    - 2*tree_wo_k$DIST2*
                                      subplot_df$DIST2[n]*
                                      cos(pi/180*(tree_wo_k$AZIMUTH2-
                                                    subplot_df$AZIMUTH2[n])))
        
        tree_wo_k_neighbour <- 
          tree_wo_k[tree_wo_k$dist_to_k < radius & tree_wo_k$dist_to_k > 0, ]   # tree_wo_k within radius (m) of k
        tree_wo_k_neighbour_Nfixer <- 
          tree_wo_k_neighbour[tree_wo_k_neighbour$Genus %in% c(genusi), ]
        tree_wo_k_neighbour_nonfixer <- 
          tree_wo_k_neighbour[!(tree_wo_k_neighbour$Genus %in% c(genusi)), ]
        
        subplot_df$NCI_unweighted_Nfixer[n]<-
          sum((tree_wo_k_neighbour_Nfixer$DIA2/100)^2/tree_wo_k_neighbour_Nfixer$dist_to_k^2, na.rm=T)
        subplot_df$NCI_unweighted_nonfixer[n]<-
          sum((tree_wo_k_neighbour_nonfixer$DIA2/100)^2/tree_wo_k_neighbour_nonfixer$dist_to_k^2, na.rm=T)
      } #end k loop through trees
      
      subplot_df$area <- 
        2*radius^2*acos(subplot_df$DIST2^2/(2*subplot_df$DIST2*radius))-
        0.5*sqrt((-subplot_df$DIST2+2*radius)
                 *(subplot_df$DIST2)
                 *(subplot_df$DIST2)
                 *(subplot_df$DIST2+2*radius))
      
      subplot_df$NCI_weighted_Nfixer <- 
        subplot_df$NCI_unweighted_Nfixer*subplot_df$area/(a)
      subplot_df$NCI_weighted_nonfixer <- 
        subplot_df$NCI_unweighted_nonfixer*subplot_df$area/(a)
      
      # output1 <- subplot_df[,c("PCN2", "TCN2", "SUBP2", "NCI_weighted_Nfixer",
      #                         "NCI_weighted_nonfixer")]
      ##output <- rbind(output,output1)
      #  output <- rbindlist(list(output, output1))
      
      for(k in 1:length(trees)){
        l <- which(FIA$ts3 == trees[k])
        n <- which(subplot_df$ts3 == trees[k])
        FIA[l,"NCI_Nfixer2"] <- subplot_df$NCI_weighted_Nfixer[n]
        FIA[l,"NCI_nonfixer2"] <- subplot_df$NCI_weighted_nonfixer[n]
      }
    } #end j loop through subplots
  } #end else
  
} #end plot loop
#}) #end profvis


#subset just t1
#select only trees where DIST1 and DIA1 are not na
#run analysis on quarters of dataset
#for PCN1, TREE1

plots <- as.numeric(levels(factor(FIA$PCN1))) #PCN2 on grow, mort, surv

output0 <- NULL
for(i in 1:length(plots)){
  print(i/length(plots)*100)
  #progress(i, length(plots))
  plot_df <- dplyr::filter(FIA, PCN1 == plots[i])  # FIA for plot
  plot_df <- dplyr::filter(plot_df, !is.na(DIST1))
  subplots <- as.numeric(levels(factor(plot_df$SUBP1))) #vector of subplots in plot
  
  if(nrow(plot_df) == 0){
    print("skip")
  }else{
    for(j in 1:length(subplots)){
      subplot_df <- dplyr::filter(plot_df, SUBP1 == subplots[j])   # FIA for subplot
      subplot_df$NCI_unweighted_Nfixer <- 0
      subplot_df$NCI_unweighted_nonfixer <- 0
      trees <- as.numeric(levels(factor(subplot_df$ts3))) # vector of trees in subplot
      
      for(k in 1:length(trees)){
        tree_wo_k <- subplot_df %>% filter(ts3 != trees[k])  # FIA for subplot excluding focaltree(k)
        n <- which(subplot_df$ts3 == trees[k])
        tree_wo_k$dist_to_k <- 
          sqrt((tree_wo_k$DIST1)^2 + 
                 (subplot_df$DIST1[n])^2 - 
                 2*tree_wo_k$DIST1*subplot_df$DIST1[n]*
                 cos(pi/180*(tree_wo_k$AZIMUTH1-
                               subplot_df$AZIMUTH1[n])))
        
        tree_wo_k_neighbour <- tree_wo_k %>% 
          filter(dist_to_k < radius & dist_to_k > 0)  # tree_wo_k within radius (m) of k
        tree_wo_k_neighbour_Nfixer <- 
          tree_wo_k_neighbour[tree_wo_k_neighbour$Genus == genusi, ]
        tree_wo_k_neighbour_nonfixer <- 
          tree_wo_k_neighbour[tree_wo_k_neighbour$Genus != genusi, ]
        
        subplot_df$NCI_unweighted_Nfixer[n] <-
          sum((tree_wo_k_neighbour_Nfixer$DIA1/100)^2/tree_wo_k_neighbour_Nfixer$dist_to_k^2, na.rm = T)
        subplot_df$NCI_unweighted_nonfixer[n] <-
          sum((tree_wo_k_neighbour_nonfixer$DIA1/100)^2/tree_wo_k_neighbour_nonfixer$dist_to_k^2, na.rm = T)
      } #end k loop through trees
      
      subplot_df$area <- 
        2*radius^2*acos(subplot_df$DIST1^2/(2*subplot_df$DIST1*radius)) -
        0.5*sqrt((-subplot_df$DIST1+2*radius)
                 *(subplot_df$DIST1)
                 *(subplot_df$DIST1)
                 *(subplot_df$DIST1+2*radius))
      subplot_df$NCI_weighted_Nfixer <- 
        subplot_df$NCI_unweighted_Nfixer*subplot_df$area/(a)
      subplot_df$NCI_weighted_nonfixer <- 
        subplot_df$NCI_unweighted_nonfixer*subplot_df$area/(a)
      
      output2 <- subplot_df[,c("PCN1", "TREE1", "PREV_TCN2", "SUBP1", 
                               "NCI_weighted_Nfixer", "NCI_weighted_nonfixer")]
      output0 <- rbind(output0, output2)
      
      
      # for(k in 1:length(trees)){
      #   l <- which(FIA$ts3==trees[k])
      #   n <- which(subplot_df$ts3==trees[k])
      #   FIA[l,"NCI_Nfixer1"] <- subplot_df$NCI_weighted_Nfixer[n]
      #   FIA[l,"NCI_nonfixer1"] <- subplot_df$NCI_weighted_nonfixer[n]
      # }
    } #end j loop through subplots
  } #end else
  
} #end i plots loop

saveRDS(output0,"output/output0_withgenus.RDS")

colnames(output0) <- c("PCN1","TREE1","PREV_TCN2","SUBP1","NCI_Nfixer1","NCI_nonfixer1")
test <- merge(FIA, output0, 
              by = "PREV_TCN2", 
              all = TRUE)

colnames(output) <- c("PCN2", "TREE2", "TCN2", "SUBP2", "NCI_Nfixer2", "NCI_nonfixer2")
test1 <- merge(test, output, 
               by = "TCN2",
               all = T)

saveRDS(test1,"output/FIA_withNCI_withgenus.RDS")