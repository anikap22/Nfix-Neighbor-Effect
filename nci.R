
############################
#
# Used to calculate N-fixer crowding, graph, and evaluate significance
#
# adapted from Sian KG & Anika Petach
# 9/26/18
#
############################

rm(list=ls())

library(dplyr)
library(data.table)
library(bit64)
library(svMisc)
library(profvis)

setwd("/Users/Anika/Documents/GradSchool/FIA_CompetitionModel/")

#----------------------------------------------------------------------------
# Load tree time series (made in join.R)
#FIA <- readRDS("data/tree_ts.RDS")
#FIA <- readRDS("output/FIA_joined_Sept18.RDS")
FIA <- readRDS("output/ts.RDS")

#FIA <- FIA[!is.na(FIA$DISTm1) | !is.na(FIA$DISTm2), ] #get rid of rows where all NA
FIA <- FIA[!is.na(FIA$DISTm),]

# Get general SPCD col
#FIA$SPCD <- FIA$SPCD1
#FIA[!is.na(FIA$SPCD2), ]$SPCD = FIA[!is.na(FIA$SPCD2), ]$SPCD2

# Merge fixer list 
fixer <- read.csv("data/raw/FIA_Fixer_spcd.csv", header = T) 
FIA <- merge(FIA, fixer, 
             by = "SPCD",
             all.x = TRUE, all.y = FALSE)
rm(fixer)

FIA[is.na(FIA$FIX),]$FIX <- 0

# New tree #:
FIA$ts3 <- seq(1, nrow(FIA), 1)

# Convert to data.table
#FIA <- as.data.table(FIA)

FIA <- tbl_dt(FIA) #lets you use data.table in dplyr funcs

# Constants
radius <- 7.3152 #(m)
a <- pi*radius^2


# remove ---------------
#subset just t2 (gets growth)
#FIA <- FIAa %>% 
#  dplyr::select(TCN2, AZIMUTH2, DIST2, DIA2, FIX, ts3, SUBP2, PCN2)

#select only trees where DIST2 and DIA2 are not na
#FIA <- FIA %>% 
#  dplyr::filter(!is.na(DIST2) & !is.na(DIA2))


# end remove -----------------

FIAcopy <- FIA
FIA <- FIAcopy


# NCI 1 ---------------------------------------------
output <- NULL
#plots <- as.numeric(levels(factor(FIA$PCN1))) #PCN2 on grow, mort, surv
plots <- as.integer64(levels(factor(FIA$PCN1))) #PCN2 on grow, mort, surv
#FIA <- FIA[FATE==2 | FATE==3,] #at time 1 only want trees that grew or died

#profvis({
  for(i in 1:length(plots)){ #1:length(plots)
    print(i/length(plots)*100)
    plot_df <- FIA[(FATE==2 | FATE==3) & (FIA$PCN1 == plots[i]), ] #at time 1 only want trees that grew or died
    #plot_df <- plot_df[!is.na(plot_df$DIST2), ]
    subplots <- as.numeric(levels(factor(plot_df$SUBP1))) #vector of subplots in plot
    
    if(nrow(plot_df) == 0){
      print("skip")
    }else{
      for(j in 1:length(subplots)){
        #subplot_df <- dplyr::filter(plot_df, SUBP2==subplots[j])   # FIA for subplot
        subplot_df <- plot_df[SUBP1 == subplots[j], ]
        #subplot_df <- subplot_df[!is.na(subplot_df$DISTm1), ] # Remove recruited trees
        subplot_df$NCI_unweighted_Nfixer <- 0
        subplot_df$NCI_unweighted_nonfixer <- 0
        trees <- as.numeric(levels(factor(subplot_df$ts3)))  # vector of trees in subplot
        
        for(k in 1:length(trees)){
          #tree_wo_k <- subplot_df %>% filter(ts3 != trees[k])  # FIA for subplot excluding focaltree(k)
          tree_wo_k <- subplot_df[subplot_df$ts3 != trees[k], ]
          n <- which(subplot_df$ts3 == trees[k])
          tree_wo_k$dist_to_k <- sqrt((tree_wo_k$DISTm)^2
                                      + (subplot_df$DISTm[n])^2
                                      - 2*tree_wo_k$DISTm*
                                        subplot_df$DISTm[n]*
                                        cos(pi/180*(tree_wo_k$AZIMUTH-
                                                      subplot_df$AZIMUTH[n])))
          
         
          #tree_wo_k_neighbour <- tree_wo_k %>% filter(dist_to_k<radius & dist_to_k>0)
          tree_wo_k_neighbour <- 
            tree_wo_k[(tree_wo_k$dist_to_k < radius) & (tree_wo_k$dist_to_k > 0), ]   # tree_wo_k within radius (m) of k
          tree_wo_k_neighbour_Nfixer <- 
            tree_wo_k_neighbour[tree_wo_k_neighbour$FIX == 1, ]
          tree_wo_k_neighbour_nonfixer <- 
            tree_wo_k_neighbour[tree_wo_k_neighbour$FIX == 0, ]
          
          subplot_df$NCI_unweighted_Nfixer[n]<-
            sum((tree_wo_k_neighbour_Nfixer$DIAcm1/100)^2/tree_wo_k_neighbour_Nfixer$dist_to_k^2, 
                na.rm = T)
          subplot_df$NCI_unweighted_nonfixer[n]<-
            sum((tree_wo_k_neighbour_nonfixer$DIAcm1/100)^2/tree_wo_k_neighbour_nonfixer$dist_to_k^2,
                na.rm = T)
        } #end k loop through trees
        
        subplot_df$area <- 
          2*radius^2*acos(subplot_df$DISTm^2/(2*subplot_df$DISTm*radius))-
          0.5*sqrt((-subplot_df$DISTm+2*radius)
                   *(subplot_df$DISTm)
                   *(subplot_df$DISTm)
                   *(subplot_df$DISTm+2*radius))
        
        subplot_df$NCI_weighted_Nfixer <- 
          subplot_df$NCI_unweighted_Nfixer*subplot_df$area/(a)
        subplot_df$NCI_weighted_nonfixer <- 
          subplot_df$NCI_unweighted_nonfixer*subplot_df$area/(a)
        
        #output1 <- subplot_df[,c("PCN2", "TCN2", "SUBP2", "NCI_weighted_Nfixer",
        #                         "NCI_weighted_nonfixer")]
        #output <- rbind(output,output1)
        #output <- rbindlist(list(output, output1))
       ## output <- rbindlist(list(output, subplot_df))
        
        for(k in 1:length(trees)){
        #l <- which(FIA$ts3==trees[k])
          n <- which(subplot_df$ts3 == trees[k])
          FIA[ts3==trees[k], NCI_weighted_Nfixer := subplot_df$NCI_weighted_Nfixer[n]]
          FIA[ts3==trees[k], NCI_weighted_nonfixer := subplot_df$NCI_weighted_nonfixer[n]]
          FIA[ts3==trees[k], area := subplot_df$area[n]]
        }
      } #end j loop through subplots
    } #end else
    
  } #end plot loop
#}) #end profvis

# NCI 2 ---------------------------------------------
output <- NULL
plots <- as.integer64(levels(factor(FIA$PCN2)))

#profvis({
for(i in 1:length(plots)){ #1:length(plots)
  print(i/length(plots)*100)
  plot_df <- FIA[(FATE==1 | FATE==2) & (FIA$PCN2 == plots[i]), ] #at time 1 only want trees that grew or died
  subplots <- as.numeric(levels(factor(plot_df$SUBP2))) #vector of subplots in plot
  
  if(nrow(plot_df) == 0){
    print("skip")
  }else{
    for(j in 1:length(subplots)){
      subplot_df <- plot_df[SUBP2 == subplots[j], ]
      subplot_df$NCI_unweighted_Nfixer <- 0
      subplot_df$NCI_unweighted_nonfixer <- 0
      trees <- as.numeric(levels(factor(subplot_df$ts3)))  # vector of trees in subplot
      
      for(k in 1:length(trees)){
        tree_wo_k <- subplot_df[subplot_df$ts3 != trees[k], ]
        n <- which(subplot_df$ts3 == trees[k])
        tree_wo_k$dist_to_k <- sqrt((tree_wo_k$DISTm)^2
                                    + (subplot_df$DISTm[n])^2
                                    - 2*tree_wo_k$DISTm*
                                      subplot_df$DISTm[n]*
                                      cos(pi/180*(tree_wo_k$AZIMUTH-
                                                    subplot_df$AZIMUTH[n])))
        
        
        tree_wo_k_neighbour <- 
          tree_wo_k[(tree_wo_k$dist_to_k < radius) & (tree_wo_k$dist_to_k > 0), ]   # tree_wo_k within radius (m) of k
        tree_wo_k_neighbour_Nfixer <- 
          tree_wo_k_neighbour[tree_wo_k_neighbour$FIX == 1, ]
        tree_wo_k_neighbour_nonfixer <- 
          tree_wo_k_neighbour[tree_wo_k_neighbour$FIX == 0, ]
        
        subplot_df$NCI_unweighted_Nfixer[n]<-
          sum((tree_wo_k_neighbour_Nfixer$DIAcm2/100)^2/tree_wo_k_neighbour_Nfixer$dist_to_k^2, 
              na.rm = T)
        subplot_df$NCI_unweighted_nonfixer[n]<-
          sum((tree_wo_k_neighbour_nonfixer$DIAcm2/100)^2/tree_wo_k_neighbour_nonfixer$dist_to_k^2,
              na.rm = T)
      } #end k loop through trees
      
      subplot_df$area <- 
        2*radius^2*acos(subplot_df$DISTm^2/(2*subplot_df$DISTm*radius))-
        0.5*sqrt((-subplot_df$DISTm+2*radius)
                 *(subplot_df$DISTm)
                 *(subplot_df$DISTm)
                 *(subplot_df$DISTm+2*radius))
      
      subplot_df$NCI_weighted_Nfixer <- 
        subplot_df$NCI_unweighted_Nfixer*subplot_df$area/(a)
      subplot_df$NCI_weighted_nonfixer <- 
        subplot_df$NCI_unweighted_nonfixer*subplot_df$area/(a)
      
      #output1 <- subplot_df[,c("PCN2", "TCN2", "SUBP2", "NCI_weighted_Nfixer",
      #                         "NCI_weighted_nonfixer")]
      #output <- rbind(output,output1)
      #output <- rbindlist(list(output, output1))
      ## output <- rbindlist(list(output, subplot_df))
      
      for(k in 1:length(trees)){
        n <- which(subplot_df$ts3 == trees[k])
        FIA[ts3==trees[k], NCI_weighted_Nfixer2 := subplot_df$NCI_weighted_Nfixer[n]]
        FIA[ts3==trees[k], NCI_weighted_nonfixer2 := subplot_df$NCI_weighted_nonfixer[n]]
        FIA[ts3==trees[k], area2 := subplot_df$area[n]]
      }
    } #end j loop through subplots
  } #end else
  
} #end plot loop
#}) #end profvis

saveRDS(FIA, "FIA_nci_2.RDS")

