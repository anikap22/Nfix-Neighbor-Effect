

library(dplyr, lib.loc = "/rigel/home/arp2195/rpackages")
library(data.table, lib.loc = "/rigel/home/arp2195/rpackages")

#setwd("/Users/Anika/Documents/GradSchool/FIA_CompetitionModel/")

#don't use sci notation
options(scipen=6)

# extract state from arguments
args = commandArgs(trailingOnly=TRUE)
state <- args[1]

#Constants
radius <- 7.3152 #(m)
a <- pi*radius^2

# Load time series
#FIA <- readRDS("output/FIA_forNCI_05_19.RDS")
filename <- paste("fia_hab_", state, ".RDS", sep="")
FIA <- readRDS(filename)

# Filter state
#FIA <- dplyr::filter(FIA, state == statelist[1])

#for PCN2, TREE2
plots <- unique(FIA$pcn2) #69464 plots

#profvis({
for(i in 1:length(plots)){
  #progress(i, length(plots))
  print(i/length(plots)*100)
  plot_df <- FIA[FIA$pcn2 == plots[i], ]
  plot_df <- plot_df[!is.na(plot_df$tcn2), ]
  subplots <- as.numeric(levels(factor(plot_df$subp2))) #vector of subplots in plot
  
  if(nrow(plot_df) == 0){
    print("skip")
  }else{
    for(j in 1:length(subplots)){
      subplot_df <- plot_df[plot_df$subp2 == subplots[j], ]
      subplot_df$NCI_unweighted_Nfixer <- 0
      subplot_df$NCI_unweighted_nonfixer <- 0
      trees <- unique(subplot_df$ts3) #vector of trees in subplot
      
      for(k in 1:length(trees)){
        tree_wo_k <- subplot_df[subplot_df$ts3 != trees[k], ]
        n <- which(subplot_df$ts3 == trees[k])
        tree_wo_k$dist_to_k <- sqrt((tree_wo_k$DISTm)^2
                                    + (subplot_df$DISTm[n])^2
                                    - 2*tree_wo_k$DISTm*
                                      subplot_df$DISTm[n]*
                                      cos(tree_wo_k$AZIMUTHr-
                                                    subplot_df$AZIMUTHr[n]))
        
        tree_wo_k_neighbour <- 
          tree_wo_k[tree_wo_k$dist_to_k < radius & tree_wo_k$dist_to_k > 0, ]   # tree_wo_k within radius (m) of k
        tree_wo_k_neighbour_Nfixer <- 
          tree_wo_k_neighbour[tree_wo_k_neighbour$FIX == 1, ]
        tree_wo_k_neighbour_nonfixer <- 
          tree_wo_k_neighbour[tree_wo_k_neighbour$FIX == 0, ]
        
        subplot_df$NCI_unweighted_Nfixer[n]<-
          sum((tree_wo_k_neighbour_Nfixer$dbh2 / 100)^2 / 
                tree_wo_k_neighbour_Nfixer$dist_to_k^2, na.rm=T)
        subplot_df$NCI_unweighted_nonfixer[n]<-
          sum((tree_wo_k_neighbour_nonfixer$dbh2 / 100)^2 / 
                tree_wo_k_neighbour_nonfixer$dist_to_k^2, na.rm=T)
      } #end k loop through trees
      
      # calc fraction of sampling radius that was sampled
      # subplot_df$area <- 
      #   2 * radius^2 * acos(subplot_df$DISTm^2 / (2 * subplot_df$DISTm * radius))-
      #   0.5 * sqrt((-subplot_df$DISTm + 2 * radius)
      #              *(subplot_df$DISTm)
      #              *(subplot_df$DISTm)
      #              *(subplot_df$DISTm + 2 * radius))
      subplot_df$area <- 
        (2 * radius^2 * acos(subplot_df$DISTm / (2 * radius))-
           0.5 * subplot_df$DISTm * sqrt(-subplot_df$DISTm^2 + 4 * radius^2))
      
      # weight NCI values by sampling area
      # replace first / with * to do original weighting ( / mirrors forest out)
      subplot_df$NCI_weighted_Nfixer <- 
        subplot_df$NCI_unweighted_Nfixer / (subplot_df$area / (a))
      
      subplot_df$NCI_weighted_nonfixer <- 
        subplot_df$NCI_unweighted_nonfixer / (subplot_df$area / (a))
      
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

plots <- unique(FIA$pcn1)

for(i in 1:length(plots)){
  print(i/length(plots)*100)
  plot_df <- FIA[FIA$pcn1 == plots[i], ]
  plot_df <- plot_df[!is.na(plot_df$tcn1), ]
  subplots <- as.numeric(levels(factor(plot_df$subp1))) #vector of subplots in plot
  
  
  if(nrow(plot_df) == 0){
    print("skip")
  }else{
    for(j in 1:length(subplots)){
      subplot_df <- plot_df[plot_df$subp1 == subplots[j], ]
      subplot_df$NCI_unweighted_Nfixer <- 0
      subplot_df$NCI_unweighted_nonfixer <- 0
      trees <- as.numeric(levels(factor(subplot_df$ts3))) # vector of trees in subplot
      
      for(k in 1:length(trees)){
        tree_wo_k <- subplot_df[subplot_df$ts3 != trees[k], ]
        n <- which(subplot_df$ts3 == trees[k])
        
        tree_wo_k$dist_to_k <- sqrt((tree_wo_k$DISTm)^2
                                    + (subplot_df$DISTm[n])^2
                                    - 2*tree_wo_k$DISTm*
                                      subplot_df$DISTm[n]*
                                      cos(tree_wo_k$AZIMUTHr-
                                                    subplot_df$AZIMUTHr[n]))
        
        tree_wo_k_neighbour <- 
          tree_wo_k[tree_wo_k$dist_to_k < radius & tree_wo_k$dist_to_k > 0, ]   # tree_wo_k within radius (m) of k
        
        tree_wo_k_neighbour_Nfixer <- 
          tree_wo_k_neighbour[tree_wo_k_neighbour$FIX == 1, ]
        
        tree_wo_k_neighbour_nonfixer <- 
          tree_wo_k_neighbour[tree_wo_k_neighbour$FIX == 0, ]
        
        subplot_df$NCI_unweighted_Nfixer[n]<-
          sum((tree_wo_k_neighbour_Nfixer$dbh1 / 100)^2 / 
                tree_wo_k_neighbour_Nfixer$dist_to_k^2, na.rm=T)
        
        subplot_df$NCI_unweighted_nonfixer[n]<-
          sum((tree_wo_k_neighbour_nonfixer$dbh1 / 100)^2 / 
                tree_wo_k_neighbour_nonfixer$dist_to_k^2, na.rm=T)
      } #end k loop through trees
      
      # calc fraction of sampling radius that was sampled
      # subplot_df$area <- 
      #   2 * radius^2 * acos(subplot_df$DISTm^2 / (2 * subplot_df$DISTm * radius))-
      #   0.5 * sqrt((-subplot_df$DISTm + 2 * radius)
      #              *(subplot_df$DISTm)
      #              *(subplot_df$DISTm)
      #              *(subplot_df$DISTm + 2 * radius))
      
      subplot_df$area <- 
        (2 * radius^2 * acos(subplot_df$DISTm / (2 * radius))-
        0.5 * subplot_df$DISTm * sqrt(-subplot_df$DISTm^2 + 4 * radius^2))
                  
      
      # weight NCI values by sampling area
      # replace first / with * to do original weighting ( / mirrors forest out)
      subplot_df$NCI_weighted_Nfixer <- 
        subplot_df$NCI_unweighted_Nfixer / (subplot_df$area / (a))
      
      subplot_df$NCI_weighted_nonfixer <- 
        subplot_df$NCI_unweighted_nonfixer / (subplot_df$area / (a))
      
      for(k in 1:length(trees)){
        l <- which(FIA$ts3 == trees[k])
        n <- which(subplot_df$ts3 == trees[k])
        FIA[l,"NCI_Nfixer1"] <- subplot_df$NCI_weighted_Nfixer[n]
        FIA[l,"NCI_nonfixer1"] <- subplot_df$NCI_weighted_nonfixer[n]
      }
    } #end j loop through subplots
  } #end else
  
} #end plot loop
#}) #end profvis

# checking why we have NAs in NCI
#temp1 <- temp %>% group_by(pcn1) %>% filter(n()>1)
#temp1 <- temp[,if (.N>1) .SD, by=.(pcn1)]


filename <- paste("nci_", state, ".RDS", sep="")
saveRDS(FIA, filename)
