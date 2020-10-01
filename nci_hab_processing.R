

############################
# Habanero prep ############
############################

setwd("/Users/Anika/Documents/GradSchool/FIA_CompetitionModel/output/")

#don't use sci notation
options(scipen=6)

# Load time series
FIA <- readRDS("FIA_forNCI_05_19.RDS")

statelist <- unique(FIA$state)

# Filter by state and save
for(i in 1:length(statelist)){
  data <- dplyr::filter(FIA, state == statelist[i])
  filename <- paste("fia_hab_", statelist[i], ".RDS", sep="")
  saveRDS(data, filename)
}



temp <- readRDS("hab/nci_MI.RDS")



#########################
# post Habanero #########
#########################

setwd("/Users/Anika/Documents/GradSchool/FIA_CompetitionModel/")

fia <- NULL
n <- length(statelist)-2
for(i in 1:n){
  filename <- paste("output/hab/nci_", statelist[i], ".RDS", sep="")
  file <- readRDS(filename)
  fia <- rbind.data.frame(fia, file)
}

saveRDS(fia, "output/nci_fromHab.RDS")
