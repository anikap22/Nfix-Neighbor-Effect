

####################
#
#
# models_ind_rec.R
#
# 1/7/19
# Anika Petach
# 
#
#####################

require(bit64)
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

setwd("/Users/Anika/Documents/GradSchool/FIA_CompetitionModel/")

############################# RECRUITMENT ########################
# Load data ------------------------------------------------

rs <- readRDS("output/nci_formodels_recr.RDS") #made in growth_prep.R

spdata <- read.csv("/Users/Anika/Documents/GradSchool/FIA_EstimateProj/ACGCA_FIAv51_from_JWL/REF_SPECIES_jwl.csv")
spdata <- spdata %>% dplyr::select(SPCD, COMMON_NAME, GENUS, SPECIES, JENKINS_SPGRPCD) #get biomass parameters if needed

rs <- merge(rs, spdata, by.x="SPCD", by.y="SPCD", all.x=T, all.y=F)

source("scripts/functions.R")


# Add whether recruited or not -------------------------------------------
rs$reci <- 0
rs[rs$FATE==1,]$reci <- 1
rs$reci <- rs$reci/rs$t #probability of recruitment per year

# Add decid vs evergreen ---------------------------------------------------
classten <- read.csv("data/raw/Jenkins10_decid0ever1.csv")
classten <- classten[, c("SPCD", "DEC0EVER1")]

rs <- merge(rs, classten, by.x="SPCD", by.y="SPCD", all.x=T)
rs[rs$JENKINS_SPGRPCD <= 5, "DEC0EVER1"] <- 1
rs[rs$JENKINS_SPGRPCD >= 6 & rs$JENKINS_SPGRPCD < 10, "DEC0EVER1"] <- 0


# Subset data --------------------------------------------------------

young <- rs %>% dplyr::filter(STDAGE < 60)
old <- rs %>% dplyr::filter(STDAGE >= 60)

f <- rs %>% dplyr::filter(FIX == 1) #fixers
nf <- rs %>% dplyr::filter(FIX == 0) #nonfixers
canopy <- rs %>% dplyr::filter(CCLCD2 < 4)
noncanopy <- rs %>% dplyr::filter(CCLCD2 >= 4)

evergreen <- rs %>% dplyr::filter(DEC0EVER1 == 1)
deciduous <- rs %>% dplyr::filter(DEC0EVER1 == 0)


# Growth: Plot binned data for individual trees -------------------------------------
# bin data for plotting
number_bins <- 50
stats_all_NCI_prop_young <- stats.bin(x = young$NCI_prop2,
                                      y = young[,"reci"],
                                      N = number_bins)
stats_N_NCI_prop_young <- stats.bin(x = young[young$FIX == 1,"NCI_prop2"],
                                    y = young[young$FIX == 1,"reci"],
                                    N = number_bins)
stats_non_NCI_prop_young <- stats.bin(x = young[young$FIX == 0,"NCI_prop2"],
                                      y = young[young$FIX == 0,"reci"],
                                      N = number_bins)

#Young:
png("output/young_recr.png", 
    width = 8, 
    height = 4, 
    units = 'in', 
    res = 300)
par(mar = c(5.1,5.1,4.1,2.1)+0.1,
    cex.main = 1.5,
    cex.lab = 1.5,
    cex.axis = 1.3,
    bty = "n")
#pdf(file="OR_gr.pdf",width=14)
#par(mar=c(5.1,5.1,4.1,2.1))
#par(mfrow=c(1,2))
plot(seq(1, number_bins, 1),
     stats_non_NCI_prop_young$stats["mean", ],
     xaxt = "n",
     xlab = "Crowding from N-Fixers (%)",
     ylab = expression(paste("Recruitment Probability (yr"^"-1",")")),
     main = "Young Trees (Stand Age < 60 Yrs)",
     ylim = c(0, 0.07),
     cex.lab = 1.5,
     col = "blue")
#mtext("c", 2, adj=5,las=1,padj=-12,cex=2)
points(seq(1, number_bins, 1),
       stats_N_NCI_prop_young$stats["mean", ],
       col = "red")
par(new = TRUE)

model_non_young <- lmer(reci ~ NCI_props +
                         NCIs +
                         NCIs^2 +
                         NCI_props*NCIs +
                         (1|PCN2), 
                        data = young[young$FIX == 0, ])
model_N_young <- lmer(reci ~ NCI_props +
                        NCIs +
                        NCIs^2 +
                        NCI_props * NCIs +
                        (1 | PCN2), 
                      data = young[young$FIX == 1, ])

plot(seq(0, 1, 0.1),
     stats_non_NCI_prop_young$stats["mean", ][[1]] +
       fixef(model_non_young)[2]*seq(0, 1, 0.1),
     type = "l",
     ylim = c(0, 0.07),
     xaxt = "n",
     yaxt = "n",
     ylab = "",
     xlab = "",
     col = "blue")

lines(seq(0, 1, 0.1),
      stats_N_NCI_prop_young$stats["mean", ][[1]] +
        fixef(model_N_young)[2]*seq(0, 1, 0.1),
      ylim = c(0, 0.07),
      xaxt = "n",
      yaxt = "n",
      ylab = "",
      xlab = "",
      col = "red")
axis(side = 1, at = seq(0, 1, 0.1), labels = seq(0, 100, 10))
legend("topleft", 
       legend = c("N-Fixer","Non-Fixer"), 
       col = c("red", "blue"),
       bty = "n",
       pch = c(1, 1))
dev.off()

orig <- stats_non_NCI_prop_young$stats["mean", ][[1]] + fixef(model_non_young)[2]*0
final <- stats_non_NCI_prop_young$stats["mean", ][[1]] + fixef(model_non_young)[2]*1
ch <- 100 * (final - orig) / orig
ch

orig <- stats_N_NCI_prop_young$stats["mean", ][[1]] + fixef(model_N_young)[2]*0
final <- stats_N_NCI_prop_young$stats["mean", ][[1]] + fixef(model_N_young)[2]*1
ch <- 100 * (final - orig) / orig
ch

#Old:
png("output/old_recr.png", 
    width = 8, 
    height = 4, 
    units = 'in', 
    res = 300)
par(mar = c(5.1,5.1,4.1,2.1)+0.1,
    cex.main = 1.5,
    cex.lab = 1.5,
    cex.axis = 1.3,
    bty = "n")
stats_all_NCI_prop_old <- stats.bin(x = old$NCI_prop2,
                                    y = old[ ,"reci"],
                                    N = number_bins)
stats_N_NCI_prop_old <- stats.bin(x = old[old$FIX == 1,"NCI_prop2"],
                                  y = old[old$FIX == 1,"reci"],
                                  N = number_bins)
stats_non_NCI_prop_old <- stats.bin(x = old[old$FIX == 0,"NCI_prop2"],
                                    y = old[old$FIX == 0,"reci"],
                                    N = number_bins)
plot(seq(1, number_bins, 1),
     stats_non_NCI_prop_old$stats["mean", ],
     xaxt = "n",
     xlab = "Crowding from N-Fixers (%)",
     ylab = expression(paste("Recruitment Probability (yr"^"-1",")")),
     main = "Old Trees (Stand Age >= 60 Yrs)",
     ylim = c(0, 0.07),
     cex.lab = 1.5,
     col = "blue")
#mtext("d", 2, adj=5,las=1,padj=-12,cex=2)
points(seq(1, number_bins, 1),
       stats_N_NCI_prop_old$stats["mean", ],
       col = "red")
par(new = TRUE)

oldnf <- old[old$FIX == 0, ]

model_non_old <- lmer(reci ~ NCI_props +
        NCIs +
        NCIs^2 +
        NCI_props * NCIs +
        (1 | PCN2), 
      data = oldnf)

oldf <- old[old$FIX == 1, ]
model_N_old <- lmer(reci ~ NCI_props +
                       NCIs +
                       NCIs^2 +
                       NCI_props * NCIs +
                       (1 | PCN2), 
                     data = oldf)

plot(seq(0, 1, 0.1),
     stats_non_NCI_prop_old$stats["mean", ][[1]] + fixef(model_non_old)[2]*seq(0, 1, 0.1),
     type = "l",
     ylim = c(0, 0.07),
     xaxt = "n",
     yaxt = "n",
     ylab = "",
     xlab = "",
     col = "blue")

lines(seq(0, 1, 0.1),
      stats_N_NCI_prop_old$stats["mean", ][[1]] + fixef(model_N_old)[2]*seq(0, 1, 0.1),
      ylim = c(0, 0.07),
      xaxt = "n",
      yaxt = "n",
      ylab = "",
      xlab = "",
      col = "red")
axis(side = 1, at = seq(0, 1, 0.1), labels = seq(0, 100, 10))
legend("topleft",
       legend = c("N-Fixer", "Non-Fixer"),
       col = c("red", "blue"),
       bty = "n",
       pch = c(1, 1))
dev.off()

orig <- stats_non_NCI_prop_old$stats["mean", ][[1]] + fixef(model_non_old)[2]*0
final <- stats_non_NCI_prop_old$stats["mean", ][[1]] + fixef(model_non_old)[2]*1
ch <- 100 * (final - orig) / orig
ch

orig <- stats_N_NCI_prop_old$stats["mean", ][[1]] + fixef(model_N_old)[2]*0
final <- stats_N_NCI_prop_old$stats["mean", ][[1]] + fixef(model_N_old)[2]*1
ch <- 100 * (final - orig) / orig
ch


#All:
number_bins <- 50
stats_all_NCI_prop_gs <- stats.bin(x = rs$NCI_prop2,
                                   y = rs[,"reci"],
                                   N = number_bins)
stats_N_NCI_prop_gs <- stats.bin(x = rs[rs$FIX == 1, "NCI_prop2"],
                                 y = rs[rs$FIX == 1, "reci"],
                                 N = number_bins)
stats_non_NCI_prop_gs <- stats.bin(x = rs[rs$FIX == 0, "NCI_prop2"],
                                   y = rs[rs$FIX == 0, "reci"],
                                   N = number_bins)

#pdf(file="OR_gr.pdf",width=14)
png("output/allforest_recr.png", width = 8, height = 4, units = 'in', res = 300)
par(mar = c(5.1,5.1,4.1,2.1)+0.1,
    cex.main = 1.5,
    cex.lab = 1.5,
    cex.axis = 1.3,
    bty = "n")
#par(mfrow=c(1,2))
#par(font.lab=2,las=1)
plot(seq(1, number_bins, 1),
     stats_non_NCI_prop_gs$stats["mean", ],
     xaxt = "n",
     xlab = "Crowding from N-Fixers (%)",
     ylab = expression(paste("Recruitment Probability (yr"^"-1",")")),
     main = "Forest",
     bty = "n",
     pch = 21,
     ylim = c(0, 0.07),
     cex = 1.5,
     col = "blue")
#mtext("c", 2, adj=5,las=1,padj=-12,cex=2)
points(seq(1, number_bins, 1),
       stats_N_NCI_prop_gs$stats["mean", ],
       col = "red",
       pch = 21,
       cex = 1.5)
par(new = TRUE)

model_non_gs <- lmer(reci ~ NCI_props +
                        NCIs +
                        NCIs^2 +
                        NCI_props * NCIs +
                        (1 | PCN2), 
                      data = rs[rs$FIX == 0, ])
                     
model_N_gs <- lmer(reci ~ NCI_props +
                      NCIs +
                      NCIs^2 +
                      NCI_props * NCIs +
                      (1 | PCN2), 
                    data = rs[rs$FIX == 1, ])

plot(seq(0, 1, 0.1),
     stats_non_NCI_prop_gs$stats["mean", ][[1]] + fixef(model_non_gs)[2]*seq(0, 1, 0.1),
     type = "l",
     ylim = c(0, 0.07),
     xaxt = "n",
     yaxt = "n",
     ylab = "",
     xlab = "",
     col = "blue",
     lwd = 2)

lines(seq(0, 1, 0.1),
      stats_N_NCI_prop_gs$stats["mean", ][[1]] + fixef(model_N_gs)[2]*seq(0, 1, 0.1),
      ylim = c(0, 0.07),
      xaxt = "n",
      yaxt = "n",
      ylab = "",
      xlab = "",
      col = "red",
      lwd = 2)
axis(side = 1, at = seq(0, 1, 0.1), labels = seq(0, 100, 10))
legend("topleft",
       legend = c("N-Fixer", "Non-Fixer"),
       col = c("red", "blue"),
       bty = "n",
       pch = c(1, 1))
dev.off()

orig <- stats_non_NCI_prop_gs$stats["mean", ][[1]] + fixef(model_non_gs)[2]*0
final <- stats_non_NCI_prop_gs$stats["mean", ][[1]] + fixef(model_non_gs)[2]*1
ch <- 100 * (final - orig) / orig
ch

orig<-stats_N_NCI_prop_gs$stats["mean",][[1]]+fixef(model_N_gs)[2]*0
final<-stats_N_NCI_prop_gs$stats["mean",][[1]]+fixef(model_N_gs)[2]*1
ch<-100*(final-orig)/orig
ch

#Canopy:
stats_all_NCI_prop_canopy <- stats.bin(x = canopy$NCI_prop2,
                                       y = canopy[,"reci"],
                                       N = number_bins)
stats_N_NCI_prop_canopy <- stats.bin(x = canopy[canopy$FIX==1,"NCI_prop2"],
                                     y = canopy[canopy$FIX==1,"reci"],
                                     N = number_bins)
stats_non_NCI_prop_canopy <- stats.bin(x = canopy[canopy$FIX==0,"NCI_prop2"],
                                       y = canopy[canopy$FIX==0,"reci"],
                                       N = number_bins)

png("output/canopy_recr.png", 
    width = 8, 
    height = 4, 
    units = 'in', 
    res = 300)
par(mar = c(5.1,5.1,4.1,2.1)+0.1,
    cex.main = 1.5,
    cex.lab = 1.5,
    cex.axis = 1.3,
    bty = "n")
plot(seq(1,number_bins,1), stats_non_NCI_prop_canopy$stats["mean",],
     xaxt = "n",
     xlab = "Crowding from N-Fixers (%)",
     ylab = expression(paste("Recruitment Probability (yr"^"-1",")")),
     main = "Canopy Trees",
     ylim = c(0, 0.07),
     cex.lab = 1.5, col = "blue")
#mtext("c", 2, adj=5,las=1,padj=-12,cex=2)
points(seq(1,number_bins,1), 
       stats_N_NCI_prop_canopy$stats["mean",], 
       col = "red")
par(new = TRUE)

model_non_canopy <- lmer(reci ~ NCI_props +
                            NCIs +
                            NCIs^2 +
                            NCI_props * NCIs +
                            (1 | PCN2), 
                          data = canopy[canopy$FIX==0,])

model_N_canopy <- lmer(reci ~ NCI_props +
                          NCIs +
                          NCIs^2 +
                          NCI_props * NCIs +
                          (1 | PCN2), 
                        data = canopy[canopy$FIX==1,])

plot(seq(0,1,0.1),stats_non_NCI_prop_canopy$stats["mean",][[1]]+fixef(model_non_canopy)[2]*seq(0,1,0.1),
     type = "l",
     ylim = c(0, 0.07),
     xaxt = "n", yaxt = "n",
     ylab = "", xlab = "",
     col = "blue")

lines(seq(0,1,0.1),stats_N_NCI_prop_canopy$stats["mean",][[1]]+fixef(model_N_canopy)[2]*seq(0,1,0.1),
      ylim = c(0, 0.07),
      xaxt = "n", yaxt = "n",
      ylab = "", xlab = "",
      col = "red")
axis(side = 1, 
     at = seq(0,1,0.1), 
     labels = seq(0,100,10))
legend("topleft", 
       legend = c("N-Fixer", "Non-Fixer"), 
       col = c("red", "blue"), 
       bty = "n",
       pch = c(1, 1))
dev.off()

orig <- stats_non_NCI_prop_canopy$stats["mean",][[1]]+fixef(model_non_canopy)[2]*0
final <- stats_non_NCI_prop_canopy$stats["mean",][[1]]+fixef(model_non_canopy)[2]*1
ch <- 100*(final-orig)/orig
ch

orig <- stats_N_NCI_prop_canopy$stats["mean",][[1]]+fixef(model_N_canopy)[2]*0
final <- stats_N_NCI_prop_canopy$stats["mean",][[1]]+fixef(model_N_canopy)[2]*1
ch <- 100*(final-orig)/orig
ch

#Non-canopy:
stats_all_NCI_prop_noncanopy <- stats.bin(x = noncanopy$NCI_prop2,
                                          y = noncanopy[,"reci"],
                                          N = number_bins)
stats_N_NCI_prop_noncanopy <- stats.bin(x = noncanopy[noncanopy$FIX==1,"NCI_prop2"],
                                        y = noncanopy[noncanopy$FIX==1,"reci"],
                                        N = number_bins)
stats_non_NCI_prop_noncanopy <- stats.bin(x = noncanopy[noncanopy$FIX==0,"NCI_prop2"],
                                          y = noncanopy[noncanopy$FIX==0,"reci"],
                                          N = number_bins)

png("output/noncanopy_recr.png", 
    width = 8, 
    height = 4, 
    units = 'in', 
    res = 300)
par(mar = c(5.1,5.1,4.1,2.1)+0.1,
    cex.main = 1.5,
    cex.lab = 1.5,
    cex.axis = 1.3,
    bty = "n")
plot(seq(1,number_bins,1), stats_non_NCI_prop_noncanopy$stats["mean",],
     xaxt = "n",
     xlab = "Crowding from N-Fixers (%)",
     ylab = expression(paste("Recruitment Probability (yr"^"-1",")")),
     main = "Noncanopy Trees",
     ylim = c(0, 0.07),
     cex.lab = 1.5,
     col = "blue")
#mtext("c", 2, adj=5,las=1,padj=-12,cex=2)
points(seq(1, number_bins, 1), 
       stats_N_NCI_prop_noncanopy$stats["mean",], 
       col = "red")
par(new = TRUE)

model_non_noncanopy <- lmer(reci ~ NCI_props + 
                              NCIs + 
                              NCIs^2 + 
                              NCI_props*NCIs +
                              (1|PCN2), 
                            data = noncanopy[noncanopy$FIX==0,])

model_N_noncanopy <- lmer(reci ~ NCI_props + 
                            NCIs + 
                            NCIs^2 + 
                            NCI_props*NCIs +
                            (1|PCN2), 
                          data = noncanopy[noncanopy$FIX==1,])

plot(seq(0, 1, 0.1),
     stats_non_NCI_prop_noncanopy$stats["mean",][[1]]+fixef(model_non_noncanopy)[2]*seq(0,1,0.1),
     type = "l",
     ylim = c(0, 0.07),
     xaxt = "n",
     yaxt = "n",
     ylab = "",
     xlab = "",
     col = "blue")
lines(seq(0, 1, 0.1),
      (stats_N_NCI_prop_noncanopy$stats["mean",][[1]]+fixef(model_N_noncanopy)[2]*seq(0,1,0.1)),
      ylim = c(0, 0.07),
      xaxt = "n",
      yaxt = "n",
      ylab = "",
      xlab = "",
      col = "red")
axis(side = 1, 
     at = seq(0, 1, 0.1), 
     labels = seq(0, 100, 10))
legend("topleft", 
       legend = c("N-Fixer", "Non-Fixer"), 
       col = c("red", "blue"), 
       bty = "n",
       pch = c(1, 1))
dev.off()

orig <- stats_non_NCI_prop_noncanopy$stats["mean",][[1]]+fixef(model_non_noncanopy)[2]*0
final <- stats_non_NCI_prop_noncanopy$stats["mean",][[1]]+fixef(model_non_noncanopy)[2]*1
ch <- 100*(final-orig)/orig
ch

orig <- stats_N_NCI_prop_noncanopy$stats["mean",][[1]]+fixef(model_N_noncanopy)[2]*0
final <- stats_N_NCI_prop_noncanopy$stats["mean",][[1]]+fixef(model_N_noncanopy)[2]*1
ch <- 100*(final-orig)/orig
ch

#Evergreen vs deciduous nonfixers:
# No evergreen fixers
stats_all_NCI_prop_evergreen <- stats.bin(x = evergreen$NCI_prop2,
                                          y = evergreen[ ,"reci"],
                                          N = number_bins)
stats_non_NCI_prop_evergreen <- stats.bin(x = evergreen[evergreen$FIX == 0,"NCI_prop2"],
                                          y = evergreen[evergreen$FIX == 0,"reci"],
                                          N = number_bins)
stats_all_NCI_prop_deciduous <- stats.bin(x = deciduous$NCI_prop2,
                                          y = deciduous[ ,"reci"],
                                          N = number_bins)
stats_non_NCI_prop_deciduous <- stats.bin(x = deciduous[deciduous$FIX == 0,"NCI_prop2"],
                                          y = deciduous[deciduous$FIX == 0,"reci"],
                                          N = number_bins)

png("output/evergreen_decid_recr.png", 
    width = 8, 
    height = 4, 
    units = 'in', 
    res = 300)
par(mar = c(5.1,5.1,4.1,2.1)+0.1,
    cex.main = 1.5,
    cex.lab = 1.5,
    cex.axis = 1.3,
    bty = "n")

plot(seq(1, number_bins, 1),
     stats_non_NCI_prop_deciduous$stats["mean", ],
     xaxt = "n",
     xlab = "Crowding from N-Fixers (%)",
     ylab = expression(paste("Relative Growth Rate (yr"^"-1",")")),
     main = "Deciduous Trees",
     ylim = c(0, 0.05),
     cex.lab = 1.5,
     col = "orange")

points(seq(1, number_bins, 1),
       stats_non_NCI_prop_evergreen$stats["mean", ],
       col = "green")

par(new = TRUE)

model_non_evergreen <- lmer(reci ~ NCI_props +
                              NCIs + 
                              NCIs^2 + 
                              NCI_props*NCIs +
                              (1|PCN2),  
                            data = evergreen[evergreen$FIX == 0, ])

model_non_deciduous <- lmer(reci ~ NCI_props + 
                              NCIs + 
                              NCIs^2 + 
                              NCI_props*NCIs +
                              (1|PCN2), 
                            data = deciduous[deciduous$FIX == 0, ])

plot(seq(0, 1, 0.1),
     stats_non_NCI_prop_deciduous$stats["mean", ][[1]] + fixef(model_non_deciduous)[2]*seq(0, 1, 0.1),
     type = "l",
     ylim = c(0, 0.05),
     xaxt = "n",
     yaxt = "n",
     ylab = "",
     xlab = "",
     col = "orange")

lines(seq(0, 1, 0.1),
      stats_non_NCI_prop_evergreen$stats["mean", ][[1]] + fixef(model_non_evergreen)[2]*seq(0, 1, 0.1),
      ylim = c(0, 0.05),
      xaxt = "n",
      yaxt = "n",
      ylab = "",
      xlab = "",
      col = "green")
axis(side = 1, at = seq(0, 1, 0.1), labels = seq(0, 100, 10))
legend("topleft",
       legend = c("Evergreen", "Deciduous"),
       col = c("green", "orange"),
       bty = "n",
       pch = c(1, 1))
dev.off()

f_pctch(model_non_evergreen)

f_pctch(model_non_deciduous)

# percent change plots ------------------------------
group <- c("forest","forest","young","young","old","old","canopy","canopy","noncanopy","noncanopy")
status <- c("N-fixer","non-fixer","N-fixer","non-fixer","N-fixer","non-fixer","N-fixer","non-fixer","N-fixer","non-fixer")
change <- c(f_pctch(model_N_gs)[[1]],
            f_pctch(model_non_gs)[[1]],
            f_pctch(model_N_young)[[1]],
            f_pctch(model_non_young)[[1]],
            f_pctch(model_N_old)[[1]],
            f_pctch(model_non_old)[[1]],
            f_pctch(model_N_canopy)[[1]],
            f_pctch(model_non_canopy)[[1]],
            f_pctch(model_N_noncanopy)[[1]],
            f_pctch(model_non_noncanopy)[[1]])
se <- c(f_pctch(model_N_gs)[[2]],
        f_pctch(model_non_gs)[[2]],
        f_pctch(model_N_young)[[2]],
        f_pctch(model_non_young)[[2]],
        f_pctch(model_N_old)[[2]],
        f_pctch(model_non_old)[[2]],
        f_pctch(model_N_canopy)[[2]],
        f_pctch(model_non_canopy)[[2]],
        f_pctch(model_N_noncanopy)[[2]],
        f_pctch(model_non_noncanopy)[[2]])
obs <- c(nobs(model_N_gs),
         nobs(model_non_gs),
         nobs(model_N_young),
         nobs(model_non_young),
         nobs(model_N_old),
         nobs(model_non_old),
         nobs(model_N_canopy),
         nobs(model_non_canopy),
         nobs(model_N_noncanopy),
         nobs(model_non_noncanopy))
changes <- as.data.frame(cbind(group,status,change,se,obs))
changes$change <- as.numeric(as.character(changes$change))
changes$se <- as.numeric(as.character(changes$se))
changes$obs <- as.numeric(as.character(changes$obs))
changes$group <- factor(changes$group, 
                        levels = c("forest","young","old","canopy","noncanopy"))
changes

library(cowplot)
library(ggthemes)
setwd("/Users/Anika/Documents/GradSchool/FIA_CompetitionModel/")
png("output/changes_recr.png", 
    width = 8, 
    height = 4, 
    units = 'in', 
    res = 300)
ggplot(data = changes, 
       aes(x = factor(group), y = change, fill = status, width = 0.65)) +
  geom_bar(stat = "identity", 
           position = position_dodge()) + 
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Group", y = "% Change", subtitle = "Recruitment") +
  geom_hline(yintercept = 0, color = "black") +
  geom_errorbar(aes(ymin = change-se, ymax = change+se), width = 0.2,
                position = position_dodge(width = 0.65)) 
  #geom_text(aes(label = paste("n = ", obs, sep="")), 
   #         position = position_dodge(width = 0.65), vjust = -1)
# theme(axis.text.x = element_text(colour="grey20",size=20,hjust=.5,vjust=.5,face="plain"),
#       axis.text.y = element_text(colour="grey20",size=20,hjust=.5,vjust=.5,face="plain"),
#       base_size=22)
dev.off()

##age and canopy side by side
group <- c("young","young","old","old")
status <- c("N-fixer","non-fixer","N-fixer","non-fixer")
change <- c(f_pctch(model_N_young)[[1]],
            f_pctch(model_non_young)[[1]],
            f_pctch(model_N_old)[[1]],
            f_pctch(model_non_old)[[1]])
se <- c(f_pctch(model_N_young)[[2]],
        f_pctch(model_non_young)[[2]],
        f_pctch(model_N_old)[[2]],
        f_pctch(model_non_old)[[2]])
changes_age <- as.data.frame(cbind(group, 
                                   status, 
                                   change,
                                   se))
changes_age$change <- as.numeric(as.character(changes_age$change))
changes_age$se <- as.numeric(as.character(changes_age$se))
changes_age$group <- factor(changes_age$group, levels = c("young","old"))
group <- c("canopy","canopy","noncanopy","noncanopy")
status <- c("N-fixer","non-fixer","N-fixer","non-fixer")
change <- c(f_pctch(model_N_canopy)[[1]],
            f_pctch(model_non_canopy)[[1]],
            f_pctch(model_N_noncanopy)[[1]],
            f_pctch(model_non_noncanopy)[[1]])
se <- c(f_pctch(model_N_canopy)[[2]],
        f_pctch(model_non_canopy)[[2]],
        f_pctch(model_N_noncanopy)[[2]],
        f_pctch(model_non_noncanopy)[[2]])
changes_can <- as.data.frame(cbind(group,status,change,se))
changes_can$change <- as.numeric(as.character(changes_can$change))
changes_can$se <- as.numeric(as.character(changes_can$se))
changes_can$group <- factor(changes_can$group, levels = c("canopy","noncanopy"))

png("output/changes1_recr.png", 
    width = 8, 
    height = 4, 
    units = 'in', 
    res = 300)
g1 <- ggplot(data = changes_age, 
             aes(x = factor(group), y = change, fill = status, width = 0.65)) +
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Group", y = "% Change", subtitle = "Recruitment") +
  geom_hline(yintercept = 0, color = "black") + 
  theme(legend.position = "none") +
  geom_errorbar(aes(ymin = change-se, ymax = change+se), width = 0.2,
                position = position_dodge(width = 0.65))

g2 <- ggplot(data = changes_can, 
             aes(x = factor(group), y = change, fill = status, width = 0.65)) +
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Group", y = "% Change") +
  geom_hline(yintercept = 0, color = "black") +
  geom_errorbar(aes(ymin = change-se, ymax = change+se), width = 0.2,
                position = position_dodge(width = 0.65))

plot_grid(g1, g2, labels=c("", ""), ncol = 2, nrow = 1)
dev.off()

group <- c("evergreen","deciduous")
status <- c("non-fixer","non-fixer")
change <- c(f_pctch(model_non_evergreen)[[1]],
            f_pctch(model_non_deciduous)[[1]])
se <- c(f_pctch(model_non_evergreen)[[2]],
        f_pctch(model_non_deciduous)[[2]])
changes_form <- as.data.frame(cbind(group, 
                                    status, 
                                    change,
                                    se))
changes_form$change <- as.numeric(as.character(changes_form$change))
changes_form$se <- as.numeric(as.character(changes_form$se))
changes_form$group <- factor(changes_form$group, levels = c("evergreen","deciduous"))

changes_form
saveRDS(changes_form, "output/changes_form_recr.RDS")


png("output/changes_form_recr.png", 
    width = 8, 
    height = 4, 
    units = 'in', 
    res = 300)
ggplot(data = changes_form, 
       aes(x = factor(group), y = change, fill = status, width = 0.65)) +
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Group", y = "% Change") +
  geom_hline(yintercept = 0, color = "black") + 
  theme(legend.position = "none") +
  geom_errorbar(aes(ymin = change-se, ymax = change+se), width = 0.2,
                position = position_dodge(width = 0.65))

#plot_grid(g1, g2, labels=c("", ""), ncol = 2, nrow = 1)
dev.off()




# mapping significance -------------------------------------

p <- fread("/Users/Anika/Documents/GradSchool/FIA_CompetitionModel/data/raw/PLOT.csv", 
           select = c("CN", "LAT", "LON"))

#gs <- merge(gs,p,by.x="PCN1.x",by.y="CN",all.x=T,all.y=F)
rs[is.na(rs$NCI_props), ]$NCI_props <- 0

p$llu <- p$LAT/100 + p$LON*100
colnames(p) <- c("PCN","lat1","lon1","llu")
p$PCN <- as.integer64(p$PCN)
pu <- p[!duplicated(p$llu), ]
# Number of plot record and unique plots
print(nrow(p)) #number of plots
print(nrow(pu)) #number of unique plots

# Get latitude and longitude bounds of AK, coterminous US, Mexico, PR, VI
latmin <- min(p$lat1, na.rm=TRUE)
latmax <- max(p$lat1, na.rm=TRUE)
lonmin <- min(p$lon1, na.rm=TRUE)
lonmax <- max(p$lon1, na.rm=TRUE)
latminint <- floor(latmin)
latmaxint <- ceiling(latmax)
lonminint <- floor(lonmin)
lonmaxint <- ceiling(lonmax)
dlat <- 1 # grid cell size
dlon <- 1 # grid cell size
lat.list <- seq(latminint,latmaxint,dlat)
lon.list <- seq(lonminint,lonmaxint,dlon)
nlat <- length(lat.list)
nlon <- length(lon.list)

## Make matrix of plot numbers
data <- rs #gs for all, young, old, canopy, under
xname <- "rs"

#sig.grid: alpha values for significance
#coef.grid: coefficient for NCI_props
#num.grid: number of N-fixers present
#nump.grid: number of plots
#riz.grid: average of ACT1vsRIZ0
#gen.grid: most common genus
#se.grid: standard error of regression for NCI_props
#pct.grid: ratio of coefficient for NCI_props to intercept

prefix <- c("coef.grid","num.grid","sig.grid","nump.grid","riz.grid","water.grid","se.grid")
grids <- paste(prefix, xname, sep=".")
coef.grid <- num.grid <- sig.grid <- nump.grid <- riz.grid <- water.grid <- se.grid <- gen.grid <- pct.grid <- array(dim=c(nlon,nlat))

for(i in 1:nlon){
  print(i/nlon*100)
  lon.a <- lon.list[i] - dlon/2
  lon.b <- lon.list[i] + dlon/2
  for(j in 1:nlat){
    lat.a <- lat.list[j] - dlat/2
    lat.b <- lat.list[j] + dlat/2
    
    plots <- p[(p$lon1>=lon.a & p$lon1<lon.b & p$lat1>=lat.a & p$lat<lat.b),]$PCN
    if(length(plots)>0){
      temp <- data[which(data$PCN2 %in% plots),
                   c("NCI_props","NCIs","reci","FIX","ACT1vsRIZ0","SPCD","PCN2","GENUS")]
      if(nrow(temp)>0){
        if(length(unique(temp$PCN2))>5 && nrow(temp[temp$reci>0,])>5){
          model.fixed2 = lmer(reci ~ NCI_props +
                                NCIs +
                                NCIs^2 +
                                NCI_props*NCIs +
                                (1|PCN2), 
                              data = temp)
          coefs <- data.frame(coef(summary(model.fixed2)))
          coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
          sig.grid[i,j] <- coefs["NCI_props","p.z"]
          #sig.grid[i,j] <- coefs[2,4]
          coef.grid[i,j] <- coefs["NCI_props","Estimate"]
          #coef.grid[i,j] <- summary(model.fixed2)$coefficients[2,1]
          num.grid[i,j] <- sum(temp$FIX, na.rm=T)
          #num.grid[i,j] <- colSums(temp)[[4]]
          nump.grid[i,j] <- length(unique(temp$PCN2))
          riz.grid[i,j] <- mean(temp$ACT1vsRIZ0, na.rm=T)
          gen.grid[i,j] <- tail(names(sort(table(temp$GENUS))), 1)
          #gen.grid[i,j] <- strtoi(names(sort(-table(temp$GENUS)))[1]) #most common genus
          se.grid[i,j] <- sqrt(diag(vcov(model.fixed2)))[[2]]
          pct.grid[i,j] <- 100*((coefs[1,"Estimate"]-(coefs[1,"Estimate"]+coefs["NCI_props","Estimate"]*1))/coefs[1,"Estimate"])
          #pct.grid[i,j] <- fixef(model.fixed2)[[2]]/fixef(model.fixed2)[[1]]
          
        }else{
          coef.grid[i,j] <- NA
          num.grid[i,j] <- NA
          riz.grid[i,j] <- NA
          se.grid[i,j] <- NA
          gen.grid[i,j] <- NA
          nump.grid[i,j] <- NA
          pct.grid[i,j] <- NA
        }
        
      }else{
        coef.grid[i,j] <- NA
        num.grid[i,j] <- NA
        riz.grid[i,j] <- NA
        se.grid[i,j] <- NA
        gen.grid[i,j] <- NA
        nump.grid[i,j] <- NA
        pct.grid[i,j] <- NA
      }
    }else{
      coef.grid[i,j] <- NA
      num.grid[i,j] <- NA
      riz.grid[i,j] <- NA
      se.grid[i,j] <- NA
      gen.grid[i,j] <- NA
      nump.grid[i,j] <- NA
      pct.grid[i,j] <- NA
      
    }
  }
}

#saveRDS(coef.grid, "/output/rs.coef.grid.RDS")
#saveRDS(se.grid, "/output/rs.se.grid.RDS")
#saveRDS(pct.grid, "/output/rs.pct.grid.RDS")

#log transform coefficients
coef.grid.t <- array(dim = c(nlon, nlat))
for(i in 1:nrow(coef.grid)){
  for(j in 1:ncol(coef.grid)){
    if(!is.na(coef.grid[i,j])){
      if(coef.grid[i,j] >= 0){
        coef.grid.t[i,j] <- sqrt(coef.grid[i,j])
      }
      else if(coef.grid[i,j] < 0){
        coef.grid.t[i,j] <- -sqrt(abs(coef.grid[i,j]))
      }
    }else{
      coef.grid.t[i,j] <- NA
    }
  }
}

#categorize coefficients
coef.grid.c <- array(dim = c(nlon, nlat))
for(i in 1:nrow(coef.grid)){
  for(j in 1:ncol(coef.grid)){
    if(!is.na(coef.grid[i,j])){
      if(coef.grid[i,j] > 0){
        coef.grid.c[i,j] <- 1
      }
      else if(coef.grid[i,j] < 0){
        coef.grid.c[i,j] <- -1
      }
    }else{
      coef.grid.c[i,j] <- NA
    }
  }
}

#remove outliers
coef.grid.noout <- array(dim = c(nlon, nlat))
for(i in 1:nrow(coef.grid)){
  for(j in 1:ncol(coef.grid)){
    if(!is.na(coef.grid[i,j])){
      if(coef.grid[i,j] > 1){
        coef.grid.noout[i,j] <- NA
      }
      else if(coef.grid[i,j]< (-1)){
        coef.grid.noout[i,j] <- NA
      }
      else{
        coef.grid.noout[i,j] <- coef.grid[i,j]
      }
    }else{
      coef.grid.noout[i,j] <- NA
    }
  }
}

#remove outliers
pct.grid.noout <- array(dim = c(nlon, nlat))
for(i in 1:nrow(pct.grid)){
  for(j in 1:ncol(pct.grid)){
    if(!is.na(pct.grid[i,j])){
      if(pct.grid[i,j] > 100){
        pct.grid.noout[i,j] <- NA
      }
      else if(pct.grid[i,j]< (-100)){
        pct.grid.noout[i,j] <- NA
      }
      else{
        pct.grid.noout[i,j] <- pct.grid[i,j]
      }
    }else{
      pct.grid.noout[i,j] <- NA
    }
  }
}

#significant only
sig.grid.c <- array(dim = c(nlon, nlat))
for(i in 1:nrow(coef.grid)){
  for(j in 1:ncol(coef.grid)){
    if(!is.na(coef.grid[i,j])){
      if(sig.grid[i,j] < 0.05){
        sig.grid.c[i,j] <- coef.grid[i,j]
        if(sig.grid.c[i,j] > 3){ #remove one weird outlier
          sig.grid.c[i,j] <- NA
        }
      }
    }else{
      sig.grid.c[i,j] <- NA
    }
  }
}

#saveRDS(sig.grid.c,"neighbor_effect_sig.RDS")
saveRDS(sig.grid.c,"rs.sig.grid.RDS")

df <- cbind(as.vector(sig.grid.c), as.vector(coef.grid), as.vector(coef.grid.t),
            as.vector(gen.grid), as.vector(riz.grid), as.vector(se.grid))
df <- as.data.frame(df)
colnames(df) <- c("siggridc","coefgrid","coefgrid.t","gengrid","rizgrid","segrid")
df_sig <- df[!is.na(sig.grid.c),]

nlev <- 64
#ticks <- 2^(0:13)
ticks <- c(min(coef.grid,na.rm=T),-2,-1,0,1,2,max(coef.grid,na.rm=T))
image.plot(lon.list, lat.list, coef.grid, nlevel=nlev, col=tim.colors(nlev),
           axis.args=list(at=ticks,labels=ticks),
           xlab="longitude",ylab="latitude",
           ylim = c(10, 60), xlim = c(-160, -50))
           #ylim=c(latminint-1,latmaxint+1),xlim=c(lonminint-2,lonmaxint+1))
title(main="Regression Slope -- All Trees")
#mtext("(a)",font=2,adj=0.01,padj=-0.5)
map('world',regions='usa',add=TRUE)

# coefficient grid with not outliers
require(RColorBrewer)
nHalf = sum(!is.na(coef.grid.noout))/2
Min = min(coef.grid.noout, na.rm=T)
Max = max(coef.grid.noout, na.rm=T)
Thresh = 0
## Make vector of colors for values below threshold
rc1 = colorRampPalette(colors = c("red", "yellow"), space="Lab")(nHalf)    
## Make vector of colors for values above threshold
rc2 = colorRampPalette(colors = c("grey", "blue"), space="Lab")(nHalf)
rampcols = c(rc1, rc2)
## In your example, this line sets the color for values between 49 and 51. 
rampcols[c(nHalf, nHalf+1)] = rgb(t(col2rgb("white")), maxColorValue=256) 
rb1 = seq(Min, Thresh, length.out=nHalf+1)
rb2 = seq(Thresh, Max, length.out=nHalf+1)[-1]
rampbreaks = c(rb1, rb2)
nlev <- 64
#ticks <- 2^(0:13)
ticks <- c(min(coef.grid.noout,na.rm=T), -2, -1, 0, 1, 2, max(coef.grid.noout,na.rm=T))
image.plot(lon.list, lat.list, coef.grid.noout, 
           nlevel = nlev, 
           col = rampcols, 
           breaks = rampbreaks,
           axis.args = list(at=ticks,labels=ticks),
           xlab = "longitude", ylab = "latitude",
           #ylim=c(14,33),xlim=c(-120,-80))
           ylim = c(10, 60), xlim = c(-160, -50))
#ylim=c(latminint-1,latmaxint+1),xlim=c(lonminint-2,lonmaxint+1))
title(main = "Recruitment Regression Slope, No outliers -- All Trees")
#mtext("(a)",font=2,adj=0.01,padj=-0.5)
map('world', regions = 'usa', add = TRUE)

# Sig grid categorized
require(RColorBrewer)
nHalf = sum(!is.na(sig.grid.c))/2
Min = min(sig.grid.c, na.rm=T)
Max = max(sig.grid.c, na.rm=T)
Thresh = 0
## Make vector of colors for values below threshold
rc1 = colorRampPalette(colors = c("red", "yellow"), space="Lab")(nHalf)    
## Make vector of colors for values above threshold
rc2 = colorRampPalette(colors = c("grey", "blue"), space="Lab")(nHalf)
rampcols = c(rc1, rc2)
## In your example, this line sets the color for values between 49 and 51. 
rampcols[c(nHalf, nHalf+1)] = rgb(t(col2rgb("white")), maxColorValue=256) 
rb1 = seq(Min, Thresh, length.out=nHalf+1)
rb2 = seq(Thresh, Max, length.out=nHalf+1)[-1]
rampbreaks = c(rb1, rb2)
nlev <- 64
ticks <- c(min(sig.grid.c,na.rm=T),0,max(sig.grid.c,na.rm=T))
#ticks <- c(-0.4,0,0.4)
image.plot(lon.list, lat.list, sig.grid.c, 
           nlevel = nlev, 
           col = rampcols, 
           breaks = rampbreaks,
           axis.args = list(at=ticks,labels=ticks),
           xlab = "longitude", ylab = "latitude",
           ylim = c(10, 60), xlim = c(-160, -50))
#ylim=c(latminint-1,latmaxint+1),xlim=c(lonminint-2,lonmaxint+1))
title(main = "Regression Slope (significant)")
#mtext("(a)",font=2,adj=0.01,padj=-0.5)
map('world', regions = 'usa', add = TRUE)

## Coefficient Grid (only significant)
require(RColorBrewer)
png("output/sigcoefs_recr.png", 
    width = 8, 
    height = 6, 
    units = 'in', 
    res = 300)
nHalf = sum(!is.na(sig.grid.c))/2
Min = min(sig.grid.c, na.rm=T)
Max = max(sig.grid.c, na.rm=T)
Thresh = 0
## Make vector of colors for values below threshold
rc1 = colorRampPalette(colors = c("red", "pink"), space="Lab")(nHalf)    
## Make vector of colors for values above threshold
rc2 = colorRampPalette(colors = c("light blue", "blue"), space="Lab")(nHalf)
rampcols = c(rc1, rc2)
## In your example, this line sets the color for values between 49 and 51. 
rampcols[c(nHalf, nHalf+1)] = rgb(t(col2rgb("white")), maxColorValue=256) 
rb1 = seq(Min, Thresh, length.out=nHalf+1)
rb2 = seq(Thresh, Max, length.out=nHalf+1)[-1]
rampbreaks = c(rb1, rb2)

nlev <- 64
ticks <- c(min(sig.grid.c, na.rm=T), 0, max(sig.grid.c, na.rm=T))
#ticks <- c(-0.4,0,0.4)
image.plot(lon.list, lat.list, sig.grid.c, 
           nlevel = nlev, 
           col = rampcols, 
           breaks = rampbreaks,
           axis.args = list(at=ticks,labels=ticks),
           xlab = "longitude",
           ylab = "latitude",
           ylim = c(23, 63), xlim = c(-140,-65))
title(main = "Regression Slope (significant) -- Recruitment")
#mtext("(a)",font=2,adj=0.01,padj=-0.5)
map('world', regions='usa', add=TRUE)
dev.off()

# trying another version of gen.grid
x <- as.vector(gen.grid)
x <- x[!is.na(x)]
nlev <- length(unique(x))
ticks <- unique(x)
library(pals)
polyc <- polychrome(nlev)
image.plot(lon.list, lat.list, gen.grid, nlevel=nlev, col=polyc,
           axis.args=list(at=ticks,labels=ticks),
           xlab="longitude",ylab="latitude",
           ylim=c(0,70),xlim=c(-140,-60))
#ylim=c(latminint-1,latmaxint+1),xlim=c(lonminint-2,lonmaxint+1))
title(main="Most Common Genus")
#mtext("(a)",font=2,adj=0.01,padj=-0.5)
map('world',regions='usa',add=TRUE)

# Pct change in recruitment rate (without outliers)
png("output/pctch_recr.png", 
    width = 8, 
    height = 6, 
    units = 'in', 
    res = 300)
require(RColorBrewer)
nHalf = sum(!is.na(pct.grid.noout))/2
Min = min(pct.grid.noout, na.rm=T)
Max = max(pct.grid.noout, na.rm=T)
Thresh = 0
## Make vector of colors for values below threshold
rc1 = colorRampPalette(colors = c("red", "grey95"), space="Lab")(nHalf)    
## Make vector of colors for values above threshold
rc2 = colorRampPalette(colors = c("grey95", "blue"), space="Lab")(nHalf)
rampcols = c(rc1, rc2)
## In your example, this line sets the color for values between 49 and 51. 
rampcols[c(nHalf, nHalf+1)] = rgb(t(col2rgb("white")), maxColorValue=256) 
rb1 = seq(Min, Thresh, length.out=nHalf+1)
rb2 = seq(Thresh, Max, length.out=nHalf+1)[-1]
rampbreaks = c(rb1, rb2)
nlev <- 64
ticks <- c(-100,0,100)
ticks <- c(min(pct.grid.noout,na.rm=T), 0, max(pct.grid.noout,na.rm=T))
ticks <- c(-80, -60, -40, -20, 0, 20, 40, 60, 80)
image.plot(lon.list, lat.list, pct.grid.noout, 
           nlevel = nlev, 
           col = rampcols, 
           breaks = rampbreaks,
           axis.args = list(at=ticks,labels=ticks),
           xlab = "longitude", ylab = "latitude",
           ylim = c(23,63), xlim = c(-140,-65))
title(main = "Percent Change in Recruitment Rate")
map('world', regions = 'usa', add = TRUE)
dev.off()

#pct change in recruitment rate (with outliers)
nHalf = sum(!is.na(pct.grid))/2
Min = min(pct.grid, na.rm=T)
Max = max(pct.grid, na.rm=T)
Thresh = 0
## Make vector of colors for values below threshold
rc1 = colorRampPalette(colors = c("red", "grey95"), space="Lab")(nHalf)    
## Make vector of colors for values above threshold
rc2 = colorRampPalette(colors = c("grey95", "blue"), space="Lab")(nHalf)
rampcols = c(rc1, rc2)
## In your example, this line sets the color for values between 49 and 51. 
rampcols[c(nHalf, nHalf+1)] = rgb(t(col2rgb("white")), maxColorValue=256) 
rb1 = seq(Min, Thresh, length.out=nHalf+1)
rb2 = seq(Thresh, Max, length.out=nHalf+1)[-1]
rampbreaks = c(rb1, rb2)
nlev <- 64
ticks <- c(min(pct.grid,na.rm=T), 0, max(pct.grid,na.rm=T))
image.plot(lon.list, lat.list, pct.grid, 
           nlevel = nlev, 
           col = rampcols, 
           breaks = rampbreaks,
           axis.args = list(at=ticks,labels=ticks),
           xlab = "longitude", ylab = "latitude",
           ylim = c(23,63), xlim = c(-140,-65))
title(main = "Percent Change in Recruitment Rate")
map('world', regions = 'usa', add = TRUE)

# Coefficient grid categorized + or -
nlev <- 2
#ticks <- 2^(0:13)
ticks <- c(min(coef.grid.c,na.rm=T), max(coef.grid.c,na.rm=T))
image.plot(lon.list, lat.list, coef.grid.c, 
           nlevel = nlev, 
           col = tim.colors(nlev),
           axis.args = list(at=ticks,labels=ticks),
           xlab = "longitude", ylab = "latitude",
           ylim = c(23,63), xlim = c(-140,-65))
title(main = "Recruitment") #Regression slope sign -- all trees
map('world', regions = 'usa', add = TRUE)

# Number of plots in each grid cell
nlev <- 64
ticks <- 2^(0:13)
image.plot(lon.list, lat.list, log(nump.grid), 
           nlevel = nlev, 
           col = tim.colors(nlev),
           axis.args = list(at=log(ticks),labels=ticks),
           xlab = "longitude", ylab = "latitude",
           ylim = c(23,63), xlim = c(-140,-65))
title(main = "Number of Plots")
map('world', regions = 'usa', add = TRUE)

# Average of rhizobial vs actinorhizal
nlev <- 64
ticks <- c(1, exp(0.5), exp(1))
image.plot(lon.list, lat.list, riz.grid, 
           nlevel = nlev, 
           col = tim.colors(nlev),
           axis.args = list(at=log(ticks),labels=log(ticks)),
           xlab = "longitude", ylab = "latitude",
           #ylim=c(14,33),xlim=c(-120,-80))
           ylim = c(23,63), xlim = c(-140,-65))
title(main = "Rhizobial 0, Actinorhizal 1 -- All Trees")
map('world', regions = 'usa', add = TRUE)

# Number of fixers in each grid cell
nlev <- 64
ticks <- 2^(0:13)
image.plot(lon.list, lat.list, log(num.grid+0.01), 
           nlevel = nlev, 
           col = tim.colors(nlev),
           axis.args = list(at=log(ticks),labels=ticks),
           xlab = "longitude", ylab = "latitude",
           ylim = c(23,63), xlim = c(-160,-65))
title(main = "Number of Fixers")
map('world', regions = 'usa', add = TRUE)

# Standard error map
png("output/semap_recr.png", 
    width = 8, 
    height = 6, 
    units = 'in', 
    res = 300)
nlev <- 64
ticks <- c(min(se.grid,na.rm=T), 0.001, 0.01, max(se.grid,na.rm=T))
image.plot(lon.list, lat.list, log(se.grid), 
           nlevel = nlev, 
           col = grey((256:0)/256),
           axis.args = list(at=log(ticks),labels=ticks),
           xlab = "longitude", ylab = "latitude",
           ylim = c(23,63), xlim = c(-160,-65))
title(main = "Standard Error on Regression Coefficient -- Recruitment")
map('world', regions = 'usa', add = TRUE)
dev.off()


#dominant genus from models_ind_mort.R


###
coef <- as.vector(coef.grid)
nplot <- as.vector(nump.grid)
nfix <- as.vector(num.grid)
rhiz <- as.vector(riz.grid)
water <- as.vector(water.grid)
se <- as.vector(se.grid)
pct <- as.vector(pct.grid)

cor.test(coef, nplot, use = "complete.obs")
cor.test(coef, nfix, use = "complete.obs")
cor.test(coef, rhiz, use = "complete.obs")
cor.test(nplot, se, use = "complete.obs")
cor.test(nfix, se, use = "complete.obs")

cor.test(coef, se, use = "complete.obs")
cor.test(pct, se, use = "complete.obs")
cor.test(coef, pct, use = "complete.obs")

plot(nplot,se, xlab="Number of Plots in Grid Cell", ylab="Standard Error", main="Standard Error vs. Number of Plots")
plot(nfix,se, xlab="Number of Fixers in Grid Cell", ylab="Standard Error", main="Standard Error vs. Number of Fixers")

# Histogram of coefficients
x <- as.vector(sig.grid.c)
x <- x[!is.na(x)]
hist(x, prob=T, breaks=30, xlab="B1 for Grid Cell", 
     main="Histogram of Regression coefficients (significant only)")

a1 <- table(coef.grid.c)
aaf <- (a1[[2]]/(a1[[1]]+a1[[2]]))*100
aaf <- round(aaf, digits=2)
aaf

x <- as.vector(sig.grid.c)
x <- x[!is.na(x)]
length(x[x>0])/length(x)

# Check correlation with USDA trait data ----------------------------------

usda <- readRDS("data/derived/usda_numeric.RDS")

mdat <- merge(rs, usda, by.x="SPCD", by.y="SPCD", all.x=T, all.y=F)

mdat <- mdat[!is.na(mdat$CNRatio),]
mdat_small <- mdat[,c(1,4,5,8,12:21,23,30,40:66)]

library(Hmisc)
res2 <- rcorr(as.matrix(mdat_small))
source("scripts/flattenCorrMatrix.R")
flat <- flattenCorrMatrix(res2$r, res2$P)
View(flat)

strong <- flat[flat$cor > 0.2, ]
strong

flat[flat$row=="reci", ]
flat[flat$row=="NCI_prop2", ]

library("PerformanceAnalytics")
my_data <- mtcars[, c(1,3,4,5,6,7)]
chart.Correlation(my_data, histogram=TRUE, pch=19)


# Probe deciduous non-fixers with USDA -----------------------------
usda <- readRDS("data/derived/usda_numeric.RDS")

gd <- merge(deciduous, usda, by.x="SPCD", by.y="SPCD", all.x=T, all.y=F)

gdnf <- gd[gd$FIX == 0,]

grmod <- function(datapass) {
  model <- lmer(reci ~ NCI_props + 
                  NCIs + 
                  NCIs^2 + 
                  NCI_props*NCIs +
                  (1|PCN2),  
                data = datapass)
  return(model)
}

model_non <- grmod(datapass = gdnf[gdnf$SLOPE <= 8, ])
lowslope <- f_pctch(model_non)
lowslope

model_non <- grmod(datapass = gdnf[gdnf$SLOPE > 8, ])
highslope <- f_pctch(model_non)
highslope

model_non <- grmod(datapass = gdnf[gdnf$ASPECT >= 90 & gdnf$ASPECT <= 270, ])
south <- f_pctch(model_non)
south

model_non <- grmod(datapass = gdnf[gdnf$ASPECT < 90 | gdnf$ASPECT > 270, ])
north <- f_pctch(model_non)
north

model_non <- grmod(datapass = gdnf[gdnf$ELEVm <= 277, ])
lowelev <- f_pctch(model_non)
lowelev

model_non <- grmod(datapass = gdnf[gdnf$ELEVm > 277, ])
highelev <- f_pctch(model_non)
highelev

model_non <- grmod(datapass = gdnf[!is.na(gdnf$CNRatio) & gdnf$CNRatio == 1,])
highcn <- f_pctch(model_non)
highcn

model_non <- grmod(datapass = gdnf[!is.na(gdnf$CNRatio) & gdnf$CNRatio > 1,])
lowcn <- f_pctch(model_non)
lowcn

model_non <- grmod(datapass = gdnf[gdnf$FoliagePorosity == 1, ])
highporosity <- f_pctch(model_non)
highporosity

model_non <- grmod(datapass = gdnf[gdnf$FoliagePorosity > 1, ])
lowporosity <- f_pctch(model_non)
lowporosity

model_non <- grmod(datapass = gdnf[!is.na(gdnf$GrowthRate) & gdnf$GrowthRate == 1,])
highgr <- f_pctch(model_non)
highgr

model_non <- grmod(datapass = gdnf[!is.na(gdnf$GrowthRate) & gdnf$GrowthRate == 2,])
lowgr <- f_pctch(model_non)
lowgr

model_non <- grmod(datapass = gdnf[!is.na(gdnf$HeightMatureFt) & gdnf$HeightMatureFt <= 75,])
tall <- f_pctch(model_non)
tall

model_non <- grmod(datapass = gdnf[!is.na(gdnf$HeightMatureFt) & gdnf$HeightMatureFt > 75,])
short <- f_pctch(model_non)
short

model_non <- grmod(datapass = gdnf[!is.na(gdnf$Allelopath) & gdnf$Allelopath == 1,])
alleloy <- f_pctch(model_non)
alleloy

model_non <- grmod(datapass = gdnf[!is.na(gdnf$Allelopath) & gdnf$Allelopath == 0,])
allelon <- f_pctch(model_non)
allelon

model_non <- grmod(datapass = gdnf[!is.na(gdnf$LeafRetention) & gdnf$LeafRetention == 1,])
leafry <- f_pctch(model_non)
leafry

model_non <- grmod(datapass = gdnf[!is.na(gdnf$LeafRetention) & gdnf$LeafRetention == 0,])
leafrn <- f_pctch(model_non)
leafrn

model_non <- grmod(datapass = gdnf[!is.na(gdnf$Lifespan) & gdnf$Lifespan == 1,])
longlife <- f_pctch(model_non)
longlife

model_non <- grmod(datapass = gdnf[!is.na(gdnf$Lifespan) & gdnf$Lifespan > 1,])
shortlife <- f_pctch(model_non)
shortlife

model_non <- grmod(datapass = gdnf[!is.na(gdnf$Lifespan) & gdnf$Lifespan == 1,])
longlife <- f_pctch(model_non)
longlife

model_non <- grmod(datapass = gdnf[!is.na(gdnf$CaCO3Tolerance) & gdnf$CaCO3Tolerance <= 2,])
highcatol <- f_pctch(model_non)
highcatol

model_non <- grmod(datapass = gdnf[!is.na(gdnf$CaCO3Tolerance) & gdnf$CaCO3Tolerance > 2,])
lowcatol <- f_pctch(model_non)
lowcatol

model_non <- grmod(datapass = gdnf[!is.na(gdnf$DroughtTolerance) & gdnf$DroughtTolerance <= 2,])
highdroughttol <- f_pctch(model_non)
highdroughttol

model_non <- grmod(datapass = gdnf[!is.na(gdnf$DroughtTolerance) & gdnf$DroughtTolerance > 2,])
lowdroughttol <- f_pctch(model_non)
lowdroughttol

model_non <- grmod(datapass = gdnf[!is.na(gdnf$FerilityReq) & gdnf$FerilityReq == 1,])
highfert <- f_pctch(model_non)
highfert

model_non <- grmod(datapass = gdnf[!is.na(gdnf$FerilityReq) & gdnf$FerilityReq > 1,])
lowfert <- f_pctch(model_non)
lowfert

model_non <- grmod(datapass = gdnf[!is.na(gdnf$MaxRootDepthIn) & gdnf$MaxRootDepthIn > 32,])
longroot <- f_pctch(model_non)
longroot

model_non <- grmod(datapass = gdnf[!is.na(gdnf$MaxRootDepthIn) & gdnf$MaxRootDepthIn <= 32,])
shortroot <- f_pctch(model_non)
shortroot

model_non <- grmod(datapass = gdnf[!is.na(gdnf$ShadeTolerance) & gdnf$ShadeTolerance == 1,])
highshadetol <- f_pctch(model_non)
highshadetol

model_non <- grmod(datapass = gdnf[!is.na(gdnf$ShadeTolerance) & gdnf$ShadeTolerance > 1,])
lowshadetol <- f_pctch(model_non)
lowshadetol

model_non <- grmod(datapass = gdnf[!is.na(gdnf$SeedAbundance) & gdnf$SeedAbundance == 1,])
highseed <- f_pctch(model_non)
highseed

model_non <- grmod(datapass = gdnf[!is.na(gdnf$SeedAbundance) & gdnf$SeedAbundance > 1,])
lowseed <- f_pctch(model_non)
lowseed

group <- c("slope","slope",
           "aspect","aspect",
           "elevation","elevation",
           "CN ratio", "CN ratio",
           "foliar porosity","foliar porosity",
           "growth rate","growth rate",
           "mature height","mature height",
           "allelopathy","allelopathy",
           "leaf retention","leaf retention",
           "life span", "life span",
           "CaCO3 tolerance", "CaCO3 tolerance",
           "drought tolerance","drought tolerance",
           "fertility req","fertility req",
           "root depth","root depth",
           "shade tolerance","shade tolerance",
           "seed abundance","seed abundance")
status <- c("high","low","high","low","high","low","high","low","high","low","high","low",
            "high","low","high","low","high","low","high","low","high","low","high","low",
            "high","low","high","low","high","low","high","low")
change <- c(highslope[[1]],lowslope[[1]],
            north[[1]],south[[1]],
            highelev[[1]],lowelev[[1]],
            highcn[[1]],lowcn[[1]],
            highporosity[[1]],lowporosity[[1]],
            highgr[[1]],lowgr[[1]],
            tall[[1]],short[[1]],
            alleloy[[1]],allelon[[1]],
            leafry[[1]],leafrn[[1]],
            longlife[[1]],shortlife[[1]],
            highcatol[[1]],lowcatol[[1]],
            highdroughttol[[1]],lowdroughttol[[1]],
            highfert[[1]],lowfert[[1]],
            longroot[[1]],shortroot[[1]],
            highshadetol[[1]],lowshadetol[[1]],
            highseed[[1]],lowseed[[1]])
se <- c(highslope[[2]],lowslope[[2]],
        north[[2]],south[[2]],
        highelev[[2]],lowelev[[2]],
        highcn[[2]],lowcn[[2]],
        highporosity[[2]],lowporosity[[2]],
        highgr[[2]],lowgr[[2]],
        tall[[2]],short[[2]],
        alleloy[[2]],allelon[[2]],
        leafry[[2]],leafrn[[2]],
        longlife[[2]],shortlife[[2]],
        highcatol[[2]],lowcatol[[2]],
        highdroughttol[[2]],lowdroughttol[[2]],
        highfert[[2]],lowfert[[2]],
        longroot[[2]],shortroot[[2]],
        highshadetol[[2]],lowshadetol[[2]],
        highseed[[2]],lowseed[[2]])
changes_decid <- as.data.frame(cbind(group, 
                                     status, 
                                     change,
                                     se))
changes_decid$change <- as.numeric(as.character(changes_decid$change))
changes_decid$se <- as.numeric(as.character(changes_decid$se))
changes_decid$group <- factor(changes_decid$group, 
                              levels = c("slope", "aspect", "elevation", "CN ratio", 
                                         "foliar porosity", "growth rate", "mature height",
                                         "allelopathy", "leaf retention", "life span", 
                                         "CaCO3 tolerance", "drought tolerance", "fertility req",
                                         "root depth", "shade tolerance", "seed abundance"))
changes_decid
saveRDS(changes_decid, "output/changes_non_decid_recr.RDS")

png("output/changes_non_decid_recr.png", 
    width = 8, 
    height = 4, 
    units = 'in', 
    res = 300)
ggplot(data = changes_decid, 
       aes(x = factor(group), y = change, fill = status, width = 0.65)) +
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Group", y = "% Change") +
  geom_hline(yintercept = 0, color = "black") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_errorbar(aes(ymin = change-se, ymax = change+se), width = 0.2,
                position = position_dodge(width = 0.65))

#plot_grid(g1, g2, labels=c("", ""), ncol = 2, nrow = 1)
dev.off()


# N dep --------------------------------------
require(raster)
Ndep <- raster("/Users/Anika/Documents/GradSchool/FIA_EstimateProj/output/N_dep.grd")

coords <- rs[,c("LON","LAT")]
dep <- raster::extract(Ndep, coords)

rs <- cbind(rs, dep)

# n dep in bins of 2
model1 <- grmod(datapass = rs[rs$dep > 0 & rs$dep < 2,])
model2 <- grmod(datapass = rs[rs$dep >= 2 & rs$dep < 4,])
model3 <- grmod(datapass = rs[rs$dep >= 4 & rs$dep < 6,])
model4 <- grmod(datapass = rs[rs$dep >= 6 & rs$dep < 8,])
model5 <- grmod(datapass = rs[rs$dep >= 8,])

N1 <- f_pctch(grmod(datapass = rs[rs$dep > 0 & rs$dep < 1,]))
N2 <- f_pctch(grmod(datapass = rs[rs$dep >= 1 & rs$dep < 2,]))
N3 <- f_pctch(grmod(datapass = rs[rs$dep >= 2 & rs$dep < 3,]))
N4 <- f_pctch(grmod(datapass = rs[rs$dep >= 3 & rs$dep < 4,]))
N5 <- f_pctch(grmod(datapass = rs[rs$dep >= 4 & rs$dep < 5,]))
N6 <- f_pctch(grmod(datapass = rs[rs$dep >= 5 & rs$dep < 6,]))
N7 <- f_pctch(grmod(datapass = rs[rs$dep >= 6 & rs$dep < 7,]))
N8 <- f_pctch(grmod(datapass = rs[rs$dep >= 7 & rs$dep < 8,]))
N9 <- f_pctch(grmod(datapass = rs[rs$dep >= 8 & rs$dep < 9,]))
N10 <- f_pctch(grmod(datapass = rs[rs$dep >= 9,]))

#n dep in quartiles
# test <- rs %>% mutate(quartile = ntile(dep, 4))
# 
# model1 <- grmod(datapass = test[test$quartile == 1, ])
# model2 <- grmod(datapass = test[test$quartile == 2, ])
# model3 <- grmod(datapass = test[test$quartile == 3, ])
# model4 <- grmod(datapass = test[test$quartile == 4, ])
# 
# N1 <- f_pctch(model1)
# N2 <- f_pctch(model2)
# N3 <- f_pctch(model3)
# N4 <- f_pctch(model4)

ndep <- rbind.data.frame(N1, N2, N3, N4, N5, N6, N7, N8, N9, N10)
colnames(ndep) <- c("est", "se", "num")
quad <- c("0-1", "1-2", "2-3", "3-4", "4-5", "5-6", "6-7", "7-8", "8-9", "9+")
ndep <- cbind(ndep, quad)

ggplot(data = ndep, 
       aes(x = factor(quad), y = est, width = 0.65)) +
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Group", y = "% Change", title = "% Change in Recruitment by N deposition") +
  geom_hline(yintercept = 0, color = "black") + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_errorbar(aes(ymin = est-se, ymax = est+se), width = 0.2,
                position = position_dodge(width = 0.65))


ndep_recr <- ndep
