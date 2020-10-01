

####################
#
#
# models_plot.R
#
# 1/17/19
# Anika Petach
# 
#
#####################

require(bit64)
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
require(fields)

setwd("/Users/Anika/Documents/GradSchool/FIA_CompetitionModel/")

source("scripts/functions.R")


# read data for models ------------------------------
FIA_plot <- readRDS("output/fia_plot_foranalysis.RDS")

#FIA_plot$survyr <- 1 - FIA_plot$mortyr
#FIA_plot$survyr <- (FIA_plot$nt1 - FIA_plot$died) / (FIA_plot$nt1 * (FIA_plot$mt2 - FIA_plot$mt1))

# model development ---------------------------------------------------
# Growth (biomass accumulation):
model_growth1 <- lmer(BA_change ~ scale(BA_Nfixer_prop_total) + 
                       scale(BA_total1) + 
                       scale(BA_total2) + 
                       scale(BA_Nfixer1) + 
                       scale(BA_nonfixer1) +
                       (1 | STATECD), 
                     data = FIA_plot)

drop1(model_growth1) #same AIC if drop BA_Nfixer1 and/or BA_nonfixer1

# Mortality:
model_mort1 <- lmer(mortyr ~ scale(BA_Nfixer_prop_total) + 
                    scale(BA_total1) + 
                    scale(BA_total2) + 
                    scale(BA_Nfixer1) + 
                    scale(BA_nonfixer1) +
                    (1 | STATECD), 
                  data = FIA_plot) 
# AIC = 11290123
drop1(model_mort1) #same AIC if drop BA_Nfixer1 and/or BA_nonfixer1

model_mort2 <- lmer(mortyr ~ scale(BA_Nfixer_prop_total) + 
                      scale(BA_total1) + 
                      scale(BA_total2) + 
                      scale(BA_Nfixer2) + 
                      scale(BA_nonfixer2) +
                      (1 | STATECD), 
                    data = FIA_plot) 
AIC(model_mort2) #AIC =  9532404, better model

model_mort3 <- lmer(mortyr ~ scale(BA_Nfixer_prop_total) + 
                      scale(BA_Nfixer_prop_total)^2 + 
                      scale(BA_total2) + 
                      scale(BA_total1) +
                      (1 | STATECD), 
                    data = FIA_plot) 
AIC(model_mort3) #AIC = 9537306, model2 still better

# Recruitment:
model_recr1 <- lmer(recyr ~ scale(BA_Nfixer_prop_total) + 
                     scale(BA_total1) + 
                     scale(BA_total2) + 
                     scale(BA_Nfixer1) + 
                     scale(BA_nonfixer1) +
                     (1 | STATECD), 
                   data = FIA_plot)
drop1(model_recr1) # AIC = 8011245

# plot BA_change vs BA_Nfixer_prop_total ----------------------------
number_bins <- 50

stats_all_NCI_prop_all <- stats.bin(x = FIA_plot$BA_Nfixer_prop_total,
                                    y = FIA_plot[,"BA_change"],
                                    N = number_bins)

model_all <- lmer(BA_change ~ scale(BA_Nfixer_prop_total) + 
                    scale(BA_total1) + 
                    scale(BA_total2) + 
                    scale(BA_Nfixer1) + 
                    scale(BA_nonfixer1) +
                    (1 | STATECD), 
                  data = FIA_plot) 
#par(mfrow=c(1,1))
#png("allplot.png", width = 10, height = 6, units = 'in', res = 300)
png("output/allforest_plot_growth.png", 
    width = 10, 
    height = 6, 
    units = 'in', 
    res = 300)
par(mar = c(5.1, 5.1, 4.1, 2.1) + 0.1,
    cex.main = 1.5,
    cex.lab = 1.5,
    cex.axis = 1.3,
    bty = "n")
plot(seq(1, number_bins, 1),
     stats_all_NCI_prop_all$stats["mean",],
     xaxt = "n",
     xlab = "Basal Area from N-Fixers (%)",
     ylab = "Change in Basal Area (%)",
     main = "Forest",
     pch = 21, 
     bg = "grey", 
     cex = 2,
     ylim = c(-2, 4),
     cex.lab = 1.5)
#mtext("d", 2, adj=5,las=1,padj=-12,cex=2)
par(new = TRUE)

plot(seq(0, 1, 0.1),
     fixef(model_all)[1]+fixef(model_all)[2]*seq(0,1,0.1),
     type = "l",
     lwd = 2,
     ylim = c(-2, 4),
     xaxt = "n",
     yaxt = "n",
     ylab = "",
     xlab = "",
     col = "black")
axis(side = 1,
     at = seq(0, 1, 0.1),
     labels = seq(0, 100, 10))
dev.off()

f_pctch_p(model_all)


# rec/year vs BA_Nfixer_prop1
stats_all_NCI_prop_all <- stats.bin(x = FIA_plot$BA_Nfixer_prop1,
                                    y = FIA_plot[,"recyr"],
                                    N = number_bins)

model_all <- lmer(recyr ~ scale(BA_Nfixer_prop_total) + 
                    scale(BA_total1) + 
                    scale(BA_total2) + 
                    scale(BA_Nfixer1) + 
                    scale(BA_nonfixer1) +
                    (1|STATECD), data=FIA_plot) #add random effect for state
#par(mfrow=c(1,1))
#setwd("/Users/Anika/Documents/GradSchool/FIA_CompetitionModel/")
#png("allplot.png", width = 10, height = 6, units = 'in', res = 300)
png("output/allforest_plot_recr.png", 
    width = 10, 
    height = 6, 
    units = 'in', 
    res = 300)
par(mar = c(5.1,5.1,4.1,2.1)+0.1,
    cex.main = 1.5,
    cex.lab = 1.5,
    cex.axis = 1.3,
    bty = "n")
plot(seq(1, number_bins, 1),
     stats_all_NCI_prop_all$stats["mean",],
     xaxt = "n",
     xlab = "Basal Area from N-Fixers (%)",
     ylab = "Recruitment (individuals/plot/yr)",
     main = "Forest",
     pch = 21, 
     bg = "grey", 
     cex = 2,
     ylim = c(-0.06, 0.06),
     cex.lab = 1.5)
#mtext("d", 2, adj=5,las=1,padj=-12,cex=2)
par(new = TRUE)

plot(seq(0, 1, 0.1),
     fixef(model_all)[1]+fixef(model_all)[2]*seq(0,1,0.1),
     type = "l",
     lwd = 2,
     ylim = c(-0.06, 0.06),
     xaxt = "n",
     yaxt = "n",
     ylab = "",
     xlab = "",
     col = "black")
axis(side = 1,
     at = seq(0, 1, 0.1),
     labels = seq(0, 100, 10))
dev.off()

f_pctch_p(model_all)

# mort/year vs BA_Nfixer_prop1
stats_all_NCI_prop_all <- stats.bin(x = FIA_plot$BA_Nfixer_prop1,
                                    y = FIA_plot[,"mortyr"],
                                    N = number_bins)

model_all <- lmer(mortyr ~ scale(BA_Nfixer_prop_total) + 
                    scale(BA_total1) + 
                    scale(BA_total2) + 
                    scale(BA_Nfixer1) + 
                    scale(BA_nonfixer1) +
                    (1 | STATECD), 
                  data = FIA_plot) #add random effect for state
#par(mfrow=c(1,1))
#setwd("/Users/Anika/Documents/GradSchool/FIA_CompetitionModel/")
#png("allplot.png", width = 10, height = 6, units = 'in', res = 300)
png("output/allforest_plot_mort.png", 
    width = 10, 
    height = 6, 
    units = 'in', 
    res = 300)
par(mar = c(5.1,5.1,4.1,2.1)+0.1,
    cex.main = 1.5,
    cex.lab = 1.5,
    cex.axis = 1.3,
    bty = "n")
plot(seq(1, number_bins, 1),
     stats_all_NCI_prop_all$stats["mean", ],
     xaxt = "n",
     xlab = "Basal Area from N-Fixers (%)",
     ylab = "Mortality (deaths/plot/year)",
     main = "Forest",
     pch = 21, 
     bg = "grey", 
     cex = 2,
     ylim = c(-0.06, 0.06),
     cex.lab = 1.5)
#mtext("d", 2, adj=5,las=1,padj=-12,cex=2)
par(new = TRUE)

plot(seq(0,1,0.1),
     fixef(model_all)[1]+fixef(model_all)[2]*seq(0,1,0.1),
     type="l",
     lwd=2,
     ylim = c(-0.06, 0.06),
     xaxt = "n",
     yaxt = "n",
     ylab = "",
     xlab = "",
     col = "black")
axis(side = 1,
     at = seq(0, 1, 0.1),
     labels = seq(0, 100, 10))
dev.off()

f_pctch_p(model_all)

# make bar charts -----------------------------------

# survival

temp <- FIA_plot
model_all <- f_surv(temp)
model_forest <- f_pctch_p(model_all)

temp <- FIA_plot %>% filter(STDAGE <= 60) #young
model_all <- f_surv(temp)
model_young <- f_pctch_p(model_all)

temp <- FIA_plot %>% filter(STDAGE > 60) #old
model_all <- f_surv(temp)
model_old <- f_pctch_p(model_all)

temp <- FIA_plot %>% filter(CARBON_LITTER <= 7.78) #low litter
model_all <- f_surv(temp)
model_lowlitter <- f_pctch_p(model_all)

temp <- FIA_plot %>% filter(CARBON_LITTER > 7.78) #high litter
model_all <- f_surv(temp)
model_highlitter <- f_pctch_p(model_all)

temp <- FIA_plot %>% filter(CARBON_SOIL_ORG <= 35.5) #low soil C
model_all <- f_surv(temp)
model_lowsoilc <- f_pctch_p(model_all)

temp <- FIA_plot %>% filter(CARBON_SOIL_ORG > 35.5) #high soil C
model_all <- f_surv(temp)
model_highsoilc <- f_pctch_p(model_all)

temp <- FIA_plot %>% filter(CARBON_UNDER_AG <= 0.96) #low understory
model_all <- f_surv(temp)
model_lowunder <- f_pctch_p(model_all)

temp <- FIA_plot %>% filter(CARBON_UNDER_AG > 0.96) #high understory
model_all <- f_surv(temp)
model_highunder <- f_pctch_p(model_all)

group <- c("forest","forest","age","age","litter", "litter", "soil C", "soil C",
           "understory", "understory")
status <- c("all","all","low","high","low","high","low","high","low","high")
change <- c(model_forest[[1]],
            0,
            model_young[[1]],
            model_old[[1]],
            model_lowlitter[[1]],
            model_highlitter[[1]],
            model_lowsoilc[[1]],
            model_highsoilc[[1]],
            model_lowunder[[1]],
            model_highunder[[1]])
se <- c(model_forest[[2]],
        0,
        model_young[[2]],
        model_old[[2]],
        model_lowlitter[[2]],
        model_highlitter[[2]],
        model_lowsoilc[[2]],
        model_highsoilc[[2]],
        model_lowunder[[2]],
        model_highunder[[2]])
obs <- c(model_forest[[3]],
         0,
         model_young[[3]],
         model_old[[3]],
         model_lowlitter[[3]],
         model_highlitter[[3]],
         model_lowsoilc[[3]],
         model_highsoilc[[3]],
         model_lowunder[[3]],
         model_highunder[[3]])
changes <- as.data.frame(cbind(group,status,change,se,obs))
changes$change <- as.numeric(as.character(changes$change))
changes$se <- as.numeric(as.character(changes$se))
changes$obs <- as.numeric(as.character(changes$obs))
changes$group <- factor(changes$group, 
                        levels = c("forest","age", "litter", "soil C", "understory"))
changes
saveRDS(changes, "output/plot_changes_surv.RDS")


# setwd("/Users/Anika/Documents/GradSchool/FIA_CompetitionModel/")
png("output/changes_plot_mort.png", 
    width = 8,
    height = 4,
    units = 'in',
    res = 300)
ggplot(data = changes, 
       aes(x = factor(group), y = change, fill = status, width = 0.65)) +
  geom_bar(stat = "identity", 
           position = position_dodge()) + 
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Group", y = "% Change", subtitle = "Mortality - Plot level") +
  geom_hline(yintercept = 0, color = "black") +
  geom_errorbar(aes(ymin = change-se, ymax = change+se), width = 0.2,
                position = position_dodge(width = 0.65)) 
#geom_text(aes(label = paste("n = ", obs, sep="")), 
#         position = position_dodge(width = 0.65), vjust = -1)
# theme(axis.text.x = element_text(colour="grey20",size=20,hjust=.5,vjust=.5,face="plain"),
#       axis.text.y = element_text(colour="grey20",size=20,hjust=.5,vjust=.5,face="plain"),
#       base_size=22)
dev.off()


# biomass accumulation

temp <- FIA_plot
model_all <- f_growth(temp)
model_forest <- f_pctch_p(model_all)

temp <- FIA_plot %>% filter(STDAGE <= 60) #young
model_all <- f_growth(temp)
model_young <- f_pctch_p(model_all)

temp <- FIA_plot %>% filter(STDAGE > 60) #old
model_all <- f_growth(temp)
model_old <- f_pctch_p(model_all)

temp <- FIA_plot %>% filter(CARBON_LITTER <= 7.78) #low litter
model_all <- f_growth(temp)
model_lowlitter <- f_pctch_p(model_all)

temp <- FIA_plot %>% filter(CARBON_LITTER > 7.78) #high litter
model_all <- f_growth(temp)
model_highlitter <- f_pctch_p(model_all)

temp <- FIA_plot %>% filter(CARBON_SOIL_ORG <= 35.5) #low soil C
model_all <- f_growth(temp)
model_lowsoilc <- f_pctch_p(model_all)

temp <- FIA_plot %>% filter(CARBON_SOIL_ORG > 35.5) #high soil C
model_all <- f_growth(temp)
model_highsoilc <- f_pctch_p(model_all)

temp <- FIA_plot %>% filter(CARBON_UNDER_AG <= 0.96) #low understory
model_all <- f_growth(temp)
model_lowunder <- f_pctch_p(model_all)

temp <- FIA_plot %>% filter(CARBON_UNDER_AG > 0.96) #high understory
model_all <- f_growth(temp)
model_highunder <- f_pctch_p(model_all)

group <- c("forest","forest","age","age","litter", "litter", "soil C", "soil C",
           "understory", "understory")
status <- c("all","all","low","high","low","high","low","high","low","high")
change <- c(model_forest[[1]],
            0,
            model_young[[1]],
            model_old[[1]],
            model_lowlitter[[1]],
            model_highlitter[[1]],
            model_lowsoilc[[1]],
            model_highsoilc[[1]],
            model_lowunder[[1]],
            model_highunder[[1]])
se <- c(model_forest[[2]],
        0,
        model_young[[2]],
        model_old[[2]],
        model_lowlitter[[2]],
        model_highlitter[[2]],
        model_lowsoilc[[2]],
        model_highsoilc[[2]],
        model_lowunder[[2]],
        model_highunder[[2]])
obs <- c(model_forest[[3]],
         0,
         model_young[[3]],
         model_old[[3]],
         model_lowlitter[[3]],
         model_highlitter[[3]],
         model_lowsoilc[[3]],
         model_highsoilc[[3]],
         model_lowunder[[3]],
         model_highunder[[3]])
changes <- as.data.frame(cbind(group,status,change,se,obs))
changes$change <- as.numeric(as.character(changes$change))
changes$se <- as.numeric(as.character(changes$se))
changes$obs <- as.numeric(as.character(changes$obs))
changes$group <- factor(changes$group, 
                        levels = c("forest", "age", "litter", "soil C", "understory"))
changes
saveRDS(changes, "output/plot_changes_biomassacc.RDS")

# setwd("/Users/Anika/Documents/GradSchool/FIA_CompetitionModel/")
png("output/changes_plot_growth.png", 
    width = 8,
    height = 4,
    units = 'in',
    res = 300)
ggplot(data = changes, 
       aes(x = factor(group), y = change, fill = status, width = 0.65)) +
  geom_bar(stat = "identity", 
           position = position_dodge()) + 
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Group", y = "% Change", subtitle = "Biomass Accumulation - Plot level") +
  geom_hline(yintercept = 0, color = "black") +
  geom_errorbar(aes(ymin = change-se, ymax = change+se), width = 0.2,
                position = position_dodge(width = 0.65)) 
#geom_text(aes(label = paste("n = ", obs, sep="")), 
#         position = position_dodge(width = 0.65), vjust = -1)
# theme(axis.text.x = element_text(colour="grey20",size=20,hjust=.5,vjust=.5,face="plain"),
#       axis.text.y = element_text(colour="grey20",size=20,hjust=.5,vjust=.5,face="plain"),
#       base_size=22)
dev.off()

# recruitment

temp <- FIA_plot
model_all <- f_recr(temp)
model_forest <- f_pctch_p(model_all)

temp <- FIA_plot %>% filter(STDAGE <= 60) #young
model_all <- f_recr(temp)
model_young <- f_pctch_p(model_all)

temp <- FIA_plot %>% filter(STDAGE > 60) #old
model_all <- f_recr(temp)
model_old <- f_pctch_p(model_all)

temp <- FIA_plot %>% filter(CARBON_LITTER <= 7.78) #low litter
model_all <- f_recr(temp)
model_lowlitter <- f_pctch_p(model_all)

temp <- FIA_plot %>% filter(CARBON_LITTER > 7.78) #high litter
model_all <- f_recr(temp)
model_highlitter <- f_pctch_p(model_all)

temp <- FIA_plot %>% filter(CARBON_SOIL_ORG <= 35.5) #low soil C
model_all <- f_recr(temp)
model_lowsoilc <- f_pctch_p(model_all)

temp <- FIA_plot %>% filter(CARBON_SOIL_ORG > 35.5) #high soil C
model_all <- f_recr(temp)
model_highsoilc <- f_pctch_p(model_all)

temp <- FIA_plot %>% filter(CARBON_UNDER_AG <= 0.96) #low understory
model_all <- f_recr(temp)
model_lowunder <- f_pctch_p(model_all)

temp <- FIA_plot %>% filter(CARBON_UNDER_AG > 0.96) #high understory
model_all <- f_recr(temp)
model_highunder <- f_pctch_p(model_all)

group <- c("forest","forest","age","age","litter", "litter", "soil C", "soil C",
           "understory", "understory")
status <- c("all","all","low","high","low","high","low","high","low","high")
change <- c(model_forest[[1]],
            0,
            model_young[[1]],
            model_old[[1]],
            model_lowlitter[[1]],
            model_highlitter[[1]],
            model_lowsoilc[[1]],
            model_highsoilc[[1]],
            model_lowunder[[1]],
            model_highunder[[1]])
se <- c(model_forest[[2]],
        0,
        model_young[[2]],
        model_old[[2]],
        model_lowlitter[[2]],
        model_highlitter[[2]],
        model_lowsoilc[[2]],
        model_highsoilc[[2]],
        model_lowunder[[2]],
        model_highunder[[2]])
obs <- c(model_forest[[3]],
         0,
         model_young[[3]],
         model_old[[3]],
         model_lowlitter[[3]],
         model_highlitter[[3]],
         model_lowsoilc[[3]],
         model_highsoilc[[3]],
         model_lowunder[[3]],
         model_highunder[[3]])
changes <- as.data.frame(cbind(group,status,change,se,obs))
changes$change <- as.numeric(as.character(changes$change))
changes$se <- as.numeric(as.character(changes$se))
changes$obs <- as.numeric(as.character(changes$obs))
changes$group <- factor(changes$group, 
                        levels = c("forest", "age", "litter", "soil C", "understory"))
changes
saveRDS(changes, "output/plot_changes_recr.RDS")


# setwd("/Users/Anika/Documents/GradSchool/FIA_CompetitionModel/")
png("output/changes_plot_recr.png", 
    width = 8,
    height = 4,
    units = 'in',
    res = 300)
ggplot(data = changes, 
       aes(x = factor(group), y = change, fill = status, width = 0.65)) +
  geom_bar(stat = "identity", 
           position = position_dodge()) + 
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Group", y = "% Change", subtitle = "Recruitment - Plot level") +
  geom_hline(yintercept = 0, color = "black") +
  geom_errorbar(aes(ymin = change-se, ymax = change+se), width = 0.2,
                position = position_dodge(width = 0.65)) 
#geom_text(aes(label = paste("n = ", obs, sep="")), 
#         position = position_dodge(width = 0.65), vjust = -1)
# theme(axis.text.x = element_text(colour="grey20",size=20,hjust=.5,vjust=.5,face="plain"),
#       axis.text.y = element_text(colour="grey20",size=20,hjust=.5,vjust=.5,face="plain"),
#       base_size=22)
dev.off()




# mapping plot level effects --------------------------------------------------------
p <- fread("/Users/Anika/Documents/GradSchool/FIA_CompetitionModel/data/raw/PLOT.csv", 
           select = c("CN", "LAT", "LON"))

gs <- merge(FIA_plot, p, by.x="pcn", by.y="CN", all.x=T, all.y=F)

# gs$BA_Nfixer_prop_totals <- scale(gs$BA_Nfixer_prop_total)
# gs$BA_total1s <- scale(gs$BA_total1)
# gs$BA_total2s <- scale(gs$BA_total2)
# gs$BA_Nfixer1s <- scale(gs$BA_Nfixer1)
# gs$BA_nonfixer1s <- scale(gs$BA_nonfixer1)

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
dlat <- 2 # grid cell size
dlon <- 2 # grid cell size
lat.list <- seq(latminint,latmaxint,dlat)
lon.list <- seq(lonminint,lonmaxint,dlon)
nlat <- length(lat.list)
nlon <- length(lon.list)

## Make matrix of plot numbers

#sig.grid: alpha values for significance
#coef.grid: coefficient for NCI_props
#num.grid: number of N-fixers present
#nump.grid: number of plots
#riz.grid: average of ACT1vsRIZ0
#gen.grid: most common genus
#se.grid: standard error of regression for NCI_props
#pct.grid: ratio of coefficient for NCI_props to intercept

#### BA Change maps
coef.grid <- num.grid <- sig.grid <- nump.grid <- riz.grid <- water.grid <- se.grid <- gen.grid <- pct.grid <- array(dim=c(nlon,nlat))

for(i in 1:nlon){
  print(i/nlon*100)
  lon.a <- lon.list[i] - dlon/2
  lon.b <- lon.list[i] + dlon/2
  for(j in 1:nlat){
    lat.a <- lat.list[j] - dlat/2
    lat.b <- lat.list[j] + dlat/2
    
    #plots <- p[(p$lon1>=lon.a & p$lon1<lon.b & p$lat1>=lat.a & p$lat1<lat.b),]$PCN
    temp <- gs[(gs$LON>=lon.a & gs$LON<lon.b & gs$LAT>=lat.a & gs$LAT<lat.b),]
    #if(length(plots)>0){
     # temp <- gs[which(gs$pcn %in% plots |gs$prev_plt_cn %in% plots), ]
      if(nrow(temp)>0){
        if(length(unique(temp$pcn))>5 && nrow(temp[temp$BA_change>0,])>5 && length(unique(temp$BA_Nfixer_prop_total))>2){
          temp$BA_Nfixer_prop_totals <- scale(temp$BA_Nfixer_prop_total)
          temp$BA_total1s <- scale(temp$BA_total1)
          temp$BA_total2s <- scale(temp$BA_total2)
          temp$BA_Nfixer1s <- scale(temp$BA_Nfixer1)
          temp$BA_nonfixer1s <- scale(temp$BA_nonfixer1)
          model.fixed2 = lmer(BA_change ~ BA_Nfixer_prop_totals + 
                                BA_total1s + 
                                BA_total2s + 
                                BA_Nfixer1s + 
                                BA_nonfixer1s +
                                (1|mt1), 
                              data = temp) #removed state b/c usually in same state
          coefs <- data.frame(coef(summary(model.fixed2)))
          coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
          sig.grid[i,j] <- coefs["BA_Nfixer_prop_totals","p.z"]
          coef.grid[i,j] <- coefs["BA_Nfixer_prop_totals","Estimate"]
          nump.grid[i,j] <- length(unique(temp$pcn))
          se.grid[i,j] <- sqrt(diag(vcov(model.fixed2)))[[2]]
          pct.grid[i,j] <- 100*((coefs[1,"Estimate"]-(coefs[1,"Estimate"]+coefs["BA_Nfixer_prop_totals","Estimate"]*1))/coefs[1,"Estimate"])

        }else{
          sig.grid[i,j] <- NA
          coef.grid[i,j] <- NA
          se.grid[i,j] <- NA
          nump.grid[i,j] <- NA
          pct.grid[i,j] <- NA
        }
        
      }else{
        sig.grid[i,j] <- NA
        coef.grid[i,j] <- NA
        se.grid[i,j] <- NA
        nump.grid[i,j] <- NA
        pct.grid[i,j] <- NA
      }
    # }else{
    #   sig.grid[i,j] <- NA
    #   coef.grid[i,j] <- NA
    #   se.grid[i,j] <- NA
    #   nump.grid[i,j] <- NA
    #   pct.grid[i,j] <- NA
    #   
    # }
  }
}

# coefficient grid 
require(RColorBrewer)
nHalf = sum(!is.na(coef.grid))/2
Min = min(coef.grid, na.rm=T)
Max = max(coef.grid, na.rm=T)
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
ticks <- c(min(coef.grid, na.rm=T), -2, -1, 0, 1, 2, max(coef.grid, na.rm=T))
image.plot(lon.list, lat.list, coef.grid, 
           nlevel = nlev, 
           col = rampcols, 
           breaks = rampbreaks,
           axis.args = list(at=ticks,labels=ticks),
           xlab = "longitude", ylab = "latitude",
           #ylim=c(14,33),xlim=c(-120,-80))
           ylim = c(10, 60), xlim = c(-160, -50))
#ylim=c(latminint-1,latmaxint+1),xlim=c(lonminint-2,lonmaxint+1))
title(main = "BA Change Regression Slope -- Plot scale")
#mtext("(a)",font=2,adj=0.01,padj=-0.5)
map('world', regions = 'usa', add = TRUE)



##### Recruitment maps
coefr.grid <- num.grid <- sig.grid <- nump.grid <- riz.grid <- water.grid <- se.grid <- gen.grid <- pct.grid <- array(dim=c(nlon,nlat))

for(i in 1:nlon){
  print(i/nlon*100)
  lon.a <- lon.list[i] - dlon/2
  lon.b <- lon.list[i] + dlon/2
  for(j in 1:nlat){
    lat.a <- lat.list[j] - dlat/2
    lat.b <- lat.list[j] + dlat/2
    
    #plots <- p[(p$lon1>=lon.a & p$lon1<lon.b & p$lat1>=lat.a & p$lat1<lat.b),]$PCN
    temp <- gs[(gs$LON>=lon.a & gs$LON<lon.b & gs$LAT>=lat.a & gs$LAT<lat.b),]
    #if(length(plots)>0){
    # temp <- gs[which(gs$pcn %in% plots |gs$prev_plt_cn %in% plots), ]
    if(nrow(temp)>0){
      if(length(unique(temp$pcn))>5 && nrow(temp[temp$recyr>0,])>5 && length(unique(temp$BA_Nfixer_prop_total))>2){
        temp$BA_Nfixer_prop_totals <- scale(temp$BA_Nfixer_prop_total)
        temp$BA_total1s <- scale(temp$BA_total1)
        temp$BA_total2s <- scale(temp$BA_total2)
        temp$BA_Nfixer1s <- scale(temp$BA_Nfixer1)
        temp$BA_nonfixer1s <- scale(temp$BA_nonfixer1)
        model.fixed2 = lmer(recyr ~ BA_Nfixer_prop_totals + 
                              BA_total1s + 
                              BA_total2s + 
                              BA_Nfixer1s + 
                              BA_nonfixer1s +
                              (1|mt1), 
                            data = temp) #removed state b/c usually in same state
        coefs <- data.frame(coef(summary(model.fixed2)))
        coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
        sig.grid[i,j] <- coefs["BA_Nfixer_prop_totals","p.z"]
        coefr.grid[i,j] <- coefs["BA_Nfixer_prop_totals","Estimate"]
        nump.grid[i,j] <- length(unique(temp$pcn))
        se.grid[i,j] <- sqrt(diag(vcov(model.fixed2)))[[2]]
        pct.grid[i,j] <- 100*((coefs[1,"Estimate"]-(coefs[1,"Estimate"]+coefs["BA_Nfixer_prop_totals","Estimate"]*1))/coefs[1,"Estimate"])
        
      }else{
        sig.grid[i,j] <- NA
        coefr.grid[i,j] <- NA
        se.grid[i,j] <- NA
        nump.grid[i,j] <- NA
        pct.grid[i,j] <- NA
      }
      
    }else{
      sig.grid[i,j] <- NA
      coefr.grid[i,j] <- NA
      se.grid[i,j] <- NA
      nump.grid[i,j] <- NA
      pct.grid[i,j] <- NA
    }
    # }else{
    #   sig.grid[i,j] <- NA
    #   coef.grid[i,j] <- NA
    #   se.grid[i,j] <- NA
    #   nump.grid[i,j] <- NA
    #   pct.grid[i,j] <- NA
    #   
    # }
  }
}

# coefficient grid 
require(RColorBrewer)
nHalf = sum(!is.na(coefr.grid))/2
Min = min(coefr.grid, na.rm=T)
Max = max(coefr.grid, na.rm=T)
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
ticks <- c(min(coefr.grid, na.rm=T), -2, -1, 0, 1, 2, max(coefr.grid, na.rm=T))
image.plot(lon.list, lat.list, coefr.grid, 
           nlevel = nlev, 
           col = rampcols, 
           breaks = rampbreaks,
           axis.args = list(at=ticks,labels=ticks),
           xlab = "longitude", ylab = "latitude",
           #ylim=c(14,33),xlim=c(-120,-80))
           ylim = c(10, 60), xlim = c(-160, -50))
#ylim=c(latminint-1,latmaxint+1),xlim=c(lonminint-2,lonmaxint+1))
title(main = "Recruitment Regression Slope -- Plot scale")
#mtext("(a)",font=2,adj=0.01,padj=-0.5)
map('world', regions = 'usa', add = TRUE)


##### Mortality maps
coefm.grid <- num.grid <- sig.grid <- nump.grid <- riz.grid <- water.grid <- se.grid <- gen.grid <- pct.grid <- array(dim=c(nlon,nlat))

for(i in 1:nlon){
  print(i/nlon*100)
  lon.a <- lon.list[i] - dlon/2
  lon.b <- lon.list[i] + dlon/2
  for(j in 1:nlat){
    lat.a <- lat.list[j] - dlat/2
    lat.b <- lat.list[j] + dlat/2
    
    temp <- gs[(gs$LON>=lon.a & gs$LON<lon.b & gs$LAT>=lat.a & gs$LAT<lat.b),]
    if(nrow(temp) > 0){
      if(length(unique(temp$pcn)) > 5 && 
         nrow(temp[temp$mortyr > 0, ]) > 5 && 
         length(unique(temp$BA_Nfixer_prop_total)) > 2 &&
         length(unique(temp$mt1)) > 2){
        temp$BA_Nfixer_prop_totals <- scale(temp$BA_Nfixer_prop_total)
        temp$BA_total1s <- scale(temp$BA_total1)
        temp$BA_total2s <- scale(temp$BA_total2)
        temp$BA_Nfixer2s <- scale(temp$BA_Nfixer2)
        temp$BA_nonfixer2s <- scale(temp$BA_nonfixer2)
        model.fixed2 = lmer(mortyr ~ BA_Nfixer_prop_totals + 
                              BA_total1s + 
                              BA_total2s + 
                             # BA_Nfixer2s + 
                              BA_nonfixer2s +
                              (1|mt1), 
                            data = temp) #removed state b/c usually in same state
        coefs <- data.frame(coef(summary(model.fixed2)))
        coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
        sig.grid[i,j] <- coefs["BA_Nfixer_prop_totals","p.z"]
        coefm.grid[i,j] <- coefs["BA_Nfixer_prop_totals","Estimate"]
        nump.grid[i,j] <- length(unique(temp$pcn))
        se.grid[i,j] <- sqrt(diag(vcov(model.fixed2)))[[2]]
        pct.grid[i,j] <- 100*((coefs[1,"Estimate"]-(coefs[1,"Estimate"]+coefs["BA_Nfixer_prop_totals","Estimate"]*1))/coefs[1,"Estimate"])
        
      }else{
        sig.grid[i,j] <- NA
        coefm.grid[i,j] <- NA
        se.grid[i,j] <- NA
        nump.grid[i,j] <- NA
        pct.grid[i,j] <- NA
      }
      
    }else{
      sig.grid[i,j] <- NA
      coefm.grid[i,j] <- NA
      se.grid[i,j] <- NA
      nump.grid[i,j] <- NA
      pct.grid[i,j] <- NA
    }
    # }else{
    #   sig.grid[i,j] <- NA
    #   coef.grid[i,j] <- NA
    #   se.grid[i,j] <- NA
    #   nump.grid[i,j] <- NA
    #   pct.grid[i,j] <- NA
    #   
    # }
  }
}

# coefficient grid 
require(RColorBrewer)
nHalf = sum(!is.na(coefm.grid))/2
Min = min(coefm.grid, na.rm=T)
Max = max(coefm.grid, na.rm=T)
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
ticks <- c(min(coefm.grid, na.rm=T), -2, -1, 0, 1, 2, max(coefm.grid, na.rm=T))
image.plot(lon.list, lat.list, coefm.grid, 
           nlevel = nlev, 
           col = rampcols, 
           breaks = rampbreaks,
           axis.args = list(at=ticks,labels=ticks),
           xlab = "longitude", ylab = "latitude",
           #ylim=c(14,33),xlim=c(-120,-80))
           ylim = c(10, 60), xlim = c(-160, -50))
#ylim=c(latminint-1,latmaxint+1),xlim=c(lonminint-2,lonmaxint+1))
title(main = "Mortality Regression Slope -- Plot scale")
#mtext("(a)",font=2,adj=0.01,padj=-0.5)
map('world', regions = 'usa', add = TRUE)



# map plot characteristics -------------------------------------
plot <- fread("data/raw/plot.csv", select = c('CN', 
                                              'LAT',
                                              'LON'))
FIA_plot <- merge(FIA_plot, plot, by.x = "pcn", by.y="CN", all.x=T, all.y=F)

# understory C
par(bty = "n")
usa <- map_data("usa")
world <- map_data("world", xlim = c(-140, -65), ylim = c(23, 63))

ggplot() + 
  geom_polygon(data = usa, aes(x = long, y = lat, group = group), 
               alpha = 0.01, 
               colour = 'black', 
               fill = NA) + 
  coord_fixed(1.3) +
  geom_point(aes(x = FIA_plot$LON, y = FIA_plot$LAT, 
                 color = FIA_plot$CARBON_UNDER_AG)) +
  geom_jitter() +
  labs(x = " ", y = " ", title = "Understory C") +
  labs(color = 'Understory C') +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20)) +
  scale_colour_gradient(low = "#e5f9d9", high = "dark green")  

cor.test(FIA_plot$CARBON_UNDER_AG, FIA_plot$STDAGE)


# understory C
require(colorRamps)
par(bty = "n")
usa <- map_data("usa")
world <- map_data("world", xlim = c(-140, -65), ylim = c(23, 63))

ggplot() + 
  geom_polygon(data = usa, aes(x = long, y = lat, group = group), 
               alpha = 0.01, 
               colour = 'black', 
               fill = NA) + 
  coord_fixed(1.3) +
  geom_point(aes(x = FIA_plot$LON, y = FIA_plot$LAT, 
                 color = FIA_plot$CARBON_SOIL_ORG)) +
  geom_jitter() +
  labs(x = " ", y = " ", title = "Organic Soil C") +
  labs(color = 'Organic Soil C') +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20)) +
  scale_colour_gradientn(colours = matlab.like(10)) 

# litter C
par(bty = "n")
usa <- map_data("usa")
world <- map_data("world", xlim = c(-140, -65), ylim = c(23, 63))

ggplot() + 
  geom_polygon(data = usa, aes(x = long, y = lat, group = group), 
               alpha = 0.01, 
               colour = 'black', 
               fill = NA) + 
  coord_fixed(1.3) +
  geom_point(aes(x = FIA_plot$LON, y = FIA_plot$LAT, 
                 color = FIA_plot$CARBON_LITTER)) +
  geom_jitter() +
  labs(x = " ", y = " ", title = "Litter C") +
  labs(color = 'Litter C') +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20)) +
  scale_colour_gradientn(colours = matlab.like(10)) 



SM1_m <- f_pctch_p(f_surv(FIA_plot[FIA_plot$sm >= 0 & FIA_plot$sm < 100,]))
SM2_m <- f_pctch_p(f_surv(FIA_plot[FIA_plot$sm >= 100 & FIA_plot$sm < 200,]))
SM3_m <- f_pctch_p(f_surv(FIA_plot[FIA_plot$sm >= 200 & FIA_plot$sm < 300,]))
SM4_m <- f_pctch_p(f_surv(FIA_plot[FIA_plot$sm >= 300 & FIA_plot$sm < 400,]))
SM5_m <- f_pctch_p(f_surv(FIA_plot[FIA_plot$sm >= 400 & FIA_plot$sm < 500,]))
SM6_m <- f_pctch_p(f_surv(FIA_plot[FIA_plot$sm >= 500 & FIA_plot$sm < 600,]))
SM7_m <- f_pctch_p(f_surv(FIA_plot[FIA_plot$sm >= 600,]))

sm <- rbind.data.frame(SM1_m, SM2_m, SM3_m, SM4_m, SM5_m, SM6_m, SM7_m)
colnames(sm) <- c("est", "se", "num")
quad <- c("0-100", "100-200", "200-300", "300-400", "400-500", "500-600", "600+")
sm <- cbind(sm, quad)
sm$type <- 'survival'

SM1_r <- f_pctch_p(f_recr(FIA_plot[FIA_plot$sm >= 0 & FIA_plot$sm < 100,]))
SM2_r <- f_pctch_p(f_recr(FIA_plot[FIA_plot$sm >= 100 & FIA_plot$sm < 200,]))
SM3_r <- f_pctch_p(f_recr(FIA_plot[FIA_plot$sm >= 200 & FIA_plot$sm < 300,]))
SM4_r <- f_pctch_p(f_recr(FIA_plot[FIA_plot$sm >= 300 & FIA_plot$sm < 400,]))
SM5_r <- f_pctch_p(f_recr(FIA_plot[FIA_plot$sm >= 400 & FIA_plot$sm < 500,]))
SM6_r <- f_pctch_p(f_recr(FIA_plot[FIA_plot$sm >= 500 & FIA_plot$sm < 600,]))
SM7_r <- f_pctch_p(f_recr(FIA_plot[FIA_plot$sm >= 600,]))

sr <- rbind.data.frame(SM1_r, SM2_r, SM3_r, SM4_r, SM5_r, SM6_r, SM7_r)
colnames(sr) <- c("est", "se", "num")
quad <- c("0-100", "100-200", "200-300", "300-400", "400-500", "500-600", "600+")
sr <- cbind(sr, quad)
sr$type <- 'recruitment'

SM1_g <- f_pctch_p(f_growth(FIA_plot[FIA_plot$sm >= 0 & FIA_plot$sm < 100,]))
SM2_g <- f_pctch_p(f_growth(FIA_plot[FIA_plot$sm >= 100 & FIA_plot$sm < 200,]))
SM3_g <- f_pctch_p(f_growth(FIA_plot[FIA_plot$sm >= 200 & FIA_plot$sm < 300,]))
SM4_g <- f_pctch_p(f_growth(FIA_plot[FIA_plot$sm >= 300 & FIA_plot$sm < 400,]))
SM5_g <- f_pctch_p(f_growth(FIA_plot[FIA_plot$sm >= 400 & FIA_plot$sm < 500,]))
SM6_g <- f_pctch_p(f_growth(FIA_plot[FIA_plot$sm >= 500 & FIA_plot$sm < 600,]))
SM7_g <- f_pctch_p(f_growth(FIA_plot[FIA_plot$sm >= 600,]))

sg <- rbind.data.frame(SM1_g, SM2_g, SM3_g, SM4_g, SM5_g, SM6_g, SM7_g)
colnames(sg) <- c("est", "se", "num")
quad <- c("0-100", "100-200", "200-300", "300-400", "400-500", "500-600", "600+")
sg <- cbind(sg, quad)
sg$type <- 'growth'

sm <- rbind(sm, sr, sg)
sm[8:21,"num"] <- NA

saveRDS(sm, "output/soilmoisture.RDS")

ggplot(data = sm, 
       aes(x = factor(quad), y = est, fill = type, width = 0.65)) +
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Soil moisture (mm/yr)", y = "% Change", 
       title = "% Change in Demographic Rate by Soil Moisture Group") +
  geom_hline(yintercept = 0, color = "black") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_errorbar(aes(ymin = est-se, ymax = est+se), width = 0.2,
                position = position_dodge(width = 0.65)) +
  geom_text(aes(label = ifelse(is.na(num), "", paste("n = ", num, sep=""))), 
            position = position_dodge(width = 0.65), vjust = -1)


D1_m <- f_pctch_p(f_mort(dstrb[dstrb$dtype %in% c('none'),]))
D2_m <- f_pctch_p(f_mort(dstrb[dstrb$dtype %in% c('insect'),]))
D3_m <- f_pctch_p(f_mort(dstrb[dstrb$dtype %in% c('disease'),]))
D4_m <- f_pctch_p(f_mort(dstrb[dstrb$dtype %in% c('fire'),]))
D5_m <- f_pctch_p(f_mort(dstrb[dstrb$dtype %in% c('animal'),]))
D6_m <- f_pctch_p(f_mort(dstrb[dstrb$dtype %in% c('weather'),]))
D7_m <- f_pctch_p(f_mort(dstrb[dstrb$dtype %in% c('vegetation'),]))
D8_m <- f_pctch_p(f_mort(dstrb[dstrb$dtype %in% c('unknown'),]))
D9_m <- f_pctch_p(f_mort(dstrb[dstrb$dtype %in% c('human'),]))
D10_m <- f_pctch_p(f_mort(dstrb[dstrb$dtype %in% c('geology'),]))

dm <- rbind.data.frame(D1_m, D2_m, D3_m, D4_m, D5_m, D6_m, D7_m, D8_m, D9_m, D10_m)
colnames(dm) <- c("est", "se", "num")
quad <- c('none','insect','disease','fire','animal','weather','vegetation','unknown',
          'human','geology')
dm <- cbind(dm, quad)
dm$resp <- 'mort'

D1_r <- f_pctch_p(f_recr(dstrb[dstrb$dtype %in% c('none'),]))
D2_r <- f_pctch_p(f_recr(dstrb[dstrb$dtype %in% c('insect'),]))
D3_r <- f_pctch_p(f_recr(dstrb[dstrb$dtype %in% c('disease'),]))
D4_r <- f_pctch_p(f_recr(dstrb[dstrb$dtype %in% c('fire'),]))
D5_r <- f_pctch_p(f_recr(dstrb[dstrb$dtype %in% c('animal'),]))
D6_r <- f_pctch_p(f_recr(dstrb[dstrb$dtype %in% c('weather'),]))
D7_r <- f_pctch_p(f_recr(dstrb[dstrb$dtype %in% c('vegetation'),]))
D8_r <- f_pctch_p(f_recr(dstrb[dstrb$dtype %in% c('unknown'),]))
D9_r <- f_pctch_p(f_recr(dstrb[dstrb$dtype %in% c('human'),]))
D10_r <- f_pctch_p(f_recr(dstrb[dstrb$dtype %in% c('geology'),]))

dr <- rbind.data.frame(D1_r, D2_r, D3_r, D4_r, D5_r, D6_r, D7_r, D8_r, D9_r, D10_r)
colnames(dr) <- c("est", "se", "num")
quad <- c('none','insect','disease','fire','animal','weather','vegetation','unknown',
          'human','geology')
dr <- cbind(dr, quad)
dr$resp <- 'recr'

D1_g <- f_pctch_p(f_growth(dstrb[dstrb$dtype %in% c('none'),]))
D2_g <- f_pctch_p(f_growth(dstrb[dstrb$dtype %in% c('insect'),]))
D3_g <- f_pctch_p(f_growth(dstrb[dstrb$dtype %in% c('disease'),]))
D4_g <- f_pctch_p(f_growth(dstrb[dstrb$dtype %in% c('fire'),]))
D5_g <- f_pctch_p(f_growth(dstrb[dstrb$dtype %in% c('animal'),]))
D6_g <- f_pctch_p(f_growth(dstrb[dstrb$dtype %in% c('weather'),]))
D7_g <- f_pctch_p(f_growth(dstrb[dstrb$dtype %in% c('vegetation'),]))
D8_g <- f_pctch_p(f_growth(dstrb[dstrb$dtype %in% c('unknown'),]))
D9_g <- f_pctch_p(f_growth(dstrb[dstrb$dtype %in% c('human'),]))
D10_g <- f_pctch_p(f_growth(dstrb[dstrb$dtype %in% c('geology'),]))

dg <- rbind.data.frame(D1_g, D2_g, D3_g, D4_g, D5_g, D6_g, D7_g, D8_g, D9_g, D10_g)
colnames(dg) <- c("est", "se", "num")
quad <- c('none','insect','disease','fire','animal','weather','vegetation','unknown',
          'human','geology')
dg <- cbind(dg, quad)
dg$resp <- 'growth'

dm <- rbind(dm, dr, dg)

dm[11:30,"num"] <- NA

ggplot(data = dm, 
       aes(x = factor(quad), y = est, fill = resp, width = 0.65)) +
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Group", y = "% Change", title = "% Change in Demographic Rate by Disturbance Type") +
  geom_hline(yintercept = 0, color = "black") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_errorbar(aes(ymin = est-se, ymax = est+se), width = 0.2,
                position = position_dodge(width = 0.65)) +
  geom_text(aes(label = ifelse(is.na(num), "", paste("n = ", num, sep=""))), 
            position = position_dodge(width = 0.65), vjust = -1)


### run function -------------------------------------
FIA_plot <- readRDS("output/FIA_plot_withgroups.RDS")



cp_s_age <- f_pctch_pg(f_grmodp(FIA_plot, "surv", FIA_plot$YOU0OLD1))
cp_s_ndep <- f_pctch_pg(f_grmodp(FIA_plot, "surv", FIA_plot$NDEP))
cp_s_sm <- f_pctch_pg(f_grmodp(FIA_plot, "surv", FIA_plot$SM))
cp_s_texture <- f_pctch_pg(f_grmodp(FIA_plot, "surv", FIA_plot$texture))
cp_s_fortypsm <- f_pctch_pg(f_grmodp(FIA_plot, "surv", FIA_plot$fortypsm))
cp_s_disturb <- f_pctch_pg(f_grmodp(FIA_plot, "surv", FIA_plot$dtype))
cp_s_ndep2 <- f_pctch_pg(f_grmodp(FIA_plot, "surv", FIA_plot$NDEP2))
cp_s_uc <- f_pctch_pg(f_grmodp(FIA_plot, "surv", FIA_plot$Cunder))

cp_g_age <- f_pctch_pg(f_grmodp(FIA_plot, "growth", FIA_plot$YOU0OLD1))
cp_g_ndep <- f_pctch_pg(f_grmodp(FIA_plot, "growth", FIA_plot$NDEP))
cp_g_sm <- f_pctch_pg(f_grmodp(FIA_plot, "growth", FIA_plot$SM))
cp_g_texture <- f_pctch_pg(f_grmodp(FIA_plot, "growth", FIA_plot$texture))
cp_g_fortypsm <- f_pctch_pg(f_grmodp(FIA_plot, "growth", FIA_plot$fortypsm))
cp_g_disturb <- f_pctch_pg(f_grmodp(FIA_plot, "growth", FIA_plot$dtype))
cp_g_ndep2 <- f_pctch_pg(f_grmodp(FIA_plot, "growth", FIA_plot$NDEP2))
cp_g_uc <- f_pctch_pg(f_grmodp(FIA_plot, "growth", FIA_plot$Cunder))

cp_r_age <- f_pctch_pg(f_grmodp(FIA_plot, "recr", FIA_plot$YOU0OLD1))
cp_r_ndep <- f_pctch_pg(f_grmodp(FIA_plot, "recr", FIA_plot$NDEP))
cp_r_sm <- f_pctch_pg(f_grmodp(FIA_plot, "recr", FIA_plot$SM))
cp_r_texture <- f_pctch_pg(f_grmodp(FIA_plot, "recr", FIA_plot$texture))
cp_r_fortypsm <- f_pctch_pg(f_grmodp(FIA_plot, "recr", FIA_plot$fortypsm))
cp_r_disturb <- f_pctch_pg(f_grmodp(FIA_plot, "recr", FIA_plot$dtype))
cp_r_ndep2 <- f_pctch_pg(f_grmodp(FIA_plot, "recr", FIA_plot$NDEP2))
cp_r_uc <- f_pctch_pg(f_grmodp(FIA_plot, "recr", FIA_plot$Cunder))

cp_s_fortypsm$est <- cp_s_fortypsm$est %>%
  mutate(grp = case_when(grps == 100 ~ "white/red/jack pine",
                         grps == 120 ~ "spruce/fir",
                         grps == 140 ~ "longleaf/slash pine",
                         grps == 150 ~ "tropical softwood",
                         grps == 160 ~ "loblolly/shortleaf pine",
                         grps == 170 ~ "other eastern softwoods",
                         grps == 180 ~ "pinyon/juniper",
                         grps == 200 ~ "douglas-fir",
                         grps == 220 ~ "ponderosa pine",
                         grps == 240 ~ "western white pine",
                         grps == 260 ~ "fir/spruce/mountain hemlock",
                         grps == 280 ~ "lodgepole pine",
                         grps == 300 ~ "hemlock/sitka spruce",
                         grps == 320 ~ "western larch",
                         grps == 340 ~ "redwood",
                         grps == 360 ~ "other western softwoods",
                         grps == 370 ~ "CA mixed conifer",
                         grps == 380 ~ "exotic softwood",
                         grps == 390 ~ "other softwood",
                         grps == 400 ~ "oak/pine",
                         grps == 500 ~ "oak/hickory",
                         grps == 600 ~ "oak/gum/cypress",
                         grps == 700 ~ "elm/ash/cottonwood",
                         grps == 800 ~ "maple/beech/birch",
                         grps == 900 ~ "aspen/birch",
                         grps == 910 ~ "alder/maple",
                         grps == 920 ~ "western oak",
                         grps == 940 ~ "tanoak/laurel",
                         grps == 960 ~ "other hardwoods",
                         grps == 970 ~ "woodland hardwoods",
                         grps == 980 ~ "tropical hardwoods",
                         grps == 988 ~ "cloud forest",
                         grps == 990 ~ "exotic hardwood"))
colnames(cp_s_fortypsm$est) <- c("ch", "grp", "grps")

cp_s_texture$est <- cp_s_texture$est %>%
  mutate(grp = case_when(grps == 1 ~ "clay",
                         grps == 2 ~ "silty clay",
                         grps == 3 ~ "sandy clay",
                         grps == 4 ~ "clay loam",
                         grps == 5 ~ "silty clay loam",
                         grps == 6 ~ "sandy clay loam",
                         grps == 7 ~ "loam",
                         grps == 8 ~ "silty loam",
                         grps == 9 ~ "sandy loam",
                         grps == 10 ~ "silt",
                         grps == 11 ~ "loamy sand",
                         grps == 12 ~ "sand"))
colnames(cp_s_texture$est) <- c("ch", "grp", "grps")

cp_s_age$est <- cp_s_age$est %>%
  mutate(grp = case_when(grps == 0 ~ "young",
                         grps == 1 ~ "old"))
colnames(cp_s_age$est) <- c("ch", "grp", "grps")


cp_g_fortypsm$est <- cp_g_fortypsm$est %>%
  mutate(grp = case_when(grps == 100 ~ "white/red/jack pine",
                         grps == 120 ~ "spruce/fir",
                         grps == 140 ~ "longleaf/slash pine",
                         grps == 150 ~ "tropical softwood",
                         grps == 160 ~ "loblolly/shortleaf pine",
                         grps == 170 ~ "other eastern softwoods",
                         grps == 180 ~ "pinyon/juniper",
                         grps == 200 ~ "douglas-fir",
                         grps == 220 ~ "ponderosa pine",
                         grps == 240 ~ "western white pine",
                         grps == 260 ~ "fir/spruce/mountain hemlock",
                         grps == 280 ~ "lodgepole pine",
                         grps == 300 ~ "hemlock/sitka spruce",
                         grps == 320 ~ "western larch",
                         grps == 340 ~ "redwood",
                         grps == 360 ~ "other western softwoods",
                         grps == 370 ~ "CA mixed conifer",
                         grps == 380 ~ "exotic softwood",
                         grps == 390 ~ "other softwood",
                         grps == 400 ~ "oak/pine",
                         grps == 500 ~ "oak/hickory",
                         grps == 600 ~ "oak/gum/cypress",
                         grps == 700 ~ "elm/ash/cottonwood",
                         grps == 800 ~ "maple/beech/birch",
                         grps == 900 ~ "aspen/birch",
                         grps == 910 ~ "alder/maple",
                         grps == 920 ~ "western oak",
                         grps == 940 ~ "tanoak/laurel",
                         grps == 960 ~ "other hardwoods",
                         grps == 970 ~ "woodland hardwoods",
                         grps == 980 ~ "tropical hardwoods",
                         grps == 988 ~ "cloud forest",
                         grps == 990 ~ "exotic hardwood"))
colnames(cp_g_fortypsm$est) <- c("ch", "grp", "grps")

cp_g_texture$est <- cp_g_texture$est %>%
  mutate(grp = case_when(grps == 1 ~ "clay",
                         grps == 2 ~ "silty clay",
                         grps == 3 ~ "sandy clay",
                         grps == 4 ~ "clay loam",
                         grps == 5 ~ "silty clay loam",
                         grps == 6 ~ "sandy clay loam",
                         grps == 7 ~ "loam",
                         grps == 8 ~ "silty loam",
                         grps == 9 ~ "sandy loam",
                         grps == 10 ~ "silt",
                         grps == 11 ~ "loamy sand",
                         grps == 12 ~ "sand"))
colnames(cp_g_texture$est) <- c("ch", "grp", "grps")

cp_g_age$est <- cp_g_age$est %>%
  mutate(grp = case_when(grps == 0 ~ "young",
                         grps == 1 ~ "old"))
colnames(cp_g_age$est) <- c("ch", "grp", "grps")



cp_r_fortypsm$est <- cp_r_fortypsm$est %>%
  mutate(grp = case_when(grps == 100 ~ "white/red/jack pine",
                         grps == 120 ~ "spruce/fir",
                         grps == 140 ~ "longleaf/slash pine",
                         grps == 150 ~ "tropical softwood",
                         grps == 160 ~ "loblolly/shortleaf pine",
                         grps == 170 ~ "other eastern softwoods",
                         grps == 180 ~ "pinyon/juniper",
                         grps == 200 ~ "douglas-fir",
                         grps == 220 ~ "ponderosa pine",
                         grps == 240 ~ "western white pine",
                         grps == 260 ~ "fir/spruce/mountain hemlock",
                         grps == 280 ~ "lodgepole pine",
                         grps == 300 ~ "hemlock/sitka spruce",
                         grps == 320 ~ "western larch",
                         grps == 340 ~ "redwood",
                         grps == 360 ~ "other western softwoods",
                         grps == 370 ~ "CA mixed conifer",
                         grps == 380 ~ "exotic softwood",
                         grps == 390 ~ "other softwood",
                         grps == 400 ~ "oak/pine",
                         grps == 500 ~ "oak/hickory",
                         grps == 600 ~ "oak/gum/cypress",
                         grps == 700 ~ "elm/ash/cottonwood",
                         grps == 800 ~ "maple/beech/birch",
                         grps == 900 ~ "aspen/birch",
                         grps == 910 ~ "alder/maple",
                         grps == 920 ~ "western oak",
                         grps == 940 ~ "tanoak/laurel",
                         grps == 960 ~ "other hardwoods",
                         grps == 970 ~ "woodland hardwoods",
                         grps == 980 ~ "tropical hardwoods",
                         grps == 988 ~ "cloud forest",
                         grps == 990 ~ "exotic hardwood"))
colnames(cp_r_fortypsm$est) <- c("ch", "grp", "grps")

cp_r_texture$est <- cp_r_texture$est %>%
  mutate(grp = case_when(grps == 1 ~ "clay",
                         grps == 2 ~ "silty clay",
                         grps == 3 ~ "sandy clay",
                         grps == 4 ~ "clay loam",
                         grps == 5 ~ "silty clay loam",
                         grps == 6 ~ "sandy clay loam",
                         grps == 7 ~ "loam",
                         grps == 8 ~ "silty loam",
                         grps == 9 ~ "sandy loam",
                         grps == 10 ~ "silt",
                         grps == 11 ~ "loamy sand",
                         grps == 12 ~ "sand"))
colnames(cp_r_texture$est) <- c("ch", "grp", "grps")

cp_r_age$est <- cp_r_age$est %>%
  mutate(grp = case_when(grps == 0 ~ "young",
                         grps == 1 ~ "old"))
colnames(cp_r_age$est) <- c("ch", "grp", "grps")

saveRDS(cp_s_age, "output/cp_s_age.RDS")
saveRDS(cp_s_ndep, "output/cp_s_ndep.RDS")
saveRDS(cp_s_sm, "output/cp_s_sm.RDS")
saveRDS(cp_s_texture, "output/cp_s_texture.RDS")
saveRDS(cp_s_fortypsm, "output/cp_s_fortypsm.RDS")
saveRDS(cp_s_disturb, "output/cp_s_disturb.RDS")

saveRDS(cp_g_age, "output/cp_g_age.RDS")
saveRDS(cp_g_ndep, "output/cp_g_ndep.RDS")
saveRDS(cp_g_sm, "output/cp_g_sm.RDS")
saveRDS(cp_g_texture, "output/cp_g_texture.RDS")
saveRDS(cp_g_fortypsm, "output/cp_g_fortypsm.RDS")
saveRDS(cp_g_disturb, "output/cp_g_disturb.RDS")

saveRDS(cp_r_age, "output/cp_r_age.RDS")
saveRDS(cp_r_ndep, "output/cp_r_ndep.RDS")
saveRDS(cp_r_sm, "output/cp_r_sm.RDS")
saveRDS(cp_r_texture, "output/cp_r_texture.RDS")
saveRDS(cp_r_fortypsm, "output/cp_r_fortypsm.RDS")
saveRDS(cp_r_disturb, "output/cp_r_disturb.RDS")

### counts
xtabs(~YOU0OLD1, FIA_plot)

### plot m
#f_plotch(% change data, demographic type, title, x label size)
require(gridExtra)
p1 <- f_plotch(cp_s_age, "surv", "Age", 20)
p2 <- f_plotch(cp_s_ndep, "surv", "N deposition", 20)
p3 <- f_plotch(cp_s_sm, "surv", "Soil Moisture", 20)
p4 <- f_plotch(cp_s_texture, "surv", "Texture", 20)
p5 <- f_plotch(cp_s_fortypsm, "surv", "Forest Type", 10)
p6 <- f_plotch(cp_s_disturb, "surv", "Disturbance", 20)
grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2)

p1 <- f_plotch(cp_g_age, "growth", "age", 20)
p2 <- f_plotch(cp_g_ndep, "growth", "N deposition", 20)
p3 <- f_plotch(cp_g_sm, "growth", "Soil Moisture", 20)
p4 <- f_plotch(cp_g_texture, "growth", "Texture", 20)
p5 <- f_plotch(cp_g_fortypsm, "growth", "Forest Type", 10)
p6 <- f_plotch(cp_g_disturb, "growth", "Disturbance", 20)
grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2)

p1 <- f_plotch(cp_r_age, "recr", "age", 20)
p2 <- f_plotch(cp_r_ndep, "recr", "N deposition", 20)
p3 <- f_plotch(cp_r_sm, "recr", "Soil Moisture", 20)
p4 <- f_plotch(cp_r_texture, "recr", "Texture", 20)
p5 <- f_plotch(cp_r_fortypsm, "recr", "Forest Type", 10)
p6 <- f_plotch(cp_r_disturb, "recr", "Disturbance", 20)
grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2)

p1 <- f_plotch(cp_s_ndep2, "surv", "N deposition", 20)
p2 <- f_plotch(cp_g_ndep2, "growth", "N deposition", 20)
p3 <- f_plotch(cp_r_ndep2, "recr", "N deposition", 20)
grid.arrange(p1, p2, p3, ncol = 3, nrow = 1)
xtabs(~NDEP2, FIA_plot)


### growth overall
mp_g_all <- lmer(BA_change ~ BA_Nfixer_prop_total + 
                   BA_total1 + 
                   BA_total2 + 
                   BA_Nfixer1 + 
                   BA_nonfixer1 +
                   (1 | STATECD), 
                 data = FIA_plot)
model <- mp_g_all
orig <- final <- ch <- NULL
#untransform coefficients
#mu <- mean(FIA_plot$BA_Nfixer_prop_total, na.rm=T)
#sdev <- sd(FIA_plot$BA_Nfixer_prop_total, na.rm=T)
#mu <- 0.01361491 #same for whole data set
#sdev <- 0.08622454
b0 <- fixef(model)[["(Intercept)"]] #fixed intercept
b1 <- fixef(model)[["BA_Nfixer_prop_total"]] #fixed NCIprop slope
#ints <- b0 - (b1) * (mu/sdev)
ints <- b0
#slopes <- (b1) / sdev
slopes <- b1
re <- cbind.data.frame(ints, slopes)
colnames(re) <- c("int", "slope")
for(i in 1:nrow(re)){
  orig[i] <- re[i,"int"]
  final[i] <- re[i,"int"]+re[i,"slope"]
  ch[i] <- 100*(final[i]-orig[i])/orig[i]
}
cp_g_all <- ch

### recr overall
mp_r_all <- lmer(recyr ~ BA_Nfixer_prop_total+ 
                   BA_total1+ 
                   BA_total2+ 
                   BA_Nfixer1 + 
                   BA_nonfixer1 +
                   (1 | STATECD), 
                 data = FIA_plot)
model <- mp_r_all
orig <- final <- ch <- NULL
#untransform coefficients
#mu <- mean(FIA_plot$BA_Nfixer_prop_total, na.rm=T)
#sdev <- sd(FIA_plot$BA_Nfixer_prop_total, na.rm=T)
#mu <- 0.01361491 #same for whole data set
#sdev <- 0.08622454
b0 <- fixef(model)[["(Intercept)"]] #fixed intercept
b1 <- fixef(model)[["BA_Nfixer_prop_total"]] #fixed NCIprop slope
#ints <- b0 - (b1) * (mu/sdev)
ints <- b0
#slopes <- (b1) / sdev
slopes <- b1
re <- cbind.data.frame(ints, slopes)
colnames(re) <- c("int", "slope")
for(i in 1:nrow(re)){
  orig[i] <- re[i,"int"]
  final[i] <- re[i,"int"]+re[i,"slope"]
  ch[i] <- 100*(final[i]-orig[i])/orig[i]
}
cp_r_all <- ch

### surv overall
mp_s_all <- lmer(survyr ~ BA_Nfixer_prop_total + 
                BA_total1 + 
                BA_total2 + 
                BA_Nfixer2 + 
                BA_nonfixer2 +
                (1 | STATECD), 
              data = FIA_plot)
model <- mp_s_all
orig <- final <- ch <- NULL
#untransform coefficients
#mu <- mean(FIA_plot$BA_Nfixer_prop_total, na.rm=T)
#sdev <- sd(FIA_plot$BA_Nfixer_prop_total, na.rm=T)
#mu <- 0.01361491 #same for whole data set
#sdev <- 0.08622454
b0 <- fixef(model)[["(Intercept)"]] #fixed intercept
b1 <- fixef(model)[["BA_Nfixer_prop_total"]] #fixed NCIprop slope
#ints <- b0 - (b1) * (mu/sdev)
ints <- b0
#slopes <- (b1) / sdev
slopes <- b1
re <- cbind.data.frame(ints, slopes)
colnames(re) <- c("int", "slope")
for(i in 1:nrow(re)){
  orig[i] <- re[i,"int"]
  final[i] <- re[i,"int"]+re[i,"slope"]
  ch[i] <- 100*(final[i]-orig[i])/orig[i]
}
cp_s_all <- ch


## for growth, re-run growth all forest just before making this plot for correct ints and slopes
number_bins <- 50
stats_mp_g_all <- stats.bin(x = FIA_plot$BA_Nfixer_prop_total,
                            y = FIA_plot[ ,"BA_change"],
                            N = number_bins)
par(mar = c(5.1,5.1,4.1,2.1)+0.1,
    cex.main = 1.5,
    cex.lab = 1.5,
    cex.axis = 1.3,
    bty = "n")
plot(seq(1, number_bins, 1),
     stats_mp_g_all$stats["mean", ],
     xaxt = "n",
     xlab = "BA from fixers (%)",
     ylab = "BA increment",
     main = "Plot Level",
     ylim = c(0, 5),
     cex.lab = 1.5,
     col = "black")
par(new = TRUE)
plot(seq(0, 1, 0.1),
     ints + slopes*seq(0, 1, 0.1),
     type = "l",
     ylim = c(0, 5),
     xaxt = "n",
     yaxt = "n",
     ylab = "",
     xlab = "",
     col = "black")
axis(side = 1, at = seq(0, 1, 0.1), labels = seq(0, 100, 10))


###
# correlations
require(psych)
corr.test(FIA_plot$BA_change, FIA_plot$sm, method="pearson")
corr.test(FIA_plot$BA_change, FIA_plot$dep, method="pearson")
corr.test(FIA_plot$BA_change, FIA_plot$LAT)
corr.test(FIA_plot$BA_change, FIA_plot$LON)

corr.test(FIA_plot$recyr, FIA_plot$sm)
corr.test(FIA_plot$recyr, FIA_plot$dep)
corr.test(FIA_plot$recyr, FIA_plot$LAT)
corr.test(FIA_plot$recyr, FIA_plot$LON)

corr.test(FIA_plot$survyr, FIA_plot$sm)
corr.test(FIA_plot$survyr, FIA_plot$dep)
corr.test(FIA_plot$survyr, FIA_plot$LAT)
corr.test(FIA_plot$survyr, FIA_plot$LON)


####
# Figs
# N dep, soil moisture, understory C for all 3 demographic rates
# cp_g_ndep, cp_r_ndep, cp_s_ndep
# cp_g_sm, cp_r_sm, cp_s_sm
# cp_g_uc, cp_r_uc, cp_s_uc

cp_g_ndep$est$type <- "growth"
cp_r_ndep$est$type <- "recruitment"
cp_s_ndep$est$type <- "survival"
cp_ndep <- rbind.data.frame(cp_g_ndep$est, cp_r_ndep$est, cp_s_ndep$est)

p1 <- ggplot(data = cp_ndep, 
               aes(x = factor(grps), y = ch, fill = type, width = 0.65)) +
  geom_bar(stat = "identity", 
           position = position_dodge()) + 
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Deposition Rate (kg N/ha/yr)", y = "% Change", title = "N deposition") +
  geom_hline(yintercept = 0, color = "black") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5, size = 20),
        text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        legend.position = "none")

cp_g_sm$est$type <- "growth"
cp_r_sm$est$type <- "recruitment"
cp_s_sm$est$type <- "survival"
cp_sm <- rbind.data.frame(cp_g_sm$est, cp_r_sm$est, cp_s_sm$est)

p2 <- ggplot(data = cp_sm, 
             aes(x = factor(grps), y = ch, fill = type, width = 0.65)) +
  geom_bar(stat = "identity", 
           position = position_dodge()) + 
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Soil Moisture (m3/m3)", y = "% Change", title = "Soil Moisture") +
  geom_hline(yintercept = 0, color = "black") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5, size = 20),
        text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        legend.position = "none")

cp_g_uc$est$type <- "growth"
cp_r_uc$est$type <- "recruitment"
cp_s_uc$est$type <- "survival"
cp_uc <- rbind.data.frame(cp_g_uc$est, cp_r_uc$est, cp_s_uc$est)

p3 <- ggplot(data = cp_uc, 
             aes(x = factor(grps), y = ch, fill = type, width = 0.65)) +
  geom_bar(stat = "identity", 
           position = position_dodge()) + 
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Understory C", y = "% Change", title = "Understory C") +
  geom_hline(yintercept = 0, color = "black") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5, size = 20),
        text = element_text(size = 20),
        axis.text.x = element_text(size = 20))

grid.arrange(p1, p2, p3, ncol = 3, nrow = 1)


### uncertainty
# http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#predictions-andor-confidence-or-prediction-intervals-on-predictions 

require(MCMCglmm)
m.bayes <- MCMCglmm(survyr ~ scale(BA_Nfixer_prop_total) + 
                      scale(BA_total1) + 
                      scale(BA_total2) + 
                      scale(BA_Nfixer2) + 
                      scale(BA_nonfixer2), 
                    random = ~STATECD + BA_Nfixer_prop_total:texture, 
                    data = FIA_plot, 
                     verbose = FALSE, pl = TRUE, nitt = 500, thin = 1, burnin = 100)

model <- lmer(survyr ~ scale(BA_Nfixer_prop_total) + 
                scale(BA_total1) + 
                scale(BA_total2) + 
                scale(BA_Nfixer2) + 
                scale(BA_nonfixer2) +
                (1 | STATECD) + 
                (0 + BA_Nfixer_prop_total | group), 
              data = datapass) 

summary(m.bayes)
plot(m.bayes)
