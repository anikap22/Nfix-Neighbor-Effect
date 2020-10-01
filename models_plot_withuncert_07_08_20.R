


####################
#
#
# models_plot_withuncert.R
#
# 4/29/19
# Anika Petach
# 
#
#####################

# https://www.r-bloggers.com/quick-and-dirty-notes-on-general-linear-mix-models/

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
require(gridExtra)
require(MuMIn)
require(svMisc)


setwd("/Users/Anika/Documents/GradSchool/FIA_CompetitionModel/")

source("scripts/functions.R")

options(scipen = 6)



FIA_plot <- readRDS("output/FIA_plot_withgroups.RDS")

FIA_plot$t <- FIA_plot$mt2 - FIA_plot$mt1

# define binary groups ---------------------------------------------
FIA_plot <- FIA_plot %>%
  mutate(NDEPg = case_when(dep < 3.64 ~ 0,
                           dep >= 3.64 ~ 1))
FIA_plot$NDEPg <- as.factor(FIA_plot$NDEPg)


FIA_plot <- FIA_plot %>%
  mutate(SMg = case_when(sm < 0.273 ~ 0,
                         sm >= 0.273 ~ 1))
FIA_plot$SMg <- as.factor(FIA_plot$SMg)

FIA_plot <- FIA_plot %>%
  mutate(TEXg = case_when(texture == 1 ~ 0, #clay
                          texture == 2 ~ 2, #silty clay
                          texture == 3 ~ 1, #sandy clay
                          texture == 4 ~ 0, #clay loam
                          texture == 5 ~ 2, #silty clay loam
                          texture == 6 ~ 1, #sandy clay loam
                          texture == 7 ~ 1, #loam
                          texture == 8 ~ 2, #silty loam
                          texture == 9 ~ 1, #sandy loam
                          texture == 10 ~ 2, #silt
                          texture == 11 ~ 1, #loamy sand
                          texture == 12 ~ 1)) #sand
FIA_plot$TEXg <- as.factor(FIA_plot$TEXg)


FIA_plot <- FIA_plot %>%
  mutate(CUNg = case_when(understoryC1 < 0.972 ~ 0,
                          understoryC1 >= 0.972 ~ 1))
FIA_plot$CUNg <- as.factor(FIA_plot$CUNg) 


FIA_plot <- FIA_plot %>%
  mutate(CLg = case_when(litterC1 < 5.113 ~ 0,
                         litterC1 >= 5.113 ~ 1))
FIA_plot$CLg <- as.factor(FIA_plot$CLg)


FIA_plot <- FIA_plot %>%
  mutate(CSOg = case_when(soilC1 < 37.283 ~ 0,
                         soilC1 >= 37.283 ~ 1))
FIA_plot$CSOg <- as.factor(FIA_plot$CSOg)

FIA_plot <- FIA_plot %>%
  mutate(MATg = case_when(MAT < 11.01099 ~ 0, #mean MAT value
                          MAT >= 11.01099 ~ 1))
FIA_plot$MATg <- as.factor(FIA_plot$MATg)

FIA_plot <- FIA_plot %>%
  mutate(MAPg = case_when(MAP < 1064.784 ~ 0, #mean MAP value
                          MAP >= 1064.784 ~ 1))
FIA_plot$MAPg <- as.factor(FIA_plot$MAPg)


# scale covariates -----------------------------------------------
# or models don't converge
#FIA_plot$BA_props <- scale(FIA_plot$BA_Nfixer_prop_total, scale=T, center=T)
#FIA_plot$BA_total1s <- scale(FIA_plot$BA_total1, scale=T, center=T)

# select rows with no NAs
FIA_plot0 <- FIA_plot %>% dplyr::select(BA_change, BA_change_nonfixer, BA_change_Nfixer, BA_Nfixer_prop_total, 
                                        BA_total1, NDEPg, SMg, TEXg, CUNg, CSOg, 
                                        MATg, MAPg, t, pcn1,
                                        STATECD, LAT, LON, YOU0OLD1, recyr, survyr, nrec,
                                        nsurv, BA_total2, BA_Nfixer1, BA_nonfixer1)
FIA_plot0 <- FIA_plot0[complete.cases(FIA_plot0), ]

FIA_plot0$BA_props <- (FIA_plot0$BA_Nfixer_prop_total - mean(FIA_plot0$BA_Nfixer_prop_total))/sd(FIA_plot0$BA_Nfixer_prop_total)
FIA_plot0$BA_total1s <- (FIA_plot0$BA_total1 - mean(FIA_plot0$BA_total1))/sd(FIA_plot0$BA_total1)
FIA_plot0$BA_total2s <- (FIA_plot0$BA_total2 - mean(FIA_plot0$BA_total2))/sd(FIA_plot0$BA_total2)
FIA_plot0$BA_Nfixer1s <- (FIA_plot0$BA_Nfixer1 - mean(FIA_plot0$BA_Nfixer1))/sd(FIA_plot0$BA_Nfixer1)


FIA_plot <- FIA_plot0


# age correlation -------------------------
test <- FIA_plot[FIA_plot$BA_Nfixer1 > 0,]
hist(FIA_plot$stdage, col=rgb(0.1,0.1,0.1,0.5), breaks=50,
     prob = T, xlim=c(0,250),
     main="plots by stand age and plots with N fixers by stand age")
hist(test$stdage, col=rgb(0.8,0.8,0.8,0.5), breaks = 20, 
     prob = T, xlim=c(0,250), add=T)
box()



# model selection ----------------------------------------------
# for BAI get BA_change
# for BAIn get BA_change_nonfixer
# for Recruitment get recyr
# for Survival get survyr
fia.s <- FIA_plot %>% dplyr::select(recyr, BA_Nfixer_prop_total, BA_total1,
                               STATECD)
fia.s <- fia.s[complete.cases(fia.s), ]

# m0 <- lmer(BA_change_nonfixer ~ BA_Nfixer_prop_total + 
#                 BA_total1 + 
#                 BA_total2 + 
#                 BA_Nfixer1 + 
#                 BA_nonfixer1 +
#                 (1 | STATECD) - 1, 
#               data = fia.s)

m0 <- lmer(BA_change ~ 
             BA_Nfixer_prop_total + 
             BA_total1 + 
             (1 | STATECD) - 1, 
           data = fia.s)

AICc(m0)

m0 <- lmer(BA_change ~ 
             BA_Nfixer_prop_total + 
             (1 | STATECD) - 1, 
           data = fia.s)
AICc(m0)

#make response variable BA_change_nonfixer
m0 <- lmer(BA_change_nonfixer ~ 
             BA_Nfixer_prop_total + 
             BA_total1 + 
             (1 | STATECD) - 1, 
           data = fia.s)
AICc(m0)

m0 <- lmer(BA_change_nonfixer ~ 
             BA_Nfixer_prop_total + 
             (1 | STATECD) - 1, 
           data = fia.s)
AICc(m0)

#make response variable survyr
m0 <- lmer(survyr ~ 
             BA_Nfixer_prop_total + 
             BA_total1 + 
             (1 | STATECD) - 1, 
           data = fia.s)
AICc(m0)

m0 <- lmer(survyr ~ 
             BA_Nfixer_prop_total + 
             (1 | STATECD) - 1, 
           data = fia.s)
AICc(m0)

#make response variable recyr
m0 <- lmer(recyr ~ 
             BA_Nfixer_prop_total + 
             BA_total1 + 
             (1 | STATECD) - 1, 
           data = fia.s)
AICc(m0)

m0 <- lmer(recyr ~ 
             BA_Nfixer_prop_total + 
             (1 | STATECD) - 1, 
           data = fia.s)
AICc(m0)

# with BA_change as response:
# Sum Sq Mean Sq NumDF  DenDF    F value Pr(>F)    
# BA_Nfixer_prop_total      0       0     1 183974 2.4480e-01 0.6208    
# BA_total1            341839  113946     3 183987 6.8285e+04 <2e-16 ***
#   BA_total2                 0       0     1 183974 2.4480e-01 0.6208    
# BA_Nfixer1           341839  113946     3 183987 6.8285e+04 <2e-16 ***
#   BA_nonfixer1         341797  170899     2 183996 1.0241e+05 <2e-16 ***

m1 <- lmer(survyr ~ BA_Nfixer_prop_total + 
             BA_total1 + 
             BA_total2 + 
             BA_Nfixer1 + 
             NDEPg +
             BA_Nfixer_prop_total*NDEPg +
             (1 | STATECD), 
           data = FIA_plot)
m2 <- lmer(survyr ~ BA_Nfixer_prop_total + 
             BA_Nfixer1 + 
             NDEPg +
             BA_Nfixer_prop_total*NDEPg +
             (1 | STATECD), 
           data = FIA_plot)
drop1(m1, test = "none")
step(m1)

m2 <- lmer(survyr ~ BA_Nfixer_prop_total + 
             BA_total1 + 
             BA_total2 + 
             (1 | STATECD) - 1, 
           data = fia.s)
drop1(m2, test = "none")

m3 <- lmer(survyr ~ BA_Nfixer_prop_total + 
             BA_total2 + 
             (1 | STATECD) - 1, 
           data = fia.s)
drop1(m3, test = "none")

# try with derdge instead
library(MuMIn)
options(na.action = na.fail)
models <- dredge(m0, beta = F, evaluate = T, rank = AICc)

sel.table <- as.data.frame(models)
sel.table


# model validation ---------------------------------------------------------

m.g2 <- lmer(BA_change ~ BA_Nfixer_prop_total + 
              BA_total1 + 
              (1 | STATECD) - 1, 
            data = FIA_plot,
            na.action = na.omit)
plot(m.g2, add.smooth=FALSE)

# growth
op <- par(mfrow = c(2,2), mar = c(4, 4, 2, 2))
plot(m.g, which = c(1), col = 1, add.smooth = FALSE, caption = " ")
plot(FIA_plot$BA_Nfixer_prop_total, resid(m.g), xlab = "BA prop", ylab = "Residuals")
plot(FIA_plot$BA_total1, resid(m.g), xlab = "BA total 1", ylab = "Residuals")
plot(FIA_plot$BA_Nfixer1, resid(m.g), xlab = "BA Nfixer 1", ylab = "Residuals")
par(op)
hist(resid(m.g))

m.r <- lmer(recyr ~ BA_Nfixer_prop_total + 
               BA_total1 + 
               (1 | STATECD) - 1, 
             data = FIA_plot,
             na.action = na.omit)
plot(m.r, add.smooth=FALSE)

m.r <- lmer(log(recyr+0.01) ~ BA_Nfixer_prop_total + 
              BA_total1 + 
              (1 | STATECD) - 1, 
            data = FIA_plot,
            na.action = na.omit)

FIA_plot$t <- FIA_plot$mt2 - FIA_plot$mt1

samp <- sample(1:nrow(FIA_plot), 1000, replace=F)
temp <- FIA_plot[samp,]


# poisson with rate data: https://rpubs.com/kaz_yos/poisson
m.r <- glmer(rec ~ BA_Nfixer_prop_total +
               BA_total1 +
               (1 | STATECD) -1,
             offset = log(t),
             family = poisson(link="log"),
             data = temp)


# recr
plot(m.r, add.smooth = FALSE)
plot(FIA_plot$BA_Nfixer_prop_total, resid(m.r), xlab = "BA prop", ylab = "Residuals")
plot(FIA_plot$BA_total1, resid(m.r), xlab = "BA total 1", ylab = "Residuals")
plot(FIA_plot$BA_Nfixer1, resid(m.r), xlab = "BA Nfixer 1", ylab = "Residuals")
hist(resid(m.r))

m.s <- lmer(survyr ~ BA_Nfixer_prop_total + 
               BA_total1 + 
               (1 | STATECD) - 1, 
             data = FIA_plot,
             na.action = na.omit)
plot(m.s, add.smooth=FALSE)

m.s <- lmer(log(survyr+0.01) ~ BA_Nfixer_prop_total + 
              BA_total1 + 
              (1 | STATECD) - 1, 
            data = FIA_plot,
            na.action = na.omit)

m.s <- glmer(surv ~ BA_Nfixer_prop_total +
               BA_total1 +
               (1 | STATECD) -1,
             offset = log(t),
             family = poisson(link="logit"),
             data = temp)

# surv
par(mfrow = c(2,2))
plot(m.s, add.smooth = FALSE)
plot(FIA_plot$BA_Nfixer_prop_total, resid(m.s), xlab = "BA prop", ylab = "Residuals")
plot(FIA_plot$BA_total1, resid(m.s), xlab = "BA total 1", ylab = "Residuals")
plot(FIA_plot$BA_Nfixer1, resid(m.s), xlab = "BA Nfixer 1", ylab = "Residuals")
hist(resid(m.s))
dev.off()

# use model - whole forest ----------------------------------


#BAI
m <- lmer(BA_change ~ BA_Nfixer_prop_total + 
            BA_total1 + 
            (1 | STATECD) - 1, 
          data = FIA_plot,
          na.action = na.omit)
r.squaredGLMM(m)
confint(m)

#initial = summary(m)$coefficients["(Intercept)","Estimate"]
#final = initial + summary(m)$coefficients["BA_Nfixer_prop_total","Estimate"]
#pch = 100*(final-initial)/initial
#pch

cat("Whole forest plot level\n", file = "output/table1.txt", append = TRUE)
cat("BAI\n", file = "output/table1.txt", append = TRUE)
out <- f_nonpar_boot(FIA_plot, "growth", "NA", 1000)
capture.output(pct_change_func_byrow_forest(out), 
               file = "output/table1.txt", append = TRUE)


#BAIn
m <- lmer(BA_change_nonfixer ~ BA_Nfixer_prop_total + 
            BA_total1 + 
            (1 | STATECD) - 1, 
          data = FIA_plot,
          na.action = na.omit)
r.squaredGLMM(m)
confint(m)

cat("BAIn\n", file = "output/table1.txt", append = TRUE)
out <- f_nonpar_boot(FIA_plot, "growthnf", "NA", 1000)
capture.output(pct_change_func_byrow_forest(out), 
               file = "output/table1.txt", append = TRUE)


#Recr
m <- lmer(recyr ~ BA_Nfixer_prop_total +
            BA_total1 +
            (1 | STATECD) - 1,
          data = FIA_plot)
r.squaredGLMM(m)
confint(m)

cat("Recruitment\n", file = "output/table1.txt", append = TRUE)
out <- f_nonpar_boot(FIA_plot, "recr", "NA", 1000)
capture.output(pct_change_func_byrow_forest(out), 
               file = "output/table1.txt", append = TRUE)

#Surv
m <- lmer(survyr ~ BA_Nfixer_prop_total +
            BA_total1 +
            (1 | STATECD) - 1,
          data = FIA_plot)
r.squaredGLMM(m)
confint(m)

cat("Survival\n", file = "output/table1.txt", append = TRUE)
out <- f_nonpar_boot(FIA_plot, "surv", "NA", 1000)
capture.output(pct_change_func_byrow_forest(out), 
               file = "output/table1.txt", append = TRUE)

#BAIf
m <- lmer(BA_change_Nfixer ~ BA_Nfixer_prop_total + 
            BA_total1 + 
            (1 | STATECD) - 1, 
          data = FIA_plot,
          na.action = na.omit)
r.squaredGLMM(m)
confint(m)

out <- f_nonpar_boot(FIA_plot, "growthf", "NA", 1000)
pct_change_func_byrow_forest(out)




# use model -----------------------------------------------------
# note: texture (TEXg) requires setting up function for 3 groups...model off individual func
##### with nonparametric bootstrap
out_ndep <- readRDS("output/boots/out_ndep_bai_1000bs.RDS")
out_age <- readRDS("output/boots/out_age_bai_1000bs.RDS")
out_sm <- readRDS("output/boots/out_sm_bai_1000bs.RDS")
out_cun <- readRDS("output/boots/out_cun_bai_1000bs.RDS")
out_cso <- readRDS("output/boots/out_cso_bai_1000bs.RDS")
out_mat <- readRDS("output/boots/out_mat_bai_1000bs.RDS")
out_map <- readRDS("output/boots/out_map_bai_1000bs.RDS")


nboot = 100
# growth
out_ndep <- f_nonpar_boot(FIA_plot,"growth","NDEPg",nboot)
pch_ndep <- pct_change_func_byrow(out_ndep)
f <- xtabs(~NDEPg, FIA_plot)
pch_ndep$Freq <- rbind(f[[1]], f[[2]]) 

out_age <- f_nonpar_boot(FIA_plot,"growth","YOU0OLD1",nboot)
pch_age <- pct_change_func_byrow(out_age)
f <- xtabs(~YOU0OLD1, FIA_plot)
pch_age$Freq <- rbind(f[[1]], f[[2]]) 

out_sm <- f_nonpar_boot(FIA_plot,"growth","SMg",nboot)
pch_sm <- pct_change_func_byrow(out_sm)
f <- xtabs(~SMg, FIA_plot)
pch_sm$Freq <- rbind(f[[1]], f[[2]]) 

out_cun <- f_nonpar_boot(FIA_plot,"growth","CUNg",nboot)
pch_cun <- pct_change_func_byrow(out_cun)
f <- xtabs(~CUNg, FIA_plot)
pch_cun$Freq <- rbind(f[[1]], f[[2]]) 

out_cso <- f_nonpar_boot(FIA_plot,"growth","CSOg",nboot)
pch_cso <- pct_change_func_byrow(out_cso)
f <- xtabs(~CSOg, FIA_plot)
pch_cso$Freq <- rbind(f[[1]], f[[2]]) 

# out_mat <- f_nonpar_boot(FIA_plot,"growth","MATg",nboot)
# pch_mat <- pct_change_func_byrow(out_mat)
# f <- xtabs(~MATg, FIA_plot)
# pch_mat$Freq <- rbind(f[[1]], f[[2]]) 
# 
# out_map <- f_nonpar_boot(FIA_plot,"growth","MAPg",nboot)
# pch_map <- pct_change_func_byrow(out_map)
# f <- xtabs(~MAPg, FIA_plot)
# pch_map$Freq <- rbind(f[[1]], f[[2]]) 

pch_ndep$var <- "N deposition"
pch_age$var <- "age"
pch_sm$var <- "soil moisture"
pch_cun$var <- "understory C"
pch_cso$var <- "soil organic C"
#pch_mat$var <- "MAT"
#pch_map$var <- "MAP"
baidat <- rbind(pch_ndep, pch_age, pch_sm, pch_cun, pch_cso, pch_mat, pch_map)
baidat$type <- "BAI"
baidat$level <- as.factor(baidat$group)
plotdat <- baidat
#plotdat <- rbind(plotdat, survdat)


saveRDS(out_ndep, "output/boots/out_ndep_bai_1000bs_07_20.RDS")
saveRDS(out_age, "output/boots/out_age_bai_1000bs_07_20.RDS")
saveRDS(out_sm, "output/boots/out_sm_bai_1000bs_07_20.RDS")
saveRDS(out_cun, "output/boots/out_cun_bai_1000bs_07_20.RDS")
saveRDS(out_cso, "output/boots/out_cso_bai_1000bs_07_20.RDS")
saveRDS(out_mat, "output/boots/out_mat_bai_1000bs_07_20.RDS")
saveRDS(out_map, "output/boots/out_map_bai_1000bs_07_20.RDS")


# recruitment
out_ndep <- f_nonpar_boot(FIA_plot,"recr","NDEPg",nboot)
pch_ndep <- pct_change_func_byrow(out_ndep)

out_age <- f_nonpar_boot(FIA_plot,"recr","YOU0OLD1",nboot)
pch_age <- pct_change_func_byrow(out_age)

out_sm <- f_nonpar_boot(FIA_plot,"recr","SMg",nboot)
pch_sm <- pct_change_func_byrow(out_sm)

out_cun <- f_nonpar_boot(FIA_plot,"recr","CUNg",nboot)
pch_cun <- pct_change_func_byrow(out_cun)

out_cso <- f_nonpar_boot(FIA_plot,"recr","CSOg",nboot)
pch_cso <- pct_change_func_byrow(out_cso)

# out_mat <- f_nonpar_boot(FIA_plot,"recr","MATg",nboot)
# pch_mat <- pct_change_func_byrow(out_mat)
# 
# out_map <- f_nonpar_boot(FIA_plot,"recr","MAPg",nboot)
# pch_map <- pct_change_func_byrow(out_map)

pch_ndep$var <- "N deposition"
pch_age$var <- "age"
pch_sm$var <- "soil moisture"
pch_cun$var <- "understory C"
pch_cso$var <- "soil organic C"
#pch_mat$var <- "MAT"
#pch_map$var <- "MAP"
recrdat <- rbind(pch_ndep, pch_age, pch_sm, pch_cun, pch_cso)
recrdat$type <- "Recruitment"
recrdat$level <- as.factor(recrdat$group)
plotdat <- rbind(plotdat, recrdat)

saveRDS(out_ndep, "output/boots/out_ndep_recr_1000bs_07_20.RDS")
saveRDS(out_age, "output/boots/out_age_recr_1000bs_07_20.RDS")
saveRDS(out_sm, "output/boots/out_sm_recr_1000bs_07_20.RDS")
saveRDS(out_cun, "output/boots/out_cun_recr_1000bs_07_20.RDS")
saveRDS(out_cso, "output/boots/out_cso_recr_1000bs_07_20.RDS")
saveRDS(out_mat, "output/boots/out_mat_recr_1000bs_07_20.RDS")
saveRDS(out_map, "output/boots/out_map_recr_1000bs_07_20.RDS")

# survival
out_ndep <- f_nonpar_boot(FIA_plot,"surv","NDEPg",nboot)
pch_ndep <- pct_change_func_byrow(out_ndep)

out_age <- f_nonpar_boot(FIA_plot,"surv","YOU0OLD1",nboot)
pch_age <- pct_change_func_byrow(out_age)

out_sm <- f_nonpar_boot(FIA_plot,"surv","SMg",nboot)
pch_sm <- pct_change_func_byrow(out_sm)

out_cun <- f_nonpar_boot(FIA_plot,"surv","CUNg",nboot)
pch_cun <- pct_change_func_byrow(out_cun)

out_cso <- f_nonpar_boot(FIA_plot,"surv","CSOg",nboot)
pch_cso <- pct_change_func_byrow(out_cso)

# out_mat <- f_nonpar_boot(FIA_plot,"surv","MATg",nboot)
# pch_mat <- pct_change_func_byrow(out_mat)
# 
# out_map <- f_nonpar_boot(FIA_plot,"surv","MAPg",nboot)
# pch_map <- pct_change_func_byrow(out_map)

pch_ndep$var <- "N deposition"
pch_age$var <- "age"
pch_sm$var <- "soil moisture"
pch_cun$var <- "understory C"
pch_cso$var <- "soil organic C"
#pch_mat$var <- "MAT"
#pch_map$var <- "MAP"
survdat <- rbind(pch_ndep, pch_age, pch_sm, pch_cun, pch_cso)
survdat$type <- "Survival"
survdat$level <- as.factor(survdat$group)
plotdat <- rbind(plotdat, survdat)

saveRDS(out_ndep, "output/boots/out_ndep_surv_1000bs_07_20.RDS")
saveRDS(out_age, "output/boots/out_age_surv_1000bs_07_20.RDS")
saveRDS(out_sm, "output/boots/out_sm_surv_1000bs_07_20.RDS")
saveRDS(out_cun, "output/boots/out_cun_surv_1000bs_07_20.RDS")
saveRDS(out_cso, "output/boots/out_cso_surv_1000bs_07_20.RDS")
saveRDS(out_mat, "output/boots/out_mat_surv_1000bs_07_20.RDS")
saveRDS(out_map, "output/boots/out_map_surv_1000bs_07_20.RDS")


# growth of non-fixers (BAI-nonfixer)
out_ndep <- f_nonpar_boot(FIA_plot,"growthnf","NDEPg",nboot)
pch_ndep <- pct_change_func_byrow(out_ndep)

out_age <- f_nonpar_boot(FIA_plot,"growthnf","YOU0OLD1",nboot)
pch_age <- pct_change_func_byrow(out_age)

out_sm <- f_nonpar_boot(FIA_plot,"growthnf","SMg",nboot)
pch_sm <- pct_change_func_byrow(out_sm)

out_cun <- f_nonpar_boot(FIA_plot,"growthnf","CUNg",nboot)
pch_cun <- pct_change_func_byrow(out_cun)

out_cso <- f_nonpar_boot(FIA_plot,"growthnf","CSOg",nboot)
pch_cso <- pct_change_func_byrow(out_cso)

# out_mat <- f_nonpar_boot(FIA_plot,"growthnf","MATg",nboot)
# pch_mat <- pct_change_func_byrow(out_mat)
# 
# out_map <- f_nonpar_boot(FIA_plot,"growthnf","MAPg",nboot)
# pch_map <- pct_change_func_byrow(out_map)

pch_ndep$var <- "N deposition"
pch_age$var <- "age"
pch_sm$var <- "soil moisture"
pch_cun$var <- "understory C"
pch_cso$var <- "soil organic C"
pch_mat$var <- "MAT"
pch_map$var <- "MAP"
baindat <- rbind(pch_ndep, pch_age, pch_sm, pch_cun, pch_cso)
baindat$type <- "BAIn"
baindat$level <- as.factor(baindat$group)
plotdat <- rbind(plotdat, baindat)

saveRDS(out_ndep, "output/boots/out_ndep_bain_1000bs.RDS")
saveRDS(out_age, "output/boots/out_age_bain_1000bs.RDS")
saveRDS(out_sm, "output/boots/out_sm_bain_1000bs.RDS")
saveRDS(out_cun, "output/boots/out_cun_bain_1000bs.RDS")
saveRDS(out_cso, "output/boots/out_cso_bain_1000bs.RDS")
#saveRDS(out_mat, "output/boots/out_mat_bain_1000bs.RDS")
#saveRDS(out_map, "output/boots/out_map_bain_1000bs.RDS")

plotdat <- plotdat %>% mutate(Freq = case_when(var == 'N deposition' & level == 0 ~ pch_ndep$Freq[1],
                                            var == 'N deposition' & level == 1 ~ pch_ndep$Freq[2],
                                            var == 'age' & level == 0 ~ pch_age$Freq[1],
                                            var == 'age' & level == 1 ~ pch_age$Freq[2],
                                            var == 'soil moisture' & level == 0 ~ pch_sm$Freq[1],
                                            var == 'soil moisture' & level == 1 ~ pch_sm$Freq[2]))


saveRDS(plotdat, "output/boots/plotdat_07_20.RDS")

require(kableExtra)
kable(plotdat) %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  save_kable(file = "output/fig2_table.html", self_contained = T)


# ----------------------------
plotdat <- readRDS("output/boots/plotdat_07_20.RDS")
plotdat <- plotdat %>% mutate(Group = case_when(group == 0 ~ "low/young", #low
                                              group == 1 ~ "high/old")) #high ,group ==2 ~ "none"
plotdat$type <- ordered(plotdat$type, levels = c("BAI", "BAIn", "Recruitment", "Survival"))
plotdat$type2 <- factor(plotdat$type, labels = c("BAI[NFE]", expression(BAI[paste(n,',',NFE)]), "R[NFE]", "S[NFE]"))
plotdat <- plotdat %>% filter(var %in% c("N deposition","stand age","soil moisture",
                                         "MAT", "MAP"))
plotdat$var <- factor(plotdat$var, levels = c("MAT", "MAP", "soil moisture", 
                                               "N deposition", "stand age"))
levels(plotdat$var)[levels(plotdat$var)=="age"] <- "stand age"

plotdat$Group <- factor(plotdat$Group, levels=c('low/young','high/old'))

plotdat <- plotdat %>% filter(var %in% c('N deposition', 'stand age', 'soil moisture'))

plot <- ggplot(data = plotdat, 
               aes(x = factor(var), y = mean, fill = Group, width = 0.65)) +
  geom_bar(stat = "identity", 
           position = position_dodge()) + 
  coord_flip() + 
  facet_wrap(~type2, labeller = label_parsed) +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "", y = "N-fixer effect (%)", subtitle = "", title = "",
       caption=" ") +
  geom_hline(yintercept = 0, color = "black") +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.65)) +
  #geom_text(aes(label=round(Freq/1000, 1)), position=position_dodge(0.9), vjust=-0.2) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5, size = 20),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, vjust = 0.1, hjust = 1, size = 20),
        axis.text.y = element_text(size = 20),
        plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0)),
        axis.line = element_line(size=1, colour = "black"),
        legend.title = element_blank()) +
  guides(fill = guide_legend(reverse = TRUE))
plot

png("output/fig3_plot_07_20.png", width = 12, height = 7, units = 'in', res = 300)
par(bty = 'n')
plot
dev.off()




plotdat <- plotdat %>% mutate(level = case_when(group == 0 ~ "low",
                                                 group == 1 ~ "high"))


plot <- ggplot(data = plotdat, 
               aes(x = factor(var), y = pchc, fill = level, width = 0.65)) +
  geom_bar(stat = "identity", 
           position = position_dodge()) + 
  coord_flip() + 
  facet_wrap(~type) +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Group", y = "% Change", subtitle = "type", title = "% change in demographic rate",
       caption=" ") +
  geom_hline(yintercept = 0, color = "black") +
  geom_errorbar(aes(ymin=pchl, ymax=pchu), width=.2,
                position=position_dodge(.9)) +
  #geom_text(aes(label=round(Freq/1000, 1)), position=position_dodge(0.9), vjust=-0.2) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5, size = 20),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, vjust = 0.1, hjust = 1, size = 20),
        plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0)),
        axis.line = element_line(size=1, colour = "black")) 



# calc p values for whether groups are different -------------------------------------
# get SD from 95% CI: https://handbook-5-1.cochrane.org/chapter_7/7_7_3_2_obtaining_standard_deviations_from_standard_errors_and.htm
# SD = sqrt(Freq)(UL-LL)/3.92
# get p-value from mean and SD: https://www.cyclismo.org/tutorial/R/pValues.html
library(BSDA)

pvalfunc <- function(data){
  sd0 <- sqrt(data$Freq[data$group==0])*(data$upper[data$group==0]-data$lower[data$group==0])/3.92
  sd1 <- sqrt(data$Freq[data$group==1])*(data$upper[data$group==1]-data$lower[data$group==1])/3.92
  
  out <- tsum.test(mean.x = data$mean[data$group==0], s.x = sd0, n.x = data$Freq[data$group==0],
                   mean.y = data$mean[data$group==1], s.y = sd1, n.y = data$Freq[data$group==1])
  return(out)
  
}

cat("T-Tests\n", file = "output/fig2_ttests.txt", append = TRUE)
cat("BAI\n", file = "output/fig2_ttests.txt", append = TRUE)
cat("N dep, age, soil moisture \n", file = "output/fig2_ttests.txt", append = TRUE)
capture.output(pvalfunc(plotdat[plotdat$type=="BAI" & plotdat$var=="N deposition", ]),
               pvalfunc(plotdat[plotdat$type=="BAI" & plotdat$var=="stand age", ]),
               pvalfunc(plotdat[plotdat$type=="BAI" & plotdat$var=="soil moisture", ]),
               file = "output/fig2_ttests.txt", append = TRUE)
cat("Recruitment\n", file = "output/fig2_ttests.txt", append = TRUE)
cat("N dep, age, soil moisture \n", file = "output/fig2_ttests.txt", append = TRUE)
capture.output(pvalfunc(plotdat[plotdat$type=="Recruitment" & plotdat$var=="N deposition", ]),
               pvalfunc(plotdat[plotdat$type=="Recruitment" & plotdat$var=="stand age", ]),
               pvalfunc(plotdat[plotdat$type=="Recruitment" & plotdat$var=="soil moisture", ]),
               file = "output/fig2_ttests.txt", append = TRUE)
cat("Survival\n", file = "output/fig2_ttests.txt", append = TRUE)
cat("N dep, age, soil moisture \n", file = "output/fig2_ttests.txt", append = TRUE)
capture.output(pvalfunc(plotdat[plotdat$type=="Survival" & plotdat$var=="N deposition", ]),
               pvalfunc(plotdat[plotdat$type=="Survival" & plotdat$var=="stand age", ]),
               pvalfunc(plotdat[plotdat$type=="Survival" & plotdat$var=="soil moisture", ]),
               file = "output/fig2_ttests.txt", append = TRUE)
cat("BAIn\n", file = "output/fig2_ttests.txt", append = TRUE)
cat("N dep, age, soil moisture \n", file = "output/fig2_ttests.txt", append = TRUE)
capture.output(pvalfunc(plotdat[plotdat$type=="BAIn" & plotdat$var=="N deposition", ]),
               pvalfunc(plotdat[plotdat$type=="BAIn" & plotdat$var=="stand age", ]),
               pvalfunc(plotdat[plotdat$type=="BAIn" & plotdat$var=="soil moisture", ]),
               file = "output/fig2_ttests.txt", append = TRUE)

# pvals: #bai, #bain, #r, #s
pvalfunc(pch_ndep) #0.01018  #0.001287  #0.6616  #0.05459
pvalfunc(pch_age) #0.1143  #0.07273  0.2778  #0.2999
pvalfunc(pch_sm) #0.00414  #0.004122  #0.004958  #0.3086
pvalfunc(pch_cun) #0.104  # 0.0005738  #0.1292  #0.9286
pvalfunc(pch_cso) #0.7982  0.7285  #0.6699  #0.7592

pvalfunc(ndep_rs)
pvalfunc(age_rs)
pvalfunc(sm_rs)
pvalfunc(cun_rs)
pvalfunc(cso_rs)

pvalfunc(ndep_ms)
pvalfunc(age_ms)
pvalfunc(sm_ms)
pvalfunc(cun_ms)
pvalfunc(cso_ms)

pvalfunc(ndep_gsn)
pvalfunc(age_gsn)
pvalfunc(sm_gsn)
pvalfunc(cun_gsn)
pvalfunc(cso_gsn)

#check weighted means
wm <- function(data){
  total <- sum(data$Freq)
  w0 <- data$Freq[data$group==0]/total
  w1 <- data$Freq[data$group==1]/total
  out <- data$pchc[data$group==0]*w0 + data$pchc[data$group==1]*w1
  return(out)
}

wm(ndep_ms)
wm(age_ms)
wm(sm_ms)
wm(cun_ms)
wm(cso_ms)






# mapping ------------------------------------------

# f_map_plot from functions.R
map_gs_plot <- f_map_plot(FIA_plot, "growth")
saveRDS(map_gs_plot, "output/map_gs_plot_07_20.RDS")
m1 <- f_mapplot(map_gs_plot$pct.grid, "BAI", "plot", "percent")

map_rs_plot <- f_map_plot(FIA_plot, "recr")
saveRDS(map_rs_plot, "output/map_rs_plot_07_20.RDS")
m2 <- f_mapplot(map_rs_plot$pct.grid, "recruitment", "plot", "percent")

map_ms_plot <- f_map_plot(FIA_plot, "surv")
saveRDS(map_ms_plot, "output/map_ms_plot_07_20.RDS")
m3 <- f_mapplot(map_ms_plot$pct.grid, "survival", "plot", "percent")

map_bain_plot <- f_map_plot(FIA_plot, "bain")
saveRDS(map_bain_plot, "output/map_bain_plot_07_20.RDS")
f_mapplot(map_bain_plot$pct.grid, "BAIn", "plot", "percent")

# Fig 1
require(fields)
source('scripts/functions.R')
map_gs_plot <- readRDS("output/map_gs_plot_07_20.RDS")
map_rs_plot <- readRDS("output/map_rs_plot_07_20.RDS")
map_ms_plot <- readRDS("output/map_ms_plot_07_20.RDS")
gs_noout <- as.vector(map_gs_plot$pct.grid)
#gs_noout <- gs_noout[gs_noout > -1000 & gs_noout < 1000] #216 long
gs_trunc <- as.matrix(map_gs_plot$pct.grid)
gs_trunc[gs_trunc > 100] <- 100
gs_trunc2 <- as.matrix(map_gs_plot$pct.grid)
gs_trunc2[gs_trunc2 > 10000] <- 10000
gs_trunc2[gs_trunc2 < -10000] <- -10000
rs_noout <- as.vector(map_rs_plot$pct.grid)
#rs_noout <- rs_noout[rs_noout > -1000 & rs_noout < 1000] #231 long
rs_trunc <- as.matrix(map_rs_plot$pct.grid)
rs_trunc[rs_trunc > 100] <- 100
ms_noout <- as.vector(map_ms_plot$pct.grid)
#ms_noout <- ms_noout[ms_noout > -1000 & ms_noout < 1000] #260 long
ms_trunc <- as.matrix(map_ms_plot$pct.grid)
ms_trunc[ms_trunc > 100] <- 100
lon.list <- readRDS("output/lon_list.RDS") #made from the first section in mapping function
lat.list <- readRDS("output/lat_list.RDS")

png("output/fig1_plot_07_20.png", width = 12, height = 7, units = 'in', res = 300)
par(mfrow=c(2,3),
    #mai = c(1, 0.1, 0.1, 0.1),
    mai = c(0.6, 0.8, 0.8, 1),
    bty = 'n')
f_mapplot(gs_trunc, "BAI", "plot-scale", "percent") #used to be map_gs_plot$pct.grid
f_mapplot(rs_trunc, "R", "plot-scale", "percent")
f_mapplot(ms_trunc, "S", "plot-scale", "percent")
hist(gs_noout, breaks=100, xlab = "% change in BAI", prob=T, 
     main=expression('plot-scale BAI'[NFE]), cex.lab = 2, cex.main = 2, cex.axis = 2,
     xlim=c(-1000,1000), yaxt='n')
#top <- round(max(density(gs_noout,na.rm=T)$y), digits=3)
#maxden = max(density(gs_noout,na.rm=T)$y)
axis(side = 2, at = c(0,0.006), labels = c(0,0.006), cex.axis=2)
abline(v = 1.1, col = "blue", lwd = 4) #value from table 1
hist(rs_noout, breaks=100, xlab = "% change in recruitment rate", prob=T, 
     main=expression('plot-scale R'[NFE]), cex.lab = 2, cex.main = 2, cex.axis = 2,
     xlim=c(-1000,1000), yaxt='n')
#top <- round(max(density(rs_noout,na.rm=T)$y), digits=3)
axis(side = 2, at = c(0,0.006), labels = c(0,0.006), cex.axis=2)
abline(v = 59.9, col = "blue", lwd = 4)
hist(ms_noout, breaks=50, xlab = "% change in survival rate", prob=T, 
     main=expression('plot-scale S'[NFE]), cex.lab = 2, cex.main = 2, cex.axis = 2,
     xlim=c(-1000,1000), yaxt='n')
#top <- round(max(density(ms_noout,na.rm=T)$y), digits=3)
axis(side = 2, at = c(0,0.1), labels = c(0,0.1), cex.axis=2)
abline(v = -7.3, col = "blue", lwd = 4)
dev.off()

### Fig S2 
#x is lon.list, y is lat.list
f_geog <- function(nfe, type){
  library(reshape2)
  require(raster)
  require(sp)
  
  temp <- nfe
  rownames(temp) <- lon.list
  colnames(temp) <- lat.list
  temp1 <- setNames(reshape2::melt(temp), c('lon', 'lat', 'values'))
  temp1 <- temp1[!is.na(temp1$values), ]
  
  if(type == "BAI"){
    ylab=expression('BAI'[NFE])
    title1 = "a"; title2 = "b"; title3 = "c"; title4 = "d"
  }else if(type == "Recr"){
    ylab=expression('R'[NFE])
    title1 = "e"; title2 = "f"; title3 = "g"; title4 = "h"
  }else if(type == "Surv"){
    ylab=expression('S'[NFE])
    title1 = "i"; title2 = "j"; title3 = "k"; title4 = "l"
  }
  
  lmlon <- lm(temp1$lon ~ temp1$values)
  summary(lmlon)
  r2 <- summary(lmlon)[8]$r.squared
  pltlon <- ggplot(data = temp1, aes(x = lon, y = values)) +
    geom_point() +
    labs(x = "Longitude", y = ylab) +
    annotate("text", x = min(temp1$lon,na.rm=T)+5, y = max(temp1$values,na.rm=T)-1000, 
             label = bquote(italic(R)^2 == .(round(r2, digits = 4)))) +
    theme_classic() +
    theme(strip.text.x = element_text(size = 18),
          strip.text.y = element_text(size = 18),
          text = element_text(size = 20)) +
    ggtitle(title2)
  
  lmlat <- lm(temp1$lat ~ temp1$values)
  summary(lmlat)
  r2 <- summary(lmlat)[8]$r.squared
  pltlat <- ggplot(data = temp1, aes(x = lat, y = values)) +
    geom_point() +
    labs(x = "Latitude", y = ylab) +
    annotate("text", x = min(temp1$lat,na.rm=T)+2, y = max(temp1$values,na.rm=T)-1000,
             label = bquote(italic(R)^2 == .(round(r2, digits = 4)))) +
    theme_classic() +
    theme(strip.text.x = element_text(size = 18),
          strip.text.y = element_text(size = 18),
          text = element_text(size = 20),
          plot.margin = unit(c(0.5,1,0.5,0.5),"cm")) + #top, right, bottom, and left
    ggtitle(title1)

  #get MAT and MAP (code snippet from Ben)
  clim <- getData("worldclim",var="bio",res=10)
  clim <- clim[[c(1,12)]]
  names(clim) <- c("MAT","MAP")
  #coords <- DATA[,c("longitude","latitude")]
  coords <- temp1[,c("lon","lat")]
  points <- SpatialPoints(coords,proj4string=clim@crs)
  climvals <- extract(clim,points)
  temp1 <- cbind(temp1,climvals)
  temp1$MAT <- (temp1$MAT/10)
  
  lmmap <- lm(temp1$MAP ~ temp1$values)
  r2 <- summary(lmmap)[8]$r.squared
  pltmap <- ggplot(data = temp1, aes(x = MAP, y = values)) +
    geom_point() +
    labs(x = "MAP", y = ylab) +
    annotate("text", x = min(temp1$MAP,na.rm=T)+250, y = max(temp1$values,na.rm=T)-1000,
             label = bquote(italic(R)^2 == .(round(r2, digits = 4)))) +
    theme_classic() +
    theme(strip.text.x = element_text(size = 18),
          strip.text.y = element_text(size = 18),
          text = element_text(size = 20),
          plot.margin = unit(c(0.5,1,0.5,0.5),"cm")) + #top, right, bottom, and left
    ggtitle(title3)
  
  lmmat <- lm(temp1$MAT ~ temp1$values)
  summary(lmmat)
  r2 <- summary(lmmat)[8]$r.squared
  pltmat <- ggplot(data = temp1, aes(x = MAT, y = values)) +
    geom_point() +
    labs(x = "MAT", y = ylab) +
    annotate("text", x = min(temp1$MAT,na.rm=T)+3, y = max(temp1$values,na.rm=T)-1000,
             label = bquote(italic(R)^2 == .(round(r2, digits = 4)))) +
    theme_classic() +
    theme(strip.text.x = element_text(size = 18),
          strip.text.y = element_text(size = 18),
          text = element_text(size = 20)) +
    ggtitle(title4)
  
  return(list(pltlat, pltlon, pltmap, pltmat))
  
}

#BAIgeog <- f_geog(gs_trunc2, "BAI") -- function not working ... instead use nfe = map_gs_plot$pct.grid etc.
BAIgeog <- f_geog(map_gs_plot$pct.grid, "BAI")
RECgeog <- f_geog(map_rs_plot$pct.grid, "Recr")
SURVgeog <- f_geog(map_ms_plot$pct.grid, "Surv")

png("output/figS2a_geog.png", width = 12, height = 7, units = 'in', res = 300)
grid.arrange(pltlat, pltlon, pltmap, pltmat, ncol=2, 
             top=textGrob('Plot-scale N-Fixer Effect (NFE)' ,gp=gpar(fontsize=20)))
dev.off()
png("output/figS2b_geog.png", width = 13, height = 7, units = 'in', res = 300)
grid.arrange(pltlat, pltlon, pltmap, pltmat, ncol=2, 
             top=textGrob("Plot-scale N-Fixer Effect (NFE)",gp=gpar(fontsize=20)))
dev.off()
png("output/figS2c_geog.png", width = 12, height = 7, units = 'in', res = 300)
grid.arrange(pltlat, pltlon, pltmap, pltmat, ncol=2, 
             top=textGrob("Plot-scale N-Fixer Effect (NFE)",gp=gpar(fontsize=20)))
dev.off()

##Plot S3 - Maps of MAT, MAP (do we want SM, N dep?)
temp <- map_gs_plot$pct.grid
rownames(temp) <- lon.list
colnames(temp) <- lat.list
temp1 <- setNames(reshape2::melt(temp), c('lon', 'lat', 'values'))
temp1 <- temp1[!is.na(temp1$values), ]
clim <- getData("worldclim",var="bio",res=10)
clim <- clim[[c(1,12)]]
names(clim) <- c("MAT","MAP")
#coords <- DATA[,c("longitude","latitude")]
coords <- temp1[,c("lon","lat")]
points <- SpatialPoints(coords,proj4string=clim@crs)
climvals <- extract(clim,points)
temp1 <- cbind(temp1,climvals)
temp1$MAT <- (temp1$MAT/10)

#gridmat <- matrix(as.vector(temp1$MAT), nrow=length(lon.list), ncol=length(lat.list))
gridmat <- matrix(as.vector(temp1$MAT), nrow=length(unique(temp1$lon)), length(unique(temp1$lat)))
geom_point(aes(x = lon, y = lat, colour = MAT),
           data = temp1)
usa <- map_data("usa")
world <- map_data("world", xlim = c(-140, -65), ylim = c(23, 63))

ggplot() + 
  geom_polygon(data = usa, aes(x = long, y = lat, group = group), 
               alpha = 0.01, 
               colour = 'black', 
               fill = NA) + 
  coord_fixed(1.3) +
  geom_tile(aes(x = temp1$lon, y = temp1$lat, fill = temp1$MAT)) +
  geom_jitter() +
  labs(x = " ", y = " ", title = "MAT") +
  labs(fill = 'MAT') +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 30)) +
  guides(colour = guide_legend(override.aes = list(size=10))) 
  scale_colour_brewer(palette = "Paired")  
  
  ggplot() + 
    geom_polygon(data = usa, aes(x = long, y = lat, group = group), 
                 alpha = 0.01, 
                 colour = 'black', 
                 fill = NA) + 
    coord_fixed(1.3) +
    geom_tile(aes(x = temp1$lon, y = temp1$lat, fill = temp1$MAP)) +
    geom_jitter() +
    labs(x = " ", y = " ", title = "MAP") +
    labs(fill = 'MAT') +
    theme(panel.background = element_rect(fill = 'white', colour = 'white'),
          plot.title = element_text(hjust = 0.5),
          text = element_text(size = 30)) +
    guides(colour = guide_legend(override.aes = list(size=10))) 
  scale_colour_brewer(palette = "Paired")  

# nlev = 64
# plot <- image.plot(lon.list, lat.list, gridmat, 
#                    nlevel = nlev,
#                    xlab = "longitude",
#                    ylab = "latitude",
#                    ylim = c(23, 53), xlim = c(-130,-65),
#                    cex.axis=2,cex.lab=2,legend.cex=1.8)
# title(main = "MAT", cex.main=2)
# maps::map('world', regions='usa', add=TRUE)
# 
# gridmap <- matrix(as.vector(temp1$MAP), nrow=length(lon.list), ncol=length(lat.list))
# nlev = 64
# plot <- image.plot(lon.list, lat.list, gridmap, 
#                    nlevel = nlev,
#                    xlab = "longitude",
#                    ylab = "latitude",
#                    ylim = c(23, 53), xlim = c(-130,-65),
#                    cex.axis=2,cex.lab=2,legend.cex=1.8)
# title(main = "MAP", cex.main=2)
# maps::map('world', regions='usa', add=TRUE)
  
plot(temp1$values, temp1$MAP)
plot(temp1$values, temp1$MAT)
plot(temp1$values, temp1$lat)
plot(temp1$values, temp1$lon)




# mantels test for spatial autocorrelation
# https://stats.idre.ucla.edu/r/faq/how-can-i-perform-a-mantel-test-in-r/

lonvec <- rep(lon.list, 78) #78 is length of lat list
latvec <- rep(lat.list, 319) #319 is lenght of lon list
library(ape)

dat <- map_bain_plot$pct.grid #rotate through map_gs_plot, map_rs_plot, map_ms_plot, map_bain_plot
datv <- as.vector(dat)

dat <- cbind.data.frame(datv, latvec, lonvec)
colnames(dat) <- c("pctch", "lat", "lon")

dat1 <- dat[complete.cases(dat),]

plot.dists <- as.matrix(dist(cbind(dat1$lon, dat1$lat)))

plot.dists.inv <- 1/plot.dists
diag(plot.dists.inv) <- 0

Moran.I(dat1$pctch, plot.dists.inv) #moran's I test (right)


# checking for simpsons paradox ------------------------------
require(Simpsons)
FIA_plot$age <- as.numeric(as.character(FIA_plot$YOU0OLD1))
FIA_plot$ndep <- as.numeric(as.character(FIA_plot$NDEPg))
FIA_plot$sm <- as.numeric(as.character(FIA_plot$SMg))
FIA_plot$cun <- as.numeric(as.character(FIA_plot$CUNg))
FIA_plot$cso <- as.numeric(as.character(FIA_plot$CSOg))
example1 = Simpsons(BA_Nfixer_prop_total, BA_change, 
                    clusterid=cso, 
                    data=FIA_plot)
coef(example1)

# plot for simpsons paradox
number_bins <- 50

all_bachange <- stats.bin(x=FIA_plot$BA_Nfixer_prop_total,
                          y=FIA_plot[,"BA_change"],
                          N=number_bins)
young_bachange <- stats.bin(x=FIA_plot$BA_Nfixer_prop_total,
                            y=FIA_plot[FIA_plot$YOU0OLD1==0,"BA_change"],
                            N=number_bins)
old_bachange <- stats.bin(x=FIA_plot$BA_Nfixer_prop_total,
                            y=FIA_plot[FIA_plot$YOU0OLD1==1,"BA_change"],
                            N=number_bins)

zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
plot(seq(1,number_bins,1),all_bachange$stats["mean",],xaxt="n",
     xlab="BA from N-Fixers (%)",
     ylab=expression(paste("BAI (m2 yr"^"-1",")")),
     ylim=c(-6,6),cex.lab=1.5,col="black")
mtext("c", 2, adj=5,las=1,padj=-12,cex=2)
points(seq(1,number_bins,1),young_bachange$stats["mean",],col="red")
points(seq(1,number_bins,1),old_bachange$stats["mean",],col="blue")

par(new=TRUE)
model_young <- lm(FIA_plot[FIA_plot$YOU0OLD1==0,"BA_change"]~FIA_plot[FIA_plot$YOU0OLD1==0,"BA_Nfixer_prop_total"])
model_old <- lm(FIA_plot[FIA_plot$YOU0OLD1==1,"BA_change"]~FIA_plot[FIA_plot$YOU0OLD1==1,"BA_Nfixer_prop_total"])
plot(seq(0,1,0.1),model_young$coefficients[1]+model_young$coefficients[2]*seq(0,1,0.1),
     type="l",ylim=c(-6,6),
     xaxt="n",yaxt="n",ylab="",xlab="",col="red")
lines(seq(0,1,0.1),model_old$coefficients[1]+model_old$coefficients[2]*seq(0,1,0.1),
      ylim=c(-6,6),
      xaxt="n",yaxt="n",ylab="",xlab="",col="blue")
axis(side=1,at=seq(0,1,0.1),labels=seq(0,100,10))
legend("topleft",legend=c("young","old"),col=c("red","blue"),bty="n",pch=c(1,1))

#lines(density(FIA_plot[FIA_plot$YOU0OLD1==0, "BA_Nfixer_prop_total"]), col="red")
#lines(density(FIA_plot[FIA_plot$YOU0OLD1==1, "BA_Nfixer_prop_total"]), col="blue")

par(mar=c(0,3,1,1))
hy <- hist(FIA_plot[FIA_plot$YOU0OLD1==0, "BA_Nfixer_prop_total"], breaks = 49, 
     prob=TRUE, plot=FALSE)
ho <- hist(FIA_plot[FIA_plot$YOU0OLD1==1, "BA_Nfixer_prop_total"], breaks = 49, 
     prob=TRUE, plot=FALSE)

top = max(c(hy$counts, ho$counts))
barplot(hy$counts, axes=FALSE, ylim=c(0, top), border = 'red', space=0)
barplot(ho$counts, axes=FALSE, ylim=c(0, top), border = 'blue', space=0, add=T)


# N fixer abundance hist ------------------------------------
# Figure SI 1
temp <- FIA_plot[FIA_plot$BA_Nfixer1>0,]
temp$BA_pct <- temp$BA_Nfixer_prop_total*100

png("output/figS1.png", width = 12, height = 7, units = 'in', res = 300)
par(bty = 'n')
hist(temp$BA_pct, breaks=30, prob=T,
     xlab="Percent basal area comprised by N-fixing trees", main = "",
     cex.lab = 1.5, cex.axis = 1.5)
dev.off()

