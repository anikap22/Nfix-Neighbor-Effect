
############################
#
# Bayesian modeling
#
# Anika Petach
# 7/8/18
#
############################
rm(list=ls())

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

############################# GROWTH ########################
# Growth start here ------------------------------------------------
# Load data
gs <- readRDS("output/nci_formodels_growth.RDS") #made in growth_prep.R

spdata <- read.csv("/Users/Anika/Documents/GradSchool/FIA_EstimateProj/ACGCA_FIAv51_from_JWL/REF_SPECIES_jwl.csv")
spdata <- spdata %>% dplyr::select(SPCD, COMMON_NAME, GENUS, SPECIES, JENKINS_SPGRPCD) #get biomass parameters if needed

gs <- merge(gs, spdata, by.x="SPCD", by.y="SPCD", all.x=T, all.y=F)

source("scripts/functions.R")

# Define evergreen vs deciduous
classten <- read.csv("data/raw/Jenkins10_decid0ever1.csv")
classten <- classten[, c("SPCD", "DEC0EVER1")]

gs <- merge(gs, classten, by.x="SPCD", by.y="SPCD", all.x=T)
gs$DEC0EVER1 <- NA
gs[gs$JENKINS_SPGRPCD <= 5, "DEC0EVER1"] <- 1
gs[gs$JENKINS_SPGRPCD >= 6 & gs$JENKINS_SPGRPCD < 10, "DEC0EVER1"] <- 0


# Check data distribution -----------------------------------------
qqp(gs$LGRp, "norm") #normal?
#ok -2 to 2, really skewed above 2
qqp(gs$LGRp, "lnorm") #lognormal?
nbinom <- fitdistr(gs$LGRp, "Negative Binomial")
qqp(gs$LGRp, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]]) #nbin?
poisson <- fitdistr(gs$LGRp, "Poisson")
qqp(gs$LGRp, "pois", poisson$estimate) #poisson?
gamma <- fitdistr(gs$LGRp, "gamma")
qqp(gs$LGRp, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]]) #gamma?



# Growth model development ---------------------------------------
## LMER model
model.int = lmer(LGRp ~ NCI_props + 
                   NCIs +
                   dbhs + 
                   NCIs^2 + 
                   NCI_props * NCIs + 
                   dbhs * NCIs +
                   (1 | PCN1) +
                   (1 | SPCD), 
                 data = gs)
summary(model.int)
MuMIn::r.squaredGLMM(model.int)
plot_model(model.int, vline = "grey")

model.int1 = lmer(LGRp ~ NCI_props +
                    NCIs +
                    dbhs +
                    NCIs^2 +
                    NCI_props * NCIs +
                    dbhs * NCIs +
                    (1 | PCN1), 
                  data = gs)
plot_model(model.int1)
MuMIn::r.squaredGLMM(model.int1)
ggCaterpillar(ranef(model.int1, condVar = TRUE), QQ = FALSE) #see function below
dotplot(ranef(model.int1, condVar = TRUE))


model.int2 = lmer(LGRp ~ NCI_props + 
                    NCIs +
                    dbhs +
                    NCIs^2 +
                    NCI_props * NCIs +
                    dbhs * NCIs +
                    (1 | SPCD), 
                  data = gs)
plot_model(model.int2)
ggCaterpillar(ranef(model.int2, condVar = TRUE), QQ = FALSE)

library(dotwhisker)
dwplot(list(plotsp = model.int,
            plot = model.int1,
            sp = model.int2))+
  theme_bw()+
  geom_vline(xintercept = 0,lty = 2)

model.int3 = lmer(LGRp ~ NCI_props + 
                    NCIs +
                    dbhs +
                    NCIs^2 +
                    NCI_props * NCIs +
                    dbhs * NCIs +
                    (1 | Genus), 
                  data = gs)

model.int4 = lmer(LGRp ~ NCI_props + 
                    NCIs +
                    dbhs +
                    NCIs^2 +
                    NCI_props * NCIs +
                    dbhs * NCIs +
                    (1 | Genus) +
                    (1|PCN1), 
                  data = gs)

model.int5 = lmer(LGRp ~ NCI_props + 
                    NCIs +
                    dbhs +
                    (1 | Genus) +
                    (1|PCN1), 
                  data = gs)

anova(model.int, model.int1, model.int2, model.int3, model.int4)

aic(model.int, model.int1, model.int2)

# Coefficient plot with different models developed
sjt.lmer(model.int, model.int1, model.int2,
         p.kr = F,
         show.aic = T,
         show.r2 = T,
         show.header = T,
         string.est = "Estimate", digits.est = 5,
         string.ci = "Conf. Int.", digits.ci = 5,
         string.p = "p-value", digits.p = 5,
         digits.se = 5, digits.summary = 5, digits.std = 5,
         string.dv = "Response",
         string.pred = "Coefficients",
         depvar.labels = c("RE=PCN&SPCD", "RE=PCN", "RE=SPCD"))

anova(model.int, model.plot.nt)


# Growth models with random slope: ------------------
model1 <- lmer(LGRp ~ NCI_props + 
                 NCIs + 
                 dbhs + 
                 (1 + NCI_prop | PCN1) +
                 (1 + NCI_prop | Genus), 
               data = f, REML = F)
summary(model1)
MuMIn::r.squaredGLMM(model1) #gives marginal and conditional R2
# marginal R2 is proportion of variance explained by fixed effects
# conditional R2 is proportion of variance explained by random effects
plot_model(model1)

ggCaterpillar(ranef(model1, condVar=TRUE), QQ=FALSE)

levels(f$Genus) <- c("Acacia", "Albizia", "Alnus", "Casuarina", "Cercocarpus",
                     "Elaeagnus", "Olneya", "Piscidia", "Prosopis", "Robinia", 
                     "Sophora", "NA")
droplevels.data.frame(f$Genus)

boxplot(LGRp ~ Genus, data=f)

model2 <- lmer(LGRp ~ NCI_props + 
                 NCIs + 
                 dbhs + 
                 (1 | PCN1) +
                 (1 + NCI_prop | GENUS), 
               data = nf, REML = F)
summary(model2)
MuMIn::r.squaredGLMM(model2) #gives marginal and conditional R2
plot_model(model2)

source("scripts/ggCaterpillar.R")
ggCaterpillar(ranef(model2, condVar=TRUE), QQ=FALSE)

ranef(model2)$SPCD
pp <- profile(model2)
confint(pp)

model3 <- lmer(LGRp ~ NCI_props + 
                 NCIs + 
                 dbhs + 
                 (1 | PCN1) +
                 (1 + NCI_prop | STDAGE), 
               data = gs, REML = F)
summary(model3)
MuMIn::r.squaredGLMM(model3) #gives marginal and conditional R2
plot_model(model3)
ggCaterpillar(ranef(model3, condVar=TRUE), QQ=FALSE)

model4 <- lmer(LGRp ~ NCI_props + 
                 NCIs +
                 dbhs +
                 NCIs^2 +
                 NCI_props * NCIs +
                 dbhs * NCIs +
                 (1 | Genus) + (0 + NCI_pros | Genus)
                 (1 | PCN1), 
               data = gs)


# Growth models on subsets of data ---------------------------------

# Subset data
young <- gs %>% dplyr::filter(STDAGE < 60)
old <- gs %>% dplyr::filter(STDAGE >= 60)

f <- gs %>% dplyr::filter(FIX == 1) #fixers
nf <- gs %>% dplyr::filter(FIX == 0) #nonfixers
canopy <- gs %>% dplyr::filter(CCLCD1 < 4)
noncanopy <- gs %>% dplyr::filter(CCLCD1 >= 4)

evergreen <- gs %>% dplyr::filter(DEC0EVER1 == 1)
deciduous <- gs %>% dplyr::filter(DEC0EVER1 == 0)

# Plot data density
plot(density(young$NCI_prop))
plot(density(old$NCI_prop))
plot(density(f$NCI_prop))
plot(density(nf$NCI_prop))
plot(density(canopy$NCI_prop))
plot(density(noncanopy$NCI_prop))
plot(density(evergreen$NCI_prop))
plot(density(deciduous$NCI_prop), col='red')

# Model with all data
model.plot.nt = lmer(LGRp ~ NCI_props + 
                       NCIs + 
                       dbhs + 
                       NCIs^2 +
                       NCI_props * NCIs +
                       dbhs * NCIs +
                       (1 | PCN1), 
                     data = gs, REML = F)
summary(model.plot.nt)
MuMIn::r.squaredGLMM(model.plot.nt) #gives marginal and conditional R2
plot_model(model.plot.nt)

# Model with nonfixers
model.plot.nf = lmer(LGRp ~ NCI_prop +
                       NCI + 
                       dbh + 
                       NCIs^2 +
                       NCI_props * NCIs +
                       dbhs * NCIs +
                       (1 | PCN1), 
                     data = nf, REML = F)
summary(model.plot.nf)
MuMIn::r.squaredGLMM(model.plot.nf)

# Model with fixers
model.plot.f = lmer(LGRp ~ NCI_prop + 
                      NCI + 
                      dbh + 
                      NCIs^2 +
                      NCI_props * NCIs +
                      dbhs * NCIs +
                      (1 | PCN1), 
                    data = f, REML = F)
summary(model.plot.f)
MuMIn::r.squaredGLMM(model.plot.f)

sjt.lmer(model.plot.nt, model.plot.nf, model.plot.f,
         p.kr = F,show.aic = T)

sjt.lmer(model.plot.nt, model.plot.nf, model.plot.f,
         p.kr = F,
         show.aic = T,
         show.r2 = T,
         show.header = T,
         string.est = "Estimate", digits.est = 5,
         string.ci ="Conf. Int.", digits.ci = 5,
         string.p = "p-value", digits.p = 5,
         digits.se = 5, digits.summary = 5,digits.std = 5,
         string.dv = "Response",
         string.pred = "Coefficients",
         depvar.labels = c("Forest GR", "Non Fixer GR", "Fixer GR"))

stargazer::stargazer(model.plot.nt, 
                     type = "text",
                     digits = 5, 
                     star.cutoffs = c(0.05, 0.01, 0.001),
                     digit.separator = "")

qqnorm(resid(model.plot.nt))
qqline(resid(model.plot.nt))
qqnorm(resid(model.plot.nf))
qqline(resid(model.plot.nf))
qqnorm(resid(model.plot.f))
qqline(resid(model.plot.f))

model.all = lmer(LGRp ~ NCI_props + 
                   NCIs + 
                   NCIs^2 + 
                   dbhs + 
                   NCI_props * NCIs^2 + 
                   dbhs * NCIs + 
                   NCI_props * NCIs + 
                   (1 | PCN1) + 
                   (1 + NCI_props | Genus), 
                 data = gs, REML = F)
summary(model.all)
MuMIn::r.squaredGLMM(model.all) #gives marginal and conditional R2
plot_model(model.all, vline.color = "red")

source("scripts/ggCaterpillar.R")
ggCaterpillar(ranef(model.all, condVar=TRUE), QQ=FALSE)

drop1(model.all)

# Growth model what drives +/- effect -------------------------------
model.plot.f = lmer(LGRp ~ NCI_props + 
                      NCIs + 
                      dbhs + 
                      (NCI_props | LAT) +
                      (NCI_props | LON) +
                      (NCI_props | ELEVm) +
                      (1 | PCN1), 
                    data = gs, REML = F)

gs$latcat <- as.factor(trunc(gs$LAT))
gs$loncat <- as.factor(trunc(gs$LON))
gs$elevcat <- as.factor(trunc(gs$ELEVm))

m.lat = lmer(LGRp ~ NCI_props + 
                      NCIs + 
                      dbhs +
                      NCIs^2 +
                      NCI_props * NCIs +
                      dbhs * NCIs +
                      (NCI_props | latcat) +
                      (1 | PCN1), 
                    data = gs, REML = F)
ggCaterpillar(ranef(m.lat, condVar=TRUE), QQ=FALSE)

m.lon = lmer(LGRp ~ NCI_props + 
               NCIs + 
               dbhs +
               NCIs^2 +
               NCI_props * NCIs +
               dbhs * NCIs +
               (NCI_props | loncat) +
               (1 | PCN1), 
             data = gs, REML = F)
ggCaterpillar(ranef(m.lon, condVar=TRUE), QQ=FALSE)

m.elev = lmer(LGRp ~ NCI_props + 
               NCIs + 
               dbhs +
               NCIs^2 +
               NCI_props * NCIs +
               dbhs * NCIs +
               (NCI_props | elevcat) +
               (1 | PCN1), 
             data = gs, REML = F)
ggCaterpillar(ranef(m.elev, condVar=TRUE), QQ=FALSE)


# Growth model eval ------------------------------------------
plot(predict(model.int), 
     gs$LGRp,
     xlab = "predicted",
     ylab = "actual")
abline(a = 0, b = 1, col = 'red')
plot(model.int)

plot(predict(model.ints), nt$LGRp)
plot(model.ints)

#############################
## Run the model.int on the subsets of data

model.int = lmer(LGRp ~ NCI_props +
                   NCIs +
                   dbhs +
                   NCIs^2 +
                   NCI_props * NCIs +
                   dbhs * NCIs +
                   (1 | PCN1) +
                   (1 | SPCD), 
                 data = gs)
summary(model.int)
MuMIn::r.squaredGLMM(model.int)
plot_model(model.int, vline = "grey")

ma <- lmer(LGRp ~ NCI_props + NCIs + dbhs + NCIs^2 + NCI_props*NCIs + dbhs*NCIs + 
             (1|PCN1.x) + (1|SPCD2), data = gs[gs$FIX==0,])
mf <- lmer(LGRp ~ NCI_props + NCIs + dbhs + NCIs^2 + NCI_props*NCIs + dbhs*NCIs +
             (1|PCN1.x) + (1|SPCD2), data = f)
mnf <- lmer(LGRp ~ NCI_props + NCIs + dbhs + NCIs^2 + NCI_props*NCIs + dbhs*NCIs +
              (1|PCN1.x) + (1|SPCD2), data = nf)
mc <- lmer(LGRp ~ NCI_props + NCIs + dbhs + NCIs^2 + NCI_props*NCIs + dbhs*NCIs +
             (1|PCN1.x) + (1|SPCD2), data = canopy)
mnc <- lmer(LGRp ~ NCI_props + NCIs + dbhs + NCIs^2 + NCI_props*NCIs + dbhs*NCIs +
              (1|PCN1.x) + (1|SPCD2), data = noncanopy)
my <- lmer(LGRp ~ NCI_props + NCIs + dbhs + NCIs^2 + NCI_props*NCIs + dbhs*NCIs +
             (1|PCN1.x) + (1|SPCD2), data = young[young$FIX==0,])
mo <- lmer(LGRp ~ NCI_props + NCIs + dbhs + NCIs^2 + NCI_props*NCIs + dbhs*NCIs +
             (1|PCN1.x) + (1|SPCD2), data = old[old$FIX==1,])


plot_model(ma) #customize from: https://www.r-bloggers.com/one-function-to-rule-them-all-visualization-of-regression-models-in-rstats-w-sjplot/

MuMIn::r.squaredGLMM(ma)
MuMIn::r.squaredGLMM(mf)
MuMIn::r.squaredGLMM(mnf)
MuMIn::r.squaredGLMM(mc)
MuMIn::r.squaredGLMM(mnc)

dwplot(list(forest = ma, fixers = mf, nonfixers = mnf, canopy = mc, understory = mnc))+
  theme_bw() +
  geom_vline(xintercept = 0, lty = 2)
multiplot(mf, mnf, mc, mnc, my, mo,
          title = "Coefficient Plot", 
          cex = 1.5,
          names = c("Fixers", "Nonfixers", "Canopy", "Noncanopy", "Young", "Old"))

stargazer::stargazer(ma, 
                     type = "text",
                     digits = 5, star.cutoffs = c(0.05, 0.01, 0.001),
                     digit.separator = "")

sjt.lmer(mf, mnf, mc, mnc, my, mo,
         p.kr = F,
         show.aic = T,
         show.r2 = T,
         show.header = T,
         string.est = "Estimate", digits.est = 5,
         string.ci = "Conf. Int.", digits.ci = 5,
         string.p = "p-value", digits.p = 5,
         digits.se = 5, digits.summary = 5, digits.std = 5,
         string.dv = "Response",
         string.pred = "Coefficients",
         depvar.labels = c("Fixer GR", "Non Fixer GR", "Canopy GR", "Understory GR", "Young GR", "Old GR"))

# Misc
plot(gs$NCI_prop, gs$LGRp)

sjp.lmer(mf, type = "re")

# Pedro's plot -------------------------------------------
# NCItot x dbh (his was pasture x perc)
forest <- rep(seq(0, 100, 1), each = 101)
prec <- rep(seq(0, 200, 2), 101)
predT1 <- exp(fixef(ma)[1] + 
                fixef(ma)[3]*scale(forest)[,1] + 
                fixef(ma)[2]*scale(log(prec+1))[,1] +
                fixef(ma)[8]*scale(forest*log(prec+1))[,1]-1)
contourplot(predT1 ~ pasture*prec, 
            cuts = 12,
            col.regions = rainbow(100),
            region = T, 
            cex.lab = 1.5, 
            ylab = "2 Days cumulative precipitation", 
            xlab = "% Pasture in watershed", 
            main = "Turbidity")



colors <- colorRampPalette(c('red', 'darkblue'), bias=1)(256) 
colors <- colorRampPalette(c('yellow', 'dark green'))(256)

ncis <- seq(0, 1, 0.1) #need 5358143 values
ncis <- runif(5358143)
predT1 <- (fixef(model_non_gs)[1] + fixef(model_non_gs)[2]*scale(ncis)[,1]
         +fixef(model_non_gs)[3]*scale(gs[gs$FIX == 0,"NCIs"])[,1]
         +fixef(model_non_gs)[4]*scale(gs[gs$FIX == 0,"dbhs"])[,1]
         +fixef(model_non_gs)[5]*scale(gs[gs$FIX == 0,"NCIs"])[,1]*scale(ncis)[,1])

predover <- data.frame(NCI_props = seq(0, 1, 0.1),
                       NCIs = rep(mean(gs[gs$FIX == 1, "NCIs"]), 11),
                       dbhs = rep(mean(gs[gs$FIX == 1, "dbhs"]), 11))
pred <- predict(model_N_gs, 
                re.form = NA, 
                newdata = predover)
predT1 <- (fixef(model_N_gs)[1] + fixef(model_N_gs)[2]*scale(ncis)[,1])

# Growth: Plot binned data for individual trees -------------------------------------



source("scripts/functions.R")

# bin data for plotting
number_bins <- 50
stats_all_NCI_prop_young <- stats.bin(x = young$NCI_prop,
                                    y = young[,"LGRp"],
                                    N = number_bins)
stats_N_NCI_prop_young <- stats.bin(x = young[young$FIX == 1,"NCI_prop"],
                                  y = young[young$FIX == 1,"LGRp"],
                                  N = number_bins)
stats_non_NCI_prop_young <- stats.bin(x = young[young$FIX == 0,"NCI_prop"],
                                    y = young[young$FIX == 0,"LGRp"],
                                    N = number_bins)

#Young:
png("output/young_growth.png", 
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
     ylab = expression(paste("Relative Growth Rate (yr"^"-1",")")),
     main = "Young Trees (Stand Age < 60 Yrs)",
     ylim = c(0, 0.05),
     cex.lab = 1.5,
     col = "blue")
#mtext("c", 2, adj=5,las=1,padj=-12,cex=2)
points(seq(1, number_bins, 1),
       stats_N_NCI_prop_young$stats["mean", ],
       col = "red")
par(new = TRUE)
# model_non_young<-lm(young[young$FIX==0,"LGRp"]~young[young$FIX==0,"NCI_prop"])
# model_N_young<-lm(young[young$FIX==1,"LGRp"]~young[young$FIX==1,"NCI_prop"])
# plot(seq(0,1,0.1),model_non_young$coefficients[1]+model_non_young$coefficients[2]*seq(0,1,0.1),
#      type="l",ylim=c(0,0.05),
#      xaxt="n",yaxt="n",ylab="",xlab="",col="blue")
# lines(seq(0,1,0.1),model_N_young$coefficients[1]+model_N_young$coefficients[2]*seq(0,1,0.1),
#       ylim=c(0,0.05),
#       xaxt="n",yaxt="n",ylab="",xlab="",col="red")
model_non_young <- lmer(LGRp ~ NCI_props +
                          NCIs +
                          dbhs +
                          NCIs^2 +
                          NCI_props * NCIs +
                          dbhs * NCIs +
                          (1 | PCN1) + (1 | GENUS), 
                        data = young[young$FIX == 0, ])
model_N_young <- lmer(LGRp ~ NCI_props +
                        NCIs +
                        dbhs +
                        NCIs^2 +
                        NCI_props * NCIs +
                        dbhs * NCIs +
                        (1 | PCN1) + (1 | GENUS), 
                      data = young[young$FIX == 1, ])
#predf<-function(x){fixef(model_non_young)[1]+colMeans(ranef(model_non_young)$SPCD2)+colMeans(ranef(model_non_young)$PCN1.x)+fixef(model_non_young)[2]*(x-mean(x))/sd(x)}
#predT2<-predf(ncis)
plot(seq(0, 1, 0.1),
     stats_non_NCI_prop_young$stats["mean", ][[1]] +
       fixef(model_non_young)[2]*seq(0, 1, 0.1),
     type = "l",
     ylim = c(0, 0.05),
     xaxt = "n",
     yaxt = "n",
     ylab = "",
     xlab = "",
     col = "blue")
#predf<-function(x){fixef(model_N_young)[1]+colMeans(ranef(model_N_young)$SPCD2)+colMeans(ranef(model_N_young)$PCN1.x)+fixef(model_N_young)[2]*(x-mean(x))/sd(x)}
#predT2<-predf(ncis)
lines(seq(0, 1, 0.1),
      stats_N_NCI_prop_young$stats["mean", ][[1]] +
        fixef(model_N_young)[2]*seq(0, 1, 0.1),
      ylim = c(0, 0.05),
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

f_pctch(model_non_young)
# orig <- stats_non_NCI_prop_young$stats["mean", ][[1]] + fixef(model_non_young)[2]*0
# final <- stats_non_NCI_prop_young$stats["mean", ][[1]] + fixef(model_non_young)[2]*1
# ch <- 100 * (final - orig) / orig
# ch

f_pctch(model_N_young)
# orig <- stats_N_NCI_prop_young$stats["mean", ][[1]] + fixef(model_N_young)[2]*0
# final <- stats_N_NCI_prop_young$stats["mean", ][[1]] + fixef(model_N_young)[2]*1
# ch <- 100 * (final - orig) / orig
# ch

#Old:
png("output/old_growth.png", 
    width = 8, 
    height = 4, 
    units = 'in', 
    res = 300)
par(mar = c(5.1,5.1,4.1,2.1)+0.1,
    cex.main = 1.5,
    cex.lab = 1.5,
    cex.axis = 1.3,
    bty = "n")
stats_all_NCI_prop_old <- stats.bin(x = old$NCI_prop,
                                  y = old[ ,"LGRp"],
                                  N = number_bins)
stats_N_NCI_prop_old <- stats.bin(x = old[old$FIX == 1,"NCI_prop"],
                                y = old[old$FIX == 1,"LGRp"],
                                N = number_bins)
stats_non_NCI_prop_old <- stats.bin(x = old[old$FIX == 0,"NCI_prop"],
                                  y = old[old$FIX == 0,"LGRp"],
                                  N = number_bins)
plot(seq(1, number_bins, 1),
     stats_non_NCI_prop_old$stats["mean", ],
     xaxt = "n",
     xlab = "Crowding from N-Fixers (%)",
     ylab = expression(paste("Relative Growth Rate (yr"^"-1",")")),
     main = "Old Trees (Stand Age >= 60 Yrs)",
     ylim = c(0, 0.05),
     cex.lab = 1.5,
     col = "blue")
#mtext("d", 2, adj=5,las=1,padj=-12,cex=2)
points(seq(1, number_bins, 1),
       stats_N_NCI_prop_old$stats["mean", ],
       col = "red")
par(new = TRUE)
# model_non_old<-lm(old[old$FIX==0,"LGRp"]~old[old$FIX==0,"NCI_prop"])
# model_N_old<-lm(old[old$FIX==1,"LGRp"]~old[old$FIX==1,"NCI_prop"])
# plot(seq(0,1,0.1),model_non_old$coefficients[1]+model_non_old$coefficients[2]*seq(0,1,0.1),
#      type="l",ylim=c(0,0.05),
#      xaxt="n",yaxt="n",ylab="",xlab="",col="blue")
# lines(seq(0,1,0.1),model_N_old$coefficients[1]+model_N_old$coefficients[2]*seq(0,1,0.1),
#       ylim=c(0,0.05),
#       xaxt="n",yaxt="n",ylab="",xlab="",col="red")
oldnf <- old[old$FIX == 0, ]
model_non_old <- lmer(LGRp ~ NCI_props +
                        NCIs +
                        dbhs +
                        NCIs^2 +
                        NCI_props * NCIs +
                        dbhs * NCIs +
                        (1 | PCN1) + (1 | GENUS), 
                      data = oldnf)
oldf <- old[old$FIX == 1, ]
model_N_old <- lmer(LGRp ~ NCI_props +
                      NCIs +
                      dbhs +
                      NCIs^2 +
                      NCI_props * NCIs +
                      dbhs * NCIs +
                      (1 | PCN1) + (1 | GENUS), 
                    data = oldf)
#predf<-function(x){fixef(model_non_old)[1]+colMeans(ranef(model_non_old)$SPCD2)+colMeans(ranef(model_non_old)$PCN1.x)+fixef(model_non_old)[2]*(x-mean(x))/sd(x)}
#predT2<-predf(ncis)
plot(seq(0, 1, 0.1),
     stats_non_NCI_prop_old$stats["mean", ][[1]] + fixef(model_non_old)[2]*seq(0, 1, 0.1),
     type = "l",
     ylim = c(0, 0.05),
     xaxt = "n",
     yaxt = "n",
     ylab = "",
     xlab = "",
     col = "blue")
#predf<-function(x){fixef(model_N_old)[1]+colMeans(ranef(model_N_old)$SPCD2)+colMeans(ranef(model_N_old)$PCN1.x)+fixef(model_N_old)[2]*(x-mean(x))/sd(x)+fixef(model_N_old)[3]*mean(scale(oldf$NCI))+fixef(model_N_old)[4]*mean(scale(oldf$dbh))+fixef(model_N_old)[5]*x*mean(scale(oldf$NCI))+fixef(model_N_old)[6]*mean(scale(oldf$NCI))*mean(scale(oldf$dbh))}
#predT2<-predf(ncis)
lines(seq(0, 1, 0.1),
      stats_N_NCI_prop_old$stats["mean", ][[1]] + fixef(model_N_old)[2]*seq(0, 1, 0.1),
      ylim = c(0, 0.05),
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

f_pctch(model_non_old)
# orig <- stats_non_NCI_prop_old$stats["mean", ][[1]] + fixef(model_non_old)[2]*0
# final <- stats_non_NCI_prop_old$stats["mean", ][[1]] + fixef(model_non_old)[2]*1
# ch <- 100 * (final - orig) / orig
# ch

f_pctch(model_N_old)
# orig <- stats_N_NCI_prop_old$stats["mean", ][[1]] + fixef(model_N_old)[2]*0
# final <- stats_N_NCI_prop_old$stats["mean", ][[1]] + fixef(model_N_old)[2]*1
# ch <- 100 * (final - orig) / orig
# ch


#All:
number_bins <- 50
stats_all_NCI_prop_gs <- stats.bin(x = gs$NCI_prop,
                                   y = gs[,"LGRp"],
                                   N = number_bins)
stats_N_NCI_prop_gs <- stats.bin(x = gs[gs$FIX == 1, "NCI_prop"],
                                 y = gs[gs$FIX == 1, "LGRp"],
                                 N = number_bins)
stats_non_NCI_prop_gs <- stats.bin(x = gs[gs$FIX == 0, "NCI_prop"],
                                   y = gs[gs$FIX == 0, "LGRp"],
                                   N = number_bins)

#pdf(file="OR_gr.pdf",width=14)
png("output/allforest_growth.png", 
    width = 8, 
    height = 4, 
    units = 'in', 
    res = 300)
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
     ylab = expression(paste("Relative Growth Rate (yr"^"-1",")")),
     main = "Forest",
     bty = "n",
     pch = 21,
     ylim = c(0, 0.05),
     cex = 1.5,
     col = "blue")
#mtext("c", 2, adj=5,las=1,padj=-12,cex=2)
points(seq(1, number_bins, 1),
       stats_N_NCI_prop_gs$stats["mean", ],
       col = "red",
       pch = 21,
       cex = 1.5)
par(new = TRUE)
#model_non_gs<-lm(gs[gs$FIX==0,"LGRp"]~gs[gs$FIX==0,"NCI_prop"])
model_non_gs <- lmer(LGRp ~ NCI_props +
                       NCIs +
                       dbhs +
                       NCIs^2 +
                       NCI_props * NCIs +
                       dbhs * NCIs +
                       (1 | PCN1) + (1 | GENUS), 
                     data = gs[gs$FIX == 0, ])
#model_N_gs<-lm(gs[gs$FIX==1,"LGRp"]~gs[gs$FIX==1,"NCI_prop"])
model_N_gs <- lmer(LGRp ~ NCI_props +
                     NCIs +
                     dbhs +
                     NCIs^2 +
                     NCI_props * NCIs +
                     dbhs * NCIs +
                     (1 | PCN1) + (1 | GENUS), 
                   data = gs[gs$FIX == 1, ])
#predf<-function(x){fixef(model_non_gs)[1]+colMeans(ranef(model_non_gs)$SPCD2)+colMeans(ranef(model_non_gs)$PCN1.x)+fixef(model_non_gs)[2]*(x-mean(x))/sd(x)}
#predT2<-predf(ncis)
plot(seq(0, 1, 0.1),
     stats_non_NCI_prop_gs$stats["mean", ][[1]] + fixef(model_non_gs)[2]*seq(0, 1, 0.1),
     type = "l",
     ylim = c(0, 0.05),
     xaxt = "n",
     yaxt = "n",
     ylab = "",
     xlab = "",
     col = "blue",
     lwd = 2)
#predf<-function(x){fixef(model_N_gs)[1]+colMeans(ranef(model_N_gs)$SPCD2)+colMeans(ranef(model_N_gs)$PCN1.x)+fixef(model_N_gs)[2]*(x-mean(x))/sd(x)}
#predT2<-predf(ncis)
lines(seq(0, 1, 0.1),
      stats_N_NCI_prop_gs$stats["mean", ][[1]] + fixef(model_N_gs)[2]*seq(0, 1, 0.1),
      ylim = c(0, 0.05),
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
       pch = c(1, 1),
       cex = c(1.5, 1.5))
dev.off()

f_pctch(model_non_gs)
# orig <- stats_non_NCI_prop_gs$stats["mean", ][[1]] + fixef(model_non_gs)[2]*0
# final <- stats_non_NCI_prop_gs$stats["mean", ][[1]] + fixef(model_non_gs)[2]*1
# ch <- 100 * (final - orig) / orig
# ch

f_pctch(model_N_gs)
# orig<-stats_N_NCI_prop_gs$stats["mean",][[1]]+fixef(model_N_gs)[2]*0
# final<-stats_N_NCI_prop_gs$stats["mean",][[1]]+fixef(model_N_gs)[2]*1
# ch<-100*(final-orig)/orig
# ch

#Canopy:
stats_all_NCI_prop_canopy <- stats.bin(x = canopy$NCI_prop,
                                       y = canopy[,"LGRp"],
                                       N = number_bins)
stats_N_NCI_prop_canopy <- stats.bin(x = canopy[canopy$FIX==1,"NCI_prop"],
                                     y = canopy[canopy$FIX==1,"LGRp"],
                                     N = number_bins)
stats_non_NCI_prop_canopy <- stats.bin(x = canopy[canopy$FIX==0,"NCI_prop"],
                                       y = canopy[canopy$FIX==0,"LGRp"],
                                       N = number_bins)

png("output/canopy_growth.png", 
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
     ylab = expression(paste("Relative Growth Rate (yr"^"-1",")")),
     main = "Canopy Trees",
     ylim = c(0,0.05),
     cex.lab = 1.5, col = "blue")

points(seq(1,number_bins,1), 
       stats_N_NCI_prop_canopy$stats["mean",], 
       col = "red")
par(new = TRUE)

model_non_canopy <- lmer(LGRp ~ NCI_props + 
                           NCIs + 
                           dbhs + 
                           NCIs^2 + 
                           NCI_props*NCIs + 
                           dbhs*NCIs +
                           (1|PCN1) + (1|GENUS), 
                         data = canopy[canopy$FIX==0,])

model_N_canopy <- lmer(LGRp ~ NCI_props + 
                         NCIs + 
                         dbhs + 
                         NCIs^2 + 
                         NCI_props*NCIs + 
                         dbhs*NCIs +
                         (1|PCN1) + (1|GENUS), 
                       data = canopy[canopy$FIX==1,])

plot(seq(0,1,0.1),stats_non_NCI_prop_canopy$stats["mean",][[1]]+fixef(model_non_canopy)[2]*seq(0,1,0.1),
     type = "l",
     ylim = c(0, 0.05),
     xaxt = "n", yaxt = "n",
     ylab = "", xlab = "",
     col = "blue")
lines(seq(0,1,0.1),stats_N_NCI_prop_canopy$stats["mean",][[1]]+fixef(model_N_canopy)[2]*seq(0,1,0.1),
      ylim = c(0, 0.05),
      xaxt = "n", yaxt = "n",
      ylab = "", xlab = "",
      col = "red")
axis(side = 1, 
     at = seq(0, 1, 0.1), 
     labels = seq(0, 100, 10))
legend("topleft", 
       legend = c("N-Fixer","Non-Fixer"), 
       col = c("red","blue"), 
       bty = "n",
       pch = c(1, 1))
dev.off()

f_pctch(model_non_canopy)
# orig <- stats_non_NCI_prop_canopy$stats["mean",][[1]]+fixef(model_non_canopy)[2]*0
# final <- stats_non_NCI_prop_canopy$stats["mean",][[1]]+fixef(model_non_canopy)[2]*1
# ch <- 100*(final-orig)/orig
# ch

f_pctch(model_N_canopy)
# orig <- stats_N_NCI_prop_canopy$stats["mean",][[1]]+fixef(model_N_canopy)[2]*0
# final <- stats_N_NCI_prop_canopy$stats["mean",][[1]]+fixef(model_N_canopy)[2]*1
# ch <- 100*(final-orig)/orig
# ch

#Non-canopy:
stats_all_NCI_prop_noncanopy <- stats.bin(x = noncanopy$NCI_prop,
                                          y = noncanopy[,"LGRp"],
                                          N = number_bins)
stats_N_NCI_prop_noncanopy <- stats.bin(x = noncanopy[noncanopy$FIX==1,"NCI_prop"],
                                        y = noncanopy[noncanopy$FIX==1,"LGRp"],
                                        N = number_bins)
stats_non_NCI_prop_noncanopy <- stats.bin(x = noncanopy[noncanopy$FIX==0,"NCI_prop"],
                                          y = noncanopy[noncanopy$FIX==0,"LGRp"],
                                          N = number_bins)

png("output/noncanopy_growth.png", 
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
     ylab = expression(paste("Relative Growth Rate (yr"^"-1",")")),
     main = "Noncanopy Trees",
     ylim = c(0,0.05),
     cex.lab = 1.5,
     col = "blue")
#mtext("c", 2, adj=5,las=1,padj=-12,cex=2)
points(seq(1, number_bins, 1), 
       stats_N_NCI_prop_noncanopy$stats["mean",], 
       col = "red")
par(new = TRUE)

model_non_noncanopy <- lmer(LGRp ~ NCI_props + 
                              NCIs + 
                              dbhs + 
                              NCIs^2 + 
                              NCI_props*NCIs +
                              dbhs*NCIs +
                              (1|PCN1) + (1|GENUS), 
                            data = noncanopy[noncanopy$FIX==0,])

model_N_noncanopy <- lmer(LGRp ~ NCI_props + 
                            NCIs + 
                            dbhs + 
                            NCIs^2 + 
                            NCI_props*NCIs +
                            dbhs*NCIs +
                            (1|PCN1) + (1|GENUS), 
                          data = noncanopy[noncanopy$FIX==1,])
# plot(seq(0,1,0.1),model_non_noncanopy$coefficients[1]+model_non_noncanopy$coefficients[2]*seq(0,1,0.1),
#      type="l",ylim=c(0,0.05),
#      xaxt="n",yaxt="n",ylab="",xlab="",col="blue")
# lines(seq(0,1,0.1),model_N_noncanopy$coefficients[1]+model_N_noncanopy$coefficients[2]*seq(0,1,0.1),
#       ylim=c(0,0.05),
#       xaxt="n",yaxt="n",ylab="",xlab="",col="red")
plot(seq(0,1,0.1),stats_non_NCI_prop_noncanopy$stats["mean",][[1]]+fixef(model_non_noncanopy)[2]*seq(0,1,0.1),
     type = "l",
     ylim = c(0,0.05),
     xaxt = "n",
     yaxt = "n",
     ylab = "",
     xlab = "",
     col = "blue")
lines(seq(0,1,0.1),(stats_N_NCI_prop_noncanopy$stats["mean",][[1]]+fixef(model_N_noncanopy)[2]*seq(0,1,0.1)),
      ylim = c(0,0.05),
      xaxt = "n",
      yaxt = "n",
      ylab = "",
      xlab = "",
      col = "red")
axis(side=1, 
     at = seq(0, 1, 0.1), 
     labels = seq(0, 100, 10))
legend("topleft", 
       legend = c("N-Fixer", "Non-Fixer"), 
       col = c("red", "blue"), 
       bty = "n",
       pch = c(1, 1))
dev.off()


f_pctch(model_non_noncanopy)

f_pctch(model_N_noncanopy)
# orig <- stats_N_NCI_prop_noncanopy$stats["mean",][[1]]+fixef(model_N_noncanopy)[2]*0
# final <- stats_N_NCI_prop_noncanopy$stats["mean",][[1]]+fixef(model_N_noncanopy)[2]*1
# ch <- 100*(final-orig)/orig
# ch

#Evergreen vs deciduous nonfixers:
# No evergreen fixers
stats_all_NCI_prop_evergreen <- stats.bin(x = evergreen$NCI_prop,
                                          y = evergreen[ ,"LGRp"],
                                          N = number_bins)
stats_non_NCI_prop_evergreen <- stats.bin(x = evergreen[evergreen$FIX == 0,"NCI_prop"],
                                          y = evergreen[evergreen$FIX == 0,"LGRp"],
                                          N = number_bins)
stats_all_NCI_prop_deciduous <- stats.bin(x = deciduous$NCI_prop,
                                          y = deciduous[ ,"LGRp"],
                                          N = number_bins)
stats_non_NCI_prop_deciduous <- stats.bin(x = deciduous[deciduous$FIX == 0,"NCI_prop"],
                                          y = deciduous[deciduous$FIX == 0,"LGRp"],
                                          N = number_bins)

png("output/evergreen_decid_growth.png", 
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

model_non_evergreen <- lmer(LGRp ~ NCI_props +
                        NCIs +
                        dbhs +
                        NCIs^2 +
                        NCI_props * NCIs +
                        dbhs * NCIs +
                        (1 | PCN1) + (1 | GENUS), 
                      data = evergreen[evergreen$FIX == 0, ])

model_non_deciduous <- lmer(LGRp ~ NCI_props +
                              NCIs +
                              dbhs +
                              NCIs^2 +
                              NCI_props * NCIs +
                              dbhs * NCIs +
                              (1 | PCN1) + (1 | GENUS), 
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


#Deciduous:
stats_all_NCI_prop_deciduous <- stats.bin(x = deciduous$NCI_prop,
                                          y = deciduous[ ,"LGRp"],
                                          N = number_bins)
stats_N_NCI_prop_deciduous <- stats.bin(x = deciduous[deciduous$FIX == 1,"NCI_prop"],
                                        y = deciduous[deciduous$FIX == 1,"LGRp"],
                                        N = number_bins)
stats_non_NCI_prop_deciduous <- stats.bin(x = deciduous[deciduous$FIX == 0,"NCI_prop"],
                                          y = deciduous[deciduous$FIX == 0,"LGRp"],
                                          N = number_bins)

png("output/deciduous_growth.png", 
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
     col = "blue")
#mtext("d", 2, adj=5,las=1,padj=-12,cex=2)
points(seq(1, number_bins, 1),
       stats_N_NCI_prop_deciduous$stats["mean", ],
       col = "red")
par(new = TRUE)

model_non_deciduous <- lmer(LGRp ~ NCI_props +
                        NCIs +
                        dbhs +
                        NCIs^2 +
                        NCI_props * NCIs +
                        dbhs * NCIs +
                        (1 | PCN1) + (1 | GENUS), 
                      data = deciduous[deciduous$FIX == 0, ])

model_N_deciduous <- lmer(LGRp ~ NCI_props +
                      NCIs +
                      dbhs +
                      NCIs^2 +
                      NCI_props * NCIs +
                      dbhs * NCIs +
                      (1 | PCN1) + (1 | GENUS), 
                    data = deciduous[deciduous$FIX == 1, ])

plot(seq(0, 1, 0.1),
     stats_non_NCI_prop_deciduous$stats["mean", ][[1]] + fixef(model_non_deciduous)[2]*seq(0, 1, 0.1),
     type = "l",
     ylim = c(0, 0.05),
     xaxt = "n",
     yaxt = "n",
     ylab = "",
     xlab = "",
     col = "blue")

lines(seq(0, 1, 0.1),
      stats_N_NCI_prop_deciduous$stats["mean", ][[1]] + fixef(model_N_deciduous)[2]*seq(0, 1, 0.1),
      ylim = c(0, 0.05),
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

f_pctch(model_non_deciduous)

f_pctch(model_N_deciduous)


# barcharts of models ----------------------------------------
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
png("output/changes_growth.png", 
    width = 8, 
    height = 4, 
    units = 'in', 
    res = 300)
ggplot(data = changes, 
       aes(x = factor(group), y = change, fill = status, width = 0.65)) +
  geom_bar(stat = "identity", 
           position = position_dodge()) + 
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Group", y = "% Change", subtitle = "Growth") +
  geom_hline(yintercept = 0, color = "black") +
  geom_errorbar(aes(ymin = change-se, ymax = change+se), width = 0.2,
                position = position_dodge(width = 0.65)) +
  geom_text(aes(label = paste("n = ", obs, sep="")), 
            position = position_dodge(width = 0.65), vjust = -1)
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

png("output/changes1_growth.png", 
    width = 8, 
    height = 4, 
    units = 'in', 
    res = 300)
g1 <- ggplot(data = changes_age, 
             aes(x = factor(group), y = change, fill = status, width = 0.65)) +
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Group", y = "% Change") +
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
saveRDS(changes_form, "output/changes_form_growth.RDS")

png("output/changes_form_growth.png", 
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

# correlation between stand age and growth form
gs$age <- NULL
gs[!is.na(gs$STDAGE) && gs$STDAGE <= 60, "age"] <- 1
gs[!is.na(gs$STDAGE) && gs$STDAGE > 60, "age"] <- 2
cor.test(gs$DEC0EVER1, gs$age, use="complete.obs")

# Probe deciduous non-fixers with USDA -----------------------------
usda <- readRDS("data/derived/usda_numeric.RDS")

gd <- merge(deciduous, usda, by.x="SPCD", by.y="SPCD", all.x=T, all.y=F)

gdnf <- gd[gd$FIX == 0,]

grmod <- function(datapass) {
  model <- lmer(LGRp ~ NCI_props +
                               NCIs +
                               dbhs +
                               NCIs^2 +
                               NCI_props * NCIs +
                               dbhs * NCIs +
                               (1 | PCN1) + (1 | GENUS.x), 
                             data = datapass)
  return(model)
}

model_non <- grmod(datapass=gdnf[gdnf$SLOPE <= 8, ])
lowslope <- f_pctch(model_non)
lowslope

model_non <- grmod(datapass=gdnf[gdnf$SLOPE > 8, ])
highslope <- f_pctch(model_non)
highslope

model_non <- grmod(datapass=gdnf[gdnf$ASPECT >= 90 & gdnf$ASPECT <= 270, ])
south <- f_pctch(model_non)
south

model_non <- grmod(datapass=gdnf[gdnf$ASPECT < 90 | gdnf$ASPECT > 270, ])
north <- f_pctch(model_non)
north

model_non <- grmod(datapass=gdnf[gdnf$ELEVm <= 277, ])
lowelev <- f_pctch(model_non)
lowelev

model_non <- grmod(datapass=gdnf[gdnf$ELEVm > 277, ])
highelev <- f_pctch(model_non)
highelev

model_non <- grmod(datapass=gdnf[!is.na(gdnf$CNRatio) & gdnf$CNRatio == 1,])
highcn <- f_pctch(model_non)
highcn

model_non <- grmod(datapass=gdnf[!is.na(gdnf$CNRatio) & gdnf$CNRatio > 1,])
lowcn <- f_pctch(model_non)
lowcn

model_non <- grmod(datapass=gdnf[gdnf$FoliagePorosity == 1, ])
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
obs <- c(highslope[[3]],lowslope[[3]],
        north[[3]],south[[3]],
        highelev[[3]],lowelev[[3]],
        highcn[[3]],lowcn[[3]],
        highporosity[[3]],lowporosity[[3]],
        highgr[[3]],lowgr[[3]],
        tall[[3]],short[[3]],
        alleloy[[3]],allelon[[3]],
        leafry[[3]],leafrn[[3]],
        longlife[[3]],shortlife[[3]],
        highcatol[[3]],lowcatol[[3]],
        highdroughttol[[3]],lowdroughttol[[3]],
        highfert[[3]],lowfert[[3]],
        longroot[[3]],shortroot[[3]],
        highshadetol[[3]],lowshadetol[[3]],
        highseed[[3]],lowseed[[3]])
changes_decid <- as.data.frame(cbind(group, 
                                    status, 
                                    change,
                                    se,
                                    obs))
changes_decid$change <- as.numeric(as.character(changes_decid$change))
changes_decid$se <- as.numeric(as.character(changes_decid$se))
changes_decid$obs <- as.numeric(as.character(changes_decid$obs))
changes_decid$group <- factor(changes_decid$group, 
                              levels = c("slope", "aspect", "elevation", "CN ratio", 
                                         "foliar porosity", "growth rate", "mature height",
                                         "allelopathy", "leaf retention", "life span", 
                                         "CaCO3 tolerance", "drought tolerance", "fertility req",
                                         "root depth", "shade tolerance", "seed abundance"))
changes_decid
saveRDS(changes_decid, "output/changes_non_decid_growth.RDS")

png("output/changes_non_decid_growth.png", 
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
  # geom_text(aes(label = paste("n = ", obs, sep="")), 
  #           position = position_dodge(width = 0.65), vjust = -1)

#plot_grid(g1, g2, labels=c("", ""), ncol = 2, nrow = 1)
dev.off()

# Probe with Glopnet data --------------------------------------
# get data, from publication Wright et al 2004. DOI: [10.1038/nature02403](http://doi.org/10.1038/nature02403)
dat.Wright <- read.csv2("/Users/Anika/Downloads/nature02403-s2.csv", 
                        header = TRUE, 
                        skip = 10, 
                        stringsAsFactors = FALSE,
                        sep = ',',
                        na.strings = c("", "NA"))
dat.Wright <- dat.Wright[,c(1:18,21)]

# convert chr to numeric
dat.Wright$log.LL <- as.numeric(as.character(dat.Wright$log.LL))
dat.Wright$log.LMA <- as.numeric(as.character(dat.Wright$log.LMA))
dat.Wright$log.Nmass <- as.numeric(as.character(dat.Wright$log.Nmass))
dat.Wright$log.Narea <- as.numeric(as.character(dat.Wright$log.Narea))
dat.Wright$log.Pmass <- as.numeric(as.character(dat.Wright$log.Pmass))
dat.Wright$log.Parea <- as.numeric(as.character(dat.Wright$log.Parea))
dat.Wright$log.Amass <- as.numeric(as.character(dat.Wright$log.Amass))
dat.Wright$log.Aarea <- as.numeric(as.character(dat.Wright$log.Aarea))
dat.Wright$log.Gs <- as.numeric(as.character(dat.Wright$log.Gs))
dat.Wright$Ca...Ci <- as.numeric(as.character(dat.Wright$Ca...Ci))

# multiple entries per species, need to average them first
# ddply(d, .(Name), summarize,  Rate1=mean(Rate1), Rate2=mean(Rate2))
require(plyr)
wright.grp <- ddply(dat.Wright, .(Species), summarize,  
                    Decid.E.green = first(Decid.E.green), 
                    log.LL = mean(log.LL, na.rm=T),
                    log.LMA = mean(log.LMA, na.rm=T),
                    log.Nmass = mean(log.Nmass, na.rm=T),
                    log.Narea = mean(log.Narea, na.rm=T),
                    log.Pmass = mean(log.Pmass, na.rm=T),
                    log.Parea = mean(log.Parea, na.rm=T),
                    log.Amass = mean(log.Amass, na.rm=T),
                    log.Aarea = mean(log.Aarea, na.rm=T),
                    log.Gs = mean(log.Gs, na.rm=T))

gs$name <- paste(gs$GENUS, gs$SPECIES, sep=" ")

spec <- intersect(unique(gs$name), unique(wright.grp$Species))
length(spec) #105 species
length(unique(gs$name)) #360 species (only ~1/4 species in wright data)

#df <- deciduous[deciduous$FIX == 0,]

wright <- merge(gs, wright.grp, by.x="name", by.y="Species", all.x=T, all.y=F)
wright1 <- wright

evergreen <- wright1 %>% dplyr::filter(DEC0EVER1 == 1)
deciduous <- wright1 %>% dplyr::filter(DEC0EVER1 == 0)

#CN NCIprop interaction term
model <- lmer(LGRp ~ NCI_props +
                log.Nmass +
                NCI_props*log.Nmass +
                NCIs +
                dbhs +
                NCIs^2 +
                NCI_props * NCIs +
                dbhs * NCIs +
                (1 | PCN1) + (1 | GENUS), 
              data = wright)
coefs <- data.frame(coef(summary(model)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

#GF, C3C4, log.LL, log.LMA, log.Nmass, log.Narea, log.Pmass, log.Parea
#log.Amass, log.Aarea, log.Gs, 

grmod <- function(datapass) {
  model <- lmer(LGRp ~ NCI_props +
                  NCIs +
                  dbhs +
                  NCIs^2 +
                  NCI_props * NCIs +
                  dbhs * NCIs +
                  (1 | PCN1) + (1 | GENUS), 
                data = datapass)
  return(model)
}

model_non <- grmod(datapass = wright[wright$log.LL <= 1, ])
lowll <- f_pctch(model_non)
lowll

model_non <- grmod(datapass = wright[wright$log.LL > 1, ])
highll <- f_pctch(model_non)
highll

model_non <- grmod(datapass = wright[wright$log.LMA <= 2, ])
lowlma <- f_pctch(model_non)
lowlma

model_non <- grmod(datapass = wright[wright$log.LMA > 2, ])
highlma <- f_pctch(model_non)
highlma

model_non <- grmod(datapass = wright[wright$log.Nmass <= 0.2, ])
lownmass <- f_pctch(model_non)
lownmass

model_non <- grmod(datapass = wright[wright$log.Nmass > 0.2, ])
highnmass <- f_pctch(model_non)
highnmass

model_non <- grmod(datapass = wright[wright$log.Narea <= 0.2, ])
lownarea <- f_pctch(model_non)
lownarea

model_non <- grmod(datapass = wright[wright$log.Narea > 0.2, ])
highnarea <- f_pctch(model_non)
highnarea

model_non <- grmod(datapass = wright[wright$log.Pmass <= -0.8, ])
lowpmass <- f_pctch(model_non)
lowpmass

model_non <- grmod(datapass = wright[wright$log.Pmass > -0.8, ])
highpmass <- f_pctch(model_non)
highpmass

model_non <- grmod(datapass = wright[wright$log.Parea <= -0.9, ])
lowparea <- f_pctch(model_non)
lowparea

model_non <- grmod(datapass = wright[wright$log.Parea > -0.9, ])
highparea <- f_pctch(model_non)
highparea

model_non <- grmod(datapass = wright[wright$log.Aarea <= 0.9, ])
lowaarea <- f_pctch(model_non)
lowaarea

model_non <- grmod(datapass = wright[wright$log.Aarea > 0.9, ])
highaarea <- f_pctch(model_non)
highaarea

model_non <- grmod(datapass = wright[wright$log.Amass <= 1.8, ])
lowamass <- f_pctch(model_non)
lowamass

model_non <- grmod(datapass = wright[wright$log.Amass > 1.8, ])
highamass <- f_pctch(model_non)
highamass

model_non <- grmod(datapass = wright[wright$log.Gs <= 2.3, ])
lowgs <- f_pctch(model_non)
lowgs

model_non <- grmod(datapass = wright[wright$log.Gs > 2.3, ])
highgs <- f_pctch(model_non)
highgs

#ll, lma, nmass, narea, pmass, parea, aarea, amass, gs

group <- c("ll", "ll", "lma", "lma", "nmass", "nmass", "narea", "narea", "pmass", 
           "pmass", "parea", "parea", "amass", "amass", "aarea", "aarea", 
           "gs", "gs")
status <- rep(c("high", "low"), 9)
change <- c(highll[[1]],lowll[[1]],
            highlma[[1]],lowlma[[1]],
            highnmass[[1]],lownmass[[1]],
            highnarea[[1]],lownarea[[1]],
            highpmass[[1]],lowpmass[[1]],
            highparea[[1]],lowparea[[1]],
            highamass[[1]],lowamass[[1]],
            highaarea[[1]],lowaarea[[1]],
            highgs[[1]],lowgs[[1]])
se <- c(highll[[2]],lowll[[2]],
        highlma[[2]],lowlma[[2]],
        highnmass[[2]],lownmass[[2]],
        highnarea[[2]],lownarea[[2]],
        highpmass[[2]],lowpmass[[2]],
        highparea[[2]],lowparea[[2]],
        highamass[[2]],lowamass[[2]],
        highaarea[[2]],lowaarea[[2]],
        highgs[[2]],lowgs[[2]])
obs <- c(highll[[3]],lowll[[3]],
         highlma[[3]],lowlma[[3]],
         highnmass[[3]],lownmass[[3]],
         highnarea[[3]],lownarea[[3]],
         highpmass[[3]],lowpmass[[3]],
         highparea[[3]],lowparea[[3]],
         highamass[[3]],lowamass[[3]],
         highaarea[[3]],lowaarea[[3]],
         highgs[[3]],lowgs[[3]])
changes_decid <- as.data.frame(cbind(group, 
                                     status, 
                                     change,
                                     se,
                                     obs))
changes_decid$change <- as.numeric(as.character(changes_decid$change))
changes_decid$se <- as.numeric(as.character(changes_decid$se))
changes_decid$obs <- as.numeric(as.character(changes_decid$obs))
changes_decid$group <- factor(changes_decid$group, 
                              levels = c("ll", "lma", "nmass", "narea", "pmass", 
                                         "parea", "amass", "aarea", "gs"))
changes_decid

png("output/changes_wright_growth.png", 
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
# geom_text(aes(label = paste("n = ", obs, sep="")), 
#           position = position_dodge(width = 0.65), vjust = -1)

#plot_grid(g1, g2, labels=c("", ""), ncol = 2, nrow = 1)
dev.off()



# high leaf retention sp ---------------------------------
leafret <- gd[gd$LeafRetention == 1, ] #leaf retention is yes
unique(leafret$SPCD)
unique(leafret$COMMON_NAME)
unique(leafret$GENUS.x)
unique(leafret[c("GENUS.x", "SPECIES.x", "JENKINS_SPGRPCD")])

#classify species in category 10 as evergreen or decid
classten <- gd[gd$JENKINS_SPGRPCD == 10, ]
unique(classten[c("GENUS.x", "SPECIES.x", "SPCD")])

# Add leaf traits as fixed effects ----------------------------
model <- lmer(LGRp ~ NCI_props +
                NCIs +
                dbhs +
                NCIs^2 +
                log.LL +
                log.LMA +
                log.Nmass +
                log.Narea +
                log.Pmass +
                (1 | PCN1) +
                (NCI_props | DEC0EVER1), 
              data = wright)
coefplot(model)

ranef(model)

library(sjPlot)
library(sjmisc)
plot_model(model, type = "re", terms = "DEC0EVER1", axis.lim = c(-0.05, 0.05))


# mapping significance -----------------------------------

p <- fread("/Users/Anika/Documents/GradSchool/FIA_CompetitionModel/data/raw/PLOT.csv", 
           select = c("CN", "LAT", "LON"))

gs <- readRDS("output/gs_withgroups.RDS")
#gs <- merge(gs,p,by.x="PCN1.x",by.y="CN",all.x=T,all.y=F)
gs[is.na(gs$NCI_props), ]$NCI_props <- 0

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
data <- gs #gs for all, young, old, canopy, under
xname <- "gs"

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


mu <- 0.008925627 #same for whole data set (gs, ms, rs)
sdev <- 0.08225147


for(i in 1:nlon){
  print(i/nlon*100)
  lon.a <- lon.list[i] - dlon/2
  lon.b <- lon.list[i] + dlon/2
  for(j in 1:nlat){
    lat.a <- lat.list[j] - dlat/2
    lat.b <- lat.list[j] + dlat/2
    
    plots <- p[(p$lon1>=lon.a & p$lon1<lon.b & p$lat1>=lat.a & p$lat<lat.b),]$PCN
    if(length(plots)>0){
      temp <- data[which(data$PCN1 %in% plots | data$PCN2 %in% plots),
                   c("NCI_props","NCIs","dbhs","LGRp","FIX","ACT1vsRIZ0","SPCD","PCN1","GENUS")]
      if(nrow(temp)>0){
        if(length(unique(temp$PCN1))>5 && length(unique(temp$GENUS))>3){
          model = lmer(LGRp ~ NCI_props +
                                NCIs +
                                dbhs +
                                NCIs^2 +
                                NCI_props*NCIs +
                                dbhs*NCIs + 
                                (1|PCN1), 
                              data = temp)
          coefs <- data.frame(coef(summary(model)))
          coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
          ints <- coefs["(Intercept)","Estimate"] - (coefs["NCI_props","Estimate"]) * (mu/sdev)
          slopes <- (coefs["NCI_props","Estimate"] + coefs["NCI_props:NCIs","Estimate"]) / sdev
          pch <- 100*((ints + slopes) - ints)/ints
          
          sig.grid[i,j] <- coefs["NCI_props","p.z"]
          #sig.grid[i,j] <- coefs[2,4]
          coef.grid[i,j] <- coefs["NCI_props","Estimate"]
          #coef.grid[i,j] <- summary(model.fixed2)$coefficients[2,1]
          num.grid[i,j] <- sum(temp$FIX, na.rm=T)
          #num.grid[i,j] <- colSums(temp)[[4]]
          nump.grid[i,j] <- length(unique(temp$PCN1))
          riz.grid[i,j] <- mean(temp$ACT1vsRIZ0, na.rm=T)
          gen.grid[i,j] <- toString(tail(names(sort(table(temp$GENUS))), 1))
          #gen.grid[i,j] <- strtoi(names(sort(-table(temp$GENUS)))[1]) #most common genus
          se.grid[i,j] <- sqrt(diag(vcov(model)))[[2]]
          pct.grid[i,j] <- pch
          
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
      else if(pct.grid[i,j] < (-100)){
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

df <- cbind(as.vector(sig.grid.c), as.vector(coef.grid), as.vector(coef.grid.t),
            as.vector(gen.grid), as.vector(riz.grid), as.vector(se.grid))
df <- as.data.frame(df)
colnames(df) <- c("siggridc","coefgrid","coefgrid.t","gengrid","rizgrid","segrid")
df_sig <- df[!is.na(sig.grid.c),]

nlev <- 64
#ticks <- 2^(0:13)
ticks <- c(min(coef.grid,na.rm=T),-2,-1,0,1,2,max(coef.grid,na.rm=T))
image.plot(lon.list, lat.list, coef.grid, 
           nlevel=nlev, 
           col=tim.colors(nlev),
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
png("output/sigcoefs_growth.png", 
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
           ylim=c(23, 63),xlim=c(-140,-60))
#ylim=c(latminint-1,latmaxint+1),xlim=c(lonminint-2,lonmaxint+1))
title(main="Most Common Genus")
#mtext("(a)",font=2,adj=0.01,padj=-0.5)
map('world',regions='usa',add=TRUE)

# Pct change in recruitment rate (without outliers)
png("output/pctch_growth.png", 
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
rc1 = colorRampPalette(colors = c("red4", "#FDE2DF"), space="Lab")(nHalf)    
## Make vector of colors for values above threshold
rc2 = colorRampPalette(colors = c("#E0DFFD", "blue4"), space="Lab")(nHalf)
rampcols = c(rc1, rc2)
## In your example, this line sets the color for values between 49 and 51. 
rampcols[c(nHalf, nHalf+1)] = rgb(t(col2rgb("white")), maxColorValue=256) 
rb1 = seq(Min, Thresh, length.out=nHalf+1)
rb2 = seq(Thresh, Max, length.out=nHalf+1)[-1]
rampbreaks = c(rb1, rb2)
nlev <- 64
ticks <- c(-100,0,100)
ticks <- c(min(pct.grid.noout,na.rm=T), 0, max(pct.grid.noout,na.rm=T))
ticks <- c(-100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100)
image.plot(lon.list, lat.list, pct.grid.noout, 
           nlevel = nlev, 
           col = rampcols, 
           breaks = rampbreaks,
           axis.args = list(at=ticks,labels=ticks),
           xlab = "longitude", ylab = "latitude",
           ylim = c(23,63), xlim = c(-140,-65))
title(main = "Percent Change in Growth Rate")
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
png("output/semap_growth.png", 
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


###
# correlations ---------------------------------------
coef <- (as.vector(coef.grid))
nplot <- as.vector(nplot.grid)
nfix <- as.vector(num.grid)
rhiz <- as.vector(riz.grid)
water <- as.vector(water.grid)
se <- as.vector(se.grid)
pct <- as.vector(pct.grid)

cor.test(coef, nplot, use="complete.obs")
cor.test(coef,nfix,use="complete.obs")
cor.test(coef,rhiz,use="complete.obs")
cor.test(nplot,se,use="complete.obs")
cor.test(nfix,se,use="complete.obs")

cor.test(coef, se, use = "complete.obs")
cor.test(pct, se, use = "complete.obs")
cor.test(coef, pct, use = "complete.obs")

plot(nplot,se, xlab="Number of Plots in Grid Cell", ylab="Standard Error", main="Standard Error vs. Number of Plots")
plot(nfix,se, xlab="Number of Fixers in Grid Cell", ylab="Standard Error", main="Standard Error vs. Number of Fixers")

x <- as.vector(sig.grid.c)
x <- x[!is.na(x)]
hist(x, prob=T, breaks=30, xlim=c(-0.1,0.3), xlab="B1 for Grid Cell", 
     main="Histogram of Regression coefficients (significant only)")

a1 <- table(coef.grid.c)
aaf <- (a1[[2]]/(a1[[1]]+a1[[2]]))*100
aaf <- round(aaf, digits=2)
aaf

x <- as.vector(sig.grid.c)
x <- x[!is.na(x)]
length(x[x>0])/length(x)

# N dep --------------------------------------
require(raster)
Ndep <- raster("/Users/Anika/Documents/GradSchool/FIA_EstimateProj/output/N_dep.grd")

coords <- gs[,c("LON","LAT")]
dep <- raster::extract(Ndep, coords)

gs <- cbind(gs, dep)

# n dep in bins of 2
model1 <- grmod(datapass = gs[gs$dep > 0 & gs$dep < 2,])
model2 <- grmod(datapass = gs[gs$dep >= 2 & gs$dep < 4,])
model3 <- grmod(datapass = gs[gs$dep >= 4 & gs$dep < 6,])
model4 <- grmod(datapass = gs[gs$dep >= 6 & gs$dep < 8,])
model5 <- grmod(datapass = gs[gs$dep >= 8,])

N1 <- f_pctch(grmod(datapass = gs[gs$dep > 0 & gs$dep < 1,]))
N2 <- f_pctch(grmod(datapass = gs[gs$dep >= 1 & gs$dep < 2,]))
N3 <- f_pctch(grmod(datapass = gs[gs$dep >= 2 & gs$dep < 3,]))
N4 <- f_pctch(grmod(datapass = gs[gs$dep >= 3 & gs$dep < 4,]))
N5 <- f_pctch(grmod(datapass = gs[gs$dep >= 4 & gs$dep < 5,]))
N6 <- f_pctch(grmod(datapass = gs[gs$dep >= 5 & gs$dep < 6,]))
N7 <- f_pctch(grmod(datapass = gs[gs$dep >= 6 & gs$dep < 7,]))
N8 <- f_pctch(grmod(datapass = gs[gs$dep >= 7 & gs$dep < 8,]))
N9 <- f_pctch(grmod(datapass = gs[gs$dep >= 8 & gs$dep < 9,]))
N10 <- f_pctch(grmod(datapass = gs[gs$dep >= 9,]))

#n dep in quartiles
test <- gs %>% mutate(quartile = ntile(dep, 4))

model1 <- grmod(datapass = test[test$quartile == 1, ])
model2 <- grmod(datapass = test[test$quartile == 2, ])
model3 <- grmod(datapass = test[test$quartile == 3, ])
model4 <- grmod(datapass = test[test$quartile == 4, ])

N1 <- f_pctch(model1)
N2 <- f_pctch(model2)
N3 <- f_pctch(model3)
N4 <- f_pctch(model4)

ndep <- rbind.data.frame(N1, N2, N3, N4, N5, N6, N7, N8, N9, N10)
colnames(ndep) <- c("est", "se", "num")
quad <- c("0-1", "1-2", "2-3", "3-4", "4-5", "5-6", "6-7", "7-8", "8-9", "9+")
ndep <- cbind(ndep, quad)

ggplot(data = ndep, 
       aes(x = factor(quad), y = est, width = 0.65)) +
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Group", y = "% Change", title = "% Change in RGR by N deposition") +
  geom_hline(yintercept = 0, color = "black") + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
geom_errorbar(aes(ymin = est-se, ymax = est+se), width = 0.2,
              position = position_dodge(width = 0.65))

ndep_gr <- ndep

# Split by forest type ------------------------------------------------------

cond <- fread("data/raw/cond.csv", select = c('PLT_CN', 
                                              'FORTYPCD'))
cond <- data.table(cond)
FIA <- merge(gs, cond, by.x="PCN1", by.y="PLT_CN", all.x=T, all.y=F)
FIA$fortypsm <- trunc(FIA$FORTYPCD/10)*10
FIA$fortypsm <- as.integer(FIA$fortypsm)
FIA[(FIA$fortypsm %in% c(510, 520)), "fortypsm"] <- 500
FIA[FIA$fortypsm %in% c(720), "fortypsm"] <- 700
FIA[FIA$fortypsm %in% c(930), "fortypsm"] <- 920
FIA$fortypsm <- as.factor(FIA$fortypsm)

model <- lmer(LGRp ~ NCI_props +
                NCIs +
                dbhs +
                NCIs^2 +
                NCI_props * NCIs +
                dbhs * NCIs +
                (1 | PCN1) + (1 + NCI_props | fortypsm), 
              data = FIA)

fortype <- ranef(model)$fortypsm

out <- NULL
for(i in 1:nrow(fortype)){
  orig <- fortype[i, 1]
  final <- fortype[i, 1] + 1*fortype[i, 2]
  pctch <- 100*(final - orig)/orig
  out <- rbind(out, pctch)
}

fortype <- cbind(fortype, out)

cnttype <- dplyr::count(FIA, fortypsm) #get number in each group

forcodes <- c(100,120,140,150,160,170,180,200,220,240,260,280,300,320,340,360,370,380,390,
              400,500,600,700,800,900,910,920,940,960,970,980,988,990,999)
fornames <- c("jack pine","spruce/fir","longleaf/slash pine","tropical softwoods",
              "loblolly/shortleaf pine","other E softwoods", "pinyon/juniper",
              "douglas-fir","ponderosa pine","W white pine","fir/spruce/mountain hemlock",
              "lodgepole pine","hemlock/sitka spruce","W larch","redwood","other W softwoods",
              "CA mixed conifer","exotic softwoods","other softwoods","oak/pine","oak/hickory",
              "oak/gum/cypress","elm/ash/cottonwood","maple/beech/birch","aspen/birch",
              "alder/maple","W oak","tanoak/laurel","other hardwoods",
              "woodland hardwoods",
              "tropical hardwoods","cloud forest","exotic hardoods","nonstocked")
forcodes <- cbind.data.frame(forcodes,fornames) #get text names for types

fortype <- merge(fortype, cnttype, by.x="row.names", by.y="fortypsm")
colnames(fortype) <- c("fortypsm", "int", "slope", "pctch", "n")
fortype <- merge(fortype, forcodes, by.x="fortypsm", by.y="forcodes", all.x=T, all.y=F)

ggplot(data = fortype, 
       aes(x = factor(fornames), y = pctch, width = 0.65)) +
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Group", y = "% Change") +
  geom_hline(yintercept = 0, color = "black") + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) 
  geom_errorbar(aes(ymin = change-se, ymax = change+se), width = 0.2,
                position = position_dodge(width = 0.65))
  
sorted <- fortype[order(fortype$pctch),] 

# split data into forest type to find dominant fixer
temp <- FIA %>% filter(FIX==1) %>% group_by(fortypsm, Genus) %>% summarize(cnt=n())
out <- temp %>% group_by(fortypsm) %>% top_n(1, cnt)
colnames(out) <- c("fortypesm","domfixer","nfixer")

fortype <- merge(fortype, out, by.x = "fortypsm", by.y="fortypesm", all.x=T)
fortype %>% arrange(nfixer)
fortype <- fortype[!is.na(fortype$nfixer), ]
fortype <- fortype[!is.na(fortype$fornames), ]
fortype <- fortype[fortype$nfixer > 10, ]

fortype %>%
  arrange(domfixer, desc(nfixer)) %>%
  mutate(fornames=factor(fornames, levels=fornames)) %>%
  ggplot(aes(x = factor(fornames), y = pctch, width = 0.65, fill = domfixer)) +
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Group", y = "% Change") +
  geom_hline(yintercept = 0, color = "black") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_text(aes(label=nfixer), vjust=1)



# nested model --------------------------------
gs <- gs %>% mutate(ndep_quart = ntile(dep, 4))

model <- lmer(LGRp ~ NCI_props +
                NCIs +
                dbhs +
                NCIs^2 +
                NCI_props * NCIs +
                dbhs * NCIs +
                (1 | PCN1) + (NCI_props | DEC0EVER1:ndep_quart), 
              data = gs)

# species list (for TRY) --------------------------------------------------------------
splist <- paste(gs$GENUS, gs$SPECIES)
splist <- unique(splist)
splist <- as.data.frame(splist)
colnames(splist) <- c("spnames")

mydata <- read.delim("data/raw/TryAccSpecies.txt")

splist <- merge(splist, mydata, by.x="spnames", by.y="AccSpeciesName", all.x=T, all.y=F)

#write.csv(splist, "output/fia_splist.csv")

# mycorrhizae -------------------------------------------------

#from Wang and Qui
wang <- read.csv("data/raw/WangQiu_AM_EM_data.csv")
wang <- wang[ ,2:4]

splist <- gs[,c("GENUS","SPECIES")]
splist$spnames <- paste(splist$GENUS, splist$SPECIES)
splist <- splist[!duplicated(splist$spnames), ]

mrlist <- merge(splist, wang, by.x="spnames", by.y="Species", all.x=T, all.y=F)

mrlist1 <- mrlist[!is.na(mrlist$Mycostatus),]
mrlist2 <- mrlist1 %>% 
  group_by(GENUS) %>% 
  summarize(myco = first(Mycostatus)) #MR for genus level (can merge back on to mrlist)

inwang <- unique(mrlist2$GENUS)
infia <- unique(gs$GENUS)

missingmr <- subset(infia, !(infia %in% inwang))

#from Phillips
phil <- read.csv("data/raw/Phillips_AM_EM_data.csv")
phil <- phil[,3:5]
colnames(phil) <- c("GENUS", "SPECIES", "Mycostatus")
phil$spnames <- paste(phil$GENUS, phil$SPECIES)

mrlist0 <- merge(mrlist, phil, by="spnames", all.x=T, all.y=F)
mrlist1 <- mrlist0[!is.na(mrlist0$Mycostatus.x) | !is.na(mrlist0$Mycostatus.y),]
mrlist1 <- mrlist1[,c(1:5,8)]
colnames(mrlist1) <- c("spnames", "GENUS", "SPECIES", "Mycostatus_wang", "References",
                       "Mycostatus_phil")
mrlist2 <- mrlist1 %>% 
  group_by(GENUS) %>% 
  summarize(myco_wang = first(Mycostatus_wang),
            myco_phil = first(Mycostatus_phil)) #MR for genus level (can merge back on to mrlist)

#use myco_phil because has value for everything 

gsmr <- merge(gs, mrlist2, by="GENUS", all.x=T, all.y=F)

model <- lmer(LGRp ~ NCI_props +
                NCIs +
                dbhs +
                NCIs^2 +
                NCI_props * NCIs +
                dbhs * NCIs +
                (1 | PCN1) + (1 + NCI_props | myco_phil), 
              data = gsmr)

# run model with groups as random effect ---------------------------
gs$group <- gs$DEC0EVER1*5 + gs$FIX*3
plyr::count(gs, "group")
gs$group <- as.factor(gs$group)
levels(gs$group)[1] <- decnon
levels(gs$group)[2]
