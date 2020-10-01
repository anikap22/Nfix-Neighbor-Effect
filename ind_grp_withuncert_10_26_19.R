



require(dplyr)
library(lme4)
require(lmerTest)
require(ggplot2)
require(GGally)
require(car)
require(tidyverse)
require(svMisc)

setwd("/Users/Anika/Documents/GradSchool/FIA_CompetitionModel/")

options(scipen = 6)

# growth ----------------------------------------------------------
gs <- readRDS("output/gs_withgroups.RDS")

gs <- gs[!is.na(gs$LGRp),]

gs$FIXg <- as.factor(gs$FIX)
gs$CANg <- as.factor(gs$NON0CAN1)
gs$DECg <- gs$Decid0Ever1
gs$DECg[gs$DECg < 0.6] <- 0
gs$DECg <- as.factor(gs$DECg)
gs <- gs %>% mutate(MRg = case_when(myco_phil == "AM" ~ 0, 
                                    myco_phil == "AM+EM" ~ 0,
                                    myco_phil == "EM" ~ 1,
                                    myco_phil == "ERM" ~ 1,
                                    myco_phil == "NM" ~ 2))
gs$MRg <- as.factor(gs$MRg)
gs <- gs %>% mutate(DEPg = case_when(dep < 3.96 ~ 0, 
                                     dep >= 3.96 ~ 1))
gs$DEPg <- as.factor(gs$DEPg)
gs <- gs %>% mutate(SMg = case_when(sm < 0.27 ~ 0, 
                                    sm >= 0.27 ~ 1))
gs$SMg <- as.factor(gs$SMg)
gs <- gs %>% mutate(CNg = case_when(Tcnratio == "low" ~ 0, 
                                    Tcnratio == "high" ~ 1))
gs$CNg <- as.factor(gs$CNg)
gs <- gs %>% mutate(ALEg = case_when(Tallelopathy == "no" ~ 0, 
                                    Tallelopathy == "yes" ~ 1))
gs$ALEg <- as.factor(gs$ALEg)
gs <- gs %>% mutate(AGEg = case_when(STDAGE < 60 ~ 0,
                                    STDAGE >= 60 ~ 1))
gs$AGEg <- as.factor(gs$AGEg)

# whole forest model
m <- lmer(LGRp ~ NCI_prop + 
            NCI +
            NCI^2+
            NCI*NCI_prop +
            NCI*dbh +
            dbh +
            (1 | pcn1),
          data = gs)
cat("Individual RGR Model:\n", file = "output/table_si3.txt", append = TRUE)
capture.output(summary(m), 
               file = "output/table_si3.txt", append = TRUE)
capture.output(r.squaredGLMM(m), 
               file = "output/table_si3.txt", append = TRUE)
capture.output(confint(m), 
               file = "output/table_si3.txt", append = TRUE)


out <- f_pch_nonpar(gs, "growth", 10, "NA")
ch_rgr <- pct_change_func_byrow_forest(out)
ch_rgr

out <- f_pch_nonpar(rs, "recr", 10, "NA")
ch_recr <- pct_change_func_byrow_forest(out)
ch_recr

out <- f_pch_nonpar(ms, "surv", 10, "NA")
ch_surv <- pct_change_func_byrow_forest(out)
ch_surv


# how many bootstraps do we need?
test <- f_pch_nonpar(data=gs, type="growth", nboot=10, group="FIXg")
m10.0 <- mean(test$group0low)
m10.1 <- mean(test$group1low)
s10.0 <- sd(test$group0low)
s10.1 <- sd(test$group1low)
test <- f_pch_nonpar(data=gs, type="growth", nboot=100, group="FIXg")
m100.0 <- mean(test$group0low)
m100.1 <- mean(test$group1low)
s100.0 <- sd(test$group0low)
s100.1 <- sd(test$group1low)
test <- f_pch_nonpar(data=gs, type="growth", nboot=1000, group="FIXg")
m1000.0 <- mean(test$group0low)
m1000.1 <- mean(test$group1low)
s1000.0 <- sd(test$group0low)
s1000.1 <- sd(test$group1low)

g1 <- cbind.data.frame(rbind(m10.1, m100.1, m1000.1), rbind(s10.1, s100.1, s1000.1), rbind(10, 100, 377))
colnames(g1) <- c("means","sd","nbootstraps")
plot(g1$nbootstraps, g1$means, ylim=c(0.01,0.02), xlab="n bootstraps", ylab="LGRp")
arrows(g1$nbootstraps, g1$means-g1$sd, g1$nbootstraps, g1$means+g1$sd, length=0.05, angle=90, code=3)

#### run non-parametric bootstrap
out_all <- f_pch_nonpar(data=gs, type="growth", nboot=100, group="NA")
ch_all <- pct_change_func(out_all)

system.time(
out_fix <- f_pch_nonpar(data=gs, type="growth", nboot=10, group="FIXg"))
ch_fix <- pct_change_func_byrow(out_fix)
f <- xtabs(~FIXg, gs) #f also called cnts
ch_fix$Freq <- rbind(f[[1]], f[[2]]) 

out_can <- f_pch_nonpar(data=gs, type="growth", nboot=10, group="CANg")
f <- xtabs(~CANg, gs)
ch_can <- pct_change_func_byrow(out_can)
(cnts[1]/nrow(gs))*inddat$mean[3]+(cnts[2]/nrow(gs))*inddat$mean[4]
ch_can$Freq <- rbind(f[[1]], f[[2]]) 

out_decid <- f_pch_nonpar(data=gs, type="growth", nboot=10, group="DECg")
f <- xtabs(~DECg, gs)
ch_decid <- pct_change_func_byrow(out_decid)
ch_decid$Freq <- rbind(f[[1]], f[[2]]) 

#out_dep <- f_pch_nonpar(data=gs, type="growth", nboot=10, group="DEPg")
#xtabs(~DEPg, gs)
#ch_dep <- pct_change_func_byrow(out_dep)

#out_sm <- f_pch_nonpar(data=gs, type="growth", nboot=10, group="SMg")
#xtabs(~SMg, gs)
#ch_sm <- pct_change_func_byrow(out_sm)

out_cn <- f_pch_nonpar(data=gs, type="growth", nboot=10, group="CNg")
f <- xtabs(~CNg, gs)
ch_cn <- pct_change_func_byrow(out_cn)
ch_cn$Freq <- rbind(f[[1]], f[[2]]) 

#out_ale <- f_pch_nonpar(data=gs, type="growth", nboot=10, group="ALEg")
#xtabs(~ALEg, gs)
#ch_ale <- pct_change_func_byrow(out_ale)

out_myco <- f_pch_nonpar(data=gs, type="growth", nboot=10, group="MRg", ngroup=3)
f <- xtabs(~MRg, gs)
ch_myco <- pct_change_func_byrow_3g(out_myco)
ch_myco$Freq <- rbind(f[[1]], f[[2]], f[[3]]) 

ch_fix$var <- "fix"
ch_can$var <- "canopy"
ch_decid$var <- "deciduousness"
ch_cn$var <- "foliar C:N"
ch_myco$var <- "mycorrhizae"
inddat <- rbind(ch_fix, ch_can, ch_decid, ch_cn, ch_myco)
inddat$type <- "RGR"

saveRDS(out_fix, "output/boots/out_fix_rgr_10bs.RDS")
saveRDS(out_can, "output/boots/out_can_rgr_10bs.RDS")
saveRDS(out_decid, "output/boots/out_decid_rgr_10bs.RDS")
saveRDS(out_cn, "output/boots/out_cn_rgr_10bs.RDS")
saveRDS(out_myco, "output/boots/out_myco_rgr_10bs.RDS")

system.time(
  out_fix <- f_pch_nonpar(data=rs, type="recr", nboot=100, group="FIXg"))
out_fix <- readRDS("output/boots/out_fix_recr_100bs.RDS")
ch_fix_r <- pct_change_func_byrow(out_fix)
f <- xtabs(~FIXg, rs)
ch_fix_r$Freq <- rbind(f[[1]], f[[2]]) 

out_can <- f_pch_nonpar(data=rs, type="recr", nboot=100, group="CANg")
out_can <- readRDS("output/boots/out_can_recr_100bs.RDS")
f <- xtabs(~CANg, rs)
ch_can_r <- pct_change_func_byrow(out_can)
ch_can_r$Freq <- rbind(f[[1]], f[[2]])

out_decid <- f_pch_nonpar(data=rs, type="recr", nboot=100, group="DECg")
out_decid <- readRDS("output/boots/out_decid_recr_100bs.RDS")
f <- xtabs(~DECg, rs)
ch_decid_r <- pct_change_func_byrow(out_decid)
ch_decid_r$Freq <- rbind(f[[1]], f[[2]]) 

out_cn <- f_pch_nonpar(data=rs, type="recr", nboot=100, group="CNg")
out_cn <- readRDS("output/boots/out_cn_recr_100bs.RDS")
f <- xtabs(~CNg, rs)
ch_cn_r <- pct_change_func_byrow(out_cn)
ch_cn_r$Freq <- rbind(f[[1]], f[[2]]) 

out_myco <- f_pch_nonpar(data=rs, type="growth", nboot=10, group="MRg", ngroup=3)
out_myco <- readRDS("output/boots/out_myco_recr_10bs.RDS")
f <- xtabs(~MRg, rs)
ch_myco_r <- pct_change_func_byrow_3g(out_myco)
ch_myco_r$Freq <- rbind(f[[1]], f[[2]], f[[3]]) 

ch_fix_r$var <- "fix"
ch_can_r$var <- "canopy"
ch_decid_r$var <- "deciduousness"
ch_cn_r$var <- "foliar C:N"
ch_myco_r$var <- "mycorrhizae"
recdat <- rbind(ch_fix_r, ch_can_r, ch_decid_r, ch_cn_r, ch_myco_r)
recdat$type <- "Recruitment"
recdat$level <- as.factor(recdat$group)
inddat <- rbind(inddat, recdat)

saveRDS(out_fix, "output/boots/out_fix_recr_100bs.RDS")
saveRDS(out_can, "output/boots/out_can_recr_100bs.RDS")
saveRDS(out_decid, "output/boots/out_decid_recr_100bs.RDS")
saveRDS(out_cn, "output/boots/out_cn_recr_100bs.RDS")
saveRDS(out_myco, "output/boots/out_myco_recr_10bs.RDS")

system.time(
  out_fix <- f_pch_nonpar(data=ms, type="surv", nboot=10, group="FIXg"))
out_fix <- readRDS("output/boots/out_fix_surv_10bs.RDS")
ch_fix_s <- pct_change_func_byrow(out_fix)
f <- xtabs(~FIXg, ms)
ch_fix_s$Freq <- rbind(f[[1]], f[[2]]) 

out_can <- f_pch_nonpar(data=ms, type="surv", nboot=10, group="CANg")
out_can <- readRDS("output/boots/out_can_surv_10bs.RDS")
f <- xtabs(~CANg, ms)
ch_can_s <- pct_change_func_byrow(out_can)
ch_can_s$Freq <- rbind(f[[1]], f[[2]]) 

out_decid <- f_pch_nonpar(data=ms, type="surv", nboot=10, group="DECg")
out_decid <- readRDS("output/boots/out_decid_surv_10bs.RDS")
f <- xtabs(~DECg, ms)
ch_decid_s <- pct_change_func_byrow(out_decid)
ch_decid_s$Freq <- rbind(f[[1]], f[[2]]) 

out_cn <- f_pch_nonpar(data=ms, type="surv", nboot=10, group="CNg")
out_cn <- readRDS('output/boots/out_cn_surv_10bs.RDS')
f <- xtabs(~CNg, ms)
ch_cn_s <- pct_change_func_byrow(out_cn)
ch_cn_s$Freq <- rbind(f[[1]], f[[2]]) 

out_myco <- f_pch_nonpar(data=ms, type="surv", nboot=10, group="MRg", ngroup=3)
out_myco <- readRDS("output/boots/out_myco_surv_10bs.RDS")
f <- xtabs(~MRg, ms)
ch_myco_s <- pct_change_func_byrow_3g(out_myco)
ch_myco_s$Freq <- rbind(f[[1]], f[[2]], f[[3]]) 

ch_fix_s$var <- "fix"
ch_can_s$var <- "canopy"
ch_decid_s$var <- "deciduousness"
ch_cn_s$var <- "foliar C:N"
ch_myco_s$var <- "mycorrhizae"
survdat <- rbind(ch_fix_s, ch_can_s, ch_decid_s, ch_cn_s, ch_myco_s)
survdat$type <- "Survival"
survdat$level <- as.factor(survdat$group)
inddat <- rbind(inddat, survdat)


saveRDS(out_fix, "output/boots/out_fix_surv_10bs.RDS")
saveRDS(out_can, "output/boots/out_can_surv_10bs.RDS")
saveRDS(out_decid, "output/boots/out_decid_surv_10bs.RDS")
saveRDS(out_cn, "output/boots/out_cn_surv_10bs.RDS")
saveRDS(out_myco, "output/boots/out_myco_surv_10bs.RDS")

saveRDS(inddat, "output/boots/inddat_10bs.RDS")

out_rgr <- f_pch_nonpar(data=gs, type="growth", nboot=2, group="NA")
ch_rgr <- pct_change_func_byrow_forest(out_rgr)

# functions -------------------------------------------------------------------

# predict using non-parametric bootstrap (resample)
# sample with replacement
# sample same number of plots as have in dataset
lmerControl(check.nobs.vs.nlev="ignore") #so when plots with too few trees picked its ok
require(MuMIn)
f_pch_nonpar <- function(data, type, nboot, group, ngroup=2){
  pred.y.0 <- matrix(ncol=2,nrow=nboot)
  pred.y.1 <- matrix(ncol=2,nrow=nboot)
  if(ngroup==3){
    pred.y.2 <- matrix(ncol=2,nrow=nboot)
  }

  for(i in 1:nboot){
    progress(i, max.value=nboot)
    
    # make sample
    inds <- sample(1:nrow(data), nrow(data), replace=T)
    dat <- data[inds,]
    
    # if no group (whole forest)
    if(group == "NA"){
      
      # run model 
      if(type == "growth"){
        keep <- c("pcn1","NCI_props","LGRp","NCIs","dbhs")
        dat <- dat %>% dplyr::select(keep)
        m <- lmer(LGRp ~ NCI_props + 
                    NCIs +
                    NCIs^2+
                    NCIs*NCI_props +
                    NCIs*dbhs +
                    dbhs +
                    (1 | pcn1),
                         data = dat)
        new_d <- data.frame(NCI_props=c(min(data$NCI_props), max(data$NCI_props)),
                            NCIs=c(mean(data$NCIs), mean(data$NCIs)),
                            dbhs=c(mean(data$dbhs), mean(data$dbhs)))
        new_d$Pred <- predict(m,new_d,re.form=~0)
        pred.y.0[i,1] <- new_d$Pred[1]
        pred.y.0[i,2] <- new_d$Pred[2]
      }else if(type == "recr" | type == "recruitment"){
        keep <- c("pcn1","NCI_props","NCIs","dbhs","recir")
        dat <- dat %>% dplyr::select(keep)
        m <- lmer(recir ~ NCI_props +
                    NCIs +
                    NCIs^2 +
                    NCIs*NCI_props +
                    dbhs +
                    (1 | pcn1),
                         data = dat)
        new_d <- data.frame(NCI_props=c(min(data$NCI_props), max(data$NCI_props)),
                            NCIs=c(mean(data$NCIs), mean(data$NCIs)),
                            dbhs=c(mean(data$dbhs), mean(data$dbhs)))
        new_d$Pred <- predict(m,new_d,re.form=~0)
        pred.y.0[i,1] <- new_d$Pred[1]
        pred.y.0[i,2] <- new_d$Pred[2]
      }else if(type == "mort" | type == "mortality" | type == "surv"){
        keep <- c("pcn1","NCI_props","NCIs","dbhs","survr")
        dat <- dat %>% dplyr::select(keep)
        m <- lmer(survr ~ NCI_props +
                    NCIs +
                    NCIs^2 +
                    dbhs +
                    NCIs*NCI_props +
                    NCIs*dbhs +
                    (1 | pcn1),
                  data = dat)
        new_d <- data.frame(NCI_props=c(min(data$NCI_props,na.rm=T), max(data$NCI_props,na.rm=T)),
                            NCIs=c(mean(data$NCIs,na.rm=T), mean(data$NCIs,na.rm=T)),
                            dbhs=c(mean(data$dbhs,na.rm=T), mean(data$dbhs,na.rm=T)))
        new_d$Pred <- predict(m,new_d,re.form=~0)
        pred.y.0[i,1] <- new_d$Pred[1]
        pred.y.0[i,2] <- new_d$Pred[2]
      }
    }else{ #if group (split into two groups)
      # run model & predict
      if(type == "growth"){
        keep <- c("pcn1","NCI_props","LGRp","NCIs","dbhs", group)
        dat <- dat %>% dplyr::select(keep)
        colnames(dat) <- c("pcn1","NCI_props","LGRp","NCIs","dbhs","group")
        m <- lmer(LGRp ~ NCI_props + 
                    NCIs +
                    NCIs^2+
                    NCIs*NCI_props +
                    NCIs*dbhs +
                    dbhs +
                    NCI_props*group +
                    (1 | pcn1),
                  data = dat)
        if(ngroup == 3){
          new_d <- data.frame(NCI_props=rep(c(min(data$NCI_props), max(data$NCI_props)),3),
                              NCIs=rep(mean(data$NCIs), 6),
                              dbhs=rep(mean(data$dbhs), 6),
                              group=as.factor(c(0,0,1,1,2,2)))
          new_d$Pred <- predict(m,new_d,re.form=~0)
          pred.y.0[i,1] <- new_d$Pred[1]
          pred.y.0[i,2] <- new_d$Pred[2]
          pred.y.1[i,1] <- new_d$Pred[3]
          pred.y.1[i,2] <- new_d$Pred[4]
          pred.y.2[i,1] <- new_d$Pred[5]
          pred.y.2[i,2] <- new_d$Pred[6]
        }else{
          new_d <- data.frame(NCI_props=rep(c(min(data$NCI_props), max(data$NCI_props)),2),
                              NCIs=rep(mean(data$NCIs), 4),
                              dbhs=rep(mean(data$dbhs), 4),
                              group=as.factor(c(0,0,1,1)))
          new_d$Pred <- predict(m,new_d,re.form=~0)
          pred.y.0[i,1] <- new_d$Pred[1]
          pred.y.0[i,2] <- new_d$Pred[2]
          pred.y.1[i,1] <- new_d$Pred[3]
          pred.y.1[i,2] <- new_d$Pred[4]
        }
      }else if(type == "recr" | type == "recruitment"){
        keep <- c("pcn1","NCI_props","recir","NCIs","dbhs", group)
        dat <- dat %>% dplyr::select(keep)
        colnames(dat) <- c("pcn1","NCI_props","recir","NCIs","dbhs","group")
        m <- lmer(recir ~ NCI_props +
                    NCIs +
                    NCIs^2 +
                    NCIs*NCI_props +
                    dbhs +
                    NCI_props*group +
                    (1 | pcn1),
                  data = dat)
        if(ngroup == 3){
          new_d <- data.frame(NCI_props=rep(c(min(data$NCI_props), max(data$NCI_props)),3),
                              NCIs=rep(mean(data$NCIs), 6),
                              dbhs=rep(mean(data$dbhs), 6),
                              group=as.factor(c(0,0,1,1,2,2)))
          new_d$Pred <- predict(m,new_d,re.form=~0)
          pred.y.0[i,1] <- new_d$Pred[1]
          pred.y.0[i,2] <- new_d$Pred[2]
          pred.y.1[i,1] <- new_d$Pred[3]
          pred.y.1[i,2] <- new_d$Pred[4]
          pred.y.2[i,1] <- new_d$Pred[5]
          pred.y.2[i,2] <- new_d$Pred[6]
        }else{
          new_d <- data.frame(NCI_props=rep(c(min(data$NCI_props,na.rm=T), max(data$NCI_props,na.rm=T)),2),
                              NCIs=rep(mean(data$NCIs,na.rm=T), 4),
                              dbhs=rep(mean(data$dbhs,na.rm=T), 4),
                              group=as.factor(c(0,0,1,1)))
          new_d$Pred <- predict(m,new_d,re.form=~0)
          pred.y.0[i,1] <- new_d$Pred[1]
          pred.y.0[i,2] <- new_d$Pred[2]
          pred.y.1[i,1] <- new_d$Pred[3]
          pred.y.1[i,2] <- new_d$Pred[4]
        }
      }else if(type == "mort" | type == "mortality" | type == "surv"){
        keep <- c("pcn1","NCI_props","survr","NCIs","dbhs", group)
        dat <- dat %>% dplyr::select(keep)
        colnames(dat) <- c("pcn1","NCI_props","survr","NCIs","dbhs","group")
        m <- lmer(survr ~ NCI_props +
                    NCIs +
                    NCIs^2 +
                    dbhs +
                    NCIs*NCI_props +
                    NCIs*dbhs +
                    NCI_props*group +
                    (1 | pcn1),
                  data = dat)
        if(ngroup == 3){
          new_d <- data.frame(NCI_props=rep(c(min(data$NCI_props,na.rm=T), 
                                              max(data$NCI_props,na.rm=T)),3),
                              NCIs=rep(mean(data$NCIs,na.rm=T), 6),
                              dbhs=rep(mean(data$dbhs,na.rm=T), 6),
                              group=as.factor(c(0,0,1,1,2,2)))
          new_d$Pred <- predict(m,new_d,re.form=~0)
          pred.y.0[i,1] <- new_d$Pred[1]
          pred.y.0[i,2] <- new_d$Pred[2]
          pred.y.1[i,1] <- new_d$Pred[3]
          pred.y.1[i,2] <- new_d$Pred[4]
          pred.y.2[i,1] <- new_d$Pred[5]
          pred.y.2[i,2] <- new_d$Pred[6]
        }else{
          new_d <- data.frame(NCI_props=rep(c(min(data$NCI_props,na.rm=T), max(data$NCI_props,na.rm=T)),2),
                              NCIs=rep(mean(data$NCIs,na.rm=T), 4),
                              dbhs=rep(mean(data$dbhs,na.rm=T), 4),
                              group=as.factor(c(0,0,1,1)))
          new_d$Pred <- predict(m,new_d,re.form=~0)
          pred.y.0[i,1] <- new_d$Pred[1]
          pred.y.0[i,2] <- new_d$Pred[2]
          pred.y.1[i,1] <- new_d$Pred[3]
          pred.y.1[i,2] <- new_d$Pred[4]
        }
        
      }
    }
  }#end for loop
    
  r2 <- r.squaredGLMM(m)
  
  if(group == "NA"){
    result <- cbind.data.frame(pred.y.0, r2)
    colnames(result) <- c("group0low", "group0high", "R2m", "R2c")
  }else{
    if(ngroup==3){
      result <- cbind.data.frame(pred.y.0, pred.y.1, pred.y.2, r2)
      colnames(result) <- c("group0low", "group0high", "group1low", "group1high", 
                            "group2low", "group2high", "R2m", "R2c")
    }else{
      result <- cbind.data.frame(pred.y.0, pred.y.1, r2)
      colnames(result) <- c("group0low", "group0high", "group1low", "group1high", "R2m", "R2c")
    }
  }
  

  
  return(result)
  
}#end function

test <- f_pch_nonpar(data=gs, type="growth", nboot=3, group="FIXg")



pct_change_func_byrow <- function(data){
  data$pch0 <- 100*(data$group0high-data$group0low)/data$group0low
  data$pch1 <- 100*(data$group1high-data$group1low)/data$group1low
  
  pch0_mean <- mean(data$pch0)
  pch1_mean <- mean(data$pch1)
  pch0_sd <- sd(data$pch0)
  pch1_sd <- sd(data$pch1)
  pch0_lower <- unname(quantile(data$pch0,.025))
  pch0_upper <- unname(quantile(data$pch0,0.975))
  pch1_lower <- unname(quantile(data$pch1,.025))
  pch1_upper <- unname(quantile(data$pch1,0.975))
  
  g0 <- cbind(pch0_mean, pch0_lower, pch0_upper, 0)
  colnames(g0) <- c("mean","lower","upper","group")
  g1 <- cbind(pch1_mean, pch1_lower, pch1_upper, 1)
  colnames(g1) <- c("mean","lower","upper","group")
  result <- rbind.data.frame(g0, g1)
  
  return(result)
}

pct_change_func_byrow_3g <- function(data){
  data$pch0 <- 100*(data$group0high-data$group0low)/data$group0low
  data$pch1 <- 100*(data$group1high-data$group1low)/data$group1low
  data$pch2 <- 100*(data$group2high-data$group2low)/data$group2low
  
  pch0_mean <- mean(data$pch0)
  pch1_mean <- mean(data$pch1)
  pch2_mean <- mean(data$pch2)
  pch0_sd <- sd(data$pch0)
  pch1_sd <- sd(data$pch1)
  pch2_sd <- sd(data$pch2)
  pch0_lower <- unname(quantile(data$pch0,.025))
  pch0_upper <- unname(quantile(data$pch0,0.975))
  pch1_lower <- unname(quantile(data$pch1,.025))
  pch1_upper <- unname(quantile(data$pch1,0.975))
  pch2_lower <- unname(quantile(data$pch2,.025))
  pch2_upper <- unname(quantile(data$pch2,0.975))
  
  g0 <- cbind(pch0_mean, pch0_lower, pch0_upper, 0)
  colnames(g0) <- c("mean","lower","upper","group")
  g1 <- cbind(pch1_mean, pch1_lower, pch1_upper, 1)
  colnames(g1) <- c("mean","lower","upper","group")
  g2 <- cbind(pch2_mean, pch2_lower, pch2_upper, 2)
  colnames(g2) <- c("mean", 'lower','upper','group')
  result <- rbind.data.frame(g0, g1, g2)
  
  
  return(result)
}

pct_change_func_byrow_forest <- function(data){
  data$pch0 <- 100*(data$group0high-data$group0low)/data$group0low

  pch0_mean <- mean(data$pch0)
  pch0_sd <- sd(data$pch0)
  pch0_lower <- unname(quantile(data$pch0,.025))
  pch0_upper <- unname(quantile(data$pch0,0.975))

  g0 <- cbind(pch0_mean, pch0_lower, pch0_upper)
  colnames(g0) <- c("mean","lower","upper")

  result <- g0
  
  return(result)
}




f_plot <- function(modelout, type, title){
  plot <- ggplot(data = modelout, 
                 aes(x = factor(group), y = pchc, width = 0.65)) +
    geom_bar(stat = "identity", 
             position = position_dodge()) + 
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Group", y = "% Change", subtitle = type, title = title,
         caption="Number in 1000s of samples") +
    geom_hline(yintercept = 0, color = "black") +
    geom_errorbar(aes(ymin=pchl, ymax=pchu), width=.2,
                  position=position_dodge(.9)) +
    geom_text(aes(label=round(Freq/1000, 1)), position=position_dodge(0.9), vjust=-0.2) +
    theme(panel.background = element_rect(fill = 'white', colour = 'white'),
          plot.title = element_text(hjust = 0.5, size = 20),
          text = element_text(size = 20),
          axis.text.x = element_text(angle = 90, hjust = 1, size = 20),
          plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0))) 
  return(plot)
}

# recr ---------------------------------

rs <- readRDS("output/rs_withgroups.RDS")

rs$FIXg <- as.factor(rs$FIX)
rs$CANg <- as.factor(rs$NON0CAN1)
rs$DECg <- rs$Decid0Ever1
rs$DECg[rs$DECg < 0.6] <- 0
rs$DECg <- as.factor(rs$DECg)
rs <- rs %>% mutate(MRg = case_when(myco_phil == "AM" ~ 0, 
                                    myco_phil == "AM+EM" ~ 0,
                                    myco_phil == "EM" ~ 1,
                                    myco_phil == "ERM" ~ 1,
                                    myco_phil == "NM" ~ 2))
rs$MRg <- as.factor(rs$MRg)
rs <- rs %>% mutate(DEPg = case_when(dep < 3.96 ~ 0, 
                                     dep >= 3.96 ~ 1))
rs$DEPg <- as.factor(rs$DEPg)
rs <- rs %>% mutate(SMg = case_when(sm < 0.27 ~ 0, 
                                    sm >= 0.27 ~ 1))
rs$SMg <- as.factor(rs$SMg)
rs <- rs %>% mutate(CNg = case_when(Tcnratio == "low" ~ 0, 
                                    Tcnratio == "high" ~ 1))
rs$CNg <- as.factor(rs$CNg)
rs <- rs %>% mutate(ALEg = case_when(Tallelopathy == "no" ~ 0, 
                                     Tallelopathy == "yes" ~ 1))
rs$ALEg <- as.factor(rs$ALEg)
rs <- rs %>% mutate(AGEg = case_when(STDAGE < 60 ~ 0,
                                     STDAGE >= 60 ~ 1))
rs$AGEg <- as.factor(rs$AGEg)

# remove outliers
rsmini <- rs[rs$NCI_total2 < 25,]
rsmini <- rsmini[!is.na(rsmini$MRg) & !is.na(rsmini$DEPg) & !is.na(rsmini$SMg) & !is.na(rsmini$CNg) & !is.na(rsmini$AGEg) & !is.na(rsmini$dbh2) & !is.na(rsmini$NCI_total2),]
#83.557% of rows are kept

rsmini$NCIs <- (rsmini$NCI_total2 - mean(rsmini$NCI_total2))/sd(rsmini$NCI_total2)
rsmini$dbhs <- (rsmini$dbh2 - mean(rsmini$dbh2))/sd(rsmini$dbh2)

rs <- rsmini

# whole forest model
m <- lmer(recir ~ NCI_props +
            NCIs +
            NCIs^2 +
            NCIs*NCI_props +
            dbhs +
            (1 | pcn1),
          data = rs)

cat("Individual Recruitment Model:\n", file = "output/table_si3.txt", append = TRUE)
capture.output(summary(m), 
               file = "output/table_si3.txt", append = TRUE)
capture.output(r.squaredGLMM(m), 
               file = "output/table_si3.txt", append = TRUE)
capture.output(confint(m), 
               file = "output/table_si3.txt", append = TRUE)


#
p1 <- f_plot(out_fix, "Recruitment", "fixer vs nonfixer")
p2 <- f_plot(out_can, "Recruitment", "canopy vs noncanopy")
p3 <- f_plot(out_decid, "Recruitment", "deciduous vs evergreen")
p4 <- f_plot(out_dep, "Recruitment", "low vs high N deposition")
p5 <- f_plot(out_sm, "Recruitment", "low vs high soil moisture")
p6 <- f_plot(out_cn, "Recruitment", "low vs high foliar C:N")
p7 <- f_plot(out_ale, "Recruitment", "nonallelopathic vs allelopathic")
p8 <- f_plot(out_mr, "Recruitment", "AM vs EM vs NM")

require(gridExtra)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, ncol = 4, nrow = 2, top = "% change")
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, ncol = 4, nrow = 2)

# surv -----------------------------------------------------------

ms <- readRDS("output/ms_withgroups.RDS")

ms$FIXg <- as.factor(ms$FIX)
ms$CANg <- as.factor(ms$NON0CAN1)
ms$DECg <- ms$Decid0Ever1
ms$DECg[ms$DECg < 0.6] <- 0
ms$DECg <- as.factor(ms$DECg)
ms <- ms %>% mutate(MRg = case_when(myco_phil == "AM" ~ 0, 
                                    myco_phil == "AM+EM" ~ 0,
                                    myco_phil == "EM" ~ 1,
                                    myco_phil == "ERM" ~ 1,
                                    myco_phil == "NM" ~ 2))
ms$MRg <- as.factor(ms$MRg)
ms <- ms %>% mutate(DEPg = case_when(dep < 3.96 ~ 0, 
                                     dep >= 3.96 ~ 1))
ms$DEPg <- as.factor(ms$DEPg)
ms <- ms %>% mutate(SMg = case_when(sm < 0.27 ~ 0, 
                                    sm >= 0.27 ~ 1))
ms$SMg <- as.factor(ms$SMg)
ms <- ms %>% mutate(CNg = case_when(Tcnratio == "low" ~ 0, 
                                    Tcnratio == "high" ~ 1))
ms$CNg <- as.factor(ms$CNg)
ms <- ms %>% mutate(ALEg = case_when(Tallelopathy == "no" ~ 0, 
                                     Tallelopathy == "yes" ~ 1))
ms$ALEg <- as.factor(ms$ALEg)

m <- lmer(survr ~ NCI_props +
            NCIs +
            NCIs^2 +
            dbhs +
            NCIs*NCI_props +
            NCIs*dbhs +
            (1 | pcn1),
          data = ms)

cat("Individual Survival Model:\n", file = "output/table_si3.txt", append = TRUE)
capture.output(summary(m), 
               file = "output/table_si3.txt", append = TRUE)
capture.output(r.squaredGLMM(m), 
               file = "output/table_si3.txt", append = TRUE)
capture.output(confint(m), 
               file = "output/table_si3.txt", append = TRUE)



p1 <- f_plot(out_fix, "Survival", "fixer vs nonfixer")
p2 <- f_plot(out_can, "Survival", "canopy vs noncanopy")
p3 <- f_plot(out_decid, "Survival", "deciduous vs evergreen")
p4 <- f_plot(out_dep, "Survival", "low vs high N deposition")
p5 <- f_plot(out_sm, "Survival", "low vs high soil moisture")
p6 <- f_plot(out_cn, "Survival", "low vs high foliar C:N")
p7 <- f_plot(out_ale, "Survival", "nonallelopathic vs allelopathic")
p8 <- f_plot(out_mr, "Survival", "AM vs EM vs NM")

require(gridExtra)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, ncol = 4, nrow = 2, top = "% change")
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, ncol = 4, nrow = 2)

# check weighted averages work out
total <- sum(out_fix$Freq)
weighted.mean(out_fix$pchc, out_fix$Freq/total) #0.1063442
weighted.mean(out_can$pchc, out_can$Freq/total) #-0.609144
weighted.mean(out_decid$pchc, out_decid$Freq/total) #-0.664473
weighted.mean(out_dep$pchc, out_dep$Freq/total) #-0.4717598
weighted.mean(out_sm$pchc, out_sm$Freq/total) #-0.6317469
weighted.mean(out_cn$pchc, out_cn$Freq/total) #-0.2014377
weighted.mean(out_ale$pchc, out_ale$Freq/total) #-0.6228808
weighted.mean(out_mr$pchc, out_mr$Freq/total) #-0.7159403

weighted.geomean <- function(x, w, ...){
  return(prod(x^w, ...)^(1/sum(w)))
}

weighted.geomean(out_fix$pchc, out_fix$Freq/total) #0.1025078
weighted.geomean(out_can$pchc, out_can$Freq/total) #NaN
weighted.geomean(out_decid$pchc, out_decid$Freq/total) #-NaN
weighted.geomean(out_dep$pchc, out_dep$Freq/total) #-NaN
weighted.geomean(out_sm$pchc, out_sm$Freq/total) #-NaN
weighted.geomean(out_cn$pchc, out_cn$Freq/total) #-NaN
weighted.geomean(out_ale$pchc, out_ale$Freq/total) #-NaN
weighted.geomean(out_mr$pchc, out_mr$Freq/total) #-NaN


# correlations among covariates ----------------------------
# SI Table 4
require(corrplot)

rs_small <- temp1 %>% dplyr::select(spcd, dbh1, NCI_total2, LAT, LON, ELEVm, STDAGE,
                                 texture, LeafRetention,
                                 FIXg, CANg, DECg, MRg, DEPg, SMg,
                                 CNg, AGEg)

rs_small <- mutate_if(rs_small, is.factor, ~ as.numeric(as.character(.x))) # convert factors to numeric
rs_small <- rs_small[complete.cases(rs_small), ]

G <- cor(rs_small)
head(round(G, 2))

corrplot(G, type="upper", tl.col="black")
colnames(G) <- c("Spcies","DBH","NCI","Latitude","Longitude","Elevation","Age","Texture",
                 "Leaf retention", "Fixer",
                 "Canopy","Deciduous","Mycorrhizae","N deposition", "Soil moisture",
                 "Foliar C:N", "Age factor")
rownames(G) <- c("Species","DBH","NCI","Latitude","Longitude","Elevation","Age","Texture",
                 "Leaf retention", "Fixer",
                 "Canopy","Deciduous","Mycorrhizae","N deposition", "Soil moisture",
                 "Foliar C:N", "Age factor")

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(G, method="color", col=col(200),  
         type="upper", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, tl.cex=1.5, #Text label color and rotation
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)


# correlation with NCI or DBH with groups


# model form (compare AICs) --------------------------------------------------
#rs, reci, PCN2
#gs, LGRp, PCN1
#ms, surv, PCN1

# with group interaction
m0 <- lmer(LGRp ~ NCI_props +
                NCIs +
                dbhs +
                NCIs^2 +
                NCI_props*FIXg + 
                (1 | pcn1) - 1, 
              data = gs,
           REML = FALSE)
predict(m0)

m.nci2 <- lmer(LGRp ~ NCI_props +
                NCIs +
                dbhs +
                NCI_props*FIXg + 
                (1 | pcn1) - 1, 
              data = gs,
              REML = FALSE)

m.dbh <- lmer(LGRp ~ NCI_props +
                NCIs +
                NCIs^2 +
                NCI_props*FIXg + 
                (1 | pcn1) - 1, 
              data = gs,
              REML = FALSE)

m.nci <- lmer(LGRp ~ NCI_props +
                dbhs +
                NCIs^2 +
                NCI_props*FIXg + 
                (1 | pcn1) - 1, 
              data = gs,
              REML = FALSE)

m.nci2dbh <- lmer(LGRp ~ NCI_props +
                NCIs +
                NCI_props*FIXg + 
                (1 | pcn1) - 1, 
              data = gs,
              REML = FALSE)

m.nciint <- lmer(LGRp ~ NCI_props +
             NCIs +
             dbhs +
             NCIs^2 +
             NCI_props*FIXg + 
             NCI_props*NCIs + 
             (1 | pcn1) - 1, 
           data = gs,
           REML = FALSE)

m.ncidbh <- lmer(LGRp ~ NCI_props +
                   NCIs +
                   dbhs +
                   NCIs^2 +
                   NCI_props*FIXg + 
                   NCI_props*NCIs + 
                   NCIs*dbhs + 
                   (1 | pcn1) - 1, 
                 data = gs,
                 REML = FALSE)

#RGR ~ NCIprop + NCI + NCI2 + DBH + NCI*DBH 
m.1 <- lmer(LGRp ~ NCI_props +
                   NCIs +
                   dbhs +
                   NCIs^2 +
                   NCI_props*FIXg + 
                   NCIs*dbhs + 
                   (1 | pcn1) - 1, 
                 data = gs,
                 REML = FALSE)

#RGR ~ NCIprop + NCI + NCI2 + DBH
m.2 <- lmer(LGRp ~ NCI_props +
              NCIs +
              dbhs +
              NCIs^2 +
              NCI_props*FIXg + 
              (1 | pcn1) - 1, 
            data = gs,
            REML = FALSE)

#NCIprop + DBH
m.3 <- lmer(LGRp ~ NCI_props +
              dbhs +
              NCI_props*FIXg + 
              (1 | pcn1) - 1, 
            data = gs,
            REML = FALSE)



anova(m0, m.nci2, m.dbh, m.nci, m.nci2dbh, m.nciint, m.ncidbh, m.1, m.2, m.3)
anova(m0, m.nci2, m.dbh, m.nci, m.nci2dbh) #if couldn't allocate mem for last model

# when using ms, AICs:
#Df       AIC       BIC  logLik  deviance      Chisq Chi Df
#m.dbh      7 -10752755 -10752664 5376385 -10752769                  
#m.nci2dbh  7 -10752755 -10752664 5376385 -10752769     0.0000      0
#m0         8 -10771041 -10770937 5385529 -10771057 18288.1324      1
#m.nci2     8 -10771041 -10770937 5385529 -10771057     0.0000      0
#m.nci      8 -10771041 -10770937 5385529 -10771057     0.0000      0
#m.nciint   9 -10771040 -10770923 5385529 -10771058     0.5131      1

# when using gs: AICs:
#Df       AIC       BIC   logLik  deviance     Chisq Chi Df
#m.dbh      7 -20794624 -20794532 10397319 -20794638                 
#m.nci2dbh  7 -20794624 -20794532 10397319 -20794638      0.00      0
#m0         8 -21092741 -21092635 10546378 -21092757 298118.63      1
#m.nci2     8 -21092741 -21092635 10546378 -21092757      0.00      0
#m.nci      8 -21092741 -21092635 10546378 -21092757      0.00      0
#m.nciint   9 -21092776 -21092657 10546397 -21092794     37.15      1

# when using rs: AICS:
# Df       AIC       BIC  logLik  deviance  Chisq Chi Df Pr(>Chisq)
# m.dbh      7 -13730294 -13730200 6865154 -13730308                         
# m.nci2dbh  7 -13730294 -13730200 6865154 -13730308      0      0          1
# m0         8 -14146129 -14146022 7073073 -14146145 415837      1     <2e-16
# m.nci2     8 -14146129 -14146022 7073073 -14146145      0      0          1
# m.nci      8 -14146129 -14146022 7073073 -14146145      0      0          1


# without group interaction
# to get AIC must unload lmerTest package

#rs, reci, PCN2
#gs, LGRp, PCN1
#ms, survi, PCN1

rs.s <- rs %>% dplyr::select(reci, NCI_props, NCIs, dbhs, pcn2)
rs.s <- na.omit(rs.s)
ms.s <- ms %>% dplyr::select(survi, NCI_props, NCIs, dbhs, pcn1)
ms.s <- na.omit(ms.s)

m0 <- lmer(survi ~ NCI_props +
             NCIs +
             dbhs +
             NCIs^2 +
             NCI_props*NCIs +
             NCIs * dbhs + 
             (1 | pcn1) - 1, 
           data = ms.s,
           REML = FALSE)

drop1(m0, test='none')

m1 <- lmer(survi ~ NCI_props +
             NCIs +
             dbhs +
             NCIs^2 +
             (1 | pcn1) - 1, 
           data = ms.s,
           REML = FALSE)

drop1(m1, test='none')

m2 <- lmer(survi ~ NCI_props +
             NCIs +
             dbhs +
             (1 | pcn1) - 1, 
           data = ms.s,
           REML = FALSE)

AIC(m2)
require(MuMIn)
AICc(m2)



# validate models (check assumptions) ----------------------------------

# survival
model <- lmer(surv ~ NCI_props +
                NCIs +
                dbhs +
                NCIs^2 +
                NCI_props*FIXg + 
                (1 | PCN1) - 1, 
              data = ms)

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(model, add.smooth = FALSE, which = 1)
E <- resid(model)
hist(E, xlab = "Residuals", main = "")
plot(ms[!is.na(ms$NCI_props),"NCI_props"], E, xlab = "NCI_props",
     ylab = "Residuals")
plot(ms[!is.na(ms$NCI_props),"FIXg"], E, xlab = "Fix",
     ylab = "Residuals")
par(op)

# growth
model <- lmer(LGRp ~ NCI_props +
                NCIs +
                dbhs +
                NCIs^2 +
                NCI_props*FIXg + 
                (1 | PCN1) - 1, 
              data = gs)
op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(model, add.smooth = FALSE, which = 1)
E <- resid(model)
hist(E, xlab = "Residuals", main = "")
plot(gs[!is.na(gs$NCI_props),"NCI_props"], E, xlab = "NCI_props",
     ylab = "Residuals")
plot(gs[!is.na(gs$NCI_props),"FIXg"], E, xlab = "Fix",
     ylab = "Residuals")
par(op)


# recr
model <- lmer(reci ~ NCI_props +
                NCIs +
                dbhs +
                NCIs^2 +
                NCI_props*FIXg + 
                (1 | PCN2) - 1, 
              data = rs)

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(model, add.smooth = FALSE, which = 1)
E <- resid(model)
hist(E, xlab = "Residuals", main = "")
plot(rs[!is.na(rs$NCI_props),"NCI_props"], E, xlab = "NCI_props",
       ylab = "Residuals")
plot(rs[!is.na(rs$NCI_props),"FIXg"], E, xlab = "Fix",
       ylab = "Residuals")
par(op)



# -----

#saveRDS(inddat, "output/inddata.RDS")
inddat <- readRDS("output/boots/inddat_10bs.RDS")

inddat <- inddat %>% mutate(level = case_when(group == 0 ~ "a", #low
                                              group == 1 ~ "b",
                                              group == 2 ~ "c"))
inddat$type <- ordered(inddat$type, levels = c("RGR", "Recruitment", "Survival"))
inddat$type2 <- factor(inddat$type, labels = c("RGR[NFE]", "r[NFE]", "s[NFE]"))

inddat <- inddat %>% filter(level != "c") #remove non-myco (too few so irrelevant)

inddat$var2 <- factor(c("non-fixer","N fixer","non-canopy","canopy","deciduous","evergreen",
                 "low foliar C:N","high foliar C:N","AM","EM",
                 "non-fixer","N fixer","non-canopy","canopy","deciduous","evergreen",
                 "low foliar C:N","high foliar C:N","AM","EM",
                 "non-fixer","N fixer","non-canopy","canopy","deciduous","evergreen",
                 "low foliar C:N","high foliar C:N","AM","EM"))
inddat$var2 <- factor(inddat$var2, levels = c("non-fixer","N fixer","non-canopy","canopy","deciduous","evergreen",
                                              "low foliar C:N","high foliar C:N","AM","EM"))
inddat$var <- factor(inddat$var, levels=c('fix',"canopy","deciduousness","foliar C:N","mycorrhizae"))

plot <- ggplot(data = inddat, 
               aes(x = factor(var), y = mean, fill = level, width = 0.65)) + #was var2
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
        plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0)),
        axis.line = element_line(size=1, colour = "black"),
        legend.position = "none") 

png("output/fig3_ind.png", width = 12, height = 7, units = 'in', res = 300)
par(bty = 'n')
plot
dev.off()




# calc p values for whether groups are different
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

pvalfunc(ch_fix) #0.00000001069
pvalfunc(ch_can) #p-value < 2.2e-16
pvalfunc(ch_decid) #0.0000004098
pvalfunc(ch_cn) #< 2.2e-16
pvalfunc(ch_myco) # 0.5365 #only comparing AM and EM (not to NM)

pvalfunc(ch_fix_r) #0.5787
pvalfunc(ch_can_r) #0.000007784
pvalfunc(ch_decid_r) #0.9567
pvalfunc(ch_cn_r) #0.002705
pvalfunc(ch_myco_r) #0.2914 (am to em)

pvalfunc(ch_fix_s) #0.1798
pvalfunc(ch_can_s) #0.5235
pvalfunc(ch_decid_s) # 0.000000002556
pvalfunc(ch_cn_s) #0.006005
pvalfunc(ch_myco_s) #0.006678


#check weighted means
wm <- function(data){
  total <- sum(data$Freq)
  w0 <- data$Freq[data$group==0]/total
  w1 <- data$Freq[data$group==1]/total
  out <- data$pchc[data$group==0]*w0 + data$pchc[data$group==1]*w1
  return(out)
}

wm(fix_ms)
wm(can_ms)
wm(dec_ms)
wm(cn_ms)
wm(lr_ms)


# mapping ------------------------------------
# map functions in functions.R

map_gs <- f_map(gs, "growth", "ind")
map_gs <- readRDS("output/map_gs_ind.RDS")
m1 <- f_mapplot(map_gs$pct.grid, "RGR", "individual", "percent")
saveRDS(map_gs, "output/map_gs_ind.RDS")

map_rs <- f_map(rs, "recr", "ind")
saveRDS(map_rs, "output/map_rs_ind.RDS")
m2 <- f_mapplot(map_rs$pct.grid, "recruitment rate", "individual", "percent change")

map_ms <- f_map(ms, "surv", "ind")
saveRDS(map_ms, "output/map_ms_ind.RDS")
m3 <- f_mapplot(map_ms$coef.grid, "survival rate", "individual", "percent change")


# Fig 1
map_gs <- readRDS("output/map_gs_ind.RDS")
map_rs <- readRDS("output/map_rs_ind.RDS")
map_ms <- readRDS("output/map_ms_ind.RDS")
lon.list <- readRDS("output/lon_list.RDS")
lat.list <- readRDS("output/lat_list.RDS")
gs_noout <- as.vector(map_gs$pct.grid)
gs_noout <- gs_noout[gs_noout > -1000 & gs_noout < 1000]
gs_trunc <- as.matrix(map_gs$pct.grid)
gs_trunc[gs_trunc > 100] <- 100
rs_noout <- as.vector(map_rs$pct.grid)
rs_noout <- rs_noout[rs_noout > -1000 & rs_noout < 1000]
rs_trunc <- as.matrix(map_rs$pct.grid)
rs_trunc[rs_trunc > 100] <- 100
ms_noout <- as.vector(map_ms$pct.grid)
ms_noout <- ms_noout[ms_noout > -1000 & ms_noout < 1000]
ms_trunc <- as.matrix(map_ms$pct.grid)
ms_trunc[ms_trunc > 100] <- 100

png("output/fig1_ind.png", width = 12, height = 7, units = 'in', res = 300)
par(mfrow=c(2,3),
    #mai = c(1, 0.1, 0.1, 0.1),
    mai = c(0.6, 0.8, 0.8, 1),
    bty = 'n')
source('scripts/functions.R')
require(fields)
f_mapplot(gs_trunc, "RGR", "individual-scale", "percent") #map_gs$pct.grid if don't want truncated values
f_mapplot(rs_trunc, "r", "individual-scale", "percent")
f_mapplot(ms_trunc, "s", "individual-scale", "percent")
hist(gs_noout, breaks=50, xlab = "% change in RGR", prob=T, xlim=c(-1000,1000),
     main=expression('individual-scale RGR'[NFE]), cex.lab = 2, cex.main = 2, cex.axis = 2,
     yaxt='n')
axis(side = 2, at = c(0,0.009), labels = c(0,0.009), cex.axis=2)
abline(v = -1.5, col = "blue", lwd = 4) #value from table 1
hist(rs_noout, breaks=50, xlab = "% change in recruitment rate", prob=T, xlim=c(-1000,1000),
     main=expression('individual-scale r'[NFE]), cex.lab = 2, cex.main = 2, cex.axis = 2,
     yaxt='n')
axis(side = 2, at = c(0,0.005), labels = c(0,0.005), cex.axis=2)
abline(v = 7.5, col = "blue", lwd = 4)
hist(ms_noout, breaks=50, xlab = "% change in survival rate", prob=T, xlim=c(-1000,1000),
     main=expression('individual-scale s'[NFE]), cex.lab = 2, cex.main = 2, cex.axis = 2,
     yaxt='n')
axis(side = 2, at = c(0,0.019), labels = c(0,0.019), cex.axis=2)
abline(v = -0.9, col = "blue", lwd = 4)
dev.off()


# mantels test for spatial autocorrelation
# https://stats.idre.ucla.edu/r/faq/how-can-i-perform-a-mantel-test-in-r/

lonvec <- rep(lon.list, 78) #78 is length of lat list
latvec <- rep(lat.list, 319) #319 is lenght of lon list
library(ade4) #for mantel test
library(ape) #for morans i test

dat <- map_ms$pct.grid #rotate through map_gs, map_rs, map_ms
datv <- as.vector(dat)

dat <- cbind.data.frame(datv, latvec, lonvec)
colnames(dat) <- c("pctch", "lat", "lon")

dat1 <- dat[complete.cases(dat),]

plot.dists <- as.matrix(dist(cbind(dat1$lon, dat1$lat)))
#data.dists <- dist(dat1$pctch)

plot.dists.inv <- 1/plot.dists
diag(plot.dists.inv) <- 0

#mantel.rtest(plot.dists, data.dists, nrepet = 999) #mantel test (wrong)
Moran.I(dat1$pctch, plot.dists.inv) #moran's I test (right)

#try ecodist or vegan packages
# require(vegan)
# geo <- cbind(latvec, lonvec)
# commdist <- vegdist(datv, method="bray")
# geodist <- vegdist(geo, method="euclidean")
# 
# plot(geodist,commdist,pch=16,cex=0.5,col="black",bty="l",xlab="Spatial distance",ylab="Species composition dissimilarity")

# model testing -------------------------------------------------------------
# notes from meeting with stats consultant:
# add more fixed effects
# transform response variable (log transform, box cox transform)
# randomly sample from full data to play around with it

# get sample to play around with
# needs multiple trees from same plot
# nees plots where NCI_prop is not 0
options(scipen = 6)


# play at first sample 100 plots
# then increase to 1000 plots
samp <- sample(1:length(unique(gs$pcn1)), 1000, replace = F)
samppcns <- unique(gs$pcn1)[samp]
temp <- gs[gs$pcn1 %in% samppcns, ]

temp <- dplyr::rename(temp, GENUS = GENUS.x)

# model improvements: http://docs.statwing.com/interpreting-residual-plots-to-improve-your-regression/

model <- glmer(sqrt(LGRp) ~ NCI_props +
        NCIs +
        dbhs +
        NCIs^2 +
        (1 | pcn1) - 1,
      data = temp,
      family = Gamma(link="log"))

model <- lmer(log(LGRp+0.01) ~ NCI_props +
                 NCIs +
                 dbhs +
                 NCIs^2 +
                 (1 | pcn1) - 1,
               data = temp)


m0 <- lmer(log(LGRp+0.01) ~ NCI_props +
                NCIs +
                dbhs +
                NCIs^2 +
                (1 | pcn1) + (1 | GENUS) - 1,
              data = temp)


m1 <- lmer(log(LGRp+0.01) ~ NCI_props +
             NCIs +
             dbhs +
             (1 | pcn1) + (1 | GENUS) - 1,
           data = temp)

m2 <- lmer(log(LGRp+0.01) ~ NCI_props +
             NCIs +
             dbhs +
             NCI_props*NCIs +
             (1 | pcn1) + (1 | GENUS) - 1,
           data = temp)

m3 <- lmer(LGRp ~ NCI_props +
             NCIs +
             dbhs +
             (1 | pcn1) + (1 | GENUS) - 1,
           data = temp)
plot(m1, add.smooth=F)
plot(m3, add.smooth=F)

m1 <- lmer(log(LGRp+0.01) ~ NCI_props +
             NCIs +
             dbhs +
             NCI_props*FIXg +
             (1 | pcn1) + (1 | GENUS) - 1,
           data = temp)

library(effects)
e <- allEffects(m1)
print(e)
plot(e)

# load package
library(sjPlot)
library(sjmisc)
library(sjlabelled)
tab_model(m1)


# survival

options(scipen = 6)

ms_no0 <- ms[which(ms$NCI_prop2 > 0),]

# play at first sample 100 plots
# then increase to 1000 plots
samp <- sample(1:length(unique(ms$pcn1)), 1000, replace = F)
samppcns <- unique(ms$pcn1)[samp]
temp <- ms[ms$pcn1 %in% samppcns, ]


m1 <- glmer(survi ~ NCI_props +
              NCIs +
              dbhs + 
              (1 | pcn1),
            offset = log(t),
            data = temp,
            family = poisson(link="log"))

m2 <- glmer(survi ~ NCI_props +
              NCIs +
              dbhs + 
              (1 | pcn1) + (1 | GENUS),
            offset = log(t),
            data = temp,
            family = poisson(link="log"))

AIC(m1, m2)


samp <- sample(1:length(unique(rs$pcn1)), 1000, replace = F)
samppcns <- unique(rs$pcn1)[samp]
temp <- rs[rs$pcn1 %in% samppcns, ]
m1 <- glmer(reci ~ NCI_props +
              NCIs +
              dbhs + 
              (1 | pcn1),
            offset = log(t),
            data = temp,
            family = poisson(link="log"))

m2 <- glmer(reci ~ NCI_props +
              NCIs +
              dbhs + 
              (1 | pcn1) + (1 | GENUS.x),
            offset = log(t),
            data = temp,
            family = poisson(link="log"))

m3 <- glmer(reci ~ NCI_props +
              NCIs +
              (1 | pcn1),
            offset = log(t),
            data = temp,
            family = poisson(link="log"))

AIC(m1, m2, m3)





# plot 10/16 --------------------------
f_plot <- function(modelout, type, title){
  plot <- ggplot(data = modelout, 
                 aes(x = factor(X), y = pchc, fill=as.factor(group), width = 0.65)) +
    geom_bar(stat = "identity", 
             position = position_dodge()) + 
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Group", y = "% Change", subtitle = type, title = title,
         caption="Number in 1000s of samples") +
    geom_hline(yintercept = 0, color = "black") +
    geom_errorbar(aes(ymin=pchl, ymax=pchu), width=.2,
                  position=position_dodge(.9)) +
    geom_text(aes(label=round(Freq/1000, 1)), position=position_dodge(0.9), vjust=-0.2) +
    theme(panel.background = element_rect(fill = 'white', colour = 'white'),
          plot.title = element_text(hjust = 0.5, size = 20),
          text = element_text(size = 20),
          axis.text.x = element_text(angle = 90, hjust = 1, size = 20),
          plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0))) 
  return(plot)
}

changes <- read.csv("data/derived/nonpar_results_10_16.csv")
changes$group <- as.factor(changes$group)

c1 <- changes[changes$level=="Ind",]
c1$type <- ordered(c1$type, levels = c("RGR", "Recruitment", "Survival"))

plot <- ggplot(data = c1, 
               aes(x = factor(var), y = pchc, fill = group, width = 0.65)) +
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

c2 <- changes[changes$level=="Plot",]
c2$type <- ordered(c2$type, levels = c("BAI", "Recruitment", "Survival", "BAIn"))

plot <- ggplot(data = c2, 
               aes(x = factor(var), y = pchc, fill = group, width = 0.65)) +
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


# map where evergreens are:
temp <- gs %>% dplyr::filter(Decid0Ever1 == 1)

library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)

usa <- map_data("usa")

ggplot() + 
  geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = NA, color = "red") + 
  coord_fixed(1.3) +
  geom_point(data = temp, aes(x = LON, y = LAT), color = "black")


# which evergreens coexist with N-fixers
#get plots with N fixers
nfixers <- gs[gs$Genus.x != "NA",]
plotswithnfix <- unique(nfixers$pcn1) #2470 plots with fixers
plotsfix <- gs[gs$pcn1 %in% plotswithnfix,]
everfix <- plotsfix %>% dplyr::filter(Decid0Ever1 == 1)
unique(everfix$Genus.y)
xtabs(~Genus.y , data=everfix) #numer of neighbors to nfixers that are each genus of evergreen
fix <- data.frame(xtabs(~Genus.y , data=everfix))
fix <- fix[fix$Freq > 0,]


#Alder
alders <- plotsfix %>% dplyr::filter(Genus.x == "Alnus")
plotswithalder <- unique(alders$pcn1)
plotsalder <- gs[gs$pcn1 %in% plotswithalder,]
everalder <- plotsalder %>% dplyr::filter(Decid0Ever1 == 1)
unique(everalder$Genus.y)
alder <- data.frame(xtabs(~Genus.y, data=everalder))
alder <- alder[alder$Freq > 0, ]



### what fraction of N-fixers are canopy trees
temp <- gs %>% filter(FIX == 1)
xtabs(~CANg, data=temp)

temp <- gs %>% filter(FIX == 1 & CANg == 0)
xtabs(~Genus.x, data=temp) #1355 Robinia

hist(temp$STDAGE, col=rgb(1,0,0,0.5), breaks=25, prob=T) #understory N-fixers
nfixers <- gs %>% filter(FIX == 1)
hist(nfixers$STDAGE, col=rgb(0,0,1,0.5), breaks=25, prob=T, add=T) #all N-fixers


### what is the ineteraction between C:N status and evergreenness
xtabs(~CNg + DECg, gs)

temp <- gs[gs$CNg == 0 & gs$DECg == 1, ]
unique(temp$gensp)
temp2 <- xtabs(~gensp, temp)
sort(temp2, decreasing = TRUE)

xtabs(~gensp + AGEg, temp)


##########################################
################ GEOGRAPHY ##################
##########################################

