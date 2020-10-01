require(data.table)
require(fields)

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

# ++++++++++++++++++++++++++++
# f_pctch
# ++++++++++++++++++++++++++++
# model : name of model run for model_N_group or model_non_group
# old orig: orig <- stats_non_NCI_prop_noncanopy$stats["mean",][[1]]+fixef(model)[2]*0
#standard error calc from: https://www.ncbi.nlm.nih.gov/books/NBK33484/
f_pctch <- function(model) {
  orig <- fixef(model)[[1]]+fixef(model)[[2]]*0
  final <- fixef(model)[[1]]+fixef(model)[[2]]*1
  ch <- 100*(final-orig)/orig
  ch
  var_orig <- (coef(summary(model))[1,"Std. Error"])^2*nobs(model)
  var_final <- (coef(summary(model))[1,"Std. Error"])^2*nobs(model) +
    (coef(summary(model))["NCI_props","Std. Error"])^2*nobs(model)
  var_ch <- (final/orig)^2*(var_orig/orig^2+var_final/final^2)
  se_ch <- sqrt(var_ch)/sqrt(nobs(model))
  se_ch
  num <- nobs(model)
  num
  return(c("est"=ch, "se"=se_ch, "num"=num))
}


# ++++++++++++++++++++++++++++
# f_pctch_g
# ++++++++++++++++++++++++++++
# model : name of model run with group as RE
# old orig: orig <- stats_non_NCI_prop_noncanopy$stats["mean",][[1]]+fixef(model)[2]*0
#standard error calc from: https://www.ncbi.nlm.nih.gov/books/NBK33484/
f_pctch_g <- function(model, nums) {
  orig <- final <- ch <- NULL
  #untransform coefficients
  #mu <- mean(gs$NCI_prop, na.rm=T)
  #sdev <- sd(gs$NCI_prop, na.rm=T)
  mu <- 0.008925627 #same for whole data set
  sdev <- 0.08225147
  b0 <- fixef(model)[["(Intercept)"]] #fixed intercept
  b1 <- fixef(model)[["NCI_props"]] #fixed NCIprop slope
  bg <- ranef(model)$group["NCI_props"] #random NCIprop slope
  b8 <- fixef(model)[["NCI_props:NCIs"]]
  ints <- b0 - (b1 + bg) * (mu/sdev)
  slopes <- (b1 + bg + b8) / sdev
  varc <- summary(model)$varcor$FIX
  sigma <- as.numeric(attr(varc, "stddev")) #SD for random effect
  grps <- rownames(ranef(model)$group)
  re <- cbind.data.frame(ints, slopes, grps, sigma)
  colnames(re) <- c("int", "slope", "grp", "sd")
  re <- merge(re, nobs, by.x="grp", by.y="grp", all.x=T)
  re$se <- re$sd/(re$n)^0.5
  for(i in 1:nrow(re)){
    orig[i] <- re[i,"int"]
    final[i] <- re[i,"int"]+re[i,"slope"]
    ch[i] <- 100*(final[i]-orig[i])/orig[i]
  }
  ch <- cbind.data.frame(ch, grps)
  # var_orig <- (coef(summary(model))[1,"Std. Error"])^2*nobs(model)
  # var_final <- (coef(summary(model))[1,"Std. Error"])^2*nobs(model) +
  #   (coef(summary(model))["NCI_props","Std. Error"])^2*nobs(model)
  # var_ch <- (final/orig)^2*(var_orig/orig^2+var_final/final^2)
  # se_ch <- sqrt(var_ch)/sqrt(nobs(model))
  # se_ch
  num <- nobs(model)
  num
  #return(c("est"=ch, "se"=se_ch, "num"=num))
  return(list("est"=ch, "num"=num))
}



# ++++++++++++++++++++++++++++
# f_pctch_pg
# ++++++++++++++++++++++++++++
# model : name of model run with group as RE
# old orig: orig <- stats_non_NCI_prop_noncanopy$stats["mean",][[1]]+fixef(model)[2]*0
#standard error calc from: https://www.ncbi.nlm.nih.gov/books/NBK33484/
f_pctch_pg <- function(model) {
  orig <- final <- ch <- NULL
  #untransform coefficients
  #mu <- mean(FIA_plot$BA_Nfixer_prop_total, na.rm=T)
  #sdev <- sd(FIA_plot$BA_Nfixer_prop_total, na.rm=T)
  #mu <- 0.01361491 #same for whole data set
  #sdev <- 0.08622454
  b0 <- fixef(model)[["(Intercept)"]] #fixed intercept
  b1 <- fixef(model)[["BA_Nfixer_prop_total"]] #fixed NCIprop slope
  bg <- ranef(model)$group["BA_Nfixer_prop_total"] #random NCIprop slope
  #ints <- b0 - (b1 + bg) * (mu/sdev)
  ints <- b0
  #slopes <- (b1 + bg) / sdev
  slopes <- b1 + bg
  grps <- rownames(ranef(model)$group)
  re <- cbind.data.frame(ints, slopes, grps)
  colnames(re) <- c("int", "slope", "grp")
  for(i in 1:nrow(re)){
    orig[i] <- re[i,"int"]
    final[i] <- re[i,"int"]+re[i,"slope"]
    ch[i] <- 100*(final[i]-orig[i])/orig[i]
  }
  ch <- cbind.data.frame(ch, grps)
  # var_orig <- (coef(summary(model))[1,"Std. Error"])^2*nobs(model)
  # var_final <- (coef(summary(model))[1,"Std. Error"])^2*nobs(model) +
  #   (coef(summary(model))["NCI_props","Std. Error"])^2*nobs(model)
  # var_ch <- (final/orig)^2*(var_orig/orig^2+var_final/final^2)
  # se_ch <- sqrt(var_ch)/sqrt(nobs(model))
  # se_ch
  num <- nobs(model)
  num
  #return(c("est"=ch, "se"=se_ch, "num"=num))
  return(list("est"=ch, "num"=num))
}


# ++++++++++++++++++++++++++++
# f_pctch_p
# ++++++++++++++++++++++++++++
# model : name of model run for model_N_group or model_non_group
# old orig: orig <- stats_non_NCI_prop_noncanopy$stats["mean",][[1]]+fixef(model)[2]*0
#standard error calc from: https://www.ncbi.nlm.nih.gov/books/NBK33484/
f_pctch_p <- function(model) {
  orig <- fixef(model)[[1]]+fixef(model)[[2]]*0
  final <- fixef(model)[[1]]+fixef(model)[[2]]*1
  ch <- 100*(final-orig)/orig
  ch
  var_orig <- (coef(summary(model))[1,"Std. Error"])^2*nobs(model)
  var_final <- (coef(summary(model))[1,"Std. Error"])^2*nobs(model) +
    (coef(summary(model))["scale(BA_Nfixer_prop_total)","Std. Error"])^2*nobs(model)
  var_ch <- (final/orig)^2*(var_orig/orig^2+var_final/final^2)
  se_ch <- sqrt(var_ch)/sqrt(nobs(model))
  se_ch
  num <- nobs(model)
  num
  return(c("est"=ch, "se"=se_ch, "num"=num))
}


# ++++++++++++++++++++++++++++
# f_recr
# ++++++++++++++++++++++++++++
# model : name of model for plot level recruitment
f_recr <- function(temp) {
  model_all <- lmer(recyr ~ scale(BA_Nfixer_prop_total) + 
                      scale(BA_total1) + 
                      scale(BA_total2) + 
                      scale(BA_Nfixer1) + 
                      scale(BA_nonfixer1) +
                      (1 | STATECD), 
                    data = temp)
  return(model_all)
}

# ++++++++++++++++++++++++++++
# f_mort
# ++++++++++++++++++++++++++++
# model : name of model for plot level mortality
f_mort <- function(temp) {
  model_all <- lmer(mortyr ~ scale(BA_Nfixer_prop_total) + 
                      scale(BA_total1) + 
                      scale(BA_total2) + 
                      scale(BA_Nfixer2) + 
                      scale(BA_nonfixer2) +
                      (1 | STATECD), 
                    data = temp) 
  return(model_all)
}

# ++++++++++++++++++++++++++++
# f_surv
# ++++++++++++++++++++++++++++
# model : name of model for plot level survival
f_surv <- function(temp) {
  model_all <- lmer(survyr ~ scale(BA_Nfixer_prop_total) + 
                      scale(BA_total1) + 
                      scale(BA_total2) + 
                      scale(BA_Nfixer2) + 
                      scale(BA_nonfixer2) +
                      (1 | STATECD), 
                    data = temp) 
  return(model_all)
}

# ++++++++++++++++++++++++++++
# f_growth
# ++++++++++++++++++++++++++++
# model : name of model for plot level growth rate
f_growth <- function(temp) {
  model_all <- lmer(BA_change ~ scale(BA_Nfixer_prop_total) + 
                      scale(BA_total1) + 
                      scale(BA_total2) + 
                      scale(BA_Nfixer1) + 
                      scale(BA_nonfixer1) +
                      (1 | STATECD), 
                    data = temp)
  return(model_all)
}

# ++++++++++++++++++++++++++++
# f_grmodp
# ++++++++++++++++++++++++++++
# model : name of model for individual level analysis with grouping as RE
f_grmodp <- function(datapass, type, group) {
  if(type == "growth"){
    model <- lmer(BA_change ~ BA_Nfixer_prop_total + 
                    BA_total1 + 
                    BA_total2 + 
                    BA_Nfixer1 + 
                    BA_nonfixer1 +
                    (1 | STATECD) + 
                    (0 + BA_Nfixer_prop_total | group), 
                  data = datapass,
                  na.action = na.omit)
  }else if(type == "recr" | type == "recruitment"){
    model <- lmer(recyr ~ BA_Nfixer_prop_total + 
                    BA_total1 + 
                    BA_total2 + 
                    BA_Nfixer1 + 
                    BA_nonfixer1 +
                    (1 | STATECD) + 
                    (0 + BA_Nfixer_prop_total | group), 
                  data = datapass,
                  na.action = na.omit)
  }else if(type == "mort" | type == "mortality" | type == "surv"){
    model <- lmer(survyr ~ BA_Nfixer_prop_total + 
                    BA_total1 + 
                    BA_total2 + 
                    BA_Nfixer2 + 
                    BA_nonfixer2 +
                    (1 | STATECD) + 
                    (0 + BA_Nfixer_prop_total | group), 
                  data = datapass,
                  na.action = na.omit) 
  }
  
  return(model)
}



# ++++++++++++++++++++++++++++
# f_group
# ++++++++++++++++++++++++++++
# model : name of model for individual level analysis with grouping as RE
f_grmod <- function(datapass, type, group) {
  if(type == "growth"){
    model <- lmer(LGRp ~ NCI_props +
                    NCIs +
                    dbhs +
                    NCIs^2 +
                    NCI_props * NCIs +
                    dbhs * NCIs +
                    (1 | PCN1) + 
                    (0 + NCI_props | group), 
                  data = datapass)
  }else if(type == "recr" | type == "recruitment"){
    model <- lmer(reci ~ NCI_props + 
                    NCIs + 
                    NCIs^2 + 
                    NCI_props*NCIs +
                    (1|PCN2) +
                    (0 + NCI_props|group),  
                  data = datapass,
                  na.action = na.omit)
  }else if(type == "survival" | type == "surv"){
    model <- lmer(surv ~ NCI_props + 
                    NCIs + 
                    NCIs^2 + 
                    NCI_props*NCIs +
                    (1|PCN1) +
                    (0 + NCI_props|group), 
                  data = datapass)
  }
  
  return(model)
}



# ++++++++++++++++++++++++++++
# f_plotch
# ++++++++++++++++++++++++++++
# model : name of model for individual level mortality
f_plotch <- function(chout, type, group, xsz) {
  plot <- ggplot(data = chout$est, 
                      aes(x = factor(grps), y = ch, width = 0.65)) +
    geom_bar(stat = "identity", 
             position = position_dodge()) + 
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Group", y = "% Change", subtitle = type, title = group) +
    geom_hline(yintercept = 0, color = "black") +
    theme(panel.background = element_rect(fill = 'white', colour = 'white'),
          plot.title = element_text(hjust = 0.5, size = 20),
          text = element_text(size = 20),
          axis.text.x = element_text(angle = 90, hjust = 1, size = xsz))
  return(plot)
}



#############################
ggCaterpillar <- function(re,QQ=TRUE,likeDotplot=TRUE) {
  require(ggplot2)
  f <- function(x) {
    pv   <- attr(x, "postVar")
    cols <- 1:(dim(pv)[1])
    se   <- unlist(lapply(cols, function(i) sqrt(pv[i, i, ])))
    ord  <- unlist(lapply(x, order)) + rep((0:(ncol(x) - 1)) * nrow(x), each=nrow(x))
    pDf  <- data.frame(y=unlist(x)[ord],
                       ci=1.96*se[ord],
                       nQQ=rep(qnorm(ppoints(nrow(x))), ncol(x)),
                       ID=factor(rep(rownames(x), ncol(x))[ord], levels=rownames(x)[ord]),
                       ind=gl(ncol(x), nrow(x), labels=names(x)))
    
    if(QQ) {  ## normal QQ-plot
      p <- ggplot(pDf, aes(nQQ, y))
      p <- p + facet_wrap(~ ind, scales="free")
      p <- p + xlab("Standard normal quantiles") + ylab("Random effect quantiles")
    } else {  ## caterpillar dotplot
      p <- ggplot(pDf, aes(ID, y)) + coord_flip()
      if(likeDotplot) {  ## imitate dotplot() -> same scales for random effects
        p <- p + facet_wrap(~ ind)
      } else {           ## different scales for random effects
        p <- p + facet_grid(ind ~ ., scales="free_y")
      }
      p <- p + xlab("Levels") + ylab("Random effects")
    }
    
    p <- p + theme(legend.position="none")
    p <- p + geom_hline(yintercept=0)
    p <- p + geom_errorbar(aes(ymin=y-ci, ymax=y+ci), width=0, colour="black")
    p <- p + geom_point(aes(size=1.2), colour="blue") 
    return(p)
  }
  
  lapply(re, f)
}


# ++++++++++++++++++++++++++++
# f_map
# ++++++++++++++++++++++++++++
# f_map : make grid and map
# d = data (gs_withgroups, rs, ms, FIA_plot), 
# type = growth, recr, surv, level = ind, plot
# 
f_map <- function(d, type, level) {
  p <- fread("/Users/Anika/Documents/GradSchool/FIA_CompetitionModel/data/raw/PLOT.csv", 
             select = c("CN", "LAT", "LON"))
  
  p$llu <- p$LAT/100 + p$LON*100
  colnames(p) <- c("PCN","lat1","lon1","llu")
  p$PCN <- as.numeric(p$PCN)
  pu <- p[!duplicated(p$llu), ]
  
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
  
  #sig.grid: alpha values for significance
  #coef.grid: coefficient for NCI_props
  #num.grid: number of N-fixers present
  #nump.grid: number of plots
  #riz.grid: average of ACT1vsRIZ0
  #gen.grid: most common genus
  #se.grid: standard error of regression for NCI_props
  #pct.grid: ratio of coefficient for NCI_props to intercept
  
  coef.grid <- num.grid <- sig.grid <- nump.grid <- riz.grid <- water.grid <- se.grid <- gen.grid <- pct.grid <- array(dim=c(nlon,nlat))
  
  if(level == "individual" | level == "ind"){
    #d[is.na(d$NCI_props), ]$NCI_props <- 0
    
    mu <- 0.008925627 #same for (gs, ms, rs)
    sdev <- 0.08225147
    NCIpropsmin <- -0.08077
    NCIpropsmax <- 17.85478
    
    for(i in 1:nlon){
      print(i/nlon*100)
      lon.a <- lon.list[i] - dlon/2
      lon.b <- lon.list[i] + dlon/2
      for(j in 1:nlat){
        lat.a <- lat.list[j] - dlat/2
        lat.b <- lat.list[j] + dlat/2
        
        #plots <- p[(p$lon1>=lon.a & p$lon1<lon.b & p$lat1>=lat.a & p$lat<lat.b),]$PCN
          if(type == "growth"){
            temp <- d[(d$LON>=lon.a & d$LON<lon.b & d$LAT>=lat.a & d$LAT<lat.b),
                         c("NCI_props","NCIs","dbhs","LGRp","FIX","ACT1vsRIZ0",
                           "spcd","pcn1","GENUS.x")]
            names(temp)[names(temp) == 'GENUS'] <- 'GENUS.x'
          }else if(type == "recr" | type == "recruitment"){
            temp <- d[(d$LON>=lon.a & d$LON<lon.b & d$LAT>=lat.a & d$LAT<lat.b),
                                 c("NCI_props","NCIs","dbhs","recir","FIX","ACT1vsRIZ0","spcd","pcn2","GENUS.x")]
            names(temp)[names(temp) == 'GENUS'] <- 'GENUS.x'
          }else if(type == "survival" | type == "surv"){
            temp <- d[which(d$LON>=lon.a & d$LON<lon.b & d$LAT>=lat.a & d$LAT<lat.b),
                         c("NCI_props","NCIs","dbhs","survr","FIX","ACT1vsRIZ0","spcd","pcn1","GENUS")]
            names(temp)[names(temp) == 'GENUS'] <- 'GENUS.x'
          } #end type choice
          if(nrow(temp)>0){
            if(length(unique(temp$pcn1))>5 | length(unique(temp$pcn2))>5 & length(unique(temp$NCI_props))>3){
              if(type == "growth"){
                model <- lmer(LGRp ~ NCI_props +
                                NCIs +
                                dbhs +
                                (1 | pcn1), 
                              data = temp)
                sd_nciprop <- sd(d$NCI_prop,na.rm=T)
                mean_nciprop <- mean(d$NCI_prop,na.rm=T)
              }else if(type == "recr" | type == "recruitment"){
                model <- lmer(recir ~ NCI_props + 
                                NCIs + 
                                dbhs +
                                (1 | pcn2),  
                              data = temp,
                              na.action = na.omit)
                sd_nciprop <- sd(d$NCI_prop2,na.rm=T)
                mean_nciprop <- mean(d$NCI_prop2,na.rm=T)
              }else if(type == "survival" | type == "surv"){
                model <- lmer(survr ~ NCI_props + 
                                NCIs + 
                                dbhs +
                                (1 | pcn1), 
                              data = temp)
                sd_nciprop <- sd(d$NCI_prop1,na.rm=T)
                mean_nciprop <- mean(d$NCI_prop1,na.rm=T)
              } #end type choice
              coefs <- data.frame(coef(summary(model)))
              coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
              #ints <- coefs["(Intercept)","Estimate"] - (coefs["NCI_props","Estimate"]) * (mu/sdev)
              #slopes <- (coefs["NCI_props","Estimate"] + coefs["NCI_props:NCIs","Estimate"]) / sdev
              ints <- coefs["(Intercept)","Estimate"]
              slope <- coefs["NCI_props","Estimate"]
              #pch <- 100*((ints + slopes) - ints)/ints
              init <- ints + NCIpropsmin*slope
              final <- ints + NCIpropsmax*slope
              pch <- 100*(final - init)/init
              
              sig.grid[i,j] <- coefs["NCI_props","p.z"]
              coef.grid[i,j] <- coefs["NCI_props","Estimate"]*sd_nciprop+mean_nciprop
              num.grid[i,j] <- sum(temp$FIX, na.rm=T)
              nump.grid[i,j] <- length(unique(temp$pcn1))
              riz.grid[i,j] <- mean(temp$ACT1vsRIZ0, na.rm=T)
              gen.grid[i,j] <- toString(tail(names(sort(table(temp$GENUS.x))), 1))
              se.grid[i,j] <- sqrt(diag(vcov(model)))[[2]]
              pct.grid[i,j] <- pch
              #pct.grid[i,j] <- NA
              
            }else{ #end length unique
              sig.grid[i,j] <- NA
              coef.grid[i,j] <- NA
              num.grid[i,j] <- NA
              riz.grid[i,j] <- NA
              se.grid[i,j] <- NA
              gen.grid[i,j] <- NA
              nump.grid[i,j] <- NA
              pct.grid[i,j] <- NA
            } #end else
            
          # }else{ #end nrow temp
          #   sig.grid[i,j] <- NA
          #   coef.grid[i,j] <- NA
          #   num.grid[i,j] <- NA
          #   riz.grid[i,j] <- NA
          #   se.grid[i,j] <- NA
          #   gen.grid[i,j] <- NA
          #   nump.grid[i,j] <- NA
          #   pct.grid[i,j] <- NA
         # } #end else
        }
          
        } #end nlat
      } #end nlon 
   # } #end nlon 
    
  }
  return(list("sig.grid" = sig.grid,
         "coef.grid" = coef.grid, 
         "num.grid" = num.grid, 
         "riz.grid" = riz.grid, 
         "se.grid" = se.grid,
         "gen.grid" = gen.grid,
         "nump.grid" = nump.grid,
         "pct.grid" = pct.grid
         ))
  
} #end function






f_map_plot <- function(d, type){
  p <- fread("/Users/Anika/Documents/GradSchool/FIA_CompetitionModel/data/raw/PLOT.csv", 
             select = c("CN", "LAT", "LON"))
  
  p$llu <- p$LAT/100 + p$LON*100
  colnames(p) <- c("PCN","lat1","lon1","llu")
  p$PCN <- as.numeric(p$PCN)
  pu <- p[!duplicated(p$llu), ]
  
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
  
  #sig.grid: alpha values for significance
  #coef.grid: coefficient for NCI_props
  #num.grid: number of N-fixers present
  #nump.grid: number of plots
  #riz.grid: average of ACT1vsRIZ0
  #gen.grid: most common genus
  #se.grid: standard error of regression for NCI_props
  #pct.grid: ratio of coefficient for NCI_props to intercept
  
  coef.grid <- num.grid <- sig.grid <- nump.grid <- riz.grid <- water.grid <- se.grid <- gen.grid <- pct.grid <- array(dim=c(nlon,nlat))
  
    mu <- 0.01361491 #same for FIA_plot
    sdev <- 0.09404602
    BApropsmin <- -0.1179
    BApropsmax <- 21.3829
    
    for(i in 1:nlon){
      print(i/nlon*100)
      lon.a <- lon.list[i] - dlon/2
      lon.b <- lon.list[i] + dlon/2
      for(j in 1:nlat){
        lat.a <- lat.list[j] - dlat/2
        lat.b <- lat.list[j] + dlat/2
        
        if(type == "growth"){
          temp <- d[(d$LON>=lon.a & d$LON<lon.b & d$LAT>=lat.a & d$LAT<lat.b),
                    c("BA_change","BA_total1","BA_Nfixer_prop_total","pcn1")]
        }else if(type == "recr" | type == "recruitment"){
          temp <- d[(d$LON>=lon.a & d$LON<lon.b & d$LAT>=lat.a & d$LAT<lat.b),
                    c("recyr","BA_total1","BA_Nfixer_prop_total","pcn1")]
        }else if(type == "surv" | type == "survival"){
          temp <- d[(d$LON>=lon.a & d$LON<lon.b & d$LAT>=lat.a & d$LAT<lat.b),
                    c("survyr","BA_total1","BA_Nfixer_prop_total","pcn1")]
        }else if(type == "bain"){
          temp <- d[(d$LON>=lon.a & d$LON<lon.b & d$LAT>=lat.a & d$LAT<lat.b),
                    c("BA_change_nonfixer","BA_total1","BA_Nfixer_prop_total","pcn1")]
        } #end type choice
        
        if(nrow(temp)>0){
          if(length(unique(temp$pcn1))>5){
            if(type == "growth"){
              model <- lm(BA_change ~ BA_Nfixer_prop_total + 
                            BA_total1, 
                          data = temp,
                          na.action = na.omit)
            }else if(type == "recr" | type == "recruitment"){
              model <- lm(recyr ~ BA_Nfixer_prop_total + 
                            BA_total1, 
                          data = temp,
                          na.action = na.omit)
            }else if(type == "mort" | type == "mortality" | type == "surv" | type == "survival"){
              model <- lm(survyr ~ BA_Nfixer_prop_total + 
                              BA_total1, 
                            data = temp,
                            na.action = na.omit) 
            }else if(type == "bain"){
              model <- lm(BA_change_nonfixer ~ BA_Nfixer_prop_total + 
                            BA_total1, 
                          data = temp,
                          na.action = na.omit)
            } #end type mort
            coefs <- data.frame(coef(summary(model)))
            coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
            #ints <- coefs["(Intercept)","Estimate"] - (coefs["BA_Nfixer_prop_total","Estimate"]) * (mu/sdev)
            #slopes <- (coefs["BA_Nfixer_prop_total","Estimate"]) / sdev
            #init <- coefs["BA_Nfixer_prop_total","Estimate"]*0 + coefs["(Intercept)","Estimate"]
            #final <- coefs["BA_Nfixer_prop_total","Estimate"]*1 + coefs["(Intercept)","Estimate"]
            new_d <- data.frame(BA_Nfixer_prop_total=c(min(temp$BA_Nfixer_prop_total,na.rm=T), max(temp$BA_Nfixer_prop_total,na.rm=T)),
                                BA_total1=c(mean(temp$BA_total1,na.rm=T), mean(temp$BA_total1,na.rm=T)),
                                BA_total2=c(mean(temp$BA_total2,na.rm=T), mean(temp$BA_total2,na.rm=T)),
                                BA_Nfixer1=c(mean(temp$BA_Nfixer1,na.rm=T), mean(temp$BA_Nfixer1,na.rm=T)))
            new_d$Pred <- predict(model,new_d,re.form=~0)
            init <- new_d$Pred[1]
            final <- new_d$Pred[2]
            pch <- 100*(init - final)/init
            
            sig.grid[i,j] <- coefs["BA_Nfixer_prop_total","p.z"]
            coef.grid[i,j] <- coefs["BA_Nfixer_prop_total","Estimate"]
            nump.grid[i,j] <- length(unique(temp$pcn1))
            se.grid[i,j] <- sqrt(diag(vcov(model)))[[2]]
            pct.grid[i,j] <- pch
            #pct.grid[i,j] <- NA
            num.grid[i,j] <- NA #always na for plot level
            riz.grid[i,j] <- NA #always na for plot level
            gen.grid[i,j] <- NA #always na for plot level
            
          }else{ #end enough plots
            coef.grid[i,j] <- NA
            se.grid[i,j] <- NA
            nump.grid[i,j] <- NA
            pct.grid[i,j] <- NA
            sig.grid[i,j] <- NA
            num.grid[i,j] <- NA #always na for plot level
            riz.grid[i,j] <- NA #always na for plot level
            gen.grid[i,j] <- NA #always na for plot level
          }
          
        }else{ #end enough rows
          coef.grid[i,j] <- NA
          se.grid[i,j] <- NA
          nump.grid[i,j] <- NA
          pct.grid[i,j] <- NA
          sig.grid[i,j] <- NA
          num.grid[i,j] <- NA #always na for plot level
          riz.grid[i,j] <- NA #always na for plot level
          gen.grid[i,j] <- NA #always na for plot level
        } #end else
        
      } #end nlon
    } #end nlat
    
    return(list("sig.grid" = sig.grid,
                "coef.grid" = coef.grid, 
                "num.grid" = num.grid, 
                "riz.grid" = riz.grid, 
                "se.grid" = se.grid,
                "gen.grid" = gen.grid,
                "nump.grid" = nump.grid,
                "pct.grid" = pct.grid
    ))
}
  



# ++++++++++++++++++++++++++++
# f_mapplot
# ++++++++++++++++++++++++++++
# model : name of model for individual level mortality
# grid is the output grid from f_map
# type is recr, surv, growth
# level is individual or plot
# class is type of grid (significant grid, coefficients, % change)
f_mapplot <- function(grid, type, level, class) {
  grid[grid > 1000] <- NA #remove outliers
  grid[grid < -1000] <- NA #remove outliers
  nHalf = sum(!is.na(grid))/2
  Min = min(grid, na.rm=T)
  Max = max(grid, na.rm=T)
  Thresh = 0
  rc2 = colorRampPalette(colors = c("#E0E0E0", "red4"), space="Lab")(nHalf)    # #FDE2DF
  rc1 = colorRampPalette(colors = c("blue4","#E0E0E0"), space="Lab")(nHalf) # #E0DFFD
  rampcols = c(rc1, rc2)
  #rampcols[c(nHalf, nHalf+1)] = rgb(t(col2rgb("grey")), maxColorValue=256) 
  o <- max(c(abs(Min), abs(Max)))
  u <- abs(o)
  #u <- 1000
  u <- 100
  l <- -abs(o)
  #l <- -900
  l <- -100
  rb1 = seq(l, Thresh, length.out=nHalf+1)
  rb2 = seq(Thresh, u, length.out=nHalf+1)[-1]
  rampbreaks = c(rb1, rb2)
  
  nlev <- 64
  
  #o <- ceiling(o)
  #o <- ceiling(o*10)/10
  n <- nchar(round(o))
  t <- 10^(n-1)
  o <- ceiling(o/t)*t
  #inc <- o/5
  #inc <- o/t
  inc <- 300
  #ticks <- seq(-o, o, inc)
  #ticks <- c(ticks, 0)
  #ticks <- seq(-o, o, by=inc) #best if not setting
  #ticks <- seq(-900, 1000, by = 300)
  ticks <- c(-100, -50, 0, 50, 100)

  plot <- image.plot(lon.list, lat.list, grid, 
                       nlevel = nlev, 
                       col = rampcols, 
                       breaks = rampbreaks,
                       axis.args = list(at=ticks,labels=ticks,cex.axis=2),
                       xlab = "longitude",
                       ylab = "latitude",
                       ylim = c(23, 53), xlim = c(-130,-65),
                     cex.axis=2,cex.lab=2,legend.cex=1.8)
  #title(main = paste("map of ", class, " for ", level, " level ", type), cex.main=1.5)
  #title(main = paste(level, type), cex.main=2)
  towrite = paste(level, type)
  title(main = bquote(.(towrite)[NFE]), cex.main=2)
  maps::map('world', regions='usa', add=TRUE)
  return(plot)
}


## individual functions
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
                    dbhs +
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
                    dbhs +
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

#test <- f_pch_nonpar(data=gs, type="growth", nboot=3, group="FIXg")



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


library(BSDA)

pvalfunc <- function(data){
  sd0 <- sqrt(data$Freq[data$group==0])*(data$upper[data$group==0]-data$lower[data$group==0])/3.92
  sd1 <- sqrt(data$Freq[data$group==1])*(data$upper[data$group==1]-data$lower[data$group==1])/3.92
  
  out <- tsum.test(mean.x = data$mean[data$group==0], s.x = sd0, n.x = data$Freq[data$group==0],
                   mean.y = data$mean[data$group==1], s.y = sd1, n.y = data$Freq[data$group==1])
  return(out)
  
}




#### Plot level functions
# functions -----------------------------------------------------------

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
        # BAI ~ BApct + BAtotal,1 + BAtotal,2 + BANfixer,1
        m <- lmer(LGRp ~ NCI_props + 
                    NCIs +
                    dbhs +
                    (1 | pcn1) - 1,
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
        # Recruitment ~ BApct + BAtotal,1 + BAtotal,2
        m <- lmer(recir ~ NCI_props +
                    NCIs +
                    dbhs +
                    (1 | pcn1) - 1,
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
        # Survival ~ BApct + BAtotal,1 + BAtotal,2 + BANfixer,1
        m <- lmer(survr ~ NCI_props +
                    NCIs +
                    dbhs +
                    (1 | pcn1) - 1,
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
                    dbhs +
                    NCI_props*group +
                    (1 | pcn1) - 1,
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
                    dbhs +
                    NCI_props*group +
                    (1 | pcn1) - 1,
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
                    dbhs +
                    NCI_props*group +
                    (1 | pcn1) - 1,
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
  
  if(ngroup==3){
    result <- cbind.data.frame(pred.y.0, pred.y.1, pred.y.2, r2)
    colnames(result) <- c("group0low", "group0high", "group1low", "group1high", 
                          "group2low", "group2high", "R2m", "R2c")
  }else{
    result <- cbind.data.frame(pred.y.0, pred.y.1, r2)
    colnames(result) <- c("group0low", "group0high", "group1low", "group1high", "R2m", "R2c")
  }
  
  return(result)
  
}#end function

#test <- f_pch_nonpar(data=gs, type="growth", nboot=3, group="FIXg")



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




# ++++++++++++++++++++++++++++
# f_nonpar_boot
# ++++++++++++++++++++++++++++
require(MuMIn)
# model : name of model for individual level analysis with grouping as RE

f_nonpar_boot <- function(data, type, group, nboot) {
  pred.y.0 <- matrix(ncol=2,nrow=nboot)
  pred.y.1 <- matrix(ncol=2,nrow=nboot)
  r2 <- matrix(ncol=2,nrow=nboot)
  
  for(i in 1:nboot){
    progress(i, max.value=nboot)
    
    # get random sample
    numplot <- nrow(data)
    inds <- sample(1:numplot, numplot, replace=T)
    dat <- data[inds,]
    
    if(group == "NA"){
      keep <- c("STATECD","BA_change","BA_change_Nfixer","recyr","survyr","BA_change_nonfixer","BA_Nfixer_prop_total","BA_total1","BA_total2","BA_Nfixer1")
      dat <- dat %>% dplyr::select(keep)
      
      if(type == "growth"){
        m <- lmer(BA_change ~ BA_Nfixer_prop_total + 
                    BA_total1 + 
                    (1 | STATECD), 
                  data = dat,
                  na.action = na.omit)
        new_d <- data.frame(BA_Nfixer_prop_total=c(min(data$BA_Nfixer_prop_total,na.rm=T), max(data$BA_Nfixer_prop_total,na.rm=T)),
                            BA_total1=c(mean(data$BA_total1,na.rm=T), mean(data$BA_total1,na.rm=T)))
        new_d$Pred <- predict(m,new_d,re.form=~0)
        pred.y.0[i,1] <- new_d$Pred[1]
        pred.y.0[i,2] <- new_d$Pred[2]
        r2[i,] <- r.squaredGLMM(m)
      }else if(type == "recr" | type == "recruitment"){
        m <- lmer(recyr ~ BA_Nfixer_prop_total +
                    BA_total1 +
                    (1 | STATECD),
                  data = dat)
        new_d <- data.frame(BA_Nfixer_prop_total=c(min(data$BA_Nfixer_prop_total,na.rm=T), max(data$BA_Nfixer_prop_total,na.rm=T)),
                            BA_total1=c(mean(data$BA_total1,na.rm=T), mean(data$BA_total1,na.rm=T)))
        new_d$Pred <- predict(m,new_d,re.form=~0)
        pred.y.0[i,1] <- new_d$Pred[1]
        pred.y.0[i,2] <- new_d$Pred[2]
        r2[i,] <- r.squaredGLMM(m)
      }else if(type == "mort" | type == "mortality" | type == "surv"){
        m <- lmer(survyr ~ BA_Nfixer_prop_total +
                    BA_total1 +
                    (1 | STATECD),
                  data = dat)
        new_d <- data.frame(BA_Nfixer_prop_total=c(min(data$BA_Nfixer_prop_total,na.rm=T), max(data$BA_Nfixer_prop_total,na.rm=T)),
                            BA_total1=c(mean(data$BA_total1,na.rm=T), mean(data$BA_total1,na.rm=T)))
        new_d$Pred <- predict(m,new_d,re.form=~0)
        pred.y.0[i,1] <- new_d$Pred[1]
        pred.y.0[i,2] <- new_d$Pred[2]
        r2[i,] <- r.squaredGLMM(m)
      }else if(type == "growthnf"){
        m <- lmer(BA_change_nonfixer ~ BA_Nfixer_prop_total + 
                    BA_total1 + 
                    (1 | STATECD), 
                  data = dat,
                  na.action = na.omit)
        new_d <- data.frame(BA_Nfixer_prop_total=c(min(data$BA_Nfixer_prop_total,na.rm=T), max(data$BA_Nfixer_prop_total,na.rm=T)),
                            BA_total1=c(mean(data$BA_total1,na.rm=T), mean(data$BA_total1,na.rm=T)))
        new_d$Pred <- predict(m,new_d,re.form=~0)
        pred.y.0[i,1] <- new_d$Pred[1]
        pred.y.0[i,2] <- new_d$Pred[2]
        r2[i,] <- r.squaredGLMM(m)
      }else if(type == "growthf"){
        m <- lmer(BA_change_Nfixer ~ BA_Nfixer_prop_total + 
                    BA_total1 + 
                    (1 | STATECD), 
                  data = dat,
                  na.action = na.omit)
        new_d <- data.frame(BA_Nfixer_prop_total=c(min(data$BA_Nfixer_prop_total,na.rm=T), max(data$BA_Nfixer_prop_total,na.rm=T)),
                            BA_total1=c(mean(data$BA_total1,na.rm=T), mean(data$BA_total1,na.rm=T)))
        new_d$Pred <- predict(m,new_d,re.form=~0)
        pred.y.0[i,1] <- new_d$Pred[1]
        pred.y.0[i,2] <- new_d$Pred[2]
        r2[i,] <- r.squaredGLMM(m)
      }
      
    }else{
      keep <- c("STATECD","BA_change","recyr","survyr","BA_change_nonfixer","BA_Nfixer_prop_total","BA_total1","BA_total2","BA_Nfixer1",group)
      dat <- dat %>% dplyr::select(keep)
      colnames(dat) <- c("STATECD","BA_change","recyr","survyr","BA_change_nonfixer","BA_Nfixer_prop_total","BA_total1","BA_total2","BA_Nfixer1", "group")
      
      if(type == "growth"){
        m <- lmer(BA_change ~ BA_Nfixer_prop_total + 
                    BA_total1 + 
                    BA_Nfixer_prop_total*group +
                    (1 | STATECD), 
                  data = dat,
                  na.action = na.omit)
        new_d <- data.frame(BA_Nfixer_prop_total=rep(c(min(data$BA_Nfixer_prop_total,na.rm=T), max(data$BA_Nfixer_prop_total,na.rm=T)),2),
                            BA_total1=rep(mean(data$BA_total1,na.rm=T), 4),
                            group=as.factor(c(0,0,1,1)))
        new_d$Pred <- predict(m,new_d,re.form=~0)
        pred.y.0[i,1] <- new_d$Pred[1]
        pred.y.0[i,2] <- new_d$Pred[2]
        pred.y.1[i,1] <- new_d$Pred[3]
        pred.y.1[i,2] <- new_d$Pred[4]
        r2[i,] <- r.squaredGLMM(m)
      }else if(type == "recr" | type == "recruitment"){
        m <- lmer(recyr ~ BA_Nfixer_prop_total +
                    BA_total1 +
                    BA_Nfixer_prop_total*group +
                    (1 | STATECD),
                  data = dat)
        new_d <- data.frame(BA_Nfixer_prop_total=rep(c(min(data$BA_Nfixer_prop_total,na.rm=T), max(data$BA_Nfixer_prop_total,na.rm=T)),2),
                            BA_total1=rep(mean(data$BA_total1,na.rm=T), 4),
                            group=as.factor(c(0,0,1,1)))
        new_d$Pred <- predict(m,new_d,re.form=~0)
        pred.y.0[i,1] <- new_d$Pred[1]
        pred.y.0[i,2] <- new_d$Pred[2]
        pred.y.1[i,1] <- new_d$Pred[3]
        pred.y.1[i,2] <- new_d$Pred[4]
        r2[i,] <- r.squaredGLMM(m)
      }else if(type == "mort" | type == "mortality" | type == "surv"){
        m <- lmer(survyr ~ BA_Nfixer_prop_total +
                    BA_total1 +
                    BA_Nfixer_prop_total*group +
                    (1 | STATECD),
                  data = dat)
        new_d <- data.frame(BA_Nfixer_prop_total=rep(c(min(data$BA_Nfixer_prop_total,na.rm=T), max(data$BA_Nfixer_prop_total,na.rm=T)),2),
                            BA_total1=rep(mean(data$BA_total1,na.rm=T), 4),
                            group=as.factor(c(0,0,1,1)))
        new_d$Pred <- predict(m,new_d,re.form=~0)
        pred.y.0[i,1] <- new_d$Pred[1]
        pred.y.0[i,2] <- new_d$Pred[2]
        pred.y.1[i,1] <- new_d$Pred[3]
        pred.y.1[i,2] <- new_d$Pred[4]
        r2[i,] <- r.squaredGLMM(m)
      }else if(type == "growthnf"){
        m <- lmer(BA_change_nonfixer ~ BA_Nfixer_prop_total + 
                    BA_total1 + 
                    BA_Nfixer_prop_total*group +
                    (1 | STATECD), 
                  data = dat,
                  na.action = na.omit)
        new_d <- data.frame(BA_Nfixer_prop_total=rep(c(min(data$BA_Nfixer_prop_total,na.rm=T), max(data$BA_Nfixer_prop_total,na.rm=T)),2),
                            BA_total1=rep(mean(data$BA_total1,na.rm=T), 4),
                            group=as.factor(c(0,0,1,1)))
        new_d$Pred <- predict(m,new_d,re.form=~0)
        pred.y.0[i,1] <- new_d$Pred[1]
        pred.y.0[i,2] <- new_d$Pred[2]
        pred.y.1[i,1] <- new_d$Pred[3]
        pred.y.1[i,2] <- new_d$Pred[4]
        r2[i,] <- r.squaredGLMM(m)
      }
    }
    
  }#end for loop
  
  if(group == "NA"){
    result <- cbind.data.frame(pred.y.0, r2)
    colnames(result) <- c("group0low", "group0high", "R2m", "R2c")
    
  }else{
    result <- cbind.data.frame(pred.y.0, pred.y.1, r2)
    colnames(result) <- c("group0low", "group0high", "group1low", "group1high", "R2m", "R2c")
    
  }
  
  return(result)
  
}


#test <- f_nonpar_boot(data=FIA_plot, type="growth", nboot=3, group="NDEPg")
#NDEPg, SMg, TEXg, CUNg, CSOg, YOU0OLD1



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

#testch <- pct_change_func_byrow(test)



f_plot <- function(modelout, type, title){
  plot <- ggplot(data = modelout, 
                 aes(x = factor(group), y = mean, width = 0.65)) +
    geom_bar(stat = "identity", 
             position = position_dodge()) + 
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Group", y = "% Change", subtitle = type, title = title,
         caption="Number in 1000s of samples") +
    geom_hline(yintercept = 0, color = "black") +
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                  position=position_dodge(.9)) +
    # geom_text(aes(label=round(Freq/1000, 1)), position=position_dodge(0.9), vjust=-0.2) +
    theme(panel.background = element_rect(fill = 'white', colour = 'white'),
          plot.title = element_text(hjust = 0.5, size = 20),
          text = element_text(size = 20),
          axis.text.x = element_text(angle = 90, hjust = 1, size = 20),
          plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0))) 
  return(plot)
}
  