

# dummy data --------------------------------------

test <- data.frame('PCN' = c(1, 2, 3, 4, 5, 6, 7, 8, 9), 'PREV_PCN' = c(2, 3, NA, 5, NA, 7, 8, 9, NA))
test$PCN <- as.integer64(test$PCN)
test$PREV_PCN <- as.integer64(test$PREV_PCN)

round1 <- test[is.na(test$PREV_PCN),]
colnames(round1) <- c('t3', 't4')

r1t3 <- round1$t3[!is.na(round1$t3)]
round2 <- filter(test, test$PREV_PCN %in% r1t3)
colnames(round2) <- c('t2', 't3')
round2 <- merge(round1, round2, by.x='t3', by.y='t3', all.x=T)

r2t2 <- round2$t2[!is.na(round2$t2)]
round3 <- filter(test, test$PREV_PCN %in% r2t2)
colnames(round3) <- c('t1', 't2')
round3 <- merge(round2, round3, by='t2', all.x=T)

r3t1 <- round3$t1[!is.na(round3$t1)]
round4 <- filter(test, test$PREV_PCN %in% r3t1)
colnames(round4) <- c('t0', 't1')
round4 <- merge(round3, round4, by='t1', all.x=T)

round4 <- round4[,c('t0','t1','t2','t3','t4')]

round4$t0 <- as.character(round4$t0)
round4$t1 <- as.character(round4$t1)
round4$t2 <- as.character(round4$t2)
round4$t3 <- as.character(round4$t3)
round4$t4 <- as.character(round4$t4)

round4a <-  (t(apply(round4, 1, function(x) (c(x[which(!is.na(x))], x[which(is.na(x))])))))
round4b <- as.integer64(round4a)
len <- length(round4b)/5
dim(round4b) <- c(len, 5)
round4b <- data.frame(round4b)
colnames(round4b) <- c('t0','t1','t2','t3','t4')

p_ts <- round4b[,c('t0', 't1')]
colnames(p_ts) <- c('t2', 't1')
p_ts

# make p_ts with real data ------------------------------------

setwd("/Users/Anika/Documents/GradSchool/FIA_CompetitionModel/")

plot <- fread("data/raw/plot.csv", select = c('CN', 
                                              'PREV_PLT_CN', 
                                              'PLOT', 
                                              'MEASYEAR', 
                                              'WATERCD', 
                                              'LAT', 
                                              'LON', 
                                              'ELEV',
                                              'DESIGNCD',
                                              'STATECD'))

# correct formatting
plot$CN <- as.integer64(plot$CN)
plot$PREV_PLT_CN <- as.integer64(plot$PREV_PLT_CN)
plot <- data.table(plot)
plot[plot$PREV_PLT_CN == 0, "PREV_PLT_CN"] <- NA
colnames(plot) <- c("PCN", "PREV_PCN", "PLOT", "MEASYEAR", 'WATERCD', 'LAT', 'LON', 'ELEV',
                    'DESIGNCD', 'STATECD')


# remove data with no year data
plot <- plot[!is.na(plot$MEASYEAR), ]

# select only PCN and PREV_PCN
plotsm <- plot[, c('PCN', 'PREV_PCN')]

round1 <- plotsm[is.na(plotsm$PREV_PCN),]
colnames(round1) <- c('t3', 't4')

r1t3 <- round1$t3[!is.na(round1$t3)]
round2 <- filter(plotsm, plotsm$PREV_PCN %in% r1t3)
colnames(round2) <- c('t2', 't3')
#round2a <- merge(round1, round2, by='t3', all.x=T)

setkey(round1, t3)
round2 <- as.data.table(round2)
setkey(round2, t3)
round2 <- round2[round1]

r2t2 <- round2$t2[!is.na(round2$t2)]
round3 <- filter(plotsm, plotsm$PREV_PCN %in% r2t2)
colnames(round3) <- c('t1', 't2')
#round3 <- merge(round2, round3, by='t2', all.x=T)

setkey(round2, t2)
round3 <- as.data.table(round3)
setkey(round3, t2)
round3 <- round3[round2]

r3t1 <- round3$t1[!is.na(round3$t1)]
round4 <- filter(plotsm, plotsm$PREV_PCN %in% r3t1)
colnames(round4) <- c('t0', 't1')
#round4 <- merge(round3, round4, by='t1', all.x=T)

setkey(round3, t1)
round4 <- as.data.table(round4)
setkey(round4, t1)
round4 <- round4[round3]

round4 <- round4[,c('t0','t1','t2','t3','t4')]

round4$t0 <- as.character(round4$t0)
round4$t1 <- as.character(round4$t1)
round4$t2 <- as.character(round4$t2)
round4$t3 <- as.character(round4$t3)
round4$t4 <- as.character(round4$t4)

round4a <-  (t(apply(round4, 1, function(x) (c(x[which(!is.na(x))], x[which(is.na(x))])))))
round4b <- as.integer64(round4a)
len <- length(round4b)/5
dim(round4b) <- c(len, 5)
round4b <- data.frame(round4b)
colnames(round4b) <- c('t0','t1','t2','t3','t4')

p_ts1 <- round4b[,c('t0', 't1')]
colnames(p_ts1) <- c('t2', 't1')
p_ts <- p_ts1[!is.na(p_ts1$t1) & !is.na(p_ts1$t2),]



# remove plot designs that don't have 24' fixed radius subplot --------------------------
p1 <- select(plot, PCN, DESIGNCD)

p2 <- merge(p_ts, p1, by.x='t1', by.y="PCN", all.x=T, all.y=F)
colnames(p2) <- c('t1', 't2', 'DESIGNCD1')
p2 <- merge(p2, p1, by.x='t2', by.y='PCN', all.x=T, all.y=F)
colnames(p2) <- c('t1', 't2', 'DESIGNCD1', 'DESIGNCD2')

keepcodes <- c(1, 111:113, 230:242, 328, 501:505, 601:603, 311:323, 301)
p3 <- p2[p2$DESIGNCD1 %in% keepcodes, ]
p_ts <- p3

rm(p1, p2, p3)
# remove island states -----------------------
p1 <- select(plot, PCN, STATECD)

p2 <- merge(p_ts, p1, by.x='t1', by.y='PCN', all.x=T, all.y=F)

removecodes <- c(15, 60, 66, 70, 72, 78)
p3 <- p2[!p2$STATECD %in% removecodes, ]
p_ts <- p3
p_ts <- select(p_ts, t1, t2, STATECD)
colnames(p_ts) <- c('PCN1', 'PCN2', 'STATECD')

rm(p1, p2, p3, plotsm, round1, round2, round3, round4, round4a, r1t3, r2t2, r3t1, p_ts1)

# add trees @t2 (grew or recruited) ------------------------
# add trees to t2 using PCN from TREE
tree <- fread("data/raw/tree.csv", select=c("CN",
                                            "PLT_CN",
                                            "PREV_TRE_CN",
                                            "INVYR",
                                            "PLOT",
                                            "SUBP",
                                            "TREE",
                                            "AZIMUTH",
                                            "DIST",
                                            "SPCD",
                                            "DIA",
                                            "TPA_UNADJ",
                                            "CCLCD",
                                            "CARBON_AG",
                                            "CARBON_BG"))
tree <- tree[!is.na(tree$DIST) & !is.na(tree$AZIMUTH), ]

# correct formatting
tree$CN <- as.integer64(tree$CN)
tree$PLT_CN <- as.integer64(tree$PLT_CN)
tree$PREV_TRE_CN <- as.integer64(tree$PREV_TRE_CN)
tree[tree$PREV_TRE_CN == 0, "PREV_TRE_CN"] <- NA
tree <- data.table(tree)
colnames(tree) <- c("TCN", "PCN", "PREV_TCN", "INVYR", "PLOT", "SUBP", "TREE",
                    "AZIMUTH", "DIST", "SPCD", "DIA", "TPA_UNADJ", "CCLCD", "CARBON_AG",
                    "CARBON_BG")

t_ts <- merge(p_ts, tree, by.x = 'PCN2', by.y = 'PCN', all.x=T, all.y=F)

colnames(t_ts) <- c('PCN2', 'PCN1', 'STATECD', 'TCN2', 'PREV_TCN2', 'INVYR2', 'PLOT2',
                    'SUBP2', 'TREE2', 'AZIMUTH2', 'DIST2', 'SPCD2', 'DIA2', 'TPA_UNADJ2',
                    'CCLCD2', 'CARBON_AG2', 'CARBON_BG2')


# add trees @t1 (grew) ------------------------------------
# use prev_tcn for trees at t1 matching those added to t2

t_ts1 <- merge(t_ts, tree, by.x = 'PREV_TCN2', by.y = "TCN", all.x=T, all.y=F)
colnames(t_ts1) <- c('TCN1', 'PCN2', 'PCN1', 'STATECD', 'TCN2', 'INVYR2', 'PLOT2',
                     'SUBP2', 'TREE2', 'AZIMUTH2', 'DIST2', 'SPCD2', 'DIA2', 'TPA_UNADJ2',
                     'CCLCD2', 'CARBON_AG2', 'CARBON_BG2', 'PCN0', 'PREV_TCN0', 'INVYR1',
                     'PLOT1', 'SUBP1', 'TREE1', 'AZIMUTH1', 'DIST1', 'SPCD1', 'DIA1',
                     'TPA_UNADJ1', 'CCLCD1', 'CARBON_AG1', 'CARBON_BG1')

t_ts1 <- data.table(t_ts1)

# add trees @t1 (died) --------------------------------------
# add to t1 that match PCN but no PREV_TCN

diedlist <- tree[(tree$PCN %in% t_ts1$PCN1) & !(tree$TCN %in% t_ts1$TCN1) & !(tree$TCN %in% t_ts1$TCN2), ]
died_ts <- merge(p_ts, diedlist, by.x="PCN1", by.y="PCN", all.x=T, all.y=F)
colnames(died_ts) <- c('PCN1', 'PCN2', 'STATECD', 'TCN1', 'PREV_TCN1', 'INVYR1',
                       'PLOT1', 'SUBP1', 'TREE1', 'AZIMUTH1', 'DIST1', 'SPCD1',
                       'DIA1', 'TPA_UNADJ1', 'CCLCD1', 'CARBON_AG1', 'CARBON_BG1')
died_ts <- died_ts[!is.na(died_ts$INVYR1), ] #remove plots with no trees

# make NA columns for time 2
namevector <- c('TCN2', 'INVYR2', 'PLOT2', 'SUBP2', 'TREE2', 'AZIMUTH2',
                'DIST2', 'SPCD2', 'DIA2', 'TPA_UNADJ2', 'CCLCD2', 'CARBON_AG2', 'CARBON_BG2')
died_ts[ ,namevector] <- NA
died_ts <- select(died_ts, -PREV_TCN1)
died_ts$TCN2 <- as.integer64(died_ts$TCN2)
died_ts$INVYR2 <- as.integer(died_ts$INVYR2)
died_ts$PLOT2 <- as.integer(died_ts$PLOT2)
died_ts$SUBP2 <- as.integer(died_ts$SUBP2)
died_ts$TREE2 <- as.integer(died_ts$TREE2)
died_ts$AZIMUTH2 <- as.integer(died_ts$AZIMUTH2)
died_ts$DIST2 <- as.numeric(died_ts$DIST2)
died_ts$SPCD2 <- as.integer(died_ts$SPCD2)
died_ts$DIA2 <- as.integer(died_ts$DIA2)
died_ts$TPA_UNADJ2 <- as.integer(died_ts$TPA_UNADJ2)
died_ts$CCLCD2 <- as.integer(died_ts$CCLCD2)
died_ts$CARBON_AG2 <- as.numeric(died_ts$CARBON_AG2)
died_ts$CARBON_BG2 <- as.numeric(died_ts$CARBON_BG2)

# algin colnames for t_ts1 with trees that grew or recruited
t_ts1 <- select(t_ts1, -c(PCN0, PREV_TCN0))

t_ts2 <- rbind.data.frame(t_ts1, died_ts)

#remove plots with no trees
t_ts <- t_ts2[!is.na(t_ts2$INVYR1) | !is.na(t_ts2$INVYR2), ]

#remove plots with trees having distance > 24 feet
temp <- t_ts[t_ts$DIST1 > 24, ]
temp1 <- unique(temp$PCN1)
t_ts <- t_ts[!(t_ts$PCN1 %in% temp1), ]

#add lat, lon, elev
p1 <- select(plot, PCN, LAT, LON, ELEV)
t_ts <- merge(t_ts, p1, by.x = 'PCN1', by.y='PCN', all.x=T, all.y=F)


# add fates
# used DIA for fate b/c some have tcn number but no tree actually measured
t_ts[is.na(DIA1) & !is.na(DIA2), FATE := 1] #recruited  DT[i, col := val] 
t_ts[!is.na(DIA1) & !is.na(DIA2), FATE := 2] #grew
t_ts[!is.na(DIA1) & is.na(DIA2), FATE := 3] #morality

hist(t_ts[,FATE], prob=T)

t_ts[FATE==1, AZIMUTH := AZIMUTH2]
t_ts[FATE==2, AZIMUTH := AZIMUTH1] #actually we want mean of DIA1 and DIA2 here
t_ts[FATE==3, AZIMUTH := AZIMUTH1]

t_ts[FATE==1, DIST := DIST2]
t_ts[FATE==2, DIST := DIST1] #actually we want mean of DIA1 and DIA2 here
t_ts[FATE==3, DIST := DIST1]

t_ts[FATE==1, SPCD := SPCD2]
t_ts[FATE==2, SPCD := SPCD1] #actually we want mean of DIA1 and DIA2 here
t_ts[FATE==3, SPCD := SPCD1]

t_ts[FATE==1, TPA_UNADJ := TPA_UNADJ2]
t_ts[FATE==2, TPA_UNADJ := TPA_UNADJ1] #actually we want mean of DIA1 and DIA2 here
t_ts[FATE==3, TPA_UNADJ := TPA_UNADJ1]


#save output
saveRDS(t_ts, "output/t_ts_new.RDS")
