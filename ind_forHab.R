

library(dplyr, lib.loc = "/rigel/home/arp2195/rpackages")
library(lme4, lib.loc = "/rigel/home/arp2195/rpackages")
library(pscl, lib.loc = "/rigel/home/arp2195/rpackages")

#don't use sci notation
options(scipen=6)

m.r <- zeroinfl(reci ~ - 1 +
                  NCI_props + 
                  NCIs +
                  dbhs +
                  1 | pcn1 + 1 | GENUS.x, 
                data = rs,
                dist = "negbin",
                offset = log(t)) #https://stats.idre.ucla.edu/r/dae/zinb/
saveRDS(m.r, "mod_recr_zi_nbin.RDS")

m.s <- zeroinfl(survi ~ - 1 +
                  NCI_props + 
                  NCIs +
                  dbhs +
                  1 | pcn1, 
                data = ms,
                dist = "negbin",
                offset = log(t)) #https://stats.idre.ucla.edu/r/dae/zinb/
saveRDS(m.s, "mod_surv_zi_nbin.RDS")