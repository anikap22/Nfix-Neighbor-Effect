
library(cowplot)
library(ggthemes)

setwd("/Users/Anika/Documents/GradSchool/FIA_CompetitionModel/")

# Figure 1, Individual bar graph -------------------------------------
#a) growth

#b) mortality
changes <- readRDS("output/changes_growth.RDS")

png("output/fig1b.png", 
    width = 8, 
    height = 4, 
    units = 'in', 
    res = 300)
ggplot(data = changes, 
       aes(x = factor(group), y = change, fill = status, width = 0.65)) +
  geom_bar(stat = "identity", 
           position = position_dodge()) + 
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Group", y = "% Change", subtitle = "Survival") +
  geom_hline(yintercept = 0, color = "black") +
  geom_errorbar(aes(ymin = change-se, ymax = change+se), width = 0.2,
                position = position_dodge(width = 0.65)) 
dev.off()

#c) recr

# Figure 2, % change maps --------------------------------------------

# Figure 3, % change by forest type ----------------------------------

# Figure 4, evergreen vs deciduous, traits ---------------------------
#a) growth
changes_form <- readRDS("output/changes_form_growth.RDS")

png("output/fig4a.png", 
    width = 8, 
    height = 4, 
    units = 'in', 
    res = 300)
ggplot(data = changes_form, 
       aes(x = factor(group), y = change, width = 0.65)) +
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Group", y = "% Change") +
  geom_hline(yintercept = 0, color = "black") + 
  theme(legend.position = "none") +
  geom_errorbar(aes(ymin = change-se, ymax = change+se), width = 0.2,
                position = position_dodge(width = 0.65))

dev.off()

#b) growth, evergreen
changes_ever <- readRDS("output/changes_non_ever_growth.RDS")
changes_ever <- changes_ever[c(7,8,17,18,23,24,25,26,29,30), ]

png("output/fig4b.png", 
    width = 8, 
    height = 4, 
    units = 'in', 
    res = 300)
ggplot(data = changes_ever, 
       aes(x = factor(group), y = change, fill = status, width = 0.65)) +
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Group", y = "% Change") +
  geom_hline(yintercept = 0, color = "black") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_errorbar(aes(ymin = change-se, ymax = change+se), width = 0.2,
                position = position_dodge(width = 0.65)) +
  scale_fill_grey(start = 0.2, end = 0.8)

dev.off()


#c) growth, decid
changes_decid <- readRDS("output/changes_non_decid_growth.RDS")
changes_decid <- changes_decid[c(7,8,17,18,23,24,25,26,29,30), ]

png("output/fig4c.png", 
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
                position = position_dodge(width = 0.65)) +
  scale_fill_grey(start = 0.2, end = 0.8)
dev.off()


#d) mort
changes_form <- readRDS("output/changes_form_mort.RDS")

png("output/fig4d.png", 
    width = 8, 
    height = 4, 
    units = 'in', 
    res = 300)
ggplot(data = changes_form, 
       aes(x = factor(group), y = change, width = 0.65)) +
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Group", y = "% Change") +
  geom_hline(yintercept = 0, color = "black") + 
  theme(legend.position = "none") +
  geom_errorbar(aes(ymin = change-se, ymax = change+se), width = 0.2,
                position = position_dodge(width = 0.65))

dev.off()


#e) mort, evergreen
changes_ever <- readRDS("output/changes_non_ever_mort.RDS")
changes_ever <- changes_ever[c(7,8,17,18,23,24,25,26,29,30), ]

png("output/fig4e.png", 
    width = 8, 
    height = 4, 
    units = 'in', 
    res = 300)
ggplot(data = changes_ever, 
       aes(x = factor(group), y = change, fill = status, width = 0.65)) +
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Group", y = "% Change") +
  geom_hline(yintercept = 0, color = "black") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_errorbar(aes(ymin = change-se, ymax = change+se), width = 0.2,
                position = position_dodge(width = 0.65)) +
  scale_fill_grey(start = 0.2, end = 0.8)
dev.off()

#f) mort, decid
changes_decid <- readRDS("output/changes_non_decid_mort.RDS")
changes_decid <- changes_decid[c(7,8,17,18,23,24,25,26,29,30), ]

png("output/fig4f.png", 
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
                position = position_dodge(width = 0.65)) +
  scale_fill_grey(start = 0.2, end = 0.8)
dev.off()

#g) recr
changes_form <- readRDS("output/changes_form_recr.RDS")

png("output/fig4g.png", 
    width = 8, 
    height = 4, 
    units = 'in', 
    res = 300)
ggplot(data = changes_form, 
       aes(x = factor(group), y = change, width = 0.65)) +
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Group", y = "% Change") +
  geom_hline(yintercept = 0, color = "black") + 
  theme(legend.position = "none") +
  geom_errorbar(aes(ymin = change-se, ymax = change+se), width = 0.2,
                position = position_dodge(width = 0.65))

dev.off()



#h) recr, evergreen
changes_ever <- readRDS("output/changes_non_ever_recr.RDS")
changes_ever <- changes_ever[c(7,8,17,18,23,24,25,26,29,30), ]

png("output/fig4h.png", 
    width = 8, 
    height = 4, 
    units = 'in', 
    res = 300)
ggplot(data = changes_ever, 
       aes(x = factor(group), y = change, fill = status, width = 0.65)) +
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Group", y = "% Change") +
  geom_hline(yintercept = 0, color = "black") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_errorbar(aes(ymin = change-se, ymax = change+se), width = 0.2,
                position = position_dodge(width = 0.65)) +
  scale_fill_grey(start = 0.2, end = 0.8)
dev.off()


#i) recr, decid
changes_decid <- readRDS("output/changes_non_decid_recr.RDS")
changes_decid <- changes_decid[c(7,8,17,18,23,24,25,26,29,30), ]

png("output/fig4i.png", 
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
                position = position_dodge(width = 0.65)) +
  scale_fill_grey(start = 0.2, end = 0.8)
dev.off()


# Figure 5, Plot level bar graph -------------------------------------

#a) biomass accumulation
changes <- readRDS("output/plot_changes_biomassacc.RDS")

png("output/fig5a.png", 
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

dev.off()

#b) survival
changes <- readRDS("output/plot_changes_surv.RDS")

png("output/fig5b.png", 
    width = 8,
    height = 4,
    units = 'in',
    res = 300)
ggplot(data = changes, 
       aes(x = factor(group), y = change, fill = status, width = 0.65)) +
  geom_bar(stat = "identity", 
           position = position_dodge()) + 
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Group", y = "% Change", subtitle = "Survival - Plot level") +
  geom_hline(yintercept = 0, color = "black") +
  geom_errorbar(aes(ymin = change-se, ymax = change+se), width = 0.2,
                position = position_dodge(width = 0.65)) 

dev.off()

#c) recruitment
changes <- readRDS("output/plot_changes_recr.RDS")

png("output/fig5c.png", 
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

dev.off()



# Figure 6, N deposition ---------------------------------------------

ndep_gr$grp <- "growth"
ndep_mort$grp <- "survival"
ndep_recr$grp <- "recruitment"

ndep <- rbind(ndep_gr, ndep_mort, ndep_recr)
saveRDS(ndep, "output/fig6_data.RDS")

ndep <- readRDS("output/fig6_data.RDS")

png("output/fig6.png", 
    width = 8, 
    height = 4, 
    units = 'in', 
    res = 300)

ggplot(data = ndep, 
       aes(x = factor(quad), y = est, fill = grp, width = 0.65)) +
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1") +
  labs(x = "N deposition (kgN/ha/yr)", 
       y = "% Change", 
       title = "% Change in demographic rate by N deposition level",
       fill = "Demographic Rate") +
  geom_hline(yintercept = 0, color = "black") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_errorbar(aes(ymin = est-se, ymax = est+se), width = 0.2,
                position = position_dodge(width = 0.65)) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

dev.off()
