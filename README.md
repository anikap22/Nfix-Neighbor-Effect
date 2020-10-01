# Nfix-Neighbor-Effect

This project used hierarchical models to examine the effect of nitrogen-fixing trees on forest demographics including: growth rate, recruitment rate, and survival rate at two scales: plot-scale and individual-scale.

Data was pre-processed, cleaned, and aggregated. Then a crowding index was calculated from the spatially explicit tree data which was carried out on the Columbia Habanero cluster for distributed computing. Modeling and analysis was carreid out and visualization followed.

## Pre-processing, cleaning, and aggregating:</br>
`join_new.R`</br>
`growth_prep.R`</br>
`mort_prep.R` </br>
`recr_prep.R` </br>
`plots_prep.R` </br>

## Crowding calculations:</br>
`ind_forHab.R`</br>
`nci_forhabanero.R` </br>
`nci_hab_processing.R` </br>
`nci.R` </br>
`nci_all_2ts_v4.R` </br>

## Further processing:</br>
`post_nci_proc.R` </br>
`smap.R` </br>
`usda_prep.R` </br>
`ind_define_groups.R`</br>


## Model development and analysis:</br>
### Individual scale: </br>
`models_ind_mort.R`</br>
`models_ind_recr.R`</br>
`models_ind_surv.R`</br>
`ind_grp_withuncert_10_26_19.R`</br>
### Plot scale: </br>
`plot_define_groups.R` </br>
`models_plot_withuncert_07_08_20.R` </br>

## Visualization:</br>
`figs.R`
