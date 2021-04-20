# standrc tests
# https://github.com/daniel-gerhard/standrc/tree/master/R
library(standrc)
library(dplyr)

options(scipen = 999)
setwd("/Users/seller/Desktop/projects/jax_arsenic_toxicology_analysis")
harmony_collection <- readRDS("temp_files/do_focus_merge_2_forstan.RData")

try({rm(standrm_robert);rm(plot_new)})
source('scripts/standrc_evals.R')

feature_subset <- harmony_collection$analytics_subset %>%
  group_by(individual, dose, dir) %>%
  summarise_at(vars(features_selected <- c(
    "multi_nucleated_cells_numberofobjects", # similar
    "infocus_intensity_cell_hoechst_mean_mean", # similar
    "infocus_nucleus_hoechst_ser_dark_1px_mean", 
    "infocus_nucleus_hoechst_ser_edge_1px_mean",
    "nuclei_numberofobjects", # similar
    "infocus_intensity_cell_mitotrackerdeepred_mean_mean", # similar
    "infocus_cell_area_mm2_mean",
    "infocus_cell_mitotrackerdeepred_ser_edge_1px_mean"
  )), funs(m = mean(., na.rm = TRUE))) %>%
    group_by(dir) %>%
    mutate(dir_id = cur_group_id()) %>%
    group_by(individual) %>%
    mutate(individual_id = cur_group_id())

prepare_data_for_stan <- function(ft, dat, log_ = FALSE, scale_ = TRUE){
  dat$response <- as.numeric(unlist(dat[,c(ft)]))
  dat <- dat[!is.na(dat$response), ] 
  if (log_) dat$response <- log(dat$response)
  if (scale_){
    dat$response <- scales::rescale(dat$response, 
                                    to = c(0.01,100), 
                                    from = range(dat$response), 
                                    na.rm = TRUE, 
                                    finite = TRUE)
  }
  return(dat)
}

###### 1 ######
feature_subset_intensity <- prepare_data_for_stan('infocus_intensity_cell_hoechst_mean_mean_m',
                                        feature_subset, 
                                        log_ = TRUE)

spm_intensity <- standrm_robert(formula=response ~ dose, data=feature_subset_intensity, 
                      fct = LL.4(fixed=c(NA,NA, NA,NA)), # slope, lower limit, upper limit, ed50
                      curveid= b + e ~ individual_id,
                      weights = nuclei_numberofobjects,
                      random=b + c + d ~ 1 | dir_id)

plot_new(spm_intensity, ndose=25, logx=TRUE, legend = FALSE) + theme(legend.position = "none")
ED(spm_intensity,respLev = 50)

###### 2 ######
# infocus_nucleus_hoechst_ser_dark_1px_mean_m
# rescale is enough

feature_subset_dark <- prepare_data_for_stan('infocus_nucleus_hoechst_ser_dark_1px_mean_m',
                                                  feature_subset, 
                                                  log_ = FALSE)

spm_dark <- standrm_robert(formula=response ~ dose, data=feature_subset_dark, 
                                fct = LL.4(fixed=c(NA,0, 100,NA)), # slope, lower limit, upper limit, ed50
                                # curveid= b ~ individual,
                                curveid= b + e ~ individual_id,
                                # random=c + d + e ~ dir/dose)
                                random=b + c + d  ~ 1 | dir_id)

plot_new(spm_dark, ndose=25, logx=TRUE, legend = FALSE) + theme(legend.position = "none")

###### 3 ######
# infocus_nucleus_hoechst_ser_edge_1px_mean_m
# rescale is enough

feature_subset_edge <- prepare_data_for_stan('infocus_nucleus_hoechst_ser_edge_1px_mean_m',
                                                  feature_subset, 
                                                  log_ = FALSE)

spm_edge <- standrm_robert(formula=response ~ dose, data = feature_subset_edge, 
                                fct = LL.4(fixed=c(NA, NA, NA, NA)), # slope, lower limit, upper limit, ed50
                                curveid= b + e ~ individual_id,
                                random=b + c + d ~ dir_id)

plot_new(spm_edge, ndose=25, logx=TRUE, legend = FALSE) + theme(legend.position = "none")
