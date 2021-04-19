# standrc tests

library(standrc)
library(dplyr)

source('scripts/standrc_evals.R')
options(scipen = 999)
# data(spinach)
harmony_collection <- readRDS("temp_files/do_focus_merge_2_forstan.RData")


hoescht_features <- data.frame("feature" = c(
  "multi_nucleated_cells_numberofobjects", # similar
  "infocus_intensity_cell_hoechst_mean_mean", # similar
  "infocus_nucleus_hoechst_ser_dark_1px_mean", 
  "infocus_nucleus_hoechst_ser_edge_1px_mean" # similar
))

mitotracker_features <- data.frame("feature" = c(
  "nuclei_numberofobjects", # similar
  "infocus_intensity_cell_mitotrackerdeepred_mean_mean", # similar
  "infocus_cell_area_mm2_mean",
  "infocus_cell_mitotrackerdeepred_ser_edge_1px_mean" # similar
))

features_selected <- c(hoescht_features$feature, 
                       mitotracker_features$feature)

feature_subset <- harmony_collection$analytics_subset %>%
  group_by(individual, dose, dir) %>%
  summarise_at(vars(c(
    features_selected
  )), funs(m = mean(., na.rm = TRUE)))

feature_subset <- feature_subset %>%
  group_by(dir) %>%
  mutate(dir_id = cur_group_id()) %>%
  group_by(individual) %>%
  mutate(individual_id = cur_group_id())


feature_subset <- feature_subset[!is.na(feature_subset$infocus_intensity_cell_mitotrackerdeepred_mean_mean_m), ] 
feature_subset$response <- feature_subset$infocus_intensity_cell_mitotrackerdeepred_mean_mean_m
#feature_subset$response <- scales::rescale(feature_subset$response, to = c(0.01,100), from = range(feature_subset$response, na.rm = TRUE, finite = TRUE))

rm(standrm_robert)
spm <- standrm_robert(formula=response ~ dose, data=feature_subset, 
                      fct=LL.4(),
                      curveid= b ~ individual_id,
                      random=c + d + e ~ dir_id,
                      chains = 3, iter = 1000, warmup = 500, thin = 10)

print(spm)
plot(spm, ndose=25, logx=TRUE, legend = FALSE) + theme(legend.position = "none")
ED(spm, respLev=50)
