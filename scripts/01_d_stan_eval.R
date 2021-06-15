# standrc tests
# https://github.com/daniel-gerhard/standrc/tree/master/R
library(standrc)
library(dplyr)

options(scipen = 999)
setwd("/Users/seller/Desktop/projects/jax_arsenic_toxicology_analysis")
harmony_collection <- readRDS("temp_files/do_focus_merge_2_forstan.RData")

try({rm(standrm_harmony);rm(plot_new)})
source('scripts/standrc_evals.R')
source('scripts/custom_tools.R')

feature_subset <- harmony_collection$analytics_subset %>%
  group_by(dir) %>%
  mutate(dir_id = cur_group_id()) %>%
  group_by(individual) %>%
  mutate(individual_id = cur_group_id())

prepare_data_for_stan <- function(ft, dat, log_ = FALSE, scale_ = TRUE){
  dat$response <- as.numeric(unlist(dat[,c(ft)]))
  dat <- dat[!is.na(dat$response), ] 
  if (log_) {
      dat$response <- log(dat$response)
  }
  if (scale_){
    dat$response <- scales::rescale(dat$response, 
                                    to = c(0.01,100), 
                                    from = range(dat$response), # add outlier reduction here?
                                    na.rm = TRUE, 
                                    finite = TRUE)
  }
  return(dat)
}


###### 1 ######
# https://github.com/daniel-gerhard/standrc/blob/master/vignettes/standrc_vignette.Rmd
feature_subset_intensity <- prepare_data_for_stan('infocus_intensity_cell_hoechst_mean_mean',
                                        feature_subset, 
                                        log_ = FALSE)


spm_intensity <- standrm_harmony(formula=response ~ dose, 
                                 data=feature_subset_intensity, 
                            fct = LL.4(), # slope, lower limit, upper limit, ed50
                            curveid= b + e ~ individual_id,
                            standrc_priors(e="normal(pe, 5)", 
                                           b="normal(pb, 1)"), 
                            random=b + c + d  ~ dir_id,
                            iter = 2000
                            # cores = 8,
                            # chains = 8
                            #iter = 3000
                            )

png(file="./output/stan_tests/spm_intensity_13.png",
    width=600, height=350)
plot_new(spm_intensity, ndose=25, logx=TRUE, legend = FALSE) + theme(legend.position = "none")
dev.off()

spm_intensity_ed50 <- ED(spm_intensity, respLev=c(50))
# print(spm_intensity) 
# fixef(spm_intensity) # almost matches fixef slopes == per mouse
# ranef(spm_intensity) # 80 random slopes, 80 random upper, 80 random lower
# VarCorr(spm_intensity) # ?
# fitted(spm_intensity) # output ?
#residual(spm_intensity) # doesn't run
# predict

###### 2 ######
# infocus_nucleus_hoechst_ser_dark_1px_mean_m
# rescale is enough
feature_subset_dark <- prepare_data_for_stan('infocus_nucleus_hoechst_ser_dark_1px_mean',
                                                  feature_subset, 
                                                  log_ = FALSE)

spm_dark <- standrm_harmony(formula=response ~ dose, data=feature_subset_dark, 
                            fct = LL.4(), # slope, lower limit, upper limit, ed50
                            curveid= b + e ~ individual_id,
                            standrc_priors(e="normal(pe, 5)", 
                                           b="normal(pb, 1)"), 
                            random=b + c + d  ~ dir_id,
                            iter = 2000
                            # cores = 8,
                            # chains = 8
                            #iter = 3000
)

png(file="./output/stan_tests/spm_dark_2.png",
    width=600, height=350)
plot_new(spm_dark, ndose=25, logx=TRUE, legend = FALSE) + theme(legend.position = "none")
dev.off()

spm_dark_ed50 <- ED(spm_dark, respLev=c(50))
###### 3 ######
# infocus_nucleus_hoechst_ser_edge_1px_mean_m
# rescale is enough

feature_subset_edge <- prepare_data_for_stan('infocus_nucleus_hoechst_ser_edge_1px_mean',
                                                  feature_subset, 
                                                  log_ = FALSE)

spm_edge <- standrm_harmony(formula=response ~ dose, data = feature_subset_edge, 
                            fct = LL.4(), # slope, lower limit, upper limit, ed50
                            curveid= b + e ~ individual_id,
                            standrc_priors(e="normal(pe, 5)", 
                                           b="normal(pb, 1)"), 
                            random=b + c + d  ~ dir_id,
                            iter = 2000
                            # cores = 8,
                            # chains = 8
                            #iter = 3000
)

png(file="./output/stan_tests/spm_edge_2.png",
    width=600, height=350)
plot_new(spm_edge, ndose=25, logx=TRUE, legend = FALSE) + theme(legend.position = "none")
dev.off()

spm_edge_ed50 <- ED(spm_edge, respLev=c(50))



######
######
######
newd <- expand.grid("dose" = spm_edge_ed50[[1]][1,])
pm <- predict_harmony(spm_edge, newdata=newd)
pred_mito
