lme4_extract <- function(formula_, name, data_, m_list){

  if (paste(as.character(format(formula(paste0(name,formula_))))) %notin% m_list$summary_stats$formula){
    #browser()
    m <- lmer(formula(paste0(name,formula_)), data = data_, REML = FALSE) 
    formula_string <- paste(as.character(format(formula(paste0(name,formula_)))), collapse = '')
    fscore <- as.data.frame(anova(m))[,4]
    # aic <- ranova(m)$AIC[length(ranova(m)$AIC)]
    # log_lik <- ranova(m)$logLik[length(ranova(m)$logLik)]
    m_list$summary_stats[nrow(m_list$summary_stats)+1,] <- c(
      formula = formula_string, # formula
      f_value = fscore #F-score
      # aic = aic, # AIC (lower better?)
      # loglik = log_lik (higher better?)
      ) 
    m_list$summary_stats<- m_list$summary_stats %>% distinct()
    m_list$models <- c(m_list$models, m)
    tryCatch({
      m_list$ranova <- ranova(m)
      m_list$anova <- anova(m)
    },
      error=function(cond) {
          stop(paste0("Failed ranova, delete " , formula_string))
      })
    
  }
  return (m_list)
}

try_plot_effects <- function(model, term, response){
  df <- harmony_collection$analytics_subset
  # look at this confint(model_results$models$nested_compare_b, method="Wald")
  effects_dose <- effects::effect(term= term, mod= model)
  effects_dose_df <- as.data.frame(effects_dose)
  
  effects_dose_df[,c(term)] <- as.numeric(effects_dose_df[,c(term)] )
  ggplot() + 
   geom_boxplot(data=df, aes_string(x = term, y = response), show.legend = FALSE, outlier.size = -1, na.rm = TRUE) +
 geom_point(data = df, aes_string('fdose', response), na.rm = TRUE) + 
geom_point(data=effects_dose_df, aes_string(term, y='fit'), color="blue") +
   geom_line(data=effects_dose_df, aes_string(term, y='fit'), color="blue") +
   geom_ribbon(data= effects_dose_df, aes_string(term,  ymin='lower', ymax='upper'), alpha= 0.3, fill="blue") +
   labs(x=term, y=response) + small_font_theme 
}

# combine everything for quick outputs
multi_process_v1 <- function(data, formula_list, predictor){
  # build object
  
  # formula_list <- sapply(formula_list, function(x) as.formula(paste(predictor,x)))
  
  model_results <- list(summary_stats = data.frame(formula = character(),
                          f_value = numeric(),
                          aic = numeric(),
                          loglik = numeric()),
                          models = list(),
                          plots_1 = list())
  
  # loop through each lme4 formula and save model
  for (formula_ in names(formula_list)){
    print(as.character(formula_))
    model_results <- lme4_extract(
      name = predictor,
      formula_ = formula_list[[formula_]],
      data_ = data,  
      m_list = model_results
      )
  }
  model_results$response <- predictor
  names(model_results$models) <- names(formula_list)
  return(model_results)
}

plot_lme4_1 <- function(model_results){
  plot_list <- list()
  # the effects::effects function is poorly written
  # need to derive these another way
  df <- harmony_collection$analytics_subset
  for(m in 1:length(model_results$models)){
    plot_list[[names(model_results$models[m])]] <- try_plot_effects(
     model = model_results$models[[m]],
     term  = "fdose",
     response = model_results$response
    )
  }
  ggarrange(plotlist=plot_list, widths = c(2,2,2), labels = names(model_results$models))
}

factor_index_rescale <- function(fdose, estimates){
      # rescale by factored index
    y <- as.numeric(levels(fdose))
    z <- findInterval(estimates, y)
    z + ((estimates-y[z])/(y[findInterval(estimates, y)+1]-y[z]))
}

multi_process_predict <- function(data, model_results, selected_model, group_pair_, EC = 50){
  # apply prediction results
  data$pred <- predict(model_results$models[[selected_model]], data)
    data <- subset(data, group_pair == group_pair_)
    # check against drc vals
    df_drc <- drc::drm(pred ~ dose, 
                data = data,
                curveid = individual,
                type = 'continuous',
                fct = drc::LL.4(names = c("slope", "min_value", "max_value", "ec_50")))
    drc_coef <- drc::ED(df,EC)
    drc_coef_df <- as.data.frame(drc::ED(df_drc,EC))
    drc_coef_df$individual <- gsub(":.*","\\1",gsub(".*e:","",row.names(drc_coef_df)))
    data_merged <- merge(data, drc_coef_df, by.x=c("individual"), by.y=c("individual"))
    data_merged$rescaled_ec50 <- factor_index_rescale(harmony_collection$analytics_subset$fdose,drc_coef_df$Estimate)#data_merged$res
  p <- ggplot(data_merged) + 
     geom_point(aes_string(x = 'fdose', y = model_results$response, shape = 'plate'), 
                color = 'orange', show.legend = TRUE, na.rm = TRUE) +
    geom_point(aes_string(x = 'fdose', y = 'pred', shape = 'group_pair'), 
               color="blue", size = 1, na.rm = TRUE, show.legend = FALSE) +
  geom_vline(aes(xintercept=rescaled_ec50)) +
    geom_text(aes(x = (rescaled_ec50+0.25), y = 1, label = paste0("EC",EC,"  ",round(Estimate,4)), group= NULL),color = '#00008b', show.legend = FALSE)
  p +
    facet_wrap(~individual, 
               labeller = labeller(group_pair_ = label_facet(data_merged$individual, paste0(" ", data_merged$Estimate))
                                   )) +
    small_font_theme +
    labs(title=paste0("Paired group: ",group_pair_,"   |   lme4 Prediction versus Raw"),
    subtitle = gsub("[[:space:]]", "", Reduce(paste, deparse((formula(model_results$models[[selected_model]])))))
    )
}