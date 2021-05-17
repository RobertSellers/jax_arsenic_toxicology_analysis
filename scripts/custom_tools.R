# rankZ transform function definition
rankZ <- function(x) {
  x <- rank(x, na.last = "keep") / (length(x) - sum(is.na(x)) + 1)
  return(qnorm(x))
}

# builds boolean vector/column for integer vs float
testInteger <- function(x){
  test <- all.equal(x, as.integer(x), check.attributes = FALSE)
  if(test == TRUE){ return(TRUE) }
  else { return(FALSE) }
}

LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))

# grep function that removes a pattern
remove_sd_vars <- function(data){
  # remove from all_plates
  before_cols <- ncol(data$all_plates)
  data$all_plates <- dplyr::select(data$all_plates, -dplyr::contains("stdev"))
  print(paste("removing",(before_cols-ncol(data$all_plates)), "/", before_cols,"columns"))
  # remove from features_df
  data$features_df <- data$features_df[!grepl("stdev", data$features_df$name_),]
  return (data)
}

impute.mean <- function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))

# ensure all NA values are the same
is.nan.data.frame <- function(x){
  do.call(cbind, lapply(x, is.nan))
}

# pretty print percent
percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}

# for use in dplyr %>% select(where(not_all_na))
not_all_na <- function(x) any(!is.na(x))

# for use in dplyr %>% select(where(not_any_na))
not_any_na <- function(x) all(!is.na(x))

# for pipes
`%notin%` <- Negate(`%in%`)

VLookup <- function(this, data, key, value) {
  m <- match(this, data[[key]])
  data[[value]][m]
}

label_facet <- function(original_var, custom_name){
  lev <- levels(as.factor(original_var))
  lab <- paste0(custom_name, ": ", lev)
  names(lab) <- lev
  return(lab)  
}

merge.all <- function(x, ..., by = "row.names") {
  L <- list(...)
  for (i in seq_along(L)) {
    x <- merge(x, L[[i]], by = by)
    rownames(x) <- x$Row.names
    x$Row.names <- NULL
  }
  return(x)
}

sample_n_groups = function(grouped_df, size, replace = FALSE, weight=NULL) {
  grp_var <- grouped_df %>% 
    groups %>%
    unlist %>% 
    as.character
  random_grp <- grouped_df %>% 
    summarise() %>% 
    sample_n(size, replace, weight) %>% 
    mutate(unique_id = 1:NROW(.))
  grouped_df %>% 
    right_join(random_grp, by=grp_var) %>% 
    group_by_(grp_var) 
}

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

# rescale_feature_df_0_100 <- function(data, covariates){
#   all_features <- query_features(data)
#   data[,c(all_features)] <- lapply(data[,c(all_features)], function(x) as.numeric(as.character(x)))
#   data[,c(all_features)] <- lapply(data[,c(all_features)], function(x) scales::rescale(x, to = c(0.001, 99.999)))
#   return (list(
#     "rescaled_data" = data[,c(covariates,all_features)],
#     "parameters" = list(
#       "max_n" = length(all_features)
#       )
#     )
#   )
# }

# ######## 2 #########
# mean_correlation_by_doseresponse <- function(data, group_id = 'mouse', remove_outliers = TRUE){
  
#   ############ extract dose pairs for loop #############
#   f <- function(x) {
#     s <- seq(2, length(x), 1)
#     paste(x[s-1], x[s], sep=",")
#   }
  
#   dose_pairs <- as.data.frame(do.call(rbind, strsplit(f(sort(unique(data$rescaled_data$dose))), ",")))
#   colnames(dose_pairs) <- c("low","high")
#   ######################################################
  
#   mean_corr_dosepair_df_list<-mean_stdev_df_list<-mean_corr_matrix_list<-list()
    
#   # loop through dose pairs
#   for (i in 1:nrow(dose_pairs)){
#     low <- dose_pairs[i, "low"]
#     high  <- dose_pairs[i, "high"]
#     print(paste("processing response change at", low,":",high, "at",group_id,"level"))
    
#     sub_df_mean_diff_per_dose_response <- data$rescaled_data %>%
#       dplyr::filter(dose %in% c(low,high)) %>% # select only values inside low and high
#       separate(mouse, c("strain","replicate")) %>% # extract specific - remove as is unnecessary
#       group_by(strain, dose) %>%  # group by strain and dose - also unnecessary code
#       summarise_each(funs(mean), -c(replicate)) %>% # mean of each strain per dose response change
#       group_by(strain) %>%# now just group by strain
#       #filter_at(vars(-c(strain, dose)), all_vars(!(abs(. - median(.)) > 2*sd(.)))) %>% # removing outliers not advisable
#       summarise_at(vars(-dose),diff) # calculate difference in mean values
#     sub_df_stdev_per_dose_response <- data$rescaled_data %>%
#       dplyr::filter(dose %in% c(low,high)) %>% # select only values inside low and high
#       separate(mouse, c("strain","replicate")) %>% # fix
#       group_by(strain) %>%  # group by strain 
#       summarise_at(vars(-c(replicate, dose)),diff) %>%# get difference per dose pair within strain group
#       group_by(strain) %>%
#       summarise_each(funs(sd))
    
#     # keep the first column 
#     names <-  sub_df_mean_diff_per_dose_response$strain
    
#     # Transpose everything other than the first column
#     mean_T <- as.data.frame(as.matrix(t(sub_df_mean_diff_per_dose_response[,-1])))
#     stdev_T <- as.data.frame(as.matrix(t(sub_df_stdev_per_dose_response[,-1])))
    
#     # Assign first column as the column names of the transposed dataframe
#     colnames(stdev_T) <- colnames(mean_T) <- names
#     res2<-Hmisc::rcorr(as.matrix(t(mean_T)), type="spearman")
#     ut <- upper.tri(res2$r)
    
#     features_a <- rownames(res2$r)[row(res2$r)[ut]]
#     features_b <- rownames(res2$r)[col(res2$r)[ut]]
#     correlations <- (res2$r)[ut]
    
#     flat_matrix_to_df <- data.frame(
#       feature_a = features_a,
#       feature_b = features_b,
#       cor  = correlations
#       )
#     colnames(flat_matrix_to_df)[3] <- paste0("cor_",low,"_",high)
    
#     # append DR matrix to sub list
#     mean_corr_matrix_list[[paste0("interval_",low,"_",high)]] <- res2$r
#     mean_stdev_df_list[[paste0("interval_",low,"_",high)]] <- as.matrix(stdev_T)
#     mean_corr_dosepair_df_list[[paste0("interval_",low,"_",high)]] <- flat_matrix_to_df
#   }
#   # return all 3 sublists
#   return (list("mean_corr_matrix_list" = mean_corr_matrix_list,
#                "mean_stdev_df_list" = mean_stdev_df_list,
#                "mean_corr_dosepair_df_list" = mean_corr_dosepair_df_list,
#                "parameters" = list(
#                   "max_n" = data$parameters$max_n
#                )
#                ))
# }

log_skew_transform <- function(response_var) {
  skew_val <- e1071::skewness(response_var)
  skew.score <- function(c, x)
    (e1071::skewness(log(x + c))) ^ 2
  if (!is.na(skew_val)) {
    print(paste0("performing log transformation with skew:", skew_val))
    best.c <- optimise(skew.score, c(0, 1), x = response_var)$minimum
    response_var <- log(response_var + best.c)
  } else{
    print("No optimization possible")
    # to be continued
    response_var <- response_var
    #df$response <- log(df$response)
  }
  return(response_var)
}

# QA function to generate a faceted boxplot per each plate pair
boxplot_pairs <- function(y, data, filename = NA){
  require(ggplot2)
  
  small_font_theme <- ggplot2::theme(
    axis.text.x = element_text(colour="grey20", size=6, angle=90, hjust=.5, vjust=.5),
    axis.text.y = element_text(colour="grey20", size=6),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 5)
  )
  if(!is.na(filename)){
    pdf(file = filename, height = 10, width = 10)
  }
  p <- ggplot(data) +
    theme(legend.position = "none") +
    geom_boxplot(aes_string(x = "individual", y = y, colour = "bin_id"), 
    show.legend = FALSE,
    outlier.size = -1) +
    ylab(y) 
  # browser()
  fr <- p +
    facet_wrap(~group_pair, scales = "free_x", 
      labeller = labeller(group_pair = label_facet(data$group_pair, "paired"))
    ) +
    small_font_theme
  if(!is.na(filename)){
    print(fr)
    dev.off()
  }else{
    fr
  }
}

group.center <- function(var,grp) {
  return(var-tapply(var,grp,mean,na.rm=T)[grp])
}

# adds binary field for plotting and pairing
# very slow -- could use update
iterate_bin_id <- function(df){
  df <- df %>%
    tibble::add_column(bin_id = NA)
  for (i in 1:nrow(df)){
    if (i != 1){
      if (df[i,]$dir == df[i-1,]$dir){
        df[i,]$bin_id <- df[i-1,]$bin_id
      }else{
        if(df[i-1,]$bin_id == 1){
          df[i,]$bin_id <- 0
        }else{
          df[i,]$bin_id <- 1
        }
      }
    }else{
      df[i,]$bin_id <- 1
    }
  }
  df$bin_id <- as.character(df$bin_id)
  return (df)
}

# batch correction lme4 auto script
batch_correct <- function(formula_,feature, data, r.eff, target_col, debug = FALSE){
  require(lme4)
  if (debug) browser()
  pheno_raw <- as.numeric(unlist(data[,c(feature)]))
  mn.pheno <- mean(pheno_raw, na.rm=TRUE)
  sd.pheno <- sd(pheno_raw, na.rm=TRUE)
  data$pheno <- (pheno_raw - mn.pheno) / sd.pheno
  fit <- lmer(formula_, data, REML=FALSE)
  random_effects <- ranef(fit)
  effects <- random_effects[[r.eff]]
  batch.ef <- effects[as.character(data[,c(target_col)][[1]]),]
  pheno_sub <- data$pheno - batch.ef
  pheno_adj <- (sd.pheno * pheno_sub) + mn.pheno
  return (pheno_adj)
}
colname_condense <- function(df){
  df <- df %>%
    rename_at(.vars = vars(ends_with("_per_well")),
              .funs = funs(sub("[_]per_well$", "", .))) %>%
    rename_at(.vars = vars(ends_with("_cells_cells")),
              .funs = funs(sub("[_]cells_cells$", "cells", .))) %>%
    rename_at(.vars = vars(contains("1_px")),
              .funs = funs(sub("1_px", "1px", .))) %>%
    rename_at(.vars = vars(contains("hoechst_33342")),
            .funs = funs(sub("hoechst_33342", "hoechst", .))) %>%
    rename_at(.vars = vars(contains("alexa_488")),
            .funs = funs(sub("alexa_488", "alexa488", .))) %>%
    rename_at(.vars = vars(contains("mito_tracker_deep_red")),
            .funs = funs(sub("mito_tracker_deep_red", "mitotrackerdeepred", .))) %>%
    rename_at(.vars = vars(contains("_number_of_objects")),
              .funs = funs(sub("_number_of_objects", "_numberofobjects", .))) %>%
    rename_at(.vars = vars(contains("std_dev")),
            .funs = funs(sub("std_dev", "stdev", .))) %>%
    rename_at(.vars = vars(contains("_pos_")),
            .funs = funs(sub("[_]pos[_]", "_positive_", .))) %>%
    rename_at(.vars = vars(contains("_neg_")),
            .funs = funs(sub("[_]neg[_]", "_negative_", .))) %>%
    rename_at(.vars = vars(contains("out_focus")),
            .funs = funs(sub("out_focus", "outfocus", .))) %>%
    rename_at(.vars = vars(contains("in_focus")),
          .funs = funs(sub("in_focus", "infocus", .))) %>%
    rename_at(.vars = vars(contains("number_of_spots")),
        .funs = funs(sub("number_of_spots", "numberofspots", .))) %>%
    rename_at(.vars = vars(contains("h2ax_positive")),
        .funs = funs(sub("h2ax_positive", "h2axpositive", .))) %>%
    rename_at(.vars = vars(contains("h2ax_negative")),
        .funs = funs(sub("h2ax_negative", "h2axnegative", .))) %>%
    distinct()
  return (df)
}

wrapit <- function(text) {
  wtext <- paste(strwrap(text,width=40),collapse=" \n ")
  return(wtext)
}