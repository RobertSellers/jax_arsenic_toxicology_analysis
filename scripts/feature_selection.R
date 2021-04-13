######## 1 #########
rescale_feature_df_0_100 <- function(data, covariates){
  all_features <- query_features(data)
  data[,c(all_features)] <- lapply(data[,c(all_features)], function(x) as.numeric(as.character(x)))
  data[,c(all_features)] <- lapply(data[,c(all_features)], function(x) scales::rescale(x, to = c(0.001, 99.999), from = range(x, na.rm = TRUE, finite = TRUE)))
  return (list(
    "rescaled_data" = data[,c(covariates,all_features)],
    "parameters" = list(
      "max_n" = length(all_features)
      )
    )
  )
}

######## 2 #########
mean_correlation_by_doseresponse <- function(data, group_id = 'mouse'){
  
  ############ extract dose pairs for loop #############
  f <- function(x) {
    s <- seq(2, length(x), 1)
    paste(x[s-1], x[s], sep=",")
  }
  
  dose_pairs <- as.data.frame(do.call(rbind, strsplit(f(sort(unique(data$rescaled_data$dose))), ",")))
  colnames(dose_pairs) <- c("low","high")
  ######################################################
  
  mean_corr_dosepair_df_list<-mean_stdev_df_list<-mean_corr_matrix_list<-list()
    
  # loop through dose pairs
  for (i in 1:nrow(dose_pairs)){
    low <- dose_pairs[i, "low"]
    high  <- dose_pairs[i, "high"]
    print(paste("processing response change at", low,":",high, "at",group_id,"level"))
    
    sub_df_mean_diff_per_dose_response <- data$rescaled_data %>%
      dplyr::filter(dose %in% c(low,high)) %>% # select only values inside low and high
      separate(mouse, c("strain","replicate")) %>% # extract specific - remove as is unnecessary
      group_by(strain, dose) %>%  # group by strain and dose - also unnecessary code
      summarise_each(funs(mean), -c(replicate)) %>% # mean of each strain per dose response change
      group_by(strain) %>%# now just group by strain
      #filter_at(vars(-c(strain, dose)), all_vars(!(abs(. - median(.)) > 2*sd(.)))) %>% # removing outliers not advisable
      summarise_at(vars(-dose),diff) # calculate difference in mean values
    sub_df_stdev_per_dose_response <- data$rescaled_data %>%
      dplyr::filter(dose %in% c(low,high)) %>% # select only values inside low and high
      separate(mouse, c("strain","replicate")) %>% # fix
      group_by(strain) %>%  # group by strain 
      summarise_at(vars(-c(replicate, dose)),diff) %>%# get difference per dose pair within strain group
      group_by(strain) %>%
      summarise_each(funs(sd))
    
    # keep the first column 
    names <-  sub_df_mean_diff_per_dose_response$strain
    
    # Transpose everything other than the first column
    mean_T <- as.data.frame(as.matrix(t(sub_df_mean_diff_per_dose_response[,-1])))
    stdev_T <- as.data.frame(as.matrix(t(sub_df_stdev_per_dose_response[,-1])))
    
    # Assign first column as the column names of the transposed dataframe
    colnames(stdev_T) <- colnames(mean_T) <- names
    res2<-Hmisc::rcorr(as.matrix(t(mean_T)), type="spearman")
    ut <- upper.tri(res2$r)
    
    features_a <- rownames(res2$r)[row(res2$r)[ut]]
    features_b <- rownames(res2$r)[col(res2$r)[ut]]
    correlations <- (res2$r)[ut]
    
    flat_matrix_to_df <- data.frame(
      feature_a = features_a,
      feature_b = features_b,
      cor  = correlations
      )
    colnames(flat_matrix_to_df)[3] <- paste0("cor_",low,"_",high)
    
    # append DR matrix to sub list
    mean_corr_matrix_list[[paste0("interval_",low,"_",high)]] <- res2$r
    mean_stdev_df_list[[paste0("interval_",low,"_",high)]] <- as.matrix(stdev_T)
    mean_corr_dosepair_df_list[[paste0("interval_",low,"_",high)]] <- flat_matrix_to_df
  }
  # return all 3 sublists
  return (list("mean_corr_matrix_list" = mean_corr_matrix_list,
               "mean_stdev_df_list" = mean_stdev_df_list,
               "mean_corr_dosepair_df_list" = mean_corr_dosepair_df_list,
               "parameters" = list(
                  "max_n" = data$parameters$max_n
               )
               ))
}
######## 3 #########

filter_exceedances <- function(list_data, n, corr_cut_off, rank_by){
  
  print(paste("top",n,"selected from highest", rank_by))
  
  ########## 3a ###########
  vals <- feature_exceedance_extraction(dr_list_object = list_data, corr_cut_off = corr_cut_off, n = n)
  vals$name <- gsub("\\..*","",rownames(vals))
  
  ########## 3b ###########
  sub_stdev <- extract_mean_stdev(list_data$mean_stdev_df_list, n=n)
  
  sub_stdev$name <- rownames(sub_stdev)
  
  intersecting_features <- intersect(sub_stdev$name,vals$name)
  remove <- c()
  keep <- c()
  
  # start @ lowest stdev 
  for (i in 1:nrow(sub_stdev)){
    name <- gsub("\\..*","",sub_stdev[i,]$name)
    if (name %in% intersecting_features){
        f_corr_df <- vals[vals[,'name']==name,]
        f_corr_df_mirror <- vals[(vals[,'col'] ==  f_corr_df$row),]
        f_corr_all <- cbind(f_corr_df,f_corr_df_mirror)[,c(4,8,3)] %>% 
          left_join(sub_stdev, by = c("name" = "name")) %>%
          left_join(sub_stdev, by = c("name.1" = "name")) %>%
          dplyr::select(dplyr::contains(c("name","correlation","sdmeans"))) %>% 
          arrange(desc(abs(correlation)))

        # loop through every correlation pair and select highest SD between the two
        for (r in 1:nrow(f_corr_all)){
            # compare first SD against the second
            # and confirm whether it is already slated for removal
            if ((f_corr_all$sdmeans.x[r] >= f_corr_all$sdmeans.y[r]) && 
                (f_corr_all$name.1[r] %notin% remove)){
              if(f_corr_all$sdmeans.x[r] == f_corr_all$sdmeans.y[r]){
                # handle identical features
                print(paste(f_corr_all$name[r],"is identical to",f_corr_all$name.1[r]))
                #browser()
              }else{
                # remove the second pair
                remove <- c(remove,f_corr_all$name.1[r])
              }
            }else{
              if(f_corr_all$name[r] %notin% remove){
                # remove the first pair
                remove <- c(remove,f_corr_all$name[r])

              }
            }
        }
    }
  }
  
  ### prepare return object ###
  keep <- sub_stdev[sub_stdev$name %notin% remove,]
  remove <- sub_stdev[sub_stdev$name %in% remove,]
  
  #### add stdev data ###
  stdevs_subset <- subset(sub_stdev, rownames(sub_stdev) %in% keep$name)
  
  print(paste("salvaging", nrow(keep),"/",n,"features"))
  return (list(
    "keep"=keep,
    "remove"=remove,
    "stdev_subset" = stdevs_subset, 
    "parameters"= list(
      "max_n" = list_data$parameters$max_n,
      "target_n" = n, 
      "corr_threshold" = corr_cut_off,
      "rank_by" =rank_by)
    )
  )
}

########## 3b, 3a-1 ###########
extract_mean_stdev <- function(stdev_list, n = 'all'){
  # reduce all matrices into a single dataframe of average values
  sd_df <- data.frame(
    sd_0_0.01=apply(stdev_list$interval_0_0.01,1, sd, na.rm = TRUE),
    sd_0.01_0.1=apply(stdev_list$interval_0.01_0.1,1, sd, na.rm = TRUE),
    sd_0.1_0.75=apply(stdev_list$interval_0.1_0.75,1, sd, na.rm = TRUE),
    sd_0.75_1=apply(stdev_list$interval_0.75_1,1, sd, na.rm = TRUE),
    sd_1_1.25=apply(stdev_list$interval_1_1.25,1, sd, na.rm = TRUE),
    sd_1.25_2=apply(stdev_list$interval_1.25_2,1, sd, na.rm = TRUE),
    sd_2_5=apply(stdev_list$interval_2_5,1, sd, na.rm = TRUE)
  )
  sd_df$sdmeans <- rowMeans(sd_df)
  return (sd_df %>% slice_max(sdmeans, n = n))
}

########## 4 ###########
extract_mean_corr <- function(corr_list){
  
  corr_reduced_matrix <- Reduce("+",  lapply(corr_list, function(x) replace(x, is.na(x), 0))) / length(lapply(corr_list, Negate(is.na)))

  #corr_reduced_df$corrmeans <- rowMeans(corr_reduced_df)
  return (corr_reduced_matrix)
  
}

########## 3a ###########
feature_exceedance_extraction <- function(dr_list_object, corr_cut_off, n){

  ########## 3a-1 ###########
  sub_stdev <- extract_mean_stdev(dr_list_object$mean_stdev_df_list, n)
  
  features <- data.frame("feature" = rownames(sub_stdev), "mean_sd" = sub_stdev$sdmeans)
  
    ########## 3a-2 ###########
    vals <- subset_corr_matrix(
      avg_corr_mat =  dr_list_object$mean_corr_matrix_list,
      avg_std_mat = sub_stdev,
      all_features = features,
      corr_cut_off = corr_cut_off
      )
    
    return(vals)
}
########## 3a-2 ###########
subset_corr_matrix <- function(avg_corr_mat, avg_std_mat, all_features, corr_cut_off){
  
    corr_df_mean_all <- as.data.frame(Reduce("+",  lapply(avg_corr_mat, function(x) replace(x, is.na(x), 0))) / length(lapply(avg_corr_mat, Negate(is.na))))
    
    # subset on all_features
    corr_df_mean_subset <- corr_df_mean_all[
      rownames(corr_df_mean_all) %in% all_features$feature, 
      colnames(corr_df_mean_all) %in% all_features$feature
      ] 
    
    m <- as.matrix(corr_df_mean_subset)
    
    correlations <- which( abs(m) > corr_cut_off, arr.ind=T )
    vals <- correlations[correlations[, 1] != correlations[, 2], ]
    vals <- as.data.frame(cbind(vals,correlation=m[vals]))
    
    print(paste("identified",length(unique(vals$row)),"/",length(all_features$feature),"with one or more dose pair exceeding",corr_cut_off ))
  return(vals)
}

########## 6 ###########
corr_plot_w_stdevs <- function(m, e_results, stds){
  
    # extract values from filter_exceedance results
    target_n <- e_results$parameters$target_n
    corrv <- e_results$parameters$corr_threshold
    salvaged <- e_results$keep
    stds <- e_results$stdev_subset
    
    par(cex.main=0.75)
    colrange <- c("grey",colorRampPalette(c("orange", "yellow"))(3),colorRampPalette(c("cyan", "blue"))(3), "grey")
    
    gplots::heatmap.2( m, 
                     Rowv=FALSE, 
                     Colv=FALSE, 
                     symm = TRUE,
                     main = paste0("corr value range < ",as.character(corrv), ", n = ", nrow(salvaged), "/", target_n),
                     dendrogram='none', 
                     cellnote=round(m,2), 
                     notecol="black",  
                     notecex = 0.75,
                     cexRow = 1 / log10(nrow(m)) - 0.3,
                     cexCol = 1 / log10(ncol(m)) - 0.3,
                     srtCol = 345,  
                     adjCol = c(0, 1),
                     labCol = paste(toupper(letters[1:nrow(m)])," | ",rownames(m)),
                     labRow =   paste(toupper(letters[1:nrow(m)])," | ",paste0("mean stdev: ", round(stds$sdmeans,2))),
                     col =  colrange,
                     trace='none', 
                     key=FALSE,
                    breaks = c(-1,-corrv,(-corrv/2),(-corrv/3),0,(corrv/2),(corrv/3),corrv,1),
                       lwid = c(0.2,20),
                       lhei = c(4,15),
                     margins = c(5,8))
    
} 