# not currently used
dependent_variable_analysis_dr <- function (df, ignore, response_var){
  # attempting https://stackoverflow.com/questions/12234248/printing-multiple-ggplots-into-a-single-pdf-multiple-plots-per-page
  # Plot separate ggplot figures in a loop.
  set.seed(11)
  # select dependent variables
  df_predictor_cols <- names(df)[!(names(df) %in% ignore)]
  
  # Make list of variable names to loop over.
  distribution_plot_list = list()
  drc_plot_list = list()
  
  for (i in 1:length(df_predictor_cols)){
    
      distribution_plot = ggplot(
        data = df, 
        aes_string(x=response_var, y=df_predictor_cols[i])) +
        geom_boxplot(width=0.1) +
        #geom_point(size=3, aes(colour=sex)) + 
        ggtitle(paste0(response_var," with\n",df_predictor_cols[i]))
      
      distribution_plot_list[[i]] = distribution_plot
      drc_plot_list[[i]] = distribution_plot
  }
  # Make plots.
  pdf_name <- "distributions.pdf"
  pdf(pdf_name)
  for (i in 1:length(df_predictor_cols)) {
    multi_page <- ggarrange(distribution_plot_list[[i]],drc_plot_list[[i]], nrow=2, ncol=1)
    print(multi_page)
  }
  dev.off()
}

testInteger <- function(x){
  test <- all.equal(x, as.integer(x), check.attributes = FALSE)
  if(test == TRUE){ return(TRUE) }
  else { return(FALSE) }
}