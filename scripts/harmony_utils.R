
# this script generates a 4 component harmony_collection class
# requires dplyr, readr, checkmate, plyr, janitor packages
require(dplyr)

# primary function
harmony_create_collection <- function(
  dir, 
  overwrite = FALSE, 
  temp_dir = NA,
  custom_rename = FALSE
  ){
  options(readr.num_columns = 0)
  # data cache handling
  # incomplete - need parameter for overwrite
  dir_name <- unlist(strsplit(dir, "/"))[length(unlist(strsplit(dir, "/")))]
  temp_data_storage <- paste0(temp_dir,"/",dir_name,'.RData')
  # look for r data
  if (file.exists(temp_data_storage) && overwrite == FALSE){
    cat(paste0("\nLoading existing plate collection from temp_files"))
    plate_collection <- readRDS(temp_data_storage)
  } else { 
    start.time <- Sys.time()
    
    if (!dir.exists(dir)){
      stop("directory parameter not found")
    }
    
    # comprehensive plate list
    temp_list <- list()
    i <- 1 # for progress 
    child_dir <- list.dirs(path = dir, full.names = TRUE, recursive = FALSE)
    # construct load list
    # takes last evaluation found
    for (child in child_dir){

      # file naming / loop through available evaluations
      evaluations <- list.dirs(path = child, full.names = TRUE, recursive = TRUE)
      file_list <- unlist(strsplit(evaluations[length(sort(evaluations))], "/"))
      child_name <- unlist(strsplit(child, "/"))[length(unlist(strsplit(child, "/")))]
      # extract evalution # from string
      e <- as.character(readr::parse_number(file_list[length(file_list)]))
      plate_obj_name <- paste0("d",sub(" ", "_", gsub('\\-', '',sub("/","__",child_name))),"_e",e)
      
      cur_path <- file.path(dir,file_list[length(file_list)-1])
      
      ###### INDEX FILE#######
      idx_df <- data.frame()
      try({
        index_file <- file.path(cur_path, "indexfile.txt")
        if(file.exists(index_file)){
          idx_df <- readr::read_tsv(index_file, progress = readr::show_progress()) %>%
            select_if(function(x) any(!is.na(x))) %>% # removes columns with all NA
            janitor::clean_names() %>% # removes spaces and bad characters for underscores
            mutate(channel_name = tolower(gsub("[[:punct:][:blank:]]", "_", channel_name)),
                   dir = plate_obj_name)
          
          # split dataframes based on stain value
          #idx_df_split <- split(idx_df, idx_df$channel_name)
        }
        else{
          #nothing right now
          warning('Index data not loaded')
        }
      })
      
      # ###### PLATE FILE#######
      plate_df <- data.frame()
      try({
        plate_file <- file.path(cur_path, paste0("Evaluation", e), "PlateResults.txt")
        
        if(file.exists(plate_file)){
          
          # handle metadata for plate header
          plate_header_raw <- readLines(plate_file,n = 8)
          split_header <- strsplit(plate_header_raw, split="\t")
          metadata_df <- as.data.frame(t(do.call(rbind,split_header[c(1:6)]))) %>%
            janitor::row_to_names(row_number = 1)
          
          plate_df <- readr::read_tsv(plate_file,skip=8) %>%
            select_if(function(x) any(!is.na(x))) %>% # removes columns with all NA
            janitor::clean_names() %>%# removes spaces and bad characters for underscores
            mutate(dir = plate_obj_name)
          
        }else{
          warning('Plate data not loaded')
        }
      })
      
      # ###### CELLS FILE#######
      cell_df <- data.frame()
      try({
        cells_file <- list.files(file.path(cur_path, paste0("Evaluation", e)), full.names = TRUE, pattern = "Objects_Population")[1]
        if(file.exists(cells_file)){
          
          # read file header
          cells_header_raw <- readLines(cells_file, n = 9)
          split_header <- strsplit(cells_header_raw, split="\t")
          metadata_df <- as.data.frame(t(do.call(rbind,split_header[c(1:6)]))) %>%
            janitor::row_to_names(row_number = 1)
          # read file data
          cell_df <- readr::read_tsv(cells_file,skip=9) %>%
            select_if(function(x) any(!is.na(x))) %>% # removes columns with all NA
            janitor::clean_names() %>%# removes spaces and bad characters for underscores
            mutate(plate_obj_name = plate_obj_name)
          
        }else{
          warning('Cells data not loaded')
          cells_header <- ''
          cell_df <- data.frame()
        }
      })

      new_plate_class <- plate_class(idx_df, plate_df, cell_df, metadata_df, dir)
      # add plate object to list
      len_list <- length(temp_list) + 1
      temp_list[[plate_obj_name]] <- new_plate_class
    }
    # analyze features and add to plate collection parameter class
    plate_collection <- plate_collection_class(temp_list, custom_rename)
    
    end.time <- Sys.time()
    time.taken <- round(end.time - start.time,2)
    print(time.taken)
    cat(paste0("\nSaving plate collection to temp_files"))
    tryCatch(
      expr = {
        saveRDS(plate_collection, paste0(temp_dir,"/",dir_name,'.RData'))
      },
      error = function(e){ 
        cat(paste("error - temp_dir not found / unsaved"))
      }
    )
  }
  return (plate_collection)
}

plate_class <- function(index, plate_df, cells_df, plate_metadata, plate_dir){

  structure(class = "plate", list(
    # attributes
    index = index,
    plate_df = plate_df,
    cells_df = cells_df,
    plate_metadata = plate_metadata,
    plate_dir = plate_dir,
    cells_size = nrow(cells_df),
    plate_size = nrow(plate_df),
    # methods
    get_plate_name = function() paste("plate_dir was", plate_dir)
  ))
  
}

# defining s3
# consult with http://adv-r.had.co.nz/OO-essentials.html if necessary
plate_collection_class <- function(list_of_plates, custom_rename){
  
  self.all_plates = do.call(plyr::rbind.fill, lapply(list_of_plates, function(x) {x$plate_df}))
  
  # bulk rename function as input parameter
  if (custom_rename) self.all_plates = colname_condense(self.all_plates)
  
  self.index = do.call(plyr::rbind.fill, lapply(list_of_plates, function(x) {x$index}))
  self.all_cells = do.call(plyr::rbind.fill, lapply(list_of_plates, function(x) {x$cells_df}))
  self.feature_df = feature_columns_analyze(self.all_plates)

  structure(class = "plate_collection", list(
    # plates = list_of_plates,
    # data methods
    all_plates = self.all_plates,
    all_index = self.index,
    all_cells = self.all_cells,
    features_df = self.feature_df
  ))
}

# s3 / methods
harmony_qa <- function(x) UseMethod("harmony_qa")
harmony_qa.plate_collection <- function(x) {
  print("Analyzing all_plates")
  start_idx = grep("timepoint", colnames(x$all_plates)) + 1
  end_idx = grep("number_of_analyzed_fields", colnames(x$all_plates)) - 1
  print(paste("Total predictor columns:",as.character(end_idx - start_idx)))
  print(paste("Total wells:",nrow(x$all_plates)))
  
  if(nrow(x$all_cells)>0){
    print("--------------------")
    print("Analyzing all_cells")
    start_idx = grep("timepoint", colnames(x$all_cells)) + 1
    end_idx = grep("plate_obj_name", colnames(x$all_cells)) - 1
    print(paste("Total predictor columns:",as.character(end_idx - start_idx)))
    print(paste("Total cells:",nrow(x$all_cells)))
  }
}

# uses generic s3 base summary function
# consider upgrading to lexical scoping - see https://rpubs.com/mrloh/oor example
summary.plate <- function(obj){
  cat("Dimensions of plate df", dim(obj$plate_df), "\n")
} 

# custom summary method
summary.plate_collection <- function(obj){

  cat("Total number of plates:", length(obj$plates),"\n")
  cat("Total wells:", nrow(obj$all_plates), "\n")
  cat("Wells per plate:", max(obj$all_plates$row*obj$all_plates$column),"\n") # get mode from each run?
  cat("Total number of individuals:", length(unique(obj$all_plates$individual)),"\n")
  cat("Total objects In Focus/Out of Focus:",sum(obj$all_plates$in_focus_number_of_objects,na.rm = TRUE)," / ", 
  sum(obj$all_plates$out_focus_number_of_objects,na.rm = TRUE),"\n")
}

# construction of summary stats per feature
# these aren't observable upon edits to the plate_df data
feature_columns_analyze <- function(df){

  features_df <- NULL
  # the data expects these two columns
  start_idx = grep("timepoint", colnames(df)) + 1
  end_idx = grep("number_of_analyzed_fields", colnames(df)) - 1
  for (feature in colnames(df[,start_idx:end_idx])){
    
    temp_feature <- df %>%
      select_(feature)
    f <- temp_feature[,1]
    
    name_ <- feature
    total_ <- nrow(temp_feature)
    na_count_ <- sum(is.na(f))
    median_ <- signif(as.numeric(median(f,na.rm=TRUE)), 3)
    max_ <- signif(as.numeric(max(f,na.rm=TRUE)), 3)
    min_ <- signif(as.numeric(min(f,na.rm=TRUE)), 3)
    integer_ <- checkmate::testInteger(f) 

    features_df <- rbind(features_df, data.frame(
      name_,
      integer_,
      total_,
      na_count_,
      median_,
      min_,
      max_
    ))
  }
  return (features_df)
}

# clean up column names
# this is a unique circumstance
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

# subset columns with predictors
# expects these specific columns
query_features <- function(data){
  start_idx = grep("timepoint", colnames(data)) + 1
  end_idx = grep("number_of_analyzed_fields", colnames(data)) - 1
  all_features <- colnames(data[,start_idx:end_idx])
  return(all_features)
}