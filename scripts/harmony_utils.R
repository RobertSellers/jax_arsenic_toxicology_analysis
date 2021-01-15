plate_class <- function(index, plate_df, cells_df, cells_header ,plate_metadata ,plate_dir){
  
  structure(class = "plate", list(
    # attributes
    index=index,
    plate_df=plate_df,
    cells_df=cells_df,
    cells_header=cells_header,
    plate_metadata=plate_metadata,
    plate_dir = plate_dir,
    # methods
    get_plate_name = function() paste("plate_dir was", plate_dir)
  ))
}

plate_collection_class <- function(list_of_plates){
  self.all_plates = do.call(rbind.fill, lapply(list_of_plates, function(x) {x$plate_df}))
  self.plate_features = NULL
  structure(class = "plate_collection", list(
    plates = list_of_plates,
    # data methods
    all_plates = self.all_plates,
    metadata = do.call(rbind.fill, lapply(list_of_plates, function(x) {x$plate_metadata})),
    # extration methods
    # selects only predictors + input batches
    select_features = function(batch_cols){
      if(missing(batch_cols)){
        stop("No batch columns selected.")
      }else{
        start_idx = grep("timepoint", colnames(self.all_plates))+1
        end_idx = grep("number_of_analyzed_fields", colnames(self.all_plates))-1
        # takes a vector of column names
        batch_cols <-names(self.all_plates)[(names(self.all_plates) %in% keep)]
        return (cbind(self.all_plates[, batch_cols], #batch in front
                      self.all_plates[,start_idx:end_idx]) #features in back
                )
      }
    },
    show_features = function(){
      start_idx = grep("timepoint", colnames(self.all_plates))+1
      end_idx = grep("number_of_analyzed_fields", colnames(self.all_plates))-1
      cat(colnames(cbind(self.all_plates[,start_idx:end_idx])))
    }
  ))
}

# let's custamize the summary, for example
# uses generic s3 base summary function
# consider upgrading to lexical scoping - see https://rpubs.com/mrloh/oor example
summary.plate <- function(obj){
  cat("Dimensions of plate df", dim(obj$plate_df))
} 

summary.plate_collection <- function(obj){
  cat("Contains", length(obj$plates), "plates")
}

harmony_read_plate_3 <- function(dir){
  
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
    file_list <- unlist(strsplit(evaluations, "/"))
    e <- as.character(sub('.*(?=.{1}$)', '',file_list[7], perl=T))
    plate_obj_name <- paste0('d',sub(" ", "_", gsub('\\-', '',sub("/","__",paste0(file_list[3])))),"_e",e)
    
    cur_path <- file.path(dir,file_list[3])
    
    # str determines where files should be if applicable
    index_file <- file.path(cur_path, "indexfile.txt")
    plate_file <- file.path(cur_path, paste0("Evaluation", e), "PlateResults.txt")
    
    # not dynamic...needs improvement
    cells_file <- file.path(cur_path, paste0("Evaluation", e), "Objects_Population - Nuclei Selected.txt")
    
    ###### INDEX FILE#######
    if(file.exists(index_file)){
      
      idx_df <- read_tsv(index_file) %>%
        select_if(function(x) any(!is.na(x))) %>% # removes columns with all NA
        clean_names() %>%# removes spaces and bad characters for underscores
        mutate(channel_name = gsub("[[:punct:][:blank:]]", "_", channel_name),
               channel_name = tolower(channel_name))
      
      # split dataframes based on stain value
      idx_df_split <- split(idx_df, idx_df$channel_name)
    }else{
      #nothing right now
      stop('Index data not loaded')
    }
    
    # ###### PLATE FILE#######
    # 
    if(file.exists(plate_file)){
      
      # handle metadata for plate header
      plate_header_raw <- readLines(plate_file,n = 8)
      split_header <- strsplit(plate_header_raw, split="\t")
      metadata_df <- as.data.frame(t(do.call(rbind,split_header[c(1:6)]))) %>%
        row_to_names(row_number = 1)
      
      plate_df <- read_tsv(plate_file,skip=8) %>%
        select_if(function(x) any(!is.na(x))) %>% # removes columns with all NA
        clean_names() %>%# removes spaces and bad characters for underscores
        mutate(dir = plate_obj_name)
      plate_df[is.nan(plate_df)] <- NA #change NaN to NA
      
    }else{
      stop('Plate data not loaded')
    }
    
    # ###### CELLS FILE#######
    if(file.exists(cells_file)){
      
      # read file header
      cells_header <- readLines(cells_file, n = 9)
      
      # read file data
      cell_df <- read_tsv(cells_file,skip=9) %>%
        select_if(function(x) any(!is.na(x))) %>% # removes columns with all NA
        clean_names() %>%# removes spaces and bad characters for underscores
        mutate(plate_obj_name = plate_obj_name)
    }else{
      warning('Cells data not loaded')
      cells_header <- ''
      cell_df <- data.frame()
    }
    
    new_plate_class <- plate_class(idx_df_split, plate_df,cell_df,cells_header,metadata_df,dir)
    # add plate object to list
    len_list <- length(temp_list) + 1
    temp_list[[plate_obj_name]] <- new_plate_class
  }
  plate_collection <- plate_collection_class(temp_list)
  return(plate_collection)
}