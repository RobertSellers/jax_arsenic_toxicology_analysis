harmony_read_plate_2 <- function(dir){

  if (!dir.exists(dir)){
    stop("directory parameter not found")
  }
  
  # comprehensive plate list
  temp_list <- list()
  
  # create a plate class to keep things organized
  setClass(Class="plate",
   representation(
      index="list",
      plate_df="data.frame",
      cells_df="data.frame",
      cells_header="character",
      plate_metadata="data.frame",
      plate_dir = "character"
    ))
  
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

    # generate new class and append to specified list by plate_obj_name
    new_plate_class <- new("plate",
                index=idx_df_split,
                plate_df=plate_df,
                cells_df=cell_df,
                cells_header=cells_header,
                plate_metadata=metadata_df,
                plate_dir = dir
            )
    
    # add plate object to list
    len_list <- length(temp_list) + 1
    temp_list[[plate_obj_name]] <- new_plate_class
  }
  return(temp_list)
}