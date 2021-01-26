plate_class <- function(index, plate_df, cells_df, plate_metadata, plate_dir){

  structure(class = "plate", list(
    # attributes
    index=index,
    plate_df=plate_df,
    cells_df=cells_df,
    plate_metadata=plate_metadata,
    plate_dir = plate_dir,
    cells_size = nrow(cells_df),
    plate_size = nrow(plate_df),
    # methods
    get_plate_name = function() paste("plate_dir was", plate_dir)
  ))
  
}

# working on qa for dataset
plate_qa <- function(df){
  #start_idx = grep("timepoint", colnames(x$all_cells)) + 1
  #end_idx = grep("plate_obj_name", colnames(x$all_cells)) - 1
  pct_na_df = as.data.frame(colMeans(is.na(df)))
  pct_na_df$field <- rownames(pct_na_df)
  rownames(pct_na_df) <- NULL
  pct_na_df<-pct_na_df[,c(2,1)] # reorderv
  colnames(pct_na_df)[2] <- "pct_na"
  return(pct_na_df)
}

# defining s3
# consult with http://adv-r.had.co.nz/OO-essentials.html if necessary
plate_collection_class <- function(list_of_plates){
  
  self.all_plates = do.call(rbind.fill, lapply(list_of_plates, function(x) {x$plate_df}))
  self.index = do.call(rbind.fill, lapply(list_of_plates, function(x) {x$index}))
  self.all_cells = do.call(rbind.fill, lapply(list_of_plates, function(x) {x$cells_df}))
  self.plates_qa = plate_qa(self.all_plates)
  
  structure(class = "plate_collection", list(
    plates = list_of_plates,
    # data methods
    all_plates = self.all_plates,
    all_index = self.index,
    all_cells = self.all_cells,
    plates_qa = self.plates_qa,
    # metadata = list(
    #   cell_data_size = nrow(self.all_cells),
    #   plate_data_size = nrow(self.all_plates)
    # ),
    # extration methods
    # selects only predictors + input batches
    select_features = function(keep){
      if(missing(keep)){
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
    analyze_indexfile = function(){
      # outputs image "completeness"
      # create index class and add this to summary?
      cbind(
        count(self.all_index,"dir")[1],
        count(self.all_index,"dir")[2]/
          max(self.all_index$channel_id)/max(self.all_index$field)/96)
    },
    effect_analysis = function(group_vars, reduce){
      start.time <- Sys.time()
      temp <- self.all_plates
      if(exists(reduce)){
        temp <- temp[temp[reduce[1]] == reduce[2],]
      }
      start_idx = grep("timepoint", colnames(self.all_plates))+1
      end_idx = grep("number_of_analyzed_fields", colnames(self.all_plates))-1
      summary_df <- self.all_plates %>% 
        #dplyr::filter_at(vars(one_of(group_vars)),any_vars(!is.na(.))) %>%
        dplyr::group_by_at(vars(one_of(group_vars))) %>% 
        dplyr::summarise_at(vars(start_idx:end_idx),list(
          mean = ~ mean(as.numeric(.), na.rm = TRUE), 
          median = ~ median(as.numeric(.), na.rm = TRUE),
          min= ~ min(as.numeric(.), na.rm = TRUE), 
          max= ~ max(as.numeric(.), na.rm = TRUE), 
          sd = ~ sd(as.numeric(.), na.rm = TRUE),
          nas = ~ sum(is.na(as.numeric(.)))
          #,mad = ~ mad(as.numeric(.),na.rm = TRUE)
          )
        )
      end.time <- Sys.time()
      time.taken <- end.time - start.time
      print(time.taken)
      return (summary_df)
    }
  ))
}

# Generic Methods for quality assurance
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

# let's custamize the summary, for example
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

harmony_create_collection <- function(dir){
  
    # look for r data
    if (file.exists(paste0("../../temp_files/",dir,'.RData'))){
      cat(paste0("\nLoading plate collection from temp_files"))
      plate_collection <- readRDS(paste0("../../temp_files/",dir,'.RData'))
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
          file_list <- unlist(strsplit(evaluations, "/"))
          e <- as.character(sub('.*(?=.{1}$)', '',file_list[7], perl=T))
          plate_obj_name <- paste0('d',sub(" ", "_", gsub('\\-', '',sub("/","__",paste0(file_list[3])))),"_e",e)
          
          cur_path <- file.path(dir,file_list[3])
          
          ###### INDEX FILE#######
          index_file <- file.path(cur_path, "indexfile.txt")
          
          if(file.exists(index_file)){
            
            idx_df <- read_tsv(index_file) %>%
              select_if(function(x) any(!is.na(x))) %>% # removes columns with all NA
              clean_names() %>% # removes spaces and bad characters for underscores
              mutate(channel_name = tolower(gsub("[[:punct:][:blank:]]", "_", channel_name)),
                    dir = plate_obj_name)
            
            # split dataframes based on stain value
            #idx_df_split <- split(idx_df, idx_df$channel_name)
          }else{
            #nothing right now
            stop('Index data not loaded')
          }
          
          # ###### PLATE FILE#######
          plate_file <- file.path(cur_path, paste0("Evaluation", e), "PlateResults.txt")
          
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
          #cells_file <- file.path(cur_path, paste0("Evaluation", e), "Objects_Population - Nuclei Selected.txt")
          cells_file <- list.files(file.path(cur_path, paste0("Evaluation", e)), full.names = TRUE, pattern = "Objects_Population")[1]
          if(file.exists(cells_file)){
            
            # read file header
            cells_header_raw <- readLines(cells_file, n = 9)
            split_header <- strsplit(cells_header_raw, split="\t")
            metadata_df <- as.data.frame(t(do.call(rbind,split_header[c(1:6)]))) %>%
              row_to_names(row_number = 1)
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
          
          new_plate_class <- plate_class(idx_df, plate_df,cell_df,metadata_df,dir)
          # add plate object to list
          len_list <- length(temp_list) + 1
          temp_list[[plate_obj_name]] <- new_plate_class
        }
        plate_collection <- plate_collection_class(temp_list)
        
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        cat(paste0("\nLoad Time:",time.taken))
        cat(paste0("\nSaving plate collection to temp_files"))
        saveRDS(plate_collection, paste0("../../temp_files/",dir,'.RData'))
    }
    return (plate_collection)
}

#(value - medianzeroval) / madzeroval*1.4826
# madadj<- function(dat) {
#   dt <- dat
#   a<- colnames(dt[,18:79])
#   b<- colnames(dt[,c(1:4, 14:16)])
#   borp<- data.table::melt(dt, id.vars = c("Individual", "Concentration", "Run","Plate","X","Y","Bounding Box","Field","in_focus - Object No in Nuclei", "Column", "Row","Position X [µm]", "Position Y [µm]", "Object No" ) ,
#                           measure.vars = list(testid =  a))
#   borp$value <- as.numeric(borp$value)
#   #borp1<- na.omit(borp)     #if needed?
#   borp1[, grp1 := .GRP, by=.(Individual,Column,Run,Plate,variable)]    
#   borp1<- borp[, grp1 := .GRP, by=.(Individual,Column,Run,Plate, variable)]    
#   ark<- borp1[,  .SD[Concentration == 0, .(medval=median(value), madval=mad(value)), by=c("Individual", "Run","Plate","Column", "variable")], by=grp1]
#   carry<- left_join(borp1,ark[,c("grp1", "medval","madval")], by = "grp1")
#   setnames(carry, tolower(names(carry)))
#   newc<- carry[,`:=`(adjval = (value-medval)/madval*1.4826), by="grp1"]
#   newmean<- newc[,`:=`(wellval = mean(adjval)), by=c("individual", "run","plate","column", "variable", "row")]
#   newmedian<- newc[,`:=`(wellval = median(adjval)), by=c("individual", "run","plate","column", "variable", "row")]
#   medtable<- data.table::dcast(newmedian, individual + run + plate + column + row ~ variable, value.var = "wellval")
#   meantable<- data.table::dcast(newmean, individual + run + plate + column + row ~ variable, value.var = "wellval")
#   medtable
#   meantable
# }
