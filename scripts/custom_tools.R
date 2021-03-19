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

# see http://rstudio-pubs-static.s3.amazonaws.com/1563_1ae2544c0e324b9bb7f6e63cf8f9e098.html
# log_skew_transform <- function(response_var = 'response') {
#   skew_val <- e1071::skewness(response_var)
#   skew.score <- function(c, x)
#     (e1071::skewness(log(x + c))) ^ 2
#   if (!is.na(skew_val)) {
#     print(paste0("performing log transformation with skew:", skew_val))
#     best.c <- optimise(skew.score, c(0, 20), x = response_var)$minimum
#     response_var <- log(response_var + best.c)
#   } else{
#     print("No optimization possible")
#     # to be continued
#     #df$response <- log(df$response)
#   }
#   return(response_var)
# }

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