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