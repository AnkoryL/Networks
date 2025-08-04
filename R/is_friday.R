is_friday <- function(date) {

  if (!inherits(date, "Date")) {
    date <- lubridate::ymd(date)
  }

  lubridate::wday(date) == 6

}
