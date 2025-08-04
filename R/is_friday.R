#' @title Is Friday?
#' @description Check which elements in a vector of dates are Fridays.
#' @param date A Date class object or a string in the form "yyyymmdd".
#' @returns A boolean vector the same length as `date`.
#' @export
#' @examples
#' # TRUE
#' is_friday("2025-03-28")
#' "2025-03-28" |> lubridate::ymd() |> is_friday()
#'
#' # FALSE
#' is_friday("2025-03-24")
#'
#' # TRUE, FALSE
#' c("2025-03-28", "2025-03-24") |> is_friday()
#'
is_friday <- function(date) {

  if (!inherits(date, "Date")) {
    date <- lubridate::ymd(date)
  }

  lubridate::wday(date) == 6

}
