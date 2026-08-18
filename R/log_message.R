#' Write message to Networks log
#'
#' @param message Text message to save
#' @noRd

log_message <- function(message) {

  log_dir <- file.path(getwd(), "logs")

  if (!dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE)
  }

  log_file <- file.path(
    log_dir,
    paste0("Networks_", Sys.Date(), ".log")
  )

  cat(
    paste0(
      "[",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      "] ",
      message,
      "\n"
    ),
    file = log_file,
    append = TRUE
  )

}
