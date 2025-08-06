user_lib <- file.path(Sys.getenv("HOME"), "R", "libs")
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))
