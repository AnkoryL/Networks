
user_lib <- file.path(Sys.getenv("HOME"), "R", "libs")
.libPaths(c(user_lib, .libPaths()))

remotes::install_github(
  "AnkoryL/Networks",
  ref = "main",
  lib = user_lib,
  force = TRUE,
  dependencies = c("Imports", "LinkingTo"),
  build_vignettes = FALSE,
  upgrade = "never"
)
