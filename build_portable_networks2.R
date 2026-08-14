create_structure <- function(main_dir = "NetworksPortable") {

  app_dir <- file.path(main_dir, "app")

  dir.create(main_dir, showWarnings = FALSE)
  dir.create(app_dir, showWarnings = FALSE)
  dir.create(file.path(app_dir, "library"), recursive = TRUE, showWarnings = FALSE)

  invisible(
    list(
      main_dir = normalizePath(main_dir),
      app_dir = normalizePath(app_dir),
      library_dir = normalizePath(file.path(app_dir, "library"))
    )
  )
}


download_rportable <- function(version, app_dir) {

  folder_name <- sprintf("portable-r-%s-win-x64", version)
  r_dir <- file.path(app_dir, folder_name)

  if (dir.exists(r_dir)) {
    message("Portable R ", version, " already exists. Skipping download.")
    return(invisible(r_dir))
  }

  filename <- sprintf("%s.zip", folder_name)

  url <- sprintf(
    "https://github.com/portable-r/portable-r-windows/releases/download/v%s/%s",
    version,
    filename
  )

  archive <- file.path(app_dir, filename)

  message("Downloading portable R ", version)

  download.file(
    url,
    archive,
    mode = "wb"
  )

  message("Extracting portable R")

  unzip(
    archive,
    exdir = app_dir
  )

  unlink(archive)

  invisible(r_dir)
}

install_packages <- function(paths) {

  options(
    timeout = 1200,
    download.file.method = "libcurl"
  )

  portable_lib <- paths$library_dir

  cran_packages <- c(
    "shinycssloaders",
    "shiny",
    "DT",
    "dplyr",
    "igraph",
    "openxlsx",
    "magrittr"
  )

  bioc_packages <- c(
    "AnnotationDbi",
    "org.Hs.eg.db",
    "org.Mm.eg.db",
    "org.Mmu.eg.db",
    "org.Dr.eg.db"
  )


  remove_base_packages_from_deps <- function(x) {
    x[!(x %in% rownames(installed.packages(priority = "base")))]
  }


  is_package_cran_or_bioconductor <- function(package) {

    package <- remove_base_packages_from_deps(package)

    out_vec <- rep("unknown", length(package))

    av_pack_cran <- available.packages(repos = findCRANmirror("web"))
    av_pack_bioconductor <- available.packages(repos = BiocManager::repositories())

    out_vec[package %in% av_pack_bioconductor[, 1]] <- "bioc"
    out_vec[package %in% av_pack_cran[, 1]] <- "cran"
    out_vec[package == "terminal_empty_dependency"] <- "cran"

    out_vec
  }


  get_package_dep <- function(package) {

    package <- remove_base_packages_from_deps(package)

    deps <- vector(
      length(package),
      mode = "list"
    )

    names(deps) <- package

    bioc_or_cran <- is_package_cran_or_bioconductor(package)

    db_cran <- available.packages(repos = findCRANmirror("web"))

    if (any(bioc_or_cran == "bioc")) {
      db_bioc <- available.packages(
        repos = BiocManager::repositories()
      )
    }

    for (i in seq_along(package)) {

      if (bioc_or_cran[i] == "cran") {

        deps[[i]] <- unlist(
          tools::package_dependencies(
            packages = package[i],
            db = db_cran,
            recursive = FALSE
          ),
          use.names = FALSE
        )

      } else if (bioc_or_cran[i] == "bioc") {

        deps[[i]] <- unlist(
          tools::package_dependencies(
            packages = package[i],
            db = db_bioc,
            recursive = FALSE
          ),
          use.names = FALSE
        )

      } else {

        stop("Package ", package[i], " was not identified as CRAN or Bioconductor.")
      }
    }

    deps <- lapply(deps,remove_base_packages_from_deps)
    deps
  }


  get_package_deps <- function(package) {

    found <- character(0)
    result <- character(0)

    current <- package

    while (length(current) > 0) {

      current <- current[!(current %in% found)]

      if (length(current) == 0) {
        break
      }

      message("Resolving dependencies for: ",paste(current, collapse = ", "))

      deps <- get_package_dep(current)

      found <- unique(c(found, current))

      dependencies <- unique( unlist(deps, use.names = FALSE))
      dependencies <- dependencies[dependencies != "terminal_empty_dependency"]

      result <- unique(c(result, current, dependencies))
      current <- dependencies[!(dependencies %in% found)]

    }

    result
  }


  message("Resolving CRAN dependencies...")

  cran_full <- get_package_deps(cran_packages)

  message("Resolved ",length(cran_full), " CRAN packages." )

  message("Resolving Bioconductor dependencies...")

  bioc_full <- get_package_deps(bioc_packages)

  message("Resolved ",length(bioc_full)," Bioconductor packages.")


  install_until_done <- function(
    packages,
    installer,
    attempts = 10
  ) {

    for (i in seq_len(attempts)) {

      installed <- rownames(
        installed.packages(
          lib.loc = portable_lib
        )
      )

      missing <- packages[!(packages %in% installed)]

      if (length(missing) == 0) {
        return(invisible(TRUE))
      }

      message("Attempt ", i, ". Remaining packages: ", paste(missing, collapse = ", "))

      try(
        installer(missing),
        silent = TRUE
      )

      Sys.sleep(10)
    }

    installed <- rownames(
      installed.packages(
        lib.loc = portable_lib
      )
    )

    missing <- packages[!(packages %in% installed)]

    if (length(missing) > 0) {
      stop("Failed to install: ",paste(missing, collapse = ", "))
    }

    invisible(TRUE)
  }


  install_until_done(
    cran_full,
    function(x) {

      install.packages(
        x,
        lib = portable_lib,
        dependencies = FALSE
      )

    }
  )


  install_until_done(
    bioc_full,
    function(x) {

      BiocManager::install(
        x,
        lib = portable_lib,
        ask = FALSE,
        update = FALSE,
        force = TRUE,
        dependencies = FALSE
      )

    }
  )


  message("Installing Networks from GitHub...")

  Sys.setenv(GITHUB_PAT = "")

  remotes::install_github(
    "AnkoryL/Networks",
    ref = "dev",
    lib = portable_lib,
    dependencies = FALSE,
    upgrade = "never"
  )

  message("Networks installation completed.")

  invisible(TRUE)
}

write_launcher <- function(paths) {

  r_version <- as.character(getRversion())

  launcher <- paste0(
    "@echo off\n\n",
    "title Networks - Gene Network Analysis\n",
    "echo.\n",
    "echo ==========================================\n",
    "echo          NETWORKS\n",
    "echo ==========================================\n",
    "echo.\n",
    "echo Starting application.\n",
    "echo.\n",
    "echo Loading R and libraries. Please wait.\n",
    "echo.\n\n",
    'set "ROOT=%~dp0"\n',
    'set "APP=%ROOT%app"\n',
    'set "R_HOME=%APP%\\portable-r-', r_version, '-win-x64"\n',
    'set "R_LIBS_USER=%APP%\\library"\n\n',
    'cd /d "%APP%"\n\n',
    '"%R_HOME%\\bin\\Rscript.exe" -e ".libPaths(c(Sys.getenv(\'R_LIBS_USER\'), .libPaths())); suppressPackageStartupMessages(Networks::run_app())"\n'
  )

  writeLines(
    launcher,
    file.path(paths$main_dir, "launcher.bat")
  )
}

main <- function() {

  paths <- create_structure()


  download_rportable(
    version = as.character(getRversion()),
    app_dir = paths$app_dir
  )


  install_packages(paths)


  write_launcher(paths)

}


main()
