example_config_path <- system.file("extdata/example/example_main_config.txt", package = "Networks")

example_wd <- setwd(example_dir)
on.exit(setwd(example_wd))

run_pipeline(example_config_path)
