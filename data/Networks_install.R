#Rscript "C:/Users/Ankory/Desktop/Networks_install.R"
Sys.setenv(LANG="en")
user_home <- file.path(Sys.getenv("HOME"))
user_lib <- file.path(Sys.getenv("HOME"), "R", "libs")
here<-try(file.create(paste0(user_home,"/","test_test")))
if (!file.exists(paste0(user_home,"/","test_test"))) {print("error 1")}


mirror_list_3<-c("https://cloud.r-project.org",
"https://mirrors.tuna.tsinghua.edu.cn",
"https://cran.uni-muenster.de/")

persistent_install_packages<-function(pkgs,...,randomize_mirror_order=FALSE,global_tries_max=10,force_reinstall=FALSE,mirror_list="https://cloud.r-project.org") {
	#The idea of this fucntion is to be as hands off and persistent as reasonable in installing the package in the enviornment of not-so-stable internet connection.
	#Well, if that is the idea - we should have an individual access to every package installation
	#Therefore, we are feeding install packages with packages 1 by 1.
	for (i in 1:length(pkgs)) {
		pkg_i<-pkgs[i]
		global_tries<-0
		is_pkg_i_installed_already<-suppressMessages(require(pkg_i,character.only=TRUE))
		if (force_reinstall & is_pkg_i_installed_already) {
		remove.packages(pkg_i,character.only=TRUE)
		successfull_install<-FALSE
		} else {
		successfull_install<-is_pkg_i_installed_already
		}
		if (randomize_mirror_order) {
			mirror_list=sample(mirror_list,size=length(mirror_list))
		}
			while ((global_tries<global_tries_max) & !(successfull_install)) {
				global_tries<-global_tries+1
				mirror_list_pos<-(global_tries %% length(mirror_list))+1 
				current_mirror<-mirror_list[mirror_list_pos]
				options(repos = c(CRAN = current_mirror))
				tryCatch({
				av_pack<-available.packages()
				# install.packages(pkg_i,...)	
				package_dep_inside<-tools::package_dependencies(packages=pkg_i,db=av_pack,recursive=FALSE)
				# print(package_dep_inside)
				toinstall<-unlist(package_dep_inside,use.names=FALSE)
					if (!is.null(toinstall) & !any(is.na(toinstall)) & !identical(character(0),toinstall)) {
						# print(toinstall)
						# if (length) {
						stat_is_installed_already=rep(FALSE,length(toinstall))
					for (j in 1:length(toinstall)) {
						tryCatch({
							stat_is_installed_already[j]<-suppressMessages(require(toinstall[j],character.only=TRUE))
						}, error=function(cond) {
							message(paste0("Error while checking for dependency ",toinstall[j], " of package ", pkg_i, " using mirror ",current_mirror))
							message("Original error message:")
							message(conditionMessage(cond))
						}, warning=function(cond) {
							# message(paste0("Warning checking for dependency ",toinstall[j], " of package ", pkg_i, " using mirror ",current_mirror))
							# message("Original warning message:")
							# message(conditionMessage(cond))
						})
					}


						if (sum(!stat_is_installed_already)!=0){
							unmet_dependencies_to_install<-toinstall[!stat_is_installed_already]
							# print(unmet_dependencies_to_install)
							print(paste("Missing dependencies",paste(unmet_dependencies_to_install,sep=" ",collapse=" "), sep=" ",collapse=NULL))
							# Here we are working on an assumption that if av_pack can be resolved, mirror is working and therefore we should use it
							suppressWarnings(persistent_install_packages(unmet_dependencies_to_install, mirror_list=current_mirror))
						}	
					}
					# print(paste0("Attempting to instal package ",pkg_i," check following conditions: stat_is_installed_already ",stat_is_installed_already, " unmet_dependencies_to_install ", unmet_dependencies_to_install ))
					if (is.na(pkg_i)) {stop("trying to install package named NA, this should not happen")}
					install.packages(pkg_i,...,dependencies=NULL)
				}, error=function(cond){
					message(paste0("Error while trying to install package ", pkg_i, " using mirror ",current_mirror))
					message("Original error message:")
					message(conditionMessage(cond))
				}, warning=function(cond) {
					message(paste0("Install packages caused a warning while installing  ", pkg_i, " using mirror ",current_mirror))
					message("Original error message:")
					message(conditionMessage(cond))
				}, finally={
						is_pkg_i_installed_already<-suppressMessages(require(pkg_i,character.only=TRUE))
						successfull_install<-is_pkg_i_installed_already
						if (successfull_install) {outcome<-"successfully"} else {outcome<-"unsuccessfully"}
						message(paste0("Package ",pkg_i, " was installed ", outcome," at try ",global_tries, " using mirror ",current_mirror))
				})
			
			}
	}
}


rprofile_path <- file.path(user_home, ".Rprofile")
if (!file.exists(rprofile_path)) {
  cat(
    ".First <- function() {\n",
    "  user_lib <- file.path(Sys.getenv('HOME'), 'R', 'libs')\n",
    "  if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)\n",
    "  .libPaths(c(user_lib, .libPaths()))\n",
    "}\n",
    file = rprofile_path
  )
}

if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))


# install.packages("BiocManager",lib=user_lib)
# persistent_install_packages(c("BiocManager","remote","dplyr"),lib=user_lib,randomize_mirror_order=TRUE,global_tries_max=10,force_reinstall=TRUE,mirror_list=mirror_list_3)
persistent_install_packages(c("dplyr"),lib=user_lib,randomize_mirror_order=TRUE,global_tries_max=10,force_reinstall=TRUE,mirror_list=mirror_list_3)

# BiocManager::install(c(
  # "AnnotationDbi",
  # "org.Hs.eg.db",
  # "org.Mm.eg.db",
  # "org.Mmu.eg.db",
  # "org.Dr.eg.db"
# ), lib = user_lib, ask = FALSE, force = TRUE)

#zip_url <- "https://github.com/AnkoryL/Networks/archive/refs/heads/main.zip"
#zip_dest <- file.path(user_home, "Networks-main.zip")
#download.file(zip_url, zip_dest, mode = "wb")
#unzip(zip_dest, exdir = user_home)
#pkg_dir <- file.path(user_home, "Networks-main")
#remotes::install_local(pkg_dir, lib = user_lib, force = TRUE, dependencies = TRUE)

# remotes::install_github("AnkoryL/Networks",
			# ref = "main",
			# lib = user_lib,
			# force = TRUE,
			# dependencies = c("Imports", "LinkingTo"),
			# build_vignettes = FALSE,
			# upgrade = "never")

