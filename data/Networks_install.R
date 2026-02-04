#Rscript "C:/Users/Ankory/Desktop/Networks_install.R"
Sys.setenv(LANG="en")
user_home <- file.path(Sys.getenv("HOME"))
user_lib <- file.path(Sys.getenv("HOME"), "R", "libs")
here<-try(file.create(paste0(user_home,"/","test_test")))
if (!file.exists(paste0(user_home,"/","test_test"))) {print("error 1")}


mirror_list_3<-c("https://cloud.r-project.org",
"https://cran.uni-muenster.de/")

persistent_install_packages<-function(pkgs,...,randomize_mirror_order=FALSE,global_tries_max=2,force_reinstall=FALSE,mirror_list="https://cloud.r-project.org",mode_cran_or_bioc=c("CRAN")) {
	#The idea of this fucntion is to be as hands off and persistent as reasonable in installing the package in the enviornment of not-so-stable internet connection.
	#Well, if that is the idea - we should have an individual access to every package installation
	#Therefore, we are feeding install packages or bioc manager with packages 1 by 1.
	
	# this part is enitrely about making sure that correct mode_cran_or_bioc is chosen, correct mirror list is chosen and bioconductor is installed (if needed)
	acceptable_modes<-tolower(as.character(c("CRAN","Bioconductor",1,2)))
	mode_cran_or_bioc<-tolower(as.character(mode_cran_or_bioc[1]))
	# print(mode_cran_or_bioc)
	# print((mode_cran_or_bioc == tolower("Bioconductor")))
	if (!(mode_cran_or_bioc %in% acceptable_modes)) {
	stop(paste("Mode ",mode_cran_or_bioc, "in not an acceptable mode - select one of the",sep=" ",collapse=" "))
	} else if (mode_cran_or_bioc %in% acceptable_modes[3:4]) { 
	mode_cran_or_bioc<-acceptable_modes[as.numeric(mode_cran_or_bioc)]
	}
	if (mode_cran_or_bioc == tolower("Bioconductor")) {
		is_bioc_manager_installed<-suppressMessages(require("BiocManager",character.only=TRUE))

		if (!is_bioc_manager_installed) {
			suppressWarnings(persistent_install_packages("BiocManager", mirror_list="https://cloud.r-project.org",mode_cran_or_bioc="CRAN",force_reinstall=FALSE))
			did_it_install_correctly<-suppressMessages(require("BiocManager",character.only=TRUE))
			if (!did_it_install_correctly) {
				stop("Biocmanager is not installed and could not be installed automatically in bioconductor mode")
			}
		
		}
		
	}
	if (mode_cran_or_bioc==tolower("Bioconductor") & any(!is.numeric(mirror_list))) {
	mirror_list<-as.character(1:16)
	} 
	# print("mode check succesfull")
	
	for (i in 1:length(pkgs)) {
		pkg_i<-pkgs[i]
		global_tries<-0
		is_pkg_i_installed_already<-suppressMessages(require(pkg_i,character.only=TRUE))
		if (force_reinstall & is_pkg_i_installed_already) {
		print(paste0("Removing package ",pkg_i, " since force reinstall is enabled"))
		remove.packages(pkg_i)
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
				# print(current_mirror)
				tryCatch({
				print("check 1")
				if (mode_cran_or_bioc==tolower("CRAN")) {
					options(repos = c(CRAN = current_mirror))
					av_pack<-available.packages()

				} 
				if (mode_cran_or_bioc==tolower("Bioconductor")) {
					print("we are here")
					chooseBioCmirror(ind=as.character(current_mirror))
					options(repos = c(CRAN = repositories()[1]))

					av_pack<-available.packages()
					print("av_pack fetch sucesfull")
				}

				print("check 2")
				# print(dim(av_pack))
				print("check 3")
				package_dep_inside<-tools::package_dependencies(packages=pkg_i,db=av_pack,recursive=FALSE)				

				# print(package_dep_inside)
				toinstall<-unlist(package_dep_inside,use.names=FALSE)
				toinstall<-toinstall[!(toinstall=="BiocManager")]
				
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
								if (mode_cran_or_bioc==tolower("CRAN")) {
									suppressWarnings(persistent_install_packages(unmet_dependencies_to_install, mirror_list=current_mirror,mode_cran_or_bioc=1))
								} else if (mode_cran_or_bioc==tolower("Bioconductor")) {
									suppressWarnings(persistent_install_packages(unmet_dependencies_to_install,mode_cran_or_bioc="bioconductor", mirror_list=current_mirror))
								}
							
						}	
					}
					# print(paste0("Attempting to instal package ",pkg_i," check following conditions: stat_is_installed_already ",stat_is_installed_already, " unmet_dependencies_to_install ", unmet_dependencies_to_install ))
					if (is.na(pkg_i)) {stop("trying to install package named NA, this should not happen")}
					if (mode_cran_or_bioc==tolower("CRAN")) {
						install.packages(pkg_i,...,dependencies=NULL)
					} else if (mode_cran_or_bioc==tolower("bioconductor")) {
						BiocManager::install(pkg_i,...,dependencies=NULL)
					}
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
# persistent_install_packages(c("dplyr"),mode_cran_or_bioc=1,lib=user_lib,randomize_mirror_order=TRUE,global_tries_max=10,force_reinstall=TRUE,mirror_list=mirror_list_3)
persistent_install_packages(c("AnnotationDbi"),mode_cran_or_bioc=2,lib=user_lib,randomize_mirror_order=FALSE,global_tries_max=10,force_reinstall=TRUE,mirror_list=1:3)

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

