#Rscript "C:/Users/Ankory/Desktop/Networks_install.R"
#Rscript "C:\ivan_work\Networks\data\Networks_install.R"
#Rscript "G:\work\Networks_test\Networks_install.R"
#Rscript "G:\work\my_r_packages\Networks\data\Networks_install.R"
Sys.setenv(LANG="en")
user_home <- file.path(Sys.getenv("HOME"))
user_lib <- file.path(Sys.getenv("HOME"), "R", "libs")
here<-try(file.create(paste0(user_home,"/","test_test")))

mirror_list_3<-c("https://cloud.r-project.org",
"https://cran.uni-muenster.de/")


options(repos = c(CRAN = mirror_list_3[1]))
if (!file.exists(paste0(user_home,"/","test_test"))) {print("error 1")}
# if (!require("BiocManager", quietly = TRUE))
try({install.packages("BiocManager",lib=user_lib)
BiocManager::install(version = "3.22",lib=user_lib,ask=FALSE)
})



persistent_install_packages<-function(pkgs,...,randomize_mirror_order=FALSE,global_tries_max=2,force_reinstall=FALSE,mirror_list_cran=NA,mirror_list_bioc=NA,mode_cran_or_bioc=c("CRAN"),lib_path=file.path(Sys.getenv("HOME"), "R", "libs")) {
	#The idea of this fucntion is to be as hands off and persistent as reasonable in installing the package in the enviornment of not-so-stable internet connection.
	#Well, if that is the idea - we should have an individual access to every package installation
	#Therefore, we are feeding install packages or bioc manager with packages 1 by 1.
	# pkg_i<-"KEGGREST"
	# BiocManager::install(pkg_i,dependencies=FALSE,lib=lib_path,force = TRUE)
	# user_lib<-lib
	# print(user_lib)
	if (is.na(mirror_list_cran)) {
	mirror_list_cran<-c("https://cloud.r-project.org","https://cran.uni-muenster.de/")
	}
	if (is.na(mirror_list_bioc)) {
	mirror_list_bioc<-as.character(1:16)
	}
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
		is_bioc_manager_installed<-suppressMessages(require("BiocManager",character.only=TRUE,lib=lib_path))
		# is_bioc_generics_installed<-suppressMessages(require("BiocGenerics",character.only=TRUE))
		# is_bioc_version_installed<-suppressMessages(require("BiocVersion",character.only=TRUE))

		if (!is_bioc_manager_installed) {
			suppressWarnings(persistent_install_packages("BiocManager", mirror_list="https://cloud.r-project.org",mode_cran_or_bioc="CRAN",force_reinstall=FALSE))
			did_it_install_correctly<-suppressMessages(require("BiocManager",character.only=TRUE,lib=lib_path))
			if (!did_it_install_correctly) {
				stop("Biocmanager is not installed and could not be installed automatically in bioconductor mode")
			}
		
		}
		
		
	}
	if (mode_cran_or_bioc==tolower("Bioconductor")) {
	mirror_list<-mirror_list_bioc
	} else if (mode_cran_or_bioc==tolower("cran")) {
	mirror_list<-mirror_list_cran
	} else {
	print(mode_cran_or_bioc)
	stop("no expected mode detected")
	}
	# print("mode check succesfull")
	mirror_list<-as.character(mirror_list)
	for (i in 1:length(pkgs)) {
		pkg_i<-pkgs[i]
		print(paste0("Attempting to install package ",pkg_i, " in mode ",mode_cran_or_bioc))
		global_tries<-0
		is_pkg_i_installed_already<-suppressMessages(require(pkg_i,character.only=TRUE,lib=lib_path))
		if (force_reinstall & is_pkg_i_installed_already) {
		print(paste0("Removing package ",pkg_i, " since force reinstall is enabled"))
		remove.packages(pkg_i,lib=lib_path)
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
				current_mirror<-as.character(current_mirror)
				tryCatch({
				options(repos = c(CRAN = mirror_list_cran[1]))
				# options(repos = c(CRAN = "@CRAN@"))
				# pkg_i2<-"KEGGREST"
				# BiocManager::install(pkg_i2,dependencies=FALSE,lib=lib_path,force = TRUE,site_repository=repositories()[1])

				av_pack_cran<-available.packages()
				package_dep_inside<-tools::package_dependencies(packages=pkg_i,db=av_pack_cran,recursive=FALSE)
				is_a_cran_pack<-pkg_i %in% av_pack_cran[,1]
				# print("we are here")
				# print(dim(av_pack_cran))
				if (mode_cran_or_bioc==tolower("Bioconductor")) {
				# suppressMessages(options(repos = c(CRAN = repositories()[1])))
				options(repos = c(CRAN = repositories()[1]))
				av_pack_bioc<-available.packages()
				is_a_bioc_pack<-pkg_i %in% av_pack_bioc[,1]
				} else {
				av_pack_bioc<-FALSE
				is_a_bioc_pack<-FALSE
				}	
				# print(dim(av_pack_bioc))
				# print(paste0(pkg_i," is a cran package:",is_a_cran_pack))
				# print(paste0(pkg_i," is a bioc package:",is_a_bioc_pack))
				# chooseBioCmirror(ind=as.character(current_mirror))
				if (is_a_cran_pack) {
				package_dep_inside<-tools::package_dependencies(packages=pkg_i,db=av_pack_cran,recursive=FALSE)
				} else if (is_a_bioc_pack) {
				package_dep_inside<-tools::package_dependencies(packages=pkg_i,db=av_pack_bioc,recursive=FALSE)
				} else {
				stop(paste0("Package ",pkg_i, " found in neither CRAN or Bioc"))
				}

			

				toinstall<-unlist(package_dep_inside,use.names=FALSE)
				# toinstall<-toinstall[!(toinstall %in% c("BiocManager")]

				
					if (!is.null(toinstall) & !any(is.na(toinstall)) & !identical(character(0),toinstall)) {
						# if (length) {
						stat_is_installed_already=rep(FALSE,length(toinstall))
					for (j in 1:length(toinstall)) {
						tryCatch({
							stat_is_installed_already[j]<-suppressMessages(require(toinstall[j],character.only=TRUE,lib=lib_path))
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
							is_depndency_a_cran_pack<-unmet_dependencies_to_install %in% av_pack_cran
							is_depndency_a_bioc_pack<-unmet_dependencies_to_install %in% av_pack_bioc
							# Here we are working on an assumption that if av_pack can be resolved, mirror is working and therefore we should use it
							for (i in 1:length(unmet_dependencies_to_install)) {
								if (is_depndency_a_cran_pack[i]) {
									suppressWarnings(persistent_install_packages(unmet_dependencies_to_install[i], mode_cran_or_bioc=1,...))
								} else if (is_depndency_a_bioc_pack[i]) {
									suppressWarnings(persistent_install_packages(unmet_dependencies_to_install[i],mode_cran_or_bioc=2,...))								
								}								
						}
						}
					}
					# print(paste0("Attempting to instal package ",pkg_i," check following conditions: stat_is_installed_already ",stat_is_installed_already, " unmet_dependencies_to_install ", unmet_dependencies_to_install ))
					if (is.na(pkg_i)) {stop("trying to install package named NA, this should not happen")}
					
				}, error=function(cond){
					message(paste0("Error while trying to install package ", pkg_i, " using mirror ",current_mirror))
					message("Original error message:")
					message(conditionMessage(cond))
				}, warning=function(cond) {
					message(paste0("Install packages caused a warning while installing  ", pkg_i, " using mirror ",current_mirror))
					message("Original warning message:")
					message(conditionMessage(cond))
				}, finally={
						is_pkg_i_installed_already<-suppressMessages(require(pkg_i,character.only=TRUE,lib=lib_path))
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
persistent_install_packages(c("AnnotationDbi"),mode_cran_or_bioc=2,lib_path=user_lib,randomize_mirror_order=FALSE,global_tries_max=3,force_reinstall=TRUE,mirror_list=1:3)

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
test_f<-function () {
		options(repos = c(CRAN = mirror_list_3[1]))
		pkg_i2<-"KEGGREST"
		BiocManager::install(pkg_i2,dependencies=FALSE,lib=user_lib,force = TRUE)
}
test_f()