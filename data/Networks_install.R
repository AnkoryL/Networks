#Rscript "C:/Users/Ankory/Desktop/Networks_install.R"
#Rscript "C:\ivan_work\Networks\data\Networks_install.R"
#Rscript "G:\work\Networks_test\Networks_install.R"
#Rscript "G:\work\my_r_packages\Networks\data\Networks_install.R"
Sys.setenv(LANG="en")
user_home <- file.path(Sys.getenv("HOME"))
user_lib <- file.path(Sys.getenv("HOME"), "R", "libs")
here<-try(file.create(paste0(user_home,"/","test_test")))
options(repos = c(CRAN = mirror_list_3[1]))
if (!file.exists(paste0(user_home,"/","test_test"))) {print("error 1")}

if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager",lib=user_lib)
BiocManager::install(version = "3.22",lib=user_lib)

}

mirror_list_3<-c("https://cloud.r-project.org",
"https://cran.uni-muenster.de/")






persistent_install_packages<-function(pkgs,...,randomize_mirror_order=FALSE,global_tries_max=2,force_reinstall=FALSE,mirror_list_cran=NA,mirror_list_bioc=NA,mode_cran_or_bioc=c("CRAN"),lib_path=file.path(Sys.getenv("HOME"), "R", "libs")) {
	#The idea of this fucntion is to be as hands off and persistent as reasonable in installing the package in the enviornment of not-so-stable internet connection.
	#Well, if that is the idea - we should have an individual access to every package installation
	#Therefore, we are feeding install packages or bioc manager with packages 1 by 1.
	acceptable_modes<-tolower(as.character(c("CRAN","Bioconductor",1,2)))
	mode_cran_or_bioc<-tolower(as.character(mode_cran_or_bioc[1]))
	if (is.na(mirror_list_bioc)) {mirror_list_bioc<-as.character(1:16)}
	if (is.na(mirror_list_cran)) {mirror_list_cran<-c("https://cloud.r-project.org",
	"https://cran.uni-muenster.de/")}

	if (!(mode_cran_or_bioc %in% acceptable_modes)) {
	stop(paste("Mode ",mode_cran_or_bioc, "in not an acceptable mode - select one of the",sep=" ",collapse=" "))
	} else if (mode_cran_or_bioc %in% acceptable_modes[3:4]) { 
	mode_cran_or_bioc<-acceptable_modes[as.numeric(mode_cran_or_bioc)]
	}
	
	for (i in 1:length(pkgs)) {
		pkg_i<-pkgs[i]
		print(paste0("Attempting to install package ",pkg_i, " in mode ",mode_cran_or_bioc))
		global_tries<-0
		is_pkg_i_installed_already<-require(pkg_i,character.only=TRUE,lib=lib_path, quietly = TRUE)
		if (force_reinstall & is_pkg_i_installed_already) {
		print(paste0("Removing package ",pkg_i, " since force reinstall is enabled"))
		remove.packages(pkg_i,lib=lib_path)
		successfull_install<-FALSE
		} else {
		successfull_install<-is_pkg_i_installed_already
		}
		if (randomize_mirror_order) {
			mirror_list_cran=sample(mirror_list_cran,size=length(mirror_list_cran))
			mirror_list_bioc=sample(mirror_list_bioc,size=length(mirror_list_bioc))
			mirror_list_cran<-c("https://cloud.r-project.org",mirror_list_cran)
			mirror_list_bioc<-c("1",mirror_list_bioc)
		}
			while ((global_tries<global_tries_max) & !(successfull_install)) {
				global_tries<-global_tries+1
				mirror_list_pos<-(global_tries %% length(mirror_list))+1 
				tryCatch({
				
				if (mode_cran_or_bioc==tolower("CRAN")) {
						message<-paste0("Runnin install packages to install ", pkg_i, " in mode ",mode_cran_or_bioc, " from mirror ", current_mirror)
						current_mirror<-mirror_list_cran[mirror_list_pos]
						install.packages(pkg_i,dependencies=TRUE,lib=lib_path,repos=current_mirror,force = FALSE)
						message<-paste0("Install packages ran sucesfully")
						print(message)

					} else if (mode_cran_or_bioc==tolower("bioconductor")) {
					# print("test")
						message<-paste0("Runnin BIOCMANAGER INSTALL to install ", pkg_i, " in mode ",mode_cran_or_bioc, " from mirror ", current_mirror)
						print(message)
						current_mirror<-mirror_list_bioc[mirror_list_pos]
						chooseBioCmirror(ind=as.character(current_mirror))
						BiocManager::install(pkgs=pkg_i,dependencies=TRUE,lib=lib_path,force = FALSE)
						message<-paste0("BIOCmanager finished sucesfully")
						print(message)
					}
				


			
					}
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
persistent_install_packages(c("dplyr"),mode_cran_or_bioc=1,lib=user_lib,randomize_mirror_order=TRUE,global_tries_max=10,force_reinstall=TRUE,mirror_list=mirror_list_3)
persistent_install_packages(c("AnnotationDbi"),mode_cran_or_bioc=2,lib_path=user_lib,randomize_mirror_order=FALSE,global_tries_max=10,force_reinstall=TRUE,mirror_list=1:3)

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