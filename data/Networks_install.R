#Rscript "C:/Users/Ankory/Desktop/Networks_install.R"
#Rscript "C:\ivan_work\Networks\data\Networks_install.R"
Sys.setenv(LANG="en")
user_home <- file.path(Sys.getenv("HOME"))
user_lib <- file.path(Sys.getenv("HOME"), "R", "libs")
here<-try(file.create(paste0(user_home,"/","test_test")))
if (!file.exists(paste0(user_home,"/","test_test"))) {print("error 1")}


# install.packages("BiocManager",lib=user_lib,repos="https://cloud.r-project.org")
# BiocManager::install(version = "3.22",lib=user_lib)
bioconductor_prefix<-"cannot open URL \\'https://bioconductor.org/"
test<-"cannot open URL 'https://bioconductor.org/packages/3.22/books/bin/windows/contrib/4.5/PACKAGES.rds': HTTP status was '404 Not Found'"
mirror_list_3<-c("https://cloud.r-project.org",
"https://cran.uni-muenster.de/")
tries<-15
# options(repos = c(CRAN = mirror_list_3[1]))
full_cran_mirror_list<-c("https://cloud.r-project.org/",
"http://mirror.fcaglp.unlp.edu.ar/CRAN/",
"https://cran.csiro.au/",
"https://mirror.aarnet.edu.au/pub/CRAN/",
"https://cran.ms.unimelb.edu.au/",
"https://cran.wu.ac.at/",
"https://www.freestatistics.org/cran/",
"https://ftp.belnet.be/mirror/CRAN/",
"https://cran-r.c3sl.ufpr.br/",
"https://vps.fmvz.usp.br/CRAN/",
"https://brieger.esalq.usp.br/CRAN/",
"https://ftp.uni-sofia.bg/CRAN/",
"https://muug.ca/mirror/cran/",
"https://mirror.csclub.uwaterloo.ca/CRAN/",
"https://cran.mirror.rafal.ca/",
"https://cran.dcc.uchile.cl/",
"https://mirrors.tuna.tsinghua.edu.cn/CRAN/",
"https://mirrors.bfsu.edu.cn/CRAN/",
"https://mirrors.pku.edu.cn/CRAN/",
"https://mirrors.ustc.edu.cn/CRAN/",
"https://mirrors.zju.edu.cn/CRAN/",
"https://mirror-hk.koddos.net/CRAN/",
"https://mirrors.qlu.edu.cn/CRAN/",
"https://mirror.lzu.edu.cn/CRAN/",
"https://mirrors.nju.edu.cn/CRAN/",
"https://mirrors.sjtug.sjtu.edu.cn/cran/",
"https://mirrors.sustech.edu.cn/CRAN/",
"https://mirrors.hust.edu.cn/CRAN/",
"https://mirrors.nwafu.edu.cn/cran/",
"https://mirror.uned.ac.cr/cran/",
"https://mirror.library.ucy.ac.cy/cran/",
"https://mirrors.dotsrc.org/cran/",
"https://cran.asia/",
"https://mirror.cedia.org.ec/CRAN/",
"https://cran.030-datenrettung.de/",
"https://pbil.univ-lyon1.fr/CRAN/",
"https://mirror.ibcp.fr/pub/CRAN/",
"https://cran.asnr.fr/",
"https://ftp.fau.de/cran/",
"https://cran.datenrettung360.de/",
"https://ftp.gwdg.de/pub/misc/cran/",
"https://mirror.dogado.de/cran/",
"https://cran.uni-muenster.de/",
"https://mirror.clientvps.com/CRAN/",
"https://mirror.kamp.de/cran/",
"https://ftp.cc.uoc.gr/mirrors/CRAN/",
"https://cran.r-project.hu/",
"https://cran.hafro.is/",
"https://cran.icts.res.in/",
"https://mirror.niser.ac.in/cran/",
"https://cran.isid.ac.in/",
"https://cran.usk.ac.id/",
"https://cran.um.ac.ir/",
"https://cran.mirror.garr.it/CRAN/",
"https://cran.stat.unipd.it/",
"https://ftp.yz.yamagata-u.ac.jp/pub/cran/",
"https://cran.itam.mx/",
"https://est.colpos.mx/",
"https://mirror.marwan.ma/cran/",
"https://mirrors.evoluso.com/CRAN/",
"https://mirror.lyrahosting.com/CRAN/",
"https://cran.stat.auckland.ac.nz/",
"https://cran.uib.no/",
"https://cran.radicaldevelop.com/",
"https://mirror.truenetwork.ru/CRAN/",
"https://mirror.maeen.sa/cran/",
"https://ftp.cixug.es/CRAN/",
"https://cran.rediris.es/",
"https://mirror.accum.se/mirror/CRAN/",
"https://stat.ethz.ch/CRAN/",
"https://mirror.metanet.ch/cran/",
"https://cran.csie.ntu.edu.tw/",
"https://www.stats.bris.ac.uk/R/",
"https://cran.ma.imperial.ac.uk/",
"https://mirror.las.iastate.edu/CRAN/",
"http://ftp.ussg.iu.edu/CRAN/",
"https://mirror.its.umich.edu/cran/",
"https://cran.wustl.edu/",
"https://archive.linux.duke.edu/cran/",
"https://cran.case.edu/",
"https://ftp.osuosl.org/pub/cran/",
"https://lib.stat.cmu.edu/R/CRAN/",
"https://cran.mirrors.hoobly.com/",
"https://mirrors.nics.utk.edu/cran/",
"https://mirror.chpc.utah.edu/pub/cran/",
"https://cran.nyuad.nyu.edu/",
"https://espejito.fder.edu.uy/cran/",
"https://mirrors.cicku.me/cran/")

mirror_selector<-function(iteration,mirror_list) {
mirror_list_pos<-(iteration %% length(iteration))+1
current_mirror<-mirror_list[mirror_list_pos]
return(current_mirror)
}
mirror_list_cran=sample(full_cran_mirror_list,size=length(full_cran_mirror_list))
mirror_list_cran<-c("https://cloud.r-project.org",mirror_list_cran)
for (i in 1:tries) {
if (require("BiocManager", quietly = TRUE,character.only=TRUE,lib=lib_path)) {break}
tryCatch({
if (!require("BiocManager", quietly = TRUE,character.only=TRUE,lib=lib_path)) {
	current_mirror<-mirror_selector(i,mirror_list_cran)
	install.packages("BiocManager",lib=user_lib,repos=current_mirror)
	# options(repos=NULL)
	chooseBioCmirror(ind=as.character(1))
	BiocManager::install(version = "3.22",lib=user_lib,site_repository=repositories()[1])
}
}, error=function(cond){
					message(paste0("Error while trying to install package ", "BiocManager"))
					message("Original error message:")
					message(conditionMessage(cond))
				}, warning=function(cond) {
					message(paste0("Install packages caused a warning while installing  ", "BiocManager"))
					message("Original warning message:")
					message(conditionMessage(cond))
				}, finally={
						is_pkg_i_installed_already<-suppressMessages(require("BiocManager",character.only=TRUE,lib=user_lib))
						successfull_install<-is_pkg_i_installed_already
						if (successfull_install) {outcome<-"successfully"} else {outcome<-"unsuccessfully"}
						message(paste0("Package ", "BiocManager", " was installed ", outcome," at try ",i))
				})
}
is_installed<-suppressMessages(require("BiocManager",character.only=TRUE,lib=user_lib))



persistent_install_packages<-function(pkgs,...,randomize_mirror_order=FALSE,global_tries_max=2,force_reinstall=FALSE,mirror_list_cran=NA,mirror_list_bioc=1,mode_cran_or_bioc=c("CRAN"),lib_path=file.path(Sys.getenv("HOME"), "R", "libs")) {
	#The idea of this fucntion is to be as hands off and persistent as reasonable in installing the package in the enviornment of not-so-stable internet connection.
	#Well, if that is the idea - we should have an individual access to every package installation
	#Therefore, we are feeding install packages or bioc manager with packages 1 by 1.
	acceptable_modes<-tolower(as.character(c("CRAN","Bioconductor",1,2)))
	mode_cran_or_bioc<-tolower(as.character(mode_cran_or_bioc[1]))
	if (is.na(mirror_list_bioc[1])) {mirror_list_bioc<-as.character(1:16)}
	if (is.na(mirror_list_cran[1])) {mirror_list_cran<-c("https://cloud.r-project.org",
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
				if (mode_cran_or_bioc==tolower("CRAN")) {
					current_mirror<-mirror_selector(global_tries,mirror_list_cran)
					
				} else if ((mode_cran_or_bioc==tolower("Bioconductor")) ){
					current_mirror<-mirror_selector(global_tries,mirror_list_bioc)
				} else {
				stop(paste0("Error mode ", mode_cran_or_bioc))
				}
				tryCatch({
				
				if (mode_cran_or_bioc==tolower("CRAN")) {
						message<-paste0("Runnin install packages to install ", pkg_i, " in mode ",mode_cran_or_bioc, " from mirror ", current_mirror)
						install.packages(pkg_i,dependencies=TRUE,lib=lib_path,repos=current_mirror,force = FALSE)
						message<-paste0("Install packages ran sucesfully")
						print(message)

					} else if (mode_cran_or_bioc==tolower("bioconductor")) {
					# print("test")
						message<-paste0("Runnin BIOCmanager INSTALL to install ", pkg_i, " in mode ",mode_cran_or_bioc, " from mirror ", current_mirror)
						print(message)
						chooseBioCmirror(ind=as.character(current_mirror))
						BiocManager::install(pkgs=pkg_i,dependencies=TRUE,lib=lib_path,force = FALSE,site_repository=repositories()[1])
						message<-paste0("BIOCmanager finished sucesfully")
						print(message)
					}
				


			
					
				}, error=function(cond){
					message(paste0("Error while trying to install package ", pkg_i, " using mirror ",current_mirror))
					message("Original error message:")
					message(conditionMessage(cond))
				}, warning=function(cond) {
					if ()
					message(paste0("Install packages caused a warning while installing  ", pkg_i, " using mirror ",current_mirror))
					message("Original warning message:")
					message(conditionMessage(cond))
					invokeRestart("muffleWarning")
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
persistent_install_packages(c("remote","dplyr"),mode_cran_or_bioc=1,lib=user_lib,randomize_mirror_order=TRUE,global_tries_max=tries,force_reinstall=FALSE,mirror_list_cran=full_cran_mirror_list)
persistent_install_packages(c("AnnotationDbi",
  "org.Hs.eg.db",
  "org.Mm.eg.db",
  "org.Mmu.eg.db",
  "org.Dr.eg.db"
),mode_cran_or_bioc=2,lib_path=user_lib,randomize_mirror_order=FALSE,global_tries_max=tries,force_reinstall=FALSE,mirror_list_bioc=1:16)

# BiocManager::install(c(
  # "AnnotationDbi",
  # "org.Hs.eg.db",
  # "org.Mm.eg.db",
  # "org.Mmu.eg.db",
  # "org.Dr.eg.db"
# ), lib = user_lib, ask = FALSE, force = TRUE)

zip_url <- "https://github.com/AnkoryL/Networks/archive/refs/heads/main.zip"
zip_dest <- file.path(user_home, "Networks-main.zip")
download.file(zip_url, zip_dest, mode = "wb")
unzip(zip_dest, exdir = user_home)
pkg_dir <- file.path(user_home, "Networks-main")
remotes::install_local(pkg_dir, lib = user_lib, force = TRUE, dependencies = TRUE)

remotes::install_github("AnkoryL/Networks",
			ref = "main",
			lib = user_lib,
			force = TRUE,
			dependencies = c("Imports", "LinkingTo"),
			build_vignettes = FALSE,
			upgrade = "never")
			
is_installed<-suppressMessages(require("Networks",character.only=TRUE,lib=user_libs))
if (is_installed) {outcome<-"successfully"} else {outcome<-"unsuccessfully"}
paste0("Package ","Networks", " was installed ", outcome)