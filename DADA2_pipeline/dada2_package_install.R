# Installation of packages requires for DADA2 pipeline
# Bolívar Aponte Rolón
# 14 June 2023


# Activate commands as needed.
# DADA2 package and associated

#install.packages("gert", repos = c(
#  ropensci = 'https://ropensci.r-universe.dev',
#  CRAN = 'https://cloud.r-project.org')) #usethis needs gert Doesn't install on the HPC cluster for some reason.
#install.packages("usethis") #Doesn't install on the HPC cluster for some reason.
#install.packages("devtools") #Doesn't install on the HPC cluster for some reason.
#devtools::install_github("benjjneb/dada2") #change the ref argument to get other versions

install.packages("Rcpp")

# if (!require("BiocManager", quietly = TRUE)){
#   install.packages("BiocManager", repo="http://cran.rstudio.com/")}
# BiocManager::install(version = "3.17") #Version 3.17
# BiocManager::install(c("dada2", "ShortRead", "Biostrings"))

if (!require("BiocManager", quietly = TRUE)){ #Another way of installing the latest version of dada2
  install.packages("BiocManager", repo="http://cran.rstudio.com/")}

BiocManager::install(version='devel', ask = FALSE) #BiocManager 3.17 (dada2 1.28.0) or developer version (dada2 1.29.0)
BiocManager::install(c("dada2", "ShortRead", "BioStrings"))

packageVersion("dada2") #checking if it is the latest version
packageVersion("ShortRead")
packageVersion("Biostrings")