# 1.1 Installing LIANA
{
# install.packages('devtools')
# library(devtools)
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE) # ignore warning from iTALK 
# BiocManager::install("ComplexHeatmap") # required for Connectome
# devtools::install_github('saezlab/OmnipathR@ff3ad88e3915747e1b557bf44ac5396f9525dd7e') # install 4.0 version of OmnipathR
# 
# # install tools
# devtools::install_github("sqjin/CellChat")  
# devtools::install_github('msraredon/Connectome', ref = 'master')   
# devtools::install_github("Coolgenome/iTALK", build_vignettes = TRUE)
# 
# # A modified version of SingleCellSignalR (SCA) that enables external resources
# devtools::install_github(repo = "CostaLab/SingleCellSignalR_v1", subdir = "SingleCellSignalR")
# 
# # Finally, install LIANA
# devtools::install_github('saezlab/liana')
# 
# # Conda env needs to be set up separately
}

# 1.2 Load Packages
{
require(tidyverse)
require(Seurat)
require(liana)
require(tibble)
require(purrr)
  
#library(pillar)

}

# 2.1 Liana Wrapper

  # get liana package
  liana_path <- system.file(package = 'liana')
  # retrieve testdata from liana package
  testdata <- 
    readRDS(file.path(liana_path, "testdata", "input", "testdata.rds"))
  
  testdata %>% glimpse
  
  
  
  
  
  
  
  

