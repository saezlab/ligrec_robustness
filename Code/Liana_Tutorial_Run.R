# 1.1 Installing LIANA
{
install.packages('tidyverse')
install.packages('Seurat')
install.packages('RobustRankAggreg')
install.packages('devtools')
library(devtools)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE) # ignore warning from iTALK
BiocManager::install("ComplexHeatmap") # required for Connectome
devtools::install_github('saezlab/OmnipathR@ff3ad88e3915747e1b557bf44ac5396f9525dd7e') # install 4.0 version of OmnipathR

# install tools
devtools::install_github("sqjin/CellChat")
devtools::install_github('msraredon/Connectome', ref = 'master')
devtools::install_github("Coolgenome/iTALK", build_vignettes = TRUE)

# A modified version of SingleCellSignalR (SCA) that enables external resources
devtools::install_github(repo = "CostaLab/SingleCellSignalR_v1", subdir = "SingleCellSignalR")

# Finally, install LIANA
devtools::install_github('saezlab/liana')

# Conda env needs to be set up separately
}

# 1.2 Updating Liana
{
# #### UPDATING LIANA ON WINDOWS: remove.packages('liana'), delete the liana folder from "C:\Users\plabu\OneDrive\Documents\R\win-library\4.1", restart computer, then:
# library(devtools)
# Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE) # ignore warning from iTALK
# devtools::install_github('saezlab/liana')
# ### and ignore any calls to update more packages on the way
}

# 2. Load Packages
{
require(tidyverse)
require(Seurat)
require(liana)
  
load(RobustRankAggreg)
}

# 3. Liana Wrapper
{
  # Load Testdata
  {
  liana_path <- system.file(package = 'liana')                                    # get liana package filepath
  testdata <- 
    readRDS(file.path(liana_path, "testdata", "input", "testdata.rds"))           # retrieve testdata from filepath

  testdata %>% glimpse                                                            # View Testdata
  }
  
  # Run wrapper on testdata for omnipath x squidppy and connectome
  liana_test <- liana_wrap(testdata,
                           method = c('cellchat', 'connectome', 'italk', 'natmi', 'sca'),
                           resource = c('OmniPath'))                              # All methods but Squidpy Work so far
  
}
  
# 4. Consensus Ranking
{
liana_aggregated_ranks <- liana_test %>% liana_aggregate()
liana_aggregated_ranks
}
