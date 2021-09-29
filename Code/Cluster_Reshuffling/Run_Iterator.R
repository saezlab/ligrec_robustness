#------------------------------------------------------------------------------#
# 0. Introduction and Goals ----------------------------------------------------
{
  # Welcome! This is the Run_Iterator.R script.
  
  # This script runs a wrapper function that analyses the robustness of 
  # predictions made by SC transcriptomic CCI inference methods in relation to
  # cluster annotation inaccuracy. If you want to understand the wrapper, go to 
  # CD_Robustness_Iterator.R, where it is defined. In short, it analyses 
  # robustness in the the following steps:
  

  
}



#------------------------------------------------------------------------------#
# 1. Setup ---------------------------------------------------------------------
{
  # In this segment, we set up all the required infrastructure to run the
  # wrapper. The default parameters for the wrapper are set in
  # CD_Robustness_Iterator.R, but custom parameters can be set in the function
  # call.
  
  # Load required Packages
  require(tidyverse)
  require(Seurat)
  require(liana)
  require(lubridate)
  
  
  # Define the functions needed to perform our analysis
  
  source("Code/Cluster_Reshuffling/CR_LIANA_Functions.R")
  source("Code/Cluster_Reshuffling/CR_Shuffler_Functions.R")
  source("Code/Cluster_Reshuffling/CR_Iterator.R")
  
  source("Code/Cluster_Reshuffling/Iterator_Meta_and_Saves.R")
  source("Code/Cluster_Reshuffling/Iterator_Top_Ranks.R")
  
  
  source("Code/Utilities/Iterator_Functions.R")
  source("Code/Utilities/User_Outputs_and_Plots.R")
  
  
  
  
}

#------------------------------------------------------------------------------#
# 2. Get cluster Reshuffling Robustness Results --------------------------------


testdata_type     <- "liana_test"  # seurat_pbmc or liana_test

# Retrieve either seurat_pbmc or liana_test data
testdata <- extract_Testdata(testdata_type = testdata_type)


# We run the wrapper with default settings and twice the standard permutations
robustness_default <- 
  wrap_cluster_Iterator(testdata      = testdata,
                        testdata_type = testdata_type,
                        methods_vector = c("call_connectome", "call_sca"))


