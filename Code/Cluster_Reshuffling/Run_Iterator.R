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
  
  # Define the iterator wrapper function (wrap_cluster_Iterator), which 
  # produces the end results. The wrapper iterates the evaluator function and 
  # collates its results.
  source("Code/Cluster_Reshuffling/CD_Robustness_Iterator.R")
  # Define general functions for data processing in the iterator
  source("Code/Cluster_Reshuffling/Iterator_Processing_Functions.R")
  # Define functions for plotting the iterator results
  source("Code/Cluster_Reshuffling/Iterator_Plotting.R")
  # Define functions for capturing metadata and saving iterator results
  source("Code/Cluster_Reshuffling/Iterator_Metadata_and_Saves.R")
  

  
  
  
  
}

#------------------------------------------------------------------------------#
# 2. Get cluster Reshuffling Robustness Results ----------------------------------





# We run the wrapper with default settings and twice the standard permutations
robustness_default <- wrap_cluster_Iterator(number_seeds = 2,
                                             methods_vector = c("call_sca"),
                                             testdata_type = "liana_test",
                                             number_ranks = list("call_sca" = 20))


