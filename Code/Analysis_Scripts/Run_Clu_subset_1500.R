#------------------------------------------------------------------------------#
# 1. Setup ---------------------------------------------------------------------
{
  # In this segment, we set up all the required infrastructure to run the
  # wrapper. The default parameters for the wrapper are set in
  # CR_Iterator.R, but custom parameters can be set in the function call.
  
  # Load required Packages
  require(tidyverse)
  require(Seurat)
  require(liana)
  require(lubridate)
  
  
  # Define the functions needed to perform our analysis
  
  # Define wrap_cluster_Iterator(), a function that performs the core analysis
  source("Code/Cluster_Reshuffling/CR_Iterator.R")
  # Define Iterator helper functions that work with the top-ranked CCI's
  source("Code/Cluster_Reshuffling/Iterator_Top_Ranks.R")
  # Define Iterator helpers that work at a meta level, saving results and so on 
  source("Code/Cluster_Reshuffling/Iterator_Meta_and_Saves.R")

  
  # Define functions for generating the many reshuffled cluster annotations
  source("Code/Cluster_Reshuffling/CR_Shuffler_Functions.R")
  # Define functions for subsetting cluster annotations, instead of reshuffling
  source("Code/Cluster_Reshuffling/CR_Subsetter_Functions.R")
  # Define functions for iterating LIANA on all the reshuffled clusters
  source("Code/Cluster_Reshuffling/CR_LIANA_Functions.R")


  # Define utility functions that help the iterator run 
  source("Code/Utilities/Iterator_Functions.R")
  # Define utiltiy functions for user-facing console outputs and plotting
  source("Code/Utilities/User_Outputs_and_Plots.R")
  # Define common functions that help with NATMI parameters
  source("Code/Utilities/LIANA_Utilities.R")
  
  
}



#------------------------------------------------------------------------------#
# 2. Get Cluster Reshuffling Robustness Results --------------------------------
{
  
  
  # First we load testdata from the data folder. 
  # We also give a label (testdata_type, choose "seurat_pbmc" or "liana_test")
  testdata_type <- "liana_test"  
  testdata      <- extract_Testdata(testdata_type = testdata_type)
  
  
  # We run the wrapper function, feeding it the testdata and the testdata label
  robustness_reshuffle_default <- 
    wrap_cluster_Iterator(testdata      = testdata,
                          testdata_type = testdata_type,
                          reshuffle_or_subset = "subset",
                          top_n = 1500)
  
  
}



