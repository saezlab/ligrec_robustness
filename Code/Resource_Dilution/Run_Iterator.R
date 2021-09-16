#------------------------------------------------------------------------------#
# 0. Introduction and Goals ----------------------------------------------------
{
  # Welcome! This is the Run_Iterator.R script.
  
  # The goal of this script is to measure the robustness of CCI inference
  # methods in LIANA when it comes to resource change. This is accomplished in
  # the following steps:
  
  # 1. Run every Cell-Cell-Interaction (CCI) prediction method in LIANA++ while
  #    using OmniPath (OP) and a user specified test data.
  
  # 2. Determine the highest ranked CCI's for each method
  #    -> these are the top ranks in the generic scenario ("0% dilution")
  
  # 3. Modify ("dilute") OP by replacing a proportion of its interactions with
  #    false interactions. These following rules apply for the dilution:
  
  #   -> Interactions that were top-ranked by any method cannot be replaced
  #   -> New interactions are created from genes found in the test data, and
  #      only replace interactions relevant to the test data
  #   -> All new interactions are unique, the diluted resource has no duplicates
  #   -> The diluted interactions don't include any genes that are communication
  #      partner to themselves.
  #   -> Depending on the the selected parameters, it can also be arranged that
  #      the topology of the resource is semi-preserved.
  
  #    In general, this portion has lots of adjustable arguments. What quality
  #    of genes should dilution occur with? Should the output be logged? Etc.
  
  # 4. Once we have OP diluted to various percentages, for example 10 %, 40 %
  #    and 90 %, we rerun all the LIANA++ methods for each proportion of
  #    dilution.
  
  # 5. We determine the top-ranked interactions for these dilution stages, and
  #    then compare the overlap between the new predictions and the old 0 %
  #    prediction. We also determine the rate of diluted interactions in the
  #    top_ranks, this is usually equal to 1 - the overlap.
  
  # 6. Since dilution has elements of randomness, we repeat the above process
  #    a few times so we can collate information from multiple samples. 1-5 are
  #    all handled by resource_Robustness(), a function we iterate here and that
  #    can be found in Resource_Robustness_Functions.R.
  
  # 7. We plot the overlap over the dilution proportion for each method,
  #    combining all the iterations into one box plot. Then we save the plots
  #    and all the results to automatically generated, descriptive file names.
  
  # This entire process is handled by the wrap_resource_Iterator() function, 
  # found in the Robustness_Iterator.R script. In this script, we simply run
  # that wrapper with varying parameters fed to wrap_resource_Iterator(). If you
  # want to understand more on how the function works, please check that script.
}



#------------------------------------------------------------------------------#
# 1. Setup ---------------------------------------------------------------------
{
  # In this segment, we set up all the required infrastructure to run this
  # script. The parameters for how this script runs are not set here, they are
  # set in Iterator_Parameters.R
  
  # Load required Packages
  require(tidyverse)
  require(Seurat)
  require(liana)
  require(lubridate)
  
  
  # Define the functions needed to perform our analysis
  
  # Define the iterator wrapper function, which produces the end results
  source("Code/Resource_Dilution/Iterator_Resource_Robustness.R")
  # Define Parameters for the following iterative robustness test
  source("Code/Resource_Dilution/Iterator_Parameters.R")
  # Define functions needed for the iterator directly itself
  source("Code/Resource_Dilution/Iterator_Functions.R")
  # Define functions that specifically reference the iterator parameters
  source("Code/Resource_Dilution/Iterator_Parameter_Dependents.R")
  
  # Define functions for testing resource robustness
  source("Code/Resource_Dilution/Resource_Robustness_Functions.R")
  # Define functions to dilute resources
  source("Code/Resource_Dilution/Resource_Dilution_Functions.R")
  # Define functions to work with top_ranked CCIs and other miscellaney
  source("Code/Resource_Dilution/Ranking_and_Misc_Functions.R")
  
  
  
  
}

#------------------------------------------------------------------------------#
# 2. Get Resource Dilution Robustness Results ----------------------------------

robustness_new <- 
  wrap_resource_Iterator(
    methods_vector = c("call_connectome", "call_sca"),
    trial_run = FALSE)


