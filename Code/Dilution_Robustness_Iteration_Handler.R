#------------------------------------------------------------------------------#
# A. Setup ---------------------------------------------------------------------
{
  # 0.1 Overview of Goals:
  {
    # The Idea with this script is to:
    # .	run all methods on a single resource (OP)
    # .	take the topmost ranked interactions for each method -> R-zero
    # .	create a modified OP resource in which the topmost ranked interactions
    #     remain, but of the remainder x% of the interactions have been diluted 
    #     and  replaced new genes not in the resource derived from the test 
    #     data, x =10,20,40% etc.
    # .	Rerun methods on diluted omnipath resource, get top ranks -> R-modified
    # .	plot percentage of R-zero in R-modified over x and investigate result
    
    
  } # end of subpoint
  
  
  # 0.2 Loading Packages and Starting Runtime
  {

    require(tidyverse)
    require(Seurat)
    require(liana)
    
    library(lubridate)
    
    # each runtime is named after the section that it marks the conclusion of
    
  } # end of subpoint
  
}



#------------------------------------------------------------------------------#
# B. Script  Parameters --------------------------------------------------------
{
  testdata_type  <- c("liana_test") # choose "liana_test" or "seurat_pbmc"
  
  # If feature_type is generic, dilution will be completed with any genes
  # in the seurat count matrix. If dilution is variable, only the variable
  # features are used for dilution.
  feature_type <- c("variable") # choose "generic" or "variable"
  
  preserve_topology <- FALSE  # TRUE = preserve_Dilute(), FALSE = random_Dilute()
  
  dilution_props <- c(seq(0.20, 0.40, 0.20)) # should be consistent between tests
  
  number_ranks   <- list("call_connectome" = 20, 
                         "call_natmi"      = 20,
                         "call_italk"      = 20,
                         "call_sca"        = 20,
                         "cellchat"        = 20)
  

  

  
  # All the methods we're using (almsot all six of liana)
  # squidpy won't be used unthetil I get it to work on windows
  methods_vector <- c('call_connectome',
                      'call_natmi', 
                      'call_italk',
                      'call_sca',
                      'cellchat')
  
  cellchat_nperms <- 10 # number of cellchat permutations, default 100
  
  run_mode <- "trial_run" # select between trial_run and real
  
  save_results <- FALSE # should results be saved?
  
  sink_output <- FALSE # If the output is sunk, all console outputs and warning
  # messages go to a txt file in Outputs folder, but you won't be able to see
  # them in the console. In essence, logs will be generated instead of console
  # outputs

}   



