#------------------------------------------------------------------------------#
# 0. Setup ---------------------------------------------------------------------
{
  # 0.1 Overview of Goals
    # The Idea with this script is to:
    # .	run all methods on a single resource (OP)
    # .	take the topmost ranked interactions for each method -> R-zero
    # .	create a modified OP resource in which the topmost ranked interactions
    #     remain, but of the remainder x% of the interactions have been diluted 
    #     and  replaced new genes not in the resource derived from the test 
    #     data, x =10,20,40% etc.
    # .	Rerun methods on diluted omnipath resource and extract the top ranked 
    #   CCI -> R-modified
    # .	plot percentage of R-zero in R-modified over x and investigate result
  
  
}







#------------------------------------------------------------------------------#
# 1. Setup ---------------------------------------------------------------------
{
  require(tidyverse)
  require(Seurat)
  require(liana)
  require(lubridate)
  
  
  # Define Parameters for the following iterative robustness test
  source("Code/Resource_Dilution_Robustness/Iterator_Parameters.R")
  # Define functions needed for the iterator itself
  source("Code/Resource_Dilution_Robustness/Iterator_Functions.R")
  # Define iterator functions that specifically reference the iterator parameters
  source("Code/Resource_Dilution_Robustness/Iterator_Parameter_Dependants.R")
  
  # Define functions for testing resource robustness
  source("Code/Resource_Dilution_Robustness/Resource_Dilution_Functions.R")
  # Define functions to dilute resources 
  source("Code/Resource_Dilution_Robustness/Ranking_and_Misc_Functions.R")
  # Define functions to work with top_ranked CCIs and other miscellanea
  source("Code/Resource_Dilution_Robustness/Resource_Robustness_Functions.R")
  
  
  # How many permutations of dilution should be performed?
  master_seed_list <- create_Params()$master_seed_list
  time_of_run      <- create_Params()$time_of_run


}   



#------------------------------------------------------------------------------#
# 2. Iterate resource_Robustness() ---------------------------------------------
{
  # resource_Robustness is an entire script that can be iterated as a function
  # There is randomness in dilution. Each master seed passed to 
  # dilute_Resource() gives us one permutation of many theoretically possible
  # dilutions. By iterating over master_seed, we can produce many permutations
  # and tally up their results. In this way, master_seed serves as an index too.
  
  # Apply resource_Robustness(), provide every argument but master_seed
  collated_robustness_results <- lapply(master_seed_list, 
                    wrap_resource_Robustness,
                    time_of_run = time_of_run)
  
                    
}



#------------------------------------------------------------------------------#
# 3. Reformatting Results --------------------------------------------------------
{
  # In this segment we extract the data from the results object, which is poorly
  # formatted by default, and put it into a more appropriate hierarchy. We then
  # extract the most relevant sublists for the rest of the analysis.

  # We name our restructured results more informatively, then extract the most
  # relevant sublists from them for the rest of the analysis
  collated_robustness_results <- 
    reformat_Results(results = collated_robustness_results)

}




#------------------------------------------------------------------------------#
# 5. Collate all iterations of  top_ranks_analysis -----------------------------
{
  collated_top_ranks_overlap <- extract_top_ranks(collated_robustness_results)
}



#------------------------------------------------------------------------------#
# 6. Plotting of Aggregate Results ---------------------------------------------
{
  # Preparing Plotting Inputs
  {
    # We format 
    alpha <- 0.4
    
    # The dilution proportion and overlap are clearer in percentage
    # NAs can't be displayed in the plot anyway and cause uneccesary warnings
    # And we rename the methods from the liana++ internal string to their 
    # official names.
    tr_overlap_for_plot <-  collated_top_ranks_overlap  %>%
      as.data.frame()                             %>%
      mutate(dilution_prop = dilution_prop * 100) %>%
      mutate(Overlap       = Overlap       * 100) %>%
      as_tibble()                                 %>%
      drop_na()                                   %>%
      mutate("Method" = recode(Method,
                               "call_connectome" = "Connectome",
                               "squidpy"         = "CellPhoneDB",
                               "call_natmi"      = "NATMI",   
                               "call_italk"      = "iTALK", 
                               "call_sca"        = "SingleCellSignalR", 
                               "cellchat"        = "CellChat")) 
     
    # automatically generate a plot description based on wrap_resource_Robustness
    # defaults and the tr_overlaps_for_plot data structure
    
    plotting_caption <- 
      auto_plot_Description(tr_overlap_for_plot, 
                            formals(wrap_resource_Robustness),
                            formals(summarise_Metadata), 
                            time_of_run)
    

  }

  
  # Generating and printing Plots
  {
    plot_line <- overlap_line_Plot(tr_overlap_for_plot)
     
    
    plot_box <- overlap_box_Plot(tr_overlap_for_plot,
                                 plotting_caption)
    
    print(plot_line)
    print(plot_box)
  }
  
  # Removing Clutter
  rm(tr_overlap_for_plot, alpha, plotting_caption)
  
}

#------------------------------------------------------------------------------#
# 4. Reformatting Runtime and Metadata -----------------------------------------
{
  # 4.1 Calculate Runtime
  {
    # By flattening we get one long seqence of time points and the portions of
    # the script they belong to
    runtime <- collated_robustness_results$runtime %>% 
      flatten() %>%
      calculate_Runtime()


  }
  
  #4.2 format metadata
  metadata <- summarise_Metadata(runtime, 
                                 time_of_run,
                                 formals(wrap_resource_Robustness),
                                 master_seed_list)
  
  rm(runtime, master_seed_list)
  
}

#------------------------------------------------------------------------------#
# 7. Saving Results ------------------------------------------------------------
{
  # Save the plot automatically to the outputs folder, if desired
  if (formals(summarise_Metadata)$save_results == TRUE) {
    
    save_Results(dilution_params = formals(wrap_resource_Robustness),
                 meta_params     = formals(summarise_Metadata),
                 time_of_run     = time_of_run)
    
  }
}
