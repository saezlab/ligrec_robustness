#------------------------------------------------------------------------------#
# 0. Introduction and Goals ----------------------------------------------------
{
  # Welcome!
  
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
  
  
  # This is the robustness iterator. It funnels data into a robustness test for
  # various dilutions (resource_Robustness()), repeats this test many times, and
  # collates the resulting information and presents it to the user. This is the
  # script that needs to be executed to obtain the resource dilution robustness
  # for each LIANA method. All the other scripts are infrastructure for this
  # script. To set the parameters for this script, go to Iterator_Params.R.
  
  
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
  
  # Define Parameters for the following iterative robustness test
  source("Code/Resource_Dilution_Robustness/Iterator_Parameters.R")
  # Define functions needed for the iterator directly itself
  source("Code/Resource_Dilution_Robustness/Iterator_Functions.R")
  # Define functions that specifically reference the iterator parameters
  source("Code/Resource_Dilution_Robustness/Iterator_Parameter_Dependants.R")
  
  # Define functions for testing resource robustness
  source("Code/Resource_Dilution_Robustness/Resource_Dilution_Functions.R")
  # Define functions to dilute resources
  source("Code/Resource_Dilution_Robustness/Ranking_and_Misc_Functions.R")
  # Define functions to work with top_ranked CCIs and other miscellanea
  source("Code/Resource_Dilution_Robustness/Resource_Robustness_Functions.R")
  
  
  # Generate parameters that need to be objects in the environment to work.
  # Check Iterator_Params.R for more information.
  master_seed_list <- create_Params()$master_seed_list
  time_of_run      <- create_Params()$time_of_run
  
  
}



#------------------------------------------------------------------------------#
# 2. Iterate resource_Robustness() ---------------------------------------------
{
  # resource_Robustness is a function that performs a single test of robustness
  # by comparing unmodified LIANA predictions with ones run on diluted OP
  # resources. Since there is randomness to resource dilution, we iterate over
  # multiple starting seeds in master_seed_list, and will collate these samples
  # later.
  
  # Apply resource_Robustness() wrapper, provide every argument but master_seed
  # To modify the defaults of the wrapper, go to Iterator_Params.R
  collated_robustness_results <- lapply(master_seed_list,
                                        wrap_resource_Robustness,
                                        time_of_run = time_of_run)
  
  
}



#------------------------------------------------------------------------------#
# 3. Reformatting Results ------------------------------------------------------
{
  # In this segment we extract the data from the results object, which is poorly
  # formatted by default, and put it into a more appropriate hierarchy. We then
  # extract the most important information from it for visualization.
  
  # We reformat the results so they are grouped more intuitively
  collated_robustness_results <-
    reformat_Results(results = collated_robustness_results)
  
  # We extract the top_ranks_overlap data and turn it into a convenient tibble
  collated_top_ranks_overlap <-
    extract_top_ranks(collated_robustness_results)
  
  
}



#------------------------------------------------------------------------------#
# 4. Plotting of Collated Results ---------------------------------------------
{
  # Here we visualize the overlap between top ranked CCI predictions as they
  # change with the dilution of OmniPath. once as a boxplot, and once as a
  # scatter/line plot.
  
  
  
  # We reformat the collated_top_ranks_overlap tibble so its more suitable for
  # plotting
  
  # The dilution proportion and overlap are clearer in percentage
  # NAs can't be displayed in the plot anyway and cause uneccesary warnings
  # Rename the methods from the LIANA++ internal string to their official name
  tr_overlap_for_plot <-  collated_top_ranks_overlap  %>%
    as.data.frame()                             %>%
    mutate(dilution_prop = dilution_prop * 100) %>% # proportion to percent
    mutate(Overlap       = Overlap       * 100) %>% # proportion to percent
    as_tibble()                                 %>%
    drop_na()                                   %>% # no NAs
    mutate("Method" = recode(Method,
                             "call_connectome" = "Connectome",
                             "squidpy"         = "CellPhoneDB",
                             "call_natmi"      = "NATMI",
                             "call_italk"      = "iTALK",
                             "call_sca"        = "SingleCellSignalR",
                             "cellchat"        = "CellChat")) # renaming
  
  # To directly be able to associate the box plot with the settings that
  # produced it, we automatically generate a plot description based on the
  # wrap_resource_Robustness defaults, the summarise_Metadata defaults, the
  # date and time of the run and the tr_overlaps_for_plot data structure
  plotting_caption <-
    auto_plot_Description(
      tr_overlap_for_plot,
      formals(wrap_resource_Robustness),
      formals(summarise_Metadata),
      time_of_run)
  
  
  
  # Generate our plots with functions.
  plot_line <- overlap_line_Plot(tr_overlap_for_plot)
  
  plot_box <- overlap_box_Plot(tr_overlap_for_plot,
                               plotting_caption)
  
  # Print out visualizations
  print(plot_line)
  print(plot_box)
  
  
  
  
  # Removing Clutter
  rm(tr_overlap_for_plot, plotting_caption)
  
  
}



#------------------------------------------------------------------------------#
# 5. Capturing Script Metadata -------------------------------------------------
{
  # In this segment, we summarize the metadata of the resource_Robustness run
  # and the iterator in general. When troubleshooting or reproducing results,
  # this information will be useful to the user.
  
  # rescource_Robustness return the points in time certain processes started,
  # we use this long list of named points in time to create a succinct and
  # informative runtime overview.
  runtime <- collated_robustness_results$runtime %>%
    flatten() %>%
    calculate_Runtime()
  
  
  # We then summarise the above information and more into a single object
  metadata <- summarise_Metadata(runtime,
                                 time_of_run,
                                 formals(wrap_resource_Robustness),
                                 master_seed_list)
  
  # Now that these objects are stored in the metadata object, we can remove
  # this clutter from the environment.
  rm(runtime, master_seed_list)
  
  
}



#------------------------------------------------------------------------------#
# 6. Saving Results ------------------------------------------------------------
{
  # In this segment we save the plots and environment to the outputs folder,
  # if it's specified in Iterator_Params.R
  if (formals(summarise_Metadata)$save_results == TRUE) {
    
    save_Results(dilution_params = formals(wrap_resource_Robustness),
                 meta_params     = formals(summarise_Metadata),
                 time_of_run     = time_of_run)
    
  }
  
  
}
