#------------------------------------------------------------------------------#
# 0. Introduction and Goals ----------------------------------------------------
{
  # This is the cluster . It funnels data into a robustness test for
  # various dilutions (resource_Robustness()), repeats this test many times, and
  # collates the resulting information and presents it to the user. The code to
  # accomplish this is presented in a  wrapper with default arguments below.
  
  # Running parameters can be left a their defaults or altered when the function
  # is called. To learn more, check the Run_Iterator.R script.
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
  
  # Define general functions for data processing in the iterator
  source("Code/Cluster_Reshuffling/Iterator_Processing_Functions.R")
  # Define functions for plotting the iterator results
  source("Code/Cluster_Reshuffling/Iterator_Plotting.R")
  # Define functions for capturing metadata and saving iterator results
  source("Code/Cluster_Reshuffling/Iterator_Metadata_and_Saves.R")
  
  source("Code/Cluster_Reshuffling/CD_Reshuffler.R")
  source("Code/Cluster_Reshuffling/CR_Iterator_Helpers.R")
  source("Code/Cluster_Reshuffling/Liana_wrapper.R")
  
  
  
  
  
  
}

#------------------------------------------------------------------------------#
# 1. Define wrap_resource_Iterator() -----------------------------------------
{
  #' Iterator
  #'
  #' @description This function iteratively evaluates the robustness of CCI
  #' inference methods by running them on the same testdata with somewhat
  #' reshuffled cluster annotations.
  #'
  
  # turn this into wrapper
  # wrap_cluster_Iterator <-
  #  function(
  number_seeds      <-
    3            # how many seeds should we iterate over
  testdata_type     <- "liana_test" # seurat_pbmc or liana_test
  mismatch_props    <- c(seq(0.40, 1.00, 0.40))
  
  number_ranks <- list(
    "call_connectome" = 20,
    "call_natmi"      = 20,
    "call_italk"      = 20,
    "call_sca"        = 20,
    "cellchat"        = 20,
    "squidpy"         = 20
  )
  
  methods_vector <- c('call_connectome' ,
                      # 'call_natmi'      ,
                      'call_italk'      ,
                      'call_sca'        #,
                      # 'cellchat'        ,
                      # 'squidpy'
    )
    
    sink_output     <- FALSE   # TRUE or FALSE
    liana_warnings  <- "divert" # TRUE, FALSE, or "divert"
    
    save_results    <- TRUE
    trial_run       <- FALSE
    
    
    
    
    cellchat_nperms <- 10      # default 100 for real data
    
    bundled_outputs <- c(
      "top_ranks_analysis",
      "runtime"          
    )
    
    master_outputs <- c(
      "collated_top_ranks_overlap",
      "plot_box",
      "plot_line",
      "collated_robustness_results",
      "metadata"
    )                                  

}    
    


#----------------------------------------------------------------------------#
# 1.1 Generate Parameters  ---------------------------------------------------
{
  # Retrieve either seurat_pbmc or liana_test data
  testdata <- extract_Testdata(testdata_type = testdata_type)
  
  
  
  # Format a named list of seeds, it contains as many seeds as the user
  # specified, from 1 to n, and each entry has an appropriate name, "Seed_n"
  # This list can be used to iterate over for every seed.
  master_seed_list <- as.list(1:number_seeds)
  
  names(master_seed_list) <-
    map(master_seed_list, function(seed) {
      # Name each element of master_seed_list appropriately
      str_glue("Seed_", seed)
      
    })
  
  
  # This is our second iterable structure, we iterate over every mismatch
  # proportion the user specified.
  
  # Convert mismatch_props to a list and name it, creating a named list
  mismatch_props <- as.list(mismatch_props)
  
  names(mismatch_props) <- map(mismatch_props, function(prop) {
    # Name every dilution proportion
    str_glue("Reshuffle_", as.character(prop * 100))
    
  })
  
  # This is our third iterable structure, we iterate over every method
  # the user specified.
  
  # Convert methods_vector to a list and name it, creating a named list
  methods_list <- as.list(methods_vector)
  
  names(methods_list) <- methods_vector
  
  
  
  # depending on the testdata type, the name of the column in the metadata
  # that has the cluster annotations is different. We need the name of that
  # column for multiple applications
  if (testdata_type == "seurat_pbmc") {
    cluster_col <- "seurat_clusters"
    
  } else if (testdata_type == "liana_test") {
    cluster_col <- "seurat_annotations"
    
  } else {
    stop("Testdata type is not recognized!")
    
  }
  
  
  # Format the Sys.time() of the run to not contain characters that are bad to
  # have in save file names. We will later use this to tag file names and
  # plots so they can be grouped according to run, and all have unique names.
  time_of_run <-  Sys.time() %>%
    as.character()    %>%
    gsub(':', '-', .) %>% # save files can't have colons
    gsub(' ', '_', .) %>% # save files shouldn't have spaces
    str_sub(1 , nchar(.) - 3) # the code never runs in under a minute, so the
  # number of seconds isn't valuable information.
  
  
  

  
  
  
}



#----------------------------------------------------------------------------#
# 1.2 Create Reshuffled Cluster Annotations ----------------------------------
{
  #
  reshuffled_clusters <- lapply(
    mismatch_props,
    wrap_Shuffler,
    master_seed_list  = master_seed_list,
    metadata          = testdata@meta.data,
    cluster_col       = cluster_col
  ) %>%
    append(list("Reshuffle_0" = testdata@meta.data), .)
  
  
  
}

#----------------------------------------------------------------------------#
# 1.3 Run LIANA --------------------------------------------------------------

# On original and reshuffled  metadata

liana_results <-
  iterate_liana_wrap(
    master_seed_list = master_seed_list,
    mismatch_props   = mismatch_props,
    testdata         = testdata,
    methods_vector   = methods_vector,
    
    liana_warnings   = liana_warnings,
    warning_logfile  = warning_logfile,
    
    cellchat_nperms  = cellchat_nperms
  )  

runtime <- liana_results$runtime
liana_results$runtime <- NULL

liana_results <- liana_results %>%
  map_depth(., .depth = 1, transpose) %>%
  map_depth(., .depth = 0, transpose)


if(liana_warnings == "divert") {
  
  # let the user know where to find the log
  cat(str_wrap(str_glue("LIANA warnings saved at ~/", 
                        warning_logfile, "."), width = 60), " \n\n")
  
  
}



#----------------------------------------------------------------------------#
# 1.4 Compare Top-Ranked Predictions -----------------------------------------

top_ranks <-
  map(methods_list, function(method) {
    
    if(method != "cellchat") {
      top_ranks_for_method <- liana_results[[method]] %>%
        map_depth(
          .,
          .depth = 2,
          get_top_ranks_clust,
          method = method,
          top_n = number_ranks[[method]]
        )
    } else if(method == "cellchat") {
      top_ranks_for_method <- liana_results[[method]] %>%
        map_depth(
          .,
          .depth = 2,
          get_top_ranks_clust,
          method = method,
          top_n = number_ranks[[method]],
          with_ties = TRUE
        )
    }

  })  %>%
  map_depth(., .depth = 3, format_top_ranks)


overlaps <- map(methods_list, function(method) {
  
  if(method != "cellchat") {
    
    overlaps_for_method <-  top_ranks[[method]] %>%
      map_depth(.,
                .depth = 2,
                rank_overlap,
                main_ranks = top_ranks[[method]]$Reshuffle_0[[1]],
                verbose = FALSE)
    
    
    
  } else if(method == "cellchat") {

    overlaps_for_method <-  top_ranks[[method]] %>%
      map_depth(.,
                .depth = 2,
                cellchat_rank_overlap,
                main_ranks = top_ranks[[method]]$Reshuffle_0[[1]],
                verbose = FALSE)
    
  }
  
  
})


# reformatting overlap as a tibble
top_ranks_overlap <- as_tibble(overlaps)         %>%
  mutate("mismatch_prop" = c(0, mismatch_props)) %>%
  unnest(cols = all_of(methods_vector))          %>%
  unnest(cols = everything())                    %>%
  relocate("mismatch_prop")                      %>%
  pivot_longer(cols = !(starts_with("mismatch_prop")), names_to = "Method") %>%
  arrange(Method)                                %>%
  rename("Overlap" = value)



# removing superfluous values
rm(overlaps)


#----------------------------------------------------------------------------#
# 1.4 Plotting of Collated Results -------------------------------------------
{
  # Here we visualize the overlap between top ranked CCI predictions as they
  # change with the dilution of OmniPath. once as a boxplot, and once as a
  # scatter/line plot.
  
  
  
  # We reformat the collated_top_ranks_overlap tibble so its more suitable for
  # plotting
  
  # The dilution proportion and overlap are clearer in percentage
  # NAs can't be displayed in the plot anyway and cause uneccesary warnings
  # Rename the methods from the LIANA++ internal string to their official name
  tr_overlap_for_plot <- top_ranks_overlap  %>%
    as.data.frame()                             %>%
    mutate(mismatch_prop = mismatch_prop * 100) %>% # proportion to percent
    mutate(Overlap       = Overlap       * 100) %>% # proportion to percent
    as_tibble()                                 %>%
    mutate(
      "Method" = recode(
        Method,
        "call_connectome" = "Connectome",
        "squidpy"         = "CellPhoneDB",
        "call_natmi"      = "NATMI",
        "call_italk"      = "iTALK",
        "call_sca"        = "SingleCellSignalR",
        "cellchat"        = "CellChat"
      )
    ) # renaming
  
  # To directly be able to associate the box plot with the settings that
  # produced it, we automatically generate a plot description
  plotting_caption <-
    clust_plot_Description(
      mismatch_props   = mismatch_props,
      trial_run        = trial_run,
      testdata_type    = testdata_type,
      master_seed_list = master_seed_list,
      number_ranks     = number_ranks,
      time_of_run      = time_of_run
    )
  
  
  
  # Generate our plots with functions.
  plot_line <- 
    overlap_line_Plot(tr_overlap_for_plot,
                      plotting_caption,
                      x_axis_var = "mismatch_prop",
                      x_axis_name = "Mismatch of Cluster Annotations [%]")
  
  plot_box <- 
    overlap_box_Plot(tr_overlap_for_plot,
                     plotting_caption,
                     x_axis_var = "mismatch_prop",
                     x_axis_name = "Mismatch of Cluster Annotations [%]")
  
  # Print out visualizations
  print(plot_line)
  print(plot_box)
  
  
  
  
  # Removing Clutter
  rm(tr_overlap_for_plot, plotting_caption)
  
  
}



#----------------------------------------------------------------------------#
# 1.5 Capturing Script Metadata ----------------------------------------------
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
  
  # Delete the old unordered runtime from the collated robustness results 
  # bundle
  collated_robustness_results$runtime <- NULL
  
  
  # We then summarise the above information and more metadata into a single 
  # object
  metadata <- summarise_Metadata(number_seeds      = number_seeds,
                                 master_seed_list  = master_seed_list,
                                 testdata_type     = testdata_type,
                                 feature_type      = feature_type, 
                                 preserve_topology = preserve_topology,    
                                 dilution_props    = dilution_props,
                                 number_ranks      = number_ranks ,
                                 methods_vector    = methods_vector,
                                 
                                 sink_output       = sink_output,    
                                 liana_warnings    = liana_warnings,
                                 
                                 cellchat_nperms   = cellchat_nperms,       
                                 bundled_outputs   = bundled_outputs,
                                 master_outputs    = master_outputs,
                                 
                                 
                                 save_results = save_results,
                                 trial_run    = trial_run,
                                 
                                 runtime      = runtime,
                                 time_of_run  = time_of_run)
  
  # Now that these objects are stored in the metadata object, we can remove
  # this clutter from the environment.
  rm(runtime, master_seed_list)
  
  
}



#----------------------------------------------------------------------------#
# 1.6 Packaging Results to return them ---------------------------------------
{
  #In order to save and return our results we package it in a succinct object
  iterator_results <- 
    list(
      "collated_top_ranks_overlap"  = collated_top_ranks_overlap,
      "plot_box"  = plot_box,
      "plot_line" = plot_line,
      "collated_robustness_results" = collated_robustness_results,
      "metadata"  = metadata
    )
  
  # Filter our results by the master_outputs the user wants to retrieve
  # Usually´all the data is requested so this step doesn't chaneg anything.
  iterator_results <- iterator_results[master_outputs]
  
  # Get rid of clutter we already summarized in different objects
  rm(collated_top_ranks_overlap, collated_robustness_results, metadata)
}



#----------------------------------------------------------------------------#
# 1.7 Saving Results ---------------------------------------------------------
{
  # In this segment we save the plots and environment to the outputs folder,
  # if it's specified in by the user.
  if (save_results == TRUE) {
    
    save_Results(plot_box  = plot_box,
                 plot_line = plot_line, 
                 iterator_results = iterator_results,
                 
                 preserve_topology  = preserve_topology,
                 testdata_type      = testdata_type,
                 feature_type       = feature_type,
                 number_ranks       = number_ranks,
                 time_of_run        = time_of_run,
                 trial_run          = trial_run)
    
  }
  
  
  # And now return our results to the user
  return(iterator_results)
  
  
}










