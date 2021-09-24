#------------------------------------------------------------------------------#
# 0. Introduction and Goals ----------------------------------------------------
{
  # This is the Iterator_Metadata_and_Saves.R script.

  # The goal of this script is do define the more meta level functions in the 
  # RD_Robustness_Iterator.R script that don't directly process data, but 
  # instead are relevant for saving files and summarizing metadata.

}



#------------------------------------------------------------------------------#
# 1. Defining Functions---------------------------------------------------------

# calculate_Runtime()
{
  #' Converts the resource_Robustness() runtime output into a convenient tibble
  #' output
  #' 
  #' @description Takes the runtime output of resource_Robustness(), which is a 
  #' named list of checkpoints in time, and creates a  output tibble that has 
  #' the time between each checkpoint in it and the time elapsed up until that
  #' checkpoint.
  #' 
  #' @param runtime The resource_Robustness() runtime output, as a list.
  #' 
  #' @return A tibble giving an overview of the runtime.

  
  calculate_Runtime <- function(runtime) {
    
    # save the names of the time-points for later
    runtime_labels <- names(runtime)
    
    # convert run time to numeric so we can perform arithmetic operations on
    # them. In this case we need it for subtractions, to calculate the duration
    # between checkpoints
    runtime_numeric <- as.numeric(runtime)
    
    # We calculate the passage of time between checkpoints in the 
    # resource_Robustness().
    # Step duration is the duration of a step between neighboring checkpoints.
    # Time elapsed is the duration between the completion of a step and the 
    # start of the script.
    
    step_duration <- c(0) # No time has passed when the script is initialized.
    time_elapsed  <- c(0) # No time has passed when the script is initialized.
    
    # starting with the second index of runtime_numeric until the last index
    for (i in 2:length(runtime_numeric)) {
      
      # subtract the preceding checkpoint from the checkpoint at i, this is the 
      # amount of time that passed between these two checkpoints
      step_duration <- c(step_duration, 
                         runtime_numeric[[i]] - runtime_numeric[[i-1]])
      
      # subtract the very first checkpoint from the checkpoint at i, this is all
      # the time that has elapsed up until now.
      time_elapsed  <- c(time_elapsed,
                         runtime_numeric[[i]] - runtime_numeric[[1]])
      
    }
    
    # Turn seconds into time periods using lubridate and round for simplicity
    # Time periods are HH:MM:SS, which is earier to understand than just values
    # in seconds.
    step_duration <- round(seconds_to_period(step_duration))
    time_elapsed  <- round(seconds_to_period(time_elapsed))
    
    
    # summarize all the runtime data in a tibble
    runtime <- runtime               %>%
      as_tibble_col()                %>%
      unnest(cols = c(value))        %>%
      rename("Start Time" = "value") %>% 
      add_column("Step Name"      = runtime_labels, .before = 1) %>%
      add_column("Step Duration"  = step_duration) %>%
      add_column("Time Elapsed"   = time_elapsed) 
    
    
    # Get rid of clutter in the environment
    rm(runtime_numeric, 
       step_duration, 
       time_elapsed, 
       runtime_labels,
       i)
    
    # return the runtime tibble
    return(runtime)
  }
}


# clust_summarise_Metadata()
{
  #' Summarizes the metadate relevant for the Robustness Iterator
  #' 
  #' @description This function summarizes all the iterator parameters, file 
  #' names that were used, runtime data, and more into one metadata object.
  #' 
  #' @inheritParams wrap_resource_Iterator
  #' 
  #' @param master_seed_list The list of seeds that resource_Robustness() 
  #' iterated over.
  #' 
  #' @param runtime The tibble runtime output of calulate_Runtime().
  #' 
  #' @param time_of_run The char tag of the time the script started being 
  #' executed.
  #' 
  #' @return Returns a compiled list of metadata, parameters and save file 
  #' locations (if files were saved to the computer).
  
  
  
  clust_summarise_Metadata <- function(master_seed_list,
                                 mismatch_props,
                                 methods_list,
                                 
                                 testdata_type,
                                 cluster_col,
                                 number_ranks,
                                 cellchat_nperms,
                                 
                                 outputs,
                                 liana_warnings,
                                 save_results,
                                 trial_run,
                                 runtime,
                                 time_of_run,
                                 
                                 warning_logfile,
                                 line_plot_png_name,
                                 box_plot_png_name,
                                 iterator_results_save_path) {
    
    # Summarize the metadata parameters
    meta_params <- list(
      "outputs"         = outputs,
      "liana_warnings"  = liana_warnings,
      "save_results"    = save_results,
      "trial_run"       = trial_run,
      "time_of_run"     = time_of_run
    )
    
    # Summarize Save files
    save_files <- list()
    
    if(liana_warnings == "divert") {
      
      save_files <- save_files %>%
        append(list("warning_logfile" = warning_logfile))
      
    }
    
    if(save_results == TRUE) {
      
      save_files <- save_files %>%
        append(list("line_plot_png_name" = line_plot_png_name,
                    "box_plot_png_name"  = box_plot_png_name,
                    "iterator_results_save_path" = iterator_results_save_path)) 
      
    }
    

    
    # summarise all the parameters from wrap_resource_Iterator()
    reshuffle_params <- list(
      "master_seed_list" = master_seed_list,
      "mismatch_props"   = mismatch_props,
      "methods_list"     = methods_list,
      
      "testdata_type"    = testdata_type,
      "cluster_col"      = cluster_col,
      "number_ranks"     = number_ranks,
      "cellchat_nperms"  = cellchat_nperms
    )
    
    # Put all the parameters in a list
    metadata <- list(
      "runtime"          = runtime,
      "reshuffle_params" = reshuffle_params,
      "meta_params"      = meta_params,
      "save_files"       = save_files,
      "sessionInfo"      = sessionInfo()
    )
  
    
    # return the metadata.
    return(metadata)
    
  } # end of function
  
  
  
}


# save_Results()
{
  #' Saves its three arguments to custom filepaths
  #'
  #' @param plot_box Takes the boxplot generated by overlap_box_Plot as an 
  #' input. Saves it to the outputs folder under a descriptive name.
  #'
  #' @param plot_line Takes the lineplot generated by overlap_line_Plot as an 
  #' input. Saves it to the outputs folder under a descriptive name.
  #' 
  #' @param iterator_results Takes the list of results from the iterator, saves 
  #' them to a descriptive file name in the outputs folder.
  #' 
  #' @param trial_run The same parameter from wrap_resource_Iterator(). Used
  #' in the file name to mark the file.
  #'
  #' @param preserve_topology The same parameter from 
  #' wrap_resource_Iterator(). Used in the file name to mark the file.
  #' 
  #' @param testdata_type The same parameter from wrap_resource_Iterator(). 
  #' Used in the file name to mark the file.
  #' 
  #' @param feature_type The same parameter from wrap_resource_Iterator(). 
  #' Used in the file name to mark the file.
  #' 
  #' @param number_ranks The same parameter from wrap_resource_Iterator(). 
  #' Used in the file name to mark the file.
  #' 
  #' @param time_of_run The char tag of the time the script started being 
  #' executed.

  
  
  clust_save_Results <- function(plot_box,
                           plot_line,
                           iterator_results,
                           
                           line_plot_png_name,
                           box_plot_png_name,
                           iterator_results_save_path) {
    
    
    
    
    # Save both plots
    ggsave(
      plot = plot_box,
      box_plot_png_name,
      height = 7.75,
      width = 8.00,
      path = "Outputs/Cluster_Dilution"
    )
    
    ggsave(
      plot = plot_line,
      line_plot_png_name,
      height = 9.00,
      width = 8.00,
      path = "Outputs/Cluster_Dilution"
    )
    
    # Save R environment and all the results within it
    save(iterator_results, file = iterator_results_save_path)
    
    # Let the user know where everything was stored.
    cat(str_wrap(str_glue("Box Plot saved at ~/Outputs/Cluster_Dilution/",
                          box_plot_png_name, "."), width = 60), "\n\n")
    
    cat(str_wrap(str_glue("Line Plot saved at ~/Outputs/Cluster_Dilution/",
                          line_plot_png_name, "."), width = 60), "\n\n")
    
    cat(str_wrap(str_glue("Iterator Results saved at ~/",
                          iterator_results_save_path, "."), width = 60), "\n\n")
  }  # end of function
}


# cluster_auto_file_Name()
{
  #' Automatically generates a file name or file path
  #' 
  #' @param prefix As a char. What should the file name start with? It could be
  #' a folder to make it a file path, such as "Outputs/", or any other tag, or 
  #' "".
  #' 
  #' @param suffix As a char. What should the file name end with? It should be a 
  #' file extension such as ".txt" or ".RData" at minimum, but it could also be
  #' more, such as "report.txt".
  #' 
  #' @param trial_run The same parameter from wrap_resource_Iterator(). Used
  #' in the file name to mark the file.
  #'
  #' @param preserve_topology The same parameter from 
  #' wrap_resource_Iterator(). Used in the file name to mark the file.
  #' 
  #' @param testdata_type The same parameter from wrap_resource_Iterator(). 
  #' Used in the file name to mark the file.
  #' 
  #' @param feature_type The same parameter from wrap_resource_Iterator(). 
  #' Used in the file name to mark the file.
  #' 
  #' @param number_ranks The same parameter from wrap_resource_Iterator(). 
  #' Used in the file name to mark the file.
  #' 
  #' @param time_of_run The char tag of the time the script started being 
  #' executed.
  #' 
  #' @return A file name that starts with the prefix, ends with the suffix and 
  #' contains a bunch of parameter tags in between. This way the user can 
  #' identify the save file by the parameters it was set up with.
  
  
  cluster_auto_file_Name <- function(prefix,
                                     suffix,
                             
                                     trial_run,
                                     testdata_type,
                                     feature_type,
                                     number_ranks,
                                     time_of_run) {
    
    # We define individual comments related to relevant parameters and then 
    # string them all together for the save file name.
    
    # If this is a trial run, mark the save files as such
    if (trial_run == FALSE) {
      
      test_run_comment <- ""
      
    } else if (trial_run == TRUE) {
      
      test_run_comment <- "TRIAL_RUN_"
      
    }
    
    
    # ´What testdata_type was extracted and used with resource_Robustness?
    testdata_comment <-
      str_glue(testdata_type, "_")
    
    
    
    # Make a comment out of the median top_n that was considered top_ranked. 
    # As a note, it could be a different number per method, but often it will
    # be the same number for every method.
    top_ranks_comment <-
      str_glue("top", median(unlist(number_ranks)), "_", )
    
    
    # Mash all the comments together with the suffix and prefix to create our
    # custom and hopefully informative file names
    auto_file_name <-
      str_glue(
        prefix,
        test_run_comment,
        testdata_comment,
        top_ranks_comment,
        time_of_run,
        suffix
      )
    
    return(auto_file_name)
  }  # end of function
  
}

