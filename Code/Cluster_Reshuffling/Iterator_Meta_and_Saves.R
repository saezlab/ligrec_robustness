#------------------------------------------------------------------------------#
# 0. Introduction and Goals ----------------------------------------------------
{}



#------------------------------------------------------------------------------#
# 1. Define Functions ----------------------------------------------------------

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

# clust_plot_Description()
{
  #' Automatically creates a verbose caption of a top ranks overlap plot
  #' 
  #' @param top_ranks_overlap As a tibble in the form of a extract_top_ranks 
  #' output, though ideally it will be preprocessed for plotting (better method 
  #' names, no NAs, etc.). This is the top_ranks_overlap that would be plotted
  #' with this caption. The function takes data from the tibble's general 
  #' structure to describe it accurately. 
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
  #' @return A verbose caption describing the parameters used to generate the 
  #' results in the plot.
  
  clust_plot_Description <- function(mismatch_props,
                                     trial_run,
                                     testdata_type,
                                     master_seed_list,
                                     number_ranks,
                                     time_of_run) {
    
    
    
    ## 
    {
      
      if (length(mismatch_props) > 1) {
        
        reshuffle_comment <- str_glue(
          "This plot was created using the ",
          testdata_type,
          " data. ",
          "The cluster annotations were reshuffled in ",
          (mismatch_props[[2]]-mismatch_props[[1]]) *100, 
          " % intervals to a maximum of ",
          max(unlist(mismatch_props))*100,
          " %. \n",
          "When an annotation was being reshuffled, ",
          "it was replaced by a random sign from ",
          "all annotations that did not match itself."
        )
        
      } else {
        reshuffle_comment <- str_glue(
          "This plot was created using the ",
          testdata_type,
          " data. ",
          "The cluster annotations were reshuffled to ",
          (mismatch_props[[1]]) *100, 
          " %. \n",
          "When an annotation was being reshuffled, ",
          "it was replaced by a random sign from ",
          "all annotations that did not match itself."
        )
      }
      
      
    }
    
    
    ## Nperms and top_ranks comment
    {
      
      top_ranks_permutations_comment <-
        str_glue(
          "The overlap was compared between the ",
          median(unlist(number_ranks)),
          " highest ranked interactions over ",
          length(master_seed_list),
          " permutations."
        )
      
    }
    
    
    ## Date and time comment
    time_comment <- str_glue("Generated at ",
                             time_of_run,
                             ".")
    
    
    ## Assemple plotting caption
    plotting_caption <-
      str_glue(
        reshuffle_comment,
        "\n\n",
        top_ranks_permutations_comment,
        "\n",
        time_comment
      )
    
    
    ## Add addendum if trial run
    if (trial_run == TRUE) {
      plotting_caption <-
        str_glue(plotting_caption, "   --   [TRIAL RUN]")
    }
    
    return(plotting_caption)
  }  # end of function
  
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

# clust_save_Results()
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
      path = "Outputs/Cluster_Reshuffling"
    )
    
    ggsave(
      plot = plot_line,
      line_plot_png_name,
      height = 9.00,
      width = 8.00,
      path = "Outputs/Cluster_Reshuffling"
    )
    
    # Save R environment and all the results within it
    save(iterator_results, file = iterator_results_save_path)
    
    # Let the user know where everything was stored.
    cat(str_wrap(str_glue("Box Plot saved at ~/Outputs/Cluster_Reshuffling/",
                          box_plot_png_name, "."), width = 60), "\n\n")
    
    cat(str_wrap(str_glue("Line Plot saved at ~/Outputs/Cluster_Reshuffling/",
                          line_plot_png_name, "."), width = 60), "\n\n")
    
    cat(str_wrap(str_glue("Iterator Results saved at ~/",
                          iterator_results_save_path, "."), width = 60), "\n\n")
  }  # end of function
}