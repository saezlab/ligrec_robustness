#------------------------------------------------------------------------------#
# 0. Setup ---------------------------------------------------------------------
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


source("Code/Ranking_Misc_Functions.R")
source("Code/Resource_Dilution_Functions.R")
source("Code/Resource_Robustness_Functions.R")
source("Code/Resource_Iterator_Functions.R")




#------------------------------------------------------------------------------#
# 1. Script  Parameters --------------------------------------------------------
{

    # How many permutations of dilution should be performed?
    master_seed_list <- as.list(c(1:2))
    
    # A tag for your results that will mark them as a test run or serious data.
    run_mode <- "trial_run" # select between "trial_run" and "real"
    
    # should the results automatically be saved? TRUE or FALSE
    save_results <- TRUE # Saved under automatically generated name in Outputs
    
    ## Process Master Seed List
    {
      # By naming each seed we can use this to label data conveniently later
      # By formatting seeds as a list we can lapply over them for an easy-to-
      # untangle results format.
      seed_names <- c()
      
      # Name each element of master_seed_list appropriately name it
      for (seed in master_seed_list) {
        seed_names <- 
          c(seed_names, str_glue("Seed_", seed))
      }
      
      names(master_seed_list) <- seed_names
      
      # Remove clutter
      rm(seed_names, seed)
    }


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
  results <- lapply(master_seed_list, 
                    wrap_resource_Robustness)
                    
}



#------------------------------------------------------------------------------#
# 3. Reformatting Results --------------------------------------------------------
{
  # In this segment we extract the data from the results object, which is poorly
  # formatted by default, and put it into a more appropriate hierarchy. We then
  # extract the most relevant sublists for the rest of the analysis.

  # We name our restructured results more informatively, then extract the most
  # relevant sublists from them for the rest of the analysis
  resource_Robustness_results <- reformat_Results(results = results)
  
  top_ranks_analysis <- resource_Robustness_results$top_ranks_analysis
  runtime            <- resource_Robustness_results$metadata
  
  
  # Remove unnecessary clutter from the environment.
  rm(results)
}



#------------------------------------------------------------------------------#
# 4. Reformatting Runtime and Metadata -----------------------------------------
{
  # 4.1 Calculate Runtime
  {
    
    
    # By flattening we get one long seqence of time points and the portions of
    # the script they belong to
    runtime <- flatten(runtime)
    
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
    
    
  }
  
  # 4.2 Reformat Metadata
  # {
  #   # Create a metadata subsection of script_params and assign it runtime
  #   script_params[["metadata"]][["runtime"]]      <- runtime
  #   
  #   # If output was sunk, store the file path to metadata
  #   if(script_params$sink_output == TRUE) {
  #     
  #     script_params$metadata[["sink_logfile"]] <- 
  #       script_params[["sink_logfile"]]
  #     
  #   }
  #   
  #   # If warnings were diverted, store the file path to metadata 
  #   if(script_params$liana_warnings == "divert") {
  #     
  #     script_params$metadata[["warning_logfile"]] <- 
  #       script_params[["warning_logfile"]]
  #     
  #   }
  #   
  #   # Remove any file paths and session info outside of metadata
  #   script_params[["sink_logfile"]]    <- NULL
  #   script_params[["warning_logfile"]] <- NULL
  #   
  #   # Remove runtime now that it's a part of script_params$metadata
  #   rm(runtime)
  # }
  
}


#------------------------------------------------------------------------------#
# 5. Aggregate top_ranks_analysis ----------------------------------------------
{
  
  where_overlap <- str_detect(names(top_ranks_analysis), "Overlap")
  
  collated_top_ranks_overlap <- top_ranks_analysis[where_overlap] %>%
    bind_rows() 
  
  
  collated_top_ranks_overlap <- collated_top_ranks_overlap %>%
    arrange(dilution_prop) %>%
    pivot_longer(cols = !(starts_with("dilution_prop")), names_to = "Method") %>%
    arrange(Method) %>%
    rename("Overlap" = value) 
  
  rm(seed_assignment, where_overlap)
  

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
                               "squidpy"    = "CellPhoneDB",
                               "call_natmi"      = "NATMI",   
                               "call_italk"      = "iTALK", 
                               "call_sca"        = "SingleCellSignalR", 
                               "cellchat"        = "CellChat")) 
     
    # Automatically assemble a file name and plot subtitle
    if (formals(wrap_resource_Robustness)$preserve_topology == FALSE) {
      
      topology_comment <- "random_Dilute()"
      
    } else if (formals(wrap_resource_Robustness)$preserve_topology == TRUE) {
      
      topology_comment <- "preserve_Dilute()"
      
    }
    
    dilution_overview <- count(tr_overlap_for_plot, dilution_prop)
    
    if (nrow(dilution_overview) > 1) {
      dilution_comment <- str_glue("The dilution occured in ", 
                                   dilution_overview$dilution_prop[2] -
                                     dilution_overview$dilution_prop[1],
                                   " % increments up to a maximum of ",
                                   max(tr_overlap_for_plot$dilution_prop),
                                   " %. ")
    } else {
      stop("Expected at least two dilution proportions in input (0, and one ",
           "more. But found only one instead, namely ",
           dilution_overview$dilution_prop)
    }
    
    if (length(unique(dilution_overview$n)) != 1) {
      stop("There should be an equal number of samples for every dilution, ",
           "but there is not.")
    }
    
    top_ranks_vector <- 
      unlist(as.list(formals(wrap_resource_Robustness)$number_ranks)[-1])
    
    permutations_overview <- tr_overlap_for_plot %>%
      filter(dilution_prop == 0) %>%
      count(Method)
    
    if (length(unique(permutations_overview$n)) != 1) {
      stop("There should be an equal number of samples for each method at , ",
           "dilution proportion 0, but there is not.")
    }
    
    
    top_ranks_permutations_comment <-
      str_glue(
        "The overlap was compared between the ",
        median(top_ranks_vector),
        " highest ranked interactions over ",
        permutations_overview$n[1],
        " permutations."
      )
    
    plotting_caption <- 
      str_glue("This plot was created using the ",
               formals(wrap_resource_Robustness)$testdata_type,
               " data. Dilution was performed using ",
               formals(wrap_resource_Robustness)$feature_type,
               " features and the ",
               topology_comment,
               " function. \n",
               dilution_comment,
               "\n\n",
               "The overlap was compared between the ",
               median(top_ranks_vector),
               " highest ranked interactions over ",
               dilution_overview$n[1],
               " permutations."
               )
               
    
    if (run_mode == "trial_run") {
      plotting_caption <- 
        str_glue(plotting_caption, "   --   [TRIAL RUN]")
    }
    
    

  }

  
  # Generating and printing Plots
  {
    plot_line <- 
      ggplot(data = tr_overlap_for_plot, aes(dilution_prop,
                                             Overlap,
                                             group = Method,
                                             color = Method)) + 
      
      
      geom_point(alpha = alpha) +
      stat_summary(alpha = 0.6,
                   fun   = mean, 
                   geom  = "line") +
      
      
      scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0,100)) +
      scale_x_continuous(breaks = seq(0, 100, 20), limits = c(0,100)) +
      
      ggtitle("Robustness of Method Predictions") +
      ylab("Overlap of Top Ranks [%]") +
      xlab("Dilution of Resource [%]") +
      labs(subtitle = "Point scatter plot.",
           caption = plotting_caption,
           color = "Method") +
      
      theme_bw() +
      
      theme(plot.caption = element_text(hjust = 0),
            legend.position = "bottom")
    
    
    
    
    
    
    plot_box <- 
      ggplot(data = tr_overlap_for_plot, aes(x = dilution_prop, 
                                             y = Overlap, 
                                             group = dilution_prop,
                                             color = Method)) + 
      geom_boxplot(outlier.shape = NA) + 
      geom_point(alpha = alpha) + 

      scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0,100)) +
      scale_x_continuous(breaks = seq(0, 100, 20)) +

      
      ggtitle("Robustness of Method Predictions") +
      ylab("Overlap of Top Ranks [%]") +
      xlab("Dilution of Resource [%]") +
      labs(subtitle = "Boxplot by Method.",
           caption = plotting_caption,
           color = "Method") +
      
       theme_bw() + 
      
      
      theme(plot.caption = element_text(hjust = 0),
            legend.position = "bottom") +     
      
      facet_wrap(~Method, nrow = 2, ncol = 3, scales = "free")
    
    
    
    print(plot_line)
    print(plot_box)
  }
  
  # Removing Clutter
  rm(tr_overlap_for_plot, alpha, plotting_caption,topology_comment)
  
}



#------------------------------------------------------------------------------#
# 7. Saving Results ------------------------------------------------------------
{
  # Save the plot automatically to the outputs folder, if desired
  if (save_results == TRUE) {
    if (run_mode == "real") {
      test_run_comment <- ""
      
    } else if (run_mode == "trial_run") {
      test_run_comment <- "TRIAL_RUN_"
      
    }
    
    if (formals(wrap_resource_Robustness)$preserve_topology == FALSE) {
      topology_comment <- "_random_topology_"
      
    } else if (formals(wrap_resource_Robustness)$preserve_topology == TRUE) {
      topology_comment <- "_preserved_topology_"
      
    }
    
    # Define the time of run to uniquely tag every save file
    time_of_run <-  Sys.time() %>%
      as.character()       %>%
      gsub(':', '-', .)    %>% 
      gsub(' ', '_', .)
    
    top_ranks_vector <- 
      unlist(as.list(formals(wrap_resource_Robustness)$number_ranks)[-1])
    
    # Generate the filepaths to save the data under
    box_plot_png_name <-
      str_glue(
        test_run_comment,
        "Boxplot_Resource_Dilution_",
        formals(wrap_resource_Robustness)$testdata_type,
        topology_comment,
        formals(wrap_resource_Robustness)$feature_type,
        "_top",
        median(top_ranks_vector),
        "_",
        time_of_run,
        ".png"
      )
    
    line_plot_png_name <-
      str_glue(
        test_run_comment,
        "Lineplot_Resource_Dilution_",
        formals(wrap_resource_Robustness)$testdata_type,
        topology_comment,
        formals(wrap_resource_Robustness)$feature_type,
        "_top",
        median(top_ranks_vector),
        "_",
        time_of_run,
        ".png"
      )
    
    
    env_save_path <- 
      str_glue(
        "Outputs/",
        test_run_comment,
        "DilutionEnv_",
        formals(wrap_resource_Robustness)$testdata_type,
        topology_comment,
        formals(wrap_resource_Robustness)$feature_type,
        "_top",
        median(top_ranks_vector),
        "_",
        time_of_run,
        ".RData"
      )
    
    
    
    # Save both plots
    ggsave(
      plot = plot_box,
      box_plot_png_name,
      height = 7.75,
      width = 8,
      path = "Outputs"
    )
    
    ggsave(
      plot = plot_line,
      line_plot_png_name,
      height = 8.5,
      width = 8,
      path = "Outputs"
    )
    
    
    # Store save locations for plots and session info in script_params metadata
    # script_params$metadata[["box_plot_png_name"]]  <- box_plot_png_name
    # script_params$metadata[["line_plot_png_name"]] <- line_plot_png_name
    # script_params$metadata[["env_save_path"]]      <- env_save_path
    # 
    # script_params$metadata[["Session_Info"]]       <- sessionInfo
    # 
    
    
    # Save R environment and all the results within it
    save.image(file = env_save_path)
    
    # Let the user know where everything was stored.
    print(str_glue("Box Plot saved at ~/Outputs/", 
                   box_plot_png_name, "."))
    
    print(str_glue("Line Plot saved at ~/Outputs/", 
                   line_plot_png_name, "."))
    
    print(str_glue("Environment saved at ~/", 
                   env_save_path, "."))
    
  }
}
