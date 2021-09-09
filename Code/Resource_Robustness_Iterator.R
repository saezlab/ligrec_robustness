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



#------------------------------------------------------------------------------#
# 1. Script  Parameters --------------------------------------------------------
{
  # 1.1 Define Script Parameters
  {
    # More information on most of these can be found in the dilution_Robustness
    # description.
    
    # Which scRNA data set do you want to use?
    testdata_type  <- c("liana_test") # choose "liana_test" or "seurat_pbmc"
    
    # Should dilution be with any genes from testdata or most variable ones?
    feature_type <- c("variable") # choose "generic" or "variable"
    
    # TRUE = preserve_Dilute(), FALSE = random_Dilute()
    preserve_topology <- FALSE # Preserve topology after dilution or not?
    
    # Which proportions should the resources that are analysed have?
    # have at least two else the plotting fails
    dilution_props <- c(seq(0.40, 1.00, 0.40)) 
    
    # How many permutations of dilution should be performed?
    master_seed_list <- as.list(c(1:2))
    
    # Which Outputs from dilution_Robustness() do you want? Choose from:
    # outputs = c("liana_results_OP", "resources_OP", "top_ranks_OP", 
    #             "top_ranks_analysis","metadata", "testdata")
    outputs = c("liana_results_OP", "resources_OP", "top_ranks_OP",
                "top_ranks_analysis","metadata", "testdata")
    
    
    # Which methods should dilution_Robustness() use? Choose from:
    # methods_vector <- c('call_connectome', 'call_squidpy', 'call_natmi', 
    #                     'call_italk', 'call_sca', 'cellchat')
    # Squidpy doesn't work on windows.
    methods_vector <- c('call_connectome', 'call_natmi',
                        'call_italk', 'call_sca', 'cellchat')
    
    # Which top n of interactions should be considered top-ranked per method?
    number_ranks   <- list("call_connectome" = 20,
                           "call_squidpy"    = 20,
                           "call_natmi"      = 20,
                           "call_italk"      = 20,
                           "call_sca"        = 20,
                           "cellchat"        = 20)
    
    cellchat_nperms <- 10 # number of cellchat permutations, default 100
    
    # A tag for your results that will mark them as a test run or serious data.
    run_mode <- "trial_run" # select between "trial_run" and "real"
    
    # should the results automatically be saved? TRUE or FALSE
    save_results <- TRUE # Saved under automatically generated name in Outputs
    
    sink_output <- FALSE # Should the entire output be redirected to a log file?
    
    # Should liana warnings be visible in output? 
    liana_warnings <- "divert" # TRUE, FALSE, or "divert" to create a log file.
  } # end of subpoint
  
  # 1.2 Process Script Parameters
  {
    
    ## Process Dilution Proportions List
    {
      # By naming each dilution proportion we can use this to label data easily
      # By formatting proportions as a list we can lapply over them conveniently
      
      # intitialize dilution names
      dilution_names <- c()
      
      # Generate a dilution name for each proportion
      for (i in dilution_props) {
        dilution_names <- 
          c(dilution_names, str_glue("OmniPath_", as.character(i*100)))
      }
      
      # Convert dilution.proprs to a list and name it, creating a named list
      dilution_props <- as.list(dilution_props)
      names(dilution_props) <- dilution_names
      
      #remove clutter
      rm(dilution_names, i)
    }
    
    
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
    
    
    ## Initialize Script_Params
    {
      # Summarize all the above parameters into one compact item
      script_params <- list("master_seed_list"  = master_seed_list,
                            "dilution_props"    = dilution_props, 
                            "number_ranks"      = number_ranks, 
                            "feature_type"      = feature_type, 
                            "methods_vector"    = methods_vector, 
                            "testdata_type"     = testdata_type,
                            "preserve_topology" = preserve_topology,
                            
                            "outputs"         = outputs,
                            "cellchat_nperms" = cellchat_nperms, 
                            "run_mode"        = run_mode, 
                            "save_results"    = save_results,
                            "sink_output"     = sink_output,
                            "liana_warnings"  = liana_warnings,
                            "sink_logfile"    = "",
                            "warning_logfile" = "")
      # Notice sink_logfile and warning_logfile, these are both variables that
      # haven't been defined yet above, but based on whether or not they will
      # be needed they will be generated and inserted into the blank spaces here
      # later.
    }

    ## Define File Paths for Logs
    {
      # If they are required, we auto generate log filepaths here. The filepaths
      # have the current time imminetnly before the lapply in them. Each lapply
      # should be assigned to a unique log this way that has all iterations in
      # it. 
      
      # Define time of run and remove characters problematic for file names
      time_of_run <-  Sys.time() %>%
        as.character()       %>%
        gsub(':', '-', .)    %>% 
        gsub(' ', '_', .)
      
      # If necessary, create a log name, store it in script params, then remove
      # the leftover clutter
      if(sink_output == TRUE) {
        # The file name includes many script_params in it to be informative and
        # unique.
        sink_logfile <- str_glue("Outputs/Logs/Complete_Log_", 
                                 script_params$run_mode,
                                 "_",
                                 testdata_type, 
                                 "_top",
                                 as.character(median(unlist(number_ranks))),
                                 "_res",
                                 as.character(length(dilution_props)),
                                 "_",
                                 feature_type,                            
                                 "_dil_at_",
                                 time_of_run,
                                 ".txt")
        
        script_params[["sink_logfile"]] <- sink_logfile
        
        rm(sink_logfile)
      }
      
      # If necessary, create a log name, store it in script params, then remove
      # the leftover clutter
      if(liana_warnings == "divert") {
        # The file name includes many script_params in it to be informative and
        # unique.
        warning_logfile <- str_glue("Outputs/Logs/LIANA_warnings_", 
                                    script_params$run_mode,
                                    "_",
                                    testdata_type, 
                                    "_top",
                                    as.character(median(unlist(number_ranks))),
                                    "_res",
                                    as.character(length(dilution_props)),
                                    "_",
                                    feature_type,                            
                                    "_dil_at_",
                                    time_of_run,
                                    ".txt")
        
        script_params[["warning_logfile"]] <- warning_logfile
        
        rm(warning_logfile)
      }
    }

    
    ## Remove Clutter
    {
      # Since all the relevant args are stored in script_params now, we can delete
      # all the random variables still scattered in the environment.
      
      rm(master_seed_list, dilution_props, number_ranks, feature_type, 
         methods_vector, testdata_type, preserve_topology, outputs, 
         cellchat_nperms, run_mode, save_results, sink_output, liana_warnings,
         time_of_run)
    }

    
    } # end of subpoint

}   



#------------------------------------------------------------------------------#
# 2. Iterate dilution_Robustness() ---------------------------------------------
{
  # dilution_Robustness is an entire script that can be iterated as a function
  # There is randomness in dilution. Each master seed passed to 
  # dilute_Resource() gives us one permutation of many theoretically possible
  # dilutions. By iterating over master_seed, we can produce many permutations
  # and tally up their results. In this way, master_seed serves as an index too.
  
  # Apply dilution_Robustness(), provide every argument but master_seed
  results <- lapply(script_params$master_seed_list, 
                    dilution_Robustness,
                    
                    testdata_type     = script_params$testdata_type,
                    feature_type      = script_params$feature_type,
                    preserve_topology = script_params$preserve_topology,
                    dilution_props    = script_params$dilution_props,
                    number_ranks      = script_params$number_ranks,
                    outputs           = script_params$outputs,
                    
                    methods_vector    = script_params$methods_vector,
                    cellchat_nperms   = script_params$cellchat_nperms,
                    sink_output       = script_params$sink_output,
                    liana_warnings    = script_params$liana_warnings,
                    
                    sink_logfile      = script_params$sink_logile,
                    warning_logfile   = script_params$warning_logfile)
}



#------------------------------------------------------------------------------#
# 3. Extracting Results --------------------------------------------------------
{
  # In this segment we extract the data from the results object, which is poorly
  # formatted by default, and put it into a more appropriate hierarchy. We then
  # extract the most relevant sublists for the rest of the analysis.
  
  # Initiate the restructured results list.
  restructured_results <- list()
  
  # The two outputs mentioned here can be formatted the same way
  # But only execute this code if that output is actually in results
  for(output in intersect(script_params$outputs, c("liana_results_OP",
                                                   "top_ranks_OP"))) {
  
  # We need the name of the method, name of dilution and seed as coordinates to
  # uniquely identify the tibbles we want to transfer, so we iterate over every 
  # combination of these three markers.
  # Using these three markers, we identify a tibble, and then copy it over in
  # a more convenient hierarchy to the restructure_results list. 
    for(method in script_params$methods_vector) {
      
      # All the dilution names plus OmniPath_0
      for(dilution_name in c("OmniPath_0", 
                             names(script_params$dilution_props))) {
        
        for (seed in names(script_params$master_seed_list)) {
          
          # For tracking purposes we tack the seed name onto the data
          name <- str_glue(dilution_name, "_",seed)
          
          # This top line is the hierarchy we are trying to achieve
          restructured_results[[output]][[method]][[dilution_name]][[name]] <- 
            results[[seed]][[output]][[method]][[dilution_name]]
          # This lower line is the hierarchy as it is per default
          
        }
      }
    }
  }
  
  
  if ("resources_OP" %in% script_params$outputs == TRUE) {
    
    # We need the name of dilution and seed as coordinates to uniquely identify 
    # the tibbles we want to transfer, so we iterate over every combination of 
    # these two markers.
    # Using these two markers, we identify a tibble, and then copy it over in
    # a more convenient hierarchy to the restructure_results list. 
      
    # All the dilution names plus OmniPath_0
    for(dilution_name in c("OmniPath_0", 
                           names(script_params$dilution_props))) {
      
      for (seed in names(script_params$master_seed_list)) {
        
        # For tracking purposes we tack the seed name onto the data
        name <- str_glue(dilution_name, "_",seed)
        
        # This top line is the hierarchy we are trying to achieve
        restructured_results[["resources_OP"]][[dilution_name]][[name]] <- 
          results[[seed]][["resources_OP"]][[dilution_name]]
        # This lower line is the hierarchy as it is per default
        
      }
    }
  }
  
  
  # This code formats top_ranks_analysis outputs, but only if they are actually
  # in results. We need the seed and the top_ranks analysis type to subset this 
  # part of results into the units of data we want to transfer, so we iterate 
  # over all combinations of these two.
  if("top_ranks_analysis" %in% script_params$outputs == TRUE) {
    
    for (analysis in names(results$Seed_1$top_ranks_analysis)) {
      
      for (seed in names(script_params$master_seed_list)) {
        
        # For tracking purposes we tack the seed name onto the data
        name <- str_glue(analysis, "_", seed)
        
        # This top line represents the list hierarchy as we would like it
        restructured_results[["top_ranks_analysis"]][[analysis]][[name]] <- 
          results[[seed]][["top_ranks_analysis"]][[analysis]]
        # and this lower line represents the hierarchy in results.
        
      }
    }
  }
  
  
  
  # Only format metadata if its actually in the results
  if("metadata" %in% script_params$outputs == TRUE) {
    
    # We iterate over every permutation seed
    for (seed in names(script_params$master_seed_list)) {
      
      # Mark the runtime with the iteration it belongs to
      name <- str_glue("Runtime", "_", seed)
      # Grab the runtime from the metadata and save it as its own sublist
      restructured_results[["runtime"]][[name]] <-      # new hierarchy
        results[[seed]][["metadata"]][["runtime"]]      # old hierarchy
      
      
    }
    
  }
  
  
  
  # Only format testdata if its actually in the results
  if("testdata" %in% script_params$outputs == TRUE) {
    
    # We only need the seed to access the testdata tibble of every seed, the 
    # hierarchy is very simple and flat here.
    for (seed in names(script_params$master_seed_list)) {
      
      # new hierarchy
      restructured_results[["testdata"]][[str_glue("Testdata", "_", seed)]] <- 
        results[[seed]][["testdata"]] # old hierarchy
    }
    
  }
  
  
  # We name our restructured results more informatively, then extract the most
  # relevant sublists from them for the rest of the analysis
  dilution_robustness_results <- restructured_results
  
  top_ranks_analysis <- dilution_robustness_results$top_ranks_analysis
  runtime            <- dilution_robustness_results$runtime
  
  
  # Remove unnecessary clutter from the environment.
  rm(analysis, dilution_name, name, method, output, seed, results, 
     restructured_results)
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
    # dilution_Robustness().
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
  {
    # Create a metadata subsection of script_params and assign it runtime
    script_params[["metadata"]][["runtime"]]      <- runtime
    
    # If output was sunk, store the file path to metadata
    if(script_params$sink_output == TRUE) {
      
      script_params$metadata[["sink_logfile"]] <- 
        script_params[["sink_logfile"]]
      
    }
    
    # If warnings were diverted, store the file path to metadata 
    if(script_params$liana_warnings == "divert") {
      
      script_params$metadata[["warning_logfile"]] <- 
        script_params[["warning_logfile"]]
      
    }
    
    # Remove any file paths and session info outside of metadata
    script_params[["sink_logfile"]]    <- NULL
    script_params[["warning_logfile"]] <- NULL
    
    # Remove runtime now that it's a part of script_params$metadata
    rm(runtime)
  }
  
}


#------------------------------------------------------------------------------#
# 5. Aggregate top_ranks_analysis ----------------------------------------------
{
  complete_top_ranks_overlap <- top_ranks_analysis$Overlap %>%
    bind_rows() 
  
  seed_assignment <- sort(rep(1:length(script_params$master_seed_list),
                              length(script_params$dilution_props) +1 ))
  
  complete_top_ranks_overlap <- complete_top_ranks_overlap %>%
    mutate("Seed" = seed_assignment) %>%
    arrange(dilution_prop) %>%
    pivot_longer(cols = script_params$methods_vector, names_to = "Method") %>%
    arrange(Method) %>%
    rename("Overlap" = value) 
  
  rm(seed_assignment)
  

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
    tr_overlap_for_plot <-  complete_top_ranks_overlap  %>%
      as.data.frame()                             %>%
      mutate(dilution_prop = dilution_prop * 100) %>%
      mutate(Overlap       = Overlap       * 100) %>%
      as_tibble()                                 %>%
      drop_na()                                   %>%
      mutate("Method" = recode(Method,
                               "call_connectome" = "Connectome",
                               "call_squidpy"    = "CellPhoneDB",
                               "call_natmi"      = "NATMI",   
                               "call_italk"      = "iTALK", 
                               "call_sca"        = "SingleCellSignalR", 
                               "cellchat"        = "CellChat")) 
     
    # Automatically assemble a file name and plot subtitle
    if (script_params$preserve_topology == FALSE) {
      
      topology_comment <- "random_Dilute()"
      
    } else if (script_params$preserve_topology == TRUE) {
      
      topology_comment <- "preserve_Dilute()"
      
    }
    
    plotting_caption <- 
      str_glue("This plot was created using the ",
               script_params$testdata_type,
               " data. Dilution was performed using ",
               script_params$feature_type,
               " features and the ",
               topology_comment,
               " function. \n",
               "The dilution occured in ",
               (script_params$dilution_props[[2]] -
                  script_params$dilution_props[[1]]) * 100,
               " % increments up to a maximum of ",
               max(tr_overlap_for_plot$dilution_prop),
               " %. \n\n",
               "The overlap was compared between the ",
               as.character(median(unlist(script_params$number_ranks))),
               " highest ranked interactions over ",
               length(script_params$master_seed_list),
               " permutations."
               )
               
    
    if (script_params$run_mode == "trial_run") {
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
  if (script_params$save_results == TRUE) {
    if (script_params$run_mode == "real") {
      test_run_comment <- ""
      
    } else if (script_params$run_mode == "trial_run") {
      test_run_comment <- "TRIAL_RUN_"
      
    }
    
    if (script_params$preserve_topology == FALSE) {
      topology_comment <- "_random_topology_"
      
    } else if (script_params$preserve_topology == TRUE) {
      topology_comment <- "_preserved_topology_"
      
    }
    
    # Define the time of run to uniquely tag every save file
    time_of_run <-  Sys.time() %>%
      as.character()       %>%
      gsub(':', '-', .)    %>% 
      gsub(' ', '_', .)
    
    # Generate the filepaths to save the data under
    box_plot_png_name <-
      str_glue(
        test_run_comment,
        "Boxplot_Resource_Dilution_",
        script_params$testdata_type,
        topology_comment,
        script_params$feature_type,
        "_top",
        as.character(median(unlist(
          script_params$number_ranks
        ))),
        "_",
        time_of_run,
        ".png"
      )
    
    line_plot_png_name <-
      str_glue(
        test_run_comment,
        "Lineplot_Resource_Dilution_",
        script_params$testdata_type,
        topology_comment,
        script_params$feature_type,
        "_top",
        as.character(median(unlist(
          script_params$number_ranks
        ))),
        "_",
        time_of_run,
        ".png"
      )
    
    
    env_save_path <- 
      str_glue(
        "Outputs/",
        test_run_comment,
        "DilutionEnv_",
        script_params$testdata_type,
        topology_comment,
        script_params$feature_type,
        "_top",
        as.character(median(unlist(
          script_params$number_ranks
        ))),
        "_",
        time_of_run,
        ".RData"
      )
    
    env_save_path
    
    
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
    script_params$metadata[["box_plot_png_name"]]  <- box_plot_png_name
    script_params$metadata[["line_plot_png_name"]] <- line_plot_png_name
    script_params$metadata[["env_save_path"]]      <- env_save_path
    
    script_params$metadata[["Session_Info"]]       <- sessionInfo
    
    # Remove clutter
    rm(box_plot_png_name, line_plot_png_name, env_save_path, test_run_comment, 
       topology_comment, time_of_run)
    
    
    # Save R environment and all the results within it
    save.image(file = script_params$metadata$env_save_path)
    
    # Let the user know where everything was stored.
    print(str_glue("Box Plot saved at ~/Outputs/", 
                   script_params$metadata$box_plot_png_name, "."))
    
    print(str_glue("Line Plot saved at ~/Outputs/", 
                   script_params$metadata$line_plot_png_name, "."))
    
    print(str_glue("Environment saved at ~/", 
                   script_params$metadata$env_save_path, "."))
    
  }
}
