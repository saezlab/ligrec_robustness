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
    dilution_props <- c(seq(0.20, 0.80, 0.20))
    
    # How many permutations of dilution should be performed?
    master_seed_list <- as.list(c(1:5))
    
    # Which Outputs from dilution_Robustness() do you want? Choose from:
    # outputs = c("liana_results_OP", "resources_OP", "top_ranks_OP", 
    #             "top_ranks_analysis","metadata", "testdata")
    outputs = c("liana_results_OP", "resources_OP", "top_ranks_OP",
                "top_ranks_analysis","metadata")
    
    
    # Which methods should dilution_Robustness() use? Choose from:
    # methods_vector <- c('call_connectome', 'call_natmi', 'call_italk', 
    #                     'call_sca', 'cellchat')
    methods_vector <- c('call_connectome', 'call_italk', 
                        'call_sca')
    # no squidpy until it works on windows.
    
    # Which top n of interactions should be considered top-ranked per method?
    number_ranks   <- list("call_connectome" = 20, 
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
      time_of_run <- Sys.time()     %>%
        as.character()       %>%
        gsub(':', '_', .)    %>% 
        gsub('-', '_', .)    %>% 
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
  
  # The three outputs mentioned here can all be formatted the same way
  # But only execute this code if that output is actually in results
  for(output in intersect(script_params$outputs, c("liana_results_OP",
                                                   "resources_OP",
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
  
  for (row in 1:(length(dilution_props)+1)) {
    
    mean_overlap <- mean(as.numeric(
      top_ranks_overlaps[row, grepl(method, names(top_ranks_overlaps))]))
    
    sd_overlap   <- sd(as.numeric(
      top_ranks_overlaps[row, grepl(method, names(top_ranks_overlaps))]))
    
    vector_mean <- c(vector_mean, mean_overlap)
    vector_sd   <- c(vector_sd,   sd_overlap)
    
  }
  
  top_ranks_aggreg_mean[[str_glue(method, "_mean")]] <- vector_mean
  top_ranks_aggreg_sd[[str_glue(method, "_sd")]]     <- vector_sd
}

agg_top_ranks_overlap_mean <- as_tibble(top_ranks_aggreg_mean) 

agg_top_ranks_overlap_sd   <- as_tibble(top_ranks_aggreg_sd) 

agg_top_ranks <- tibble(agg_top_ranks_overlap_mean,
                        agg_top_ranks_overlap_sd,
                        .name_repair = c("universal"))

agg_top_ranks <- agg_top_ranks %>%
  mutate(dilution_prop = c(0, dilution_props))  %>%
  unnest(cols = c(dilution_prop))               %>%
  relocate("dilution_prop")


#------------------------------------------------------------------------------#
# 4. Visualizing the results -------------------------------------------------
{ 

  
  # 4.1 Plotting, labeling and saving top_ranks_overlap 
  {
    
    
    # The plot is better in percent than proportion
    tr_overlap_for_plot <-  agg_top_ranks * 100
    
    # Automatically assemble a file name and plot subtitle
    if (preserve_topology == FALSE) {
      
      topology_comment <- "using random_Dilute()"
      
    } else if (preserve_topology == TRUE) {
      
      topology_comment <- "using preserve_Dilute()"
      
    }
    
    plotting_subtitle <- str_glue(feature_type,
                                  " dilution, top ",
                                  as.character(median(unlist(number_ranks))),
                                  " ranks, ",
                                  testdata_type,
                                  " data, ",
                                  run_mode,
                                  " results, ",
                                  topology_comment,
                                  ", ",
                                  length(master_seed_list),
                                  " permutations")
    
    plot_png_name     <- str_glue(run_mode,
                                  "_",
                                  testdata_type, 
                                  "_top",
                                  as.character(median(unlist(number_ranks))),
                                  "_res",
                                  as.character(length(dilution_props)),
                                  "_",
                                  feature_type,
                                  "_dil_on_",
                                  as.character(Sys.Date()),
                                  ".png")
    
    # Plot top_ranks_overlap with lines and points at each value
    overlap_plot <-  ggplot(data = tr_overlap_for_plot) + 
      geom_line(mapping = aes(dilution_prop, 
                              call_connectome_mean, 
                              color =  "Connectome")) +
      
      geom_line(mapping = aes(dilution_prop, 
                              call_natmi_mean, 
                              color = "NATMI")) + 
      
      geom_line(mapping = aes(dilution_prop,
                              call_italk_mean, 
                              color = "iTALK")) +
      
      geom_line(mapping = aes(dilution_prop, 
                              call_sca_mean, 
                              color = "SCA")) +
      
      geom_line(mapping = aes(dilution_prop, 
                              cellchat_mean, 
                              color = "CellChat")) +
      
      
      
      geom_point(mapping = aes(dilution_prop, 
                               call_connectome_mean, 
                               color =  "Connectome")) +
      
      geom_point(mapping = aes(dilution_prop,
                               call_natmi_mean, 
                               color = "NATMI")) +
      
      geom_point(mapping = aes(dilution_prop, 
                               call_italk_mean,
                               color = "iTALK")) +
      
      geom_point(mapping = aes(dilution_prop, 
                               call_sca_mean, 
                               color = "SCA")) +
      
      geom_point(mapping = aes(dilution_prop, 
                               cellchat_mean, 
                               color = "CellChat")) +
      

      geom_errorbar(aes(ymin = call_connectome_mean - call_connectome_sd,
                        ymax = call_connectome_mean + call_connectome_sd,
                        x = dilution_prop), 
                    width=.2,
                    position=position_dodge(0.05)) + 
      
      geom_errorbar(aes(ymin = call_natmi_mean - call_natmi_sd,
                        ymax = call_natmi_mean + call_natmi_sd,
                        x = dilution_prop), 
                    width=.2,
                    position=position_dodge(0.05)) + 
      
      geom_errorbar(aes(ymin = call_italk_mean - call_italk_sd,
                        ymax = call_italk_mean + call_italk_sd,
                        x = dilution_prop), 
                    width=.2,
                    position=position_dodge(0.05)) + 
      
      geom_errorbar(aes(ymin = call_sca_mean - call_sca_sd,
                        ymax = call_sca_mean + call_sca_sd,
                        x = dilution_prop), 
                    width=.2,
                    position=position_dodge(0.05)) + 
      
      geom_errorbar(aes(ymin = cellchat_mean - cellchat_sd,
                        ymax = cellchat_mean + cellchat_sd,
                        x = dilution_prop), 
                    width=.2,
                    position=position_dodge(0.05)) + 
      
      
      
      # Show full breadth of 100-0 percent overlap
      ylim(0, 100) +
      
      ggtitle("Robustness of Method Predictions") +
      ylab("Overlap of Top Ranks [%]") +
      xlab("Dilution of Resource [%]") +
      labs(subtitle = plotting_subtitle,
           color = "Method")
    
    
    
    # Print the Plot
    print(overlap_plot)
    
    # Save the plot automatically to the outputs folder, if desired
    if (save_results) {
      
      ggsave(plot_png_name, 
             height = 5, width = 8, 
             path = "Outputs")
      
      print(str_glue("Plot saved at ~/Outputs/", plot_png_name, "."))
      
    }
    
    # Remove unnecessary variables
    rm(tr_overlap_for_plot, overlap_plot, plotting_subtitle, topology_comment)
    
    
    
    
  } # end of subpoint
  
 }

