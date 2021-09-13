summarise_Metadata <- function(runtime,
                               time_of_run,
                               dilution_params,
                               master_seed_list,
                               save_results = TRUE,
                               trial_run    = TRUE) {
  
  meta_params <- list("time_of_run"  = time_of_run,
                      "save_results" = save_results,
                      "trial_run"    = trial_run)
  
  
  metadata <- list(
    "runtime"         = runtime,
    "dilution_params" = append(dilution_params, 
                               master_seed_list),
    "meta_params"     = meta_params,
    "sessionInfo"     = sessionInfo()
  )
  
  
  if (save_results == TRUE) {
    # Generate the filepaths to save the data under
    metadata[["box_plot_png_name"]] <-
      auto_file_Name(
        prefix = "Boxplot_Resource_Dilution_",
        suffix = ".png",
        dilution_params = formals(wrap_resource_Robustness),
        meta_params     = formals(summarise_Metadata),
        time_of_run     = time_of_run
      )
    
    metadata[["line_plot_png_name"]] <-
      auto_file_Name(
        prefix = "Lineplot_Resource_Dilution_",
        suffix = ".png",
        dilution_params = formals(wrap_resource_Robustness),
        meta_params     = formals(summarise_Metadata),
        time_of_run     = time_of_run
      )
    
    metadata[["env_save_path"]] <-
      auto_file_Name(
        prefix = "Outputs/DilutionEnv_",
        suffix = ".RData",
        dilution_params = formals(wrap_resource_Robustness),
        meta_params     = formals(summarise_Metadata),
        time_of_run     = time_of_run
      )
  }

  
  if(dilution_params$sink_output == TRUE) {
    # The file name includes many script_params in it to be informative and
    # unique.
    metadata[["sink_logfile"]] <- auto_file_Name(prefix = "Outputs/Logs/Complete_Log_",
                                   suffix =  ".txt",
                                   dilution_params = formals(wrap_resource_Robustness),
                                   meta_params = formals(summarise_Metadata),
                                   time_of_run = time_of_run)
    
    
  }
  
  # If necessary, create a log name, store it in script params, then remove
  # the leftover clutter
  if(dilution_params$liana_warnings == "divert") {
    # The file name includes many script_params in it to be informative and
    # unique.
    metadata[["warning_logfile"]] <- auto_file_Name(prefix = "Outputs/Logs/LIANA_warnings_",
                                      suffix =  ".txt",
                                      dilution_params = formals(wrap_resource_Robustness),
                                      meta_params = formals(summarise_Metadata),
                                      time_of_run = time_of_run)
    
    
  }
  
  
  return(metadata)
  
}



# wrap_resource_Robustness()
{
  # Which scRNA data set do you want to use?
  
  
  # Should dilution be with any genes from testdata or most variable ones?
  
  
  # Preserve topology after dilution or not?
  # TRUE = preserve_Dilute(), FALSE = random_Dilute()
  
  
  # Which proportions should the resources that are analysed have?
  # Must be at least one integer.
  
  
  # Which Outputs from resource_Robustness() do you want? Choose from:
  
  
  # Which methods should resource_Robustness() use? Choose from:
  # Squidpy doesn't work on windows.
  
  # Which top n of interactions should be considered top-ranked per method?
  
  
  # number of cellchat permutations, default 100
  
  # Should the entire output be redirected to a log file?
  
  # Should liana warnings be visible in output? 
  
}

wrap_resource_Robustness <- 
  function(master_seed,
           
           testdata_type     = "liana_test", # "liana_test" or "seurat_pbmc", only as a string, not a vector of a string
           feature_type      = "variable",   # "generic" or "variable", only as a string, not a vector of a string
           preserve_topology = FALSE,           # TRUE or FALSE
           dilution_props    = c(seq(0.40, 1.00, 0.40)),
           
           outputs = c(
             "liana_results_OP"  ,
             "resources_OP"      ,
             "top_ranks_OP"      ,
             "top_ranks_analysis",
             "runtime"          ,
             "testdata"
           ),
           
           methods_vector = c(
             'call_connectome' ,
             #'squidpy'         ,
             'call_natmi'      ,
             'call_italk'      ,
             'call_sca'        ,
             'cellchat'
           ), 
           
           number_ranks = list(
             "call_connectome" = 20,
             "squidpy"         = 20,
             "call_natmi"      = 20,
             "call_italk"      = 20,
             "call_sca"        = 20,
             "cellchat"        = 20
           ),
           
           cellchat_nperms = 10,      # default 100 for real data
           sink_output     = FALSE,   # TRUE or FALSE
           liana_warnings  = "divert", # TRUE, FALSE, or "divert"
           
           time_of_run
           
  ) {
    
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
    
    ## Define File Paths for Logs
    {
      sink_logfile <- ""
      warning_logfile <- ""
      
      
      # If they are required, we auto generate log filepaths here. The filepaths
      # have the current time imminetnly before the lapply in them. Each lapply
      # should be assigned to a unique log this way that has all iterations in
      # it. 
      

      # If necessary, create a log name, store it in script params, then remove
      # the leftover clutter
      if(sink_output == TRUE) {
        # The file name includes many script_params in it to be informative and
        # unique.
        sink_logfile <- auto_file_Name(prefix = "Outputs/Logs/Complete_Log_",
                                       suffix =  ".txt",
                                       dilution_params = formals(wrap_resource_Robustness),
                                       meta_params = formals(summarise_Metadata),
                                       time_of_run = time_of_run)
        
        
      }
      
      # If necessary, create a log name, store it in script params, then remove
      # the leftover clutter
      if(liana_warnings == "divert") {
        # The file name includes many script_params in it to be informative and
        # unique.
        warning_logfile <- auto_file_Name(prefix = "Outputs/Logs/LIANA_warnings_",
                                          suffix =  ".txt",
                                          dilution_params = formals(wrap_resource_Robustness),
                                          meta_params = formals(summarise_Metadata),
                                          time_of_run = time_of_run)
          
        
      }
    }
    
    
    
    resource_Robustness(
      master_seed       = master_seed, 
      
      testdata_type     = testdata_type,
      feature_type      = feature_type,
      preserve_topology = preserve_topology,
      dilution_props    = dilution_props,
      number_ranks      = number_ranks,
      outputs           = outputs,
      
      methods_vector    = methods_vector,
      cellchat_nperms   = cellchat_nperms,
      sink_output       = sink_output,
      liana_warnings    = liana_warnings,
      
      sink_logfile      = sink_logile,
      warning_logfile   = warning_logfile
    )
    
    
    
  }




auto_file_Name <- function(prefix, 
                           suffix, 
                           dilution_params,
                           meta_params,
                           time_of_run) {
  
  if (meta_params$trial_run == FALSE) {
    test_run_comment <- ""
    
  } else if (meta_params$trial_run == TRUE) {
    test_run_comment <- "TRIAL_RUN_"
    
  }
  
  
  
  if (dilution_params$preserve_topology == FALSE) {
    topology_comment <- "random_topology_"
    
  } else if (dilution_params$preserve_topology == TRUE) {
    topology_comment <- "preserved_topology_"
    
  }
  
  
  
  testdata_comment <-     
    str_glue(dilution_params$testdata_type, "_")
  
  
  
  feature_type_comment <-
    str_glue(dilution_params$feature_type, "_")
  
  
  
  top_ranks_vector <- 
    unlist(as.list(dilution_params$number_ranks)[-1])
  
  top_ranks_comment <- 
    str_glue("top", median(top_ranks_vector), "_",)
  
  
  
  auto_file_path <- 
    str_glue(
      prefix,
      test_run_comment,
      testdata_comment,
      topology_comment,
      feature_type_comment,
      top_ranks_comment, 
      time_of_run,
      suffix
    )
  
  return(auto_file_path)
}


auto_plot_Description <- function(top_ranks_overlap, 
                                  dilution_params, 
                                  meta_params,
                                  time_of_run) {
  
  ## General comment, on testdata type, feature_type and topology 
  {
    if (dilution_params$preserve_topology == FALSE) {
      
      topology_comment <- "random_Dilute()"
      
    } else if (dilution_params$preserve_topology == TRUE) {
      
      topology_comment <- "preserve_Dilute()"
      
    }
    
    
    general_comment <- str_glue("This plot was created using the ",
                                dilution_params$testdata_type,
                                " data. Dilution was performed using ",
                                dilution_params$feature_type,
                                " features and the ",
                                topology_comment,
                                " function. ")
    
    rm(topology_comment)
  }
  
  
  
  ## Dilution comment, on proportions 
  {
    dilution_overview <- count(top_ranks_overlap, 
                               dilution_prop, 
                               run_mode = "real")
    
    
    dilution_comment <- str_glue("The dilution occured in ", 
                                 dilution_overview$dilution_prop[2] -
                                   dilution_overview$dilution_prop[1],
                                 " % increments up to a maximum of ",
                                 max(top_ranks_overlap$dilution_prop),
                                 " %. ")
    
    if (nrow(dilution_overview) < 1) {
      stop("Expected at least two dilution proportions in input (0, and one ",
           "more. But found only one instead, namely ",
           dilution_overview$dilution_prop)
    }
    
    if (length(unique(dilution_overview$n)) != 1) {
      stop("There should be an equal number of samples for every dilution, ",
           "but there is not.")
    }
    
    
    rm(dilution_overview)
    
  }
  
  
  ## Nperms and top_ranks comment
  {
    top_ranks_vector <- 
      unlist(as.list(dilution_params$number_ranks)[-1])
    
    permutations_overview <- top_ranks_overlap %>%
      filter(dilution_prop == 0) %>%
      count(Method)
    
    
    top_ranks_permutations_comment <-
      str_glue(
        "The overlap was compared between the ",
        median(top_ranks_vector),
        " highest ranked interactions over ",
        permutations_overview$n[1],
        " permutations."
      )
    
    if (length(unique(permutations_overview$n)) != 1) {
      stop("There should be an equal number of samples for each method at , ",
           "dilution proportion 0, but there is not.")
    }
    
    rm(permutations_overview, top_ranks_vector)
  }
  
  
  ## Date and time comment
  time_comment <- str_glue("Generated at ",
                           time_of_run,
                           ".")
  
  
  ## Assemple plotting caption
  plotting_caption <- 
    str_glue(general_comment,
             "\n",
             dilution_comment,
             "\n\n",
             top_ranks_permutations_comment,
             "\n",
             time_comment
    )
  
  
  ## Add addendum if trial run
  if (meta_params$trial_run == TRUE) {
    plotting_caption <- 
      str_glue(plotting_caption, "   --   [TRIAL RUN]")
  }
  
  return(plotting_caption)
}






save_Results <- function(dilution_params,
                         meta_params, 
                         time_of_run) {
  
  if (meta_params$trial_run == FALSE) {
    test_run_comment <- ""
    
  } else if (meta_params$trial_run == TRUE) {
    test_run_comment <- "TRIAL_RUN_"
    
  }
  
  if (dilution_params$preserve_topology == FALSE) {
    topology_comment <- "_random_topology_"
    
  } else if (dilution_params$preserve_topology == TRUE) {
    topology_comment <- "_preserved_topology_"
    
  }
  
  top_ranks_vector <- 
    unlist(as.list(dilution_params$number_ranks)[-1])
  
  # Generate the filepaths to save the data under
  box_plot_png_name <-
    auto_file_Name(
      prefix = "Boxplot_Resource_Dilution_",
      suffix = ".png",
      dilution_params = dilution_params,
      meta_params     = meta_params,
      time_of_run     = time_of_run
    )
  
  line_plot_png_name <-
    auto_file_Name(
      prefix = "Lineplot_Resource_Dilution_",
      suffix = ".png",
      dilution_params = dilution_params,
      meta_params     = meta_params,
      time_of_run     = time_of_run
    )
  
  
  env_save_path <- auto_file_Name(
    prefix = "Outputs/DilutionEnv_",
    suffix = ".RData",
    dilution_params = dilution_params,
    meta_params     = meta_params,
    time_of_run     = time_of_run
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
