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
