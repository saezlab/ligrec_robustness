  auto_file_Name <- function(prefix,
                             suffix,
                             
                             trial_run,
                             preserve_topology,
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
    
    
    # How was the topology of diluted resources handled? Make an appropiate 
    # comment.
    if (preserve_topology == FALSE) {
      topology_comment <- "random_topology_"
      
    } else if (preserve_topology == TRUE) {
      topology_comment <- "preserved_topology_"
      
    }
    
    
    # ´What testdata_type was extracted and used with resource_Robustness?
    testdata_comment <-
      str_glue(testdata_type, "_")
    
    
    # What feature_type was dilution specified to occur with?
    feature_type_comment <-
      str_glue(feature_type, "_")
    
    
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
        topology_comment,
        feature_type_comment,
        top_ranks_comment,
        time_of_run,
        suffix
      )
    
    return(auto_file_name)
  }  # end of function
  
}


  
  auto_plot_Description <- function(top_ranks_overlap,
                                    
                                    trial_run,
                                    preserve_topology,
                                    testdata_type,
                                    feature_type,
                                    number_ranks,
                                    time_of_run) {
    
    ## General comment, on testdata type, feature_type and topology
    {
      if (preserve_topology == FALSE) {
        topology_comment <- "random_Dilute()"
        
      } else if (preserve_topology == TRUE) {
        topology_comment <- "preserve_Dilute()"
        
      }
      
      
      general_comment <-
        str_glue(
          "This plot was created using the ",
          testdata_type,
          " data. Dilution was performed using ",
          feature_type,
          " features and the ",
          topology_comment,
          " function. "
        )
      
      rm(topology_comment)
    }
    
    
    
    ## Dilution comment, on proportions
    {
      dilution_overview <- count(top_ranks_overlap,
                                 dilution_prop,
                                 run_mode = "real")
      
      
      dilution_comment <- str_glue(
        "The dilution occured in ",
        dilution_overview$dilution_prop[2] -
          dilution_overview$dilution_prop[1],
        " % increments up to a maximum of ",
        max(top_ranks_overlap$dilution_prop),
        " %. "
      )
      
      if (nrow(dilution_overview) < 1) {
        stop(
          "Expected at least two dilution proportions in input (0, and one ",
          "more. But found only one instead, namely ",
          dilution_overview$dilution_prop
        )
      }
      
      if (length(unique(dilution_overview$n)) != 1) {
        stop(
          "There should be an equal number of samples for every dilution, ",
          "but there is not."
        )
      }
      
      
      rm(dilution_overview)
      
    }
    
    
    ## Nperms and top_ranks comment
    {
      top_ranks_vector <- unlist(number_ranks)
      
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
        stop(
          "There should be an equal number of samples for each method at , ",
          "dilution proportion 0, but there is not."
        )
      }
      
      rm(permutations_overview, top_ranks_vector)
    }
    
    
    ## Date and time comment
    time_comment <- str_glue("Generated at ",
                             time_of_run,
                             ".")
    
    
    ## Assemple plotting caption
    plotting_caption <-
      str_glue(
        general_comment,
        "\n",
        dilution_comment,
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
