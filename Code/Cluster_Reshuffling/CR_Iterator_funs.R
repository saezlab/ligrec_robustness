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


wrap_Shuffler <- function(master_seed_list,
                          mismatch_prop,
                          metadata,
                          cluster_key) {
  
  print_Title(str_glue("Creating ",
                       mismatch_prop*100, 
                       " % mismatched cluster annotations."))
  
  reshuffled_metadatas <- lapply(master_seed_list,
                                 shuffle_Clusters,
                                 mismatch_prop = mismatch_prop,
                                 metadata      = metadata)
  
  
  return(reshuffled_metadatas)
  
}


shuffle_Clusters <- function(master_seed,
                             mismatch_prop,
                             metadata) {
  
  set.seed(master_seed)
  
  
  metadata_old <- metadata %>%
    rownames_to_column(var = "Bar_Code") %>%
    as_tibble()
  
  metadata_old$cluster_key <- metadata_old$cluster_key %>%
    as.numeric()
  
  meta_clusters_dilute <- slice_sample(metadata_old, prop = mismatch_prop) %>%
    select(c("Bar_Code", all_of("cluster_key")))
  
  
  
  
  meta_clusters_dilute$cluster_key <- 
    map(meta_clusters_dilute$cluster_key, function(annotation) {
      
      replacement_candidates <-
        filter(metadata_old, 
               metadata_old$cluster_key != annotation)$cluster_key
      
      new_annotation <- sample(replacement_candidates, 1)
      
      return(new_annotation)
      
    }) %>%
    unlist() 
  
  
  
  metadata_new <-
    left_join(metadata_old, 
              meta_clusters_dilute, 
              by = "Bar_Code",
              suffix = c("_old", ""))
  
  helper_index <- which(is.na(metadata_new$cluster_key))
  
  metadata_new$cluster_key[helper_index] <-
    metadata_new[["cluster_key_old"]][helper_index]
  
  metadata_new <- metadata_new %>%
    mutate("Mismatched" = 
             ifelse(.$cluster_key_old == .$cluster_key, FALSE, TRUE)) %>%
    select(-"cluster_key_old") %>%
    as.data.frame()
  
  rownames(metadata_new) <- metadata_new$Bar_Code
  
  metadata_new <- metadata_new %>%
    select(-"Bar_Code")
  
  
  
  metadata_new$cluster_key <- metadata_new$cluster_key %>%
    as.factor()
  
  print(
    str_glue("Cluster annotations reshuffled. ", 
             round((sum(metadata_new$Mismatched) / nrow(metadata_new)) *100, 2),
             " % mismatch to the original annotation."))
  
  
  
  
  return(metadata_new)
}


iterate_liana_wrap <- function(master_seed_list,
                               mismatch_props,
                               reshuffled_clusters,
                               testdata,
                               methods_vector,
                               liana_warnings,
                               warning_logfile,
                               ...) {
  
  
  
  runtime <- list("Start Iterations" = Sys.time())
  
  liana_results <- map(mismatch_props, function(mismatch_prop) {
    
    print_Title(str_glue("LIANA for ", 
                         mismatch_prop*100, 
                         " % mismatched cluster annotations"), 
                super = TRUE)
    
    liana_results_mismatch <- map(master_seed_list, function(seed) {
      
      print_Title(str_glue(mismatch_prop*100,
                           " % Iteration ",
                           seed,
                           " ",
                           as.character(Sys.time())))
      
      reshuffled_testdata <- testdata
      mismatch_name <- names(which(mismatch_props == mismatch_prop))
      
      reshuffled_testdata@meta.data <- 
        reshuffled_clusters[[mismatch_name]][[seed]]
      
      Idents(reshuffled_testdata) <- 
        reshuffled_testdata@meta.data$cluster_key
      
      liana_results_mismatch_seed <- 
        liana_with_warnings(liana_warnings  = liana_warnings,
                            methods_vector  = methods_vector,
                            testdata        = reshuffled_testdata,
                            warning_logfile = warning_logfile,
                            ...)
      
      return(liana_results_mismatch_seed)
      
    })
    
    return(liana_results_mismatch)
    
  })
  
  runtime[["Shuffled Clusters"]] <- Sys.time()
  
  print_Title("LIANA with default annotations.", super = TRUE)
  
  original_results <- liana_with_warnings(testdata        = testdata,
                                          methods_vector  = methods_vector,
                                          liana_warnings  = liana_warnings,
                                          warning_logfile = warning_logfile,
                                          ...)
  
  original_results <- 
    map(master_seed_list, function(seed) {return(original_results)})
  
  runtime[["Default Clusters"]] <- Sys.time()
  
  
  
  complete_liana_results <- original_results %>%
    list("Reshuffle_0" = .)  %>%
    append(., liana_results) %>%
    append(., list("runtime" = runtime))
  
  return(complete_liana_results)
  
}


liana_with_warnings <- function(liana_warnings,
                                methods_vector,
                                testdata,
                                warning_logfile,
                                ...) {
  
  
  # Generate Undiluted liana results by running wrapper function
  # Omnipath x the methods vector, on the selected data
  
  # NATMI results are contaminated with results from earlier runs if you
  # don't specify a special output folder for the results to go in
  natmi_output <-  Sys.time() %>%
    as.character()       %>%
    gsub(':', '-', .)    %>% 
    gsub(' ', '_', .)    %>%
    str_glue('Test_', .)
  
  # The if statements give the user control over how warnings are handled
  if (liana_warnings == TRUE) {
    
    liana_results <- 
      liana_wrap(testdata, 
                 method = methods_vector, 
                 resource = c('OmniPath'), 
                 call_natmi.params = list(output_dir = natmi_output),
                 ...)
    
    
  } else if (liana_warnings == "divert") {
    
    divert_Warnings(
      {    
        liana_results <- 
          liana_wrap(testdata, 
                     method = methods_vector, 
                     resource = c('OmniPath'), 
                     call_natmi.params = list(output_dir = natmi_output),
                     ...)
        
      }, logFile = warning_logfile)
    
  } else if (liana_warnings == FALSE) {
    
    suppressWarnings(
      {    
        liana_results <- 
          liana_wrap(testdata, 
                     method = methods_vector, 
                     resource = c('OmniPath'), 
                     call_natmi.params = list(output_dir = natmi_output),
                     ...)
        
      })
    
  }
  
  # If only one method is fed to liana_wrap(), the output data structure is
  # different. So here we convert it to the same data structure that would 
  # exist if multiple methods had been called, so the code below works
  # properly for single method runs too.
  if (length(methods_vector) == 1) {
    
    liana_results <- list(liana_results)
    
    names(liana_results) <- methods_vector
    
  }
  
  return(liana_results)
  
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


# get_top_ranks_clust()
{
  #' Get the top n ranked items of a method from the tibble liana wrapper or
  #' call_x results
  #'
  #' @param data_set The tibble output by the liana wrapper function or call_x
  #' functions.
  #' @param method The method for which you would like to extract the top ranked
  #' interactions, as a string.
  #' @param top_n The number of items to return, returns items ranked #1 to #n.
  #' @param with_ties TRUE or FALSE. An argument passed to slice_min and 
  #' slice_max respectively. When slicing the top_n, sometimes a cutoff occurs
  #' between equally ranked items. Should the items also be included, expanding 
  #' the selection (TRUE), or should preserving top_n be prioritized (FALSE). 
  #'
  #' @return Returns the tibble input cut down to the top n highest ranked
  #' interactions.
  
  get_top_ranks_clust <- function(data_set, method, top_n, with_ties = FALSE) {
    # generate a list that describes how to rank in each method to get what the
    # method considers best
    rank_spec_list <-
      list(
        "cellchat"        = list(method_score = "pval",
                                 descending_order =  FALSE),
        
        "call_connectome" = list(method_score = "weight_sc",
                                 
                                 descending_order =  TRUE),
        
        "call_italk"      = list(method_score = "logfc_comb",
                                 descending_order =  TRUE),
        
        "call_natmi"      = list(method_score = "edge_specificity",
                                 descending_order =  TRUE),
        
        "call_sca"        = list(method_score = "LRscore",
                                 descending_order =  TRUE),
        
        "squidpy"         = list(method_score = "pvalue",
                                 descending_order =  FALSE)
      )
    
    
    
    # If its a p-value method, the code will go into this if statement
    if (rank_spec_list[[method]]$descending_order == FALSE) {
      topranks <-
        slice_min(
          data_set,
          n = top_n,
          order_by = !!sym(rank_spec_list[[method]]$method_score),
          with_ties = with_ties
        )
      
      # order_by is assigned the criterion dictated by rank_spec_list
      
      
    } else {
      # if it's not one of the p-value methods
      
      topranks <-
        slice_max(
          data_set,
          n = top_n,
          order_by = !!sym(rank_spec_list[[method]]$method_score),
          with_ties = with_ties
        )
      
      # order_by is assigned the criterion dictated by rank_spec_list
      
      
    }
    
    return(topranks)
    
  } #end of function
  
  
  
}


format_top_ranks <- function(top_ranks) {
  
  # Format the top_ranks data frame for future processing steps.
  top_ranks <- top_ranks %>%
    unite("LR_Pair",
          c(ligand, receptor),
          remove = FALSE,
          sep = "_") %>%
    relocate("LR_Pair", .after = last_col()) %>%
    unite("LR_ID",
          c(source, target, ligand, receptor),
          remove = FALSE,
          sep = "_") %>%
    relocate("LR_ID", .after = last_col())
  
  
  return(top_ranks)
  
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


