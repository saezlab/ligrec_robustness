liana_with_warnings <- function(liana_warnings,
                                methods_vector,
                                cellchat_nperms,
                                testdata,
                                warning_logfile) {
  
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
                 expr_prop = 0,
                 cellchat.params   = list(nboot = cellchat_nperms, 
                                          expr_prop = 0,
                                          thresh = 1),
                 call_natmi.params = list(output_dir = natmi_output))
    
    
  } else if (liana_warnings == "divert") {
    
    divert_Warnings(
      {    
        liana_results <- 
          liana_wrap(testdata, 
                     method = methods_vector, 
                     resource = c('OmniPath'), 
                     expr_prop = 0,
                     cellchat.params   = list(nboot      = cellchat_nperms, 
                                              expr_prop  = 0,
                                              thresh     = 1),
                     call_natmi.params = list(output_dir = natmi_output))
        
      }, logFile = warning_logfile)
    
  } else if (liana_warnings == FALSE) {
    
    suppressWarnings(
      {    
        liana_results <- 
          liana_wrap(testdata, 
                     method = methods_vector, 
                     resource = c('OmniPath'), 
                     expr_prop = 0,
                     cellchat.params   = list(nboot = cellchat_nperms, 
                                              expr_prop = 0,
                                              thresh = 1),
                     call_natmi.params = list(output_dir = natmi_output))
        
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















iterate_liana_wrap <- function(master_seed_list,
                               mismatch_props,
                               reshuffled_clusters,
                               testdata,
                               methods_vector,
                               liana_warnings,
                               warning_logfile,
                               cellchat_nperms) {
  
  print_Title("LIANA with default annotations.", super = TRUE)
  
  runtime <- list("Start Iterations" = Sys.time())
  
  original_results <- liana_with_warnings(testdata        = testdata,
                                          methods_vector  = methods_vector,
                                          liana_warnings  = liana_warnings,
                                          warning_logfile = warning_logfile,
                                          cellchat_nperms = cellchat_nperms)
  
  original_results <- 
    map(master_seed_list, function(seed) {return(original_results)})
  
  runtime[["Default Clusters"]] <- Sys.time()
  
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
        reshuffled_testdata@meta.data[[cluster_col]]
      
      liana_results_mismatch_seed <- 
        liana_with_warnings(liana_warnings  = liana_warnings,
                            methods_vector  = methods_vector,
                            cellchat_nperms = cellchat_nperms,
                            testdata        = reshuffled_testdata,
                            warning_logfile = warning_logfile)
      
      return(liana_results_mismatch_seed)
      
    })
    
    return(liana_results_mismatch)
    
  })
  
  runtime[["Shuffled Clusters"]] <- Sys.time()
  
  complete_liana_results <- original_results %>%
    list("Reshuffle_0" = .)  %>%
    append(., liana_results) %>%
    append(., list("runtime" = runtime))
  
  return(complete_liana_results)
  
}