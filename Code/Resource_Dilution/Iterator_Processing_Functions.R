# wrap_resource_Robustness()
{
  #' Wrapper function for resource_Robustnes()
  #' 
  #' @description This function feeds its arguments to resource_Robustness. In 
  #' addition, the dilution_proportions are converted from a simple vector to a 
  #' named list, and filepaths for logs are automatically generated if needed. 
  #' 
  #' @param master_seed  At some stages when diluting  a resource, randomness
  #' is at play. Setting a master seed will ensure that dilute_Resource(), the 
  #' function that dilutes the resources, runs the same way every time. Each 
  #' instance of randomness within the scope of dilute_Resource() will run 
  #' reproducibly.
  #' 
  #' @param testdata A preprocessed seurat object of test data that LIANA++ can 
  #' run on.
  #' 
  #' @param feature_type Should dilution occur with all the genes profiled in
  #' testdata (choose "generic") or with the most variable features (choose 
  #' "variable")? Only takes a string as default, not a vector of a string.
  #' 
  #' @param preserve_topology When diluting, two methods are implemented.
  #' random_Dilute makes only a small effort to preserve the topology of the 
  #' rows that are diluted from the resource, preserve_Dilute makes a far 
  #' greater effort. Choose TRUE for preserve_Dilute and FALSE for 
  #' random_Dilute. Only takes a boolean as default, not a vector of a boolean.
  #' 
  #' @param dilution_props A sequence of numerics (0-1) that indicate which 
  #' proportions to dilute the resource with. For example c(0.1, 0.2, 0.3) would
  #' compare the top ranked CCI's using undiluted OmniPath compared to OmniPath 
  #' with 10 % of its rows diluted, undiluted vs 20 % diluted, and undiluted vs
  #' 30 % diluted.
  #' 
  #' @param bundled_outputs Which outputs of the calculation would you like to 
  #' return? Outputs are returned bundled up in a list. By default, only the 
  #' analysis of top_ranks and runtime data is returned, but more information 
  #' can be returned. Construct an atomic vector using all or some of 
  #' "liana_results_OP", "resources_OP", "top_ranks_OP", "top_ranks_analysis",
  #' "runtime", and "testdata" to tell the script which outputs should go in the
  #' output bundle.
  #' 
  #' "liana_results_OP": All the predicted CCI's from LIANA++ for each method 
  #' and dilution stage.
  #' 
  #' "resources_OP": All the resources that were used are returned.
  #' 
  #' "top_ranks_OP": All the highest ranked CCI's from LIANA++ for each method 
  #' and dilution stage.
  #' 
  #' "top_ranks_analysis": How the diluted top predictions related to the 
  #' undiluted predictions in terms of overlap, mismatch, proportion of fake 
  #' interactions in top_ranks, and proportion of mismatch caused by fake 
  #' interactions. 
  #' 
  #' "runtime": What was the runtime of this function was like.
  #' 
  #' "testdata": What was the testdata that was used in this run. Since that's 
  #' also a user supplied Seurat, there is almost never a reason to return this.
  #' 
  #' @param methods_vector Which methods should the function run? Choose from
  #' "call_connectome", "squidpy", "call_natmi", "call_italk", "call_sca" and
  #' "cellchat". Supply the argument in the form of e.g. 
  #' c("call_conncectome, "call_italk"). Only takes a vector of chars as 
  #' default.
  #' 
  #' @param number_ranks A named list. Each item is named after a method and is 
  #' equal to the number of top interactions considered relevant for that 
  #' method. Only takes a list named after methods and containing numerics as
  #' default.
  #'  
  #' @param cellchat_nperms Cellchat is one of the slower methods, for test runs
  #' it may be useful to set this parameter to 10 to speed up the analysis.
  #'
  #' @param sink_otuput TRUE or FALSE. Should the function save a full log of 
  #' the Console Output as a log to the log folder? Warnings and messages will 
  #' not be visible in the console output if this is enabled, so unless there is
  #' a reason why such a record is necessary this option is not recommended.
  #' 
  #' @param liana_warnings Either TRUE, FALSE or "divert". Should the warnings 
  #' from liana, which are often repetitive and unhelpful, be either suppressed
  #' or alternatively diverted to the log folder? When these types of warning 
  #' are left in, they can often displace valuable warnings. Be careful with 
  #' this setting, as suppressing warnings is obviously risky.
  #' 
  #' @param time_of_run What time was this script run at? Used as an argument 
  #' for auto_file_Name() to uniquely tag logs so that thee user can associate 
  #' the logs with the correct plots and data, and so that a new log never 
  #' overwrites an old one.
  #' 
  #' @return A list that contains all the data specified in outputs. At minimum,
  #' the top_ranks_analysis list will be contained, which contains various 
  #' comparisons of the top_ranks predicted at dilution_prop = 0 and the 
  #' various higher dilution proportions.
  
  wrap_resource_Robustness <-
    function(master_seed,
             testdata,
             feature_type,
             preserve_topology, 
             dilution_props,
             number_ranks,
             methods_vector,
             
             bundled_outputs,
             cellchat_nperms,
             
             sink_output,     
             liana_warnings,
             trial_run,
             time_of_run,
             testdata_type) {
      
      
      ## Process Dilution Proportions List
      {
        # By naming each dilution proportion we can use this to label data 
        # easily. By formatting proportions as a list we can lapply over them 
        # conveniently
        
        # Convert dilution.props to a list and name it, creating a named list
        dilution_props <- as.list(dilution_props)
        
        names(dilution_props) <- map(dilution_props, function(prop) {
          
          # Name every dilution proportion
          str_glue("OmniPath_", as.character(prop * 100))
          
        }) 
      }
      
      
      ## Define File Paths for Logs
      {
        # If they are required, we auto generate log filepaths here. The 
        # filepaths have the current time imminently before the lapply in them. 
        # Each run of the iterator should be assigned to a unique log this way 
        # that has all iterations in it.
        
        # Initiate empty filepaths
        sink_logfile <- ""
        warning_logfile <- ""
        
        
        # If necessary, create a log name, store it in script params, then 
        # removethe leftover clutter
        if (sink_output == TRUE) {
          # The file name includes many script_params in it to be informative and
          # unique. RD stands for Resource Dilution.
          sink_logfile <-
            auto_file_Name(
              prefix = "Outputs/Resource_Dilution/Logs/Complete_Log_RD_",
              suffix =  ".txt",
              
              preserve_topology  = preserve_topology,
              testdata_type      = testdata_type,
              feature_type       = feature_type,
              number_ranks       = number_ranks,
              time_of_run        = time_of_run,
              trial_run          = trial_run)
          
          
        }
        
        # If necessary, create a log name, store it in script params, then 
        # remove the leftover clutter
        if (liana_warnings == "divert") {
          # The file name includes many script_params in it to be informative and
          # unique. RD stands for Resource Dilution.
          warning_logfile <-
            auto_file_Name(
              prefix = "Outputs/Resource_Dilution/Logs/LIANA_warnings_RD_",
              suffix =  ".txt",
              
              preserve_topology  = preserve_topology,
              testdata_type      = testdata_type,
              feature_type       = feature_type,
              number_ranks       = number_ranks,
              time_of_run        = time_of_run,
              trial_run          = trial_run)
          
          
        }
      }
      
      
      ## Run resource robustness with the wrapper defaults and new arguments.
      resource_Robustness(
        master_seed       = master_seed,
        testdata          = testdata ,
        
        feature_type      = feature_type,
        preserve_topology = preserve_topology,
        dilution_props    = dilution_props,
        number_ranks      = number_ranks,
        bundled_outputs   = bundled_outputs,
        
        methods_vector    = methods_vector,
        cellchat_nperms   = cellchat_nperms,
        sink_output       = sink_output,
        liana_warnings    = liana_warnings,
        
        sink_logfile      = sink_logile,
        warning_logfile   = warning_logfile
      )
      
      
      
    } # end of function
  
  
}

# extract_Testdata()
{
  #' Helper function that gets a specific seurat object from the outputs folder
  #' 
  #' @param testdata_type As a string. Which testdata should be retrieved? 
  #' Either "seurat_pbmc" or "liana_test". Seurat_pbmc is the data set used in 
  #' the seurat tutorial, while liana_test is the testdata that comes with 
  #' LIANA++, and is a small subset of seurat_pbmc. Only takes a string as 
  #' default, not a vector of a string.
  #' 
  #' @return A seurat object loaded from the outputs folder or liana package.
  
  
  extract_Testdata <- function(testdata_type) {
    
    # Get seurat or liana test data
    if (testdata_type == "seurat_pbmc") {
      
      # Read testdata from outputs
      testdata <- readRDS(file = "Data/pbmc3k_final.rds")     
      
    } else if (testdata_type == "liana_test") {
      
      # Where is the liana testdata located?
      liana_path <- system.file(package = 'liana')       
      # Read the testdata from its location.
      testdata <- 
        readRDS(file.path(liana_path, "testdata", "input", "testdata.rds"))   
      
      # removing superfluous values
      rm(liana_path)
      
      
    } else {
      
      # error if its not one of the supported data sets
      stop("Testdata name not recognized!")
      
    }
    
    
    # Return the seurat.
    return(testdata)
    
  } # end of function
}

extract_top_ranks <- function(results) {
  
  top_ranks_analysis <- results$top_ranks_analysis
  
  where_overlap <- str_detect(names(top_ranks_analysis), "Overlap")
  
  collated_top_ranks_overlap <- top_ranks_analysis[where_overlap] %>%
    bind_rows() 
  
  
  collated_top_ranks_overlap <- collated_top_ranks_overlap %>%
    arrange(dilution_prop) %>%
    pivot_longer(cols = !(starts_with("dilution_prop")), names_to = "Method") %>%
    arrange(Method) %>%
    rename("Overlap" = value) 
  
  rm(where_overlap)
  
  return(collated_top_ranks_overlap)
}

















flatten_names <- function(three_tier_list, depth) {
  
  new_three_tier <-
    map_depth(three_tier_list, depth, function(two_tier_list) {
      two_tier_list %>%
        map2(names(.), function(one_tier_list, one_tier_list_name)
          rename_list(one_tier_list, one_tier_list_name)) %>%
        flatten()
    })
  
  return(new_three_tier)
  
}

rename_list <- function(list_element, str){
  
  new_list <- setNames(list_element,
                       str_glue("{str}_{names(list_element)}"))
  
  return(new_list)
  
}





reformat_Results <- function(results) {
  
  # At most, there are six outputs. They fall into three pairs of two that are
  # each formatted the same way. Liana results and top ranks are formatted
  # the same way, resources and top_ranks analysis are formatted the same way,
  # and runtime and testdata is formatted the same way.
  
  
  # We start by transposing results
  results <- transpose(results)
  
  
  # This is the segment of the results containing runtime and testdata
  segment_runtime_test <-
    results[names(results) %in% intersect(names(results),
                                          c("runtime",
                                            "testdata"))]
  
  # This is the segment of the results containing resources and analysis
  segment_resources_analysis <-
    results[names(results) %in% intersect(names(results),
                                          c("resources_OP",
                                            "top_ranks_analysis"))]
  
  # This is the segment of the results containing ranks and results
  segment_results_ranks <-
    results[names(results) %in% intersect(names(results),
                                          c("liana_results_OP",
                                            "top_ranks_OP"))]
  
  # Runtime_test is already corrrectly formatted.
  
  # Format resources_analysis, we transpose at a deeper level
  segment_resources_analysis <- segment_resources_analysis %>%
    map_depth(.depth = 1, transpose) %>%
    flatten_names(depth = 1)
  
  # Format results_ranks
  segment_results_ranks <- segment_results_ranks %>%
    map_depth(.depth = 1, transpose) %>%
    map_depth(.depth = 2, transpose) %>%
    flatten_names(depth = 2)
    
  restructured_results <- list(segment_results_ranks, 
                               segment_resources_analysis,
                               segment_runtime_test) %>%
    flatten()

   return(restructured_results) 
  
}






