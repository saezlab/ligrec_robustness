#------------------------------------------------------------------------------#
# 0. Introduction and Goals ----------------------------------------------------
{
  
  # This is the Iterator_Processing_Functions.R script. These function work as 
  # miscellaneous intermediary steps in the process of running the iterator. 
  # Information on how the iterator runs can be found in 
  # RD_Robustness_Iterator.R.
  
}



#------------------------------------------------------------------------------#
# 1. Defining Functions --------------------------------------------------------


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
  #' random_Dilute makes only a small effort to preserve the topology of the rows
  #' that are diluted from resource, preserve_Dilute makes a far greater effort.
  #' Choose TRUE for preserve_Dilute and FALSE for random_Dilute.
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
  #' c("call_conncectome, "call_italk"). This argument is also very impactful for 
  #' the runtime of this function. For example on widnows, natmi and cellchat are
  #' the slowest methods.
  #' 
  #' @param number_ranks A named list. Each item is named after a method and is 
  #' equal to the number of top interactions considered relevant for that method. 
  #' 
  #' For example, item one on the list may be called call_connectome and be equal 
  #' to 500. This would signal to the function that for call_connectome, the top 
  #' 500 CCI's are considered relevant and that these 500 are the ones that are to
  #' be compared between the dilutions. Usually you consider the same number of
  #' top_ranks for each method relevant.
  #'  
  #' @param cellchat_nperms Cellchat is one of the slower methods, for test runs
  #' it may be useful to set this parameter to 10 to speed up the analysis. 
  #' Unless you're using this function for a trial run, this should be left at its
  #' default.
  #'
  #' @param sink_otuput TRUE or FALSE. Should resource_Robustness() save a full 
  #' log of the Console Output to the log folder? Warnings and messages will not 
  #' be visible in the console output if this is enabled, so unless there is a 
  #' reason why such a record is necessary this option is not recommended.
  #' 
  #' @param liana_warnings Either TRUE, FALSE or "divert". A less extreme 
  #' alternative to sinking the output. Should the warnings from LIANA++, which 
  #' are often repetitive and unhelpful, be either suppressed or alternatively 
  #' diverted to the log folder? When these types of warning are left in, they can
  #' often displace valuable warnings. Be careful with this  setting, as 
  #' suppressing warnings is obviously risky.
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


# reformat_Results()
{
  #' Reformats the outputs from lapplying resource_Robustness over the master
  #' seed list into a better organized hierarchy
  #' 
  #' @description The standard outputs from lapplying resource_Robustness() are 
  #' highly nested and unintuitive, this being a product of how they are 
  #' created. This function uses transpositions, flattenings and renaming to 
  #' create a flatter and easier to understand result object from iterating 
  #' resource_Robustness(), which can be used for further data extraction more 
  #' easily.
  #' 
  #' @param results The output from lapplying resource_Robustness over the 
  #' master seed list.
  #' 
  #' @return A better structured version of the input.
  
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
    
    # Recombine our results again and return them.
    restructured_results <- list(segment_results_ranks, 
                                 segment_resources_analysis,
                                 segment_runtime_test) %>%
      flatten()
    
    return(restructured_results) 
    
  }
}


# flatten_names()
{
  #' A helper function for restructure results, a variant on flatten() that 
  #' preserves names
  #' 
  #' @description This function takes a list with minimum three layers (a 
  #' sublist and a sub-sub list) and renames the lowest elements to include the 
  #' name of their parent list one step up in the hierarchy, and then flattens 
  #' that level of the list. In short, it removes one level of hierarchy from 
  #' the bottom, but makes sure that the names of the items indicate what 
  #' sublist they used to be a part of. 
  #' 
  #' This helps avoid situations where the lowest tier items in a list all have 
  #' the same names, the list gets flattened and you can't tell where the items
  #' came from anymore. 
  #' 
  #' @param high_tier_list The list to flatten the ends of.
  #' 
  #' @param depth The depth of the lowest items in the list. 
  #' 
  #' @return The input list with it's lowest level of hierarchy removed, and the
  #' lowest level items in it renamed to show whqat sublist they used to be in.
  
  
  flatten_names <- function(high_tier_list, depth) {
    
    new_high_tier <-
      # map to the user specififed depth
      map_depth(high_tier_list, depth, function(two_tier_list) {
        
        # Once at the second lowest tier, grab the name of the lists here
        # and the lists themselves.
        # We then map to the lowest tier, with the names in hand and rename the 
        # lowest tier elements.
        two_tier_list %>%
          map2(names(.), function(one_tier_list, one_tier_list_name)
            rename_list(one_tier_list, one_tier_list_name)) %>%
          flatten()
        
      })
    
    # return our new list
    return(new_high_tier)
    
  }
}


# rename_list()
{
  #' Helper function for reformat_Results, takes a list element and adds a tag 
  #' to its name
  #' 
  #' @param list_element A list element to rename.
  #' 
  #' @param str The tag you want to slap onto the list elements name.
  #' 
  #' @description Usually this is passed the lowest elements in a nested list
  #' as list_elements and then the name of the sub-list they are a part of as
  #' str.
  #' 
  #' @return A renamed list_element.
  
  rename_list <- function(list_element, str){
    
    new_list <- setNames(list_element,
                         str_glue("{str}_{names(list_element)}"))
    
    return(new_list)
    
  }
}


# extract_top_ranks()
{
  #' Extracts top_rank_overlaps from the restructured resource_Robustness() 
  #' output and formats and collates it into one convenient tibble
  #' 
  #' @param results The restructured resource_Robustness() output, i.e. the 
  #' output from reformat_Results().
  #' 
  #' @return A tibble that collates all the top_ranks overlap data from the 
  #' input.
  
  extract_top_ranks <- function(results) {
    
    # get just the top_ranks anylsis information
    top_ranks_analysis <- results$top_ranks_analysis
    
    # narrow it down to top_ranks_overlaps tibbles by searching for entries that 
    # contain the word "overlap"
    where_overlap <- str_detect(names(top_ranks_analysis), "Overlap")
    
    # bind all these separate tibbles into one 
    collated_top_ranks_overlap <- top_ranks_analysis[where_overlap] %>%
      bind_rows() 
    
    # reorganize the tibble to be more usable and easily understandable 
    collated_top_ranks_overlap <- collated_top_ranks_overlap %>%
      arrange(dilution_prop) %>%
      pivot_longer(cols = !(starts_with("dilution_prop")), names_to = "Method") %>%
      arrange(Method) %>%
      rename("Overlap" = value) 
    
    # remove unnecessary clutter
    rm(where_overlap)
    
    # returnt eh collated top_ranks_overlap information
    return(collated_top_ranks_overlap)
  }
}


# auto_plot_Description()
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

