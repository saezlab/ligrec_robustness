#------------------------------------------------------------------------------#
# 0. Introduction and Goals ----------------------------------------------------
{
  
  # This is the Iterator_Processing_Functions.R script. These function work as 
  # intermediary steps in the process of running the iterator. Information on
  # how the iterator runs can be found in RD_Robustness_Iterator.R.
  
}



#------------------------------------------------------------------------------#
# 1. Defining Functions --------------------------------------------------------




# extract_Testdata()
{
  #' Helper function that gets a specific seurat object from the outputs folder
  #' 
  #' @param testdata_type As a string. Which testdata should be retrieved? 
  #' Either "seurat_pbmc" or "liana_test". Seurat_pbmc is the data set used in 
  #' the seurat tutorial, while liana_test is the testdata that comes with 
  #' LIANA++, and is a small subset of seurat_pbmc. 
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

