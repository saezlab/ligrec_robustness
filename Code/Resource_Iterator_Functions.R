

  
  
  
  
  
  
  
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
  # and metadata and testdata is formatted the same way.
  
  
  # We start by transposing results
  results <- transpose(results)
  
  
  # This is the segment of the results containing metadata and testdata
  segment_meta_test <-
    results[names(results) %in% intersect(names(results),
                                          c("metadata",
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
  
  # Format  meta_test
  segment_meta_test[["metadata"]]  <- segment_meta_test[["metadata"]] %>%
    map_depth(.depth = 0, transpose) %>%
    flatten_names(depth = 0)
  
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
                               segment_meta_test) %>%
    flatten()

   return(restructured_results) 
  
}