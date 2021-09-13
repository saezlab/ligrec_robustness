

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

calculate_Runtime <- function(runtime) {

  # save the names of the time-points for later
  runtime_labels <- names(runtime)
  
  # convert run time to numeric so we can perform arithmetic operations on
  # them. In this case we need it for subtractions, to calculate the duration
  # between checkpoints
  runtime_numeric <- as.numeric(runtime)
  
  # We calculate the passage of time between checkpoints in the 
  # resource_Robustness().
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
  
  return(runtime)
}


overlap_line_Plot <- function(top_ranks_overlap) {
  
  plot_line <- 
    ggplot(data = top_ranks_overlap, aes(dilution_prop,
                                           Overlap,
                                           group = Method,
                                           color = Method)) + 
    
    
    geom_point(alpha = alpha) +
    stat_summary(alpha = 0.6,
                 fun   = mean, 
                 geom  = "line") +
    
    
    scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0,100)) +
    scale_x_continuous(breaks = seq(0, 100, 20), limits = c(0,100)) +
    
    ggtitle("Robustness of Method Predictions") +
    ylab("Overlap of Top Ranks [%]") +
    xlab("Dilution of Resource [%]") +
    labs(subtitle = "Point scatter plot.",
         caption = plotting_caption,
         color = "Method") +
    
    theme_bw() +
    
    theme(plot.caption = element_text(hjust = 0),
          legend.position = "bottom")
  
  return(plot_line)
  
}

overlap_box_Plot <- function(top_ranks_overlap, plotting_caption) {
  
  plot_box <- 
    ggplot(data = tr_overlap_for_plot, aes(x = dilution_prop, 
                                           y = Overlap, 
                                           group = dilution_prop,
                                           color = Method)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_point(alpha = alpha) + 
    
    scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0,100)) +
    scale_x_continuous(breaks = seq(0, 100, 20)) +
    
    
    ggtitle("Robustness of Method Predictions") +
    ylab("Overlap of Top Ranks [%]") +
    xlab("Dilution of Resource [%]") +
    labs(subtitle = "Boxplot by Method.",
         caption = plotting_caption,
         color = "Method") +
    
    theme_bw() + 
    
    
    theme(plot.caption = element_text(hjust = 0),
          legend.position = "bottom") +     
    
    facet_wrap(~Method, nrow = 2, ncol = 3, scales = "free")
  
  return(plot_box)
  
}


