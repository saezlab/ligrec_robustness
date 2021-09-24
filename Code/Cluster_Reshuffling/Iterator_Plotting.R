#------------------------------------------------------------------------------#
# 0. Introduction and Goals ----------------------------------------------------
{
  # This is the Iterator_Plotting.R script. It defines functions relevant for 
  # the plots created in the iterator.
}



#------------------------------------------------------------------------------#
# 1. Define Plotting Functions  ------------------------------------------------

# overlap_line_Plot()
{
  
  #' Draw a line/scatter plot from the collated top_ranks_overlap data
  #' 
  #' @param top_ranks_overlap As a tibble, formatted like the collated_top_ranks
  #' overlap data.
  #' 
  #' @param plotting_caption A caption for your plot, as a string. Potentially
  #' could be an auto_plot_Description() output.
  #' 
  #' @return A ggplot of the top_ranks_overlap as a line and scatter plot, 
  #' complete with a caption.
  
  
  overlap_line_Plot <- function(top_ranks_overlap, plotting_caption) {
    
    plot_line <- 
      ggplot(data = top_ranks_overlap, aes(dilution_prop,
                                           Overlap,
                                           group = Method,  # for colors
                                           color = Method)) + 
      
      
      # add transparency because some points may overlap
      geom_point(alpha = 0.4) +
      stat_summary(alpha = 0.6,
                   fun   = mean, 
                   geom  = "line") + # add a mean trendline per color
      
      
      # make sure the scale is the same on both axes, since both are in percent
      scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0,100)) +
      scale_x_continuous(breaks = seq(0, 100, 20), limits = c(0,100)) +
      
      # add text
      ggtitle("Robustness of Method Predictions") +
      ylab("Overlap of Top Ranks [%]") +
      xlab("Dilution of Resource [%]") +
      labs(subtitle = "Line / point scatter plot.",
           caption = plotting_caption,
           color = "Method") +
      
      # remove grey background so that points are bteer visible
      theme_bw() +
      
      # shift caption and legend position
      theme(plot.caption = element_text(hjust = 0),
            legend.position = "bottom")
    
    
    # return the plot
    return(plot_line)
    
  }
  
}


# overlap_box_Plot()
{
  
  #' Draw multiple boxplots plot from the collated top_ranks_overlap data
  #' 
  #' @param top_ranks_overlap As a tibble, formatted like the collated_top_ranks
  #' overlap data.
  #' 
  #' @param plotting_caption A caption for your plot, as a string. Potentially
  #' could be an auto_plot_Description() output.
  #' 
  #' @return A ggplot of the top_ranks_overlap as a box plot, complete with a 
  #' caption.

  overlap_box_Plot <- function(top_ranks_overlap, plotting_caption) {
    
    plot_box <- 
      ggplot(data = top_ranks_overlap, aes(x = dilution_prop, 
                                           y = Overlap, 
                                           group = dilution_prop, # for faceting 
                                           color = Method)) +     # ...and colors
      
      # Add boxplots, no outliers because we add points in a second
      geom_boxplot(outlier.shape = NA) + 
      geom_point(alpha = 0.4) +  # plot semi transparent plots over them so
      # you can see where the distirbution comes from
      
      
      # make sure the scale is the same on both axes, since both are in percent
      scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0,100)) +
      scale_x_continuous(breaks = seq(0, 100, 20)) +
      
      
      # add text
      ggtitle("Robustness of Method Predictions") +
      ylab("Overlap of Top Ranks [%]") +
      xlab("Dilution of Resource [%]") +
      labs(subtitle = "Boxplot by Method.",
           caption = plotting_caption,
           color = "Method") +
      
      # remove gray background so the transparent points are more visible
      theme_bw() + 
      
      # adjust legend and caption position
      theme(plot.caption = element_text(hjust = 0),
            legend.position = "bottom") +     
      
      # facet wrap for each method, so they don't overlap and are easily 
      # comparable between subplots
      facet_wrap(~Method, nrow = 2, ncol = 3, scales = "free")
    
    
    
    # return the plot we crafted
    return(plot_box)
    
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
  
  clust_plot_Description <- function(mismatch_props,
                                     trial_run,
                                     testdata_type,
                                     master_seed_list,
                                     number_ranks,
                                     time_of_run) {
    
    ## General comment, on testdata type, feature_type and topology
    {
      
      
      general_comment <-
        str_glue(
          "This plot was created using the ",
          testdata_type,
          " data.  "
        )

    }
    
    
    
    ## Dilution comment, on proportions
    {
      
      if (length(mismatch_props) > 1) {
        
        reshuffle_comment <- str_glue(
          "The cluster annotations were reshuffled in ",
          (mismatch_props[[2]]-mismatch_props[[1]]) *100, 
          " % intervals to a maximum of ",
          max(unlist(mismatch_props))*100,
          " %. When an annotation was being reshuffled, it was replaced by a
        random sign from all annotations that did not match itself."
        )
        
      } else {
        reshuffle_comment <- str_glue(
          "The cluster annotations were reshuffled to ",
          (mismatch_props[[1]]) *100, 
          " %. When an annotation was being reshuffled, it was replaced by a
        random sign from all annotations that did not match itself."
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
        general_comment,
        "\n",
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

