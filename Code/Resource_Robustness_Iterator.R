#------------------------------------------------------------------------------#
# 0. Setup ---------------------------------------------------------------------
{
  # 0.1 Overview of Goals:
  {
    # The Idea with this script is to:
    # .	run all methods on a single resource (OP)
    # .	take the topmost ranked interactions for each method -> R-zero
    # .	create a modified OP resource in which the topmost ranked interactions
    #     remain, but of the remainder x% of the interactions have been diluted 
    #     and  replaced new genes not in the resource derived from the test 
    #     data, x =10,20,40% etc.
    # .	Rerun methods on diluted omnipath resource and extract the top ranked 
    #   CCI -> R-modified
    # .	plot percentage of R-zero in R-modified over x and investigate result
    
    
  } # end of subpoint
  
  
  # 0.2 Loading Packages and Starting Runtime
  {

    require(tidyverse)
    require(Seurat)
    require(liana)
    require(lubridate)
    
  } # end of subpoint
  
  
}


source("Code/Ranking_Misc_Functions.R")
source("Code/Resource_Dilution_Functions.R")
source("Code/Resource_Robustness_Functions.R")
source("Code/Resource_Iterator_Functions.R")




#------------------------------------------------------------------------------#
# 1. Master Seed list --------------------------------------------------------
{
  # Define script params
  {
    source("Code/Iterator_Parameter_Functions.R")
    # Define the time of run to uniquely tag every save file
    time_of_run <-  Sys.time() %>%
      as.character()       %>%
      gsub(':', '-', .)    %>% 
      gsub(' ', '_', .)
    
  }
  
    # How many permutations of dilution should be performed?
    master_seed_list <- as.list(c(1:2))
    
    
    ## Process Master Seed List
    {
      # By naming each seed we can use this to label data conveniently later
      # By formatting seeds as a list we can lapply over them for an easy-to-
      # untangle results format.
      seed_names <- c()
      
      # Name each element of master_seed_list appropriately name it
      for (seed in master_seed_list) {
        seed_names <- 
          c(seed_names, str_glue("Seed_", seed))
      }
      
      names(master_seed_list) <- seed_names
      
      # Remove clutter
      rm(seed_names, seed)
    }


}   



#------------------------------------------------------------------------------#
# 2. Iterate resource_Robustness() ---------------------------------------------
{
  # resource_Robustness is an entire script that can be iterated as a function
  # There is randomness in dilution. Each master seed passed to 
  # dilute_Resource() gives us one permutation of many theoretically possible
  # dilutions. By iterating over master_seed, we can produce many permutations
  # and tally up their results. In this way, master_seed serves as an index too.
  
  # Apply resource_Robustness(), provide every argument but master_seed
  collated_robustness_results <- lapply(master_seed_list, 
                    wrap_resource_Robustness,
                    time_of_run = time_of_run)
  
                    
}



#------------------------------------------------------------------------------#
# 3. Reformatting Results --------------------------------------------------------
{
  # In this segment we extract the data from the results object, which is poorly
  # formatted by default, and put it into a more appropriate hierarchy. We then
  # extract the most relevant sublists for the rest of the analysis.

  # We name our restructured results more informatively, then extract the most
  # relevant sublists from them for the rest of the analysis
  collated_robustness_results <- 
    reformat_Results(results = collated_robustness_results)

}




#------------------------------------------------------------------------------#
# 5. Collate all iterations of  top_ranks_analysis -----------------------------
{
  collated_top_ranks_overlap <- extract_top_ranks(collated_robustness_results)
}



#------------------------------------------------------------------------------#
# 6. Plotting of Aggregate Results ---------------------------------------------
{
  # Preparing Plotting Inputs
  {
    # We format 
    alpha <- 0.4
    
    # The dilution proportion and overlap are clearer in percentage
    # NAs can't be displayed in the plot anyway and cause uneccesary warnings
    # And we rename the methods from the liana++ internal string to their 
    # official names.
    tr_overlap_for_plot <-  collated_top_ranks_overlap  %>%
      as.data.frame()                             %>%
      mutate(dilution_prop = dilution_prop * 100) %>%
      mutate(Overlap       = Overlap       * 100) %>%
      as_tibble()                                 %>%
      drop_na()                                   %>%
      mutate("Method" = recode(Method,
                               "call_connectome" = "Connectome",
                               "squidpy"         = "CellPhoneDB",
                               "call_natmi"      = "NATMI",   
                               "call_italk"      = "iTALK", 
                               "call_sca"        = "SingleCellSignalR", 
                               "cellchat"        = "CellChat")) 
     
    # automatically generate a plot description based on wrap_resource_Robustness
    # defaults and the tr_overlaps_for_plot data structure
    
    plotting_caption <- 
      auto_plot_Description(tr_overlap_for_plot, 
                            formals(wrap_resource_Robustness),
                            formals(summarise_Metadata), 
                            time_of_run)
    

  }

  
  # Generating and printing Plots
  {
    plot_line <- 
      ggplot(data = tr_overlap_for_plot, aes(dilution_prop,
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
    
    
    
    print(plot_line)
    print(plot_box)
  }
  
  # Removing Clutter
  rm(tr_overlap_for_plot, alpha, plotting_caption)
  
}

#------------------------------------------------------------------------------#
# 4. Reformatting Runtime and Metadata -----------------------------------------
{
  # 4.1 Calculate Runtime
  {
    # By flattening we get one long seqence of time points and the portions of
    # the script they belong to
    runtime <- collated_robustness_results$runtime %>% 
      flatten() %>%
      calculate_Runtime()


  }
  
  #4.2 format metadata
  metadata <- summarise_Metadata(runtime, 
                                 time_of_run,
                                 formals(wrap_resource_Robustness),
                                 master_seed_list)
  
  rm(runtime, master_seed_list)
  
}

#------------------------------------------------------------------------------#
# 7. Saving Results ------------------------------------------------------------
{
  # Save the plot automatically to the outputs folder, if desired
  if (save_results == TRUE) {
    if (run_mode == "real") {
      test_run_comment <- ""
      
    } else if (run_mode == "trial_run") {
      test_run_comment <- "TRIAL_RUN_"
      
    }
    
    if (formals(wrap_resource_Robustness)$preserve_topology == FALSE) {
      topology_comment <- "_random_topology_"
      
    } else if (formals(wrap_resource_Robustness)$preserve_topology == TRUE) {
      topology_comment <- "_preserved_topology_"
      
    }
    
    # Define the time of run to uniquely tag every save file
    time_of_run <-  Sys.time() %>%
      as.character()       %>%
      gsub(':', '-', .)    %>% 
      gsub(' ', '_', .)
    
    top_ranks_vector <- 
      unlist(as.list(formals(wrap_resource_Robustness)$number_ranks)[-1])
    
    # Generate the filepaths to save the data under
    box_plot_png_name <-
      str_glue(
        test_run_comment,
        "Boxplot_Resource_Dilution_",
        formals(wrap_resource_Robustness)$testdata_type,
        topology_comment,
        formals(wrap_resource_Robustness)$feature_type,
        "_top",
        median(top_ranks_vector),
        "_",
        time_of_run,
        ".png"
      )
    
    line_plot_png_name <-
      str_glue(
        test_run_comment,
        "Lineplot_Resource_Dilution_",
        formals(wrap_resource_Robustness)$testdata_type,
        topology_comment,
        formals(wrap_resource_Robustness)$feature_type,
        "_top",
        median(top_ranks_vector),
        "_",
        time_of_run,
        ".png"
      )
    
    
    env_save_path <- 
      str_glue(
        "Outputs/",
        test_run_comment,
        "DilutionEnv_",
        formals(wrap_resource_Robustness)$testdata_type,
        topology_comment,
        formals(wrap_resource_Robustness)$feature_type,
        "_top",
        median(top_ranks_vector),
        "_",
        time_of_run,
        ".RData"
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
    
    
    # Store save locations for plots and session info in script_params metadata
    # script_params$metadata[["box_plot_png_name"]]  <- box_plot_png_name
    # script_params$metadata[["line_plot_png_name"]] <- line_plot_png_name
    # script_params$metadata[["env_save_path"]]      <- env_save_path
    # 
    # script_params$metadata[["Session_Info"]]       <- sessionInfo
    # 
    
    
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
}
