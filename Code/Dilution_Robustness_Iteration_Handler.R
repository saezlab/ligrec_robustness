#------------------------------------------------------------------------------#
# A. Setup ---------------------------------------------------------------------
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
    # .	Rerun methods on diluted omnipath resource, get top ranks -> R-modified
    # .	plot percentage of R-zero in R-modified over x and investigate result
    
    
  } # end of subpoint
  
  
  # 0.2 Loading Packages and Starting Runtime
  {

    require(tidyverse)
    require(Seurat)
    require(liana)
    
    library(lubridate)
    
    # each runtime is named after the section that it marks the conclusion of
    
  } # end of subpoint
  
}



#------------------------------------------------------------------------------#
# B. Script  Parameters --------------------------------------------------------
{
  testdata_type  <- c("liana_test") # choose "liana_test" or "seurat_pbmc"
  
  # If feature_type is generic, dilution will be completed with any genes
  # in the seurat count matrix. If dilution is variable, only the variable
  # features are used for dilution.
  feature_type <- c("variable") # choose "generic" or "variable"
  
  preserve_topology <- FALSE  # TRUE = preserve_Dilute(), FALSE = random_Dilute()
  
  dilution_props <- c(seq(0.10, 0.10, 0.10)) # should be consistent between tests
  
  number_ranks   <- list("call_connectome" = 30, 
                         "call_natmi"      = 30,
                         "call_italk"      = 30,
                         "call_sca"        = 30,
                         "cellchat"        = 30)
  
  # Format a master_seed_list that provides a different master_seed for each
  # iteration of dilution_Robustness()
  master_seed_list <- as.list(c(1:4))
  
  seed_names <- c()
  
  # Name each element of master_seed_list to have a name on results later
  for (seed in master_seed_list) {
    seed_names <- 
      c(seed_names, str_glue("Seed_", as.character(seed)))
  }
  
  names(master_seed_list) <- seed_names
  
  # Remove clutter
  rm(seed_names, seed)
  

  # Define Outpputs
  outputs = c(
              "liana_results_OP",
              "resources_OP",
              "top_ranks_OP",
              "top_ranks_analysis",
              "metadata", 
              "testdata"
              )
  

  # All the methods we're using (almsot all six of liana)
  # squidpy won't be used unthetil I get it to work on windows
  methods_vector <- c(
                      'call_connectome',
                      #'call_natmi', 
                      'call_italk',
                      'call_sca' #,
                      #'cellchat'
                      )
  
  
  cellchat_nperms <- 10 # number of cellchat permutations, default 100
  
  run_mode <- "trial_run" # select between trial_run and real
  
  save_results <- FALSE # should results be saved?
  
  sink_output <- FALSE # If the output is sunk, all console outputs and warning
  # messages go to a txt file in Outputs folder, but you won't be able to see
  # them in the console. In essence, logs will be generated instead of console
  # outputs
  
  liana_warnings <- "divert"
  
  
  # define dilution proportions
  # dilution props is a user defined sequence in the setup section
  dilution_names <- c()
  
  for (i in dilution_props) {
    dilution_names <- 
      c(dilution_names, str_glue("OmniPath_", as.character(i*100)))
  }
  
  dilution_props <- as.list(dilution_props)
  names(dilution_props) <- dilution_names
  
  rm(dilution_names, i)
  
  
  # Summarize all the above parameters into one compact item
  script_params <- list("master_seed_list"  = master_seed_list,
                        "dilution_props"    = dilution_props, 
                        "number_ranks"      = number_ranks, 
                        "feature_type"      = feature_type, 
                        "methods_vector"    = methods_vector, 
                        "testdata_type"     = testdata_type,
                        "preserve_topology" = preserve_topology,
                        
                        "outputs"         = outputs,
                        "cellchat_nperms" = cellchat_nperms, 
                        "run_mode"        = run_mode, 
                        "save_results"    = save_results,
                        "sink_output"     = sink_output,
                        "Session_Info"    = sessionInfo(),
                        "liana_warnings"  = liana_warnings)
  
  # Remove the parameters we just summarized
  rm(master_seed_list, dilution_props, number_ranks, feature_type, 
     methods_vector, testdata_type, preserve_topology, outputs, cellchat_nperms,
     run_mode, save_results, sink_output, liana_warnings)
  

}   



#------------------------------------------------------------------------------#
# C. Iterate dilution_Robustness() ---------------------------------------------

# dilution_Robustness is an entire script that can be iterated as a function
# There is randomness in dilution. Each master seed passed
# to dilute_Resource() gives us one permutation of many theoretically possible
# dilutions. By iterating over master_seed, we can produce many permutations
# and aggregate their results.

# Apply dilution_Robustness(), provide every argument but master_seed
results <- lapply(script_params$master_seed_list, 
                  dilution_Robustness,
                  
                  testdata_type     = script_params$testdata_type,
                  feature_type      = script_params$feature_type,
                  preserve_topology = script_params$preserve_topology,
                  dilution_props    = script_params$dilution_props,
                  number_ranks      = script_params$number_ranks,
                  outputs           = script_params$outputs,
                   
                  methods_vector    = script_params$methods_vector,
                  cellchat_nperms   = script_params$cellchat_nperms,
                  sink_output       = script_params$sink_output,
                  liana_warnings    = script_params$liana_warnings)



# Remove uneccesary Parameters
rm(number_ranks, cellchat_nperms, feature_type,
   outputs, preserve_topology, run_mode, save_results, sink_output,
   testdata_type)


# Reorganize Outputs for Top_ranks_overlaps
top_ranks_overlaps <- list()

for (seed in names(master_seed_list)) {
  top_ranks_overlaps[[seed]] <- results[[seed]]$top_ranks_analysis$Overlap
  top_ranks_overlaps[[seed]] <- top_ranks_overlaps[[seed]] %>%
    as.data.frame()
}

top_ranks_overlaps <- bind_cols(top_ranks_overlaps)

top_ranks_aggreg_mean <- list()
top_ranks_aggreg_sd   <- list()

for (method in methods_vector) {
  
  vector_sd   <- c()
  vector_mean <- c()
  
  for (row in 1:(length(dilution_props)+1)) {
    
    mean_overlap <- mean(as.numeric(
      top_ranks_overlaps[row, grepl(method, names(top_ranks_overlaps))]))
    
    sd_overlap   <- sd(as.numeric(
      top_ranks_overlaps[row, grepl(method, names(top_ranks_overlaps))]))
    
    vector_mean <- c(vector_mean, mean_overlap)
    vector_sd   <- c(vector_sd,   sd_overlap)
    
  }
  
  top_ranks_aggreg_mean[[str_glue(method, "_mean")]] <- vector_mean
  top_ranks_aggreg_sd[[str_glue(method, "_sd")]]     <- vector_sd
}

agg_top_ranks_overlap_mean <- as_tibble(top_ranks_aggreg_mean) 

agg_top_ranks_overlap_sd   <- as_tibble(top_ranks_aggreg_sd) 

agg_top_ranks <- tibble(agg_top_ranks_overlap_mean,
                        agg_top_ranks_overlap_sd,
                        .name_repair = c("universal"))

agg_top_ranks <- agg_top_ranks %>%
  mutate(dilution_prop = c(0, dilution_props))  %>%
  unnest(cols = c(dilution_prop))               %>%
  relocate("dilution_prop")


#------------------------------------------------------------------------------#
# 4. Visualizing the results -------------------------------------------------
{ 

  
  # 4.1 Plotting, labeling and saving top_ranks_overlap 
  {
    
    
    # The plot is better in percent than proportion
    tr_overlap_for_plot <-  agg_top_ranks * 100
    
    # Automatically assemble a file name and plot subtitle
    if (preserve_topology == FALSE) {
      
      topology_comment <- "using random_Dilute()"
      
    } else if (preserve_topology == TRUE) {
      
      topology_comment <- "using preserve_Dilute()"
      
    }
    
    plotting_subtitle <- str_glue(feature_type,
                                  " dilution, top ",
                                  as.character(median(unlist(number_ranks))),
                                  " ranks, ",
                                  testdata_type,
                                  " data, ",
                                  run_mode,
                                  " results, ",
                                  topology_comment,
                                  ", ",
                                  length(master_seed_list),
                                  " permutations")
    
    plot_png_name     <- str_glue(run_mode,
                                  "_",
                                  testdata_type, 
                                  "_top",
                                  as.character(median(unlist(number_ranks))),
                                  "_res",
                                  as.character(length(dilution_props)),
                                  "_",
                                  feature_type,
                                  "_dil_on_",
                                  as.character(Sys.Date()),
                                  ".png")
    
    # Plot top_ranks_overlap with lines and points at each value
    overlap_plot <-  ggplot(data = tr_overlap_for_plot) + 
      geom_line(mapping = aes(dilution_prop, 
                              call_connectome_mean, 
                              color =  "Connectome")) +
      
      geom_line(mapping = aes(dilution_prop, 
                              call_natmi_mean, 
                              color = "NATMI")) + 
      
      geom_line(mapping = aes(dilution_prop,
                              call_italk_mean, 
                              color = "iTALK")) +
      
      geom_line(mapping = aes(dilution_prop, 
                              call_sca_mean, 
                              color = "SCA")) +
      
      geom_line(mapping = aes(dilution_prop, 
                              cellchat_mean, 
                              color = "CellChat")) +
      
      
      
      geom_point(mapping = aes(dilution_prop, 
                               call_connectome_mean, 
                               color =  "Connectome")) +
      
      geom_point(mapping = aes(dilution_prop,
                               call_natmi_mean, 
                               color = "NATMI")) +
      
      geom_point(mapping = aes(dilution_prop, 
                               call_italk_mean,
                               color = "iTALK")) +
      
      geom_point(mapping = aes(dilution_prop, 
                               call_sca_mean, 
                               color = "SCA")) +
      
      geom_point(mapping = aes(dilution_prop, 
                               cellchat_mean, 
                               color = "CellChat")) +
      

      geom_errorbar(aes(ymin = call_connectome_mean - call_connectome_sd,
                        ymax = call_connectome_mean + call_connectome_sd,
                        x = dilution_prop), 
                    width=.2,
                    position=position_dodge(0.05)) + 
      
      geom_errorbar(aes(ymin = call_natmi_mean - call_natmi_sd,
                        ymax = call_natmi_mean + call_natmi_sd,
                        x = dilution_prop), 
                    width=.2,
                    position=position_dodge(0.05)) + 
      
      geom_errorbar(aes(ymin = call_italk_mean - call_italk_sd,
                        ymax = call_italk_mean + call_italk_sd,
                        x = dilution_prop), 
                    width=.2,
                    position=position_dodge(0.05)) + 
      
      geom_errorbar(aes(ymin = call_sca_mean - call_sca_sd,
                        ymax = call_sca_mean + call_sca_sd,
                        x = dilution_prop), 
                    width=.2,
                    position=position_dodge(0.05)) + 
      
      geom_errorbar(aes(ymin = cellchat_mean - cellchat_sd,
                        ymax = cellchat_mean + cellchat_sd,
                        x = dilution_prop), 
                    width=.2,
                    position=position_dodge(0.05)) + 
      
      
      
      # Show full breadth of 100-0 percent overlap
      ylim(0, 100) +
      
      ggtitle("Robustness of Method Predictions") +
      ylab("Overlap of Top Ranks [%]") +
      xlab("Dilution of Resource [%]") +
      labs(subtitle = plotting_subtitle,
           color = "Method")
    
    
    
    # Print the Plot
    print(overlap_plot)
    
    # Save the plot automatically to the outputs folder, if desired
    if (save_results) {
      
      ggsave(plot_png_name, 
             height = 5, width = 8, 
             path = "Outputs")
      
      print(str_glue("Plot saved at ~/Outputs/", plot_png_name, "."))
      
    }
    
    # Remove unnecessary variables
    rm(tr_overlap_for_plot, overlap_plot, plotting_subtitle, topology_comment)
    
    
    
    
  } # end of subpoint
  
 }

