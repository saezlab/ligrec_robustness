#------------------------------------------------------------------------------#
# A. Setup ---------------------------------------------------------------------
{
  # 0.1 Overview of Goals:
  {
  # The Idea with this script is to:
  # .	run all methods on a single resource (OP)
  # .	take the topmost ranked interactions for each method -> R-zero
  # .	create a modified OP resource in which the topmost ranked interactions
  #     remain, but of the remainder x% of the interactions have been removed 
  #     and  replaced with entirely random pairs of genes derived from the test 
  #     data that do not exist in the resource, x =10,20,40% etc.
  # .	Rerun methods on modified omnipath resource, get top ranks -> R-modified
  # .	plot percentage of R-zero in R-modified over x and investigate result
    
    
  } # end of subpoint
  
  # 0.2 Script Structure:
  {
    # We load preprocessed seurat pbmc data and apply each method combined with 
    # undiluted OmniPath to it.
    # We define a function that can get the top ranked CCI from each method
    # output and run it for the undiluted results.
    # We restructure OmniPath, define a function that dilutes Resources and 
    # restructure our data to be used at various dilutions.
  
    # We have base OmniPath and the top ranks for each method. These basic 
    # inputs allow us to run dilute_Resource(). We create diluted resources for
    # each method.
    # We get the top ranks foreach method at each new dilution.
    # We compare the percentage overlap of the undiluted top_ranks to the 
    # top_ranks at each stage of dilution (for each method).
  
    # We visualize the results.
    
    
  } # end of subpoint
    
  # 0.3 Loading Packages and Starting Runtime
  {
    runtime <- list("Script Start" = Sys.time())
    
    require(tidyverse)
    require(Seurat)
    require(liana)
    
    library(lubridate)
    
    runtime[["Load Packages"]] <- Sys.time() 
    # each runtime is named after the section that it marks the conclusion of
    
  } # end of subpoint
  
}



#------------------------------------------------------------------------------#
# B. Set top_n, dilution props, testdata type ----------------------------------
{
  dilution_props <- c(seq(0.10, 0.3, 0.20)) # should be consistent between tests
  
  number_ranks   <- list("call_connectome" = 20, 
                         "call_natmi"      = 20,
                         "call_italk"      = 20,
                         "call_sca"        = 20,
                         "cellchat"        = 20)
  
  testdata_type  <- c("liana_test") # choose "liana_test" or "seurat_pbmc"
 
  feature_type <- c("variable") # choose "generic" or "variable"
  # If feature_type is generic, dilution will be completed with any genes
  # in the seurat count matrix. If dilution is variable, only the variable
  # features are used for dilution.
  
  # All the methods we're using (almsot all six of liana)
  # squidpy won't be used unthetil I get it to work on windows
  methods_vector <- c('call_connectome',
                      'call_natmi', 
                      'call_italk',
                      'call_sca',
                      'cellchat')
  
  cellchat_nperms <- 10 # number of cellchat permutations, default 100
  
  run_mode <- "trial_run" # select between trial_run and real
  
  save_results <- TRUE # should results be saved?
  
  
}   


  
#------------------------------------------------------------------------------#
# C. Preparing necessary inputs to dilute Resources ----------------------------
{
  # 1. Running LIANA wrapper
  {
    
  # Get seurat or liana test data
  if (testdata_type == "seurat_pbmc") {
    
    testdata <- readRDS(file = "Data/pbmc3k_final.rds")     
    
    
    
  } else if (testdata_type == "liana_test") {
    
    liana_path <- system.file(package = 'liana')       
    testdata <- 
      readRDS(file.path(liana_path, "testdata", "input", "testdata.rds"))   
    
    rm(liana_path)
    
    
    
  } else {
    
    warning("Testdata name not recognized!")
    
  }
  
 
  # Generate Undiluted liana results by running wrapper function
  # Omnipath x the methods vector, on the selected data

  # NATMI results are contaminated with results from earlier runs if you
  # don't specify a special output folder for the results to go in
    
  natmi_output <-  Sys.time()           %>%
                   as.character()       %>%
                   gsub(':', '_', .)    %>% 
                   gsub('-', '_', .)    %>% 
                   gsub(' ', '_', .)    %>%
                   str_glue('Test_', .)
    
  liana_results_OP_0 <- 
    liana_wrap(testdata, 
               method = methods_vector, 
               resource = c('OmniPath'), 
               expr_prop = 0,
               cellchat.params   = list(nboot = cellchat_nperms, 
                                        expr_prop = 0,
                                        thresh = 1),
               call_natmi.params = list(output_dir = natmi_output))
  
  
  
  # Update runtime
  runtime[["First Liana"]] <- Sys.time()
  
  
  # Remove superfluous variables
  rm(natmi_output)
  
  } # end of subpoint
  
  # 2. Get highest ranked interactions for undiluted conditions 
  {
    
    
  # Apply get_top_n_ranks for each method's results on OP_0 (i.e. undiluted)
  top_ranks_OP_0 <- list()
  
  for (method in methods_vector){
    
    top_ranks_OP_0[[method]] <- get_top_n_ranks(data_set = 
                                                liana_results_OP_0[[method]],
                                                top_n = number_ranks[[method]],
                                                method = method)
    
  }
    
  # an interesting note here is that NATMI produces far more unique LR Pairs in 
  # its top rankings than other methods. Cellchat identifies less for example, 
  # filling its top 1000 (or n) by repeating those interactions between many
  # combinations of source cell clusters and target cell clusters.
  # NATMI repeats its interactions less, providing more unique LR pairs than
  # cellchat. In exchange it only finds these interacting pairs more rarely 
  # between cell clusters.
    
    
    
  } # end of subpoint
  
  # 3. Prepare to Modify OmniPath with random genes
  {
  
    
  # Format OmniPath_0 to be easier to work with and to pair it down to the 
  # columns relevant for the methods. Also add the isRandom column which 
  # indicates whether an interaction has been randomly generated and the LR_Pair
  # columns, which helps identify individual interactions
  OmniPath_0 <- select_resource(c('OmniPath'))[["OmniPath"]] %>%
                select(source_genesymbol,
                       target_genesymbol,
                       is_directed,
                       is_stimulation,
                       consensus_stimulation,
                       is_inhibition,
                       consensus_inhibition,
                       category_intercell_source,
                       category_intercell_target,
                       genesymbol_intercell_source,
                       genesymbol_intercell_target,
                       entity_type_intercell_target,
                       sources,
                       references,
                       entity_type_intercell_source,
                       entity_type_intercell_target) %>%
                mutate(isRandom = FALSE) %>%
                unite("LR_Pair", 
                      c(source_genesymbol, target_genesymbol), 
                      remove = FALSE, 
                      sep = "_") %>%
                relocate("LR_Pair", .after = last_col())
    
    # Filter OmniPath to only include interactions between genes which are 
    # represented in the data. 
    # This has no impact on the results, since the removed interactions can't be
    # evaluated by the methods as the necessary genes are missing.
    # The advantage here is that dilution later replaces genes from the resource
    # with genes in the data set
    # If we consider genes in the resource that are also represented in the data
    # 'hits', then we are diluting the resource by inserting hits.
    # Making sure OP only had hits to begin with ensures we dilute hits with 
    # other hits, a fairer comparison than the alternative, which would be 
    # diluting hits and non-hits from OP with hits from the data.
    gene_names <- rownames(testdata@assays$RNA@data)
    
    OmniPath_0 <- OmniPath_0 %>%
      filter(source_genesymbol %in% gene_names) %>%
      filter(target_genesymbol %in% gene_names)
    
    # removing superfluous values
    rm(gene_names)
    
    
    
  } # end of subpoint
    
  # 4. Listify Data sets in the face of dilution steps
  {
    
    
  # Since we are about to perform the same analysis in the steps above but 
  # multiplied by each dilution step, we will turn our data sets into named 
  # lists sorted by method and sub categorized by dilution.
  
    
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
  
  #relist all our data into three concise named lists of named lists
  resources_OP <- list("call_connectome" = list(OmniPath_0 = OmniPath_0),
                       "call_natmi"      = list(OmniPath_0 = OmniPath_0),
                       "call_italk"      = list(OmniPath_0 = OmniPath_0),
                       "call_sca"        = list(OmniPath_0 = OmniPath_0),
                       "cellchat"        = list(OmniPath_0 = OmniPath_0))
  
  
  
  liana_results_OP <- 
    list("call_connectome" = list(OmniPath_0 = liana_results_OP_0$call_connectome),
         "call_natmi"      = list(OmniPath_0 = liana_results_OP_0$call_natmi),
         "call_italk"      = list(OmniPath_0 = liana_results_OP_0$call_italk),
         "call_sca"        = list(OmniPath_0 = liana_results_OP_0$call_sca),
         "cellchat"        = list(OmniPath_0 = liana_results_OP_0$cellchat))

  
  
  
  top_ranks_OP <- 
    list("call_connectome" = list(OmniPath_0 = top_ranks_OP_0$call_connectome),
         "call_natmi"      = list(OmniPath_0 = top_ranks_OP_0$call_natmi),
         "call_italk"      = list(OmniPath_0 = top_ranks_OP_0$call_italk),
         "call_sca"        = list(OmniPath_0 = top_ranks_OP_0$call_sca),
         "cellchat"        = list(OmniPath_0 = top_ranks_OP_0$cellchat))



  
  
  
  # remove old data frames
  rm(OmniPath_0, liana_results_OP_0, top_ranks_OP_0)
  
  
  
  } # end of subpoint
  
}



#------------------------------------------------------------------------------#
# D. Diluting Resources --------------------------------------------------------
{
  # 5. Generate diluted Resources for all methods
  {
    
    
  # Initiating a list of all dilutions
  dilutions_OP <- list()
  
  # Iterate over every method, lapply over every dilution proportion
  for (method in methods_vector){
    
    dilutions_OP[[method]] <- 
      lapply(dilution_props, dilute_Resource, 
             resource = resources_OP[[method]]$OmniPath_0, 
             top_rank_df = top_ranks_OP[[method]]$OmniPath_0, 
             data_set = testdata,
             dilution_type = feature_type)
    
  }
  
  
  
  # Merge OP_0 with the rest of the dilutions, could use mapply but its less 
  # consistent
  for (method in methods_vector) {
    for (dilution in names(dilution_props)) {
      
      resources_OP[[method]][[dilution]] <- dilutions_OP[[method]][[dilution]]
      
    }
  }
  
  
  # Remove uneccesary Variables
  rm(dilutions_OP, method, dilution)
  
  
  
  } # end of subpoint
  
}    
  


#------------------------------------------------------------------------------#
# E. Rerunning Liana and comparing  top ranks ----------------------------------
{
  # 6. Reapply individual methods with diluted resources
  {
  # results still growing somehow, don't know why (natmi and others)
  # in some cases the defaults assosciated with liana_wrap are explicitly applied
    
  # Initialize a list for liana results using diluted resources
  liana_dilutions_OP <- list()
    
  runtime[["Resource Dilution"]] <- Sys.time()
  
  # lapply liana wrap accross the diluted resources for every method 
  
  # NATMI results are contaminated with results from earlier runs if you
  # don't specify a special output folder for the results to go in
  
  natmi_output <-  Sys.time()           %>%
                   as.character()       %>%
                   gsub(':', '_', .)    %>% 
                   gsub('-', '_', .)    %>% 
                   gsub(' ', '_', .)    %>%
                   str_glue('Test_', .)
  
  
  
  for (method in methods_vector) {
    
 
    
    liana_dilutions_OP[[method]] <-
      lapply(resources_OP[[method]][-1], 
             liana_wrap,
             seurat_object = testdata,
             method = method,
             resource = c('custom'),
             expr_prop = 0,
             cellchat.params = list(nboot = cellchat_nperms, 
                                    expr_prop = 0,
                                    thresh = 1),
             call_natmi.params = list(output_dir = natmi_output))
    
    runtime[[str_glue(str_to_title(method), " rerun")]] <- Sys.time()
    
  }

  
  
  
  # Merge with undiluted results, could use mapply but its less consistent
  for (method in methods_vector) {
    for (dilution in names(dilution_props)) {
      
      liana_results_OP[[method]][[dilution]] <- 
        liana_dilutions_OP[[method]][[dilution]][[method]]
      
    }
  }

  
  # Remove uneccesary Variables
  rm(liana_dilutions_OP, method, dilution, natmi_output)

  
  
  } # end of subpoint
  
  # 7. Get top_n_ranks for each method and dilution
  {
    
    
  # lapply get_top_n_ranks over the dilution stages and save results in 
  # top_dilutions list
  top_dilutions_OP <- list()
  

  top_dilutions_OP[["call_connectome"]] <- 
    lapply(liana_results_OP$call_connectome[-1], get_top_n_ranks, 
           method = "call_connectome", top_n = number_ranks$call_connectome)
  
  top_dilutions_OP[["call_natmi"]] <-
    lapply(liana_results_OP$call_natmi[-1], get_top_n_ranks,
           method = "call_natmi", top_n = number_ranks$call_natmi)
  
  top_dilutions_OP[["call_italk"]] <- 
    lapply(liana_results_OP$call_italk[-1], get_top_n_ranks, 
           method = "call_italk", top_n = number_ranks$call_italk)  
  
  top_dilutions_OP[["call_sca"]] <- 
    lapply(liana_results_OP$call_sca[-1], get_top_n_ranks, 
           method = "call_sca", top_n = number_ranks$call_sca)  
  
  top_dilutions_OP[["cellchat"]] <- 
    lapply(liana_results_OP$cellchat[-1], get_top_n_ranks, 
           method = "cellchat", top_n = number_ranks$cellchat)  
  
  

  



  
  
  
  
  # Merge with undiluted results, could use mapply but its less consistent
  for (method in methods_vector) {
    for (dilution in names(dilution_props)) {
      
      top_ranks_OP[[method]][[dilution]] <- 
        top_dilutions_OP[[method]][[dilution]]
      
    }
  }
  
  # Remove superfluous values
  rm(top_dilutions_OP, method, dilution)
  
  
  
  } # end of subpoint

  #8. Evaluate how many of the top 200 interactions overlap between the original
  #   and the dilutions
  {
    
    
  # format top_ranks to have an ID that marks each specific interaction (LR and 
  # the source and target cell)
  top_ranks_OP$call_connectome <- 
      lapply(top_ranks_OP$call_connectome, unite, col = "LR_ID", 
             c(source, target, ligand, receptor), remove = FALSE)
    
    top_ranks_OP$call_natmi <-
      lapply(top_ranks_OP$call_natmi, unite, col = "LR_ID",
             c(source, target, ligand, receptor), remove = FALSE)
    
    
    top_ranks_OP$call_italk <- 
      lapply(top_ranks_OP$call_italk, unite, col = "LR_ID", 
             c(source, target, ligand, receptor), remove = FALSE)
    
    
    top_ranks_OP$call_sca <- 
      lapply(top_ranks_OP$call_sca, unite, col = "LR_ID", 
             c(source, target, ligand, receptor), remove = FALSE)
    
    
  top_ranks_OP$cellchat <- 
    lapply(top_ranks_OP$cellchat, unite, col = "LR_ID", 
            c(source, target, ligand, receptor), remove = FALSE)
  

  
  
  # add a column to see if an interaction is fake
  for (method in methods_vector) {
    for (dilution in c("OmniPath_0", names(dilution_props))) {
      if( !(is_null(top_ranks_OP[[method]][[dilution]]))) {
          
          top_ranks_OP[[method]][[dilution]] <- 
            top_ranks_OP[[method]][[dilution]] %>%
            mutate(isRandom = 
                     !(LR_Pair %in% resources_OP[[method]]$OmniPath_0$LR_Pair))
          
          
      } else {
        
        warning("One of the top_rank tibbles is missing! Moving on.")
        
      }
    }
  }
  
  # remove superfluous values
  rm(method,dilution)
  
  
  # lapply rank_overlap over the top rank tibbles, comparing the dilutions to 
  # the OP_0 at each stage.
  overlaps <- list()
  
  overlaps[['call_connectome']] <- 
    lapply(top_ranks_OP$call_connectome, 
           rank_overlap, 
           main_ranks = top_ranks_OP$call_connectome$OmniPath_0)
  
  
  overlaps[['call_natmi']]      <- 
    lapply(top_ranks_OP$call_natmi, 
           rank_overlap,
           main_ranks = top_ranks_OP$call_natmi$OmniPath_0)
  
  
  overlaps[['call_italk']]    <- 
    lapply(top_ranks_OP$call_italk, 
           rank_overlap, 
           main_ranks = top_ranks_OP$call_italk$OmniPath_0)
  
  
  overlaps[['call_sca']]       <- 
    lapply(top_ranks_OP$call_sca, 
           rank_overlap,
           main_ranks = top_ranks_OP$call_sca$OmniPath_0)
  
  
  overlaps[['cellchat']]      <- 
    lapply(top_ranks_OP$cellchat, 
           rank_overlap,
           main_ranks = top_ranks_OP$cellchat$OmniPath_0)
  
  
  
  # add NAs to the end of the overlaps where dilution wasn't possible
  # this way all the overlaps have the same length for tibble construction
    overlaps <- 
    lapply(overlaps, 
           function(x) { c(x, rep(NA, length(dilution_props)+1-length(x)))})
    
  
 
    
    
    
  
  # reformatting overlap as a tibble
  top_ranks_overlap <- tibble("call_connectome" = overlaps$call_connectome,
                             "call_natmi"      = overlaps$call_natmi,
                             "call_italk"      = overlaps$call_italk,
                             "call_sca"        = overlaps$call_sca,
                             "cellchat"        = overlaps$cellchat)  %>%
                       unnest(cols = all_of(methods_vector))         %>%
                       mutate(dilution_prop = c(0, dilution_props))  %>%
                       unnest(cols = c(dilution_prop))               %>%
                       relocate("dilution_prop")
  
  # removing superfluous values
  rm(overlaps)
  
  } # end of subpoint
  
}



#------------------------------------------------------------------------------#
# F. Visualizing the results ---------------------------------------------------
{ 
  # 9. Plotting, labeling and saving top_ranks_overlap 
  {
    
    
  # The plot is better in percent than proportion
  tr_overlap_for_plot <- top_ranks_overlap * 100
  
  # Automatically assemble a file name and plot subtitle
  plotting_subtitle <- str_glue(feature_type,
                                " dilution, top ",
                                as.character(median(unlist(number_ranks))),
                                " ranks, ",
                                testdata_type,
                                " data, ",
                                run_mode,
                                " results")
  
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
                                            call_connectome, 
                                            color =  "Connectome")) +
                    
                    geom_line(mapping = aes(dilution_prop, 
                                            call_natmi, 
                                            color = "NATMI")) + 
                    
                    geom_line(mapping = aes(dilution_prop,
                                            call_italk, 
                                            color = "iTALK")) +
                    
                    geom_line(mapping = aes(dilution_prop, 
                                            call_sca, 
                                            color = "SCA")) +
                    
                    geom_line(mapping = aes(dilution_prop, 
                                            cellchat, 
                                            color = "CellChat")) +
                    
                    
                    
                    geom_point(mapping = aes(dilution_prop, 
                                             call_connectome, 
                                             color =  "Connectome")) +
                    
                    geom_point(mapping = aes(dilution_prop,
                                             call_natmi, 
                                             color = "NATMI")) +
                    
                    geom_point(mapping = aes(dilution_prop, 
                                             call_italk,
                                             color = "iTALK")) +
                    
                    geom_point(mapping = aes(dilution_prop, 
                                             call_sca, 
                                             color = "SCA")) +
                    
                    geom_point(mapping = aes(dilution_prop, 
                                             cellchat, 
                                             color = "CellChat")) +
                    
                    
                    
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
  rm(tr_overlap_for_plot, overlap_plot, plotting_subtitle)
  
  
  } # end of subpoint
  
}



#------------------------------------------------------------------------------#
# G. Saving the results --------------------------------------------------------
{
  # 10. Calculating run time of script
  {
    
  # stop the stopwatch
  runtime[["Script Epilogue"]] <- Sys.time()
  
  # save the names of the time-points for later
  runtime_labels <- names(runtime)
  
  # convert run time for subtractions, to calculate durations between checkpoints
  runtime_numeric <- as.numeric(runtime)
  
  # We calculate the passage of time between checkpoints in the script, 
  # Step duration is the duration of a step between neighboring checkpoints
  # Time elapsed is the duration between the completion of a step and the 
  # start of the script.
  
  step_duration <- c(0) # No time has passed when the script is initialized.
  time_elapsed  <- c(0) # No time has passed when the script is initialized.
  
  for (i in seq(2, length(runtime_numeric), 1)) {
    
    step_duration <- c(step_duration, 
                       runtime_numeric[[i]] - runtime_numeric[[i-1]])
    
    time_elapsed  <- c(time_elapsed,
                       runtime_numeric[[i]] - runtime_numeric[[1]])
    
  }
  
  # Turn seconds into time periods and round for simplicity
  step_duration <- round(seconds_to_period(step_duration))
  time_elapsed  <- round(seconds_to_period(time_elapsed))
  
 
  # summarize all the runtime data in a tibble
  runtime <- runtime               %>%
    as_tibble_col()                %>%
    unnest(cols = c(value))        %>%
    rename("Start Time" = "value") %>% 
    add_column("Step Name" = runtime_labels, .before = 1) %>%
    add_column("Step Duration"   = step_duration)         %>%
    add_column("Time Elapsed"    = time_elapsed) 
  
  # remove uneccesary variables
  rm(runtime_numeric, 
     step_duration, 
     time_elapsed, 
     runtime_labels,
     i)
  
  
  
  
  } # end of subpoint
  
  # 11. Tidying R environment and saving to Outputs under custom name
  {
    
  # Automatically generate environment save file name
  env_save_path <- str_glue("Outputs/DilutionEnv_", 
                            run_mode,
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
                            ".RData")
  
  # Summarizing important but somewhat scattered meta data
  script_params <- list("dilution_props"  = dilution_props, 
                        "number_ranks"    = number_ranks, 
                        "runtime"         = runtime, 
                        "cellchat_nperms" = cellchat_nperms, 
                        "feature_type"    = feature_type, 
                        "methods_vector"  = methods_vector, 
                        "run_mode"        = run_mode, 
                        "save_results"    = save_results,
                        "testdata_type"   = testdata_type,
                        "save_names"      = list("plot_png_name" = plot_png_name,
                                                 "env_save_path" = env_save_path))
  # Removing now-superfluous meta data
  rm(dilution_props, number_ranks, runtime, cellchat_nperms, feature_type, 
     methods_vector, run_mode, save_results, testdata_type, plot_png_name, 
     env_save_path)
  
  # Save R environment and all the results within it
  if (script_params$save_results) {
    
    save.image(file = script_params$save_names$env_save_path)
    
    print(str_glue("Environment saved at ~/", 
                   script_params$save_names$env_save_path,
                   "."))
  
  }
  
  } # end of subpoint
  
}