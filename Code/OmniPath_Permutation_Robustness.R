#------------------------------------------------------------------------------#
# A. Setup ---------------------------------------------------------------------
{
  # 0.1 Overview of Goals:  ---
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
    
  # 0.3 Loading Packages
  {
    require(tidyverse)
    require(Seurat)
    require(liana)
  } # end of subpoint
  
}



#------------------------------------------------------------------------------#
# B. Set top_n, dilution props, testdata type ----------------------------------
{
  dilution_props <- c(seq(0.20, 0.8, 0.20)) # should be consistent between tests
  
  number_ranks   <- list("connectome" = 1000, 
                         "cellchat"   = 1000,
                         "italk"      = 1000,
                         "sca"        = 1000,
                         "natmi"      = 1000,
                         "squidpy"    = 1000)
  
  testdata_type  <- c("liana_test") # choose "liana_test" or "seurat_pbmc"
 
  feature_type <- c("variable") # choose "generic" or "variable"
  # If feature_type is generic, dilution will be completed with any genes
  # in the seurat count matrix. If dilution is variable, only the variable
  # features are used for dilution.
  
  
  
}   


  
#------------------------------------------------------------------------------#
# C. Preparing necessary inputs to dilute Resources ----------------------------
{
  runtime <- list("start" = Sys.time())
  
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
  # Omnipath x (cellchat, connectome, italk, sca) on testdata
  # no natmi because it doesn't work right and its being reimplemented
  # squidpy won't be used until I get it to work on windows
  liana_results_OP_0 <- liana_wrap(testdata,
                                   method   = c('connectome', 'cellchat',
                                                'italk', 'sca'),
                                   resource = c('OmniPath'))
  
  
  
  } # end of subpoint
  
  # 2. Get highest ranked interactions for undiluted conditions 
  {
    
    
  # Apply get_top_n_ranks for each method's results on OP_0 (i.e. undiluted)
  top_ranks_OP_0 <- list()
  
  for(method in c('connectome', 'cellchat', 'italk', 'sca')){
    
    top_ranks_OP_0[[method]] <- get_top_n_ranks(data_set = 
                                                liana_results_OP_0[[method]],
                                              top_n = number_ranks[[method]],
                                              method = method)
    
  }
    
    
    
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
  
  for(i in dilution_props) {
    dilution_names <- 
      c(dilution_names, str_glue("OmniPath_", as.character(i*100)))
  }
  
  dilution_props <- as.list(dilution_props)
  names(dilution_props) <- dilution_names
  
  rm(dilution_names, i)
  
  #relist all our data into three concise named lists of named lists
  resources_OP <- list("connectome" = list(OmniPath_0 = OmniPath_0),
                         "cellchat" = list(OmniPath_0 = OmniPath_0),
                         "italk"    = list(OmniPath_0 = OmniPath_0),
                        #"natmi"    = list(OmniPath_0 = OmniPath_0),
                         "sca"      = list(OmniPath_0 = OmniPath_0))
  
  
  
  liana_results_OP <- list("connectome" = 
                             list(OmniPath_0 = liana_results_OP_0$connectome),
                           "cellchat" = 
                             list(OmniPath_0 = liana_results_OP_0$cellchat),
                           "italk" = 
                             list(OmniPath_0 =  liana_results_OP_0$italk),
                           #"natmi"     = 
                           # list(OmniPath_0 = liana_results_OP_0$natmi),
                           "sca" = 
                             list(OmniPath_0 = liana_results_OP_0$sca))

  
  
  
  top_ranks_OP <- list("connectome" = list(OmniPath_0 = top_ranks_OP_0$connectome),
                       "cellchat"   = list(OmniPath_0 = top_ranks_OP_0$cellchat),
                       "italk"      = list(OmniPath_0 = top_ranks_OP_0$italk),
                       #"natmi"     = list(OmniPath_0 = top_ranks_OP_0$natmi),
                       "sca"        = list(OmniPath_0 = top_ranks_OP_0$sca))
  
  
  
  # remove old data frames
  rm(OmniPath_0, liana_results_OP_0, top_ranks_OP_0)
  
  
  
  } # end of subpoint
  
}



#------------------------------------------------------------------------------#
# D. Diluting Resources and comparing Top Ranks --------------------------------
{
  # 5. Generate diluted Resources for all methods
  {
    
    
  # Initiating a list of all dilutions
  dilutions_OP <- list()
  
  # Iterate over every method, lapply over every dilution proportion
  for(method in c('connectome', 'cellchat', 'italk', 'sca')){
    
    dilutions_OP[[method]] <- 
      lapply(dilution_props, dilute_Resource, 
             resource = resources_OP[[method]]$OmniPath_0, 
             top_rank_df = top_ranks_OP[[method]]$OmniPath_0, 
             data_set = testdata,
             dilution_type = feature_type)
    
  }
  
  
  
  # Merge OP_0 with the rest of the dilutions, could use mapply but its less 
  # consistent
  for (method in c('connectome', 'cellchat', 'italk', 'sca')) {
    for (dilution in names(dilution_props)) {
      
      resources_OP[[method]][[dilution]] <- dilutions_OP[[method]][[dilution]]
      
    }
  }
  
  
  # Remove uneccesary Variables
  rm(dilutions_OP, method, dilution)
    
  
    
  } # end of subpoint
  
  # 6. Reapply individual methods with diluted resources
  {
    
    
  # Initialize a list for liana results using diluted resources
  liana_dilutions_OP <- list()
    
  runtime[["methods_start"]] <- Sys.time()
    
  
  # Lapply call_x functions from liana for each method over every dilution.
  liana_dilutions_OP[["connectome"]] <- 
    lapply(resources_OP$connectome[-1], call_connectome, 
           seurat_object = testdata)
  
  runtime[["connectome_end"]] <- Sys.time()
  
  
  liana_dilutions_OP[["cellchat"]] <- 
    lapply(resources_OP$cellchat[-1], call_cellchat, 
           seurat_object = testdata, thresh = 1)
  
  runtime[["cellchat_end"]] <- Sys.time()
  
  
  liana_dilutions_OP[["italk"]] <- 
    lapply(resources_OP$italk[-1], call_italk, 
           seurat_object = testdata)
  
  runtime[["italk_end"]] <- Sys.time()
  
  # automatically iterates over list because of hurdles of conda env
  #liana_dilutions_OP[["natmi"]] <- 
  #  call_natmi(seurat_object = testdata, op_resource = resources_OP$natmi[-1]) 
  
  #runtime[["natmi_end"]] <- Sys.time()
  
  
  liana_dilutions_OP[["sca"]] <- 
    lapply(resources_OP$sca[-1], call_sca, 
           seurat_object = testdata, s.score = 0, logFC = 0.5849625) 
  
  runtime[["sca_end"]] <- Sys.time()
  # for some reason the number of results with sca still goes up, we're unsure why
  
  
  
  
  
  
  # Merge with undiluted results, could use mapply but its less consistent
  for (method in c('connectome', 'cellchat', 'italk', 'sca')) {
    for (dilution in names(dilution_props)) {
      
      liana_results_OP[[method]][[dilution]] <- 
        liana_dilutions_OP[[method]][[dilution]]
      
    }
  }

  
  # Remove uneccesary Variables
  rm(liana_dilutions_OP, method, dilution)

  
  
  } # end of subpoint
  
  # 7. Get top_n_ranks for each method and dilution
  {
    
    
  # lapply get_top_n_ranks over the dilution stages and save results in 
  # top_dilutions list
  top_dilutions_OP <- list()
  

  top_dilutions_OP[["connectome"]] <- 
    lapply(liana_results_OP$connectome[-1], get_top_n_ranks, 
           method = "connectome", top_n = number_ranks$connectome)
  
  
  top_dilutions_OP[["cellchat"]] <- 
    lapply(liana_results_OP$cellchat[-1], get_top_n_ranks, 
           method = "cellchat", top_n = number_ranks$cellchat)  
  
  
  top_dilutions_OP[["italk"]] <- 
    lapply(liana_results_OP$italk[-1], get_top_n_ranks, 
           method = "italk", top_n = number_ranks$italk)  
  
  # top_dilutions_OP[["natmi"]] <- 
  # lapply(liana_results_OP$natmi[-1], get_top_n_ranks, 
  #        method = "natmi", top_n = number_ranks$natmi) 
  
  top_dilutions_OP[["sca"]] <- 
    lapply(liana_results_OP$sca[-1], get_top_n_ranks, 
           method = "sca", top_n = number_ranks$sca)  
  
  
  
  
  # Merge with undiluted results, could use mapply but its less consistent
  for (method in c('connectome', 'cellchat', 'italk', 'sca')) {
    for (dilution in names(dilution_props)) {
      
      top_ranks_OP[[method]][[dilution]] <- top_dilutions_OP[[method]][[dilution]]
      
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
  top_ranks_OP$connectome <- 
      lapply(top_ranks_OP$connectome, unite, col = "LR_ID", 
             c(source, target, ligand, receptor), remove = FALSE)
    
    
  top_ranks_OP$cellchat <- 
    lapply(top_ranks_OP$cellchat, unite, col = "LR_ID", 
            c(source, target, ligand, receptor), remove = FALSE)
  
  
  top_ranks_OP$italk <- 
    lapply(top_ranks_OP$italk, unite, col = "LR_ID", 
           c(source, target, ligand, receptor), remove = FALSE)
  
  
  #top_ranks_OP$natmi <-
  # lapply(top_ranks_OP$natmi, unite, col = "LR_ID",
  #        c(source, target, ligand, receptor), remove = FALSE)
  
  
  top_ranks_OP$sca <- 
    lapply(top_ranks_OP$sca, unite, col = "LR_ID", 
           c(source, target, ligand, receptor), remove = FALSE)
  
  
  
  
  
  
  # add a column to see if an interaction is fake
  for (method in c('connectome', 'cellchat', 'italk', 'sca')) {
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
  overlap_connectome <- lapply(top_ranks_OP$connectome, rank_overlap, 
                               main_ranks = top_ranks_OP$connectome$OmniPath_0)
  
  
  overlap_cellchat   <- lapply(top_ranks_OP$cellchat, rank_overlap, 
                               main_ranks = top_ranks_OP$cellchat$OmniPath_0)
  
  
  overlap_italk      <- lapply(top_ranks_OP$italk, rank_overlap, 
                               main_ranks = top_ranks_OP$italk$OmniPath_0)
  
  
  overlap_sca        <- lapply(top_ranks_OP$sca, rank_overlap,
                               main_ranks = top_ranks_OP$sca$OmniPath_0)
  
  
  #overlap_natmi    <- lapply(top_ranks_OP$natmi, rank_overlap, 
  #                            main_ranks = top_ranks_OP$natmi$OmniPath_0)
  
  
  
  # reformatting overlap as a tibble
  top_rank_overlap <- tibble("connectome" = overlap_connectome,
                             "cellchat"   = overlap_cellchat,
                             "italk"      = overlap_italk,
                             "sca"        = overlap_sca) %>%
                      unnest(cols = c(connectome, cellchat, italk, sca)) %>%
                      mutate(dilution_prop = c(0, dilution_props)) %>%
                      unnest(cols = c(dilution_prop)) %>%
                      relocate("dilution_prop")
  
  # removing superfluous values
  rm(overlap_connectome, overlap_cellchat, overlap_italk, overlap_sca)
  
  } # end of subpoint
  
}



#------------------------------------------------------------------------------#
# E. Visualizing the results ---------------------------------------------------
  
  # The plot is better in percent than proportion
  top_rank_overlap_plot <- top_rank_overlap * 100
  
  # Automatically assemble a file name and plot subtitle
  plotting_subtitle <- str_glue(feature_type,
                                " dilution, top ",
                                as.character(median(unlist(number_ranks))),
                                " ranks, ",
                                testdata_type,
                                " data")
  
  plot_png_name <- str_glue(testdata_type, 
                            "_top",
                            as.character(median(unlist(number_ranks))),
                            "_",
                            feature_type,
                            "_",
                            as.character(Sys.Date()),
                            ".png")
  
  # Plot top_rank_overlap with lines and points at each value
  ggplot(data = top_rank_overlap_plot) + 
    geom_line(mapping = aes(dilution_prop, connectome, color =  "Connectome")) +
    geom_line(mapping = aes(dilution_prop, cellchat, color = "CellChat")) +
    geom_line(mapping = aes(dilution_prop, italk, color = "iTALK")) +
    geom_line(mapping = aes(dilution_prop, sca, color = "SCA")) +
    
    geom_point(mapping = aes(dilution_prop, connectome, color =  "Connectome")) +
    geom_point(mapping = aes(dilution_prop, cellchat, color = "CellChat")) +
    geom_point(mapping = aes(dilution_prop, italk, color = "iTALK")) +
    geom_point(mapping = aes(dilution_prop, sca, color = "SCA")) +
    
    # Show full breadth of 100-0 percent overlap
    ylim(0, 100) +
    
    ggtitle("Robustness of Method Predictions") +
    ylab("Overlap of Top Ranks [%]") +
    xlab("Dilution of Resource [%]") +
    labs(subtitle = plotting_subtitle,
         color = "Method")
  
  # Save the plot automatically to the outputs folder
  ggsave(plot_png_name, 
         height = 5, width = 8, 
         path = "Outputs")
  
  # Remove unnecessary variables
  rm(top_rank_overlap_plot)
  
#------------------------------------------------------------------------------#
# F. Saving the results --------------------------------------------------------
  
  
  
  runtime[["end"]] <- Sys.time()
  
  runtime_numeric <- as.numeric(runtime)
  
  seconds_elapsed <- c(0, 
                       runtime_numeric[[2]] - runtime_numeric [[1]],
                       runtime_numeric[[3]] - runtime_numeric [[2]],
                       runtime_numeric[[4]] - runtime_numeric [[3]],
                       runtime_numeric[[5]] - runtime_numeric [[4]],
                       runtime_numeric[[6]] - runtime_numeric [[5]],
                       runtime_numeric[[7]] - runtime_numeric [[6]])

  minutes_elapsed <- seconds_elapsed / 60
  hours_elapsed <- minutes_elapsed / 60
  
  seconds_elapsed <- round(seconds_elapsed, 2)
  minutes_elapsed <- round(minutes_elapsed, 2)
  hours_elapsed <- round(hours_elapsed, 2)
  
  runtime <- runtime %>%
    as_tibble_col() %>%
    unnest(cols = c(value)) %>%
    add_column(seconds_elapsed) %>%
    add_column(minutes_elapsed) %>%
    add_column(hours_elapsed)
  
  rm(runtime_numeric, seconds_elapsed, minutes_elapsed, hours_elapsed)
  
  env_save_path <- str_glue("Outputs/DilutionEnv_", 
                            testdata_type, 
                            "_top",
                            as.character(median(unlist(number_ranks))),
                            "_",
                            feature_type,                            
                            "_",
                            as.character(Sys.Date()),
                            ".RData")
  
  save.image(file = env_save_path)
  
