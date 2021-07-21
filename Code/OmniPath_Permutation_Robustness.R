# 0. Setup
{
  # 0.1 Overview of Goals: 
  {
  # The Idea with this script is to:
  # .	run all methods on a single resource (OP)
  # .	take the topmost ranked interactions for each method -> R-zero
  # .	create a modified OP resource in which the topmost ranked interactions
  #     remain, but of the remainder x% of the interactions have been removed and 
  #     replaced with entirely random pairs of genes derived from the test data 
  #     that do not exist in the resource, x =10,20,40% etc.
  # .	Rerun methods on modified omnipath resource, get top ranks -> R-modified
  # .	plot percentage of R-zero in R-modified over x and investigate result
  }
  
  # 0.2 Script Structure:
  {
    # We load preprocessed seurat data and apply each method combined with undiluted OmniPath to it.
    # We define a function that can get the top ranked CCI from each method output and run it for the undiluted results.
    # We restructure OmniPath, define a function that dilutes Resources and restructure our data to be used at various dilutions.
  
    # We have base OmniPath and the top ranks for each method. These basic inputs allow us to run dilute_Resource(). We create diluted resources for each method.
    # We get the top ranks foreach method at each new dilution.
    # We compare the percentage overlap of the undiluted top_ranks to the top_ranks at each stage of dilution (for each method).
  
    # We visualize the results.
  }
    
  # 0.2 Loading Packages
  {
    require(tidyverse)
    require(Seurat)
    require(liana)
    
    load(RobustRankAggreg)
  }
}

#------------------------------------------------------------------------------#
# A. Preparing necessary inputs to dilute Resources
{
  # 1. Running LIANA wrapper
    # 1.1 Get test data
    {
        liana_path <- system.file(package = 'liana')                                    # get liana package filepath
        testdata <- 
          readRDS(file.path(liana_path, "testdata", "input", "testdata.rds"))           # retrieve testdata from filepath
        
        rm(liana_path)
}
    
    # 1.2 Generate Undiluted liana results
    { 
    # Run wrapper on testdata for omnipath x (cellchat, connectome, italk, sca) no natmi because it doesn't work right and its being reimplemented
    #squidpy won't be used until I get it to work on windows
    liana_results_OP_0 <- liana_wrap(testdata,
                               method = c('connectome', 'cellchat', 'italk', 'sca'),
                               resource = c('OmniPath'))
  }
  
  # 2. Get highest ranked interactions for undiluted conditions 
    # 2.2 Apply get_top_n_ranks for each method's results on OP_0 (i.e. undiluted)
    {
    top_ranks_OP_0 <- list("connectome" = get_top_n_ranks(dat = liana_results_OP_0$connectome, top_n = 200, met = "connectome"),
                          "cellchat" = get_top_n_ranks(dat = liana_results_OP_0$cellchat, top_n = 30, met = "cellchat"), ## doesn't produce that many more interactions
                          "italk" = get_top_n_ranks(dat = liana_results_OP_0$italk, top_n = 100, met = "italk"), ## doesn't produce that many more interactions
                        #  "natmi" = get_top_n_ranks(dat = liana_results_OP_0$natmi, top_n = 200, met = "natmi"),
                          "sca" = get_top_n_ranks(dat = liana_results_OP_0$sca, top_n = 80, met = "sca")) ## produces variable result nrows at higher dilutions, some of them are low
    }
  
  
  # 3. Prepare to Modify OmniPath with random genes
    # 3.1 Format OmniPath_0 to be easier to work with and to pair it down to the columns relevant for the methods
    {
    OmniPath_0 <- select_resource(c('OmniPath'))
    OmniPath_0 <- OmniPath_0[["OmniPath"]]
    
    OmniPath_0 <- select(OmniPath_0, 
      source_genesymbol,
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
      entity_type_intercell_target)
    
    OmniPath_0 <- mutate(OmniPath_0,
                         isRandom = FALSE)
    
    OmniPath_0 <- unite(OmniPath_0, "LR_Pair", c(source_genesymbol, target_genesymbol), remove = FALSE, sep = "_")
    OmniPath_0 <-  relocate(OmniPath_0, "LR_Pair", .after = last_col())
}
    
  # 4. Listify Data sets in the face of dilution steps
  {
  # Since we are about to perform the same analysis in the steps above but multiplied by each dilution step, we will turn our data sets into named lists sorted by method and sub categorized by dilution.
  
  # define dilution proportions
  dilution_props <-c(seq(0.2, 0.6, 0.2))
  dilution_names <- c()
  
  for(i in dilution_props) {
    dilution_names <- c(dilution_names, str_glue("OmniPath_", as.character(i*100)))
  }
  
  dilution_props <- as.list(dilution_props)
  names(dilution_props) <- dilution_names
  
  rm(dilution_names, i)
  
  #relist all our data into three concise named lists of named lists
  resources_OP <- list("connectome" = list(OmniPath_0 = OmniPath_0),
                         "cellchat" = list(OmniPath_0 = OmniPath_0),
                         "italk" = list(OmniPath_0 = OmniPath_0),
                     #    "natmi" = list(OmniPath_0 = OmniPath_0),
                         "sca" = list(OmniPath_0 = OmniPath_0))
    
  liana_results_OP <- list("connectome" = list(OmniPath_0 = liana_results_OP_0$connectome),
                           "cellchat" = list(OmniPath_0 = liana_results_OP_0$cellchat),
                           "italk" = list(OmniPath_0 = liana_results_OP_0$italk),
                      #     "natmi" = list(OmniPath_0 = liana_results_OP_0$natmi),
                           "sca" = list(OmniPath_0 = liana_results_OP_0$sca))
  
  top_ranks_OP <- list("connectome" = list(OmniPath_0 = top_ranks_OP_0$connectome),
                       "cellchat" = list(OmniPath_0 = top_ranks_OP_0$cellchat),
                       "italk" = list(OmniPath_0 = top_ranks_OP_0$italk),
                   #    "natmi" = list(OmniPath_0 = top_ranks_OP_0$natmi),
                       "sca" = list(OmniPath_0 = top_ranks_OP_0$sca))
  
  # remove old data frames
  rm(OmniPath_0, liana_results_OP_0, top_ranks_OP_0)
}
}


#------------------------------------------------------------------------------#
# B. Diluting Resources and comparing Top Ranks
  
  # 5. Generate diluted Resources for all methods
  {
    
  # Initiating a list of all dilutions
  dilutions_OP <- list()
  
  # Iterate over every method, lapply over every dilution proportion
  for(i in c('connectome', 'cellchat', 'italk', 'sca')){
  dilutions_OP[[i]] <- lapply(dilution_props, 
                       dilute_Resource, resource = resources_OP[[i]]$OmniPath_0, 
                                       top_rank_df = top_ranks_OP[[i]]$OmniPath_0, 
                                       data_set = testdata)
  }
  
  # Merge OP_0 with the rest of the dilutions
  resources_OP <- mapply(c, resources_OP, dilutions_OP, SIMPLIFY=FALSE)
  
  # Remove uneccesary Variables
  rm(dilutions_OP, i)
      
  }
  
  # 6. Reapply individual call functions from Liana  with diluted resources and store top ranks
  {
 
     
  # Initialize a list for liana results using diluted resources
  liana_dilutions_OP <- list()
  
  # Lapply call functions from liana over every dilution. Different method every line
  liana_dilutions_OP[["connectome"]] <- lapply(resources_OP$connectome[-1], call_connectome, seurat_object = testdata)
  liana_dilutions_OP[["cellchat"]] <- lapply(resources_OP$cellchat[-1], call_cellchat, seurat_object = testdata)
  liana_dilutions_OP[["italk"]] <- lapply(resources_OP$italk[-1], call_italk, seurat_object = testdata)
  #liana_dilutions_OP[["natmi"]] <- call_natmi(seurat_object = testdata, op_resource = resources_OP$natmi[-1]) # automatically iterates over list because of hurdles of conda env
  liana_dilutions_OP[["sca"]] <- lapply(resources_OP$sca[-1], call_sca, seurat_object = testdata) 
  
  # Merge with undiluted results, could use mapply but its less consistent
  for (i in c('connectome', 'cellchat', 'italk', 'sca')) {
    for (j in names(dilution_props)) {
      liana_results_OP[[i]][[j]] <- liana_dilutions_OP[[i]][[j]]
    }
  }

  
  # Remove uneccesary Variables
  rm(liana_dilutions_OP)

  
  } 
  
  #7. Get top_n_ranks and evaluate ranks

  top_dilutions_OP <- list()
  

  top_dilutions_OP[["connectome"]] <- lapply(liana_results_OP$connectome[-1], 
                                get_top_n_ranks, met = "connectome", top_n = 200)
  top_dilutions_OP[["cellchat"]] <- lapply(liana_results_OP$cellchat[-1], 
                                              get_top_n_ranks, met = "cellchat", top_n = 30)  
  top_dilutions_OP[["italk"]] <- lapply(liana_results_OP$italk[-1], 
                                             get_top_n_ranks, met = "italk", top_n = 100)  
  # top_dilutions_OP[["natmi"]] <- lapply(liana_results_OP$natmi[-1], 
  #                                            get_top_n_ranks, met = "natmi", top_n = 200)  
  top_dilutions_OP[["sca"]] <- lapply(liana_results_OP$sca[-1], 
                                             get_top_n_ranks, met = "sca", top_n = 200)  
  
  # Merge with undiluted results, could use mapply but its less consistent
  for (i in c('connectome', 'cellchat', 'italk', 'sca')) {
    for (j in names(dilution_props)) {
      top_ranks_OP[[i]][[j]] <- top_dilutions_OP[[i]][[j]]
    }
  }
  
  rm(top_dilutions_OP)
  
  
  #8. Evaluate how many of the top 200 interactions overlap
  
  # format top_ranks to have an ID that marks each specific interaction
  top_ranks_OP$connectome <- lapply(top_ranks_OP$connectome, unite, col = "LR_ID", c(source, target, ligand, receptor), remove = FALSE)
  top_ranks_OP$cellchat <- lapply(top_ranks_OP$cellchat, unite, col = "LR_ID", c(source, target, ligand, receptor), remove = FALSE)
  top_ranks_OP$italk <- lapply(top_ranks_OP$italk, unite, col = "LR_ID", c(source, target, ligand, receptor), remove = FALSE)
  #top_ranks_OP$natmi <- lapply(top_ranks_OP$natmi, unite, col = "LR_ID", c(source, target, ligand, receptor), remove = FALSE)
  top_ranks_OP$sca <- lapply(top_ranks_OP$sca, unite, col = "LR_ID", c(source, target, ligand, receptor), remove = FALSE)
  
  
  # add a column to see if an interaction is fake
  
  for (i in c('connectome', 'cellchat', 'italk', 'sca')) {
    for (j in c("OmniPath_0", names(dilution_props))) {
      if( !(is_null(top_ranks_OP[[i]][[j]]))) {
      top_ranks_OP[[i]][[j]] <- mutate(top_ranks_OP[[i]][[j]], isRandom = 
                                         !(top_ranks_OP[[i]][[j]]$LR_Pair %in% resources_OP[[i]]$OmniPath_0$LR_Pair))
      }
    }
  }
  
  
  rm(i,j)
  
  
  
  
  
