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

--------------------------------------------------------------------------------
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
    # Run wrapper on testdata for omnipath x (cellchat, connectome, italk, natmi, sca)
    #squidpy won't be used until I get it to work on windows
    liana_results_OP_0 <- liana_wrap(testdata,
                               method = c('connectome', 'cellchat', 'italk', 'natmi', 'sca'),
                               resource = c('OmniPath'))
  }
  
  
  # 2. Get highest ranked interactions forundiluted conditionsDefine get_top_n_ranks  
    #2.1 Define get_top_n_ranks
    {
    #' Get the top n ranked items of a method from the liana wrapper results
    #'
    #' @param dat The list of tibbles output by the liana wrapper function.
    #' @param met The method for which you would like to extract the top ranked interactions, as a string.
    #' @param top_n The number of items to return, returns items ranked #1 to #n.
    #'
    #' @return Returns a tibble with all the columns for this method in the liana wrapper output but only including the rows of the top n interactions (tied values at the boundary line are cut off, no ties).
  
    get_top_n_ranks <- function(dat, met, top_n) {
      
      # generate a list that describes how to rank in each method to get what the method considers best
      rank_spec_list <- list("cellchat" = list(met_score = "pval",descending_order =  FALSE),
                             "connectome" = list(met_score = "weight_sc", descending_order =  TRUE),
                             "italk" = list(met_score = "weight_comb", descending_order =  TRUE),
                             "natmi" = list(met_score = "edge_specificity", descending_order =  TRUE),
                             "sca" = list(met_score = "LRscore", descending_order =  TRUE),
                             "squidpy" = list(met_score = "pvalue", descending_order =  FALSE))
      
      # If its a p-value method, the code will go into this if statement
      if(rank_spec_list[[met]]$descending_order == FALSE) {
         topranks <- slice_min(dat[[met]], 
                               n = top_n, 
                               order_by = !!sym(rank_spec_list[[met]]$met_score), # order by the criterion dictated by rank_spec_list
                               with_ties = FALSE)
                                       
       } else { # if it's not one of the p-value methods
         topranks <- slice_max(dat[[met]], 
                               n = top_n, 
                               order_by = !!sym(rank_spec_list[[met]]$met_score), # order by the criterion dictated by rank_spec_list
                               with_ties = FALSE)
      }
      
      # Format the top_ranks data frame for future processing steps.
      topranks <- topranks %>% 
        unite("LR_Pair", 3:4, remove = FALSE, sep = "_") %>%
        relocate("LR_Pair", .after = last_col())
      
      return(topranks)
    }
  }
    
    # 2.2 Apply get_top_n_ranks for each method's results on OP_0 (i.e. undiluted)
    {
    top_ranks_OP_0 <- list("connectome" = get_top_n_ranks(dat = liana_results_OP_0, top_n = 200, met = "connectome"),
                          "cellchat" = get_top_n_ranks(dat = liana_results_OP_0, top_n = 30, met = "cellchat"),
                          "italk" = get_top_n_ranks(dat = liana_results_OP_0, top_n = 200, met = "italk"),
                          "natmi" = get_top_n_ranks(dat = liana_results_OP_0, top_n = 200, met = "natmi"),
                          "sca" = get_top_n_ranks(dat = liana_results_OP_0, top_n = 200, met = "sca"))
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
      entity_type_intercell_target
    )
    
    OmniPath_0 <- mutate(OmniPath_0,
                         isRandom = FALSE)
    
    OmniPath_0 <- unite(OmniPath_0, "LR_Pair", 1:2, remove = FALSE, sep = "_")
    OmniPath_0 <-  relocate(OmniPath_0, "LR_Pair", .after = last_col())
}
    
    # 3.2 Format the dilute_Resource Function to make diluted OP resources for each method
    {
    #' Dilutes a resource with randomly generated interactions derived from specific genes
    #'
    #' @param resource The resource (as a tibble) which you would like to falsify / dilute with random gene relationships.
    #' @param data_set The data set (as a Seurat object) from which you will draw genes. Diluting the resource with genes that are actually found in the data set the resource will be used with guarantees that the random relationships will be relevant to the method. This function excludes dilution with genes already present in the base resource and by extension genes in the top ranking.
    #' @param top_rank_list As a list of tibbles. Each tibble in the list is named after the method that generated it and contains the top ranked rows of the liana_wrapper output (using the undiluted resource). Since we are comparing how resource dilution affects top ranked interactions, this list of top rankings using the undiluted resource ensures that the top interactions can still be caught in the same way. It's the effect of surrounding diluted noise that we'll be picking up on.
    #' @param dilution_prop As a number between 0-1. The proportion of rows of the resource to dilute. Top ranked rows can't be diluted. If attaining the requested dilution proportion requires overwriting top ranked interactions, the function throws an error and returns nothing instead.
    #' @param used_method As a string. The method used to generate the top ranked list.
    #' 
    #' @return Returns a tibble that can be used as a resource for the liana call_method functions but has a certain (marked) percentage of it replaced with random nonsensical interactions.
  
    dilute_Resource <- function(resource, top_rank_list, dilution_prop, data_set, used_method){
    
      # Generate a list of gene names that will be relevant by getting them from the data_set
      gene_name_list  <- as.list(testdata@assays$RNA@var.features)
      
      # Remove items already in OmniPath to ensure nonsense relationships, and to not mess with the top_ranks
      gene_name_list <- gene_name_list[!(gene_name_list %in% resource$source_genesymbol)]
      gene_name_list <- gene_name_list[!(gene_name_list %in% resource$target_genesymbol)]
      
      # Separate top_rank parts of resource from non top rank parts of resource, we only want to dilute the latter
      resource_top <- resource %>% 
        filter(LR_Pair %in% top_rank_list[used_method]$LR_Pair)
      
      resource_bottom <- resource %>% 
        filter(!(LR_Pair %in% top_rank_list[used_method]$LR_Pair))
      
      # Determine how many rows of resource_bottom to dilute so that the overall dilution_prop is met
      dilution_number <- round(nrow(resource)*dilution_prop)
      
      # Additional Warning message and break if the dilution prop can't be met
      if(dilution_number > nrow(resource_bottom)) {
        warning("Requested dilution proportion not attainable without overwriting 
                the given top-ranked interacions. Returning nothing instead.")
        return()
      }
      
      # Select the candidates for dilution from resource_bottom
      set.seed(123)
      resource_dilute <- slice_sample(resource_bottom, n = dilution_number)
      # Delete the dilution candidates from resource_bottom so we have an even split into dilution candidates and non-diluted candidates
      resource_bottom <- anti_join(resource_bottom, resource_dilute)
      
      
      # Sample random gene names from the gene list into the dilution candidates, creating random nonsensical relationships
      set.seed(123)
      resource_dilute$source_genesymbol <- as.character(sample(gene_name_list, size = nrow(resource_dilute), replace = TRUE))
      set.seed(123)
      resource_dilute$target_genesymbol <- as.character(sample(gene_name_list, size = nrow(resource_dilute), replace = TRUE))
      
      # Update the LR-Pair column with the new random "interaction partners", and mark the random interactions as such.
      resource_dilute <- select(resource_dilute, -LR_Pair)
      resource_dilute <- unite(resource_dilute, "LR_Pair", 1:2, remove = FALSE, sep = "_")
      resource_dilute <- mutate(resource_dilute, isRandom = TRUE)
      resource_dilute <-  relocate(resource_dilute, "LR_Pair", .after = last_col())
      
      # The new resource has top ranked interactions, non-top rank but still real interactions, and diluted random interactions.
      new_resource <- bind_rows(resource_top, resource_bottom, resource_dilute)
      
      return(new_resource)
    }
}
    
  
  # 4. Listify Data sets in the face of dilution steps
  {
  # Since we are about to perform the same analysis in the steps above but multiplied by each dilution step, we will turn our data sets into named lists sorted by method and sub categorized by dilution.
  
  # define dilution proportions
  dilution_props <-c(seq(0.1, 0.6, 0.1))
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
                         "natmi" = list(OmniPath_0 = OmniPath_0),
                         "sca" = list(OmniPath_0 = OmniPath_0))
    
  liana_results_OP <- list("connectome" = list(OmniPath_0 = liana_results_OP_0$connectome),
                           "cellchat" = list(OmniPath_0 = liana_results_OP_0$cellchat),
                           "italk" = list(OmniPath_0 = liana_results_OP_0$italk),
                           "natmi" = list(OmniPath_0 = liana_results_OP_0$natmi),
                           "sca" = list(OmniPath_0 = liana_results_OP_0$sca))
  
  top_ranks_OP <- list("connectome" = list(OmniPath_0 = top_ranks_OP_0$connectome),
                       "cellchat" = list(OmniPath_0 = top_ranks_OP_0$cellchat),
                       "italk" = list(OmniPath_0 = top_ranks_OP_0$italk),
                       "natmi" = list(OmniPath_0 = top_ranks_OP_0$natmi),
                       "sca" = list(OmniPath_0 = top_ranks_OP_0$sca))
  
  # remove old data frames
  rm(OmniPath_0, liana_results_OP_0, top_ranks_OP_0)
}
}


--------------------------------------------------------------------------------
# B. Diluting Resources and comparing Top Ranks
{
  # 5. Generate diluted Resources for all methods
  {
    # Initiating a list of all dilutions
    dilutions_OP <- list()
    
    # Iterate over every method, lapply over every dilution proportion
    for(i in c('connectome', 'cellchat', 'italk', 'natmi', 'sca')){
    dilutions_OP[[i]] <- lapply(dilution_props, 
                         dilute_Resource, resource = resources_OP[[i]]$OmniPath_0, 
                                         top_rank_list = top_ranks_OP[[i]], 
                                         data_set = testdata,
                                         used_method = i)
    }
    
    # Merge OP_0 with the rest of the dilutions
    resources_OP <- mapply(c, resources_OP, dilutions_OP, SIMPLIFY=FALSE)
    
    # Remove uneccesary Variables
    rm(dilutions_OP, i)
      
  }
  
  # 6. Reapply individual call functions from Liana  with diluted resources and store top ranks
  
  for(i in seq(1, length(dilution_props), 1)) {
    dil_resource_name <- names(OmniPath_dilutions$connectome[i+1])
    liana_results_OP[[dil_resource_name]] <- list("connectome" = call_connectome(seurat_object = testdata, op_resource = OmniPath_dilutions$connectome[[i+1]]))
                                                      #,
                                                      # "cellchat" = call_cellchat(seurat_object = testdata, op_resource = OmniPath_dilutions$cellchat[[i]]),
                                                      # "italk" = call_italk(seurat_object = testdata, op_resource = OmniPath_dilutions$italk[[i]]),
                                                      # "natmi"= call_natmi(seurat_object = testdata, op_resource = OmniPath_dilutions$natmi[[i]]),
                                                      # "sca"= call_sca(seurat_object = testdata, op_resource = OmniPath_dilutions$sca[[i]]))
                                                      
                                                  
    rm(dil_resource_name)
  }
  
  rm(i, liana_results_OP_0)
  
  ## repeat this code for every method
  top_ranks_connectome_OP <- list(OmniPath_0 = top_ranks_OP_0$connectome)
  
  for(i in seq(1, length(dilution_props), 1)) {
    dil_resource_name <- names(OmniPath_dilutions$connectome[i+1])
    top_ranks_list <- get_top_n_ranks(dat = liana_results_OP[[i+1]], met = "connectome", top_n = 200)
    
    top_ranks_list <- unite(top_ranks_list, "LR_Pair", 3:4, remove = FALSE, sep = "_")
    top_ranks_list <-  relocate(top_ranks_list, "LR_Pair", .after = last_col())
    
    top_ranks_connectome_OP[[dil_resource_name]] <- top_ranks_list
    rm(dil_resource_name, top_ranks_list)
  }
  
  top_ranks_OP <- list("connectome" = top_ranks_connectome_OP)
  
  rm(i, top_ranks_connectome_OP)
  
  ## delete uneccesary values
  rm(top_ranks_OP_0)
  
  ## restructure liana_results_OP, REPEAT for every method
  liana_results_OP_2 <- list("connectome" = list())
  
  for(i in seq(1, length(dilution_props)+1, 1)) {
    liana_results_OP_2$connectome[[names(liana_results_OP[i])]] <- liana_results_OP[[i]]$connectome
  }
  
  liana_results_OP <- liana_results_OP_2
  rm(i, liana_results_OP_2)

}
