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
  
  # 0.2 Loading Packages
  {
    require(tidyverse)
    require(Seurat)
    require(liana)
    
    load(RobustRankAggreg)
  }
}

# 1. Running LIANA wrapper

  # 1.1 Get test data
  {
      liana_path <- system.file(package = 'liana')                                    # get liana package filepath
      testdata <- 
        readRDS(file.path(liana_path, "testdata", "input", "testdata.rds"))           # retrieve testdata from filepath
}
  
  # 1.2 Run wrapper on testdata for omnipath x connectome
  { 
  # In future it should be OmniPath x (cellchat, connectome, italk, natmi, sca), squidpy won't be used until I get it to work on windows
  liana_results_OP_0 <- liana_wrap(testdata,
                             method = c('connectome'),
                             resource = c('OmniPath'))
  }

# 2. Extract top x ranked interactions for a given method  
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
  
   return(topranks)
  }
  
  # Apply get_top_n_ranks for each method's results on OP_0. Other emthods currently commented out
  top_ranks_OP_0 <- list("connectome" = get_top_n_ranks(dat = liana_results_OP_0, top_n = 200, met = "connectome")
                       #  ,
                       #  "cellchat" = get_top_n_ranks(dat = liana_results_OP_0, top_n = 30, met = "cellchat"),
                       #  "italk" = get_top_n_ranks(dat = liana_results_OP_0, top_n = 200, met = "italk"),
                       #  "natmi" = get_top_n_ranks(dat = liana_results_OP_0, top_n = 200, met = "natmi"),
                       #  "sca" = get_top_n_ranks(dat = liana_results_OP_0, top_n = 200, met = "sca"))
  )
  
  # Format the top_ranks data frames for future processing steps. Other methods would need to be added here.
  top_ranks_OP_0$connectome <- unite(top_ranks_OP_0$connectome, "LR_Pair", 3:4, remove = FALSE, sep = "_")
  top_ranks_OP_0$connectome <-  relocate(top_ranks_OP_0$connectome, "LR_Pair", .after = last_col())
}

# 3. Modifying OmniPath with random genes

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
  #' @param top_rank_list As a list of tibbles. Each tibble in the list is named after the method that generated it and contains the top ranked rows of the liana_wrapper output. Since we are comparing how resource dilution affects top ranked interactions, this list of top rankings using the undiluted resource ensures that the top interactions can still be caught in the same way. It's the effect of surrounding diluted noise that we'll be picking up on.
  #' @param dilution_prop As a number between 0-1. The proportion of rows of the resource to dilute. Top ranked rows can't be diluted. If attaining the requested dilution proportion requires overwriting top ranked interactions, the function throws an error and returns nothing instead.
  #' @param used_method As a string. The method used to generae the top ranked list.
  #' 
  #' @return Returns a tibble that can be used as a resource for the liana wrapper function but has a certain (marked) percentage of it replaced with random nonsensical interactions.

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
    # Additional Warning message and break if the dilution prop can't be met
    dilution_number <- round(nrow(resource)*dilution_prop)
    
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
  
  # 3.3 Apply dilute_Resource for connectome at various dilution stages
  {
    
    diluted_resources_OP <- list()
    dilution_props <- as.list(seq(0.1, 0.6, 0.1))
    
    for(i in dilution_props) {
      diluted_resources_OP[[str_glue("OmniPath_", as.character(i*100))]] <- dilute_Resource(resource = OmniPath_0, 
                                                                           top_rank_list = liana_results_OP_0, 
                                                                           dilution_prop = i, 
                                                                           data_set = testdata, 
                                                                           used_method = "connectome")
    }
    
  }

# 4. Reapply individual call functions from Liana  with diluted resources and store top ranks

liana_results_OP_diluted <- list()

for(i in seq(1, length(diluted_resources_OP), 1)) {
  dil_resource_name <- names(diluted_resources_OP[i])
  liana_results_OP_diluted[[dil_resource_name]] <- list("connectome" = call_connectome(seurat_object = testdata, op_resource = diluted_resources_OP[[i]]))
                                                    #,
                                                    # "cellchat" = call_cellchat(seurat_object = testdata, op_resource = diluted_resources_OP[[i]]),
                                                    # "italk" = call_italk(seurat_object = testdata, op_resource = diluted_resources_OP[[i]]),
                                                    # "natmi"= call_natmi(seurat_object = testdata, op_resource = diluted_resources_OP[[i]]),
                                                    # "sca"= call_sca(seurat_object = testdata, op_resource = diluted_resources_OP[[i]]))
                                                    
                                                
  rm(dil_resource_name)
}

top_ranks_connectome_OP <- list(OmniPath_0 = top_ranks_OP_0$connectome)

for(i in seq(1, length(diluted_resources_OP), 1)) {
  dil_resource_name <- names(diluted_resources_OP[i])
  top_ranks_connectome_OP[[dil_resource_name]] <- get_top_n_ranks(dat = liana_results_OP_diluted[[i]], met = "connectome", top_n = 200)
  rm(dil_resource_name)
}




get_top_n_ranks(dat = liana_results_OP_diluted[[i]], met = "connectome", top_n = 200)

top_ranks_OP_0$connectome <- unite(top_ranks_OP_0$connectome, "LR_Pair", 3:4, remove = FALSE, sep = "_")
top_ranks_OP_0$connectome <-  relocate(top_ranks_OP_0$connectome, "LR_Pair", .after = last_col())




















