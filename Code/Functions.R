#------------------------------------------------------------------------------#
# A. For OmniPath_Permutation_Robustness
  # get_top_n_ranks()
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
        unite("LR_Pair", c(ligand, receptor), remove = FALSE, sep = "_") %>%
        relocate("LR_Pair", .after = last_col())
      
      return(topranks)
    }
  }

  # dilute_Resource()
  {
    
    
    
  #' Dilutes a resource with randomly generated interactions derived from specific genes
  #'
  #' @param resource The resource (as a tibble) which you would like to falsify / dilute with random gene relationships.
  #' @param data_set The data set (as a Seurat object) from which you will draw genes. Diluting the resource with genes that are actually found in the data set the resource will be used with guarantees that the random relationships will be relevant to the method. This function excludes dilution with genes already present in the base resource and by extension genes in the top ranking.
  #' @param top_rank_df As a tibble. Ideally the output of get_top_n_ranks. Should contain the top ranked rows of the liana_wrapper output for the method you are diluting for (using the undiluted resource). Since we are comparing how resource dilution affects top ranked interactions, this list of top rankings using the undiluted resource ensures that the top interactions can still be caught in the same way. It's the effect of surrounding diluted noise that we'll be picking up on.
  #' @param dilution_prop As a number between 0-1. The proportion of rows of the resource to dilute. Top ranked rows can't be diluted. If attaining the requested dilution proportion requires overwriting top ranked interactions, the function throws an error and returns nothing instead.
  #' 
  #' @return Returns a tibble that can be used as a resource for the liana call_method functions but has a certain (marked) percentage of it replaced with random nonsensical interactions.
  
  # Format the dilute_Resource Function to make diluted OP resources for each method
  dilute_Resource <- function(resource_to_dil, top_rank_df, dilution_prop, data_set){
    
    # Generate a list of gene names that will be relevant by getting them from the data_set
    gene_name_list  <- as.list(rownames(data_set@assays$RNA@data))
    
    # Remove items already in OmniPath to ensure nonsense relationships, and to not mess with the top_ranks
    gene_name_list <- gene_name_list[!(gene_name_list %in% resource_to_dil$source_genesymbol)]
    gene_name_list <- gene_name_list[!(gene_name_list %in% resource_to_dil$target_genesymbol)]
    
    # Separate top_rank parts of resource from non top rank parts of resource, we only want to dilute the latter
    resource_top <- resource_to_dil %>% 
      filter(LR_Pair %in% top_rank_df$LR_Pair)
    
    resource_bottom <- resource_to_dil %>% 
      filter(!(LR_Pair %in% top_rank_df$LR_Pair))
    
    # Determine how many rows of resource_bottom to dilute so that the overall dilution_prop is met
    dilution_number <- round(nrow(resource_to_dil)*dilution_prop)
    
    # Additional Warning message and break if the dilution prop can't be met
    if(dilution_number > nrow(resource_bottom)) {
      warning("Requested dilution proportion not attainable without overwriting 
                the given top-ranked interacions. Returning nothing instead.")
      return()
    }
    
    # Select the candidates for dilution from resource_bottom
    set.seed(90)
    resource_dilute <- slice_sample(resource_bottom, n = dilution_number)
    # Delete the dilution candidates from resource_bottom so we have an even split into dilution candidates and non-diluted candidates
    resource_bottom <- anti_join(resource_bottom, resource_dilute)
    
    
    # Sample random gene names from the gene list into the dilution candidates, creating random nonsensical relationships
    set.seed(91)
    resource_dilute$source_genesymbol <- as.character(sample(gene_name_list, size = nrow(resource_dilute), replace = TRUE))
    set.seed(92)
    resource_dilute$target_genesymbol <- as.character(sample(gene_name_list, size = nrow(resource_dilute), replace = TRUE))
    
    # Update the LR-Pair column with the new random "interaction partners", and mark the random interactions as such.
    resource_dilute <- select(resource_dilute, -LR_Pair)
    resource_dilute <- unite(resource_dilute, "LR_Pair", c(source_genesymbol, target_genesymbol), remove = FALSE, sep = "_")
    resource_dilute <- mutate(resource_dilute, isRandom = TRUE)
    resource_dilute <-  relocate(resource_dilute, "LR_Pair", .after = last_col())
    
    # The new resource has top ranked interactions, non-top rank but still real interactions, and diluted random interactions.
    new_resource <- bind_rows(resource_top, resource_bottom, resource_dilute)

    #return output
    return(new_resource)
  }
    
    
    
  }






