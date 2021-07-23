#------------------------------------------------------------------------------#
# A. For OmniPath_Permutation_Robustness
  # get_top_n_ranks()
  {
    
    
  #' Get the top n ranked items of a method from the tibble liana wrapper or 
  #' call_x results
  #'
  #' @param data_set The tibble output by the liana wrapper function or call_x 
  #' functions.
  #' @param method The method for which you would like to extract the top ranked
  #' interactions, as a string.
  #' @param top_n The number of items to return, returns items ranked #1 to #n.
  #'
  #' @return Returns the tibble input cut down to the top n highest ranked 
  #' interactions.
  
  get_top_n_ranks <- function(data_set, method, top_n) {
    
    # generate a list that describes how to rank in each method to get what the 
    # method considers best
    rank_spec_list <- 
      list("cellchat"   = list(method_score = "pval",
                               descending_order =  FALSE),
           
           "connectome" = list(method_score = "weight_sc",
                               descending_order =  TRUE),
           
           "italk"      = list(method_score = "weight_comb",
                               descending_order =  TRUE),
           
           "natmi"      = list(method_score = "edge_specificity",
                               descending_order =  TRUE),
           
           "sca"        = list(method_score = "LRscore",
                               descending_order =  TRUE),
           
           "squidpy"    = list(method_score = "pvalue",
                               descending_order =  FALSE))
    
    
    
    # If its a p-value method, the code will go into this if statement
    if(rank_spec_list[[method]]$descending_order == FALSE) {
    
        topranks <- 
          slice_min(data_set, 
                    n = top_n, 
                    order_by = !!sym(rank_spec_list[[method]]$method_score), 
                    with_ties = FALSE)
        
        # order_by is assigned the criterion dictated by rank_spec_list
        
        
    } else { # if it's not one of the p-value methods
    
        topranks <- 
          slice_max(data_set, 
                    n = top_n, 
                    order_by = !!sym(rank_spec_list[[method]]$method_score), 
                    with_ties = FALSE)
        
        # order_by is assigned the criterion dictated by rank_spec_list
        
        
    }
    
    # Format the top_ranks data frame for future processing steps.
    topranks <- topranks %>% 
      unite("LR_Pair", c(ligand, receptor), remove = FALSE, sep = "_") %>%
      relocate("LR_Pair", .after = last_col())
    
    
    
    return(topranks)
    
  } #end of function
    
    
    
  }

  
  # dilute_Resource()
  {
    
    
  #' Dilutes a resource with randomly generated interactions derived from 
  #' specific genes
  #'
  #' @param resource The resource (as a tibble) which you would like to 
  #' falsify / dilute with random gene relationships.
  #' @param data_set The data set (as a Seurat object) from which you will draw 
  #' genes. Diluting the resource with genes that are actually found in the data
  #' set the resource will be used with guarantees that the random relationships
  #' will be relevant to the method. This function excludes dilution with genes 
  #' already present in the base resource and by extension genes in the top 
  #' ranking.
  #' @param top_rank_df As a tibble. Ideally the output of get_top_n_ranks. 
  #' Should contain the top ranked rows of the liana_wrapper output for the 
  #' method you are diluting for (using the undiluted resource). Since we are 
  #' comparing how resource dilution affects top ranked interactions, this list 
  #' of top rankings using the undiluted resource ensures that the top 
  #' interactions can still be caught in the same way. It's the effect of 
  #' surrounding diluted noise that we'll be picking up on.
  #' @param dilution_prop As a number between 0-1. The proportion of rows of the
  #' resource to dilute. Top ranked rows can't be diluted. If attaining the 
  #' requested dilution proportion requires overwriting top ranked 
  #' interactions, the function throws an error and returns nothing instead.
  #' @param dilution_type Choose "generic" or "variable". Generic dilutes with 
  #' generic genes from the count matrix, variable dilutes from variable features.
  #' 
  #' @return Returns a tibble that can be used as a resource for the liana 
  #' call_method functions but has a certain (marked) percentage of it replaced 
  #' with random nonsensical interactions.
  
  # Format the dilute_Resource Function to make diluted OP resources for each method
  dilute_Resource <- 
      function(resource, top_rank_df, dilution_prop, data_set, dilution_type){
        
    # Generate a list of gene names that will be relevant by getting them from 
    # the data_set. Depending on what type of dilution is requested, we pull our
    # genes from the variable or normal section of the seurat object.
        
    if (dilution_type == "generic") {
    
    gene_name_list  <- as.list(rownames(data_set@assays$RNA@data))
    
    
    } else if (dilution_type == "variable") {
      
      gene_name_list <- as.list(data_set@assays$RNA@var.features)
      
      
    } else {
      
      warning("Type of dilution was not set properly. Returning null.")
      return()
      
    }
    
    
    # Remove gene names already present in OmniPath to ensure that each diluted 
    # relationship is not a top ranked CCI and is definetly new to OmniPath
    gene_name_list <- gene_name_list %>% 
      discard(~ .x %in% resource$source_genesymbol) %>%
      discard(~ .x %in% resource$target_genesymbol)
    
    
    
    # Separate top_rank parts of resource from non top rank parts of resource, 
    # we only want to dilute the latter
    resource_top <- resource %>% 
      filter(LR_Pair %in% top_rank_df$LR_Pair)
    
    resource_bottom <- resource %>% 
      filter(!(LR_Pair %in% top_rank_df$LR_Pair))
    
    
    
    # Determine how many rows of resource_bottom to dilute so that the overall 
    # dilution_prop is met
    dilution_number <- round(nrow(resource)*dilution_prop)
    
    # Additional Warning message and break if the dilution prop can't be met
    if(dilution_number > nrow(resource_bottom)) {
      
      warning("Requested dilution proportion not attainable without overwriting 
                the given top-ranked interacions. Returning nothing instead.")
      
      return()
      
      
    }
    
    # Select the candidates for dilution from resource_bottom
    set.seed(90)
    resource_dilute <- slice_sample(resource_bottom, n = dilution_number)
    
    # Delete the dilution candidates from resource_bottom so we have an even 
    # split into dilution candidates and non-diluted candidates
    resource_bottom <- anti_join(resource_bottom, resource_dilute)
    
    
    # Sample random gene names from the gene list into the dilution candidates, 
    # creating random nonsensical relationships
    set.seed(91)
    resource_dilute$source_genesymbol <-
      as.character(sample(gene_name_list,
                          size = nrow(resource_dilute),
                          replace = TRUE))
    
    
    set.seed(92)
    resource_dilute$target_genesymbol <-
      as.character(sample(gene_name_list, 
                          size = nrow(resource_dilute), 
                          replace = TRUE))
    
    # Update the LR-Pair column with the new random "interaction partners", and 
    # mark the random interactions as such.
    resource_dilute <- resource_dilute %>%
      select(-LR_Pair) %>%
      unite("LR_Pair", c(source_genesymbol, target_genesymbol), remove = FALSE) %>%
      mutate(isRandom = TRUE) %>%
      relocate("LR_Pair", .after = last_col())
    
    
    
    # The new resource has top ranked interactions, non-top rank but still real 
    # interactions, and diluted random interactions.
    new_resource <- bind_rows(resource_top, resource_bottom, resource_dilute)

    
    #return output
    return(new_resource)
    
  } #end of function
    
    
    
  }

  
  # rank_overlap()
  {
    
    
  #' Takes get_n_top_ranks outputs that have an LR_ID and determines their overlap
  #'
  #' @param main_ranks A tibble of 
  #' @param met The method for which you would like to extract the top ranked 
  #' interactions, as a string.
  #' @param verbose Should the function describe the overlap to you or not?
  #'
  #' @return The overlap (0-1) between the two input frames in contents of the 
  #' LR_ID column, as well as an optional print statement that gives more detail.
  
  rank_overlap <- function(main_ranks, comparison_ranks, verbose = TRUE) {
    
    # calculate overlap between LR_IDs
    overlap <- sum(comparison_ranks$LR_ID %in% main_ranks$LR_ID)
    percentage_overlap <- overlap/nrow(main_ranks)
    
    
    # describe the output to the user
    if(verbose == TRUE) {
      print(str_glue("The main ranking and the comparison ranking have ",
                     as.character(overlap),
                     " LR_IDs in common. Which corrsponds to a ",
                     as.character(percentage_overlap*100),
                     "% overlap."))        
    }

    
    # put out a warning if the rankings are not of the same length, in this case 
    # the overlap percentage is only based on the main_ranks, which might catch 
    # the user off-guard
    if(nrow(main_ranks) != nrow(comparison_ranks)) {
      warning("Rankings are not of same length. 
              Percentage based on main_ranks argument.")
    }
    
    
    return(percentage_overlap)
    
  } #end of function

    
  }
  



