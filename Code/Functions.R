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
           
           "call_connectome" = list(method_score = "weight_sc",
                               descending_order =  TRUE),
           
           "call_italk"      = list(method_score = "logfc_comb",
                               descending_order =  TRUE),
           
           "call_natmi"      = list(method_score = "edge_specificity",
                               descending_order =  TRUE),
           
           "call_sca"        = list(method_score = "LRscore",
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
  #' 
  #' @param data_set The data set (as a Seurat object) from which you will draw 
  #' genes for dilution. Diluting the resource with genes that are actually 
  #' found in the data set the resource will be used with guarantees that the 
  #' random relationships will be relevant to the method. This function excludes
  #' dilution with genes already present in the base resource and by extension 
  #' genes in the top ranking.
  #' 
  #' @param top_rank_df As a tibble. Ideally the output of get_top_n_ranks. 
  #' Should contain the top ranked rows of the liana_wrapper output for the 
  #' method you are diluting for (using the undiluted resource). Since we are 
  #' comparing how resource dilution affects top ranked interactions, this list 
  #' of top rankings using the undiluted resource ensures that the top 
  #' interactions can still be caught in the same way. It's the effect of 
  #' surrounding diluted noise that we'll be picking up on.
  #' 
  #' @param dilution_prop As a number between 0-1. The proportion of rows of the
  #' resource to dilute. Top ranked rows can't be diluted. If attaining the 
  #' requested dilution proportion requires overwriting top ranked 
  #' interactions, the function throws an error and returns nothing instead.
  #' 
  #' @param dilution_feature_type Choose "generic" or "variable". Generic dilutes with 
  #' generic genes from the count matrix, variable dilutes from variable 
  #' features.
  #' 
  #' @param preserve_topology Either TRUE or FALSE. If TRUE, a dilution method 
  #' is used where among the interactions that are to be diluted, each 
  #' instance of a given gene is replaced with a fake gene counterpart. This way
  #' the topology among those diluted interactions remains the same, while
  #' all the original gene names are replaced. This preserves the overall 
  #' network topology of the resource better, but not perfectly. The topology is
  #' not preserved perfectly because the topology of genes that appeared in both
  #' the undiluted and diluted parts is split by the renaming of the diluted 
  #' part.
  #' 
  #' If FALSE, a method of dilution is chosen that randomly assigns interaction
  #' pairs to the parts of the resource meant for dilution. This alters the
  #' network topology, but the number of edges and number of unique edges is
  #' preserved. Additionally, among the diluted interactions,  no source can be 
  #' target to itself and there is no overlap between the sources and targets 
  #' (each added gene is either a source in all of its interactions or a target,
  #' but not both). This mimics similar constraints in OmniPath, where for 
  #' example, source and target overlap is uncommon (but possible). The topology
  #' received for this method is highly dependent on the number of dilution 
  #' genes provided, but in turn it will also function with far less genes than 
  #' if the topology is to be preserved.
  #' 
  #' @param verbose Set to TRUE by default. Produces a detailed output that 
  #' summarizes important elements of the resource and how they changed through
  #' dilution. The output can be used to double check the parameters that were
  #' used and whether or not the dilution went correctly. 
  #' 
  #' Normally, the number of unique edges should not change, all edges should be 
  #' unique, the percentage of edges MARKED diluted should be equal to the 
  #' percentage that ARE diluted and should be close to the given dilution 
  #' proportion. The overlap in diluted edges should be 0 if 
  #' preserve_topology == FALSE, otherwise it should be close to the source and 
  #' target overlap before dilution. No sources should be targets to themselves,
  #' ever.
  #' 
  #' @return Returns a tibble that can be used as a custom OmniPath resource for
  #' liana_wrap and the call_method functions but has a certain (marked) 
  #' percentage of it replaced with random nonsensical interactions.
  
    
  dilute_Resource <- 
      function(resource, top_rank_df, dilution_prop, data_set, 
               dilution_feature_type, preserve_topology, verbose = TRUE) {
  
        
        
    # Separate top_rank parts of resource from non top rank parts of resource, 
    # we only want to dilute the latter
    resource_top <- resource %>% 
      filter(LR_Pair %in% top_rank_df$LR_Pair)
    
    resource_bottom <- resource %>% 
      filter(!(LR_Pair %in% top_rank_df$LR_Pair))
    
    # Determine how many rows of resource_bottom to dilute so that the overall 
    # dilution_prop is met.
    dilution_number <- round(nrow(resource)*dilution_prop)
    
    # Additional Warning message and break if the dilution prop can't be met
    # without touching top-ranked interactions which need to be protected.
    if(dilution_number > nrow(resource_bottom)) {
      
      warning(str_glue("REQUESTED DILUTION PROPORTION NOT ATTAINABLE WITHOUT ",
                       "OVERWRITING THE GIVEN TOP-RANKED INTERACIONS. ",
                       "RETURNING NOTHING INSTEAD."))
      
      return()
      
      
    }
    
    # Select the candidates for dilution from resource_bottom
    set.seed(90)
    resource_dilute <- slice_sample(resource_bottom, n = dilution_number)
    
    # Delete the dilution candidates from resource_bottom so we have an even 
    # split into dilution candidates and non-diluted candidates
    resource_bottom <- anti_join(resource_bottom, resource_dilute)
    
    
    
    
    
    # Generate a list of gene names that will be relevant by getting them from 
    # the data_set. Depending on what type of dilution is requested, we pull our
    # genes from the variable or normal section of the seurat object.
    
    if (dilution_feature_type == "generic") {
      
      gene_name_list  <- as.list(rownames(data_set@assays$RNA@data))
      
      
    } else if (dilution_feature_type == "variable") {
      
      gene_name_list <- as.list(data_set@assays$RNA@var.features)
      
      
    } else {
      
      warning("TYPE OF DILUTION WAS NOT SET PROPERLY. RETURNING NULL.")
      return()
      
    }
    
    
    # Remove gene names already present in OmniPath to ensure that each diluted 
    # relationship is not a top ranked CCI and is definetly new to OmniPath
    gene_name_list <- gene_name_list %>% 
      discard(~ .x %in% resource$source_genesymbol) %>%
      discard(~ .x %in% resource$target_genesymbol)
    
    
    # Does not attempt to fully preserve dilution candidate or resource topology
    # but also doesn't randomly add interactions either. There will be no
    # duplicates, no sources that are their own target, and no source and target
    # overlap.
    if (preserve_topology == FALSE) {
      # When generating combinations from two distinct lists, the highest number
      # of interactions are achieved when both are of as equal size as possible.
      # Here we divide the given gene_name_list into two such sublists, randomly
      # in order to avoid ordering bias in gene_name_list. We name them source 
      # and target and will keep them distinct, as in real life gene-products 
      # that are both ligands and receptors are rare. This restriction makes our
      # dilutions more true to impersonating the overall structure of omnipath, 
      # while still being nonsense.
      set.seed(91)
      source_gene_list <- sample(gene_name_list, 
                                 size = ceiling(length(gene_name_list)/2))
      
      target_gene_list <- gene_name_list %>% 
        discard(~ .x %in% source_gene_list)
      
      
      max_possible_interactions <- 
        length(source_gene_list) * length(target_gene_list)
      
      
      # If we need to create more source/target interactions than is possible in
      # the best case scenario given the gene_name_lists we have, we throw an 
      # error and don't even try to compute it.
      if(max_possible_interactions < dilution_number) {
        warning(str_glue("THE NUMBER OF GENES GIVEN CANT CREATE AS MANY ",
                         "RANDOM INTERACTIONS AS THE DILUTION PROP REQUIRES. ",
                         "RETURNING NULL."))
        return()
      }
      
      
      # ST: Stands for Source and Target
      diluted_ST_interactions <- expand.grid(source = source_gene_list, 
                                             target = target_gene_list)
      
      diluted_ST_interactions$source <- unlist(diluted_ST_interactions$source)
      diluted_ST_interactions$target <- unlist(diluted_ST_interactions$target)
      
      diluted_ST_interactions        <- tibble(diluted_ST_interactions)
      
      
      
      set.seed(92)
      resource_dilute[c("source_genesymbol", "target_genesymbol")] <-
        slice_sample(diluted_ST_interactions, n = dilution_number)
      
      
      rm(source_gene_list, 
         target_gene_list, 
         diluted_ST_interactions, 
         max_possible_interactions)
      
      
    }
    
    # Fully preserves dilution candidates topology by  switching gene names.
    # Almost preserves full resource topology, however if there are genes that
    # were amongst dilution candidates and the rest of the resource, their
    # topology will be altered.

    if (preserve_topology == TRUE)  {
      
      
      resource_genes <- unlist(unique(c(resource_dilute$source_genesymbol, 
                                        resource_dilute$target_genesymbol)))
      
      if(length(gene_name_list) < length(resource_genes)) {
        warning(str_glue("NOT ENOUGH GENES FROM DATA SET PROVIDED TO PRESERVE ",
                         "TOPOLOGY. RETURNING NULL."))
        return()
      }
      
      set.seed(93)
      dilution_genes <- unlist(sample(gene_name_list, 
                                      size = length(resource_genes)))
      
      
      dilution_mapping <- tibble(resource_genes, dilution_genes)
      
      
      
      for(i in 1:nrow(resource_dilute)) {
        
        gene_to_dilute  <- resource_dilute$source_genesymbol[i]
        
        GTD_index <- which(dilution_mapping$resource_genes == gene_to_dilute)
        
        resource_dilute$source_genesymbol[i] <-
          dilution_mapping[GTD_index, "dilution_genes"]
        
        
      }
      
      
      for(i in 1:nrow(resource_dilute)) {
        
        gene_to_dilute  <- resource_dilute$target_genesymbol[i]
        
        GTD_index <- which(dilution_mapping$resource_genes == gene_to_dilute)
        
        resource_dilute$target_genesymbol[i] <-
          dilution_mapping[GTD_index, "dilution_genes"]
        
        
      }
      
      resource_dilute <- resource_dilute %>%
        unnest(cols = c(source_genesymbol,  target_genesymbol))
      
      
      rm(i, 
         gene_to_dilute, 
         GTD_index, 
         resource_genes, 
         dilution_genes, 
         dilution_mapping)
      
    }
    
    
    # Update the LR-Pair column with the new random "interaction partners", and 
    # mark the random interactions as such.
    resource_dilute <- resource_dilute %>%
      select(-LR_Pair)                 %>%
      unite("LR_Pair", 
            c(source_genesymbol, target_genesymbol), 
            remove = FALSE)            %>%
      mutate(isRandom = TRUE)          %>%
      relocate("LR_Pair", .after = last_col())
    
    
    # The new resource has top ranked interactions, non-top rank but still real 
    # interactions, and diluted random interactions.
    new_resource <- bind_rows(resource_top, resource_bottom, resource_dilute)
    
     
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    {
    # calculate what proportion of edges is marked as diluted
    prop_marked_diluted <- sum(new_resource$isRandom)*100/nrow(new_resource)
    
    prop_diluted <- sum(!(new_resource$LR_Pair %in%
                            resource$LR_Pair))   *100   /nrow(new_resource)
    
    # calculate percentage overlap between source & target inresource
    before_overlap <- 
      sum(resource$source_genesymbol %in% 
            resource$target_genesymbol)          *100   /nrow(resource)
    
    after_overlap <- 
      sum(new_resource$source_genesymbol %in% 
            new_resource$target_genesymbol)      *100   /nrow(new_resource)
    
    diluted_overlap <- 
      sum(resource_dilute$source_genesymbol %in% 
            resource_dilute$target_genesymbol)   *100   /nrow(resource_dilute)
    
    edges_okay <- all(sapply(list(nrow(resource), 
                                  nrow(new_resource),
                                  length(unique(resource$LR_Pair)),
                                  length(unique(new_resource$LR_Pair))), 
                             function(x) x == nrow(resource)))
    
    if(edges_okay == FALSE) {
      warning(str_glue("DILUTION ERROR: THE NUMBER OF EDGES OR NUMBER OF ",
                       "UNIQUE EDGES IS NOT CONSTANT."))
    }
    
    
    
    
    
    
    dilution_marker_okay <- prop_marked_diluted == prop_diluted
    dilution_prop_okay <- abs(dilution_prop*100 - prop_diluted) > 5
    
    if(dilution_marker_okay == FALSE) {
      warning(str_glue("DILUTION ERROR: THE PROPORTION OF INTERACTIONS ",
                       "MARKED AS DILUTED (isRandom) DOES NOT CORRESPOND TO ",
                       "THE PROPORTION OF EDGES THAT ARE DILUTED. SOMETHING ",
                       "IS WRONG WITH THE isRandom MARKINGS."))
    }
    
    if(dilution_prop_okay) {
      warning(str_glue("DILUTION WARNING: THE ACTUAL DILUTION PROPORTION IS ",
                       "MORE THAN 5 PERCENTAGE POINTS DIFFERENT TO THE USER ",
                       "SPECIFIED DILUTION PROPORTION. SOMETHING IS WORNG ",
                       "WITH THE NUMBER OF DILUTED ROWS."))
    }
    
    
    
    
    
    
    if(preserve_topology == TRUE) {
      
      topology_okay <- abs(before_overlap - after_overlap) <= 10
      
      if(topology_okay == FALSE) {
        warning(str_glue("DILUTION WARNING: THE SOURCE AND TARGET  OVERLAP ",
                         "CHANGED MORE THAN 10 PERCENTAGE POINTS THROUGH ",
                         "DILUTION. THE TOPOLOGY HAS CHANGED DRASTICALLY."))
      }
      
    } else if(preserve_topology == FALSE) {
      
      topology_okay <- diluted_overlap == 0
      
      if(topology_okay == FALSE) {
        warning(str_glue("DILUTION ERROR: THE SOURCE AND TARGET  OVERLAP ",
                         "SHOULD BE 0 FOR THIS DILUTION TYPE AND ISN'T. THE ",
                         "TOPOLOGY IS NOT AS IT SHOULD BE."))
      }
    }
    
    
    
    
    
    no_self_targets <- all(sapply(list(sum(resource$source_genesymbol == 
                                             resource$target_genesymbol), 
                                       sum(new_resource$source_genesymbol == 
                                             new_resource$target_genesymbol)), 
                                  function(x) x == 0))
    
    if(no_self_targets == FALSE) {
      warning(str_glue("DILUTION ERROR: THE ORIGINAL RESOURCE AND/OR THE ",
                       "DILUTED RESOURCE CONTAIN(S) INTERACTIONS WITH ",
                       "SOURCES THAT ARE TARGETS TO THEMSELVES."))
    }
    }
      
      if (verbose == TRUE) {
      # define output divider so that 80 characters per line aren't exceeded 
      # below
      output_divider <- 
        "------------------------------------------------------------
        ------------------------------------------------------------"
      
      
      #a bunch of info on the new and old data set
      print(str_glue(""))
      print(str_glue(""))
      print(str_glue(as.character(dilution_prop*100), "% Dilution Complete"))
      print(str_glue(""))
      
      if(preserve_topology == TRUE) {
        print(str_glue("Topology semi-preserved."))
        print(str_glue(""))
      } else if (preserve_topology == FALSE) {
        print(str_glue("Topology unpreserved, but: 
                         - No duplicate interactions.
                         - No sources that are their own target.
                         - No overlap between diluted sources and targets."))
        print(str_glue(""))
      }
      
      if(dilution_feature_type == "generic") {
        print(str_glue("Diluted with a list of all genes in the input Seurat."))
      } else if (dilution_feature_type == "variable") {
        print(str_glue("Diluted with the variable features in the input ",
                       "Seurat."))
      }
      
      
      print(str_glue(""))
      print(str_glue(output_divider))
      print(str_glue("Number of edges before:                           ", 
                     nrow(resource)))
      print(str_glue("Number of edges after:                            ", 
                     nrow(new_resource)))
      print(str_glue(""))
      
      print(str_glue("Number of unique edges before:                    ", 
                     length(unique(resource$LR_Pair))))
      print(str_glue("Number of unique edges after:                     ", 
                     length(unique(new_resource$LR_Pair))))
      print(str_glue(""))
      print(str_glue(output_divider))
      
      
      
      
      print(str_glue("Number of edges marked  as diluted (after):       ", 
                     sum(new_resource$isRandom)
      ))
      print(str_glue(""))
      
      print(str_glue("Proportion of edges marked as diluted:            ", 
                     round(prop_marked_diluted, 2),
                     " %"))
      print(str_glue("Proportion of edges actually diluted:             ", 
                     round(prop_marked_diluted, 2),
                     " %"))
      print(str_glue(""))
      print(str_glue(output_divider))
      
      
      
      
      print(str_glue("Source and Target overlap before:                 ", 
                     round(before_overlap, 2),
                     " %"))
      print(str_glue("Source and Target overlap after:                  ", 
                     round(after_overlap, 2),
                     " %"))
      print(str_glue(""))
      print(str_glue("Source and target overlap in diluted edges:       ", 
                     round(diluted_overlap, 2),
                     " %"))
      print(str_glue(""))
      print(str_glue(output_divider))
      
      
      
      
      print(str_glue("Sources that are Targets to themselves, before:   ", 
                     sum(resource$source_genesymbol == 
                           resource$target_genesymbol)))
      print(str_glue("Sources that are Targets to themselves, after:    ", 
                     sum(new_resource$source_genesymbol == 
                           new_resource$target_genesymbol)))
      
      print(str_glue(""))
      print(str_glue(""))
  
    
    
    
    

    }
      
      rm(dilution_marker_okay,
         dilution_prop_okay,
         edges_okay,
         no_self_targets,
         topology_okay)
      
      
      rm(before_overlap,
         after_overlap,
         diluted_overlap,
         prop_marked_diluted,
         prop_diluted,
         output_divider)
    
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
    #return output
    return(new_resource)
  } #end of function
    
    
    
  }

  
  # rank_overlap()
  {
    
    
  #' Takes get_n_top_ranks outputs that have an LR_ID and determines their overlap
  #'
  #' @param main_ranks A tibble of of top ranked interactions
  #' @param comparison_ranks A tibble of top ranked interactions
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
  
  # prop_isRandom()
  {
    #' Takes top_rank tibbles and checks what proportion are diluted interactions
    #'
    #' @param top_rank_df A tibble of top ranked interactions. get_top_n_ranks 
    #' output. Requires isRandom column
    #'
    #' @return The proportion of interactions within the top_rank_df that are 
    #' diluted interactions.
    
    prop_isRandom <- function(top_rank_df) {
      
     FP_rate <- sum(top_rank_df$isRandom) / nrow(top_rank_df)
     
     return(FP_rate)
      
    } #end of function
    
  }
  



