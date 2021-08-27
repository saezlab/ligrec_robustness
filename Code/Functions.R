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
    
    
  #' Dilutes a subset of a resource while analysing its output
  #' 
  #' @description This function acts as a handler for the entire process of 
  #' diluting a resource. Diluting a resource is the process of replacing 
  #' interactions of a resource with new interactions not from the resource.
  #' These new interactions are obviously not supported by the literature
  #' the resource is based on, they represent spurious "fake" interactions. 
  #' From the viewpoint of the resource, they are unverified or irrelevant to 
  #' CCI.
  #' 
  #' Adding these interactions can model the addition of lower quality noise
  #' to a resource or mimic the switching from one resource to another. A 
  #' certain core remains the same while the rest is exchanged. When we rerun
  #' liana methods with a diluted resource, the degree of similarity of its 
  #' predictions gives insight into the robustness of the method as the
  #' resource changes. In particular,gradual changes in resource can be 
  #' achieved that aren't possible when a resource is swapped wholesale.
  #' 
  #' This function first formats an input resource so that only a proportion 
  #' (dilution_prop) of it is diluted. Optionally, a list of interactions 
  #' (top_rank_list) can be supplied, when the resource is diluted these 
  #' interactions will be untouched. Commonly this is used so that top_ranked 
  #' interactions aren't altered by dilution, so that only the false positive 
  #' rate is measured. Ideally, the input resource should only include
  #' interactions that are also possible in data_set. Excess rows only
  #' increase the chances that dilution un-uniformly affects the rows that are
  #' relevant. 
  #' 
  #' Once the subset of the resource that will be diluted is known, a list of
  #' genes to dilute with is necessary. Ideally , these genes should come from
  #' the data set you will later use in LIANA with the diluted resource you
  #' generate. Since the resource only contains relevant genes, it is only
  #' accurate to replace them with genes that will also be relevant. To do this,
  #' the method runs  extract_unconflicting_Genes() on data_set. It will extract
  #' genes of the given feature_type.
  #' 
  #' With these prerequesites, the method calls either random_Dilute() or
  #' preserve_Dilut() to perform the actual dilution of the subset using the 
  #' list of genes from extract_unconflicting_Genes().
  #' 
  #' Finally, the method analyses the diluted output resource in detail and
  #' outputs warnings if anything is off about it. A detailed output can also
  #' be achieved with the verbose argument.
  #'
  #' @param resource The resource (as a tibble) which you would like to 
  #' falsify / dilute with new gene relationships. Ideally, the resource 
  #' should only include interactions that are also possible in the data set.
  #' Dilution should uniformly affect the relevant parts of a resource after
  #' all, and not the rest. 
  #' 
  #' @param top_rank_list Optional. Gene interactions that should not be 
  #' diluted, given as a list of chars or vector of chars. Every item of the 
  #' list or vector should be formatted as an LR_Pair (source + _ + target, 
  #' e.g "ITGB2_ICAM2"). Often this will be the LR_Pair column of a top_ranks 
  #' data frame. For example, if the output of get_top_n_ranks() was 
  #' top_ranks_df, top_rank_list = top_ranks_df$LR_Pair would work.
  #' 
  #' @param dilution_prop As a number between 0-1. The proportion of rows of the
  #' resource to dilute. Interactions in top_rank_list can't be diluted. If 
  #' attaining the requested dilution proportion requires overwriting top ranked 
  #' interactions, the function throws an error and returns nothing instead.
  #' 
  #' @param preserve_topology Either TRUE or FALSE. If FALSE, uses 
  #' random_Dilute() to dilute the resource. If TRUE,uses preserve_Dilute() to 
  #' dilute the resource.
  #' 
  #' @param data_set A parameter for extract_unconflicting_Genes(). 
  #' Specifically, the data_set the gene names for dilution are extracted from.
  #' 
  #' @param feature_type A parameter for extract_unconflicting_Genes(). More 
  #' specifically, the type of genes to extract from data_set. More information
  #' can be found in the extract_unconflicting_Genes() documentation.
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
  #' @param master_seed  At some stages when diluting  a resource, randomness 
  #' is at play. There are four stages where a random sample is taken in this
  #' function (and its subfunctions). Once when resource is subset to produce
  #' resource_dilute, twice in random_Dilute() and once in preserve_Dilute().
  #' Using master_seed, the function will produce a unique seed on each of
  #' these occasions, making sure the function runs the same every time. The 
  #' master seed must be an integer, and is floored into one in the function.
  #' 
  #' @return Returns a tibble that can be used as a custom OmniPath resource for
  #' liana_wrap and the call_method functions but has a certain (marked) 
  #' percentage of it diluted with  interactions.
  
    
  dilute_Resource <- function(resource, 
                              top_rank_list = list(),
                              dilution_prop, 
                              preserve_topology, 
                              data_set, 
                              feature_type,
                              verbose = TRUE, 
                              master_seed = 900)         {
    
    ## 0. Sanitize master_seed input
    master_seed <- floor(master_seed)
    
    ## 1. Splitting the Resource    
    {
    # We divide the resource into three separate subsections.
    #   1. Interactions that are in top_rank_list, i.e. ones that should not be 
    #      diluted, named resource_top
    #   2. Interactions that are free to be diluted, named resource_bottom
    #   3. The specific interactions from resource_bottom that will be diluted
    #      to most closely match dilution_prop
      
      
    # Separate top-ranked parts of resource, create resource_top 
    resource_top <- resource %>% 
      filter(LR_Pair %in% top_rank_list)
    
    # Everything that is not top_ranked is open for dilution
    resource_bottom <- resource %>% 
      filter(!(LR_Pair %in% top_rank_list))
    
    # Determine how many rows of resource_bottom to dilute so that the overall 
    # dilution_prop is met.
    dilution_number <- round(nrow(resource)*dilution_prop)
    
    # Additional Warning message and break if the dilution prop can't be met
    # without diluting top-ranked interactions (which this function doesn't do).
    if(dilution_number > nrow(resource_bottom)) {
      
      warning(str_glue("REQUESTED DILUTION PROPORTION NOT ATTAINABLE WITHOUT ",
                       "OVERWRITING THE GIVEN TOP-RANKED INTERACIONS. ",
                       "RETURNING NOTHING INSTEAD."))
      
      return()
      
      
    }
    
    
    # If we dilute dilution_number of rows dilution_prop will be met
    # Copy dilution_number rows at random from resource_bottom...
    set.seed(seed = master_seed + 0.1)
    resource_dilute <- slice_sample(resource_bottom, n = dilution_number)
    
    # ... and delete them from resource bottom. 
    resource_bottom <- anti_join(resource_bottom, resource_dilute)
    
    # Currently, resource_top + resource_bottom + resource_dilute == resource
    # Once resource_dilute has been properly diluted with the methods below
    # resource_top + resource_bottom + resource_dilute* == new_resource
    # The number of rows of new_resource that are diluted are equal to
    # dilution_prop.
    
    } # end of subpoint
    
    
    ## 2. Preparing a Gene Name List
    {
    # As described above, we need a gene_name_list to create diluted 
    # interactions with. The extract_unconflicting_Genes() function extracts 
    # genes from the data_set and removes any overlap to existing genes in the
    # resource.
    
    gene_name_list <- 
        extract_unconflicting_Genes(data_set       = data_set, 
                                    feature_type   = feature_type,
                                    conflict_genes = 
                                      c(resource$source_genesymbol,
                                        resource$target_genesymbol))
    
    # With the dilution candidates in resource_dilute and dilution genes ready, 
    # we move on to the actual dilution methods.
      
    } # end of subpoint
    
        
    ## 3. Implementing Dilution Method
    {
    # There are two approaches to dilution. One is random and that follows 
    # certain topological rules but makes no further effort to match the 
    # topology of OmniPath (yes, these rules are specific to OmniPath and not 
    # generalization to other resources). The other preserves the topology of 
    # resource_dilute exactly, and thus semi-preserves the topology of resource
    # in new_resource. More details in the documentation.
      
    # random_Dilute(). There will be no duplicates, no sources that are their 
    # own target, and no source and target overlap.
    if (preserve_topology == FALSE) {
      
      resource_dilute <- random_Dilute(resource_dilute = resource_dilute,
                                       gene_name_list  = gene_name_list,
                                       master_seed     = master_seed)
      
    }
    
    # preserve_Dilute(). The topology of resource_dilute will be unaltered. 
    if (preserve_topology == TRUE)  {
      
      resource_dilute <- preserve_Dilute(resource_dilute = resource_dilute, 
                                         gene_name_list  = gene_name_list,
                                         master_seed     = master_seed)
      
    }
    
    
    # Update the LR-Pair column with the new random "interaction partners".
    resource_dilute <- resource_dilute %>%
      select(-LR_Pair)                 %>%
      unite("LR_Pair", 
            c(source_genesymbol, target_genesymbol), 
            remove = FALSE)            %>%
      relocate("LR_Pair", .after = last_col())
    
    
    # The new resource has top ranked interactions, non-top rank but still real 
    # interactions, and diluted random interactions.
    new_resource <- bind_rows(resource_top, resource_bottom, resource_dilute)
    
    } # end of subpoint
    
    
    ## 4. Output Analysis and Warnings Assessment
    {
    # At this point, the output is ready, but to make sure nothing has gone
    # wrong in the dilution process, we calculate multiple metrics to capture
    # the topology of the result.
      
    # We start by calculating some more complicated expressions to save 
    # typespace elsewhere.
    
    # calculate what proportion of edges is marked as diluted
    prop_marked_diluted <- sum(new_resource$isRandom)*100/nrow(new_resource)
    
    # calculate what proportion of edges have an LR_Pair unique to the new
    # resource
    prop_diluted <- sum(!(new_resource$LR_Pair %in%
                            resource$LR_Pair))   *100   /nrow(new_resource)
    
    # calculate percentage of source & target overlap in the old resource, 
    # new resource, and in the specific rows that were diluted.
    before_overlap <- 
      sum(resource$source_genesymbol %in% 
            resource$target_genesymbol)          *100   /nrow(resource)
    
    after_overlap <- 
      sum(new_resource$source_genesymbol %in% 
            new_resource$target_genesymbol)      *100   /nrow(new_resource)
    
    diluted_overlap <- 
      sum(resource_dilute$source_genesymbol %in% 
            resource_dilute$target_genesymbol)   *100   /nrow(resource_dilute)
    
    
    # Having our complex parameters sorted, we now start a series of evaluations
    # as to whether or not the topology looks healthy. Guidelines for this
    # are also in the documentation of the verbose argument above.
    
    # Store our evlauation results in a compact format
    warning_logic <- list()

    # The number of edges should be constant, and every edge should be unique
    warning_logic["edges_okay"]       <- 
      all(sapply(list(nrow(resource), 
                      nrow(new_resource),
                      length(unique(resource$LR_Pair)),
                      length(unique(new_resource$LR_Pair))), 
                 function(x) x == nrow(resource)))
    
    # The number of rows marked with isRandom == TRUE should be equal to the
    # number of rows with LR_Pairs foreign to the original resource
    warning_logic["dil_marker_okay"]  <- 
      prop_marked_diluted == prop_diluted
    
    # The proportion of rows diluted shouldn't be too different from the user
    # defined proportion
    warning_logic["dil_prop_okay"]    <- 
      abs(dilution_prop*100 - prop_diluted) < 5
    
    # If random_Dilute() was used, the source and target overlap should be 0
    warning_logic["r_topology_okay"]  <- 
      diluted_overlap == 0
    
    # If preserve_Dilute() was used, the source and target overlap should be
    # similar to the original input
    warning_logic["p__topology_okay"] <- 
      abs(before_overlap - after_overlap) <= 10
    
    # There should never be sources that are targets to themselves in an 
    # interaction
    warning_logic["no_self_targets"]  <- 
      all(sapply(list(sum(resource$source_genesymbol == 
                            resource$target_genesymbol),  
                      sum(new_resource$source_genesymbol ==  
                            new_resource$target_genesymbol)),  
                 function(x) x == 0))
    
    
    
    
    # In this next if statement we integrate the above knowledge, is there an 
    # issue in the data? If so, is it mild (a warning) or severe (an error).
    
    # We only want to consider r_topology_okay if random_Dilute() was used
    if (preserve_topology == FALSE) {
      
      # The dilution proportion being off is concerning but not necessarily
      # indicative that the dilution went wrong. We want to distinguish this
      # issue as a warning, not an error.
      warning_logic["warning_thrown"]   <- !(all(warning_logic$dil_prop_okay))
      
      # Any if these issues are serious and indicate the dilution failed.
      warning_logic["error_thrown"]     <- !(all(warning_logic$edges_okay,  
                                                 warning_logic$dil_marker_okay, 
                                                 warning_logic$r_topology_okay,  
                                                 warning_logic$no_self_targets))
      
    # We only want to consider p_topology_okay if preserve_Dilute() was used  
    } else if (preserve_topology == TRUE) {
      
      # As above, these two issues are concerning when they occur, but only
      # at the level of a warning. The dilution may still be fine.
      warning_logic["warning_thrown"]   <- !(all(warning_logic$dil_prop_okay,  
                                                 warning_logic$p_topology_okay))
      
      # As above, these issues are serious and indicate the dilution failed.
      warning_logic["error_thrown"]     <- !(all(warning_logic$edges_okay,  
                                                 warning_logic$dil_marker_okay, 
                                                 warning_logic$no_self_targets))
      
    }
    
    
    
    } # end of subpoint
    
      
    ## 5. Optional Verbose Output
    {
    # In this segment we print a bunch of information on key parameters of the
    # the inputs and outputs. 
      
    if (verbose == TRUE) {
    # we briefly define this divider here so the code below is more readable
    output_divider <- 
      "------------------------------------------------------------
      ------------------------------------------------------------"
    
    
    # spacing outputs helps with overview.
    print(str_glue(""))
    print(str_glue(""))
    
    # If an error occurred, this is the highest priority information. Thus this
    # is what is brought to the users attention first. In this way, the output 
    # immediately warns the user of any problems and lets the user know which 
    # output is associated with which warnings. 
    if (warning_logic$error_thrown == TRUE) {
      
      print(str_glue(as.character(dilution_prop*100), 
                     "% Dilution --- ERROR OCCURRED"))
      print(str_glue("Please check ERROR below."))
      
    
    # A warning is second priority. If there was an error that is more relevant
    # than a warning, but if there was no error, a warning has the next priority
    } else if (warning_logic$warning_thrown == TRUE) {
      
      print(str_glue(as.character(dilution_prop*100), 
                     "% Dilution --- WARNING OCCURRED"))
      print(str_glue("Please check WARNING below."))
    
      
    # Only if there were no errors or warnings do we label a dilution successful
    } else {
      
      print(str_glue(as.character(dilution_prop*100), 
                     "% Dilution Successful"))
      
    }
    
    print(str_glue(""))
    
    # The output reminds the user what the input parameters for dilution type
    # were, and labels this output. When multiple outputs are produced by 
    # iterating this function multiple times and the user is in doubt, they can
    # hopefully piece together which dilution this is the output for.
    if(preserve_topology == TRUE) {
      
      # remind the user what "semi-preserved" means.
      print(str_glue("Topology semi-preserved:
                       - Diluted interactions have identical topology.
                       - Genes NOT in dilution have identical topology.
                       - Genes that existed in both have split topology."))
      
    } else if (preserve_topology == FALSE) {
      
      # remind the user of the topological constraints in random_Dilute()
      print(str_glue("Topology unpreserved, but: 
                       - No duplicate interactions.
                       - No sources that are their own target.
                       - No overlap between diluted sources and targets."))
      
    }
    
    # space
    print(str_glue(""))
    
    # What did the feature type mean again?
    if(feature_type == "generic") {
      
      print(str_glue("Diluted with a list of all genes in the input Seurat."))
      
    } else if (feature_type == "variable") {
      
      print(str_glue("Diluted with the variable features in the input ",
                     "Seurat."))
    
    }
    
    # Results for the analysis of the edge number and uniqueness
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
    
    
    
    # Analysis of number of edges with LR_Pairs foreign to original resource
    # as well as the ones marked diluted. Also gives the achieves dilution
    # proportion compared the user defined one.
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
    
    
    
    # Limited analysis of topology in source + target overlap
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
    
    
    
    # There should never be sources that are targets to themselves within a 
    # single interaction.
    print(str_glue("Sources that are Targets to themselves, before:   ", 
                   sum(resource$source_genesymbol == 
                         resource$target_genesymbol)))
    print(str_glue("Sources that are Targets to themselves, after:    ", 
                   sum(new_resource$source_genesymbol == 
                         new_resource$target_genesymbol)))
    
    print(str_glue(""))
    print(str_glue(""))
    
    }
    
    # remove the cluttered parameters for the outputs and warnings
    rm(before_overlap,
       after_overlap,
       diluted_overlap,
       prop_marked_diluted,
       prop_diluted,
       output_divider)
 
    } # end of subpoint
    
    
    ## 6. Print Warnings
    {
    # We print the warnings here even though we've had all the knowledge on 
    # which warnings are needed for a while. We print them here because it is
    # more intuitive for the warnings to be right under the verbose output.
      
    if(warning_logic$edges_okay == FALSE) {
      warning(str_glue("DILUTION ERROR: THE INPPUT AND OUTPUT SHOULD ONLY ",
                       "INCLUDE UNIQUE INTERACTIONS AND SHOULD HAVE THE ", 
                       "SAME NUMBER OF INTERACTIONS. THIS IS NOT THE CASE."))
      
    }
    
    
    if(warning_logic$dil_marker_okay == FALSE) {
      warning(str_glue("DILUTION ERROR: THE PROPORTION OF INTERACTIONS ",
                       "MARKED AS DILUTED (isRandom) DOES NOT CORRESPOND TO ",
                       "THE PROPORTION OF INTERACTIONS THAT ACTUALLY ARE DILUTED. ",
                       "SOMETHING IS WRONG WITH THE isRandom MARKINGS."))
    }
    
    
    if(warning_logic$dil_prop_okay == FALSE) {
      warning(str_glue("DILUTION WARNING: THE ACTUAL DILUTION PROPORTION IS ",
                       "MORE THAN 5 PERCENTAGE POINTS DIFFERENT TO THE USER ",
                       "SPECIFIED DILUTION PROPORTION. SOMETHING IS WORNG ",
                       "WITH THE NUMBER OF DILUTED ROWS."))
    }
    
    
    # Only consider r_topology_okay if random_Dilute() was used
    if(preserve_topology == FALSE) {
      
      if(warning_logic$r_topology_okay == FALSE) {
        warning(str_glue("DILUTION ERROR: THE SOURCE AND TARGET  OVERLAP ",
                         "SHOULD BE 0 FOR THIS DILUTION TYPE AND ISN'T. THE ",
                         "TOPOLOGY IS NOT AS IT SHOULD BE."))
      }
      
      
    # Only consider p_topology_okay if preserve_Dilute() was used  
    } else if(preserve_topology == TRUE) {
      
      if(warning_logic$p__topology_okay == FALSE) {
        warning(str_glue("DILUTION WARNING: THE SOURCE AND TARGET  OVERLAP ",
                         "CHANGED MORE THAN 10 PERCENTAGE POINTS THROUGH ",
                         "DILUTION. THE TOPOLOGY HAS CHANGED MORE ",
                         "DRASTICALLY THAN EXPECTED FOR THIS METHOD."))
      }
      
      
    }
    
    
    if(warning_logic$no_self_targets == FALSE) {
      warning(str_glue("DILUTION ERROR: THE ORIGINAL RESOURCE AND/OR THE ",
                       "DILUTED RESOURCE CONTAINS INTERACTIONS WITH ",
                       "SOURCES THAT ARE TARGETS TO THEMSELVES."))
    }
    
    rm(warning_logic)
      
    
    } # end of subpoint
    
    
    ## 7. Return Output
    return(new_resource)
    
    
  } #end of function
    
    
    
  }

  # extract_unconflicting_Genes()
  {
    
  #' Extract and filter genes
  #' 
  #' @description Extracts a subset of genes from an input Seurat, and filters 
  #' them to not include a given list of gene names.
  #' 
  #' @param data_set The data set (as a Seurat object) from which you will draw 
  #' genes from. Must be one using RNA assay data.
  #' 
  #' @param feature_type Choose "generic" or "variable" (as a char). 
  #' 
  #' If generic, genes will be extracted from the row names of the RNA count 
  #' matrix of the input, i.e. all the genes profiled in the data_set will be 
  #' extracted. 
  #' 
  #' If variable extraction will be from the variable features of the data_set, 
  #' which have to be set up beforehand.
  #' 
  #' @param conflict_genes A list of genes that should not be in the output, as
  #' a list.
  #' 
  #' @return Returns a list of genes from the data_set that does not include any
  #' members of the conflict_genes input.
  
  
  extract_unconflicting_Genes <- function (data_set, 
                                           feature_type, 
                                           conflict_genes) {
    
    # Depending on what type of dilution is requested, we pull our genes from 
    # the variable or normal section of the seurat object.
    
    if (feature_type == "generic") {
      
      gene_name_list  <- as.list(rownames(data_set@assays$RNA@data))
      
      
    } else if (feature_type == "variable") {
      
      gene_name_list <- as.list(data_set@assays$RNA@var.features)
      
      
    } else {
      
      # Throw an error if feature_type is not "generic" or "variable".
      warning("FEATURE TYPE FOR DILUTION WAS NOT SET PROPERLY. RETURNING NULL.")
      return()
      
    }
    
    # Remove any blacklisted genes from gene_name_list
    gene_name_list <- gene_name_list %>% 
      discard(~ .x %in% conflict_genes) 
    
    
    return(gene_name_list)
    
  } # end of function
    
    
  }

  # random_Dilute()
  {
    
  #' Dilutes entire resources with a random method.
  #' 
  #' @description Given a resource replaces all source_genesymbols and 
  #' target_genesymbols with random relationships made from pairwise 
  #' combinations of gene_name_list. In the process topological rules that mimic 
  #' OmniPath are followed, namely:
  #' 
  #' 1. No duplicate interactions created. All new interactions are unique.
  #' 2. In OmniPath, genes that are sources in some relationships and targets in 
  #' others is rare. This is termed source and target overlap. All diluted 
  #' interactions will have 0 % source and target overlap. A given gene will
  #' either always be a source, or always be a target.
  #' 3. No sources that are targets to themselves.
  #' 
  #' Beyond these three rules interaction topology is random. This method 
  #' requires less genes to work with than other methods, and if a solution in 
  #' the above constraints is possible, this method will always find it. 
  #' 
  #' @param resource_dilute The resource (as a tibble) which you would like to 
  #' falsify / dilute with random gene relationships. Must have a 
  #' source_genesymbol and a target_genesymbol column. The method will introduce
  #' a column that marks all interactions in the diluted output as random. No 
  #' other columns are modified in any way (LR_Pairs not updated for example).
  #' 
  #' @param gene_name_list A list of gene names from which to generate random 
  #' pairwise relationships.
  #' 
  #' @param master_seed This method has two instances of randomness. By 
  #' supplying a master_seed the function will always run the same way. The same
  #' master seed can be used here as is used elsewhere; every time randomness
  #' is employed the master_seed is used to calculate a globally unique seed for
  #' that instance of usage.  The master seed must be an integer, and is floored
  #' into one in the function.
  #' 
  #' @return Returns a tibble that can be used as a custom OmniPath resource for
  #' liana_wrap and the call_method functions but has is entirely made up of 
  #' random interactions, which are all marked.
  
  
  random_Dilute <- function (resource_dilute, 
                             gene_name_list, 
                             master_seed = 900) {
    
    # Sanitize master_seed input
    master_seed <- floor(master_seed)
    
    # This method works by generating random interactions from gene_name_list.
    # Because OmniPath has a low overlap of sources and targets, i.e. few 
    # gene products that are both ligands and receptor, we split 
    # gene_name_list into sources and targets first. Since the highest number
    # of possible interactions in this scenario would come from equally sized
    # source and target lists, we split gene_name_list evenly.
    
    # Random sample for source_gene_list to avoid any ordering bias of 
    # gene_name_list
    set.seed(seed = master_seed + 0.2)
    source_gene_list <- sample(gene_name_list, 
                               size = ceiling(length(gene_name_list)/2))
    
    # target_gene_list is what remains of gene_name_list when source_gene_list
    # is removed
    target_gene_list <- gene_name_list %>% 
      discard(~ .x %in% source_gene_list)
    
    
    
    
    # If we need to create more source/target interactions than is possible 
    # given the gene_name_lists we have, we throw an error and don't even try 
    # to compute it.
    max_possible_interactions <- 
      length(source_gene_list) * length(target_gene_list)
    
    
    if(max_possible_interactions < nrow(resource_dilute)) {
      warning(str_glue("THE NUMBER OF GENES GIVEN CANT CREATE AS MANY ",
                       "RANDOM INTERACTIONS AS THE DILUTION PROP REQUIRES. ",
                       "RETURNING NULL."))
      return()
    }
    
    
    
    # We compute every possible unique interaction given our source and target
    # list and format the results. ST: Stands for Source and Target
    diluted_ST_interactions        <- expand.grid(source = source_gene_list, 
                                                  target = target_gene_list)
    
    diluted_ST_interactions$source <- unlist(diluted_ST_interactions$source)
    diluted_ST_interactions$target <- unlist(diluted_ST_interactions$target)
    
    diluted_ST_interactions        <- tibble(diluted_ST_interactions)
    
    
    
    # We impute interactions into resource_dilute . 
    
    # Because we are sampling from every possible unique interaction using the
    # setup where the maximum possible number of interactions is created, we 
    # always find a solution if a solution with our criteria is possible. 
    
    # We sample without replacing because we don't want duplicates.
    set.seed(seed = master_seed + 0.3)
    resource_dilute[c("source_genesymbol", "target_genesymbol")] <-
      slice_sample(diluted_ST_interactions, n = nrow(resource_dilute))
    
    
    # Mark every interaction as random, because every interaction is random.
    resource_dilute <- resource_dilute %>%
      mutate(isRandom = TRUE)          
    
    # Return Output
    return(resource_dilute)
    
  } # end of function
    
  
  }
  
  # preserve_Dilute()
  {
    
  #' Dilutes entire resources while semi-preserving topology
  #' 
  #' @description Determines every unique gene in the given resource, then 
  #' remaps every one of them on a 1:1 basis with genes provided as an argument.
  #' In the process, existing verified interactions from the resource are wholly
  #' replaced ("diluted") with false interactions. The resource is diluted. At 
  #' the same time, since the process simply renames the genes, the resource
  #' topology is 100 % preserved. However, if the output was just a subset of
  #' a resource as in dilute_Resource(), the overall topology is only 
  #' semi-preserved.
  #' 
  #' @param resource_dilute The resource (as a tibble) which you would like to 
  #' falsify / dilute with random gene relationships. Must have a 
  #' source_genesymbol and a target_genesymbol column. The method will introduce
  #' a column that marks all interactions in the diluted output as random. No 
  #' other columns are modified in any way (LR_Pairs not updated for example).
  #' 
  #' @param gene_name_list A list of gene names to be used for dilution. The 
  #' original gene names from resource_dilute will be wholly replaced with the
  #' gene names in this list.
  #' 
  #' @param master_seed This method has one instance of randomness. By supplying 
  #' master_seed the function will always run the same way. The same master seed
  #' can be used here as is used elsewhere; every time randomness is employed 
  #' the master_seed is used to calculate a globally unique seed for that 
  #' instance of usage. The master seed must be an integer, and is floored into
  #' one in the function.
  #'
  #' @return Returns a tibble that can be used as a custom OmniPath resource for
  #' liana_wrap and the call_method functions but has is entirely made up of 
  #' diluted interactions, which are all marked.
  
  preserve_Dilute <- function (resource_dilute, 
                               gene_name_list,
                               master_seed) {
    
    # Sanitize master_seed input
    master_seed <- floor(master_seed)
    
    # This method creates a dictionary, where every real gene in
    # resource_dilute gets a unique fake counterpart from gene_name_list.
    # Using the dictionary, we go through resource_dilute and swap each real 
    # gene for its counterpart. As such, the resource_dilute topology is 
    # preserved.
    
    # In order to create the dictionary, we extract all the unique gene names
    # in resource dilute
    resource_genes <- unlist(unique(c(resource_dilute$source_genesymbol, 
                                      resource_dilute$target_genesymbol)))
    
    # If gene_name_list isn't long enough to provide a 1:1 dictionary for
    # resource_genes we abort the process.
    if(length(gene_name_list) < length(resource_genes)) {
      warning(str_glue("NOT ENOUGH GENES FROM DATA SET PROVIDED TO PRESERVE ",
                       "TOPOLOGY. RETURNING NULL."))
      return()
    }
    
    # We sample to avoid ordering bias in gene_name_list. We sample as many
    # genes as are in resource_genes to create our 1:1 dictionary.
    set.seed(seed = master_seed + 0.4)
    dilution_genes <- unlist(sample(gene_name_list, 
                                    size = length(resource_genes)))
    
    # Define dictionary with resource_genes and dilution_genes
    dilution_dict <- tibble(resource_genes, dilution_genes)
    
    
    
    
    # Now that we have our dictionary, we iterate through every 
    # source_genesymbol in resource_dilute, look up its counterpart in the 
    # dictionary and replace source_genesymbol with the counterpart
    for(i in 1:nrow(resource_dilute)) {
      
      # Whats the name of the ource_genesymbol?
      gene_to_dilute  <- resource_dilute$source_genesymbol[i]
      
      # Where do I find it in the dictionary?
      GTD_index <- which(dilution_dict$resource_genes == gene_to_dilute)
      
      # Replace source_genesymbol with the counterpart at the row index
      # the source_genesymbol is at in the dictionary
      resource_dilute$source_genesymbol[i] <-
        dilution_dict[GTD_index, "dilution_genes"]
      
      
    }
    
    # Repeat the aboove process for target_genesymbol
    for(i in 1:nrow(resource_dilute)) {
      
      gene_to_dilute  <- resource_dilute$target_genesymbol[i]
      
      GTD_index <- which(dilution_dict$resource_genes == gene_to_dilute)
      
      resource_dilute$target_genesymbol[i] <-
        dilution_dict[GTD_index, "dilution_genes"]
      
      
    }
    
    # The above for loops nest resource_dilute, which we undo here. We also mark
    # all the interactions as random, because they are.
    resource_dilute <- resource_dilute %>%
      unnest(cols = c(source_genesymbol,  target_genesymbol)) %>%
      mutate(isRandom = TRUE)    
    
    # Return Output
    return(resource_dilute)
    
  } # end of function
    
    
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
  

  # print_Title()
  {
    #' Prints a nice title to the console
    #'
    #' @param title A char of what the title should read
    #' @param super TRUE or FALSE. Should this title be extra large?
    #'
    #' @return Prints a title into the console.
  
    
    print_Title <- function(title, super = FALSE) {
      
      divider <- str_glue("|=======================================",
                          "=======================================|")
      
      if(super == FALSE) {
        
        print(str_glue(divider))
        print(str_glue("   ",title))
        print(str_glue(divider))
        print(str_glue(""))
        
      } else if (super == TRUE) {
        
        print(str_glue(divider))
        print(str_glue(divider))
        print(str_glue(""))
        print(str_glue("   ",title))
        print(str_glue(""))
        print(str_glue(divider))
        print(str_glue(divider))
        print(str_glue(""))
        print(str_glue(""))
        
      }

      
      
    } #end of function
    
  }
  

