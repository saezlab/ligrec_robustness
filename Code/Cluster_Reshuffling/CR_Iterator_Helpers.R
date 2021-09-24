#------------------------------------------------------------------------------#
# A. For Resource Robustness

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
  #' @param with_ties TRUE or FALSE. An argument passed to slice_min and 
  #' slice_max respectively. When slicing the top_n, sometimes a cutoff occurs
  #' between equally ranked items. Should the items also be included, expanding 
  #' the selection (TRUE), or should preserving top_n be prioritized (FALSE). 
  #'
  #' @return Returns the tibble input cut down to the top n highest ranked
  #' interactions.
  
  get_top_ranks_clust <- function(data_set, method, top_n, with_ties = FALSE) {
    # generate a list that describes how to rank in each method to get what the
    # method considers best
    rank_spec_list <-
      list(
        "cellchat"        = list(method_score = "pval",
                                 descending_order =  FALSE),
        
        "call_connectome" = list(method_score = "weight_sc",
                                 
                                 descending_order =  TRUE),
        
        "call_italk"      = list(method_score = "logfc_comb",
                                 descending_order =  TRUE),
        
        "call_natmi"      = list(method_score = "edge_specificity",
                                 descending_order =  TRUE),
        
        "call_sca"        = list(method_score = "LRscore",
                                 descending_order =  TRUE),
        
        "squidpy"         = list(method_score = "pvalue",
                                 descending_order =  FALSE)
      )
    
    
    
    # If its a p-value method, the code will go into this if statement
    if (rank_spec_list[[method]]$descending_order == FALSE) {
      topranks <-
        slice_min(
          data_set,
          n = top_n,
          order_by = !!sym(rank_spec_list[[method]]$method_score),
          with_ties = with_ties
        )
      
      # order_by is assigned the criterion dictated by rank_spec_list
      
      
    } else {
      # if it's not one of the p-value methods
      
      topranks <-
        slice_max(
          data_set,
          n = top_n,
          order_by = !!sym(rank_spec_list[[method]]$method_score),
          with_ties = with_ties
        )
      
      # order_by is assigned the criterion dictated by rank_spec_list
      
      
    }
    
    return(topranks)
    
  } #end of function
# rank_overlap()
{
  #' Takes get_n_top_ranks outputs that have an LR_ID and determines their 
  #' overlap
  #'
  #' @param main_ranks A tibble of of top ranked interactions
  #' @param comparison_ranks A tibble of top ranked interactions
  #' @param verbose Should the function describe the overlap to you or not?
  #'
  #' @return The overlap (0-1) between the two input frames in contents of the
  #' LR_ID column, as well as an optional print statement that gives more 
  #' detail.
  
  rank_overlap <-
    function(main_ranks, comparison_ranks, verbose = TRUE) {
      # calculate overlap between LR_IDs
      overlap <- sum(comparison_ranks$LR_ID %in% main_ranks$LR_ID)
      percentage_overlap <- overlap / nrow(main_ranks)
      
      
      # describe the output to the user
      if (verbose == TRUE) {
        print(str_glue(""))
        cat(
          str_wrap(
            str_glue(
              "The main ranking and the comparison ranking have ",
              as.character(overlap),
              " LR_IDs in common. Which corresponds to a ",
              as.character(round(percentage_overlap * 100, 2)),
              "% overlap."), 
            width = 60), "\n")
      }
      
      
      # put out a warning if the rankings are not of the same length, in this 
      # case the overlap percentage is only based on the main_ranks, which might
      # catch the user off-guard
      if (nrow(main_ranks) != nrow(comparison_ranks)) {
        warning("Rankings are not of same length.
              Percentage based on main_ranks argument.")
      }
      
      
      return(percentage_overlap)
      
    } #end of function
  
  
}


# cellchat_rank_overlap()
{
  #' Takes get_n_top_ranks outputs that have an LR_ID and determines their 
  #' overlap, but handles the exceptions cellchat has
  #' 
  #' @description Cellchat top_ranks are often overflowing with interactions 
  #' that have p-values of zero. If you cut without ties in this scenario you
  #' ignore a great deal of interactions that have the same ranking as your
  #' candidates. There is also an ordering bias in the data, such that you
  #' always receive the same 500 interactions. So you get incomplete and biased
  #' top_ranks. 
  #' 
  #' To combat this we cut with ties, but that makes our overlap calculations 
  #' difficult. Percentage-overlap is an intuitive metric but if the two items
  #' being compared are not the same size the comparison is directional. A x B 
  #' is not equal to B x A. 
  #' 
  #' In order to preserve the intuitive metric but create a bidirectional 
  #' comparison, we estimate the overlap one would get if one randomly 
  #' subsampled B to be the same size as A, and then calculated the overlap.
  #' 
  #' If A is entirely in B, the average overlap between A and b is P(overlap) *
  #' length(A) / length(A). P(overlap) is (A in B) / length(B), since we sample
  #' values from B uniformly. Thus, the average overlap of A and b is 
  #' (A in B) / length(B). 
  #'
  #' @param main_ranks A tibble of of top ranked interactions
  #' @param comparison_ranks A tibble of top ranked interactions
  #' @param verbose Should the function describe the overlap to you or not?
  #'
  #' @return The overlap (0-1) between the two input frames in contents of the
  #' LR_ID column, as well as an optional print statement that gives more 
  #' detail.
  
  cellchat_rank_overlap <-function(main_ranks, 
                                   comparison_ranks, 
                                   verbose = TRUE) {
    
    # calculate overlap between LR_IDs
    overlap <- sum(comparison_ranks$LR_ID %in% main_ranks$LR_ID)
    
    # If the size of both rankings is the same, calculating overlap is simple
    if (nrow(main_ranks) == nrow(comparison_ranks)) {
      
      overlap <- sum(comparison_ranks$LR_ID %in% main_ranks$LR_ID)
      percentage_overlap <- overlap / nrow(main_ranks)
      
    }
    
    
    # Otherwise, we determine the average overlap one would achieve when 
    # downsampling B to the size of A, given that A in entirely in B.
    if(overlap == nrow(main_ranks)) {
      
      percentage_overlap <- overlap / nrow(comparison_ranks)
      
    } else {
      
      # This case scenario shouldn't come up if top_ranks were preserved.
      stop("The original cellchat top_ranks could not be found in the diluted
           cellchat top_ranks. This means that cellchat_rank_overlap() can't 
           calculate what you want it to.")
      
    }
    
    
    # describe the output to the user
    if (verbose == TRUE) {
      print(str_glue(""))
      cat(
        str_wrap(
          str_glue(
            "The main ranking and the comparison ranking have ",
            as.character(overlap),
            " LR_IDs in common. Which corresponds to a ",
            as.character(round(percentage_overlap * 100, 2)),
            "% overlap."), 
          width = 60), "\n")
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
    divider <- str_glue(
      "|=======================================",
      "=======================================|"
    )
    
    if (super == FALSE) {
      print(str_glue(""))
      print(str_glue(""))
      print(str_glue(divider))
      print(str_glue(""))
      print(str_glue("   ", title))
      print(str_glue(""))
      print(str_glue(divider))
      print(str_glue(""))
      print(str_glue(""))
      print(str_glue(""))
      
    } else if (super == TRUE) {
      print(str_glue(""))
      print(str_glue(""))
      print(str_glue(""))
      print(str_glue(divider))
      print(str_glue(divider))
      print(str_glue(""))
      print(str_glue(""))
      print(str_glue("   ", title))
      print(str_glue(""))
      print(str_glue(""))
      print(str_glue(divider))
      print(str_glue(divider))
      print(str_glue(""))
      print(str_glue(""))
      print(str_glue(""))
      print(str_glue(""))
      
    }
    
    
    
  } #end of function
  
}


# divert_Warnings()
{
  #' Replaces warnings of code with a logfile
  #' 
  #' @description A function within which you can execute code. Warnings
  #' generated by the code are suppressed from the console and instead logged
  #' elsewhere.
  #'
  #' @param code Some expression to execute that produces warnings
  #' @param logFile The location you would like to log the warning at.
  #'
  #' @return Prints a title into the console.
  
  #Use a calling Handler to supress warning
  
  
  divert_Warnings <- function(code, logFile) {
    
    withCallingHandlers(code, 
                        warning = function(w) 
                        {
                          # cat warnings to a specified logFile
                          cat(as.character(Sys.time()), # time of warning
                              "\n",
                              as.character(w), # the warning as it would print
                              "\n",
                              file = logFile, 
                              append = TRUE)
                          
                          # don't print warning in console
                          invokeRestart("muffleWarning")
                        }
    )
    
  } # end of function
  
}

