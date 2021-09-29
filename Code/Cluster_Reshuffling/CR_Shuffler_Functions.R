#------------------------------------------------------------------------------#
# 0. Introduction and Goals ----------------------------------------------------
{
}



#------------------------------------------------------------------------------#
# 1. Define Functions ----------------------------------------------------------

# shuffle_Clusters()
{
  #' Reshuffle Cluster annotations in a Seurat
  #' 
  #' @description This function takes a metadata table from a Seurat Object and
  #' partially reshuffles the cluster annotations of the cells in the table.
  #'  
  #' N.B. that whenever an annotation is replaced, the replacement is drawn from
  #' the distribution of annotations already in the table that are not equal to 
  #' the original annotation. This means that if n % of annotations are 
  #' reshuffled, there will also be exactly n % mismatch between the old and new 
  #' cluster annotations.
  #' 
  #' 
  #' @param master_seed As an integer. Cluster reshuffling has inherent 
  #' randomness to it, by supplying a seed you ensure that this function runs 
  #' reproducibly.
  #' 
  #' @param mismatch_prop As a fraction. What proportion of the cluster 
  #' annotations should be reshuffled? Also, what proportion of cluster 
  #' annotations in the new metadata should differ from the old? 
  #' (Both of these proportions are the same.)
  #' 
  #' @param metadata The metadata df of a Seurat who's cell cluster annotations
  #' you want to reshuffle. Note that the metadata should have numeric cluster
  #' labels stored as a factor in a column called "cluster_key" of the metadata.
  #' The input metadata should also have the cell bar codes as its rownames.
  #' 
  #' 
  #' @return A new metadata data frame that has partially reshuffled/mismatched
  #' cluster annotations in the "cluster_key" column.
  
  
  shuffle_Clusters <- function(master_seed,
                               mismatch_prop,
                               metadata) {
    
    # There is randomness in reshuffling, so set the seed for reproducibility
    set.seed(master_seed)
    
    # We format the input df as a tibble and store that as metadata_old
    metadata_old <- metadata %>%
      rownames_to_column(var = "Bar_Code") %>%
      as_tibble()
    
    # Convert the cluster_key column froma  factor to an integer
    metadata_old$cluster_key <- metadata_old$cluster_key %>%
      as.numeric()
    
    
    
    # We select the rows we will reshuffle from the metadata. We keep the bar 
    # code so we can remerge the old and reshuffled metadata later.
    meta_clusters_dilute <-
      slice_sample(metadata_old, prop = mismatch_prop) %>%
      select(c("Bar_Code", all_of("cluster_key")))
    
    
    
    
    meta_clusters_dilute$cluster_key <-
      map(meta_clusters_dilute$cluster_key, function(annotation) {
        replacement_candidates <-
          filter(metadata_old,
                 metadata_old$cluster_key != annotation)$cluster_key
        
        new_annotation <- sample(replacement_candidates, 1)
        
        return(new_annotation)
        
      }) %>%
      unlist()
    
    
    
    metadata_new <-
      left_join(
        metadata_old,
        meta_clusters_dilute,
        by = "Bar_Code",
        suffix = c("_old", "")
      )
    
    helper_index <- which(is.na(metadata_new$cluster_key))
    
    metadata_new$cluster_key[helper_index] <-
      metadata_new[["cluster_key_old"]][helper_index]
    
    metadata_new <- metadata_new %>%
      mutate("Mismatched" =
               ifelse(.$cluster_key_old == .$cluster_key, FALSE, TRUE)) %>%
      select(-"cluster_key_old") %>%
      as.data.frame()
    
    rownames(metadata_new) <- metadata_new$Bar_Code
    
    metadata_new <- metadata_new %>%
      select(-"Bar_Code")
    
    
    
    metadata_new$cluster_key <- metadata_new$cluster_key %>%
      as.factor()
    
    print(
      str_glue(
        "Cluster annotations reshuffled. ",
        round((
          sum(metadata_new$Mismatched) / nrow(metadata_new)
        ) * 100, 2),
        " % mismatch to the original annotation."
      )
    )
    
    
    
    
    return(metadata_new)
  } # end of function
}


# wrap_Shuffler()
{
  #' shuffle_Clusters Wrapper
  wrap_Shuffler <- function(master_seed_list,
                            mismatch_prop,
                            metadata,
                            cluster_key) {
    print_Title(str_glue(
      "Creating ",
      mismatch_prop * 100,
      " % mismatched cluster annotations."
    ))
    
    reshuffled_metadatas <- lapply(
      master_seed_list,
      shuffle_Clusters,
      mismatch_prop = mismatch_prop,
      metadata      = metadata
    )
    
    
    return(reshuffled_metadatas)
    
  } # end of function
}
