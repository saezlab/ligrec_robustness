shuffle_Clusters <- function(master_seed,
                             mismatch_prop,
                             metadata) {
  
  set.seed(master_seed)
  
  
  metadata_old <- metadata %>%
    rownames_to_column(var = "Bar_Code") %>%
    as_tibble()
  
  metadata_old$cluster_key <- metadata_old$cluster_key %>%
    as.numeric()
  
  meta_clusters_dilute <- slice_sample(metadata_old, prop = mismatch_prop) %>%
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
  
  meta_clusters_dilute <- meta_clusters_dilute
  
  
  metadata_new <-
    left_join(metadata_old, 
              meta_clusters_dilute, 
              by = "Bar_Code",
              suffix = c("_old", ""))
  
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
    str_glue("Cluster annotations reshuffled. ", 
             round((sum(metadata_new$Mismatched) / nrow(metadata_new)) *100, 2),
             " % mismatch to the original annotation."))
  

        
  
  return(metadata_new)
}


wrap_Shuffler <- function(master_seed_list,
                          mismatch_prop,
                          metadata,
                          cluster_key) {
  
  print_Title(str_glue("Creating ",
                       mismatch_prop*100, 
                       " % mismatched cluster annotations."))
  
  reshuffled_metadatas <- lapply(master_seed_list,
                                 shuffle_Clusters,
                                 mismatch_prop = mismatch_prop,
                                 metadata      = metadata)
  
  
  return(reshuffled_metadatas)
  
}