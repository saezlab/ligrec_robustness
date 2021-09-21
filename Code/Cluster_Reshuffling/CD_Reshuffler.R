shuffle_Clusters <- function(master_seed,
                             mismatch_prop,
                             metadata,
                             cluster_col) {
  
  set.seed(master_seed)
  
  metadata_old <- metadata %>%
    rownames_to_column(var = "Bar_Code") %>%
    as_tibble()
  
  metadata_old[[cluster_col]] <- metadata_old[[cluster_col]] %>%
    as.character()
  
  meta_clusters_dilute <- slice_sample(metadata_old, prop = mismatch_prop) %>%
    select(c("Bar_Code", all_of(cluster_col)))
  
  
  
  
  meta_clusters_dilute[[cluster_col]] <- 
    map(meta_clusters_dilute[[cluster_col]], function(annotation) {
      
      replacement_candidates <-
        filter(metadata_old, 
               metadata_old[[cluster_col]] != annotation)[[cluster_col]]
      
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
  
  helper_index <- which(is.na(metadata_new[[cluster_col]]))
  
  metadata_new[[cluster_col]][helper_index] <-
    metadata_new[[str_glue(cluster_col, "_old")]][helper_index]
  
  metadata_new <- metadata_new %>%
    mutate("Mismatched" = 
             ifelse(.[[str_glue(cluster_col, "_old")]] ==
                      . [[cluster_col]], FALSE, TRUE)) %>%
    select(-str_glue(cluster_col, "_old")) %>%
    as.data.frame()
  
  rownames(metadata_new) <- metadata_new$Bar_Code
  
  metadata_new <- metadata_new %>%
    select(-"Bar_Code")
  
  metadata_new[[cluster_col]] <- metadata_new[[cluster_col]] %>%
    as.factor()
  
  
  print(str_glue("Cluster annotations reshuffled. There is ", 
                 round((sum(metadata_new$Mismatched) / nrow(metadata_new)) *100, 
                       2),
                 " % mismatch to the original annotation."))
  

        
  
  return(metadata_new)
}


wrap_Shuffler <- function(master_seed,
                          reshuffle_props,
                          metadata,
                          cluster_col) {
  
  print(str_glue(""))
  print(str_glue(""))
  print(str_glue("Iteration ", master_seed))
  print(str_glue(""))
  
  reshuffled_metadatas <- lapply(reshuffle_props,
                                 shuffle_Clusters,
                                 master_seed = master_seed,
                                 metadata    = metadata,
                                 cluster_col = cluster_col)
  
  metadata_list <- list("Reshuffle_0" = metadata) %>%
    append(reshuffled_metadatas)
  
  return(metadata_list)
  
}