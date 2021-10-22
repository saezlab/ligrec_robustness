# 11. NATMI produces more unique LR Pairs in the top ranking than other methods ----
# Other methods produce less unique LR Pairs and repeat them among more cluster combinations
# NATMI has less repeats
{
  
  # Load required Packages
  require(tidyverse)
  require(Seurat)
  require(liana)
  require(lubridate)
  
  # source extract_Testdata function
  source("Code/Utilities/Iterator_Functions.R")
  # source clust_get_top_Ranks function
  source("Code/Cluster_Reshuffling/Iterator_Top_Ranks.R")
  
  # First we load testdata from the data folder. 
  # We also give a label (testdata_type, choose "seurat_pbmc" or "liana_test")
  testdata_type <- "seurat_pbmc"  
  testdata      <- extract_Testdata(testdata_type = testdata_type)
  
  # We format the metadata so that the cluster annotations in metadata and 
  # the seurat idents are the same. This is a prerequisite for squidpy.
  testdata@meta.data <- testdata@meta.data %>%
    mutate("cluster_key" = as.factor(as.numeric((Idents(testdata)))))
  
  Idents(testdata) <-  testdata@meta.data$cluster_key
  
  
  # NATMI results are contaminated with results from earlier runs if you
  # don't specify a special output folder for the results to go in. Here we
  # define an output folder unique to this usage of this function.
  natmi_output <-  Sys.time() %>%
    as.character()       %>%
    gsub(':', '-', .)    %>%
    gsub(' ', '_', .)    %>%
    str_glue('Test_topology_', .)
  
  
  # Use the unique tag to create unique filepaths
  natmi_params <-
    list(output_dir = natmi_output,
         expr_file  = str_glue(natmi_output, "_expr_matrix.csv"),
         meta_file  = str_glue(natmi_output, "_metadata.csv"),
         reso_name  = str_glue(natmi_output, "Resource"))
  
  
  methods_vector <- c("call_connectome",
                      "call_natmi",
                      "call_italk",
                      "call_sca",
                      "cellchat",
                      "squidpy"
                    )
  
  # Convert methods_vector to a list and name it
  methods_list        <- as.list(methods_vector)
  names(methods_list) <- methods_vector
  

  liana_results <- 
    liana_wrap(testdata, 
               methods_vector,
               resource = c("OmniPath"),
               call_natmi.params = natmi_params)
  
  
  top_ranks <-
    # We extract the top ranked interactions of all our LIANA++ results
    map(methods_list, function(method) {

      top_ranks_for_method <- liana_results[[method]] %>%
        clust_get_top_ranks(method = method,
                            top_n = 500,
                            with_ties = TRUE)
      
    })  %>%

    map_depth(., .depth = 1, format_top_ranks)
  
  
  
  top_rank_edges <- 
    tibble("methods"             = methods_vector,
           "number_unique_edges" = map(top_ranks, function(top_ranks_tib) {
             
             number_unique_edges <- top_ranks_tib$LR_Pair %>%
               unique() %>%
               length()
             
           }) %>% 
             unlist()) %>%
    arrange(desc(number_unique_edges))
           

  
  
  
  
  
  topology <- map(methods_list, function(method) {
    
    # Investigate topology with plots
    topology <- top_ranks[[method]]$LR_Pair       %>%
      table()                                     %>%
      as_tibble(.name_repair = "universal")       %>% 
      rename("LR_Pairs" = ".", "Frequency" = "n") %>%
      arrange(., desc(Frequency)) %>%
      mutate("Method" = method)
    
    
    return(topology)
    
    
  }) %>% bind_rows()
  
  
  topology_plot <- ggplot(data = topology, aes(x = Frequency,
                                               group = Frequency,
                                               fill = Method)) +
    geom_histogram(binwidth = 1) +
    
    ggtitle("Distribution of Top-Ranked Interactions by Method") +
    labs(subtitle = "Histogramm of how many LR-Pairs are repeated how many times.") +
    ylab("Number of LR-Pairs") +
    xlab("Number of Repeats") +
    
    theme(legend.position = "bottom",) +  
    
    facet_wrap(~Method, nrow = 2, ncol = 3, scales = "free")
  
  
  print(topology_plot)
  
  # Save both plots
  ggsave(
    plot = topology_plot,
    "Topology_Plot.png",
    height = 7.75,
    width = 8.00,
    path = "Outputs"
  )
  
  results <- list(top_rank_edges = top_rank_edges,
                  topology = topology,
                  topology_plot = topology_plot)
  
  save(top_rank_edges, file = "Outputs/Top_Rank_Edges_supp_topology.RData")

}