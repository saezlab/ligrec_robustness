# 0.1 Overview: 
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

# 1.1 Running LIANA wrapper
# We use five methods with OmniPath on the liana test data (squidpy is not used as it's currently not functioning for me)

  # 1.1.1 Get test data
  {
      liana_path <- system.file(package = 'liana')                                    # get liana package filepath
      testdata <- 
        readRDS(file.path(liana_path, "testdata", "input", "testdata.rds"))           # retrieve testdata from filepath
  }
  
  # 1.1.2 Run wrapper on testdata for omnipath x (cellchat, connectome, italk, natmi, sca), (starting out with just connectome for now)
  {  
  liana_results_OP_0 <- liana_wrap(testdata,
                             method = c('connectome'),
                             resource = c('OmniPath'))
  }

# 1.2 Extract top x ranked interactions for a given method method  
{
  #' Get the top n ranked items of a method from the liana wrapper results
  #'
  #' @param dat The list of tibbles output by the liana wrapper function.
  #' @param met The method for which you would like to extract the top ranked interactions, as a string.
  #' @param top_n The number of items to return, returns items ranked #1 to #n.
  #'
  #' @return Returns a tibble with all the columns for this method in the liana wrapper output but only including the rows of the top n interactions (tied values at the boundary line are cut off, no ties).

  get_top_n_ranks <- function(dat, met, top_n) {
  
    rank_spec_list <- list("cellchat" = list(met_score = "pval",descending_order =  FALSE),
                           "connectome" = list(met_score = "weight_sc", descending_order =  TRUE),
                           "italk" = list(met_score = "weight_comb", descending_order =  TRUE),
                           "natmi" = list(met_score = "edge_specificity", descending_order =  TRUE),
                           "sca" = list(met_score = "LRscore", descending_order =  TRUE),
                           "squidpy" = list(met_score = "pvalue", descending_order =  FALSE))
  
    if(rank_spec_list[[met]]$descending_order == FALSE) {
       topranks <- slice_min(dat[[met]], 
                             n = top_n, 
                             order_by = !!sym(rank_spec_list[[met]]$met_score), 
                             with_ties = FALSE)
                                     
     } else {
       topranks <- slice_max(dat[[met]], 
                             n = top_n, 
                             order_by = !!sym(rank_spec_list[[met]]$met_score), 
                             with_ties = FALSE)
    }
  
   return(topranks)
  }
  
  
  top_ranks_OP_0 <- list("connectome" = get_top_n_ranks(dat = liana_results_OP_0, top_n = 200, met = "connectome")
                       #  ,
                       #  "cellchat" = get_top_n_ranks(dat = liana_results_OP_0, top_n = 30, met = "cellchat"),
                       #  "italk" = get_top_n_ranks(dat = liana_results_OP_0, top_n = 200, met = "italk"),
                       #  "natmi" = get_top_n_ranks(dat = liana_results_OP_0, top_n = 200, met = "natmi"),
                       #  "sca" = get_top_n_ranks(dat = liana_results_OP_0, top_n = 200, met = "sca"))
  )
}

# 2.1 Modifying OmniPath with random genes

top_ranks_OP_0$connectome <- unite(top_ranks_OP_0$connectome, "LR_Pair", 3:4, remove = FALSE, sep = "_")

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

OmniPath_0 <- unite(OmniPath_0, "LR_Pair", 1:2, remove = FALSE, sep = "_")



dilute_Resource <- function(resource, top_rank_list, permutation_prop, data_set, used_method){

    gene_name_list  <- as.list(testdata@assays$RNA@var.features)
  gene_name_list <- gene_name_list[!(gene_name_list %in% resource$source_genesymbol)]
  gene_name_list <- gene_name_list[!(gene_name_list %in% resource$target_genesymbol)]
  
  resource_top <- resource %>% 
    filter(LR_Pair %in% top_rank_list[used_method]$LR_Pair)
  
 resource_bottom <- resource %>% 
    filter(!(LR_Pair %in% top_rank_list[used_method]$LR_Pair))
  
  set.seed(123)
 resource_dilute <- slice_sample(resource_bottom, prop = permutation_prop)
  resource_bottom <- anti_join(resource_bottom, resource_dilute)
  
  set.seed(123)
  resource_dilute$source_genesymbol <- as.character(sample(gene_name_list, size = nrow(resource_dilute)))
  set.seed(123)
  resource_dilute$target_genesymbol <- as.character(sample(gene_name_list, size = nrow(resource_dilute)))
  
  resource_dilute <- select(resource_dilute, -LR_Pair)
  resource_dilute <- unite(resource_dilute, "LR_Pair", 1:2, remove = FALSE, sep = "_")
  
  
  new_resource <- bind_rows(resource_top, resource_bottom, resource_dilute)
  
  return(new_resource)
}


OmniPath_10 <- dilute_Resource(resource = OmniPath_0, 
                top_rank_list = liana_results_OP_0, 
                permutation_prop = 0.1, 
                data_set = testdata, 
                used_method = "connectome")

































# OmniPath_0_top <- OmniPath_0 %>% 
#   filter(LR_Pair %in% top_ranks_OP_0$connectome$LR_Pair)
# 
# OmniPath_0_bottom <- OmniPath_0 %>% 
#   filter(!(LR_Pair %in% top_ranks_OP_0$connectome$LR_Pair))
# 
# 
# gene_name_list  <- as.list(testdata@assays$RNA@var.features)
# gene_name_list <- gene_name_list[!(gene_name_list %in% OmniPath_0$source_genesymbol)]
# gene_name_list <- gene_name_list[!(gene_name_list %in% OmniPath_0$target_genesymbol)]
# 
# OmniPath_0_dilute <- slice_sample(OmniPath_0_bottom, prop = permutation_prop)
# OmniPath_0_bottom <- anti_join(OmniPath_0_bottom, OmniPath_0_dilute)
# 
# OmniPath_0_dilute$source_genesymbol <- sample(gene_name_list, size = nrow(OmniPath_0_dilute))
# OmniPath_0_dilute$target_genesymbol <- sample(gene_name_list, size = nrow(OmniPath_0_dilute))
# 
# OmniPath_0_dilute <- select(OmniPath_0_dilute, -LR_Pair)
# OmniPath_0_dilute <- unite(OmniPath_0_dilute, "LR_Pair", 1:2, remove = FALSE, sep = "_")
# 
# 
# OmniPath_10 <- bind_rows(OmniPath_0_top, OmniPath_0_bottom, OmniPath_0_dilute)






