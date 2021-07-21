
# 1.
### Does the percentage of OmniPath genes in the results rise as OmniPath becomes more diluted?
# run after liana results after dilution are in env
gene_names_results <- unique(c(liana_results_OP$connectome$OmniPath_0$ligand, liana_results_OP$connectome$OmniPath_0$receptor))
gene_names_resource <- unique(c(resources_OP$connectome$OmniPath_0$source_genesymbol, resources_OP$connectome$OmniPath_0$target_genesymbol))
percentage_resource_in_results <- sum(gene_names_resource %in% gene_names_results)/length(gene_names_resource)
print("OmniPath_0")
print(percentage_resource_in_results)


for (i in names(dilution_props)) {
print(i)
gene_names_results <- unique(c(liana_results_OP$connectome[[i]]$ligand, liana_results_OP$connectome[[i]]$receptor))
gene_names_resource <- unique(c(resources_OP$connectome[[i]]$source_genesymbol, resources_OP$connectome[[i]]$target_genesymbol))

percentage_resource_in_results <- sum(gene_names_resource %in% gene_names_results)/length(gene_names_resource)
print(percentage_resource_in_results)

}

rm(i, gene_names_resource, gene_names_results, percentage_resource_in_results)



# 2. 
## does the number of interactions rise when lapplying over multiple OP resources regardless of dilution?
# run after dilutions are in env
test_list <- list(OmniPath_0 = resources_OP$connectome$OmniPath_0, 
                  OmniPath_0 = resources_OP$connectome$OmniPath_0, 
                  OmniPath_0 = resources_OP$connectome$OmniPath_0, 
                  OmniPath_0 = resources_OP$connectome$OmniPath_0)



test_result_list <- lapply(test_list, call_connectome, seurat_object = testdata)
lapply(test_result_list, all.equal, test_result_list$OmniPath_0)
## the results is the same all four times, and have the same number of rows


# 3.
## if we filter OP_0 to be just hits and then dilute hits with more hits the number of rows, reflecting the number of interactions, shouldn't rise
# run after dilutions are in env
gene_names <- rownames(testdata@assays$RNA@data)
OmniPath_filter_0 <- resources_OP$connectome$OmniPath_0 %>%
                        filter(source_genesymbol %in% gene_names) %>%
                        filter(target_genesymbol %in% gene_names)
                
filter_connectome_0 <- call_connectome(seurat_object = testdata, op_resource = OmniPath_filter_0)
top_ranks_filter <- get_top_n_ranks(data_set =filter_connectome_0, method = "connectome", top_n = 200)

OmniPath_filter_60 <- dilute_Resource(resource = OmniPath_filter_0, top_rank_df = top_ranks_filter, dilution_prop = 0.6, data_set = testdata)
filter_connectome_60 <- call_connectome(seurat_object = testdata, op_resource = OmniPath_filter_60)
# the number of rows of filter_connectome_0 and filter_connectome_60 are the same


# as a test, we compare our connectome at filter_0 results with the connectome results of the standard
normal_connectome <- call_connectome(seurat_object = testdata, op_resource = resources_OP$connectome$OmniPath_0) %>%
  arrange_at(vars(everything()))
filter_connectome_0 <- call_connectome(seurat_object = testdata, op_resource = OmniPath_filter_0) %>%
  arrange_at(vars(everything()))

all.equal(normal_connectome, filter_connectome_0)
# the two results are identical







#4.
## Checking if a filtered (but undiluted) OP resource produces the same as the standard non-filtered one
## filtered means it only includes genes that are also present in the data you're going to use it with
## since gene interactions for genes not in the data should have no impact, this should make no difference.
## run when dilutionsare in env

gene_names <- rownames(testdata@assays$RNA@data)
OmniPath_filter <- resources_OP$connectome$OmniPath_0 %>%
  filter(source_genesymbol %in% gene_names) %>%
  filter(target_genesymbol %in% gene_names)

filter_connectome <- call_connectome(op_resource = OmniPath_filter, seurat_object = testdata) %>%
                      arrange_at(vars(everything()))

all.equal(arrange_at(liana_results_OP$connectome$OmniPath_0, vars(everything())), filter_connectome)


filter_italk <- call_italk(op_resource = OmniPath_filter, seurat_object = testdata) %>%
  arrange_at(vars(everything()))

all.equal(arrange_at(liana_results_OP$italk$OmniPath_0, vars(everything())), filter_italk)





OmniPath_filter_OR <- resources_OP$connectome$OmniPath_0 %>%
  filter((source_genesymbol %in% gene_names) | (target_genesymbol %in% gene_names)) 


filter_cellchat <- call_cellchat(op_resource = OmniPath_filter_OR, seurat_object = testdata) %>%
  arrange_at(vars(everything()))

all.equal(arrange_at(liana_results_OP$cellchat$OmniPath_0, vars(everything())), filter_cellchat)



filter_sca <- call_sca(op_resource = OmniPath_filter_OR, seurat_object = testdata) %>%
  arrange_at(vars(everything()))

all.equal(arrange_at(liana_results_OP$sca$OmniPath_0, vars(everything())), filter_sca)


## works for italk and connectome but not the others, not even using an OR instead of AND makes a difference