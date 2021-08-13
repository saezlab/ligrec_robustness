
# 1.
### Does the percentage of OmniPath genes in the results rise as OmniPath becomes more diluted? ----
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
## does the number of interactions rise when lapplying over multiple OP resources regardless of dilution? ----
# run after dilutions are in env
test_list <- list(OmniPath_0 = resources_OP$connectome$OmniPath_0, 
                  OmniPath_0 = resources_OP$connectome$OmniPath_0, 
                  OmniPath_0 = resources_OP$connectome$OmniPath_0, 
                  OmniPath_0 = resources_OP$connectome$OmniPath_0)



test_result_list <- lapply(test_list, call_connectome, seurat_object = testdata)
lapply(test_result_list, all.equal, test_result_list$OmniPath_0)
## the results is the same all four times, and have the same number of rows


# 3. does the number of interactions in connectome rise when diluting only hits? ----
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
## Checking if a filtered (but undiluted) OP resource produces the same as the standard non-filtered one ----
## filtered means it only includes genes that are also present in the data you're going to use it with
## since gene interactions for genes not in the data should have no impact, this should make no difference.
## run when dilutions are in env

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

## works for italk and connectome when comparing to liana_wrap but not the others, but this is because liana_wrap and call_x don't work the same somehow
#### cellchat and sca are a bit different. I have to use the call function here but then the result is the same. See more in section below
filter_cellchat <- call_cellchat(op_resource = OmniPath_filter, seurat_object = testdata) %>%
  arrange_at(vars(everything()))
filter_cellchat_call <- call_cellchat(op_resource = select_resource(c("OmniPath"))[[1]], seurat_object = testdata) %>%
  arrange_at(vars(everything()))

all.equal(filter_cellchat_call, filter_cellchat)


filter_sca <- call_sca(op_resource = OmniPath_filter, seurat_object = testdata) %>%
  arrange_at(vars(everything()))
filter_sca_call <- call_sca(op_resource = select_resource(c("OmniPath"))[[1]], seurat_object = testdata) %>%
  arrange_at(vars(everything()))

all.equal(filter_sca_call, filter_sca)




#5.
#### Call_X and liana wrap don't produce the same results ----
#### even when same resource, same method, same data
## not the same output somehow
filter_sca_call <- call_sca(op_resource = select_resource(c("OmniPath"))[[1]], seurat_object = testdata)
filter_sca_wrap <- liana_wrap(seurat_object = testdata, method = c("sca"), resource = c("OmniPath"))[[1]]
## not the same output somehow
filter_cellchat_call <- call_cellchat(op_resource = select_resource(c("OmniPath"))[[1]], seurat_object = testdata)
filter_cellchat_wrap <- liana_wrap(seurat_object = testdata, method = c("cellchat"), resource = c("OmniPath"))[[1]]




# 6. Plot and save results for dilution, top 250 and top 1000 ranks -------------------------

load("C:/Users/plabu/OneDrive/Documents/GitHub/ligrec_robustness/Data/env_after_1000.RData")

top_rank_overlap_1000 <- top_rank_overlap_1000 %>%
  unnest(dilution_prop) %>%
  add_row(dilution_prop = 0, 
          connectome = 1, 
          cellchat = 1, 
          italk = 1, 
          sca = 1, 
          .before = 1) %>%
  "*"(100)



ggplot(data = top_rank_overlap_1000) + 
  geom_line(mapping = aes(dilution_prop, connectome, color =  "Connectome")) +
  geom_line(mapping = aes(dilution_prop, cellchat, color = "CellChat")) +
  geom_line(mapping = aes(dilution_prop, italk, color = "iTALK")) +
  geom_line(mapping = aes(dilution_prop, sca, color = "SCA")) +
  
  geom_point(mapping = aes(dilution_prop, connectome, color =  "Connectome")) +
  geom_point(mapping = aes(dilution_prop, cellchat, color = "CellChat")) +
  geom_point(mapping = aes(dilution_prop, italk, color = "iTALK")) +
  geom_point(mapping = aes(dilution_prop, sca, color = "SCA")) +
  
  ylim(0, 100) +
  
  ggtitle("Robustness of Method Predictions") +
  ylab("Overlap of Top Ranks [%]") +
  xlab("Dilution of Resource [%]") +
  labs(subtitle = "Generic dilution, top 1000 ranks",
       color = "Method")

ggsave("top_rank_overlap_1000_plot.png", 
       height = 5, width = 8, 
       path = "Outputs")

load("C:/Users/plabu/OneDrive/Documents/GitHub/ligrec_robustness/Data/env_after_250_high_res.RData")

top_rank_overlap_250 <- top_rank_overlap_250 %>%
  unnest(dilution_prop) %>%
  add_row(dilution_prop = 0, 
          connectome = 1, 
          cellchat = 1, 
          italk = 1, 
          sca = 1, 
          .before = 1) %>%
  "*"(100)



ggplot(data = top_rank_overlap_250) + 
  geom_line(mapping = aes(dilution_prop, connectome, color =  "Connectome")) +
  geom_line(mapping = aes(dilution_prop, cellchat, color = "CellChat")) +
  geom_line(mapping = aes(dilution_prop, italk, color = "iTALK")) +
  geom_line(mapping = aes(dilution_prop, sca, color = "SCA")) +
  
  geom_point(mapping = aes(dilution_prop, connectome, color =  "Connectome")) +
  geom_point(mapping = aes(dilution_prop, cellchat, color = "CellChat")) +
  geom_point(mapping = aes(dilution_prop, italk, color = "iTALK")) +
  geom_point(mapping = aes(dilution_prop, sca, color = "SCA")) +
  
  ylim(0, 100) +
  
  ggtitle("Robustness of Method Predictions") +
  ylab("Overlap of Top Ranks [%]") +
  xlab("Dilution of Resource [%]") +
  labs(subtitle = "Generic dilution, top 250 ranks",
       color = "Method")

ggsave("top_rank_overlap_250_plot.png", 
       height = 5, width = 8, 
       path = "Outputs")



















## what is the unification of top ranks ----


unification_top_ranks <- top_ranks_OP_0$call_connectome[1:4] %>%
  add_row(top_ranks_OP_0$call_natmi[1:4]) %>%
  add_row(top_ranks_OP_0$call_italk[1:4]) %>%
  add_row(top_ranks_OP_0$call_sca[1:4]) %>%
  add_row(top_ranks_OP_0$cellchat[1:4])

















## are tow rsult outputs idetical when ordered? ----
all.equal(arrange_at(t, vars(everything())), 
          arrange_at(liana_results_OP$call_connectome$OmniPath_0, vars(everything()))
          )








## Investigating the relationship of OP source and target genesymbols ----
#get OP_resource and construct LR Pairs
op <- select_resource(c('OmniPath'))[["OmniPath"]] %>%
  select(source_genesymbol,
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
         entity_type_intercell_target) %>%
  mutate(isRandom = FALSE) %>%
  unite("LR_Pair", 
        c(source_genesymbol, target_genesymbol), 
        remove = FALSE, 
        sep = "_") %>%
  relocate("LR_Pair", .after = last_col())



# Are there duplicate Ligand-receptor pairs?
LR_Pair_list <- op$LR_Pair
length(unique(LR_Pair_list)) - length(LR_Pair_list) # No, every LR pair is unique

# Are there genes that are both source and target?
source_list <- op$source_genesymbol
target_list <- op$target_genesymbol

length(intersect(source_list, target_list)) # yes there are this many genes that occur as both source and target
length(unique(c(source_list, target_list))) # there are this many unique genesymbols in OP
length(intersect(source_list, target_list)) / length(unique(c(source_list, target_list))) # this proportion of genes in OP occur as both sources and targets


# Are any of these gene symbols source and target to themselves?
sum(source_list == target_list) # no, even though some genes appear in both lists they never appear twice in the same row, i.e. in relationship to themselves




# is thres = 1 actually a standard parameter of liana_wrap? ----
test_thing_no_thresh <- 
  liana_wrap(testdata, 
             method = c('cellchat'), 
             resource = c('OmniPath'), 
             expr_prop = 0,
             cellchat.params = list(nboot = cellchat_nperms, 
                                    expr_prop = 0),
             call_natmi.params = list(output_dir = natmi_output))

test_thing <- 
  liana_wrap(testdata, 
             method = c('cellchat'), 
             resource = c('OmniPath'), 
             expr_prop = 0,
             cellchat.params = list(nboot = cellchat_nperms, 
                                    expr_prop = 0,
                                    thresh = 1),
             call_natmi.params = list(output_dir = natmi_output))




# 11. NATMI produces more unique LR Pairs in the top ranking than other methods ----
# Other methods produce less unique LR Pairs and repeat them among more cluster combinations
# NATMI has less repeats

length(unique(top_ranks_OP$call_natmi$OmniPath_0$LR_Pair))

# Compare to the others
length(unique(top_ranks_OP$call_connectome$OmniPath_0$LR_Pair))
length(unique(top_ranks_OP$call_italk$OmniPath_0$LR_Pair))
length(unique(top_ranks_OP$call_sca$OmniPath_0$LR_Pair))
length(unique(top_ranks_OP$cellchat$OmniPath_0$LR_Pair))