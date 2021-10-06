require(liana)
require(Seurat)
require(tidyverse)
require(lubridate)


testdata <- readRDS("Data/pbmc3k_final.rds")

liana_results <- liana_wrap(testdata,
                            methods = c('call_connectome',
                                        'call_natmi',
                                        'call_italk',
                                        'call_sca',
                                        'cellchat',
                                        'squidpy'),
                            resource = c("OmniPath"))


test_plot <- ggplot(liana_results$cellchat,
                    aes(x = pval)) +
  geom_histogram(bins = 20)

ggsave(plot = test_plot, 
       "cluster_test.png", 
       height = 8.00,
       width = 8.00,
       path = "Outputs")

print("Hello world!")

warning("uh oh...")

save(liana_results, "Outputs/cluster_test.RData")