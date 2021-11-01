# Load packages
library(tidyverse)

# Load data and name it
{
}


# Create a function for plotting and saving plots
final_plot_fun <- function(plot_data,
                           title,
                           plot_caption,
                           x_variable,
                           x_label,
                           png_name) {
  
  
  final_plot <- 
    ggplot(data = plot_data, aes(x = .data[[x_variable]],
                                 y = Overlap,
                                 group = .data[[x_variable]],
                                 color = Method)) +
    
    # Add boxplots, no outliers because we add points in a second
    geom_boxplot(outlier.shape = NA) + 
    geom_point(alpha = 0.4) +  # plot semi transparent plots over them so
    # you can see where the distirbution comes from
    
    # add text
    ggtitle(title) +
    ylab("Overlap of Top Ranks [%]") +
    xlab(x_label) +
    labs(caption = plot_caption,
         color = "Method") +
    
    # remove gray background so the transparent points are more visible
    theme_bw() + 
    
    # make sure the scale is the same on both axes, since both are in percent
    scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, 100)) +
    scale_x_continuous(breaks = seq(0, 60, 10), limits = c(-2, 62)) +
    
    # adjust legend and caption position
    theme(plot.caption = element_text(hjust = 0),
          legend.position = "bottom",
          plot.margin = unit(c(0.5,1,0.5,0),"cm")) +    
    
    # facet wrap for each method, so they don't overlap and are easily 
    # comparable between subplots
    facet_wrap(~Method, nrow = 1, ncol = 6)
  
  print(final_plot)
  
  
  # plot
  ggsave(
    plot = final_plot,
    png_name,
    height = 14,
    width  = 20,
    path   = "Outputs",
    units  = "cm"
  )
  
  
  return(final_plot)
  
}

# Plot relevant data
final_plot_modified_resource <- 
  final_plot_fun(plot_data = modified_resource,
                 title        = "Indiscriminate Resource Dilution Robustness",
                 plot_caption = "This plot shows the robustness results for default resource dilution and a modified baseline.",
                 x_variable   = "dilution_prop",
                 x_label      = "Dilution of Resource [%]",
                 png_name     = "Final_Default_Dilution_Modify.png")

final_plot_conserved_resource <- 
  final_plot_fun(plot_data = conserved_resource,
                 title        = "Discriminant Resource Dilution Robustness",
                 plot_caption = "This plot shows the robustness results for default resource dilution and a conserved baseline.",
                 x_variable   = "dilution_prop",
                 x_label      = "Dilution of Resource [%]",
                 png_name     = "Final_Default_Dilution_Conserve.png")

final_plot_subset_clusters <- 
  final_plot_fun(plot_data = subset_clusters,
                 title        = "Cluster Subsetting Robustness",
                 plot_caption = "This plot shows the robustness results for default cluster subsetting.",
                 x_variable   = "Mismatch",
                 x_label      = "Proportion of Cells Removed [%]",
                 png_name     = "Final_Default_Cluster_Subsetting.png")

final_plot_reshuffled_clusters <- 
  final_plot_fun(plot_data = reshuffled_clusters,
                 title        = "Cluster Reshuffling Robustness",
                 plot_caption = "This plot shows the robustness results for default cluster reshuffling.",
                 x_variable   = "Mismatch",
                 x_label      = "Proportion of Cluster Annotation Mismatch [%]",
                 png_name     = "Final_Default_Cluster_Reshuffling.png")



## top ranks

# subset
final_plot_subset_1000 <- 
  final_plot_fun(plot_data = subset_1000,
                 title        = "Cluster Subsetting with 1000 Significant Interactions",
                 plot_caption = "This plot shows the robustness results for 1000 cluster subsetting.",
                 x_variable   = "Mismatch",
                 x_label      = "Proportion of Cells Removed [%]",
                 png_name     = "Final_1000_Cluster_Subsetting.png")

final_plot_subset_1500 <- 
  final_plot_fun(plot_data = subset_1500,
                 title        = "Cluster Subsetting with 1500 Significant Interactions",
                 plot_caption = "This plot shows the robustness results for 1500 cluster subsetting.",
                 x_variable   = "Mismatch",
                 x_label      = "Proportion of Cells Removed [%]",
                 png_name     = "Final_1500_Cluster_Subsetting.png")

# reshuffle
final_plot_reshuffle_1000 <- 
  final_plot_fun(plot_data = reshuffle_1000,
                 title        = "Cluster Reshuffling with 1000 Significant Interactions",
                 plot_caption = "This plot shows the robustness results for 1000 cluster reshuffling",
                 x_variable   = "Mismatch",
                 x_label      = "Proportion of Cluster Annotation Mismatch [%]",
                 png_name     = "Final_1000_Cluster_Reshuffling.png")

final_plot_reshuffle_1500 <- 
  final_plot_fun(plot_data = reshuffle_1500,
                 title        = "Cluster Reshuffling with 1500 Significant Interactions",
                 plot_caption = "This plot shows the robustness results for 1500 cluster reshuffling",
                 x_variable   = "Mismatch",
                 x_label      = "Proportion of Cluster Annotation Mismatch [%]",
                 png_name     = "Final_1500_Cluster_Reshuffling.png")

# dilution
final_plot_dilution_1000 <- 
  final_plot_fun(plot_data = dilute_1000,
                 title        = "Indiscriminate Resource Dilution with 1000 Significant Interactions",
                 plot_caption = "This plot shows the robustness results for 1000 resource dilution and a modified baseline.",
                 x_variable   = "dilution_prop",
                 x_label      = "Dilution of Resource [%]",
                 png_name     = "Final_1000_Dilution_Modify.png")

final_plot_dilution_1500 <- 
  final_plot_fun(plot_data = dilute_1500,
                 title        = "Indiscriminate Resource Dilution with 1500 Significant Interactions",
                 plot_caption = "This plot shows the robustness results for 1500 resource dilution and a modified baseline.",
                 x_variable   = "dilution_prop",
                 x_label      = "Dilution of Resource [%]",
                 png_name     = "Final_1500_Dilution_Modify.png")


## variant analysis
final_plot_dilution_preserve <- 
  final_plot_fun(plot_data = preserve_dilution,
                 title        = "Discriminant Dilution with Preserved Resource Topology",
                 plot_caption = "This plot shows the robustness results for preserve resource dilution and a modified baseline.",
                 x_variable   = "dilution_prop",
                 x_label      = "Dilution of Resource [%]",
                 png_name     = "Final_preserve_Dilution.png")

final_plot_dilution_generic <- 
  final_plot_fun(plot_data = generic_dilution,
                 title        = "Discriminant Resource Dilution with Generic Genes",
                 plot_caption = "This plot shows the robustness results for generic resource dilution and a modified baseline.",
                 x_variable   = "dilution_prop",
                 x_label      = "Dilution of Resource [%]",
                 png_name     = "Final_generic_Dilution.png")



## calculations and plotting for topology analysis

topology <- results$top_rank_proportions

# recode proper method names
topology <- topology %>%
  mutate("methods" = recode(methods,
                           "call_connectome" = "Connectome",
                           "squidpy"         = "CellPhoneDB",
                           "call_natmi"      = "NATMI",
                           "call_italk"      = "LogFC Product",
                           "call_sca"        = "SingleCellSignalR",
                           "cellchat"        = "CellChat")) %>%
  rename("Method" = methods)
  

# consistent color scheme with other plots despite the reordering
library(RColorBrewer)
myColors <- brewer.pal(5,"Set1")
names(myColors) <- levels(topology$Method)
colScale <- scale_colour_manual(name = "Method",values = myColors)



# create the bar plot
topology_plot <- 
  ggplot(data = topology, aes(x = reorder(Method, number_unique_edges),
                               y = number_unique_edges,
                              fill = Method)) +
  
  # Add boxplots, no outliers because we add points in a second
  geom_bar(stat = "identity") +  
  
  # add text
  ggtitle("Degreeness of Method Predictions") +
  xlab("") +
  ylab("Unique Interactions per Significant Interaction") +

  
  # adjust legend and caption position
  theme(plot.caption = element_text(hjust = 0),
        legend.position = "bottom",
        plot.margin = unit(c(0.5,1,0.5,0),"cm")) +  
  coord_flip() +
  colScale



# save the bar plot
ggsave(
  plot = topology_plot,
  "Final_Topology_Plot.png",
  height = 20,
  width  = 20,
  path   = "Outputs",
  units  = "cm"
)
