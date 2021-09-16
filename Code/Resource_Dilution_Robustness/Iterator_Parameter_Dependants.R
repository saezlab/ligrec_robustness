save_Results <- function(dilution_params,
                         meta_params, 
                         testdata_type,
                         time_of_run) {
  
  if (meta_params$trial_run == FALSE) {
    test_run_comment <- ""
    
  } else if (meta_params$trial_run == TRUE) {
    test_run_comment <- "TRIAL_RUN_"
    
  }
  
  if (dilution_params$preserve_topology == FALSE) {
    topology_comment <- "_random_topology_"
    
  } else if (dilution_params$preserve_topology == TRUE) {
    topology_comment <- "_preserved_topology_"
    
  }
  
  top_ranks_vector <- 
    unlist(as.list(dilution_params$number_ranks)[-1])
  
  # Generate the filepaths to save the data under
  box_plot_png_name <-
    auto_file_Name(
      prefix = "Boxplot_Resource_Dilution_",
      suffix = ".png",
      dilution_params = dilution_params,
      meta_params     = meta_params,
      testdata_type   = testdata_type,
      time_of_run     = time_of_run
    )
  
  line_plot_png_name <-
    auto_file_Name(
      prefix = "Lineplot_Resource_Dilution_",
      suffix = ".png",
      dilution_params = dilution_params,
      meta_params     = meta_params,
      testdata_type   = testdata_type,
      time_of_run     = time_of_run
    )
  
  
  env_save_path <- auto_file_Name(
    prefix = "Outputs/Resource_Dilution/DilutionEnv_",
    suffix = ".RData",
    dilution_params = dilution_params,
    meta_params     = meta_params,
    testdata_type   = testdata_type,
    time_of_run     = time_of_run
  )
  
  
  
  # Save both plots
  ggsave(
    plot = plot_box,
    box_plot_png_name,
    height = 7.75,
    width = 8.00,
    path = "Outputs/Resource_Dilution"
  )
  
  ggsave(
    plot = plot_line,
    line_plot_png_name,
    height = 9.00,
    width = 8.00,
    path = "Outputs/Resource_Dilution"
  )
  
  
  # Save R environment and all the results within it
  save.image(file = env_save_path)
  
  # Let the user know where everything was stored.
  print(str_glue("Box Plot saved at ~/Outputs/", 
                 box_plot_png_name, "."))
  
  print(str_glue("Line Plot saved at ~/Outputs/", 
                 line_plot_png_name, "."))
  
  print(str_glue("Environment saved at ~/", 
                 env_save_path, "."))
  
}
