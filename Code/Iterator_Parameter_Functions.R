summarise_Metadata <- function(runtime,
                               time_of_run,
                               dilution_params,
                               master_seed_list,
                               save_results = TRUE,
                               trial_run    = TRUE) {
  
  meta_params <- list("time_of_run"  = time_of_run,
                      "save_results" = save_results,
                      "trial_run"    = trial_run)
  
  
  metadata <- list(
    "runtime"         = runtime,
    "dilution_params" = append(dilution_params, 
                               master_seed_list),
    "meta_params"     = meta_params,
    "sessionInfo"     = sessionInfo()
  )
  
  
  if (save_results == TRUE) {
    # Generate the filepaths to save the data under
    metadata[["box_plot_png_name"]] <-
      auto_file_Name(
        prefix = "Boxplot_Resource_Dilution_",
        suffix = ".png",
        dilution_params = formals(wrap_resource_Robustness),
        meta_params     = formals(summarise_Metadata),
        time_of_run     = time_of_run
      )
    
    metadata[["line_plot_png_name"]] <-
      auto_file_Name(
        prefix = "Lineplot_Resource_Dilution_",
        suffix = ".png",
        dilution_params = formals(wrap_resource_Robustness),
        meta_params     = formals(summarise_Metadata),
        time_of_run     = time_of_run
      )
    
    metadata[["env_save_path"]] <-
      auto_file_Name(
        prefix = "Outputs/DilutionEnv_",
        suffix = ".RData",
        dilution_params = formals(wrap_resource_Robustness),
        meta_params     = formals(summarise_Metadata),
        time_of_run     = time_of_run
      )
  }

  
  if(dilution_params$sink_output == TRUE) {
    # The file name includes many script_params in it to be informative and
    # unique.
    metadata[["sink_logfile"]] <- auto_file_Name(prefix = "Outputs/Logs/Complete_Log_",
                                   suffix =  ".txt",
                                   dilution_params = formals(wrap_resource_Robustness),
                                   meta_params = formals(summarise_Metadata),
                                   time_of_run = time_of_run)
    
    
  }
  
  # If necessary, create a log name, store it in script params, then remove
  # the leftover clutter
  if(dilution_params$liana_warnings == "divert") {
    # The file name includes many script_params in it to be informative and
    # unique.
    metadata[["warning_logfile"]] <- auto_file_Name(prefix = "Outputs/Logs/LIANA_warnings_",
                                      suffix =  ".txt",
                                      dilution_params = formals(wrap_resource_Robustness),
                                      meta_params = formals(summarise_Metadata),
                                      time_of_run = time_of_run)
    
    
  }
  
  
  return(metadata)
  
}
