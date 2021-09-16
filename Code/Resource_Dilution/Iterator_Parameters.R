#------------------------------------------------------------------------------#
# 0. Introduction and Goals ----------------------------------------------------
{
  # This is the Iterator_Params.R script.
  
  # The goal of this script is to define three central functions, the default 
  # arguments of these functions define all of the parameters needed to run the  
  # robustness iterator ("Robustness_Iterator.R"). In this manner, the iterator 
  # can always simply use formals() on the functions defined here when it needs
  # to know the iterator parameters that were specified. 
  
  # formals(fun) is a base R function that returns a list of arguments for the 
  # function you passed it. For each argument, the default value of that 
  # argument is returned. Define the iterator parameters in the defaults of the
  # functions here, and the iterator will run with them. When a function in the
  # iterator wants to know what parameters were used for the iterator, it can 
  # use formals(fun) to refer to your inputs here. Because formals returns the
  # default as "language" objects, it's important that the defaults are set in 
  # the manner specified in the parameter description. For example, individual 
  # values should be passed directly, not in a c() vector. If you are getting
  # errors in these functions or the ones in the Iterator_Paramter_Dependents.R
  # script, it could be because the defaults are not set properly.

  # One example would be the auto_file_Name() function. It needs to know the 
  # parameters used in the iterator in order to be able to create informative
  # and accurate file names. Since the iterator parameters are all defined 
  # in the default arguments of the functions below, auto_file_Name() can use 
  # formals(wrap_resource_Robustness), formals(summarise_Metadata) etc. to get 
  # all the information it needs. This approach allows the user to set all 
  # their iterator parameters in this script to their needs, while avoiding 
  # having those parameters saved in editable objects in the environment. 

  # One function below is an exception, create_Params(). It handles the script 
  # parameters that can't be saved as function defaults (namely 
  # master_seed_list and time_of_run). Unfortunately, these have to be handled
  # as editable environment objects for practicality's sake. Still, the way in
  # which these objects are created is still under the control of the function 
  # defaults of create_Params().
  
  # Finally, all the functions that need to know Iterator_Params.R to function 
  # have been grouped. They are all defined in Iterator_Parameter_Dependents.R. 
  # Every single usage of those functions is restricted to this script and the
  # Robustness_Iterator.R script.
  
  # In summary, please use the argument defaults of the functions below to
  # input the script parameters you would like to run the 
  # Robustness_Iterator.R with. The functions each explain what each argument
  # will impact, and how.

  
}



#------------------------------------------------------------------------------#
# 1. Defining Functions---------------------------------------------------------

# wrap_resource_Robustness()
{
  #' Wrapper function for resource_Robustnes()
  #' 
  #' @description This function feeds its arguments to resource_Robustness. In 
  #' addition, the dilution_proportions are converted from a simple vector to a 
  #' named list, and filepaths for logs are automatically generated if needed. 
  #' 
  #' @param master_seed  At some stages when diluting  a resource, randomness
  #' is at play. Setting a master seed will ensure that dilute_Resource(), the 
  #' function that dilutes the resources, runs the same way every time. Each 
  #' instance of randomness within the scope of dilute_Resource() will run 
  #' reproducibly.
  #' 
  #' @param testdata A preprocessed seurat object of test data that LIANA++ can 
  #' run on.
  #' 
  #' @param feature_type Should dilution occur with all the genes profiled in
  #' testdata (choose "generic") or with the most variable features (choose 
  #' "variable")? Only takes a string as default, not a vector of a string.
  #' 
  #' @param preserve_topology When diluting, two methods are implemented.
  #' random_Dilute makes only a small effort to preserve the topology of the 
  #' rows that are diluted from the resource, preserve_Dilute makes a far 
  #' greater effort. Choose TRUE for preserve_Dilute and FALSE for 
  #' random_Dilute. Only takes a boolean as default, not a vector of a boolean.
  #' 
  #' @param dilution_props A sequence of numerics (0-1) that indicate which 
  #' proportions to dilute the resource with. For example c(0.1, 0.2, 0.3) would
  #' compare the top ranked CCI's using undiluted OmniPath compared to OmniPath 
  #' with 10 % of its rows diluted, undiluted vs 20 % diluted, and undiluted vs
  #' 30 % diluted.
  #' 
  #' @param bundled_outputs Which outputs of the calculation would you like to 
  #' return? Outputs are returned bundled up in a list. By default, only the 
  #' analysis of top_ranks and runtime data is returned, but more information 
  #' can be returned. Construct an atomic vector using all or some of 
  #' "liana_results_OP", "resources_OP", "top_ranks_OP", "top_ranks_analysis",
  #' "runtime", and "testdata" to tell the script which outputs should go in the
  #' output bundle.
  #' 
  #' "liana_results_OP": All the predicted CCI's from LIANA++ for each method 
  #' and dilution stage.
  #' 
  #' "resources_OP": All the resources that were used are returned.
  #' 
  #' "top_ranks_OP": All the highest ranked CCI's from LIANA++ for each method 
  #' and dilution stage.
  #' 
  #' "top_ranks_analysis": How the diluted top predictions related to the 
  #' undiluted predictions in terms of overlap, mismatch, proportion of fake 
  #' interactions in top_ranks, and proportion of mismatch caused by fake 
  #' interactions. 
  #' 
  #' "runtime": What was the runtime of this function was like.
  #' 
  #' "testdata": What was the testdata that was used in this run. Since that's 
  #' also a user supplied Seurat, there is almost never a reason to return this.
  #' 
  #' @param methods_vector Which methods should the function run? Choose from
  #' "call_connectome", "squidpy", "call_natmi", "call_italk", "call_sca" and
  #' "cellchat". Supply the argument in the form of e.g. 
  #' c("call_conncectome, "call_italk"). Only takes a vector of chars as 
  #' default.
  #' 
  #' @param number_ranks A named list. Each item is named after a method and is 
  #' equal to the number of top interactions considered relevant for that 
  #' method. Only takes a list named after methods and containing numerics as
  #' default.
  #'  
  #' @param cellchat_nperms Cellchat is one of the slower methods, for test runs
  #' it may be useful to set this parameter to 10 to speed up the analysis.
  #'
  #' @param sink_otuput TRUE or FALSE. Should the function save a full log of 
  #' the Console Output as a log to the log folder? Warnings and messages will 
  #' not be visible in the console output if this is enabled, so unless there is
  #' a reason why such a record is necessary this option is not recommended.
  #' 
  #' @param liana_warnings Either TRUE, FALSE or "divert". Should the warnings 
  #' from liana, which are often repetitive and unhelpful, be either suppressed
  #' or alternatively diverted to the log folder? When these types of warning 
  #' are left in, they can often displace valuable warnings. Be careful with 
  #' this setting, as suppressing warnings is obviously risky.
  #' 
  #' @param time_of_run What time was this script run at? Used as an argument 
  #' for auto_file_Name() to uniquely tag logs so that thee user can associate 
  #' the logs with the correct plots and data, and so that a new log never 
  #' overwrites an old one.
  #' 
  #' @return A list that contains all the data specified in outputs. At minimum,
  #' the top_ranks_analysis list will be contained, which contains various 
  #' comparisons of the top_ranks predicted at dilution_prop = 0 and the 
  #' various higher dilution proportions.
  
  wrap_resource_Robustness <-
    function(master_seed,
             testdata,
             feature_type,
             preserve_topology, 
             dilution_props,
             number_ranks,
             methods_vector,
             
             bundled_outputs,
             cellchat_nperms,
             
             sink_output,     
             liana_warnings,
             trial_run,
             time_of_run,
             testdata_type) {
      
      
      ## Process Dilution Proportions List
      {
        # By naming each dilution proportion we can use this to label data 
        # easily. By formatting proportions as a list we can lapply over them 
        # conveniently
        
        # Convert dilution.props to a list and name it, creating a named list
        dilution_props <- as.list(dilution_props)
        
        names(dilution_props) <- map(dilution_props, function(prop) {
          
          # Name every dilution proportion
          str_glue("OmniPath_", as.character(prop * 100))
          
        }) 
      }

                                     
      ## Define File Paths for Logs
      {
        # If they are required, we auto generate log filepaths here. The 
        # filepaths have the current time imminently before the lapply in them. 
        # Each run of the iterator should be assigned to a unique log this way 
        # that has all iterations in it.
        
        # Initiate empty filepaths
        sink_logfile <- ""
        warning_logfile <- ""
        
        
        # If necessary, create a log name, store it in script params, then 
        # removethe leftover clutter
        if (sink_output == TRUE) {
        # The file name includes many script_params in it to be informative and
        # unique.
          sink_logfile <-
            auto_file_Name(
              prefix = "Outputs/Resource_Dilution/Logs/Complete_Log_",
              suffix =  ".txt",
              
              preserve_topology  = preserve_topology,
              testdata_type      = testdata_type,
              feature_type       = feature_type,
              number_ranks       = number_ranks,
              time_of_run        = time_of_run,
              trial_run          = trial_run)
          
          
        }
        
        # If necessary, create a log name, store it in script params, then 
        # remove the leftover clutter
        if (liana_warnings == "divert") {
        # The file name includes many script_params in it to be informative and
        # unique.
          warning_logfile <-
            auto_file_Name(
              prefix = "Outputs/Resource_Dilution/Logs/LIANA_warnings_",
              suffix =  ".txt",
              
              preserve_topology  = preserve_topology,
              testdata_type      = testdata_type,
              feature_type       = feature_type,
              number_ranks       = number_ranks,
              time_of_run        = time_of_run,
              trial_run          = trial_run)
          
          
        }
      }
      
      
      ## Run resource robustness with the wrapper defaults and new arguments.
      resource_Robustness(
        master_seed       = master_seed,
        testdata          = testdata ,
        
        feature_type      = feature_type,
        preserve_topology = preserve_topology,
        dilution_props    = dilution_props,
        number_ranks      = number_ranks,
        bundled_outputs   = bundled_outputs,
        
        methods_vector    = methods_vector,
        cellchat_nperms   = cellchat_nperms,
        sink_output       = sink_output,
        liana_warnings    = liana_warnings,
        
        sink_logfile      = sink_logile,
        warning_logfile   = warning_logfile
      )
      
      
      
    } # end of function
  
  
}

# extract_Testdata()
{
  #' Helper function that gets a specific seurat object from the outputs folder
  #' 
  #' @param testdata_type As a string. Which testdata should be retrieved? 
  #' Either "seurat_pbmc" or "liana_test". Seurat_pbmc is the data set used in 
  #' the seurat tutorial, while liana_test is the testdata that comes with 
  #' LIANA++, and is a small subset of seurat_pbmc. Only takes a string as 
  #' default, not a vector of a string.
  #' 
  #' @return A seurat object loaded from the outputs folder or liana package.
  
  
  extract_Testdata <- function(testdata_type) {
    
    # Get seurat or liana test data
    if (testdata_type == "seurat_pbmc") {
      
      # Read testdata from outputs
      testdata <- readRDS(file = "Data/pbmc3k_final.rds")     
      
    } else if (testdata_type == "liana_test") {
      
      # Where is the liana testdata located?
      liana_path <- system.file(package = 'liana')       
      # Read the testdata from its location.
      testdata <- 
        readRDS(file.path(liana_path, "testdata", "input", "testdata.rds"))   
      
      # removing superfluous values
      rm(liana_path)
      
      
    } else {
      
      # error if its not one of the supported data sets
      stop("Testdata name not recognized!")
      
    }
    
    
    # Return the seurat.
    return(testdata)
    
  } # end of function
}

# summarise_Metadata()
{
  #' Summarizes the metadate relevant for the Robustness Iterator
  #' 
  #' @description This function summarizes all the iterator parameters and used
  #' file names into one metadata object.
  #' 
  #' @param runtime The output of the calculate_Runtime() function, showing
  #' the time spent on various steps in resource_Robustness.
  #' 
  #' @param time_of_run What time was this script run at? Used as an argument 
  #' for auto_file_Name() to generate the same file names that were used 
  #' earlier. Also stored in the metadata object. Should always be passed the 
  #' char tag generated by create_Params().
  #' 
  #' @param dilution_params What arguments were used in dilution robustness?
  #' Should always be passed formals(wrap_resource_Robustness()).
  #' 
  #' @param testdata_type What type of testdata was used in the iterator? Should
  #' always be passed formals(extract_Testdata)$testdata_type.
  #' 
  #' @param master_seed_list What master_seed_list was the iterator run with?
  #' Should always be passed the seed list output from create_Params().
  #' 
  #' @param save_results Should the results from the iterator be saved 
  #' (including plots)? The save_Results() function will later refer to this.
  #' This value is also stored in the metadata. Only takes a boolean as a
  #' default, not a vector of a boolean.
  #' 
  #' @param trial_run Is this a trial run of the iterator or serious results?
  #' Takes a boolean. If this is a trial run, the save file names, logs and plot
  #' captions will reflect this. Only takes a boolean as a default, not a vector
  #' of a boolean.
  #' 
  #' @return Returns a list of metadata and parameters.
  
  
  
  summarise_Metadata <- function(number_seeds,
                                 master_seed_list,
                                 testdata_type,
                                 feature_type, 
                                 preserve_topology,    
                                 dilution_props,
                                 number_ranks,
                                 methods_vector,
                                
                                 sink_output,    
                                 liana_warnings,
                                
                                 cellchat_nperms,       
                                 bundled_outputs,
                                 master_outputs,
                                

                                 save_results,
                                 trial_run,
                            
                                 runtime,
                                 time_of_run) {
    
    # Summarize the metadata parameters
    meta_params <- list(
      "time_of_run"  = time_of_run,
      "save_results" = save_results,
      "trial_run"    = trial_run
    )
    
    dilution_params <- list(
      "number_seeds"      = number_seeds,
      "master_seed_list"  = master_seed_list,
      "testdata_type"     = testdata_type,
      "feature_type"      = feature_type, 
      "preserve_topology" = preserve_topology,    
      "dilution_props"    = dilution_props,
      "number_ranks"      = number_ranks ,
      "methods_vector"    = methods_vector,
      
      "sink_output"       = sink_output,    
      "liana_warnings"    = liana_warnings,
      
      "cellchat_nperms"   = cellchat_nperms,       
      "bundled_outputs"   = bundled_outputs,
      "master_outputs"    = master_outputs
    )
    
    # Put all the parameters in a list
    metadata <- list(
      "runtime"         = runtime,
      "dilution_params" = dilution_params,
      "meta_params"     = meta_params,
      "sessionInfo"     = sessionInfo()
    )
    
    
    # If the results were saved, tack the file names that the saves are under 
    # onto the end of metadata.
    if (save_results == TRUE) {
      # Generate the filepaths to save the data under
      metadata[["box_plot_png_name"]] <-
        auto_file_Name(
          prefix = "Boxplot_Resource_Dilution_",
          suffix = ".png",
          
          preserve_topology  = preserve_topology,
          testdata_type      = testdata_type,
          feature_type       = feature_type,
          number_ranks       = number_ranks,
          time_of_run        = time_of_run,
          trial_run          = trial_run)
      
      metadata[["line_plot_png_name"]] <-
        auto_file_Name(
          prefix = "Lineplot_Resource_Dilution_",
          suffix = ".png",
          
          preserve_topology  = preserve_topology,
          testdata_type      = testdata_type,
          feature_type       = feature_type,
          number_ranks       = number_ranks,
          time_of_run        = time_of_run,
          trial_run          = trial_run)
      
      metadata[["iterator_results_save_path"]] <-
        auto_file_Name(
          prefix = "Outputs/Resource_Dilution/Iterator_Results_",
          suffix = ".RData",
          
          preserve_topology  = preserve_topology,
          testdata_type      = testdata_type,
          feature_type       = feature_type,
          number_ranks       = number_ranks,
          time_of_run        = time_of_run,
          trial_run          = trial_run)
    }
    
    
    # If the resource_Robustness() output was sunk and logged, append the 
    # file name of the log to the metadata.
    if (sink_output == TRUE) {
      # The file name includes many script_params in it to be informative and
      # unique.
      metadata[["sink_logfile"]] <-
        auto_file_Name(
          prefix = "Outputs/Resource_Dilution/Logs/Complete_Log_",
          suffix = ".txt",
          
          preserve_topology  = preserve_topology,
          testdata_type      = testdata_type,
          feature_type       = feature_type,
          number_ranks       = number_ranks,
          time_of_run        = time_of_run,
          trial_run          = trial_run)
      
      
    }
    
    # If a warnings log was created for resource_Robustness(), append the file
    # name of the log to the metadata.
    if (liana_warnings == "divert") {
      # The file name includes many script_params in it to be informative and
      # unique.
      metadata[["warning_logfile"]] <-
        auto_file_Name(
          prefix = "Outputs/Resource_Dilution/Logs/LIANA_warnings_",
          suffix = ".txt",
          
          preserve_topology  = preserve_topology,
          testdata_type      = testdata_type,
          feature_type       = feature_type,
          number_ranks       = number_ranks,
          time_of_run        = time_of_run,
          trial_run          = trial_run)
      
      
    }
    
    
    # return the metadata.
    return(metadata)
    
  } # end of function
  
  
  
}
