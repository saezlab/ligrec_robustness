#------------------------------------------------------------------------------#
# 0. Introduction and Goals ----------------------------------------------------
{
  # The goal of this script is to define three central functions, the default 
  # arguments of which define all of the parameters needed to run the the 
  # robustness iterator ("Robustness_Iterator.R"). In this manner, the iterator 
  # can always simply use formals() on the functions defined here when it needs
  # to know the iterator parameters that were specified. 

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
  # as editable environment objects for practicalities sake. Still, the way in
  # which these objects are created is still under the control of the function 
  # defaults of create_Params().
  
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
  #' @description This function feeds its defaults to resource_Robustness, 
  #' allowing other functions to determine under which exact parameters the 
  #' robustness test occurred.
  #' 
  #' In addition, the dilution_proportions are converted from a simple vector
  #' to a named list, and filepaths for logs are automatically generated if 
  #' needed. 
  #' 
  #' All if these arguments are then funneled into resource_Robustness().
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
  #' "variable")? Only takes a string as input, not a vector of a string.
  #' 
  #' @param preserve_topology When diluting, two methods are implemented.
  #' random_Dilute makes only a small effort to preserve the topology of the 
  #' rows that are diluted from the resource, preserve_Dilute makes a far 
  #' greater effort. Choose TRUE for preserve_Dilute and FALSE for 
  #' random_Dilute.
  #' 
  #' @param dilution_props A sequence of numerics (0-1) that indicate which 
  #' proportions to dilute the resource with. For example c(0.1, 0.2, 0.3) would
  #' compare the top ranked CCI's using undiluted OmniPath compared to OmniPath 
  #' with 10 % of its rows diluted, undiluted vs 20 % diluted, and undiluted vs
  #' 30 % diluted.
  #' 
  #' @param outputs Which outputs of the calculation would you like to return? 
  #' By default, all the method results, resources used, top ranked CCIs, 
  #' analysis of top ranks, script parameters and the chosen testdata are 
  #' returned in a list. Construct an atomic vector using all or some of 
  #' "liana_results_OP", "resources_OP", "top_ranks_OP", "top_ranks_analysis",
  #' "metadata", and "testdata" to tell the script which outputs should go in 
  #' the returned list. It's probably best to run this once with testdata to 
  #' better understand which list element holds which data, and then pair it 
  #' down to what is needed.
  #' 
  #' @param methods_vector Which methods should the function run? Choose from
  #' "call_connectome", "squidpy", "call_natmi", "call_italk", "call_sca" and
  #' "cellchat". Supply the argument in the form of e.g. 
  #' c("call_conncectome, "call_italk").
  #' 
  #' @param number_ranks A named list. Each item is named after a method and is 
  #' equal to the number of top interactions considered relevant for that method.
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

             feature_type      = "variable", # "generic" or "variable"
             preserve_topology = FALSE,      # TRUE or FALSE
             dilution_props    = c(seq(0.40, 1.00, 0.40)),
             
             outputs = c(
               "liana_results_OP"  ,
               "resources_OP"      ,
               "top_ranks_OP"      ,
               "top_ranks_analysis",
               "runtime"          ,
               "testdata"
             ),
             
             methods_vector = c('call_connectome' ,
                                #'squidpy'         ,
                                'call_natmi'      ,
                                'call_italk'      ,
                                'call_sca'        ,
                                'cellchat'),
             
             number_ranks = list(
               "call_connectome" = 20,
               "squidpy"         = 20,
               "call_natmi"      = 20,
               "call_italk"      = 20,
               "call_sca"        = 20,
               "cellchat"        = 20
             ),
             
             cellchat_nperms = 10,       # default 100 for real data
             sink_output     = FALSE,    # TRUE or FALSE
             liana_warnings  = "divert", # TRUE, FALSE, or "divert"
             
             time_of_run) {
      
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
        
        
        # If necessary, create a log name, store it in script params, then remove
        # the leftover clutter
        if (sink_output == TRUE) {
          # The file name includes many script_params in it to be informative and
          # unique.
          sink_logfile <-
            auto_file_Name(
              prefix = "Outputs/Resource_Dilution/Logs/Complete_Log_",
              suffix =  ".txt",
              dilution_params = formals(wrap_resource_Robustness),
              meta_params     = formals(summarise_Metadata),
              testdata_type   = formals(extract_Testdata)$testdata_type,
              time_of_run     = time_of_run
            )
          
          
        }
        
        # If necessary, create a log name, store it in script params, then remove
        # the leftover clutter
        if (liana_warnings == "divert") {
          # The file name includes many script_params in it to be informative and
          # unique.
          warning_logfile <-
            auto_file_Name(
              prefix = "Outputs/Resource_Dilution/Logs/LIANA_warnings_",
              suffix =  ".txt",
              dilution_params = formals(wrap_resource_Robustness),
              meta_params     = formals(summarise_Metadata),
              testdata_type   = formals(extract_Testdata)$testdata_type,
              time_of_run     = time_of_run
            )
          
          
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
        outputs           = outputs,
        
        methods_vector    = methods_vector,
        cellchat_nperms   = cellchat_nperms,
        sink_output       = sink_output,
        liana_warnings    = liana_warnings,
        
        sink_logfile      = sink_logile,
        warning_logfile   = warning_logfile
      )
      
      
      
    }
  
  
}

# extract_Testdata()
{
  #' Helper function that gets a specific seurat object from the outputs folder
  #' 
  #' @param testdata_type As a string. Which testdata should be retrieved? 
  #' Either "seurat_pbmc" or "liana_test". Seurat_pbmc is the data set used in 
  #' the seurat tutorial, while liana_test is the testdata that comes with 
  #' LIANA++, and is a small subset of seurat_pbmc.
  #' 
  #' @return A seurat object loaded from the outputs folder or liana apckage.
  
  
  extract_Testdata <- function(testdata_type = "liana_test") {
    
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
    
  }
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
  #' This value is also stored in the metadata.
  #' 
  #' @param trial_run Is this a trial run of the iterator or serious results?
  #' Takes a boolean. If this is a trial run, the save file names, logs and plot
  #' captions will reflect this.
  #' 
  #' @return Returns a list of metadata and parameters.
  
  
  
  summarise_Metadata <- function(runtime,
                                 time_of_run,
                                 dilution_params,
                                 testdata_type,
                                 master_seed_list,
                                 save_results = TRUE,
                                 trial_run    = TRUE) {
    
    # Summarize the metadata parameters
    meta_params <- list(
      "time_of_run"  = time_of_run,
      "save_results" = save_results,
      "trial_run"    = trial_run
    )
    
    # Put all the parameters in a list
    metadata <- list(
      "runtime"         = runtime,
      "dilution_params" = append(dilution_params,
                                 master_seed_list),
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
          dilution_params = formals(wrap_resource_Robustness),
          meta_params     = formals(summarise_Metadata),
          testdata_type   = formals(extract_Testdata)$testdata_type,
          time_of_run     = time_of_run
        )
      
      metadata[["line_plot_png_name"]] <-
        auto_file_Name(
          prefix = "Lineplot_Resource_Dilution_",
          suffix = ".png",
          dilution_params = formals(wrap_resource_Robustness),
          meta_params     = formals(summarise_Metadata),
          testdata_type   = formals(extract_Testdata)$testdata_type,
          time_of_run     = time_of_run
        )
      
      metadata[["env_save_path"]] <-
        auto_file_Name(
          prefix = "Outputs/DilutionEnv_",
          suffix = ".RData",
          dilution_params = formals(wrap_resource_Robustness),
          meta_params     = formals(summarise_Metadata),
          testdata_type   = formals(extract_Testdata)$testdata_type,
          time_of_run     = time_of_run
        )
    }
    
    
    # If the resource_Robustness() output was sunk and logged, append the 
    # file name of the log to the metadata.
    if (dilution_params$sink_output == TRUE) {
      # The file name includes many script_params in it to be informative and
      # unique.
      metadata[["sink_logfile"]] <-
        auto_file_Name(
          prefix = "Outputs/Logs/Complete_Log_",
          suffix =  ".txt",
          dilution_params = formals(wrap_resource_Robustness),
          meta_params     = formals(summarise_Metadata),
          testdata_type   = formals(extract_Testdata)$testdata_type,
          time_of_run     = time_of_run
        )
      
      
    }
    
    # If a warnings log was created for resource_Robustness(), append the file
    # name of the log to the metadata.
    if (dilution_params$liana_warnings == "divert") {
      # The file name includes many script_params in it to be informative and
      # unique.
      metadata[["warning_logfile"]] <-
        auto_file_Name(
          prefix = "Outputs/Logs/LIANA_warnings_",
          suffix =  ".txt",
          dilution_params = formals(wrap_resource_Robustness),
          meta_params     = formals(summarise_Metadata),
          testdata_type   = formals(extract_Testdata)$testdata_type,
          time_of_run     = time_of_run
        )
      
      
    }
    
    
    # return the metadata.
    return(metadata)
    
  }
  
  
  
}

# create_Params()
{
  #' Produces parameters for Robustness_Iterator.R that´need to be stored as
  #' objects in the environment
  #'
  #' @description Most of the time, when a function in the iterator needs 
  #' information about how the iterator was run, its enough to pass it the 
  #' formals() of the relevant function. However, there are some iterator 
  #' parameters that can't be stored as function defaults for one reason or 
  #' another. This function handles these cases by creating the relevant 
  #' parameters as objects in the environment, where they can be referenced.
  #' 
  #' @param master_seeds As a vector of numerics. How many times should 
  #' resource_Robustness be iterated? resource_Robustness is run once for each
  #' numeric in this vector, while using the numeric as the master_seed argument
  #' for resource_Robustness. This would be tricky to store as a function
  #' argument default because you would need to manually format and name it each
  #' time you changed its lenghth.
  #' 
  #' @return This function returns a formatted named list of seeds derived from
  #' the input master_seeds.
  #' 
  #' It also eturns a char called time_of_run. This is the Sys.time() formatted 
  #' as a char tag that can be appended to file save names, marking them all 
  #' uniquely, and making the relationship between logs, plots and .RData files 
  #' directly apparent (all the files that were created at the same time go
  #' together). This couldn't be stored in a functions argument defaults because
  #' it needs to be dynamically generated, unless the user is supposed to type
  #' in the Sys.time() at every run.
  #' 
  #' In sum, this function simply loads parameters into the environment that 
  #' would be too inconvenient to call with formals().
  
  
  create_Params <- function(master_seeds = 1:2) {
    
    # Format Sys.time() to not contain characters that are bad to have in
    # save file names.
    time_of_run <-  Sys.time() %>%
      as.character()    %>% 
      gsub(':', '-', .) %>% # save files can't have colons
      gsub(' ', '_', .)     # save files shouldn't have spaces
    
    
    
    
    # Format master_seeds as a list, this makes it easy to lapply over it later 
    master_seed_list <- as.list(master_seeds)
    
    # By naming each seed we can use this to label data conveniently as the
    # lapply iteration is happening
    names(master_seed_list) <-
      map(master_seed_list, function(seed) {
        
        # Name each element of master_seed_list appropriately
        str_glue("Seed_", seed)
        
      })
    
    
    
    
    # Summarize and return our two parameters.
    parameters <- list("master_seed_list" = master_seed_list,
                       "time_of_run"       = time_of_run)
    
    return(parameters)
  }
  
  
}