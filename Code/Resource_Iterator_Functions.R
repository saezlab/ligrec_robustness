
wrap_resource_Robustness <- 
  function(master_seed,
           
           testdata_type     = c("liana_test"), # "liana_test" or "seurat_pbmc"
           feature_type      = c("variable"),   # "generic" or "variable"
           preserve_topology = FALSE,           # TRUE or FALSE
           dilution_props    = c(seq(0.40, 1.00, 0.40)),
           
           outputs = c(
             "liana_results_OP"  ,
             "resources_OP"      ,
             "top_ranks_OP"      ,
             "top_ranks_analysis",
             "metadata"          #,
           #  "testdata"
           ),
           
           methods_vector = c(
              'call_connectome' ,
           #  'call_squidpy'    ,
           #  'call_natmi'      ,
              'call_italk'      ,
              'call_sca'        #,
           #  'cellchat'
           ), 
           
           number_ranks = list(
             "call_connectome" = 20,
             "call_squidpy"    = 20,
             "call_natmi"      = 20,
             "call_italk"      = 20,
             "call_sca"        = 20,
             "cellchat"        = 20
           ),
           
           cellchat_nperms = 10,      # default 100
           sink_output     = FALSE,   # TRUE or FALSE
           liana_warnings  = "divert" # TRUE, FALSE, or "divert"
           
           ) {
    
  print(master_seed)
  
  ## Process Dilution Proportions List
  {
    # By naming each dilution proportion we can use this to label data easily
    # By formatting proportions as a list we can lapply over them conveniently
    
    # intitialize dilution names
    dilution_names <- c()
    
    # Generate a dilution name for each proportion
    for (i in dilution_props) {
      dilution_names <- 
        c(dilution_names, str_glue("OmniPath_", as.character(i*100)))
    }
    
    # Convert dilution.proprs to a list and name it, creating a named list
    dilution_props <- as.list(dilution_props)
    names(dilution_props) <- dilution_names
    
    #remove clutter
    rm(dilution_names, i)
  }
  
  ## Define File Paths for Logs
  {
    sink_logfile <- ""
    warning_logfile <- ""
    
    
    # If they are required, we auto generate log filepaths here. The filepaths
    # have the current time imminetnly before the lapply in them. Each lapply
    # should be assigned to a unique log this way that has all iterations in
    # it. 
    
    # Define time of run and remove characters problematic for file names
    time_of_run <-  Sys.time() %>%
      as.character()       %>%
      gsub(':', '-', .)    %>% 
      gsub(' ', '_', .)
    
    # If necessary, create a log name, store it in script params, then remove
    # the leftover clutter
    if(sink_output == TRUE) {
      # The file name includes many script_params in it to be informative and
      # unique.
      sink_logfile <- str_glue(
        "Outputs/Logs/Complete_Log_",
        run_mode,
        "_",
        testdata_type,
        "_top",
        as.character(median(unlist(number_ranks))),
        "_res",
        as.character(length(dilution_props)),
        "_",
        feature_type,
        "_dil_at_",
        time_of_run,
        ".txt"
      )
      

    }
    
    # If necessary, create a log name, store it in script params, then remove
    # the leftover clutter
    if(liana_warnings == "divert") {
      # The file name includes many script_params in it to be informative and
      # unique.
      warning_logfile <- str_glue(
        "Outputs/Logs/LIANA_warnings_",
        run_mode,
        "_",
        testdata_type,
        "_top",
        as.character(median(unlist(number_ranks))),
        "_res",
        as.character(length(dilution_props)),
        "_",
        feature_type,
        "_dil_at_",
        time_of_run,
        ".txt"
      )
      
    }
  }
  
  
    
  resource_Robustness(
    master_seed       = master_seed, 
    
    testdata_type     = testdata_type,
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






# wrap_resource_Robustness()


# Which scRNA data set do you want to use?


# Should dilution be with any genes from testdata or most variable ones?


# Preserve topology after dilution or not?
# TRUE = preserve_Dilute(), FALSE = random_Dilute()


# Which proportions should the resources that are analysed have?
# have at least two else the plotting fails


# Which Outputs from resource_Robustness() do you want? Choose from:


# Which methods should resource_Robustness() use? Choose from:
# Squidpy doesn't work on windows.

# Which top n of interactions should be considered top-ranked per method?


# number of cellchat permutations, default 100

# Should the entire output be redirected to a log file?

# Should liana warnings be visible in output? 
