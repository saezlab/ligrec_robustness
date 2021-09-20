#------------------------------------------------------------------------------#
# 0. Introduction and Goals ----------------------------------------------------
{
  # This is the robustness iterator. It funnels data into a robustness test for
  # various dilutions (resource_Robustness()), repeats this test many times, and
  # collates the resulting information and presents it to the user. The code to
  # accomplish this is presented in a  wrapper with default arguments below.
  
  # Running parameters can be left a their defaults or altered when the function
  # is called. To learn more, check the Run_Iterator.R script.
}




#------------------------------------------------------------------------------#
# 1. Define wrap_robustness_Iterator() -----------------------------------------
{
  
#' resource_Robustness() Iterator
#' 
#' @description This function iterates the resource_Robustness() function in
#' RD_Robustness_Evaluator.R. It compares CCI predictiuions between methods 
#' running on undiluted resources and resources diluted to varying degrees.
#' 
#' This function feeds resource_Robusntess() the defaults it needs, runs it many 
#' times and collates the results. 
#' 
#' @return A list object that contains the resource_Robustness() results 
#' separated by iteration, the collated top_rank_overlap data, two plots 
#' describing that data and metadata statistics. 
#' 
#' 
#' 
#' @param number_seeds As an integer. How many iterations of the robustness 
#' evaluator should be run? This is the most important value for a 
#' representative analysis, but also directly influences the length it takes for
#'  the function to execute.
#'  
#' @param testdata_type As a string. Which testdata should be used in the 
#' analysis? Either "seurat_pbmc" or "liana_test". Seurat_pbmc is the data set 
#' used in the seurat tutorial, while liana_test is the testdata that comes with 
#' LIANA++, and is a small subset of seurat_pbmc. 
#' 
#' @param feature_type Should dilution occur with all the genes profiled in
#' testdata (choose "generic") or with the most variable features (choose 
#' "variable")? As a string.
#' 
#' @param preserve_topology When diluting, two methods are implemented.
#' random_Dilute makes only a small effort to preserve the topology of the rows
#' that are diluted from resource, preserve_Dilute makes a far greater effort.
#' Choose TRUE for preserve_Dilute and FALSE for random_Dilute.
#' 
#' @param dilution_props A sequence of doubles (0-1) that indicate which 
#' proportions to dilute the resource with. For example c(0.1, 0.2, 0.3) would
#' compare the top ranked CCI's using undiluted OmniPath compared to OmniPath 
#' with 10 % of its rows diluted, undiluted vs 20 % diluted, and undiluted vs
#' 30 % diluted.
#' 
#' @param number_ranks A named list. Each item is named after a method and is 
#' equal to the number of top interactions considered relevant for that method. 
#' 
#' For example, item one on the list may be called call_connectome and be equal 
#' to 500. This would signal to the function that for call_connectome, the top 
#' 500 CCI's are considered relevant and that these 500 are the ones that are to
#' be compared between the dilutions. Usually you consider the same number of
#' top_ranks for each method relevant.
#' 
#' @param methods_vector Which methods should the function run? Choose from
#' "call_connectome", "squidpy", "call_natmi", "call_italk", "call_sca" and
#' "cellchat". Supply the argument in the form of e.g. 
#' c("call_conncectome, "call_italk"). This argument is also very impactful for 
#' the runtime of this function. For example on widnows, natmi and cellchat are
#' the slowest methods.
#' 
#' @param sink_otuput TRUE or FALSE. Should resource_Robustness() save a full 
#' log of the Console Output to the log folder? Warnings and messages will not 
#' be visible in the console output if this is enabled, so unless there is a 
#' reason why such a record is necessary this option is not recommended.
#' 
#' @param liana_warnings Either TRUE, FALSE or "divert". A less extreme 
#' alternative to sinking the output. Should the warnings from LIANA++, which 
#' are often repetitive and unhelpful, be either suppressed or alternatively 
#' diverted to the log folder? When these types of warning are left in, they can
#' often displace valuable warnings. Be careful with this  setting, as 
#' suppressing warnings is obviously risky.
#' 
#' @param save_results Either TRUE or FALSE. Should the plots and iterator
#' results objects be saved to the outputs folder with an automatically
#' generated file name? Recommended for keeping the results safe.
#' 
#' @param trial_run Is this a trial run of the iterator or serious results?
#' Takes a boolean. If this is a trial run, the save file names, logs and plot
#' captions will reflect this. Only takes a boolean as a default, not a vector
#' of a boolean.
#' 
#' @param cellchat_nperms Cellchat is one of the slower methods, for test runs
#' it may be useful to set this parameter to 10 to speed up the analysis. 
#' Unless you're using this function for a trial run, this should be left at its
#' default.
#' 
#' @param bundled_outputs Which outputs of resource_Robustness() would you like 
#' to return? Outputs are returned bundled up in a list. By default, only the 
#' analysis of top_ranks and runtime data is returned, but more infomation can 
#' be returned, if desired (though usually uneccesary). If the 
#' top_ranks_analysis and runtime aren't returned, the script wont work.
#' 
#' To specify your outputs, construct an atomic vector using all or some of 
#' "liana_results_OP", "resources_OP", "top_ranks_OP", "top_ranks_analysis", 
#' "runtime", and "testdata" to tell the script which outputs should go in the 
#' output bundle.
#' 
#' "liana_results_OP": All the predicted CCI's from LIANA++ for each method and
#' dilution stage.
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
#' @param master_outputs Which of the iterator results do you want in the final
#' output? By default, all the data generated is returned in a list, but in
#' unusual scenarios you could subset this list with a vector passed to this 
#' argument. Each of the strings below is the name of a piece of data that can
#' be retruned by adding it to the vector you pass here.
#' 
#' "collated_top_ranks_overlap": The combined information on the overlap of 
#' predicted top ranks between dilution is returned.
#' 
#' "plot_box" and "plot_line": The two different plots that describe the 
#' collated_top_ranks_overlap data.
#' 
#' "collated_robustness_results": The bundle of information that was the
#'  resource_Robusntess() output, though more sensibly structured. Check the 
#'  bundled_outputs parameter.


wrap_robustness_Iterator <- 
  function(
    number_seeds      = 10,            # how many seeds should we iterate over
    testdata_type     = "seurat_pbmc", # seurat_pbmc or liana_test
    feature_type      = "variable",    # "generic" or "variable"
    preserve_topology = FALSE,         # TRUE or FALSE
    dilution_props    = c(seq(0.05, 1.00, 0.05)),
    
    number_ranks = list(
      "call_connectome" = 500,
      "call_natmi"      = 500,
      "call_italk"      = 500,
      "call_sca"        = 500,
      "cellchat"        = 500,
      "squidpy"         = 500
    ),
    
    methods_vector = c('call_connectome' ,
                       'call_natmi'      ,
                       'call_italk'      ,
                       'call_sca'        ,
                       'cellchat'        ,
                       'squidpy'
    ),
    
    sink_output     = FALSE,    # TRUE or FALSE
    liana_warnings  = "divert", # TRUE, FALSE, or "divert"
    
    save_results    = TRUE,
    trial_run       = FALSE,
    
    
    
    
    cellchat_nperms = 100,       # default 100 for real data
    
    bundled_outputs = c(
      "top_ranks_analysis",
      "runtime"          
    ),
    
    master_outputs = c(
      "collated_top_ranks_overlap",
      "plot_box",
      "plot_line",
      "collated_robustness_results",
      "metadata"
    )                                  
  ) {
    
    
  
  
  #----------------------------------------------------------------------------#
  # 1.1 Generate Parameters  --------------------------------------------------- 
  {
    # Retrieve either seurat_pbmc or liana_test data
    testdata <- extract_Testdata(testdata_type = testdata_type)
    
    
    
    # Format a named list of seeds, it contains as many seeds as the user 
    # specified, from 1 to n, and each entry has an appropriate name, "Seed_n"
    master_seed_list <- as.list(1:number_seeds)
    
    names(master_seed_list) <-
      map(master_seed_list, function(seed) {
        
        # Name each element of master_seed_list appropriately
        str_glue("Seed_", seed)
        
      })
    
    # Format the Sys.time() of the run to not contain characters that are bad to
    # have in save file names. We will later use this to tag file names and
    # plots so they can be grouped according to run, and all have unique names.
    time_of_run <-  Sys.time() %>%
      as.character()    %>% 
      gsub(':', '-', .) %>% # save files can't have colons
      gsub(' ', '_', .) %>% # save files shouldn't have spaces
      str_sub(1 , nchar(.) - 3) # the code never runs in under a minute, so the 
    # number of seconds isn't valuable information.
    
  }  
  
  
  
  #----------------------------------------------------------------------------#
  # 1.2 Iterate resource_Robustness() ------------------------------------------
  {
    # resource_Robustness is a function that performs a single test of 
    # robustness by comparing unmodified LIANA predictions with ones run on 
    # diluted OP resources. Since there is randomness to resource dilution, we 
    # iterate over multiple starting seeds in master_seed_list, and will collate 
    # these samples later.
    
    # Apply resource_Robustness() wrapper, provide every argument but 
    # master_seed.
    # To modify the defaults of the wrapper, go to Iterator_Params.R
    collated_robustness_results <- lapply(master_seed_list,
                                          wrap_resource_Robustness,
                                          
                                          testdata          = testdata,
                                          feature_type      = feature_type,
                                          preserve_topology = preserve_topology, 
                                          dilution_props    = dilution_props,
                                          number_ranks      = number_ranks,
                                          methods_vector    = methods_vector,
                                          
                                          bundled_outputs = bundled_outputs,
                                          cellchat_nperms = cellchat_nperms, 
                                          
                                          sink_output     = sink_output,     
                                          liana_warnings  = liana_warnings,
                                          trial_run       = trial_run,
                                          time_of_run     = time_of_run,
                                          testdata_type   = testdata_type)
    
    # We don't need the testdata after this point.
    rm(testdata)
    
    
  }
  
  
  
  #----------------------------------------------------------------------------#
  # 1.3 Reformatting Results ---------------------------------------------------
  {
    # In this segment we extract the data from the results object, which is 
    # poorly formatted by default, and put it into a more appropriate hierarchy.
    # We then extract the most important information from it for visualization.
    
    # We reformat the results so they are grouped more intuitively
    collated_robustness_results <-
      reformat_Results(results = collated_robustness_results)
    
    # We extract the top_ranks_overlap data and turn it into a convenient tibble
    collated_top_ranks_overlap <-
      extract_top_ranks(collated_robustness_results)
    
    
  }
  
  
  
  #----------------------------------------------------------------------------#
  # 1.4 Plotting of Collated Results -------------------------------------------
  {
    # Here we visualize the overlap between top ranked CCI predictions as they
    # change with the dilution of OmniPath. once as a boxplot, and once as a
    # scatter/line plot.
    
    
    
    # We reformat the collated_top_ranks_overlap tibble so its more suitable for
    # plotting
    
    # The dilution proportion and overlap are clearer in percentage
    # NAs can't be displayed in the plot anyway and cause uneccesary warnings
    # Rename the methods from the LIANA++ internal string to their official name
    tr_overlap_for_plot <-  collated_top_ranks_overlap  %>%
      as.data.frame()                             %>%
      mutate(dilution_prop = dilution_prop * 100) %>% # proportion to percent
      mutate(Overlap       = Overlap       * 100) %>% # proportion to percent
      as_tibble()                                 %>%
      drop_na()                                   %>% # no NAs
      mutate("Method" = recode(Method,
                               "call_connectome" = "Connectome",
                               "squidpy"         = "CellPhoneDB",
                               "call_natmi"      = "NATMI",
                               "call_italk"      = "iTALK",
                               "call_sca"        = "SingleCellSignalR",
                               "cellchat"        = "CellChat")) # renaming
    
    # To directly be able to associate the box plot with the settings that
    # produced it, we automatically generate a plot description
    plotting_caption <-
      auto_plot_Description(
        tr_overlap_for_plot,
        
        preserve_topology  = preserve_topology,
        testdata_type      = testdata_type,
        feature_type       = feature_type,
        number_ranks       = number_ranks,
        time_of_run        = time_of_run,
        trial_run          = trial_run)
    
    
    
    
    
    
    # Generate our plots with functions.
    plot_line <- overlap_line_Plot(tr_overlap_for_plot,
                                   plotting_caption)
    
    plot_box <- overlap_box_Plot(tr_overlap_for_plot,
                                 plotting_caption)
    
    # Print out visualizations
    print(plot_line)
    print(plot_box)
    
    
    
    
    # Removing Clutter
    rm(tr_overlap_for_plot, plotting_caption)
    
    
  }
  
  
  
  #----------------------------------------------------------------------------#
  # 1.5 Capturing Script Metadata ----------------------------------------------
  {
    # In this segment, we summarize the metadata of the resource_Robustness run
    # and the iterator in general. When troubleshooting or reproducing results,
    # this information will be useful to the user.
    
    # rescource_Robustness return the points in time certain processes started,
    # we use this long list of named points in time to create a succinct and
    # informative runtime overview.
    runtime <- collated_robustness_results$runtime %>%
      flatten() %>%
      calculate_Runtime()
    
    # Delete the old unordered runtime from the collated robustness results 
    # bundle
    collated_robustness_results$runtime <- NULL
    
    
    # We then summarise the above information and more metadata into a single 
    # object
    metadata <- summarise_Metadata(number_seeds      = number_seeds,
                                   master_seed_list  = master_seed_list,
                                   testdata_type     = testdata_type,
                                   feature_type      = feature_type, 
                                   preserve_topology = preserve_topology,    
                                   dilution_props    = dilution_props,
                                   number_ranks      = number_ranks ,
                                   methods_vector    = methods_vector,
                                   
                                   sink_output       = sink_output,    
                                   liana_warnings    = liana_warnings,
                                   
                                   cellchat_nperms   = cellchat_nperms,       
                                   bundled_outputs   = bundled_outputs,
                                   master_outputs    = master_outputs,
                                   
                                   
                                   save_results = save_results,
                                   trial_run    = trial_run,
                                   
                                   runtime      = runtime,
                                   time_of_run  = time_of_run)
    
    # Now that these objects are stored in the metadata object, we can remove
    # this clutter from the environment.
    rm(runtime, master_seed_list)
    
    
  }
  
  
  
  #----------------------------------------------------------------------------#
  # 1.6 Packaging Results to return them ---------------------------------------
  {
    #In order to save and return our results we package it in a succinct object
    iterator_results <- 
      list(
        "collated_top_ranks_overlap"  = collated_top_ranks_overlap,
        "plot_box"  = plot_box,
        "plot_line" = plot_line,
        "collated_robustness_results" = collated_robustness_results,
        "metadata"  = metadata
      )
    
    # Filter our results by the master_outputs the user wants to retrieve
    # Usually´all the data is requested so this step doesn't chaneg anything.
    iterator_results <- iterator_results[master_outputs]
    
    # Get rid of clutter we already summarized in different objects
    rm(collated_top_ranks_overlap, collated_robustness_results, metadata)
  }
  
  
  
  #----------------------------------------------------------------------------#
  # 1.7 Saving Results ---------------------------------------------------------
  {
    # In this segment we save the plots and environment to the outputs folder,
    # if it's specified in by the user.
    if (save_results == TRUE) {
      
      save_Results(plot_box  = plot_box,
                   plot_line = plot_line, 
                   iterator_results = iterator_results,
                   
                   preserve_topology  = preserve_topology,
                   testdata_type      = testdata_type,
                   feature_type       = feature_type,
                   number_ranks       = number_ranks,
                   time_of_run        = time_of_run,
                   trial_run          = trial_run)
      
    }
    
    
    # And now return our results to the user
    return(iterator_results)
    
    
  }
  
  
  
} # end of function

  
  
}



