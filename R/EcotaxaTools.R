
#' EcotaxaTools
#'
#' A routine to compute biovolumes and save updated csv files. Save also summary tables for each taxonomic group and figures to check the data.
#' 1. Choose a directory and create a new one inside to save results.
#' 2. Use function check_metadata to allow metadata edit.
#' 3. Compute biovolume and BSS tables with chosen metadata.
#' 4. Save results and summary tables (AB or BV for each taxa by sample).
#' 5. Graphical output : metadata, project and samples.
#'
#' @return Biovolumes, dataframes and summary tables, figures to check the data.
#' @export
#'
#' @examples EcotaxaTools()
EcotaxaTools <- function() {
  options(dplyr.summarise.inform = FALSE)
  # Choose your instrument
  instru <- dlg_list(c("PlanktoScope", "FlowCam","ZooScan","IFCB"), title="Instrument")$res

  # Files 
  start <- dlg_message("Please select the data directory. All the tsv files in this directory will be processed.", type="okcancel")$res
  if (start=="cancel") stop("Canceled.")

if (rstudioapi::isAvailable()) {
  mainpath <- rstudioapi::selectDirectory()
} else {
  if (Sys.info()["sysname"] == "Windows") {
    mainpath <- svDialogs::dlg_dir()
  } else {
    mainpath <- svDialogs::dlg_dir_choose()
  }
}

if (!is.null(mainpath) && mainpath != "") {
  path <- list.files(mainpath, pattern = "\\.tsv$", full.names = TRUE) 
} else {
  message("No valid directory selected.")
}

  # Create a new directory for created files
  if (!file.exists(file.path(mainpath, paste0("EcoTaxa_data_analysis_",format(Sys.time(), "%Y-%m-%d %H%M"))))) {
    dir.create(file.path(mainpath, paste0("EcoTaxa_data_analysis_",format(Sys.time(), "%Y-%m-%d %H%M"))))
  }
  output <- file.path(mainpath, paste0("EcoTaxa_data_analysis_",format(Sys.time(), "%Y-%m-%d %H%M")))

  # COMPUTE DATA
  # ------------------------------------------------------------------------------

  final_dataset <- ecotaxa_tools(path, output, instru)
 
  # GRAPHICAL OUPUTS
  # ------------------------------------------------------------------------------
  # Create directory
  if (!file.exists(file.path(output,"graph"))) {
    dir.create(file.path(output,"graph"))
  }
  path.graph <- file.path(output,"graph")

  # for the original metadata (saved by check_metadata) - Commented as not very informative
  #pdf(file.path(path.graph, "original_metadata.pdf"),width = 10, paper="a4r")
  #graph.metadata(read_csv2(file.path(output, "metadata", "original_metadata.csv")))
  #dev.off()

  #### for the edited metadata (using the returned 'processed_metadata')
  
  graph.metadata(final_dataset,path.graph)
  
  #### for the project - Create directory + Plot graphics
  if (!file.exists(file.path(path.graph,"global raw analysis"))) {
    dir.create(file.path(path.graph,"global raw analysis"))
  }
  path.graph_project <- file.path(path.graph,"global raw analysis")
  print("Creating graphical output for the whole project")
  
  graph.project(final_dataset, path.graph_project)
  
  #### for each sample - Create directory + Plot graphic per sample
  if (!file.exists(file.path(path.graph,"raw analysis per sample"))) {
    dir.create(file.path(path.graph,"raw analysis per sample"))
  }
  path.graph_sample <- file.path(path.graph,"raw analysis per sample")
  
  for (i in unique(final_dataset$sample_id)) {
    final_dataset %>% filter(sample_id==i) %>% graph.sample() %>%
      ggsave(filename=file.path(path.graph_sample, paste0(i,".jpg")),
             width=297, height=210, units = "mm")
  }
  
  # MULTIVARIATE ANALYSIS
  # ------------------------------------------------------------------------------
   # Ask for multivariate analysis
   yesno <- dlg_message("Do you want to continue the pipeline with multivariate analysis (still in progress)", type="yesno")$res
   
   if(yesno=="yes") {
     
     #### Create  output directory + start analysis
     if (!file.exists(file.path(output,"multivariates analysis"))) {
       dir.create(file.path(output,"multivariates analysis"))
     }
     path.analysis <- file.path(output,"multivariates analysis")
     print("Creating graphical output for multivariates analysis")
     
     multivariates_analysis(final_dataset, path.analysis)
     }else{
       print("End of the script.")
       }
  print("End of the script.")
}
