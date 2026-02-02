
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
  # Check metadata
  processed_metadata <- check_metadata(path, output, instru)

  # Compute biovolumes and BSS summary (warning : not normalized by size class)
  yesno <- dlg_message("IMPORTANT: Do you want to use the edited metadata ? If not you can select the original metadata or another metadata table", type="yesno")$res

  if(yesno=="yes") {
    metadata <- processed_metadata
    bss <- lapply(path, function(x) BSS_table(compute_bv(x, output, metadata, instru))) %>% bind_rows()
  } else {
    metadata <- file.choose() %>% read_csv2()
    bss <- lapply(path, function(x) BSS_table(compute_bv(x, output, metadata, instru))) %>% bind_rows()
  }


  # SUMMARY DATA AND SAVING TABLES
  # ------------------------------------------------------------------------------
  # Create directory
  if (!file.exists(file.path(output,"summary"))) {
    dir.create(file.path(output,"summary"))
  }
  path.summary <- file.path(output,"summary")

  # Global summary table without size class
  res <- bss %>% group_by(sample_id, object_annotation_category, type) %>%
    summarize(AB = sum(AB, na.rm=T),
              BV= sum(BV, na.rm=T))

  # Loop to save tables for Biovolume (Elli/Plain/Riddled)
  for (bv.type in unique(res$type)) {
    t <- res %>% filter(type==bv.type) %>% select(-type, -AB) %>%
      pivot_wider(names_from = sample_id,
                  values_from = BV)
    t$Total <- rowSums(t[-1], na.rm=T) # sum by sp
    tot <- t[1,] %>% mutate(object_annotation_category="Total") # sum by sample
    tot[-1] <- t(colSums(t[-1], na.rm=T))
    t <- rbind(t, tot)
    write_csv2(t, file.path(path.summary, paste0("Biovolume_",bv.type,".csv")))
  }

  # Summary table for Abundance (same for elli, plain, etc.)
  t <- res %>% filter(type=="plain") %>% select(-type, -BV) %>%
    pivot_wider(names_from = sample_id,
                values_from = AB)
  t$Total <- rowSums(t[-1], na.rm=T) # sum by sp
  tot <- t[1,] %>% mutate(object_annotation_category="Total") # sum by sample
  tot[-1] <- t(colSums(t[-1], na.rm=T))
  t <- rbind(t, tot)
  write_csv2(t, file.path(path.summary, "Abundance.csv"))

  # For the taxonomy
  taxo <- add.taxo(unique(bss$object_annotation_hierarchy)) %>% add.trophiclvl(., output)

  # Replacing all the NA in case the original metadata was selected
  final_metadata <- metadata
  final_metadata[is.na(final_metadata)] <- 1

  # Saving tables
  write_csv2(final_metadata, file.path(path.summary, "metadata_used.csv"))
  write_csv2(bss, file.path(path.summary, "BSS.csv"))
  write_csv2(res, file.path(path.summary, "summary_all.csv"))
  write_csv2(taxo, file.path(path.summary, "taxo.csv"))


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
  
  graph.metadata(final_metadata)
  
 
  #### for the project - Create directory + Plot graphics
  if (!file.exists(file.path(path.graph,"global raw analysis"))) {
    dir.create(file.path(path.graph,"global raw analysis"))
  }
  path.graph_project <- file.path(path.graph,"global raw analysis")
                  
  graph.project(bss, final_metadata, taxo)
  
  #### for each sample - Create directory + Plot graphic per sample
  if (!file.exists(file.path(path.graph,"raw analysis per sample"))) {
    dir.create(file.path(path.graph,"raw analysis per sample"))
  }
  path.graph_sample <- file.path(path.graph,"raw analysis per sample")
  
  for (i in unique(bss$sample_id)) {
    bss %>% filter(sample_id==i) %>% graph.sample(final_metadata, taxo) %>%
      ggsave(filename=file.path(path.graph, paste0(i,".jpg")),
             width=297, height=210, units = "mm")
  }



  print("End of the script.")
}
