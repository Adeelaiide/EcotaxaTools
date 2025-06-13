
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
  # Definition de l'instrument
  instru <- dlg_list(c("PlanktoScope", "FlowCam","ZooScan","IFCB", "UVP"), title="Instrument")$res

  # Fichiers (faire option selection seulement certains fichiers?)
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
  check_metadata(path, output)

  
  # Set common color palette
  plankton_groups_colors <- c("#709699", #cyanobacteria
               "#F2F2F2", #detritus
               "#FAFABA", #other
               "#C9C4FF", #ciliophora
               "#B5D493", #dinoflagellata
               "#FAD4EB", #rhizaria
               "#1A9C75", #bacillariophyta
               "#E3E699", #dictyochophyceae
               "#EDB022", #crustacea
               "#E88F6B", #copepoda
               "#B0D6E0", #chaetognatha
               "#3B638A", #tunicata
               "#6999C7", #cnidaria
               "#C68181", #mollusca
               "#668F3B", #coccolithophyceae
               "#FFD3CF", #other_unidentified
               "#0073BD" #plastics
               )
  
  names(plankton_groups_colors)<- c("cyanobacteria","detritus","other","ciliophora","dinoflagellata",
                     "rhizaria","bacillariophyta","dictyochophyceae","crustacea",
                     "copepoda","chaetognatha","tunicata","cnidaria","mollusca",
                     "coccolithophyceae","other_unidentified","plastics")
  plankton_groups_colScale <- scale_colour_manual(name = "taxonomic group",values = plankton_groups_colors)
  plankton_groups_colFill <- scale_fill_manual(name = "taxonomic group",values = plankton_groups_colors)

  # Compute biovolumes and BSS summary (warning : not normalized by size class)
  yesno <- dlg_message("IMPORTANT: Do you want to select the edited metadata or another metadata table ? If not the original metadata will be used.", type="yesno")$res

  if(yesno=="yes") {
    metadata <- file.choose() %>% read_csv2()
    bss <- lapply(path, function(x) BSS_table(compute_bv(x, output, metadata))) %>% bind_rows()
  } else {
    bss <- lapply(path, function(x) BSS_table(compute_bv(x, output))) %>% bind_rows()
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

  # Unite taxonomiques
  taxo <- add.taxo(unique(bss$object_annotation_hierarchy)) %>% add.trophiclvl(., output)

  # Search for metadata if not loaded already
  if(!exists("metadata")){
    metadata <- read_csv2(file.path(output, "metadata", "original_metadata.csv"))
  }
  metadata[is.na(metadata)] <- 1

  # Saving tables
  write_csv2(metadata, file.path(path.summary, "metadata_used.csv"))
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

  # for the original metadata
  pdf(file.path(path.graph, "original_metadata.pdf"), paper="a4r")
  graph.metadata(read_csv2(file.path(output, "metadata", "original_metadata.csv")))
  dev.off()

  # for the edited metadata
  pdf(file.path(path.graph, "metadata.pdf"), paper="a4r")
  graph.metadata(metadata)
  dev.off()
  # system(paste0('open "', file.path(path.graph, "metadata.pdf"), '"'))

  # for the project
  pdf(file.path(path.graph, "graph_project.pdf"), paper="a4r")
  graph.project(bss, metadata, taxo)
  dev.off()

  # for each sample
  for (i in unique(bss$sample_id)) {
    bss %>% filter(sample_id==i) %>% graph.sample(metadata, taxo) %>%
      ggsave(filename=file.path(path.graph, paste0(i,".jpg")),
             width=297, height=210, units = "mm")
  }



  print("End of the script.")
}
