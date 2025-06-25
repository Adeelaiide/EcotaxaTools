#' check_metadata
#'
#' Check the metadata of all samples in a project and allows to edit the metadata in a shiny app.
#' The metadata (original and edited) are saved in the output directory.
#'
#' @param path path to the directory containing all the tsv files
#' @param output path to the output directory where the results are saved
#'
#' @return return a metadata summary of the directory and save the results
#' @export
#'
#' @examples check_metadata(path=project.directory, output=where to save results)
check_metadata <- function(path, output) {

  if (!file.exists(file.path(output,"metadata"))) {
    dir.create(file.path(output,"metadata"))
  }

  meta_file <- function(x) {
    metadata <- read_tsv(x, col_types = list(object_time=col_time(),
                                            object_annotation_time=col_time())) %>%
      group_by(object_label, acq_id) %>% mutate(ghost_id=1:n()) %>% ungroup %>%
      mutate(unique_id = paste(acq_id,sample_operator,ghost_id,
                               object_date,object_time,
                               object_lat,object_lon,sep="_")) %>%
      group_by(unique_id) %>%
      mutate(number_object = n(),
             percentValidated = sum(object_annotation_status=="validated", na.rm=T)/n()*100) %>%
      ungroup() %>%
      mutate(sample_total_volume = ifelse("sample_total_volume" %in% colnames(.), sample_total_volume, NA),
             sample_concentrated_sample_volume = ifelse("sample_concentrated_sample_volume" %in% colnames(.), sample_concentrated_sample_volume, NA),
             sample_dilution_factor = ifelse("sample_dilution_factor" %in% colnames(.), sample_dilution_factor, NA)) %>%
      mutate(sample_dilution_factor = as.numeric(gsub(",", ".",sample_dilution_factor))) %>%
      select(sample_id,
             acq_id,
             unique_id,
             ghost_id,
             object_date,
             object_time,
             object_lat,
             object_lon,
             sample_operator,
             percentValidated,
             number_object,
             acq_nb_frame,
             sample_total_volume,
             sample_concentrated_sample_volume,
             acq_celltype,
             acq_imaged_volume,
             process_pixel,
             sample_dilution_factor) %>%
      distinct() %>% group_by(sample_id) %>% mutate(ghost_id=1:n()) %>% ungroup()

    return(metadata)
  }

  metadata <- do.call("rbind", lapply(path, meta_file))
  
  # Save original
  write_csv2(metadata, file.path(output,"metadata","original_metadata.csv"))
  print("Original metadata saved.")

  # DATA EDIT
  check <- metadata %>% group_by(sample_id) %>% summarize(nb=n())
  check <- max(check$nb, na.rm=T)
  if (check>1){
    dlg_message("Warning ! Some of your samples have more than one acquisition. Be sure that they belong to the same sample. If not, please create a different tsv file for the supplementary acquisition and restart the process.", type="ok")
  }

  time <- format(Sys.time(), "%d-%m-%Y_%H%M")
  dlg_message("You can now edit metadata (click on SYNCHRONIZE and DONE button to update edition). The original .tsv files will not be edited. NA will be replaced by 1. Do not change the unique_id.", type="ok")

  metadata$object_date <- as.character(metadata$object_date)
  metadata$object_time <- as.character(metadata$object_time)
  edited_metadata <- data_edit(metadata, viewer="pane")

  print("Data editing completed.")

  # Arrange by the *edited* date and time before creating sample_num
  edited_metadata <- arrange(edited_metadata, object_date, object_time)

  # Create sample_num based on the *edited and arranged* metadata: Each unique sample_id will get a unique sequential number ordered by the EDITED date and time
  metadata <- edited_metadata %>%
    mutate(sample_num = as.numeric(factor(sample_id, levels = unique(sample_id))))

  # Save the *final* edited metadata
  write_csv2(metadata, file.path(output, "metadata",
                                        paste0("edited_metadata_",
                                               time,
                                               ".csv")))
  
  print("Edited metadata saved.")

  return(metadata)
 
}
