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
check_metadata <- function(path, output, instru) {
  # Columns to hide for the data edit but needed for compute_bv and Instrument_Specific_processing
  cols_to_hide <- c("object_area_exc", "object_feret", "object_area", "object_major", "object_minor")

  if (!file.exists(file.path(output,"metadata"))) {
    dir.create(file.path(output,"metadata"))
  }

  # the helper function for reading base data is used here
  if (instru == "PlanktoScope") {
    print("You chose PlanktoScope. Applying PlanktoScope specific processing...")    
    metadata <- do.call("rbind", lapply(path, function(p) {
      read_base_metadata_file(p) %>% transform_planktoscope_data()
    }))
  } else if (instru == "FlowCam") {
    print("You chose FlowCam. Applying FlowCam specific processing...")
    metadata <- do.call("rbind", lapply(path, function(p) {
      read_base_metadata_file(p) %>% transform_flowcam_data()
    }))
  } else if (instru == "ZooScan") {
    print("You chose ZooScan. Applying ZooScan specific processing...")
    metadata <- do.call("rbind", lapply(path, function(p) {
      read_base_metadata_file(p) %>% transform_zooscan_data()
    }))
  } else if (instru == "IFCB") {
    print("You chose IFCB. Applying IFCB specific processing...")
    metadata <- do.call("rbind", lapply(path, function(p) {
      read_base_metadata_file(p) %>% transform_ifcb_data()
    }))
  } else {
    stop("Error: Invalid instrument specified. Choose 'PlanktoScope', 'FlowCam', 'ZooScan', or 'IFCB'.") 
  }

  # --- Rest of the function (common steps for all instruments) ---
  
  # Save original
  write_csv2(metadata, file.path(output,"metadata","original_metadata.csv"))
  print("Original metadata saved.")

   # Extract and save the hidden columns, along with the key identifier
  edited_metadata_for_viewer <- metadata %>%
    select(-all_of(cols_to_hide)) %>%
    distinct()

  # DATA EDIT
  check <- metadata %>% group_by(sample_id) %>% summarize(nb=n())
  check <- max(check$nb, na.rm=T)
  if (check>1){
    dlg_message("Warning ! Some of your samples have more than one acquisition. Be sure that they belong to the same sample. If not, please create a different tsv file for the supplementary acquisition and restart the process.", type="ok")
  }

  time <- format(Sys.time(), "%d-%m-%Y_%H%M")
  dlg_message("You can now edit metadata (click on SYNCHRONIZE and DONE button to update edition). The original .tsv files will not be edited. NA will be replaced by 1. Do not change the unique_id.", type="ok")

  edited_metadata_for_viewer$object_date <- as.character(edited_metadata_for_viewer$object_date)
  edited_metadata_for_viewer$object_time <- as.character(edited_metadata_for_viewer$object_time)
  edited_metadata <- data_edit(edited_metadata_for_viewer, viewer="pane")

  print("Data editing completed.")

  # Arrange by the *edited* date and time before creating sample_num
  edited_metadata <- arrange(edited_metadata, object_date, object_time) %>%
    left_join(hidden_data, by = "unique_id")

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



