# --- HELPER FUNCTIONS FOR compute_bv ---

# Helper function for Planktoscope-specific processing
process_planktoscope_data <- function(data, metadata) {
 # Planktoscope-specific unit conversions
  data <- mutate(metadata,
                 sample_total_volume = sample_total_volume / 1000,
                 sample_concentrated_sample_volume = sample_concentrated_sample_volume / 1000000,
                 acq_imaged_volume = acq_imaged_volume / 1000000,
                 acq_celltype = acq_celltype/1000/1000,
                 pixelsize = process_pixel / 1000)

  # Planktoscope-specific conver.uniqueID formula
  data <- mutate(data, conver.uniqueID = (sample_concentrated_sample_volume) /
                   (acq_imaged_volume * sample_total_volume * sample_dilution_factor))
  
  # Planktoscope-specific conver.sample calculation
  vimgsample <- data %>% select(sample_id, unique_id, acq_imaged_volume) %>% distinct() %>%
    group_by(sample_id) %>% summarize(sample_imaged_volume = sum(acq_imaged_volume, na.rm = TRUE))
  data <- merge(data, vimgsample, "sample_id", all.x = TRUE)
  data <- mutate(data, conver.sample = (sample_concentrated_sample_volume) /
                   (sample_imaged_volume * sample_total_volume * sample_dilution_factor))
  
  return(data)
}

# Helper function for FlowCam-specific processing
process_flowcam_data <- function(data, metadata) {
 # FlowCam-specific unit conversions
  data <- mutate(metadata,
                 sample_conc_vol_ml = sample_conc_vol_ml / 1000000,
                 acq_fluid_volume_imaged = acq_fluid_volume_imaged / 1000000,
                 acq_celltype = acq_celltype/1000/1000,
                 pixelsize = process_pixel / 1000) 
  
 # FlowCam-specific conver.uniqueID formula
  data <- mutate(data, conver.uniqueID = (sample_conc_vol_ml) /
                   (acq_fluid_volume_imaged * sample_initial_col_vol_m3 * sample_volconc))
  
  # FlowCam-specific conver.sample calculation
  vimgsample <- data %>% select(sample_id, unique_id, acq_fluid_volume_imaged) %>% distinct() %>%
    group_by(sample_id) %>% summarize(sample_imaged_volume = sum(acq_fluid_volume_imaged, na.rm = TRUE))
  data <- merge(data, vimgsample, "sample_id", all.x = TRUE)
  data <- mutate(data, conver.sample = (sample_conc_vol_ml * sample_volconc) /
                   (sample_imaged_volume * sample_initial_col_vol_m3))
  
  return(data)
}

