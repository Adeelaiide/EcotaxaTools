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

# Helper function for ZooScan-specific processing
process_zooscan_data <- function(data, metadata) {
 # ZooScan-specific unit conversions
  data <- mutate(metadata,
                pixelsize = unique(process_particle_pixel_size_mm),
                perimferet = object_feret * pixelsize) 
  
 # ZooScan-specific conver.uniqueID formula
  data <- mutate(data, conver.uniqueID = unique(acq_sub_part) / unique(sample_tot_vol))
  
  # ZooScan-specific conver.sample calculation
  vimgsample <- data %>% select(sample_id, unique_id, acq_sub_part) %>% distinct() %>%
    group_by(sample_id) %>% summarize(sample_imaged_volume = sum(acq_sub_part, na.rm = TRUE))
  data <- merge(data, vimgsample, "sample_id", all.x = TRUE)
  data <- mutate(data, conver.sample = unique(acq_sub_part) / unique(sample_tot_vol))
  
  return(data)
}

# Helper function for IFCB-specific processing
process_ifcb_data <- function(data, metadata) {
 # IFCB-specific unit conversions and variable's names conversions
  data <- mutate(metadata,
                 vol = unique(acq_volume_sampled) / 1000000, 
                 pixelsize = (1 / unique(acq_resolution_pixel_per_micron)) / 1000,
                 summedbiovolume = object_summed_biovolume * (pixelsize^3),
                 summedarea = object_summed_surface_area * (pixelsize^2),
                 object_major = object_major_axis_length,
                 object_minor = object_minor_axis_length,
                 object_area = object_surface_area) 
  
 # IFCB-specific conver.uniqueID formula
  data <- mutate(data, conver.uniqueID = 1 / vol)
  
  # IFCB-specific conver.sample calculation
  vimgsample <- data %>% select(sample_id, unique_id, vol) %>% distinct() %>%
    group_by(sample_id) %>% summarize(sample_imaged_volume = sum(vol, na.rm = TRUE))
  data <- merge(data, vimgsample, "sample_id", all.x = TRUE)
  data <- mutate(data, conver.sample =  1 / vol)
  
  return(data)
}


