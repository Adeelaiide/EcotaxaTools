# --- HELPER FUNCTIONS FOR compute_bv ---

# Helper function for Planktoscope-specific processing
process_planktoscope_data <- function(data, metadata) {
  # Planktoscope-specific initial mutate cols and unique_id creation
  data <- data %>%
    group_by(object_label, acq_id) %>%
    mutate(ghost_id = 1:n()) %>%
    ungroup() %>%
    mutate(unique_id = paste(acq_id, sample_operator, ghost_id,
                             object_date, object_time,
                             object_lat, object_lon, sep = "_")) %>%
    # Ensure all columns exist, creating NA columns if they are missing
    mutate(sample_total_volume = ifelse("sample_total_volume" %in% colnames(.), sample_total_volume, NA),
           sample_concentrated_sample_volume = ifelse("sample_concentrated_sample_volume" %in% colnames(.), sample_concentrated_sample_volume, NA),
           sample_dilution_factor = ifelse("sample_dilution_factor" %in% colnames(.), sample_dilution_factor, NA)) %>%
    mutate(sample_dilution_factor = as.numeric(gsub(",", ".", sample_dilution_factor)))
  
  # Metadata update (if metadata is provided)
  if (!is.null(metadata)) {
    if ("object_time" %in% colnames(metadata) && !inherits(metadata$object_time, "hms")) {
      metadata$object_time <- hms::as_hms(metadata$object_time)
    }
    if ("object_date" %in% colnames(metadata) && !inherits(metadata$object_date, "Date")) {
      metadata$object_date <- lubridate::as_date(metadata$object_date)
    }
    for (i in unique(data$unique_id)) {
      meta_row <- metadata[metadata$unique_id == i, ]
      if (nrow(meta_row) > 0) {
        update_cols <- intersect(names(meta_row), names(data))
        data[data$unique_id == i, update_cols] <- meta_row[1, update_cols]
      }
    }
  }
  
  # Planktoscope-specific unit conversions
  data <- mutate(data,
                 sample_total_volume = sample_total_volume / 1000,
                 sample_concentrated_sample_volume = sample_concentrated_sample_volume / 1000000,
                 acq_celltype = acq_celltype / 1000000,
                 acq_imaged_volume = acq_imaged_volume / 1000000,
                 pixelsize = process_pixel / 1000)
  
  # Planktoscope-specific metadata check columns
 metadata_cols <- c(
    "unique_id", "object_date", "object_time", "object_lat", "object_lon",
    "acq_id", "sample_id", "sample_operator", "sample_total_volume",
    "sample_concentrated_sample_volume", "acq_celltype", "acq_imaged_volume",
    "pixelsize", "sample_dilution_factor"
  )
  metadata_check <- select(data, all_of(metadata_cols)) %>% distinct()
  
  # Replace NA by 1 
  for (i in unique(metadata_check$unique_id)) {
    if (sum(is.na(metadata_check[metadata_check$unique_id == i, ])) > 0) {
      pb <- colnames(metadata_check[metadata_check$unique_id == i, ])[is.na(metadata_check[metadata_check$unique_id == i, ])]
      print(paste0("Warning! [sample: ", unique(data$sample_id), "] These metadata do not have value. Default is 1: ", paste(pb, collapse = ", ")))
      for (col_name_to_set_one in pb) {
        if (col_name_to_set_one %in% colnames(data)) {
          data[[col_name_to_set_one]][is.na(data[[col_name_to_set_one]]) & data$unique_id == i] <- 1
        }
      }
      metadata_check[is.na(metadata_check) & metadata_check$unique_id == i] <- 1
    }
  }
  
  
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
  # FlowCam-specific initial mutate cols and unique_id creation
  data <- data %>%
    group_by(object_id, acq_id) %>% 
    mutate(ghost_id = 1:n()) %>%
    ungroup() %>%
    mutate(unique_id = paste(acq_id, ghost_id,
                             object_date, object_time,
                             object_lat, object_lon, sep = "_")) %>%
    # Ensure all columns exist, creating NA columns if they are missing
    mutate(sample_initial_col_vol_m3 = ifelse("sample_initial_col_vol_m3" %in% colnames(.), sample_initial_col_vol_m3, NA),
           sample_conc_vol_ml = ifelse("sample_conc_vol_ml" %in% colnames(.), sample_conc_vol_ml, NA),
           sample_volconc = ifelse("sample_volconc" %in% colnames(.), sample_volconc, NA)) %>%
    mutate(sample_volconc = as.numeric(gsub(",", ".", sample_volconc))) 
  
 # Metadata update (if metadata is provided)
  if (!is.null(metadata)) {
    if ("object_time" %in% colnames(metadata) && !inherits(metadata$object_time, "hms")) {
      metadata$object_time <- hms::as_hms(metadata$object_time)
    }
    if ("object_date" %in% colnames(metadata) && !inherits(metadata$object_date, "Date")) {
      metadata$object_date <- lubridate::as_date(metadata$object_date)
    }
    for (i in unique(data$unique_id)) {
      meta_row <- metadata[metadata$unique_id == i, ]
      if (nrow(meta_row) > 0) {
        update_cols <- intersect(names(meta_row), names(data))
        data[data$unique_id == i, update_cols] <- meta_row[1, update_cols]
      }
    }
  }
  
  # FlowCam-specific unit conversions
  data <- mutate(data,
                 sample_conc_vol_ml = sample_conc_vol_ml / 1000000,
                 acq_celltype = acq_celltype / 1000000,
                 acq_fluid_volume_imaged = acq_fluid_volume_imaged / 1000000,
                 pixelsize = process_pixel / 1000) 
  
  # FlowCam-specific metadata check columns
  metadata_cols <- c(
    "unique_id", "object_date", "object_time", "object_lat", "object_lon",
    "acq_id", "sample_id", "sample_initial_col_vol_m3", "sample_conc_vol_ml",
    "acq_celltype", "acq_fluid_volume_imaged", "pixelsize", "sample_volconc"
  )
  metadata_check <- select(data, all_of(metadata_cols)) %>% distinct()
  
  # Replace NA by 1 
  for (i in unique(metadata_check$unique_id)) {
    if (sum(is.na(metadata_check[metadata_check$unique_id == i, ])) > 0) {
      pb <- colnames(metadata_check[metadata_check$unique_id == i, ])[is.na(metadata_check[metadata_check$unique_id == i, ])]
      print(paste0("Warning! [sample: ", unique(data$sample_id), "] These metadata do not have value. Default is 1: ", paste(pb, collapse = ", ")))
      for (col_name_to_set_one in pb) {
        if (col_name_to_set_one %in% colnames(data)) {
          data[[col_name_to_set_one]][is.na(data[[col_name_to_set_one]]) & data$unique_id == i] <- 1
        }
      }
      metadata_check[is.na(metadata_check) & metadata_check$unique_id == i] <- 1
    }
  }
  
  # FlowCam-specific conver.uniqueID formula
  data <- mutate(data, conver.uniqueID = (sample_conc_vol_ml) /
                   (acq_fluid_volume_imaged * sample_initial_col_vol_m3 * sample_volconc))
  
  # FlowCam-specific conver.sample calculation
   vimgsample <- data %>% select(sample_id, unique_id, acq_fluid_volume_imaged) %>% distinct() %>%
    group_by(sample_id) %>% summarize(sample_imaged_volume = sum(acq_fluid_volume_imaged, na.rm = TRUE))
  data <- merge(data, vimgsample, "sample_id", all.x = TRUE)
  data <- mutate(data, conver.sample =  (sample_conc_vol_ml) /
                   (acq_fluid_volume_imaged * sample_initial_col_vol_m3 * sample_volconc))
  
  return(data)
}

# Helper function for ZooScan-specific processing
process_zooscan_data <- function(data, metadata) {
    # ZooScan-specific initial mutate cols and unique_id creation
  data <- data %>%
    group_by(object_id, acq_id) %>% 
  mutate(ghost_id = 1:n()) %>%
    ungroup() %>%
    mutate(unique_id = paste(acq_id, sample_scan_operator, ghost_id,
                             object_date, object_time,
                             object_lat, object_lon, acq_sub_part, sep = "_")) %>%
  # Ensure all columns exist, creating NA columns if they are missing
   mutate(sample_tot_vol = ifelse("sample_tot_vol" %in% colnames(.), sample_tot_vol, NA),
          acq_sub_part = ifelse("acq_sub_part" %in% colnames(.), acq_sub_part, NA),
          process_particle_pixel_size_mm = ifelse("process_particle_pixel_size_mm" %in% colnames(.), process_particle_pixel_size_mm, NA)) %>%
  
  # Metadata update (if metadata is provided)
  if (!is.null(metadata)) {
    if ("object_time" %in% colnames(metadata) && !inherits(metadata$object_time, "hms")) {
      metadata$object_time <- hms::as_hms(metadata$object_time)
    }
    if ("object_date" %in% colnames(metadata) && !inherits(metadata$object_date, "Date")) {
      metadata$object_date <- lubridate::as_date(metadata$object_date)
    }
    for (i in unique(data$unique_id)) {
      meta_row <- metadata[metadata$unique_id == i, ]
      if (nrow(meta_row) > 0) {
        update_cols <- intersect(names(meta_row), names(data))
        data[data$unique_id == i, update_cols] <- meta_row[1, update_cols]
      }
    }
  }
  
  # ZooScan-specific unit conversions
  data <- mutate(metadata,
                pixelsize = unique(process_particle_pixel_size_mm) #,
                #perimferet = object_feret * pixelsize
                )
  
  # ZooScan-specific metadata check columns
 metadata_cols <- c(
    "unique_id", "object_date", "object_time", "object_lat", "object_lon",
    "acq_id", "sample_id", "sample_scan_operator", "acq_sub_part", "pixelsize", "sample_tot_vol")
  metadata_check <- select(data, all_of(metadata_cols)) %>% distinct()
  
  # Replace NA by 1 
  for (i in unique(metadata_check$unique_id)) {
    if (sum(is.na(metadata_check[metadata_check$unique_id == i, ])) > 0) {
      pb <- colnames(metadata_check[metadata_check$unique_id == i, ])[is.na(metadata_check[metadata_check$unique_id == i, ])]
      print(paste0("Warning! [sample: ", unique(data$sample_id), "] These metadata do not have value. Default is 1: ", paste(pb, collapse = ", ")))
      for (col_name_to_set_one in pb) {
        if (col_name_to_set_one %in% colnames(data)) {
          data[[col_name_to_set_one]][is.na(data[[col_name_to_set_one]]) & data$unique_id == i] <- 1
        }
      }
      metadata_check[is.na(metadata_check) & metadata_check$unique_id == i] <- 1
    }
  }
  
  
# ZooScan-specific conver.uniqueID formula
  data <- mutate(data, conver.uniqueID = (acq_sub_part) / (sample_tot_vol))
  
# ZooScan-specific conver.sample calculation
  vimgsample <- data %>% select(sample_id, unique_id, acq_sub_part) %>% distinct() %>%
    group_by(sample_id) %>% summarize(sample_imaged_volume = sum(acq_sub_part, na.rm = TRUE))
  data <- merge(data, vimgsample, "sample_id", all.x = TRUE)
  data <- mutate(data, conver.sample = (acq_sub_part) / (sample_tot_vol))
  
  return(data)
}

# Helper function for IFCB-specific processing
#process_ifcb_data <- function(data, metadata) {
 # IFCB-specific unit conversions and variable's names conversions
 # data <- mutate(metadata,
                 # vol = unique(acq_volume_sampled) / 1000000, 
                #  pixelsize = (1 / unique(acq_resolution_pixel_per_micron)) / 1000,
               #   summedbiovolume = object_summed_biovolume * (pixelsize^3),
                #  summedarea = object_summed_surface_area * (pixelsize^2),
                #  object_major = object_major_axis_length,
                #  object_minor = object_minor_axis_length,
                #  object_area = object_surface_area) 
  
 # IFCB-specific conver.uniqueID formula
  #data <- mutate(data, conver.uniqueID = 1 / vol)
  
  # IFCB-specific conver.sample calculation
 # vimgsample <- data %>% select(sample_id, unique_id, vol) %>% distinct() %>%
   # group_by(sample_id) %>% summarize(sample_imaged_volume = sum(vol, na.rm = TRUE))
  #data <- merge(data, vimgsample, "sample_id", all.x = TRUE)
 # data <- mutate(data, conver.sample =  1 / vol)
  
  #return(data)
#}



