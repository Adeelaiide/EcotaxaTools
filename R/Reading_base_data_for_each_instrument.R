# --- HELPER FUNCTIONS FOR check_metadata ---

#Helper Function for Reading Base Data
#This function handles the initial reading and common column type conversions
read_base_metadata_file <- function(file_path) {
  read_tsv(file_path, col_types = list(object_time = col_time(),
                                       object_date = col_date(),
                                       object_annotation_time = col_time()))
}

#Instrument-Specific Transformation Functions

# Function for PlanktoScope specific data transformation
transform_planktoscope_data <- function(df) {
  df %>%
    group_by(object_id, acq_id) %>% 
    mutate(unique_id = paste(acq_id, sample_operator, object_date, object_time,
                             object_lat, object_lon, sep = "_")) %>%
    group_by(unique_id) %>%
    mutate(number_object = n(),
           percentValidated = sum(object_annotation_status == "validated", na.rm = TRUE) / n() * 100) %>%
    ungroup() %>%
  # Ensure all columns exist, creating NA columns if they are missing
    mutate(sample_total_volume = ifelse("sample_total_volume" %in% colnames(.), sample_total_volume, NA),
           sample_concentrated_sample_volume = ifelse("sample_concentrated_sample_volume" %in% colnames(.), sample_concentrated_sample_volume, NA),
           sample_dilution_factor = ifelse("sample_dilution_factor" %in% colnames(.), sample_dilution_factor, NA),
           acq_imaged_volume = ifelse("acq_imaged_volume" %in% colnames(.), acq_imaged_volume, NA),
           acq_celltype = ifelse("acq_celltype" %in% colnames(.), acq_celltype, NA),
           process_pixel = ifelse("process_pixel" %in% colnames(.), process_pixel, NA)) %>%
  
    mutate(sample_dilution_factor = as.numeric(gsub(",", ".", sample_dilution_factor))) %>%
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
           acq_minimum_mesh,
           acq_maximum_mesh,
           sample_total_volume,
           sample_concentrated_sample_volume,
           acq_celltype,
           acq_imaged_volume,
           process_pixel,
           sample_dilution_factor,
          object_area,
         object_major,
         object_minor) %>%
    distinct() %>%
    group_by(sample_id) %>% mutate(ghost_id=1:n()) %>% ungroup() 
}

# Function for FlowCam specific data transformation
transform_flowcam_data <- function(df) {
  df %>%
    group_by(object_id, acq_id) %>% 
    mutate(unique_id = paste(acq_id, object_date, object_time,
                             object_lat, object_lon, sep = "_")) %>%
    group_by(unique_id) %>%
    mutate(number_object = n(),
           percentValidated = sum(object_annotation_status == "validated", na.rm = TRUE) / n() * 100) %>%
    ungroup() %>%
  # Ensure all columns exist, creating NA columns if they are missing
    mutate(sample_initial_col_vol_m3 = ifelse("sample_initial_col_vol_m3" %in% colnames(.), sample_initial_col_vol_m3, NA),
           sample_conc_vol_ml = ifelse("sample_conc_vol_ml" %in% colnames(.), sample_conc_vol_ml, NA),
           sample_volconc_temp = ifelse("sample_volconc" %in% colnames(.), sample_volconc, NA_real_),
           acq_fluid_volume_imaged = ifelse("acq_fluid_volume_imaged" %in% colnames(.), acq_fluid_volume_imaged, NA),
           acq_celltype = ifelse("acq_celltype" %in% colnames(.), acq_celltype, NA),
           process_pixel = ifelse("process_pixel" %in% colnames(.), process_pixel, NA),
           sample_volconc = as.numeric(gsub(",", ".", sample_volconc_temp)), 
           acq_celltype = parse_number(acq_celltype)) %>%
    select(-sample_volconc_temp) %>%
    select(sample_id,
           acq_id,
           unique_id,
           ghost_id,
           object_date,
           object_time,
           object_lat,
           object_lon,
           percentValidated,
           number_object,
           acq_raw_image_total,
           acq_celltype,
           acq_min_esd,
           acq_max_esd,
           sample_initial_col_vol_m3,
           sample_conc_vol_ml,
           acq_fluid_volume_imaged,
           process_pixel,
           sample_volconc,
          object_area,
         object_major,
         object_minor) %>%
    distinct() %>%
    group_by(sample_id) %>% mutate(ghost_id=1:n()) %>% ungroup() 
}

# Function for Zooscan specific data transformation
transform_zooscan_data <- function(df) {
  df %>%
    group_by(object_id, acq_id) %>% 
    mutate(unique_id = paste(acq_id, object_id, sample_scan_operator, 
                             object_date, object_time,
                             object_lat, object_lon, sep = "_")) %>%
   group_by(unique_id) %>%
    mutate(number_object = n(),
            percentValidated = sum(object_annotation_status == "validated", na.rm = TRUE) / n() * 100) %>%
  ungroup() %>%
  # Ensure all columns exist, creating NA columns if they are missing
   mutate(sample_tot_vol = ifelse("sample_tot_vol" %in% colnames(.), sample_tot_vol, NA),
           acq_sub_part = ifelse("acq_sub_part" %in% colnames(.), acq_sub_part, NA),
          object_feret = ifelse("object_feret" %in% colnames(.), object_feret, NA),
           object_area = ifelse("object_area" %in% colnames(.), object_area, NA),
           object_major = ifelse("object_major" %in% colnames(.), object_major, NA),
           object_minor = ifelse("object_minor" %in% colnames(.), object_minor, NA)) %>%
  select(sample_id,
         unique_id,
         ghost_id,
         number_object,
         sample_scan_operator,
         sample_barcode,
         sample_tot_vol,
         object_date,
         object_time,
         object_lat,
         object_lon,
         acq_id,
         acq_min_mesh,
         acq_max_mesh,
         acq_sub_part,
         process_particle_pixel_size_mm,
         object_feret,
         object_area,
         object_major,
         object_minor,
         percentValidated) %>%
   distinct() %>%
    group_by(sample_id) %>% mutate(ghost_id=1:n()) %>% ungroup() 
}

# Function for IFCB specific data transformation
transform_ifcb_data <- function(df) {
  df %>% 
  group_by(object_id, acq_id) %>% 
    mutate(unique_id = paste(acq_id, object_date, object_time,
                             object_lat, object_lon, sep = "_")) %>%
  group_by(unique_id) %>%
    mutate(number_object = n(), 
           percentValidated = sum(object_annotation_status == "validated") / n() * 100) %>%
  ungroup() %>%
  # Ensure all columns exist, creating NA columns if they are missing
  mutate(object_lat_end = ifelse("object_lat_end" %in% colnames(.), object_lat_end, NA_real_),
         object_lon_end = ifelse("object_lon_end" %in% colnames(.), object_lon_end, NA_real_),
         acq_volume_sampled = ifelse("acq_volume_sampled" %in% colnames(.), acq_volume_sampled, NA),
         acq_resolution_pixel_per_micron = ifelse("acq_resolution_pixel_per_micron" %in% colnames(.), acq_resolution_pixel_per_micron, NA),
         object_major_axis_length = ifelse("object_major_axis_length" %in% colnames(.), object_major_axis_length, NA),
         object_minor_axis_length = ifelse("object_minor_axis_length" %in% colnames(.), object_minor_axis_length, NA),
         object_surface_area = ifelse("object_surface_area" %in% colnames(.), object_surface_area, NA),
         object_summed_biovolume = ifelse("object_summed_biovolume" %in% colnames(.), object_summed_biovolume, NA),
         object_summed_surface_area = ifelse("object_summed_surface_area" %in% colnames(.), object_summed_surface_area, NA)) %>%
  select(sample_id,
      acq_id,
      object_date,
      object_time,
      object_lat,
      object_lon,
      object_lat_end,
      object_lon_end,
      acq_volume_sampled,
      acq_resolution_pixel_per_micron,
      object_major_axis_length,
      object_minor_axis_length,
      object_surface_area,
      object_summed_biovolume,
      object_summed_surface_area,
      percentValidated) %>%
   distinct() %>%
    group_by(sample_id) %>% mutate(ghost_id=1:n()) %>% ungroup() 
}











