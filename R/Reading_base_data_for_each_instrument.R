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
  group_by(sample_id, acq_id) %>% 
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
  mutate(sample_dilution_factor = as.numeric(gsub(",", ".", sample_dilution_factor)))  %>%
  select(sample_id,
           acq_id,
           unique_id, 
           object_date,
           object_time,
           object_lat,
           object_lon,
           sample_operator, 
           percentValidated,
           number_object,
           acq_minimum_mesh,
           acq_maximum_mesh,
           sample_total_volume,
           sample_concentrated_sample_volume,
           acq_celltype,
           acq_imaged_volume,
           process_pixel,
           sample_dilution_factor) %>%
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
           sample_volconc,
           acq_fluid_volume_imaged,
           process_pixel) %>%
    distinct() %>%
    group_by(sample_id) %>% mutate(ghost_id=1:n()) %>% ungroup() 
}

# Function for Zooscan specific data transformation
transform_zooscan_data <- function(df) {
  df %>%
    group_by(sample_id, acq_id) %>%
    mutate(unique_id = paste(acq_id, sample_scan_operator, 
                             object_date, object_time,
                             object_lat, object_lon, acq_sub_part, sep = "_")) %>%
    # Ensure all columns exist, creating NA columns if they are missing
    mutate(sample_tot_vol = ifelse("sample_tot_vol" %in% colnames(.), sample_tot_vol, NA),
           acq_sub_part = ifelse("acq_sub_part" %in% colnames(.), acq_sub_part, NA),
           process_particle_pixel_size_mm = ifelse("process_particle_pixel_size_mm" %in% colnames(.), process_particle_pixel_size_mm, NA)) %>%
   group_by(unique_id) %>%
    mutate(number_object = n(),
            percentValidated = sum(object_annotation_status == "validated", na.rm = TRUE) / n() * 100) %>%
  ungroup() %>%
    select(sample_id,
           acq_id,
           unique_id,
           object_date,
           object_time,
           object_lat,
           object_lon,
           sample_scan_operator, 
           percentValidated,
           number_object,
           acq_min_mesh,
           acq_max_mesh,
           acq_sub_part,
           sample_barcode,
           sample_tot_vol,
           process_particle_pixel_size_mm) %>%
           distinct() %>%
 group_by(sample_id) %>% mutate(ghost_id=1:n()) %>% ungroup() 
}


















