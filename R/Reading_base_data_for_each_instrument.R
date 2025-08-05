# --- HELPER FUNCTIONS FOR check_metadata ---

#Helper Function for Reading Base Data
#This function handles the initial reading and common column type conversions
read_base_metadata_file <- function(file_path) {
  read_tsv(file_path, col_types = list(object_time = col_time(),
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
    mutate(sample_total_volume = ifelse("sample_total_volume" %in% colnames(.), sample_total_volume, NA),
           sample_concentrated_sample_volume = ifelse("sample_concentrated_sample_volume" %in% colnames(.), sample_concentrated_sample_volume, NA),
           sample_dilution_factor = ifelse("sample_dilution_factor" %in% colnames(.), sample_dilution_factor, NA)) %>%
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
    mutate(sample_initial_col_vol_m3 = ifelse("sample_initial_col_vol_m3" %in% colnames(.), sample_initial_col_vol_m3, NA),
    sample_conc_vol_ml = ifelse("sample_conc_vol_ml" %in% colnames(.), sample_conc_vol_ml, NA),
    sample_volconc_temp = ifelse("sample_volconc" %in% colnames(.), sample_volconc, NA_real_), 
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
           acq_max_esd
           sample_initial_col_vol_m3,
           sample_conc_vol_ml,
           acq_fluid_volume_imaged,
           process_pixel,
           sample_volconc) %>%
    distinct() %>%
    group_by(sample_id) %>% mutate(ghost_id=1:n()) %>% ungroup() 
}

# Function for Zooscan specific data transformation
transform_zooscan_data <- function(df) {
  df %>%
    group_by(acq_id) %>% 
  mutate(pixelsize = unique(process_particle_pixel_size_mm),
      percentValidated = sum(object_annotation_status == "validated", na.rm = TRUE) / n() * 100,
      major = object_major * pixelsize,
      minor = object_minor * pixelsize,
      area_exc = object_area_exc * (pixelsize^2),
      area = object_area * (pixelsize^2),
      perimferet = object_feret * pixelsize,
      ESD = 2 * (((object_area * (pixelsize^2)) / pi)^0.5),
      conver = unique(acq_sub_part) / unique(sample_tot_vol)) %>%
  ungroup() %>%
  mutate(unique_id = paste(acq_id, object_id, sample_scan_operator, 
                             object_date, object_time,
                             object_lat, object_lon, sep = "_")) %>%
  select(sample_id,
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
      process_img_resolution,
      percentValidated,
      major,
      minor,
      area_exc,
      area,
      perimferet,
      ESD,
      conver) %>%
   distinct() %>%
    group_by(sample_id) %>% mutate(ghost_id=1:n()) %>% ungroup() 
}

# Function for IFCB specific data transformation
transform_ifcb_data <- function(df) {
  df %>%
  mutate(vol = unique(acq_volume_sampled) / 1000000, 
      pixelsize = (1 / unique(acq_resolution_pixel_per_micron)) / 1000, 
      percentValidated = 100 * sum(object_annotation_status == "validated") / n() ) %>%
  mutate(object_lat_end = ifelse("object_lat_end" %in% colnames(.), object_lat_end, NA_real_),
      object_lon_end = ifelse("object_lon_end" %in% colnames(.), object_lon_end, NA_real_)) %>%
   mutate(major = object_major_axis_length * pixelsize,
      minor = object_minor_axis_length * pixelsize,
      area = object_surface_area * (pixelsize^2),
      ESD = 2 * (((object_surface_area * (pixelsize^2)) / pi)^0.5),
      summedbiovolume = object_summed_biovolume * (pixelsize^3),
      summedarea = object_summed_surface_area * (pixelsize^2),
      conver = 1 / vol) %>%
  group_by(acq_id) %>% 
    mutate(unique_id = paste(acq_id, object_date, object_time,
                             object_lat, object_lon, sep = "_")) %>%
   select(sample_id,
      acq_id,
      object_date,
      object_time,
      object_lat,
      object_lon,
      object_lat_end,
      object_lon_end,
      vol,
      percentValidated,
      major,
      minor,
      area,
      ESD,
      summedbiovolume,
      summedarea,
      conver) %>%
   distinct() %>%
    group_by(sample_id) %>% mutate(ghost_id=1:n()) %>% ungroup() 
}



