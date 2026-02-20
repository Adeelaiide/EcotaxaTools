#' ecotaxa_tools
#'
#' This function merges all tsv exported from Ecotaxa to a single database combining metadata,
#'  raw data, and quantitative descriptors (abundance, biovolume, etc.).
#' It returns a processed dataset with Abundance, Biovolume, size spectra and taxonomic units 
#' computes for each objects. 
#'
#' @param path list of path for exported tsv files
#' @param output path to save dataset and figures 
#' @param instru instrument chosen for the analysis
#'
#' @return A processed dataset called final_dataset
#' @export
#'
#' @examples ecotaxa_tools(path= path,output= output = "PlanktoScope")

ecotaxa_tools <- function(path,output,instru){
  
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
  
  # Create a colum with the trophic categories
  taxo <- taxo %>% mutate(Trophic_lvl = case_when( Value == 1 ~ "Phototrophs", Value == 1.5 ~ "Mixotrophs",
                                                   Value == 2 ~ "Grazers", Value == 2.5 ~ "Omnivorous",
                                                   Value == 3 ~ "Predators", Value == 3.5 ~ "Unknown",
                                                   Value == -1 ~ "None"))
  
  # Replacing all the NA in case the original metadata was selected
  final_metadata <- metadata
  final_metadata[is.na(final_metadata)] <- 1
  
  # Create final processed table for the user - Compute the ESD + Delete unnecessary metadata
  final_dataset <- merge(final_metadata, bss, all.x=T) %>% 
    mutate(ESD=bv_to_esdum(max)) %>% 
    merge(taxo[,c("object_annotation_hierarchy","Type","Sub_type","Value","Trophic_lvl")], all.x = T) 
  
  
  # Saving tables
  write_csv2(final_dataset, file.path(path.summary, "Processed database.csv"))
  write_csv2(final_metadata, file.path(path.summary, "metadata_used.csv"))
  write_csv2(bss, file.path(path.summary, "BSS.csv"))
  write_csv2(res, file.path(path.summary, "summary_all.csv"))
  write_csv2(taxo, file.path(path.summary, "taxo.csv"))
  
  return(final_dataset)
}