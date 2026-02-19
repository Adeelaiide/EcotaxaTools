#' select_taxa_to_analyse
#'
#' Allow user to keep or remove taxonomical groups of interest for analysis.
#'
#' 
#' @param data dataframe generated during multivariates_analysis.R. It is Species X Samples
#' with the first four columnes with sample_ID,sample_num, latitude and longitude
#' @param path file output to save the selection.
#'
#' @return A dataframe filtered with the selected taxa
#' @export
#'
#' @examples select_taxa_to_analyse(data= analyse_dataset, path= output)

select_taxa_to_analyse <- function(data,path){

yesno <- dlg_message(paste(c("Do you want to delete some groups or taxa (non-living for exemple) before the analysis ?\n",
                              "If no, you can also import an existing database or ignore them.\n")),
                      type="yesno")$res

if(yesno=="yes") {
  
  replace <- data.frame(Groups = colnames(data[-c(1:4)]), Status = "Keep")
  print("Fill the value 'Keep' or 'Remove' in the Data Editor panel. Click on Synchronise and Done to save the changes")
  replace <- data_edit(replace,
                       col_options = list(Status = c("Keep","Remove")), viewer="pane")
  
  ## Check if all groups have been given a value
  replace$Status[replace$Status == ""] <- NA
  
  if(any(is.na(replace$Status))){
    dlg_message("Status is missing in some groups. They are keep by precaution",type = "ok")
    replace$Status[is.na(replace$Status)] <- "Keep"}
  
  ## Save the species selected 
  keep_species <- replace %>%
    filter(Status == "Keep") %>%
    pull(Groups)
  
  ## Stop the script if all groups are removed
  if(length(keep_species) == 0){
    dlg_message("All groups were removed. No analysis can be performed", type = "ok")
    stop()}
  
  ## Filtering the database with our selection 
  data_filtered <- data %>%
    select(1:4, all_of(keep_species))
  
}else{
  
  yesno2 <- dlg_message("Do you want to import an existing taxa selection? If no, everything will be kept for the analysis.",
                        type="yesno")$res
  
  if(yesno2=="yes") {
    
    replace <- file.choose() %>% read_csv2()
    
    ## Check if all groups have been given a value
    replace$Status[replace$Status == ""] <- NA
    
    if(any(is.na(replace$Status))){
      dlg_message("Status is missing in some groups. They are keep by precaution",type = "ok")
      replace$Status[is.na(replace$Status)] <- "Keep"}
    
    ## Save the species selected 
    keep_species <- replace %>%
      filter(Status == "Keep") %>%
      pull(Groups)
    
    ## Stop the script if all groups are removed
    if(length(keep_species) == 0){
      dlg_message("All groups were removed. No analysis can be performed", type = "ok")
      stop()}
    
    ## Filtering the database with our selection 
    data_filtered <- data %>%
      select(1:4, all_of(keep_species))
  }else{
    data_filtered<-data
    }
}

write_csv(replace, file.path(path,"Taxa selected for multivariate analysis.csv"))
return(data_filtered)
}