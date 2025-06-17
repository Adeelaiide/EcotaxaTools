

#' add.trophiclvl
#'
#' Add trophic levels and OTU from a database constructed within the package, or external.
#' object_annotation_hierarchy2 = object annotation hierarchy from ecotaxa with corrected delim >.
#' Type = OTU main type "Living, non-living, etc."
#' Sub_Type = OTU sub type "Diatoms, Copepods, etc."
#' Value = Trophic level "Unknown = 3.5, Predators = 3, Omnivorous = 2.5, Grazer = 2, Mixotrophe = 1.5, Phototroph = 1 and None = -1"
#'
#' @param taxo OTU table generated with "add_taxo" containing trophic levels and main OTU groups.
#' @param output file output to save the new zoo table if edited.
#'
#' @return A taxonomic table with also trophic levels.
#' @export
#'
#' @examples add.trophiclvl(taxo=table from add.taxo, output=where to save the table)
add.trophiclvl <- function(taxo, output){
  # load the original database of OTU and make it "clean"
  zo <- trophic_affiliation_of_organisms
  zoo <- merge(taxo, zo, all.x=T)
  zoo$Type[zoo$n1=="temporary"] <- "temporary"
  zoo$Sub_type[zoo$n1=="temporary"] <- "temporary"
  zoo$Type[zoo$n1=="not-living"] <- "non_living"
  zoo$Sub_type[zoo$n1=="not-living"] <- zoo$n2[zoo$n1=="not-living"]
  zoo$Value[zoo$Type=="non_living"] <- -1

  # find  objects who are not in the otu database
  liste.choix <- zo %>% mutate(Category=paste0(zo$Type,">",zo$Sub_type)) %>%
    select(-object_annotation_hierarchy2, -Value) %>% distinct()
  liste.value <- data.frame(Value=c(-1,1,1.5,2,2.5,3,3.5),
                            Trophic=c("None","Phototroph","Mixotroph","Grazer","Omnivorous","Predator","Unknown"))
  non <- zoo$object_annotation_hierarchy2[is.na(zoo$Type)]

  if(sum(is.na(zoo$Type)>0)) {
    yesno2 <- dlg_message(paste(c("These taxa do not exist in the database. Would you like to create them?\n",
                                  "If no, you can also import an existing database or ignore them.\n",non)),
                          type="yesno")$res
    if(yesno2=="yes") {
      # You can edit the new taxa in the database
      replace <- data.frame(object_annotation_hierarchy2=non, Category="temporary>temporary", Trophic=NA)
      replace <- data_edit(replace,
                           col_options = list(Category = c(liste.choix$Category),
                                              Trophic = c(liste.value$Trophic)), viewer="pane")

      replace <- replace %>% separate(Category, into=c("Type","Sub_type")) %>%
        merge(liste.value, "Trophic") %>% select(-Trophic)
      zo <- bind_rows(zo, replace) %>% as.data.frame()
      write_csv2(zo, file.path(output,"summary", "trophic_affiliation_of_organisms.csv"))
      # Restart the process and ignore if NA
      zoo <- merge(taxo, zo, all.x=T)
      zoo$Type[zoo$n1=="temporary"] <- "temporary"
      zoo$Sub_type[zoo$n1=="temporary"] <- "temporary"
      zoo$Type[zoo$n1=="not_living"] <- "non_living"
      zoo$Sub_type[zoo$n1=="not_living"] <- zoo$n2[zoo$n1=="not_living"]
      zoo$Value[zoo$Type=="non_living"] <- -1
      zoo$Type[is.na(zoo$Type)] <- "temporary"
      zoo$Sub_type[is.na(zoo$Type)] <- "temporary"

    } else {
      yesno3 <- dlg_message("Do you want to import an existing database ? If no, they will be classified as temporary.",
                            type="yesno")$res
      if(yesno3=="yes"){
        # You can also import an already edited table
        zo <- file.choose() %>% read_csv2()
        # Restart the process and ignore if NA
        zoo <- merge(taxo, zo, all.x=T)
        zoo$Type[zoo$n1=="temporary"] <- "temporary"
        zoo$Sub_type[zoo$n1=="temporary"] <- "temporary"
        zoo$Type[zoo$n1=="not-living"] <- "non_living"
        zoo$Sub_type[zoo$n1=="not-living"] <- zoo$n2[zoo$n1=="not-living"]
        zoo$Value[zoo$Type=="non_living"] <- -1
        zoo$Type[is.na(zoo$Type)] <- "temporary"
        zoo$Sub_type[is.na(zoo$Type)] <- "temporary"
      } else{
        # You can ignore them
        zoo$Type[is.na(zoo$Type)] <- "temporary"
        zoo$Sub_type[is.na(zoo$Type)] <- "temporary"
      }
    }
  }
  return(zoo)
}
