


#' add.taxo
#'
#' Return a taxonomic table with different levels per column based on annotation hierarchy from EcoTaxa.
#' The taxa are separated by the delim ">".
#'
#' @param object_hierarchy_vector Vector of object annotation hierarchy from EcoTaxa.
#'
#' @return Return a taxonomic table.
#' @export
#'
#' @examples add.taxo(object_hierarchy_vector)
add.taxo <- function(object_hierarchy_vector) {
  # 1. Create initial data frame
  taxo <- data.frame(object_annotation_hierarchy = object_hierarchy_vector)
  
  # 2. Dynamically determine max taxonomic depth to name columns cleanly
  max_depth <- max(stringr::str_count(object_hierarchy_vector, ">")) + 1
  col_names <- paste0("n", 1:max_depth)
  
  # 3. Process using native, vectorized tidyverse operations
  taxo <- taxo %>%
    separate_wider_delim(
      cols = object_annotation_hierarchy,
      delim = ">",
      names = col_names,
      too_few = "align_start",
      cols_remove = FALSE
    ) %>%
    tidyr::fill(all_of(col_names), .direction = "down") %>%
    mutate(
      object_annotation_hierarchy2 = stringr::str_replace_all(object_annotation_hierarchy, "<", ">")
    )
  
  return(taxo)
}
