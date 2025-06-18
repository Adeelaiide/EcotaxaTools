


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
add.taxo <- function(object_hierarchy_vector){

  # rep.last is a function that allow to copy the last taxonomic level available to the lower levels if they are not provided.
  rep.last <- function(x) {
    for (i in 2:length(x)) {if(is.na(x[i])) x[i] <- x[i-1]}
    return(x)}

  # separate by the delim ">" into columns
  taxo <- data.frame(n=object_hierarchy_vector) %>%
    separate_wider_delim(n, delim=">", names_sep = "", too_few = "align_start", cols_remove=F) %>%
    apply(1, rep.last) %>% t() %>% as.data.frame() %>%
    rename(object_annotation_hierarchy=nn) %>%
    mutate(object_annotation_hierarchy2 = gsub("<",">",object_annotation_hierarchy))

  return(taxo)
}
