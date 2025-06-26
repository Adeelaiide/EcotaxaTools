#' NBSS.plot
#'
#' NBSS plot function.
#'
#' @param x bss data.table produced by "BSS_table"
#' @param taxon taxon you want to plot. Choose "total" if you want to see the global NBSS.
#' @param bv.type elli, plain or riddled biovolume. default is "elli"
#' @param samples samples you want to plot.
#'
#' @return A NBSS ggplot
#' @export
#'
#' @examples NBSS.plot(x=bss table of project, taxon="all", bv.type="elli", samples="all")
NBSS.plot <- function(x, taxon, bv.type, samples="all"){

  # Define min/max size
  x$max <- bv_to_esdum(x$max)
  mini <- min(x$max)
  maxi <- max(x$max)

  # Define colors according to common taxonomy vector
  N <- length(taxon)
 plankton_groups_colors <- c("#709699", #cyanobacteria
               "#F2F2F2", #detritus
               "#FAFABA", #other
               "#C9C4FF", #ciliophora
               "#B5D493", #dinoflagellata
               "#FAD4EB", #rhizaria
               "#1A9C75", #bacillariophyta
               "#E3E699", #dictyochophyceae
               "#EDB022", #crustacea
               "#E88F6B", #copepoda
               "#B0D6E0", #chaetognatha
               "#3B638A", #tunicata
               "#6999C7", #cnidaria
               "#C68181", #mollusca
               "#668F3B", #coccolithophyceae
               "#FFD3CF", #other_unidentified
               "#0073BD" #plastics
               )

  names(plankton_groups_colors)<- c("cyanobacteria","detritus","other","ciliophora","dinoflagellata",
                     "rhizaria","bacillariophyta","dictyochophyceae","crustacea",
                     "copepoda","chaetognatha","tunicata","cnidaria","mollusca",
                     "coccolithophyceae","other_unidentified","plastics")
  plankton_groups_colScale <- scale_colour_manual(name = "taxonomic group",values = plankton_groups_colors)
  plankton_groups_colFill <- scale_fill_manual(name = "taxonomic group",values = plankton_groups_colors)

  if(samples=="all") samples <- unique(x$sample_num)

  if(taxon=="total") {
    g <- x %>% filter(type==bv.type, sample_id %in% samples) %>% group_by(sample_id, max, norm) %>%
      summarise(BV=sum(BV, na.rm=T)) %>%
      ggplot(aes(x=max, y=BV/norm)) +
      geom_line() +
      geom_point() +
      scale_x_log10("Size (mm)", limits=c(mini,maxi)) +
      scale_y_log10("NBSS (mm\u00b3.mm\u207B\u00b3.m\u207B\u00b3)", labels=trans_format('log10',math_format(10^.x))) +
      theme_minimal()
  } else {
    g <- x %>% filter(type==bv.type, object_annotation_category %in% taxon, sample_id %in% samples) %>%
      ggplot(aes(x=max, y=BV/norm, color=taxon)) +
      geom_line() +
      plankton_groups_colFill +
      geom_point() +
      scale_x_log10("Size (mm)", limits=c(mini,maxi)) +
      scale_y_log10("NBSS (mm\u00b3.mm\u207B\u00b3.m\u207B\u00b3)", labels=trans_format('log10',math_format(10^.x))) +
      theme_minimal()
  }

  if (length(samples)>1) g <- g + facet_wrap(~sample_id)

  print(g)
}
