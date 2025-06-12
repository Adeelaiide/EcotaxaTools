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
  myColors <- colorRampPalette(brewer.pal(8, "Set2"))(N)

  if(samples=="all") samples <- unique(x$sample_num)

  if(taxon=="total") {
    g <- x %>% filter(type==bv.type, sample_num %in% samples) %>% group_by(sample_num, max, norm) %>%
      summarise(BV=sum(BV, na.rm=T)) %>%
      ggplot(aes(x=max, y=BV/norm)) +
      geom_line() +
      geom_point() +
      scale_x_log10("Size (mm)", limits=c(mini,maxi)) +
      scale_y_log10("NBSS (mm\u00b3.mm\u207B\u00b3.m\u207B\u00b3)", labels=trans_format('log10',math_format(10^.x))) +
      theme_minimal()
  } else {
    g <- x %>% filter(type==bv.type, object_annotation_category %in% taxon, sample_num %in% samples) %>%
      ggplot(aes(x=max, y=BV/norm, color=taxon)) +
      geom_line() +
      scale_color_manual(values = myColors) +
      geom_point() +
      scale_x_log10("Size (mm)", limits=c(mini,maxi)) +
      scale_y_log10("NBSS (mm\u00b3.mm\u207B\u00b3.m\u207B\u00b3)", labels=trans_format('log10',math_format(10^.x))) +
      theme_minimal()
  }

  if (length(samples)>1) g <- g + facet_wrap(~sample_id)

  print(g)
}
