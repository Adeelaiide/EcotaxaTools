#' graph.metadata
#'
#' A set of graphics to resume the project metadata.
#' Return a map of the samples, date and time, and metadata comparison between all samples.
#'
#' @param metadata metadata table generated by "check_metadata"
#'
#' @return A set of graphics to resume the project metadata.
#' @export
#'
#' @examples graph.metadata(metadata.file from check_metadata)

# ordonner par datetime

graph.metadata <- function(metadata) {
  # 1. MAP
  ex = 3
  latmin <- min(metadata$object_lat, na.rm=T)-ex
  lonmin <- min(metadata$object_lon, na.rm=T)-ex
  latmax <- max(metadata$object_lat, na.rm=T)+ex
  lonmax <- max(metadata$object_lon, na.rm=T)+ex

  if(latmin<(-90)) latmin <- (-90)
  if(lonmin<(-180)) lonmin <- -180
  if(latmax>90) latmax <- 90
  if(lonmax>180) lonmax <- 180

  metadata$time <- as.POSIXct(paste(metadata$object_date, metadata$object_time))

  sf_use_s2(FALSE)

  worldmap <- ne_countries(scale = 'medium', type = 'map_units', returnclass = 'sf') %>%
    st_crop(xmin=lonmin, xmax=lonmax, ymax=latmax, ymin=latmin)
  meta.point <- st_as_sf(metadata, coords=c("object_lon","object_lat"), crs=st_crs(worldmap))

  print(ggplot() +
          geom_sf(data = worldmap, color=NA, fill="gray54") +
          geom_sf(data = meta.point, size=1, aes(color=time)) +
          ggtitle("Sampling map") +
          theme_bw())

  sf_use_s2(TRUE)

  # 2. DATE and TIME
  print(ggplot(metadata, aes(x=time, y=reorder(sample_num, time, decreasing=T), color=as.factor(ghost_id))) +
          geom_point(position=position_dodge(width=0.3), size=3) +
          labs(color="Acq. number") +
          xlab(NULL) +
          ylab(NULL) +
          scale_y_discrete() +
          theme_minimal() +
          ggtitle("Number of acquisition per sample") +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                legend.position = "right"))

  # 3. METADATA
  metadata.long <- pivot_longer(metadata, 10:18) %>% arrange(time)
  id <- unique(metadata$sample_num)
  nb <- length(id)
  n <- ceiling(nb/10)
  a = 1

  for (i in 1:n) {
    b = a + 9
    if(b>nb) b = nb
    temp <- metadata.long %>% filter(sample_num %in% id[a:b])
    temp$sample_num <- gsub("_", " ", temp$sample_num )
    print(ggplot(temp, aes(x=reorder(sample_num,time), fill=as.factor(ghost_id), y=value)) +
            geom_bar(stat="identity", position=position_dodge()) +
            scale_fill_brewer(palette="Paired") +
            scale_x_discrete(labels = label_wrap(10)) +
            labs(fill="Acq. number") +
            xlab(NULL) +
            ylab(NULL) +
            theme_bw() +
            ggtitle("Metadata values per acquisition and per sample") +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 8),
                  legend.position = "right") +
            facet_wrap(~name, scales="free_y", ncol=1))
    a = b + 1 # conditions
    if(a>nb) a = nb
  }
}
