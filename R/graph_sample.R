

#' graph.sample
#'
#' Return ggplot graphics to summarise one sample.
#' 1. Text summary
#' 2. Relative AB and BV
#' 3. Trophic pyramid
#' 4. BSS and NBSS
#' 5. Relative NBSS
#' 6. Map
#'
#' @param x BSS table for one sample
#' @param final_metadata metadata table generated by "check_metadata" and used for the analysis
#' @param bv.type elli, plain or riddled biovolume
#' @param living.only TRUE by default
#' @param taxo taxonomic table generated with "add_trophiclvl" and "add_taxo"
#'
#' @return A set of graphics to summarise the sample.
#' @export
#'
#' @examples graph.sample(x=bss of one sample, final_metadata, trophic_affiliation_of_organisms.csv, bv.type="elli", living.only=TRUE)

graph.sample <- function(x, final_metadata, taxo, bv.type="elli", living.only=T) {
 
  # Select type of biovolume
  x <- filter(x, type==bv.type)

  # Unite taxonomiques
  x <- merge(x, taxo, "object_annotation_hierarchy", all.x=T)

  # Filter living
  if (living.only==T) {
    x <- filter(x, n1=="living")
    taxo <- filter(taxo, n1=="living")
  }

  # Set common color palette
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
 
 # ------------------------------------------------------------------------------
  # Resume
  div <- x %>% group_by(object_annotation_category) %>%
    summarise(AB=sum(AB, na.rm=T)) %>% select(AB) %>% vegan::diversity()
  text = paste0(
    "\n\nTotal abundance (ind.m-3):\n",
    round(sum(x$AB, na.rm=T),2),
    "\nTotal biovolume (mm3.mm-3):\n",
    round(sum(x$BV, na.rm=T),2),
    "\nShannon Index:\n",
    round(div,2))
  p1 <- ggplot() +
    annotate("text", x = 1, y=10, size=4, label = text) +
    theme_void()

  # Relative abundance
  tot <- sum(sum(x$AB, na.rm=T))

  p2 <- x %>% group_by(Sub_type) %>%
    summarise(per=sum(AB, na.rm=T)/tot*100) %>%
    ggplot(aes(x="", y=per, fill=Sub_type)) +
    geom_bar(stat="identity", width=1, size=0.15, color="black") +
    plankton_groups_colFill +
    coord_polar("y", start=0) +
    ggtitle("Relative abundance") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, vjust = -4,size = 10,face = "bold"))

  # Relative biovolume
  tot <- sum(sum(x$BV, na.rm=T))

  p3 <- x %>% group_by(Sub_type) %>%
    summarise(per=sum(BV, na.rm=T)/tot*100) %>%
    ggplot(aes(x="", y=per, fill=Sub_type)) +
    geom_bar(stat="identity", width=0.25, size=0.15, color="black") +
    plankton_groups_colFill +
    coord_polar("y", start=0) +
    ggtitle("Relative biovolume") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, vjust = -4,size = 10,face = "bold"))

  # NBSS
  p4 <- x %>% group_by(max) %>% summarise(BV=sum(BV/norm, na.rm=T)) %>%
    ggplot(aes(x=bv_to_esdum(max), y=BV)) +
     geom_point(size=2, fill="lightgrey", colour="black",shape=21) +
    scale_x_log10("Size on ESD (\u00b5m)") +
    scale_y_log10("NBSS (mm\u00b3.mm\u207B\u00b3.m\u207B\u00b3)", labels=trans_format('log10',math_format(10^.x))) +
    theme_bw() +
    ggtitle("NBSS") +
    theme(plot.title = element_text(hjust = 0.5, size = 10,face = "bold"), legend.text = element_text(size = 0.5))

  # BSS
  p5 <- x %>% group_by(max) %>% summarise(BV=sum(BV, na.rm=T)) %>%
    ggplot(aes(x=bv_to_esdum(max), y=BV)) +
    geom_point(size=2, fill="lightgrey", colour="black",shape=21) +
    scale_x_log10("Size on ESD (\u00b5m)") +
    scale_y_log10("BSS (mm\u00b3.m\u207B\u00b3)", labels=trans_format('log10',math_format(10^.x))) +
    theme_bw() +
    ggtitle("BSS") +
    theme(plot.title = element_text(hjust = 0.5, size = 10,face = "bold"), legend.text = element_text(size = 0.5))

  # BV compo/size class
  p6 <- x %>% group_by(Sub_type, max, class) %>%
    group_by(class) %>% mutate(per = BV/sum(BV, na.rm=T)*100) %>%
    group_by(class, max, Sub_type) %>% summarise(per=sum(per, na.rm=T)) %>%
    ggplot(aes(x=bv_to_esdum(max), y=per, fill=Sub_type)) +
    geom_col(position = "fill", width = 0.02) +
    plankton_groups_colFill +
    ylab("Biovolume (%)") +
    scale_x_log10("Size on ESD (\u00b5m)") +
    ggtitle("Biovolume composition per size class") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 10,face = "bold"), legend.text = element_text(size = 0.5))

  # Trophic pyramid p7
  # Create a colum with the trophic categories
  x <- x %>% mutate(categorie = case_when(
      Value == 1 ~ "Phototrophs",
      Value == 1.5 ~ "Mixotrophs",
      Value == 2 ~ "Grazers",
      Value == 2.5 ~ "Omnivorous",
      Value == 3 ~ "Predators",
      Value == 3.5 ~ "Unknown",
      Value == -1 ~ "None"))

#Create color map for each trophic category
  color_map_troph <- c( 
  "Unknown" = "#E6E6E6",
  "Predators" = "#C75426", 
  "Omnivorous" = "#FFEBA8", 
  "Grazers" = "#4D8ABA", 
  "Mixotrophs" = "#A0C487",
  "Phototrophs" = "#66B064",
  "None" = "#555555")

 # Calculus of the required variables for geom_rect: xmin, ymin, xmax, and ymax
 plot_data <- x %>%
  mutate(trophic_level_num = Value,half_width = log(BV + 1), 
    xmin = -half_width,
    xmax = half_width,
    ymin = trophic_level_num - 0.25,
    ymax = trophic_level_num + 0.25,
    trophic_group_name = factor(categorie, levels = names(color_map_troph))) %>%
  arrange(trophic_level_num)

p7<-ggplot(plot_data) +
  geom_rect(aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax, fill = trophic_group_name)) +
  scale_fill_manual(values = color_map_troph, name = "Trophic groups") +
  scale_y_continuous(breaks = plot_data$trophic_level_num, labels = as.character(plot_data$trophic_group_name)) +
  labs(x = "Log Biovolume +1 (mm³⋅m⁻³)", y = NULL) +
  ggtitle("Trophic pyramid") + 
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 10,face = "bold"), legend.text = element_text(size = 0.5))

  # Map
   ex = 10
  latmin <- min(final_metadata$object_lat, na.rm=T)-ex
  lonmin <- min(final_metadata$object_lon, na.rm=T)-ex
  latmax <- max(final_metadata$object_lat, na.rm=T)+ex
  lonmax <- max(final_metadata$object_lon, na.rm=T)+ex

  if(latmin<(-90)) latmin <- (-90)
  if(lonmin<(-180)) lonmin <- -180
  if(latmax>90) latmax <- 90
  if(lonmax>180) lonmax <- 180

  final_metadata$time <- as.POSIXct(paste(final_metadata$object_date, final_metadata$object_time))

  sf_use_s2(FALSE)

 bbox_area <- st_bbox(c(xmin = lonmin, ymin = latmin, xmax = lonmax, ymax = latmax), crs = 4326)
 
  worldmap <- ne_countries(scale = 'medium', type = 'map_units', returnclass = 'sf') %>%
              st_filter(st_as_sfc(bbox_area))
 #st_crop(xmin=lonmin, xmax=lonmax, ymax=latmax, ymin=latmin)
  meta.x <- filter(final_metadata, sample_id==unique(x$sample_id))
  meta.point <- st_as_sf(meta.x, coords=c("object_lon","object_lat"), crs=st_crs(worldmap))

p8 <-ggplot() +
          geom_sf(data = worldmap, color=NA, fill="gray54") +
          geom_sf(data = meta.point, size=1, color="red") +
          coord_sf(xlim = c(lonmin, lonmax), ylim = c(latmin, latmax), crs = st_crs(worldmap), expand = FALSE) +
          ggtitle("Sampling map") +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5, size = 10,face = "bold"))

  sf_use_s2(TRUE)
 
  ptot <- ggarrange(p1, p2, p3, p7, p5, p4, p6, p8,
                 common.legend = T, legend="bottom",
                 ncol=4, nrow=2) %>%
    annotate_figure(top=unique(x$sample_id), fig.lab.face="bold")

  # sf_use_s2(TRUE)

  return(ptot)
}
