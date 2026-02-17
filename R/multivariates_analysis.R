## Multivariate analysis 
## To add : Data to exclude 

multivariates_analysis <- function(data, path, bv.type="elli", living.only=T) {}

# Set commun parameters for the maps------------------------------------------
ex = 7
latmin <- min(data$object_lat, na.rm=T)-ex
lonmin <- min(data$object_lon, na.rm=T)-ex
latmax <- max(data$object_lat, na.rm=T)+ex
lonmax <- max(data$object_lon, na.rm=T)+ex

if(latmin<(-90)) latmin <- (-90)
if(lonmin<(-180)) lonmin <- -180
if(latmax>90) latmax <- 90
if(lonmax>180) lonmax <- 180

bbox_area <- st_bbox(c(xmin = lonmin, ymin = latmin, xmax = lonmax, ymax = latmax), crs = 4326)

worldmap <- ne_countries(scale = 'medium', type = 'map_units', returnclass = 'sf') %>%
  st_filter(st_as_sfc(bbox_area))

###-----------------------------------------------------------------------------
dataset <- filter(data, type=="elli")

## Choose data to analyse
list_data <- c("Abundance per taxo",
               "Abundance per plankton group",
               "Abundance per trophic group",
               "Biovolume per taxo",
               "Biovolume per plankton group",
               "Biovolume per trophic group")

choice <- dlg_list(list_data,
                   title = "What data should be used for analysis?",
                   multiple = FALSE)$res

data <- switch(choice,
               
               "Abundance per taxo" = dataset %>% group_by(sample_id, object_lat, object_lon, object_annotation_category) %>% 
                 summarise(value = sum(AB, na.rm = TRUE), .groups = "keep") %>% spread(object_annotation_category,value,fill = 0),
               
               "Abundance per plankton group" = dataset %>% group_by(sample_id, object_lat, object_lon, Sub_type) %>% 
                 summarise(value = sum(AB, na.rm = TRUE), .groups = "keep") %>% spread(Sub_type,value,fill = 0),
               
               "Abundance per trophic group" = dataset %>% group_by(sample_id, object_lat, object_lon, Value) %>% 
                 summarise(value = sum(AB, na.rm = TRUE), .groups = "keep") %>% spread(Value,value,fill = 0),
               
               "Biovolume per taxo" = dataset %>% group_by(sample_id, object_lat, object_lon, object_annotation_category) %>% 
                 summarise(value = sum(BV, na.rm = TRUE), .groups = "keep") %>% spread(object_annotation_category,value,fill = 0),
               
               "Biovolume per plankton group" = dataset %>% group_by(sample_id, object_lat, object_lon, Sub_type) %>% 
                 summarise(value = sum(BV, na.rm = TRUE), .groups = "keep") %>% spread(Sub_type,value,fill = 0),
               
               "Biovolume per trophic group" = dataset %>% group_by(sample_id, object_lat, object_lon, Value) %>% 
                 summarise(value = sum(BV, na.rm = TRUE), .groups = "keep") %>% spread(Value,value,fill = 0)
)
               
## Choose normalisation

norm_list <- c("No normalisation",
               "Log x+1",
               "Double cube root",
               "Relative abundance",
               "Hellinger transform",
               "Centered and scaled",
               "Box-Cox")

norm_choice <- dlg_list(norm_list,
                        title = "With what normalisation?",
                        multiple = FALSE)$res

data_stat <- switch(norm_choice,
                    
                    "No normalisation" = data[,-c(1:3)],
                    
                    "Log x+1" = log1p(data[,-c(1:3)]),
                    
                    "Double cube root" = data[,-c(1:3)]^(1/4),
                    
                    "Relative abundance" = decostand(data[,-c(1:3)], method = "total") * 100,
                    
                    "Hellinger transform" = decostand(data[,-c(1:3)], method = "hellinger"),
                      
                    "Centered and scaled" = scale(data[,-c(1:3)]),
                    
                    # "Box-Cox" = {
                    #   library(MASS)
                    #   bc <- apply(data$value + 1e-7, 2, function(x) {
                    #     bc <- boxcox(x ~ 1, plotit = FALSE)
                    #     lambda <- bc$x[which.max(bc$y)]
                    #     if (lambda == 0) log(x) else (x^lambda - 1)/lambda
                    #   })
                    #   scale(bc)
                    # }
)

##---------------------------------------------------------
### Perfom PCA analysis on transformed data

## To add : Highlight grapgh, save graph in dedicated path, filter if too many species (Zoé)

pca_res <- rda(data_stat, scale = FALSE)

#Retrieve three first axis values for sites and variables 
scores_sites  <- scores(pca_res, choices=c(1,2,3), display = "sites")
scores_vars   <- scores(pca_res, choices=c(1,2,3), display = "species")
var_explained <- summary(pca_res)$cont$importance[2,1:3] * 100

pca_df <- as.data.frame(scores_sites)
vars_df <- as.data.frame(scores_vars)

## To obtain better clarity on the graph, we decided to keep only the top 30% variable
## if the number of variables exceed 15. 

vars_df$Contribution <- sqrt(vars_df$PC1^2 + vars_df$PC2^2) #Compute contribution 
if(nrow(vars_df) > 15){  
  n_keep <- ceiling(nrow(vars_df) * 0.30) # Number to keep = top 30%
  vars_df <- vars_df[order(-vars_df$Contribution), ][1:n_keep, ]
}

## Color in RGB the difference between sites 
# each PCA axis is attributed to a color (R,G,B) and value is scaled
# Therefore each station as a unique color. Color similarity indicates relative similarity

rgb_vals <- rescale(scores_sites[,1:3])
pca_df$R <- rgb_vals[,1]
pca_df$G <- rgb_vals[,2]
pca_df$B <- rgb_vals[,3]

# Simple PCA plot
ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(size = 3) +
  labs(x = paste0("PCA 1 (", round(var_explained[1],1), "%) variance"),
       y = paste0("PCA 2 (", round(var_explained[2],1), "%) variance")) +
  theme_classic()

# Colored PCA
ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(size = 3,aes(colour = rgb(R,G,B))) +
  scale_colour_identity() +
  labs(x = paste0("PCA 1 (", round(var_explained[1],1), "%) variance"),
       y = paste0("PCA 2 (", round(var_explained[2],1), "%) variance")) +
  theme_classic()


# Biplot
ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(size = 3, aes(colour = rgb(R,G,B))) +
  scale_colour_identity() +
  geom_segment(data = vars_df,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2,"cm"))) +
  geom_text(data = vars_df,
            aes(PC1, PC2, label = rownames(vars_df)),
            size = 3) +
  labs(x = paste0("PCA 1 (", round(var_explained[1],1), "%) variance"),
       y = paste0("PCA 2 (", round(var_explained[2],1), "%) variance")) +
  theme_classic()



##-------------------------------------------------
# PCoA

dist_choice <- dlg_list(c("Euclidean - recommanded for Hellinger normalisation","Bray-Curtis - recommanded for others normalisation"),
                        title = "Distance method for PCoA")$res
# Compute PCoA
pcoa_res <- if(dist_choice=="Euclidean"){
  capscale(data_stat ~ 1, distance = "euclidiean")
} else {
  capscale(data_stat ~ 1, distance = "bray")
}

scores_sites  <- scores(pcoa_res, choices=c(1,2,3), display="sites")
pcoa_df <- as.data.frame(scores_sites)

var_explained <- summary(pcoa_res)$cont$importance[2,1:3] * 100

fit <- envfit(pcoa_res, data_stat, permutations = 999)
scores_vars <- scores(fit, display="vectors")
vars_df <- as.data.frame(scores_vars)

## To obtain better clarity on the graph, we decided to keep only the top 30% variable
## if the number of variables exceed 15. 
vars_df$Contribution <- sqrt(vars_df$MDS1^2 + vars_df$MDS2^2)
if(nrow(vars_df) > 15){  
  n_keep <- ceiling(nrow(vars_df) * 0.30) # Number to keep = top 30%
  vars_df <- vars_df[order(-vars_df$Contribution), ][1:n_keep, ]
}

## Color in RGB the difference between sites 
rgb_vals <- rescale(scores_pcoa)

pcoa_df$R <- rgb_vals[,1]
pcoa_df$G <- rgb_vals[,2]
pcoa_df$B <- rgb_vals[,3]

# Plot the PCoA
ggplot(pcoa_df, aes(MDS1, MDS2)) +
  geom_point(size=3, aes(colour = rgb(R,G,B))) +
  scale_colour_identity() +
  labs(x = paste0("PCoA 1 (", round(explained[1],1), "%) variance"),
       y = paste0("PCoA 2 (", round(explained[2],1), "%) variance")) +
  theme_classic()

# Biplot
ggplot(pcoa_df, aes(MDS1, MDS2)) +
  geom_point(size = 3, aes(colour = rgb(R,G,B))) +
  scale_colour_identity() +
  geom_segment(data = vars_df,
               aes(x = 0, y = 0, xend = MDS1, yend = MDS2),
               arrow = arrow(length = unit(0.2,"cm"))) +
  geom_text(data = vars_df,
            aes(MDS1, MDS2, label = rownames(vars_df)),
            size = 3) +
  labs(x = paste0("PCoA 1 (", round(var_explained[1],1), "%) variance"),
       y = paste0("PCoA 2 (", round(var_explained[2],1), "%) variance")) +
  theme_classic()

###--------------------------------
# Hierarchical Cluster analysis 
# We use the average distance between clusters 

## To add : Better selection of cut

hc <- hclust(vegdist(data_stat, method="euclidean"), method="complete")
cut <- dendrocut_auto(hc)
clusters <- cutree(hc, h=cut)

plot(hc)
abline(h=cut)


# coupe automatique possible :
k <- 4
clusters <- cutree(hc, h=0.6)

pca_df$cluster <- factor(clusters)

ggplot(pca_df, aes(PC1, PC2, colour=cluster)) +
  scale_colour_identity() +
  geom_point(size=3) +
  theme_minimal()

##-----------------------------------------------------------------
#MAP

pca_map<- pca_df %>% mutate(sample_ID = data$sample_id,lon = data$object_lon, lat = data$object_lat) %>% 
  st_as_sf(coords=c("lon","lat"), crs=st_crs(worldmap))

ggplot(pca_map, aes(lon, lat)) +
  geom_point(aes(colour = cluster), size=3) +
  scale_colour_identity() +
  coord_quickmap() +
  theme_minimal()


#BV map
ggplot() +
        geom_sf(data = worldmap, color=NA, fill="gray54") +
        geom_sf(data = pca_map, size=3, aes(color= cluster)) +
        coord_sf(xlim = c(lonmin, lonmax), ylim = c(latmin, latmax), crs = st_crs(worldmap), expand = FALSE) +
        scale_colour_identity()
        geom_text_repel(data = sample.point,aes(X, Y, label = sample_num),
                        size = 4,max.overlaps = Inf, box.padding = 0.15, 
                        point.padding = 0.15, min.segment.length = 0.3, seed = 42) +
        scale_color_viridis_c(option = "H") +
        labs(color = "Biovolume (mm3.m-3)", x = NULL, y = NULL) +
        
        ggtitle("Biovolume of the living per sample") +
        theme_bw() +
        theme(axis.text = element_text(size = 10), plot.title = element_text(hjust = 0.5, face = "bold")))
ggsave(filename= file.path(path, "Map of biovolume.png"),
       width=297, height=210, units = "mm")


