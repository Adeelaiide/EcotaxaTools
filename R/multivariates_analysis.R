## Multivariate analysis 
## To add : Data to exclude 

multivariates_analysis <- function(data,path, bv.type="elli", living.only=T) {

# Set commun parameters for the maps------------------------------------------
ex = 3
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

### Data selection and process (filtering + normalisation)
###-----------------------------------------------------------------------------
dataset <- filter(data, type=="elli") #Use Ellipsoidal BV

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
               
               "Abundance per taxo" = dataset %>% group_by(sample_id, sample_num, object_lat, object_lon, object_annotation_category) %>% 
                 summarise(value = sum(AB, na.rm = TRUE), .groups = "keep") %>% spread(object_annotation_category,value,fill = 0) %>% arrange(sample_num),
               
               "Abundance per plankton group" = dataset %>% group_by(sample_id, sample_num, object_lat, object_lon, Sub_type) %>% 
                 summarise(value = sum(AB, na.rm = TRUE), .groups = "keep") %>% spread(Sub_type,value,fill = 0) %>% arrange(sample_num),
               
               "Abundance per trophic group" = dataset %>% group_by(sample_id, sample_num, object_lat, object_lon, Value) %>% 
                 summarise(value = sum(AB, na.rm = TRUE), .groups = "keep") %>% spread(Value,value,fill = 0) %>% arrange(sample_num),
               
               "Biovolume per taxo" = dataset %>% group_by(sample_id, sample_num, object_lat, object_lon, object_annotation_category) %>% 
                 summarise(value = sum(BV, na.rm = TRUE), .groups = "keep") %>% spread(object_annotation_category,value,fill = 0) %>% arrange(sample_num),
               
               "Biovolume per plankton group" = dataset %>% group_by(sample_id, sample_num, object_lat, object_lon, Sub_type) %>% 
                 summarise(value = sum(BV, na.rm = TRUE), .groups = "keep") %>% spread(Sub_type,value,fill = 0) %>% arrange(sample_num),
               
               "Biovolume per trophic group" = dataset %>% group_by(sample_id, sample_num, object_lat, object_lon, Value) %>% 
                 summarise(value = sum(BV, na.rm = TRUE), .groups = "keep") %>% spread(Value,value,fill = 0) %>% arrange(sample_num)
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
                    
                    "No normalisation" = data[,-c(1:4)],
                    
                    "Log x+1" = log1p(data[,-c(1:4)]),
                    
                    "Double cube root" = data[,-c(1:4)]^(1/4),
                    
                    "Relative abundance" = decostand(data[,-c(1:4)], method = "total") * 100,
                    
                    "Hellinger transform" = decostand(data[,-c(1:4)], method = "hellinger"),
                      
                    "Centered and scaled" = scale(data[,-c(1:4)]),
                    
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
print("Start PCA analysis...")

## To add : save graph in dedicated path

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

print("...Save graphics")
# Simple PCA plot
ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(size = 4) +
  labs(x = paste0("PCA 1 (", round(var_explained[1],1), "%) variance"),
       y = paste0("PCA 2 (", round(var_explained[2],1), "%) variance")) +
  theme_classic()
ggsave(filename= file.path(path, paste0("PCA on ",choice," ",norm_choice,".png")),
       width=297, height=210, units = "mm")

# Colored PCA
ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(size = 4,aes(colour = rgb(R,G,B))) +
  scale_colour_identity() +
  labs(x = paste0("PCA 1 (", round(var_explained[1],1), "%) variance"),
       y = paste0("PCA 2 (", round(var_explained[2],1), "%) variance")) +
  theme_classic()
ggsave(filename= file.path(path, paste0("PCA with RGB colors on ",choice," ",norm_choice,".png")),
       width=297, height=210, units = "mm")

# Biplot
ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(size = 4, aes(colour = rgb(R,G,B))) +
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
ggsave(filename= file.path(path, paste0("PCA biplot with RGB colors on ",choice," ",norm_choice,".png")),
       width=297, height=210, units = "mm")


##-------------------------------------------------
# PCoA
print("Start PCoA analysis...")

dist_choice <- dlg_list(c("Euclidean - recommanded for Hellinger normalisation","Bray-Curtis - recommanded for others normalisation"),
                        title = "Distance method for PCoA")$res
# Compute PCoA
pcoa_res <- if(dist_choice=="Euclidean - recommanded for Hellinger normalisation"){
  capscale(data_stat ~ 1, distance = "euclidean")
} else {
  capscale(data_stat ~ 1, distance = "bray")
}

scores_sites  <- scores(pcoa_res, choices=c(1,2,3), display="sites")
pcoa_df <- as.data.frame(scores_sites)

var_explained <- summary(pcoa_res)$cont$importance[2,1:3] * 100

fit <- envfit(pcoa_res, data_stat, permutations = 999)
scores_vars <- scores(fit, display="vectors")
vars_pcoa_df <- as.data.frame(scores_vars)

## To obtain better clarity on the graph, we decided to keep only the top 30% variable
## if the number of variables exceed 15. 
vars_pcoa_df$Contribution <- sqrt(vars_pcoa_df$MDS1^2 + vars_pcoa_df$MDS2^2)
if(nrow(vars_pcoa_df) > 15){  
  n_keep <- ceiling(nrow(vars_pcoa_df) * 0.30) # Number to keep = top 30%
  vars_pcoa_df <- vars_pcoa_df[order(-vars_pcoa_df$Contribution), ][1:n_keep, ]
}

## Color in RGB the difference between sites 
rgb_vals <- rescale(scores_sites)

pcoa_df$R <- rgb_vals[,1]
pcoa_df$G <- rgb_vals[,2]
pcoa_df$B <- rgb_vals[,3]


print("...Save graphics")
# Plot the PCoA
ggplot(pcoa_df, aes(MDS1, MDS2)) +
  geom_point(size=4, aes(colour = rgb(R,G,B))) +
  scale_colour_identity() +
  labs(x = paste0("PCoA 1 (", round(var_explained[1],1), "%) variance"),
       y = paste0("PCoA 2 (", round(var_explained[2],1), "%) variance")) +
  theme_classic()
ggsave(filename= file.path(path, paste0("PCoA with RGB colors on ",choice," ",norm_choice,".png")),
       width=297, height=210, units = "mm")

# Biplot
ggplot(pcoa_df, aes(MDS1, MDS2)) +
  geom_point(size = 4, aes(colour = rgb(R,G,B))) +
  scale_colour_identity() +
  geom_segment(data = vars_pcoa_df,
               aes(x = 0, y = 0, xend = MDS1, yend = MDS2),
               arrow = arrow(length = unit(0.2,"cm"))) +
  geom_text(data = vars_pcoa_df,
            aes(MDS1, MDS2, label = rownames(vars_pcoa_df)),
            size = 3) +
  labs(x = paste0("PCoA 1 (", round(var_explained[1],1), "%) variance"),
       y = paste0("PCoA 2 (", round(var_explained[2],1), "%) variance")) +
  theme_classic()
ggsave(filename= file.path(path, paste0("PCoA biplot with RGB colors on ",choice," ",norm_choice,".png")),
       width=297, height=210, units = "mm")

###--------------------------------
# Hierarchical Cluster analysis 
# We use the average distance between clusters 
print("Start Hierarchical cluster analysis...")
print("Clustering is computed on Bray-Curtis distance matrix with Ward method")

###This repeat loop allows the user to manually choose the number of cluster. We set 4 cluster as initial cut
cut = 4
repeat {
  
  # Hierarchical clustering on Bray-Curtis distance and complete method
  hc <- hclust(vegdist(data_stat, method = "bray"), method = "ward.D2")
  
  # Plotting dendrogram
  plot(hc, main = paste("Dendrogram -", cut, "clusters"),xlab = "")
  rect.hclust(hc, k = cut, border = "red")
  
  # Dialog box Y/N
  choix <- dlg_message(paste("Do you want to keep", cut, "clusters ?"), type = "yesno")$res
  
  if (choix == "yes") {
    #Save the dendrogram
    png(filename = file.path(path, paste0("Dendrogram with ",cut," clusters.png")),
        width = 1400, height = 1200, res = 120)
    plot(hc, main = paste("Dendrogram -", cut, "clusters"), xlab = "")
    rect.hclust(hc, k = cut, border = "red")
    dev.off()
    break
    }
  
  if (choix == "no") {
    
    new_cut <- dlg_input("How many clusters do you want ?",default = cut)$res
    
    # Check if value is numerical
    if (!is.na(as.numeric(new_cut)) && as.numeric(new_cut) > 1) {
      cut <- as.numeric(new_cut)
    } else {
      dlg_message( "Veuillez entrer un nombre valide supérieur à 1.",type = "ok")
    }
  }
}

#Extract cluster value
clusters <- cutree(hc, k=cut)
pca_df$cluster <- factor(clusters)
pcoa_df$cluster <- factor(clusters)

print("...Plot multivariate analysis with clusters")

## Plot PCA, PCoA with cluster values 

ggplot(pca_df, aes(PC1, PC2, colour=cluster)) +
  scale_colour_identity() +
  geom_point(size=4) +
  labs(x = paste0("PCA 1 (", round(var_explained[1],1), "%) variance"),
       y = paste0("PCA 2 (", round(var_explained[2],1), "%) variance")) +
  theme_classic()
ggsave(filename= file.path(path, paste0("PCA with clusters colors on ",choice," ",norm_choice,".png")),
       width=297, height=210, units = "mm")

ggplot(pcoa_df, aes(MDS1, MDS2, colour=cluster)) +
  scale_colour_identity() +
  geom_point(size=4) +
  labs(x = paste0("PCoA 1 (", round(var_explained[1],1), "%) variance"),
       y = paste0("PCoA 2 (", round(var_explained[2],1), "%) variance")) +
  theme_classic()
ggsave(filename= file.path(path, paste0("PCoA with clusters colors on ",choice," ",norm_choice,".png")),
       width=297, height=210, units = "mm")

# Biplot PCA
ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(size = 4, aes(colour=cluster)) +
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
ggsave(filename= file.path(path, paste0("PCA biplot with clusters colors on ",choice," ",norm_choice,".png")),
       width=297, height=210, units = "mm")

# Biplot PCoA
ggplot(pcoa_df, aes(MDS1, MDS2)) +
  geom_point(size = 4, aes(colour=cluster)) +
  scale_colour_identity() +
  geom_segment(data = vars_pcoa_df,
               aes(x = 0, y = 0, xend = MDS1, yend = MDS2),
               arrow = arrow(length = unit(0.2,"cm"))) +
  geom_text(data = vars_pcoa_df,
            aes(MDS1, MDS2, label = rownames(vars_pcoa_df)),
            size = 3) +
  labs(x = paste0("PCoA 1 (", round(var_explained[1],1), "%) variance"),
       y = paste0("PCoA 2 (", round(var_explained[2],1), "%) variance")) +
  theme_classic()
ggsave(filename= file.path(path, paste0("PCoA biplot with clusters colors on ",choice," ",norm_choice,".png")),
       width=297, height=210, units = "mm")

##-----------------------------------------------------------------
#MAP

print("Plotting maps")
pca_map<- pca_df %>% mutate(sample_ID = data$sample_id,sample_num = data$sample_num,lon = data$object_lon, lat = data$object_lat) %>% 
  st_as_sf(coords=c("lon","lat"), crs=st_crs(worldmap),remove = FALSE) 

pcoa_map<- pcoa_df %>% mutate(sample_ID = data$sample_id,sample_num = data$sample_num,lon = data$object_lon, lat = data$object_lat) %>% 
  st_as_sf(coords=c("lon","lat"), crs=st_crs(worldmap),remove = FALSE)

## Cluster map
ggplot() +
        geom_sf(data = worldmap, color=NA, fill="gray54") +
        geom_sf(data = pca_map, size=4, aes(color= cluster)) +
        coord_sf(xlim = c(lonmin, lonmax), ylim = c(latmin, latmax), crs = st_crs(worldmap), expand = FALSE) +
        scale_colour_identity()+
        theme_bw()+
        ggtitle("Map of clusters") +
        theme(axis.text = element_text(size = 10), plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(filename= file.path(path, paste0("Cluster map on ",choice," ",norm_choice,".png")),
       width=297, height=210, units = "mm")

## Cluster map with sample number
ggplot() +
  geom_sf(data = worldmap, color=NA, fill="gray54") +
  geom_sf_text(data = pca_map, size=4, aes(label = sample_num, color= cluster, fontface = "bold")) +
  coord_sf(xlim = c(lonmin, lonmax), ylim = c(latmin, latmax), crs = st_crs(worldmap), expand = FALSE) +
  scale_colour_identity()+
  theme_bw()+
  xlab(NULL)+ylab(NULL)+
  ggtitle("Map of clusters with sample number (chronological order)") +
  theme(axis.text = element_text(size = 10), plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(filename= file.path(path, paste0("Cluster map on ",choice," ",norm_choice,".png")),
       width=297, height=210, units = "mm")

## Map with PCA RGB colors
ggplot() +
  geom_sf(data = worldmap, color=NA, fill="gray54") +
  geom_sf(data = pca_map, size=4, aes(color= rgb(R,G,B))) +
  coord_sf(xlim = c(lonmin, lonmax), ylim = c(latmin, latmax), crs = st_crs(worldmap), expand = FALSE) +
  scale_colour_identity()+
  theme_bw()+
  ggtitle("Map of sampling point colorised with PCA values") +
  theme(axis.text = element_text(size = 10), plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(filename= file.path(path, paste0("PCA with RGB colors map on ",choice," ",norm_choice,".png")),
       width=297, height=210, units = "mm")

## Map with PCoA RGB colors
ggplot() +
  geom_sf(data = worldmap, color=NA, fill="gray54") +
  geom_sf(data = pcoa_map, size=4, aes(color= rgb(R,G,B))) +
  coord_sf(xlim = c(lonmin, lonmax), ylim = c(latmin, latmax), crs = st_crs(worldmap), expand = FALSE) +
  scale_colour_identity()+
  theme_bw()+
  ggtitle("Map of sampling point colorised with PCoA values") +
  theme(axis.text = element_text(size = 10), plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(filename= file.path(path, paste0("PCoA with RGB colors map on ",choice," ",norm_choice,".png")),
       width=297, height=210, units = "mm")
}
