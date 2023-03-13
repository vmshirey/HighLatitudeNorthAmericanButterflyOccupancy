# VERSION 2: Restricted to latitudinal bands that intersect with each species' range rather
# than components of each species' range.
northern_dx <- list() # upper quartile
southern_dx <- list() # lower quartile
core_dx <- list() # interquartile range
total_dx <- list()

# Grab conversion from WGS to Albers for the latitudinal range of points
conv_test <- data.frame(x=c(-100,-100,-100,-100, -100),
                        y=c(45, 55, 65, 75, 90)) %>%
  sf::st_as_sf(coords=c("x", "y"), remove=FALSE)
st_crs(conv_test) <- "wgs84"
conv_test <- conv_test %>% sf::st_transform(crs=crs_1)

for(i in 1:length(raster_occ_dx)){
  
  # Print what species is being worked on in the iteration
  message(paste("Processing shift for:", sp_traits$species[i], "..."))
  
  # Grab the species-specific shift data
  my_rast_sf <- as.data.frame(raster_occ_dx[[i]], xy=TRUE) %>%
    na.omit()
  
  # Calculate and save off the northern, core, and southern shift data
  my_lat_groups <- Hmisc::cut2(my_rast_sf$y, cuts=c(1767664,
                                                    2876653,
                                                    4483202))
  my_rast_sf$group <- my_lat_groups
  my_lat_groups <- levels(unique(my_lat_groups))
  
  northern_dx[[i]] <- c(SPID=i,
                        mean=mean(dplyr::filter(my_rast_sf, 
                                                group==my_lat_groups[3])$layer),
                        sd=sd(dplyr::filter(my_rast_sf, 
                                            group==my_lat_groups[3])$layer),
                        se=sd(dplyr::filter(my_rast_sf, 
                                            group==my_lat_groups[3])$layer)/
                          sqrt(nrow(dplyr::filter(my_rast_sf, 
                                                  group==my_lat_groups[3]))))
  southern_dx[[i]] <- c(SPID=i,
                        mean=mean(dplyr::filter(my_rast_sf, 
                                                group==my_lat_groups[1])$layer),
                        sd=sd(dplyr::filter(my_rast_sf, 
                                            group==my_lat_groups[1])$layer),
                        se=sd(dplyr::filter(my_rast_sf, 
                                            group==my_lat_groups[1])$layer)/
                          sqrt(nrow(dplyr::filter(my_rast_sf, 
                                                  group==my_lat_groups[1])))) 
  core_dx[[i]] <- c(SPID=i,
                    mean=mean(dplyr::filter(my_rast_sf, 
                                            group==my_lat_groups[2])$layer),
                    sd=sd(dplyr::filter(my_rast_sf, 
                                        group==my_lat_groups[2])$layer),
                    se=sd(dplyr::filter(my_rast_sf, 
                                        group==my_lat_groups[2])$layer)/
                      sqrt(nrow(dplyr::filter(my_rast_sf, 
                                              group==my_lat_groups[2]))))
  
  total_dx[[i]] <- c(SPID=i,
                     mean=mean(my_rast_sf$layer),
                     sd=sd(my_rast_sf$layer),
                     se=sd(my_rast_sf$layer)/
                       sqrt(nrow(my_rast_sf))) 
}
