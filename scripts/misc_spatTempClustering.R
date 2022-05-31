# Load libraries
library(tidyverse); library(data.table); library(sf)
library(cowplot); library(mapdata); library(lubridate);
library(reclin); library(sp); library(rgdal); library(geosphere)

# Load the occurrence data using the analysis
occ <- readRDS("../output/occ_sf.rds") %>%
  dplyr::filter(between(year, 1970, 2019)) %>% 
  data.table()

n_obs_year <- table(occ$year) 
single_obs_dates <- names(n_obs_year[n_obs_year == 1])
unique_dates <- names(n_obs_year[n_obs_year > 1])

cluster_lists <- list()
sin_obs_data <- occ[occ$year %in% single_obs_dates]
sin_obs_data[, "cluster" := paste0(year, "-", 1)]
cluster_lists[[1]] <- sin_obs_data
counter <- 2

for(date_use in unique_dates){
  cur_date <- occ[year == date_use]
  
  lat_lon <- cur_date[,.(decimalLongitude, decimalLatitude)]
  
  xy <- SpatialPointsDataFrame(
    lat_lon, data.frame(ID=seq(1:nrow(lat_lon))),
    proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
  
  # use the distm function to generate a geodesic distance matrix in meters
  mdist <- distm(xy)
  
  # cluster all points using a hierarchical clustering approach
  hc <- hclust(as.dist(mdist), method="complete")
  
  # define the distance threshold
  d=1000
  
  cur_date[, "cluster" := paste0(year, "-", cutree(hc, h=d))]
  
  cluster_lists[[counter]] <- cur_date
  
  counter <- counter+1
}

all_clusters <- rbindlist(cluster_lists)
size_clusters <- table(all_clusters$cluster) %>% table()

