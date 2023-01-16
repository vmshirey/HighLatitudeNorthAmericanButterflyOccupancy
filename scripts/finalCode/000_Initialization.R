# Load required libraries
library(tidyverse)
library(data.table)
library(sf)
library(raster)
library(terra)
library(ape)
library(exactextractr)
library(cowplot)
library(ggridges)

# Set a global seed for reproducibility
set.seed(04262022)

# Use legacy sf geometry
sf_use_s2(FALSE)

# Source my custom helper functions
source("../../../../000_DataResources/scripts/helperFunctions_shirey.R")
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# Set project CRS to North American Equal Area Albers Conic
crs_1 <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"

# Load the basemap and reproject to the project reference system (NAEAC)
basemap_wgs <- st_read("../../data/shapefile/ne_10m_land.shp") %>%
  st_crop(xmin=-180, xmax=-50,
          ymin=45, ymax=80)
wgs_crop <- st_bbox(basemap_wgs)
basemap <- basemap_wgs %>% 
  st_transform(crs_1) %>%
  st_make_valid()