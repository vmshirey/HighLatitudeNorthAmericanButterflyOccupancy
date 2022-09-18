# Load required libraries
library(tidyverse); library(data.table); library(sf); library(sp); library(rgdal) 
library(geosphere); library(nimble); library(raster); library(lubridate); library(ggpubr) 
library(cowplot); library(scatterpie); library(ggforce); library(exactextractr); library(gridExtra)
library(grid); library(parallel); library(coda); library(mcmcr); library(purrr); library(taxotools) 
library(taxize); library(terra); library(ggnewscale); library(pals); library(colorspace) 
library(jagsUI); library(biscale); library(gganimate); library(ggforce); library(pbapply) 
library(egg); library(HDInterval); library(GGally); library(MCMCvis); library(grid) 
library(ggrepel); library(ape); library(rotl); library(ggtree); library(tidytree);
library(ggridges); library(png); library(tidytree); library(gridExtra); library(gridtext)
library(sfheaders); library(phylosignal)

# Set a global seed for reproducibility
set.seed(04262022)

# Use legacy sf geometry
sf_use_s2(FALSE)

# Source my custom helper functions
source("../../../000_DataResources/scripts/helperFunctions_shirey.R")

# Set project CRS to North American Equal Area Albers Conic
crs_1 <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"

# Load the base map and reproject to the project reference system (NAEAC)
basemap <- st_read("../data/shapefile/ne_10m_land.shp") %>%
  st_crop(xmin=-180, xmax=-50,
          ymin=45, ymax=80)

wgs_crop <- st_bbox(basemap)

basemap <- basemap %>% 
  st_transform(crs_1) %>%
  st_make_valid()