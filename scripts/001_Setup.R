# Load required libraries
library(tidyverse); library(data.table); library(sf);
library(sp); library(rgdal); library(geosphere);
library(nimble);library(raster); library(lubridate);
library(ggpubr); library(cowplot); library(scatterpie);
library(ggforce); library(exactextractr); library(gridExtra);
library(grid); library(parallel); library(coda); library(mcmcr);
library(purrr); library(taxotools); library(taxize);
library(terra); library(ggnewscale); library(pals);
library(colorspace); library(jagsUI); library(biscale)
library(gganimate); library(nimble); library(ggforce);
library(pbapply); library(egg); library(HDInterval);
library(GGally); library(MCMCvis); library(grid); library(ggrepel);
library(fasterize); library(ape); library(rotl); library(ggtree); library(tidytree)

# Set a global seed for reproducibility
set.seed(04262022)

# Use legacy sf geometry
sf_use_s2(FALSE)

# Source my custom helper functions
source("../../../000_DataResources/scripts/helperFunctions_shirey.R")

# Set project CRS to North American Equal Area Albers Conic
crs_1 <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"

# Load the base map and reproject to the project reference system
basemap <- st_read("../data/shapefile/ne_10m_land.shp") %>%
  st_crop(xmin=-180, xmax=-50,
          ymin=45, ymax=80)

wgs_crop <- st_bbox(basemap)

basemap <- basemap %>% 
  st_transform(crs_1) %>%
  st_make_valid()

# Read in the species list for species to include in this analysis, remove from this list species
# that do not have a corresponding range map.
sp_list <- read.csv("../data/taxa/species_list.csv") %>%
  dplyr::filter(!is.na(rangeMapName), rangeMapName!="", rangeMapName!=" ")

# Read in the occurrence data from all three databases (GBIF, iDigBio, and SCAN). Perform some
# prefiltering as well as name harmonization with the species list.
# Import and wrangle GBIF data
gbif <- fread("../data/occurrence/gbif_occ.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t", quote="") %>%
  dplyr::mutate(species=str_replace(species, "Clossiana", "Boloria"),
                species=str_replace(species, "Speyeria", "Argynnis")) %>%
  dplyr::select(species, year, decimalLongitude, decimalLatitude, basisOfRecord)

# Import and wrangle iDigBio data
idig <- fread("../data/occurrence/idig_occ.csv", header=TRUE, stringsAsFactors=FALSE, sep=",") %>%
  dplyr::mutate(`dwc:specificEpithet`=word(`dwc:specificEpithet`, -1)) %>%
  dplyr::mutate(species=str_to_sentence(paste(`dwc:genus`, `dwc:specificEpithet`))) %>%
  dplyr::mutate(species=str_replace(species, "Clossiana", "Boloria"),
                species=str_replace(species, "Speyeria", "Argynnis")) %>%
  dplyr::select(species, year=`dwc:year`, decimalLongitude=`dwc:decimalLongitude`,
                decimalLatitude=`dwc:decimalLatitude`, basisOfRecord='dwc:basisOfRecord')

# Import and wrangle SCAN data
scan <- fread("../data/occurrence/scan_occ.csv", header=TRUE, stringsAsFactors=FALSE, sep=",") %>%
  dplyr::mutate(species=str_to_sentence(paste(genus, specificEpithet))) %>%
  dplyr::mutate(species=str_replace(species, "Clossiana", "Boloria"),
                species=str_replace(species, "Speyeria", "Argynnis")) %>%
  dplyr::select(species, year, decimalLongitude, decimalLatitude, basisOfRecord) %>%
  dplyr::filter(species!="", species!=" ", !is.na(species))

# Merge the occurrences into a single data frame
occ <- rbind(gbif, idig, scan)
occ <- occ %>%
  dplyr::filter(species %in% sp_list$verbatimName) %>%
  dplyr::distinct(species, year, decimalLongitude, decimalLatitude, .keep_all=TRUE) %>%
  dplyr::mutate(basisOfRecord=case_when(basisOfRecord %in% c("HUMAN_OBSERVATION", "humanobservation", "HumanObservation",
                                                             "observation", "OBSERVATION", "Saw 1", "Ruth ct. 1", "Saw 1 female",
                                                             "Saw 2", "Saw 3", "Saw 4", "Saw 7", "Saw few", "Saw many") ~ "Observation",
                                        basisOfRecord %in% c("MATERIAL_SAMPLE", "preserved specimen", "Preserved Specimen", "PRESERVED_SPECIMEN",
                                                             "PreservedSpecimen", "Caught 1", "Caught 2") ~ "Preserved Specimen",
                                        basisOfRecord %in% c("MACHINE_OBSERVATION") ~ "Observation")) %>%
  dplyr::filter(!is.na(basisOfRecord)) %>%
  dplyr::mutate(id=row_number()) %>%
  arrange(species)

# Filter the occurrence data such that there is a range map for it, and there are over 825
# occurrences. Save the occurrences to a .rds file. Rework a list of species kept for the
# analysis.
occ_join <- occ %>%
  left_join(sp_list, by=c("species"="verbatimName")) %>%
  dplyr::mutate(species=lamasListName) %>%
  dplyr::select(species, year, decimalLongitude, decimalLatitude, basisOfRecord, rangeMapName)

sp_list_thresh <- occ_join %>%
  group_by(species) %>%
  dplyr::mutate(n=n()) %>%
  ungroup() %>%
  dplyr::select(species, n) %>%
  dplyr::filter(n >= 500) %>%
  unique() %>%
  arrange(n)

occ_join <- occ_join %>%
  dplyr::filter(species %in% sp_list_thresh$species)

saveRDS(occ_join, "../output/finalOccurrences.rds")

sp_kept <- sp_list_thresh$species %>% as.data.frame()
colnames(sp_kept) <- c("species")
sp_kept <- sp_kept %>%
  arrange(species) %>%
  dplyr::mutate(SPID=row_number()) %>%
  inner_join(dplyr::select(sp_list, lamasListName, rangeMapName), 
             by=c("species"="lamasListName")) %>%
  unique()

# Convert the occurrences to spatial data and save as a .rds file.
occ_sf <- st_as_sf(occ_join, 
                   coords=c("decimalLongitude", "decimalLatitude"),
                   crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", 
                   remove=FALSE) %>%
  st_transform(crs_1)

occ_sf <- left_join(occ_sf, sp_kept) %>%
  dplyr::select(SPID, year, decimalLongitude, decimalLatitude, basisOfRecord, rangeMapName, geometry)

saveRDS(occ_sf, "../output/occ_sf.rds")

message(paste("Modeling for", nrow(sp_kept), "unique species..."))
message("Obtaining phylogentic tree from Open Tree of Life")

# Pull a phylogenetic tree from Open Tree of Life for the modeled taxa.
# Resolve names to best match (need to fix some manually, see below)
tree_taxa <- tnrs_match_names(sp_kept$species)

tree_taxa[105,]=tnrs_match_names("Papilio canadensis")
tree_taxa[106,]=tnrs_match_names("Papilio eurymedon")
tree_taxa[107,]=tnrs_match_names("Papilio rutulus")

# Pull the phylogenetic tree from OTL, clean the names a bit
my_tree <- tol_induced_subtree(ott_ids=tree_taxa$ott_id)
my_tree$tip.label <- stringr::word(my_tree$tip.label, 1,2, sep="_") %>%
  str_replace("_", " ")

# Join the tree names to the overall Lamas names list for consistency
tree_taxa[105,]$search_string <- "Pterourus canadensis"
tree_taxa[106,]$search_string <- "Pterourus eurymedon"
tree_taxa[107,]$search_string <- "Pterourus rutulus"
tree_taxa$search_string <- str_to_sentence(tree_taxa$search_string)
sp_kept <- sp_kept %>% left_join(tree_taxa, by=c("species"="search_string"))
sp_kept <- sp_kept[match(my_tree$tip.label, sp_kept$unique_name),]
my_tree$tip.label <- sp_kept$species

# Visualize the phylogenetic tree (as a sanity check before proceeding)
tree_plot <- ggtree(my_tree, layout="circular")+
  geom_hilight(node=119, fill="#eedd88", alpha=0.6)+ # Pieridae
  geom_hilight(node=182, fill="#99ddff", alpha=0.6)+ # Lycaenidae
  geom_hilight(node=134, fill="#ee8866", alpha=0.6)+ # Nymphalidae
  geom_hilight(node=220, fill="#ffaabb", alpha=0.6)+ # Papilionidae
  geom_hilight(node=205, fill="#44bb99", alpha=0.6)+ # Hesperiidae
  geom_tippoint(color="black")+
  geom_tiplab(size=3, color="black", hjust=-0.05)

ggsave2("../output/supplemental/tree.png", tree_plot, dpi=350, height=12, width=14)

# Create a grid over the base map region, save these to .rds files.
grid_50 <- st_make_grid(basemap, cellsize=50*1000) %>% st_sf() %>%
  st_intersection(basemap) %>% 
  dplyr::mutate(GID=row_number()) %>% st_make_valid()

grid_100 <- st_make_grid(basemap, cellsize=100*1000) %>% st_sf() %>%
  st_intersection(basemap) %>% 
  dplyr::mutate(GID=row_number()) %>% st_make_valid()

grid_200 <- st_make_grid(basemap, cellsize=200*1000) %>% st_sf() %>%
  st_intersection(basemap) %>% 
  dplyr::mutate(GID=row_number()) %>% st_make_valid()

saveRDS(grid_50, "../output/data/grid_50.rds")
saveRDS(grid_100, "../output/data/grid_100.rds")
saveRDS(grid_200, "../output/data/grid_200.rds")

# Calculate the total land-surface area for each grid cell using the basemap, save these to
# .rds files.
Area50_grid <- st_intersection(grid_50, basemap) %>%
  dplyr::mutate(area=st_area(.)) %>%
  st_drop_geometry() %>%
  pull(area)/1e6 %>%
  round(digits=2)
Area50_grid <- round(Area50_grid, 2)

Area100_grid <- st_intersection(grid_100, basemap) %>%
  dplyr::mutate(area=st_area(.)) %>%
  st_drop_geometry() %>%
  pull(area)/1e6 %>%
  round(digits=2)
Area100_grid <- round(Area100_grid, 2)

Area200_grid <- st_intersection(grid_200, basemap) %>%
  dplyr::mutate(area=st_area(.)) %>%
  st_drop_geometry() %>%
  pull(area)/1e6
Area200_grid <- round(Area200_grid, 2)

saveRDS(Area50_grid, "../output/data/grid_50_area.rds")
saveRDS(Area100_grid, "../output/data/grid_100_area.rds")
saveRDS(Area200_grid, "../output/data/grid_200_area.rds")

# Read in the range data for North American species, keeping only species used in the analysis.
range <- st_read("../data/range/na_rm_albers.shp")
range <- range %>%
  dplyr::filter(binomial %in% sp_kept$rangeMapName) %>%
  arrange(binomial)
st_crs(range) <- crs_1

st_write(range, "../output/data/range_kept.shp")

# Read in Bioclim data for average annual temperature (BIO1) and precipitation (BIO12). Using
# these values, extract the mean into the species range for a sense of climatic adaptation.
# Write this to a .csv file to merge with trait data from LepTraits v1.0.
BIO1 <- raster("../data/climate/BIO1_2.5min.tif")
BIO1 <- projectRaster(BIO1, crs=crs_1)
BIO12 <- raster("../data/climate/BIO12_2.5min.tif")
BIO12 <- projectRaster(BIO12, crs=crs_1)

basemap_full <- st_read("../data/shapefile/ne_10m_land.shp") %>%
  st_crop(xmin=-180, xmax=-10,
          ymin=0, ymax=90)

basemap_full <- basemap_full %>% 
  st_transform(crs_1) %>%
  st_make_valid()

BIO1 <- raster::crop(BIO1, extent(basemap_full))
BIO12 <- raster::crop(BIO12, extent(basemap_full))

ave_temp <- exact_extract(BIO1, range, fun="mean")
ave_precip <- exact_extract(BIO12, range, fun="mean")

range$ave_temp <- ave_temp
range$ave_precip <- ave_precip

range <- range %>%
  group_by(binomial) %>%
  dplyr::mutate(ave_temp2 = mean(ave_temp, na.rm=TRUE),
                ave_precip2 = mean(ave_precip, na.rm=TRUE))

range2 <- range %>% st_drop_geometry() %>%
  dplyr::select(binomial, ave_temp2, ave_precip2) %>%
  unique()

write.csv(range2, "../data/taxa/species_range_clim.csv")

# Join the occurrence data to the grids and save into .rds files.
grid_50_occur <- occ_sf %>%
  st_join(grid_50) %>% dplyr::filter(!is.na(GID))

grid_100_occur <-  occ_sf %>%
  st_join(grid_100) %>% dplyr::filter(!is.na(GID))

grid_200_occur <-  occ_sf %>%
  st_join(grid_200) %>% dplyr::filter(!is.na(GID))

saveRDS(grid_50_occur, "../output/data/grid_50_occur.rds")
saveRDS(grid_100_occur, "../output/data/grid_100_occur.rds")
saveRDS(grid_200_occur, "../output/data/grid_200_occur.rds")

# Intersect the ranges with the grid cells to constrain analysis to grid cells that fall
# within species ranges. Write these to .rds files for later use. To make this not take
# up all of your time, convert the ranges to rasters and then sample from there.
range_ras <- fasterize::fasterize(sf=range,
                                  raster=BIO1,
                                  by="binomial")

grid_50_range <- exactextractr::exact_extract(range_ras, grid_50) %>% 
  unique()
grid_list <- list()
for(i in 1:nrow(grid_50)){
  grid_list[[i]] <- grid_50_range[[i]] %>%
    summarise(across(1:115, sum, na.rm=TRUE)) >= 1 %>%
    as.numeric()
}
grid_50_range <- do.call(rbind, grid_list)

grid_100_range <- exactextractr::exact_extract(range_ras, grid_100) %>% 
  unique()
grid_list <- list()
for(i in 1:nrow(grid_100)){
  grid_list[[i]] <- grid_100_range[[i]] %>%
    summarise(across(1:115, sum, na.rm=TRUE)) >= 1 %>%
    as.numeric()
}
grid_100_range <- do.call(rbind, grid_list)


grid_200_range <- exactextractr::exact_extract(range_ras, grid_200) %>% 
  unique()
grid_list <- list()
for(i in 1:nrow(grid_200)){
  grid_list[[i]] <- grid_200_range[[i]] %>%
    summarise(across(1:115, sum, na.rm=TRUE)) >= 1 %>%
    as.numeric()
}
grid_200_range <- do.call(rbind, grid_list)

saveRDS(grid_50_range, "../output/data/grid_50_range.rds")
saveRDS(grid_100_range, "../output/data/grid_100_range.rds")
saveRDS(grid_200_range, "../output/data/grid_200_range.rds")

message("Finished pre-processing all data, ready for model preperations...")
