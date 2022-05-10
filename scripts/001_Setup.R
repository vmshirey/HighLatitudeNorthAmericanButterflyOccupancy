# Read in the species list for species to include in this analysis, remove from this list species
# that do not have a corresponding range map.
sp_list <- read.csv("../data/taxa/species_list.csv") %>%
  dplyr::filter(!is.na(rangeMapName), rangeMapName!="", rangeMapName!=" ")

# Read in the occurrence data from all three databases (GBIF, iDigBio, and SCAN). Perform some
# prefiltering as well as name harmonization with the species list.
# Import and wrangle GBIF data
gbif <- fread("../data/occurrence/gbif_occ.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t", quote="") %>%
  dplyr::mutate(species=str_replace(species, "Clossiana", "Boloria")) %>%
  dplyr::select(species, year, decimalLongitude, decimalLatitude, basisOfRecord)

# Import and wrangle iDigBio data
idig <- fread("../data/occurrence/idig_occ.csv", header=TRUE, stringsAsFactors=FALSE, sep=",") %>%
  dplyr::mutate(`dwc:specificEpithet`=word(`dwc:specificEpithet`, -1)) %>%
  dplyr::mutate(species=str_to_sentence(paste(`dwc:genus`, `dwc:specificEpithet`))) %>%
  dplyr::mutate(species=str_replace(species, "Clossiana", "Boloria")) %>%
  dplyr::select(species, year=`dwc:year`, decimalLongitude=`dwc:decimalLongitude`,
                decimalLatitude=`dwc:decimalLatitude`, basisOfRecord='dwc:basisOfRecord')

# Import and wrangle SCAN data
scan <- fread("../data/occurrence/scan_occ.csv", header=TRUE, stringsAsFactors=FALSE, sep=",") %>%
  dplyr::mutate(species=str_to_sentence(paste(genus, specificEpithet))) %>%
  dplyr::mutate(species=str_replace(species, "Clossiana", "Boloria")) %>%
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
  dplyr::filter(n >= 800) %>%
  unique() %>%
  arrange(n)

sp_kept <- sp_list_thresh$species %>% as.data.frame()
colnames(sp_kept) <- c("species")

# Load the tree used in Earl et al. 2020
my_tree <- ape::read.tree("../data/taxa/SupDryad_treepl.tre")
my_tree$tip.label <- stringr::word(my_tree$tip.label, 3, 4, sep="_") %>%
  str_replace("_", " ")

# Prune the tree to select taxa
taxa_drop <- my_tree$tip.label[my_tree$tip.label %!in% sp_kept$species]
my_tree_prune <- ape::drop.tip(my_tree, taxa_drop, root.edge=1)

sp_kept <- sp_kept[match(my_tree_prune$tip.label, sp_kept$species),] %>% as.data.frame()
colnames(sp_kept) <- c("species")
my_tree_prune$tip.label <- sp_kept$species

my_tree <- my_tree_prune

saveRDS(my_tree, "../output/tree_topology.rds")

my_vcv <- vcv(my_tree)
my_vcv <- my_vcv[sort(sp_kept$species), sort(sp_kept$species)]

saveRDS(my_vcv, "../output/tree_vcv.rds")

# Visualize the phylogenetic tree (as a sanity check before proceeding)
tree_plot <- ggtree(my_tree_prune, layout="circular", branch.length="edge.length")+
  geom_text2(aes(subset=!isTip, label=node))+
  #geom_hilight(node=120, fill="#44bb99", alpha=0.6)+ # Hesperiidae
  #geom_hilight(node=204, fill="#99ddff", alpha=0.6)+ # Lycaenidae
  #geom_hilight(node=136, fill="#eedd88", alpha=0.6)+ # Pieridae
  #geom_hilight(node=151, fill="#ffaabb", alpha=0.6)+ # Nymphalidae
  #geom_hilight(node=226, fill="#ee8866", alpha=0.6)+ # Papilionidae
  geom_tippoint(color="black")+
  geom_tiplab(size=3, color="black", hjust=-0.05)
tree_plot

ggsave2("../output/supplemental/tree.png", tree_plot, dpi=350, height=12, width=14)

occ_join <- occ_join %>%
  dplyr::filter(species %in% sp_kept$species)
saveRDS(occ_join, "../output/finalOccurrences.rds")

sp_kept <- sp_kept %>%
  dplyr::filter(species %!in% c("Pieris marginalis", "Agriades optilete")) %>%
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

# st_write(range, "../output/data/range_kept.shp", append=FALSE)

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

range2 <- range2 %>% arrange(binomial) %>%
  left_join(sp_kept, by=c("binomial"="rangeMapName")) %>%
  arrange(species)

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

range <- range %>% left_join(sp_kept, by=c("binomial"="rangeMapName")) %>%
  arrange(species)

range_ras <- fasterize::fasterize(sf=range,
                                  raster=BIO1,
                                  by="binomial")

grid_50_range <- exactextractr::exact_extract(range_ras, grid_50) %>% 
  unique()
grid_list <- list()
for(i in 1:nrow(grid_50)){
  grid_list[[i]] <- grid_50_range[[i]] %>%
    summarise(across(1:nrow(sp_kept), sum, na.rm=TRUE)) >= 1 %>%
    as.numeric()
}
grid_50_range <- do.call(rbind, grid_list)

grid_100_range <- exactextractr::exact_extract(range_ras, grid_100) %>% 
  unique()
grid_list <- list()
for(i in 1:nrow(grid_100)){
  grid_list[[i]] <- grid_100_range[[i]] %>%
    summarise(across(1:nrow(sp_kept), sum, na.rm=TRUE)) >= 1 %>%
    as.numeric()
}
grid_100_range <- do.call(rbind, grid_list)


grid_200_range <- exactextractr::exact_extract(range_ras, grid_200) %>% 
  unique()
grid_list <- list()
for(i in 1:nrow(grid_200)){
  grid_list[[i]] <- grid_200_range[[i]] %>%
    summarise(across(1:nrow(sp_kept), sum, na.rm=TRUE)) >= 1 %>%
    as.numeric()
}
grid_200_range <- do.call(rbind, grid_list)

saveRDS(grid_50_range, "../output/data/grid_50_range.rds")
saveRDS(grid_100_range, "../output/data/grid_100_range.rds")
saveRDS(grid_200_range, "../output/data/grid_200_range.rds")

# Load in climate raster data for minimum temperature and precipitation
temp_his <- stack(list.files(path="../data/climate/minTemp/historical", pattern="*.tif", full.names=TRUE))
temp <- stack(list.files(path="../data/climate/minTemp/", pattern="*.tif", full.names=TRUE))
precip <- stack(list.files(path="../data/climate/precip", pattern="*.tif", full.names=TRUE))

# Update NA values
NAvalue(temp_his) <- -90
NAvalue(temp) <- -90
NAvalue(precip) <- -90

# Resample the historical raster to include the same extent and resolution as the modern
# temperature rasters
temp_his <- raster::shift(temp_his, dy=0.5, dx=-360.5) %>%
  raster::crop(extent(-170,-55,45,75))
temp <- raster::shift(temp, dy=0, dx=-360)
temp_his <- raster::resample(temp_his, temp)

# Stack the temperature rasters
temp <- stack(temp_his, temp)

# Reproject the raster to the project CRS
temp <- projectRaster(temp, crs=crs_1)
precip <- projectRaster(precip, crs=crs_1)

# Clip both rasters to only terrestrial records
temp <- raster::mask(temp, basemap)
precip <- raster::mask(precip, basemap)

# Average climate raster data into grids
grid_50_temp <- exactextractr::exact_extract(temp, grid_50, fun="mean") %>%
  dplyr::mutate(GID=row_number())
grid_100_temp <- exactextractr::exact_extract(temp, grid_100, fun="mean") %>%
  dplyr::mutate(GID=row_number())
grid_200_temp <- exactextractr::exact_extract(temp, grid_200, fun="mean") %>%
  dplyr::mutate(GID=row_number())

grid_50_precip <- exactextractr::exact_extract(precip, grid_50, fun="mean") %>%
  dplyr::mutate(GID=row_number())
grid_100_precip <- exactextractr::exact_extract(precip, grid_100, fun="mean") %>%
  dplyr::mutate(GID=row_number())
grid_200_precip <- exactextractr::exact_extract(precip, grid_200, fun="mean") %>%
  dplyr::mutate(GID=row_number())

# Save the climate rasters to .rds files for further use as covariates in the 
# occupancy detection model
saveRDS(grid_50_temp, "../output/data/temp50.rds")
saveRDS(grid_100_temp, "../output/data/temp100.rds")
saveRDS(grid_200_temp, "../output/data/temp200.rds")
saveRDS(grid_50_precip, "../output/data/precip50.rds")
saveRDS(grid_100_precip, "../output/data/precip100.rds")
saveRDS(grid_200_precip, "../output/data/precip200.rds")

message("Finished pre-processing all data, ready for model preperations...")
