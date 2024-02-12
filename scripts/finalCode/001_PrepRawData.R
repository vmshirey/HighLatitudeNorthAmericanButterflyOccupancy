# Read in the species list for species to include in this analysis, remove from this list species
# that do not have a corresponding range map.
sp_list <- read.csv("../../data/taxa/species_list.csv") %>%
  dplyr::filter(!is.na(rangeMapName), rangeMapName!="", rangeMapName!=" ")

# Read in the occurrence data from all three databases (GBIF, iDigBio, and SCAN). Perform some
# prefiltering as well as name harmonization with the species list.
# Import and wrangle GBIF data
gbif <- fread("../../data/occurrence/gbif_occ.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t", quote="") %>%
  dplyr::mutate(species=str_replace(species, "Clossiana", "Boloria")) %>%
  dplyr::select(species, year, decimalLongitude, decimalLatitude, basisOfRecord, coordinateUncertaintyInMeters)

# Import and wrangle iDigBio data
idig <- fread("../../data/occurrence/idig_occ.csv", header=TRUE, stringsAsFactors=FALSE, sep=",") %>%
  dplyr::mutate(`dwc:specificEpithet`=word(`dwc:specificEpithet`, -1)) %>%
  dplyr::mutate(species=str_to_sentence(paste(`dwc:genus`, `dwc:specificEpithet`))) %>%
  dplyr::mutate(species=str_replace(species, "Clossiana", "Boloria")) %>%
  dplyr::select(species, year=`dwc:year`, decimalLongitude=`dwc:decimalLongitude`,
                decimalLatitude=`dwc:decimalLatitude`, basisOfRecord='dwc:basisOfRecord', 
                coordinateUncertaintyInMeters=`dwc:coordinateUncertaintyInMeters`)

# Import and wrangle SCAN data
scan <- fread("../../data/occurrence/scan_occ.csv", header=TRUE, stringsAsFactors=FALSE, sep=",") %>%
  dplyr::mutate(species=str_to_sentence(paste(genus, specificEpithet))) %>%
  dplyr::mutate(species=str_replace(species, "Clossiana", "Boloria")) %>%
  dplyr::select(species, year, decimalLongitude, decimalLatitude, basisOfRecord, coordinateUncertaintyInMeters) %>%
  dplyr::filter(species!="", species!=" ", !is.na(species))

# Merge the occurrences into a single data frame
occ <- rbind(gbif, idig, scan)
occ <- occ %>%
  dplyr::filter(species %in% sp_list$verbatimName) %>%
  dplyr::filter(coordinateUncertaintyInMeters < 25000) %>%
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

# Filter the occurrence data such that there is a range map for it, and there are over 800
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

sp_kept <- sp_list_thresh$species %>% 
  as.data.frame() %>%
  dplyr::filter(. %!in% c("Pieris marginalis", "Agriades optilete", "Celestrina ladon", "Danaus plexippus",
                          "Boloria alaskensis", "Lethe anthedon", "Papilio polyxenes", "Phyciodes cocyta",
                          "Megisto cymela", "Vanessa cardui", "Vanessa atalanta", "Chlosyne nycteis",
                          "Heraclides cresphontes", "Polygonia comma", "Satyrium calanus", "Colias scudderii",
                          "Colias tyche", "Erebia disa"))
colnames(sp_kept) <- c("species")

# Load the tree used in Earl et al. 2020
my_tree <- ape::read.tree("../../data/taxa/SupDryad_treepl.tre")
my_tree$tip.label <- stringr::word(my_tree$tip.label, 3, 4, sep="_") %>%
  str_replace("_", " ")

# Prune the tree to select taxa
taxa_drop <- my_tree$tip.label[my_tree$tip.label %!in% sp_kept$species]
my_tree_prune <- ape::drop.tip(my_tree, taxa_drop, root.edge=1)

sp_kept <- sp_kept[match(my_tree_prune$tip.label, sp_kept$species),] %>% as.data.frame()
colnames(sp_kept) <- c("species")
my_tree_prune$tip.label <- sp_kept$species

my_tree <- my_tree_prune

saveRDS(my_tree, "../../output/tree_topology.rds")

my_vcv <- vcv(my_tree)
my_vcv <- my_vcv[sort(sp_kept$species), sort(sp_kept$species)]

saveRDS(my_vcv, "../../output/tree_vcv.rds")

occ_join <- occ_join %>%
  dplyr::filter(species %in% sp_kept$species)
saveRDS(occ_join, "../../output/finalOccurrences.rds")

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
  dplyr::select(SPID, year, decimalLongitude, decimalLatitude, 
                basisOfRecord, rangeMapName, geometry)

saveRDS(occ_sf, "../../output/occ_sf.rds")
message(paste("Modeling for", nrow(sp_kept), "unique species..."))

# Create a grid over the base map region, save these to .rds files.
grid_100 <- st_make_grid(basemap, cellsize=100*1000) %>% st_sf() %>%
  st_intersection(basemap) %>% 
  dplyr::mutate(GID=row_number()) %>% st_make_valid()

grid_200 <- st_make_grid(basemap, cellsize=200*1000) %>% st_sf() %>%
  st_intersection(basemap) %>% 
  dplyr::mutate(GID=row_number()) %>% st_make_valid()

saveRDS(grid_100, "../../output/data/grid_100.rds")
saveRDS(grid_200, "../../output/data/grid_200.rds")

# Calculate the total land-surface area for each grid cell using the basemap, save these to
# .rds files.
Area100_grid <- st_intersection(grid_100, basemap) %>%
  dplyr::mutate(area=st_area(.)) %>%
  st_drop_geometry() %>%
  pull(area)/1e6 %>%
  round(digits=2)
Area100_grid <- round(Area100_grid, 2) %>% scale()

Area200_grid <- st_intersection(grid_200, basemap) %>%
  dplyr::mutate(area=st_area(.)) %>%
  st_drop_geometry() %>%
  pull(area)/1e6
Area200_grid <- round(Area200_grid, 2) %>% scale()

saveRDS(Area100_grid, "../../output/data/grid_100_area.rds")
saveRDS(Area200_grid, "../../output/data/grid_200_area.rds")

# Read in the range data for North American species, keeping only species used in the analysis.
range <- st_read("../../data/range/na_rm_albers.shp")
range <- range %>%
  dplyr::filter(binomial %in% sp_kept$rangeMapName) %>%
  arrange(binomial)

range_buff <- st_read("../../data/range/na_rn_albers_buff.shp")
range_buff <- range_buff %>%
  dplyr::filter(binomial %in% sp_kept$rangeMapName) %>%
  arrange(binomial)

range$area <- st_area(range)
range_area <- range %>%
  dplyr::group_by(binomial) %>%
  dplyr::mutate(areaSum=sum(area)) %>%
  ungroup() %>% dplyr::select(binomial, areaSum) %>%
  st_drop_geometry() %>% unique()

# Read in Bioclim data for average annual temperature (BIO1) and precipitation (BIO12). Using
# these values, extract the mean into the species range for a sense of climatic adaptation.
# Write this to a .csv file to merge with trait data from LepTraits v1.0.
BIO1 <- raster("../../data/climate/BIO1_2.5min.tif")
crs(BIO1) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
BIO1 <- BIO1 %>% raster::crop(extent(-180, -10, 0, 90))
BIO1 <- BIO1 %>% raster::projectRaster(crs=crs_1)

BIO12 <- raster("../../data/climate/BIO12_2.5min.tif")
crs(BIO12) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
BIO12 <- BIO12 %>% raster::crop(extent(-180, -10, 0, 90))
BIO12 <- BIO12 %>% raster::projectRaster(crs=crs_1)

basemap_full <- st_read("../../data/shapefile/ne_10m_land.shp") %>%
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
                var_temp2 = diff(range(ave_temp, na.rm=TRUE)),
                ave_precip2 = mean(ave_precip, na.rm=TRUE),
                var_precip2 = diff(range(ave_precip, na.rm=TRUE)))

range2 <- range %>% st_drop_geometry() %>%
  dplyr::select(binomial, ave_temp2, var_temp2, ave_precip2, var_precip2) %>%
  unique()

range2 <- range2 %>% arrange(binomial) %>%
  left_join(sp_kept, by=c("binomial"="rangeMapName")) %>%
  arrange(species)

write.csv(range2, "../../data/taxa/species_range_clim.csv")

# Join the occurrence data to the grids and save into .rds files.
grid_100_occur <-  occ_sf %>%
  st_join(grid_100) %>%
  dplyr::select(SPID, year, GID) %>%
  dplyr::filter(!is.na(GID)) %>%
  st_drop_geometry() %>%
  unique()

grid_200_occur <-  occ_sf %>%
  st_join(grid_200) %>%
  dplyr::select(SPID, year, GID) %>%
  dplyr::filter(!is.na(GID)) %>%
  st_drop_geometry() %>%
  unique()

saveRDS(grid_100_occur, "../../output/data/grid_100_occur.rds")
saveRDS(grid_200_occur, "../../output/data/grid_200_occur.rds")

# Intersect the ranges with the grid cells to constrain analysis to grid cells that fall
# within species ranges. Write these to .rds files for later use. To make this not take
# up all of your time, convert the ranges to rasters and then sample from there.
range_ras <- fasterize::fasterize(sf=range_buff,
                                  raster=BIO1,
                                  by="binomial")

grid_100_range <- exactextractr::exact_extract(range_ras, grid_100)
grid_list <- list()
for(i in 1:nrow(grid_100)){
  grid_list[[i]] <- grid_100_range[[i]] %>%
    summarise(across(1:nrow(sp_kept), sum, na.rm=TRUE)) >= 1 %>%
    as.numeric()
}
grid_100_range <- do.call(rbind, grid_list)

ggplot()+
  geom_sf(grid_100, mapping=aes(), fill=NA)+
  geom_sf(dplyr::filter(range_buff, binomial=="Parnassius_clodius"), 
          mapping=aes(), fill="cadetblue", alpha=0.5)+
  geom_sf(dplyr::filter(grid_100[grid_100_range[,67],]), 
          mapping=aes(), fill="cadetblue")+
  theme_map()

grid_200_range <- exactextractr::exact_extract(range_ras, grid_200)
grid_list <- list()
for(i in 1:nrow(grid_200)){
  grid_list[[i]] <- grid_200_range[[i]] %>%
    summarise(across(1:nrow(sp_kept), sum, na.rm=TRUE)) >= 1 %>%
    as.numeric()
}
grid_200_range <- do.call(rbind, grid_list)

saveRDS(grid_100_range, "../../output/data/grid_100_range.rds")
saveRDS(grid_200_range, "../../output/data/grid_200_range.rds")

# Load in climate raster data for minimum temperature and precipitation
temp_his <- raster::stack(list.files(path="../../data/climate/minTemp/historical", pattern="*.tif", full.names=TRUE))
temp <- raster::stack(list.files(path="../../data/climate/minTemp/", pattern="*.tif", full.names=TRUE))
precip <- raster::stack(list.files(path="../../data/climate/precip", pattern="*.tif", full.names=TRUE))

NAvalue(temp_his) <- -46
NAvalue(temp) <- -46
NAvalue(precip) <- -1
temp[temp < -46] <- NA
temp[temp > 50] <- NA

# Resample the historical raster to include the same extent and resolution as the modern
# temperature rasters
temp_his <- raster::shift(temp_his, dy=0.5, dx=-360.5) %>%
  raster::crop(extent(-170,-55,45,75))
temp <- raster::shift(temp, dy=0, dx=-360)
temp_his <- raster::resample(temp_his, temp)

# Stack the temperature rasters
temp <- stack(temp_his, temp)

# Reproject the raster to the project CRS
temp <- raster::projectRaster(temp, crs=crs_1)
precip <- raster::projectRaster(precip, crs=crs_1)

# Clip both rasters to only terrestrial records
temp <- raster::mask(temp, basemap)
precip <- raster::mask(precip, basemap)
precip[precip < 0] <- NA
precip[precip > 10000] <- NA

# Average climate raster data into grids
grid_100_temp <- exactextractr::exact_extract(temp, grid_100, fun="mean") %>%
  dplyr::mutate(GID=row_number())
grid_200_temp <- exactextractr::exact_extract(temp, grid_200, fun="mean") %>%
  dplyr::mutate(GID=row_number())

grid_100_precip <- exactextractr::exact_extract(precip, grid_100, fun="mean") %>%
  dplyr::mutate(GID=row_number())
grid_200_precip <- exactextractr::exact_extract(precip, grid_200, fun="mean") %>%
  dplyr::mutate(GID=row_number())

# Save the climate rasters to .rds files for further use as covariates in the
# occupancy detection model
saveRDS(grid_100_temp, "../../output/data/temp100.rds")
saveRDS(grid_200_temp, "../../output/data/temp200.rds")
saveRDS(grid_100_precip, "../../output/data/precip100.rds")
saveRDS(grid_200_precip, "../../output/data/precip200.rds")

message("Finished pre-processing all data, ready for model preperations...")

# Create a summary figure showing the occurrence data, type of data, and climate change
# FIGURE ONE
temp_70s <- calc(subset(temp, subset=grep(paste(c(1970, 1971, 1972, 1973, 1974,
                                               1975, 1976, 1977, 1978, 1979), collapse="|"),
                                         names(temp), value=TRUE)),
                fun=mean)
temp_10s <- calc(subset(temp, subset=grep(paste(c(2010, 2011, 2012, 2013, 2014,
                                        2015, 2016, 2017, 2018, 2019), collapse="|"),
                                names(temp), value=TRUE)),
       fun=mean)

temp_dx <- temp_10s - temp_70s
temp_dx <- raster::mask(temp_dx, sf::as_Spatial(basemap))
temp_dx <- as(temp_dx, "SpatialPixelsDataFrame")
temp_dx <- as.data.frame(temp_dx)
colnames(temp_dx) <- c("value", "x", "y")

occ_sf_plot <- occ_sf %>%
  dplyr::slice_sample(n=5000) %>%
  sf::st_intersection(basemap)

occ_sum <- occ_sf %>%
  dplyr::group_by(year, basisOfRecord) %>%
  dplyr::mutate(freq=n()) %>%
  dplyr::ungroup() %>%
  dplyr::select(year, basisOfRecord, freq) %>%
  sf::st_drop_geometry() %>%
  unique()

FIGURE_ONE_A <- ggplot()+
  geom_tile(temp_dx,
          mapping=aes(x=x, y=y, fill=value))+
  colorspace::scale_fill_continuous_diverging(palette="Blue-Red",
                                              name="Change in Average\nAnnualMin. Temp [C]",
                                              limits=c(-3,4))+
  geom_sf(basemap,
          mapping=aes(),
          fill=NA)+
  theme_map()+
  theme(legend.position=c(0.05, 0.3),
        legend.direction="horizontal")

FIGURE_ONE_B <- ggplot()+
  geom_sf(dplyr::slice_sample(occ_sf_plot, n=5000),
          mapping=aes(color=basisOfRecord))+
  colorspace::scale_color_discrete_qualitative(name="Type of Record")+
  geom_sf(basemap,
          mapping=aes(),
          fill=NA)+
  theme_map()+
  theme(legend.position=c(0.00, 0.3),
        legend.direction="horizontal")

FIGURE_ONE_C <- ggplot()+
  geom_point(occ_sum,
             mapping=aes(x=year, y=freq, color=basisOfRecord))+
  geom_line(occ_sum,
            mapping=aes(x=year, y=freq, color=basisOfRecord))+
  colorspace::scale_color_discrete_qualitative(name="Type of Record")+
  scale_x_continuous(limits=c(1970,2019))+
  scale_y_log10(limits=c(10,15000))+
  labs(x="Year of Collection/Observation",
       y="Number of\nOccurrences")+
  theme_cowplot()+
  theme(legend.position="none")

FIGURE_ONE_B_LEGEND <- get_legend(FIGURE_ONE_B)
FIGURE_ONE_A_LEGEND <- get_legend(FIGURE_ONE_A)

FIGURE_ONE_TOP <- cowplot::plot_grid(FIGURE_ONE_A+theme(legend.position="none"), 
                                     FIGURE_ONE_B+theme(legend.position="none"),
                   ncol=2, labels=c("(a)", "(b)"))
FIGURE_ONE_TOP2 <- cowplot::plot_grid(FIGURE_ONE_TOP,
                                      cowplot::plot_grid(NULL,
                                                         FIGURE_ONE_A_LEGEND,
                                                         FIGURE_ONE_B_LEGEND,
                                                         NULL, ncol=4,
                                                         rel_widths=c(0.15, 0.35, 0.25, 0.25)),
                                      nrow=2,
                                      rel_heights=c(1,0.1))
FIGURE_ONE <- cowplot::plot_grid(FIGURE_ONE_TOP2,
                                 FIGURE_ONE_C,
                                 nrow=2, labels=c("", "(c)"))+
  theme(plot.background=element_rect(fill="white", color="white"))
ggsave2("../../figures/main/FIGURE_001.png", FIGURE_ONE, dpi=400, height=6, width=12)




