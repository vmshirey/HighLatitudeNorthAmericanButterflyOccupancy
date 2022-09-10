# Load libraries
library(colorspace); library(ggdist); library(hrbrthemes)
source("000_init.R")
source("002_PrepData.R")

# Read in the basemap
basemap <- st_read("../../../Resources/ShapefilesAndRasters/StateProvince/ne_10m_admin_1_states_provinces.shp") %>%
  st_crop(xmin=-180, xmax=-50,
          ymin=45, ymax=80)

basemap <- basemap %>% 
  st_transform(crs_1) %>%
  st_make_valid()

# Load species trait data
my_traits <- read.csv("../data/taxa/species_range_clim.csv")
my_traits <- my_traits %>%
  dplyr::mutate(specificEpithet=word(species, 2)) %>%
  dplyr::mutate(speciesCode=toupper(paste0(substr(species, 1, 3),
                                           substr(specificEpithet, 1, 3))))

# Load data that went into the model
my_data_100_1 <- make.data(scale=100, imputeThres=1) 
my_data_200_1 <- make.data(scale=200, imputeThres=1)

# 100 by 100 km analysis
my_res_100_1 <- readRDS("../output/samples/my_res_100_1_LONGRUN.rds")
my_sum_100_1 <- MCMCvis::MCMCsummary(my_res_100_1)

# 200 by 200 km analysis
my_res_200_1 <- readRDS("../output/samples/my_res_200_1_LONGRUN_Intercept.RDS")
my_sum_200_1 <- MCMCvis::MCMCsummary(my_res_200_1)
my_sum_200_1 <- my_sum_200_1 %>% rownames_to_column("Parameter")

####################################################################################################
# Create Figure One                                                                                #
####################################################################################################

# Read in the occurrence data
my_occ_dat <- readRDS("../output/finalOccurrences.rds") %>%
  sf::st_as_sf(coords=c("decimalLongitude", "decimalLatitude"))
st_crs(my_occ_dat) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
my_occ_dat <- my_occ_dat %>% st_crop(xmin=-180, xmax=-50,
        ymin=45, ymax=80)

my_occ_dat <- my_occ_dat %>% 
  st_transform(crs_1) %>%
  st_make_valid()

# Read in the gridded occurrence data and grid
my_grid <- readRDS("../output/data/grid_100.rds")
my_occ_grid <- readRDS("../output/data/grid_100_occur.rds") %>%
  group_by(GID) %>%
  dplyr::mutate(detection_freq=n()) %>%
  select(GID, detection_freq) %>% ungroup() %>% unique() %>%
  inner_join(my_grid) %>% st_as_sf()

# Get a frequency of occurrence records per grid cell
my_grid <- my_grid %>% sf::st_join(my_occ_dat) %>%
  dplyr::select(GID, species, year, geometry) %>%
  dplyr::group_by(GID) %>%
  dplyr::mutate(n=n()) %>% ungroup() %>%
  dplyr::select(GID, n, geometry) %>% st_as_sf()

ggplot()+
  geom_sf(my_grid, mapping=aes(fill=n), color=NA)+
  geom_sf(basemap, mapping=aes(), fill=NA)


####################################################################################################
# Create Figure Two                                                                                #
####################################################################################################

# Pull the preciperature parameter estimates per species
sims_matrix <- as.matrix(my_res_200_1)
my_temp_samp <- sims_matrix[,grepl("psi.beta.temp", colnames(sims_matrix))]

# Extract summary statistics
my_temp_mean <- apply(my_temp_samp,2,mean)
my_temp_lower <- apply(my_temp_samp,2,
                       quantile, probs=c(0.025))
my_temp_upper <- apply(my_temp_samp,2,
                       quantile, probs=c(0.975))

# Extract community statistics
my_temp_mean_comm <- mean(my_temp_samp)
my_temp_int_comm <- quantile(my_temp_samp, probs=c(0.025,0.975))

# Merge the temperature parameter estimates with the species trait data
my_temp <- data.frame(mean=my_temp_mean, lower=my_temp_lower, upper=my_temp_upper) %>%
  cbind(my_traits) %>%
  arrange(mean) %>%
  dplyr::mutate(ordering=row_number()) %>%
  dplyr::mutate(trend=ifelse(sign(lower)==-1 & sign(upper)==1,0,
                             ifelse(sign(lower)==-1 & sign(upper)==-1,-1,1))) %>%
  dplyr::mutate(temperatureClass=cut(ave_temp2, 4))

# Create sample-style dataframe (for distribution plotting)
my_temp_samp_df <- my_temp_samp
colnames(my_temp_samp_df) <- my_traits$speciesCode
my_temp_samp_df <- melt(my_temp_samp_df)
my_temp_samp_df <- my_temp_samp_df %>% left_join(my_temp, by=c("Var2"="speciesCode"))

# Plot the temperature parameter estimates by species
tempEffectsPlot <- ggplot()+
  geom_rect(mapping=aes(xmin=my_temp_int_comm[1], xmax=my_temp_int_comm[2], 
                        ymin=-Inf, ymax=Inf), fill="grey92")+
  geom_vline(xintercept=my_temp_mean_comm, linetype=1, color="grey72")+
  geom_vline(xintercept=0, linetype=2, color="black")+
  geom_pointrange(my_temp, 
                  mapping=aes(y=ordering,
                              x=mean, 
                              xmin=lower, 
                              xmax=upper, 
                              group=SPID,
                              color=ave_temp2), alpha=0.9)+
  geom_point(my_temp,
             mapping=aes(y=ordering,
                         x=mean), pch=21, fill=NA, color="black", alpha=0.9)+
  scale_y_continuous(breaks=c(1:nrow(my_traits)), labels=my_temp$speciesCode)+
  scale_color_continuous_divergingx(mid=mean(my_temp$ave_temp2),
                                    palette="temps", name="Average Annual Range-\nwide Temperature (°C)",
                                    guide = guide_colorbar(
                                      direction = "horizontal",
                                      title.position = "top"
                                    ))+
  labs(y="Species Code", x="Effect of Rising Minimum Temperature (°C) on\nOccupancy Probability")+
  theme_cowplot()+
  theme(legend.position="top",
        axis.text.y=element_text(size=9),
        legend.box="vertical",
        legend.key.width=unit(1, 'cm'),
        legend.background=element_rect(fill=alpha("white", 0.9)),
        plot.background=element_rect(fill="white"))
ggsave2("../figures/tempEffectsPlot.png", tempEffectsPlot, dpi=400, height=11, width=6)

# Pull the precipitation parameter estimates per species
sims_matrix <- as.matrix(my_res_200_1)
my_precip_samp <- sims_matrix[,grepl("psi.beta.precip", colnames(sims_matrix))]
my_precip_mean <- apply(my_precip_samp,2,mean)
my_precip_lower <- apply(my_precip_samp,2,
                       quantile, probs=c(0.025))
my_precip_upper <- apply(my_precip_samp,2,
                       quantile, probs=c(0.975))

# Extract community statistics
my_precip_mean_comm <- mean(my_precip_samp)
my_precip_int_comm <- quantile(my_precip_samp, probs=c(0.025,0.975))

# Merge the precipitation parameter estimates with the species trait data
my_precip <- data.frame(mean=my_precip_mean, lower=my_precip_lower, upper=my_precip_upper) %>%
  cbind(my_traits) %>%
  arrange(mean) %>%
  dplyr::mutate(ordering=row_number()) %>%
  dplyr::mutate(trend=ifelse(sign(lower)==-1 & sign(upper)==1,0,
                             ifelse(sign(lower)==-1 & sign(upper)==-1,-1,1))) %>%
  dplyr::mutate(preciperatureClass=cut(ave_precip2, 4))

# Plot the precipitation parameter estimates by species
precipEffectsPlot <- ggplot()+
  geom_rect(mapping=aes(xmin=my_precip_int_comm[1], xmax=my_precip_int_comm[2], 
                        ymin=-Inf, ymax=Inf), fill="grey92")+
  geom_vline(xintercept=my_precip_mean_comm, linetype=1, color="grey72")+
  geom_vline(xintercept=0, linetype=2, color="black")+
  geom_pointrange(my_precip, 
                  mapping=aes(y=ordering,
                              x=mean, 
                              xmin=lower, 
                              xmax=upper, 
                              group=SPID,
                              color=ave_precip2))+
  geom_point(my_precip,
             mapping=aes(y=ordering,
                         x=mean), pch=21, fill=NA, color="black")+
  scale_y_continuous(breaks=c(1:nrow(my_traits)), labels=my_precip$speciesCode)+
  scale_color_continuous_divergingx(mid=mean(my_precip$ave_precip2),
                                    palette="Earth", name="Average Annual Range-\nwide Precipitation (cm)",
                                    guide = guide_colorbar(
                                      direction = "horizontal",
                                      title.position = "top"
                                    ))+
  labs(y="Species Code", x="Estimated Effect of Rising Precipitation (cm) on\nOccupancy Probability")+
  theme_cowplot()+
  theme(legend.position="top",
        axis.text.y=element_text(size=9),
        legend.box="vertical",
        legend.key.width=unit(1, 'cm'),
        legend.background=element_rect(fill=alpha("white", 0.9)),
        plot.background=element_rect(fill="white"))
ggsave2("../figures/precipEffectsPlot.png", precipEffectsPlot, dpi=400, height=11, width=5)

figure_two <- cowplot::plot_grid(tempEffectsPlot, precipEffectsPlot, ncol=2,
                                 labels=c("(a)", "(b)"))
ggsave2("../figures/figure_two.png", figure_two, dpi=400, height=11, width=10)

####################################################################################################
# Create Figure Three                                                                              #
####################################################################################################

sims_matrix <- as.matrix(my_res_200_1)

# Main panel figure (overall species occupancy shifts from OI 1 to OI 10)
for(i in 1:nrow(my_traits)){
  sp_sites <- my_data_200_1$my.info$range.list[i,] %>% na.omit()
  
  sp_occ_OI1 <- sims_matrix[,"mu.psi.0"]
}





