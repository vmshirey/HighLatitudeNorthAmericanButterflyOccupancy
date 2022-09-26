# Load libraries
library(colorspace); library(ggdist); library(hrbrthemes)
library(bayestestR); library(brmstools); library(beepr)

# Source additional scripts
source("000_Init.R")
source("002_PrepData.R")
source("misc_helperFunctions.R")

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
my_res_100_1 <- readRDS("../output/samples/my_res_100_1_LONGRUN_Intercept.rds")
my_sum_100_1 <- MCMCvis::MCMCsummary(my_res_100_1)
my_sum_100_1 <- my_sum_100_1 %>% rownames_to_column("Parameter")
sims_matrix_100 <- as.matrix(my_res_100_1)

# 200 by 200 km analysis
my_res_200_1 <- readRDS("../output/samples/my_res_200_1_LONGRUN_Intercept.RDS")
my_sum_200_1 <- MCMCvis::MCMCsummary(my_res_200_1)
my_sum_200_1 <- my_sum_200_1 %>% rownames_to_column("Parameter")
sims_matrix_200 <- as.matrix(my_res_200_1)

####################################################################################################
# Grab Overall Occupancy Shifts                                                                           
####################################################################################################

# 100 x 100 kilometer analysis

# 200 x 200 kilometer analysis
occShift_200_core <- compute_occ_shift(my_data_200_1, sims_matrix_200, my_traits,
                                       site_type="core", scale="200")

ggplot()+
  geom_pointrange(occShift_200_core,
                  mapping=aes(x=order, y=site_occDx, ymin=site_occDx_lower, ymax=site_occDx_upper))+
  theme_cowplot()


####################################################################################################
# Create Figure One                                                                                #
####################################################################################################
# Read in the occurrence data
my_occ_dat <- readRDS("../output/finalOccurrences.rds")

# Get count of records per species by occupancy interval
my_occ_counts <- my_occ_dat %>%
  dplyr::arrange(species) %>%
  dplyr::group_by(species) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::select(species, n) %>% ungroup() %>% unique()


####################################################################################################
# FIGURE TWO                                                                            
####################################################################################################

plot_figure_two(sims_matrix_100, my_traits, scale="100", makeCommOnly=TRUE)
plot_figure_two(sims_matrix_200, my_traits, scale="200")

####################################################################################################
# FIGURE THREE
####################################################################################################

plot_figure_three(sims_matrix_100, my_data_100, sims_matrix_200, my_data_200, my_traits, scale="100")
plot_figure_three(sims_matrix_100, my_data_100, sims_matrix_200, my_data_200, my_traits, scale="200")

####################################################################################################
# Figure Four                                                                                      #
####################################################################################################

figure_4_100 <- plot_figure_four(NULL, my_data_100_1, scale="100", SPIDs=c(5,49,17))

####################################################################################################
# Figure Five                                                                                      #
####################################################################################################

plot_figure_five(occ_dx_100, "100")
plot_figure_five(occ_dx_200, "200")


####################################################################################################
# Figure Six                                                                                       #
####################################################################################################
temp_occ_dx <- temp_occ_dx %>%
  dplyr::mutate(agreement=ifelse(sign(meanOcc_mean_100)==sign(ForisterNABASlope),
                                 "Match", "Non-match"))

figure_6 <- ggplot()+
  geom_vline(xintercept=0, linetype=2, alpha=0.25)+
  geom_hline(yintercept=0, linetype=2, alpha=0.25)+
  geom_point(temp_occ_dx,
             mapping=aes(x=meanOcc_mean_100, y=ForisterNABASlope, 
                         color=trend, shape=agreement),
             size=2)+
  geom_label_repel(dplyr::filter(temp_occ_dx, 
                                 abs(meanOcc_mean_100)>=0.02 | abs(ForisterNABASlope)>0.05),
                   mapping=aes(x=meanOcc_mean_100, y=ForisterNABASlope, 
                               color=trend, label=speciesCode),
                   size=2)+
  scale_color_manual(labels=c("Decreasing", "Stable", "Increasing"),
                     values=c("#f8a29e", "black", "#7fbff5"),
                     guide=guide_legend(direction="Vertical",
                                        title.position="top"),
                     name="Average Site Trend")+
  scale_x_continuous(labels=scales::percent)+
  labs(y="Forister et al. 2021 NABA Slope", x="Occupancy Probability Shift")+
  theme_cowplot()+
  theme(axis.text=element_text(size=14),
        legend.position="none",
        legend.box="vertical",
        legend.key.width=unit(1, 'cm'),
        plot.background=element_rect(fill="white"))
ggsave("../figures/FIGURE_6.png", figure_6, dpi=400, height=5, width=5)
