# Load libraries
library(colorspace); library(ggdist); library(hrbrthemes)
library(bayestestR); library(brmstools); library(beepr)

# Source additional scripts
source("000_Init.R")
source("002_PrepData.R")
source("misc_helperFunctions.R")
sf_use_s2(FALSE)

# Load basemap
basemap <- sf::st_read("../../../000_DataResources/shapefiles/land/ne_50m_land.shp") %>%
  st_crop(xmin=-180, xmax=-20,
          ymin=44, ymax=90) %>%
  st_transform(crs_1)
lakes <- sf::st_read("../../../000_DataResources/shapefiles/lakes/ne_50m_lakes.shp") %>%
  st_crop(xmin=-180, xmax=-20,
          ymin=44, ymax=90) %>%
  st_transform(crs_1)

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
my_res_200_1 <- readRDS("../output/samples/my_res_200_1_LONGRUN_Intercept_PCA.RDS")
my_sum_200_1 <- MCMCvis::MCMCsummary(my_res_200_1)
my_sum_200_1 <- my_sum_200_1 %>% rownames_to_column("Parameter")
sims_matrix_200 <- as.matrix(my_res_200_1)

####################################################################################################
# Grab Overall Occupancy Shifts                                                                           
####################################################################################################

# 100 x 100 kilometer analysis
occShift_100_all <- compute_occ_shift(my_data_100_1, sims_matrix_100, my_traits,
                                      site_type="all", scale="100") %>%
  dplyr::arrange(mean) %>%
  dplyr::mutate(order=row_number(),
                trend=sign(mean)) %>%
  dplyr::mutate(trend=ifelse(between(0, lower, upper), 0, trend)) %>%
  inner_join(my_traits, by="SPID")


# 200 x 200 kilometer analysis
occShift_200_all <- compute_occ_shift(my_data_200_1, sims_matrix_200, my_traits,
                                      site_type="all", scale="200") %>%
  dplyr::arrange(mean) %>%
  dplyr::mutate(order=row_number(),
                trend=sign(mean)) %>%
  dplyr::mutate(trend=ifelse(between(0, lower, upper), 0, trend)) %>%
  inner_join(my_traits, by="SPID")

occShift_200_north <- compute_occ_shift(my_data_200_1, sims_matrix_200, my_traits,
                                      site_type="north", scale="200") %>%
  dplyr::arrange(mean) %>%
  dplyr::mutate(order=row_number(),
                trend=sign(mean)) %>%
  dplyr::mutate(trend=ifelse(between(0, lower, upper), 0, trend)) %>%
  inner_join(my_traits, by="SPID")

occShift_200_south <- compute_occ_shift(my_data_200_1, sims_matrix_200, my_traits,
                                      site_type="south", scale="200") %>%
  dplyr::arrange(mean) %>%
  dplyr::mutate(order=row_number(),
                trend=sign(mean)) %>%
  dplyr::mutate(trend=ifelse(between(0, lower, upper), 0, trend)) %>%
  inner_join(my_traits, by="SPID")

####################################################################################################
# Create Figure One                                                                                #
####################################################################################################

# Figure One A: Map of all occurrence data used in this research:
my_grid <- readRDS("../output/data/grid_100.rds") %>%
  dplyr::filter(GID %in% my_data_100_1$my.info$kept.sites)

# Read in the occurrence data
my_occ_dat <- readRDS("../output/finalOccurrences.rds")
my_occ_dat_sf <- sf::st_as_sf(my_occ_dat, coords=c("decimalLongitude", "decimalLatitude"))
st_crs(my_occ_dat_sf) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
my_occ_dat_sf <- my_occ_dat_sf %>% st_transform(crs_1) %>%
  st_intersection(my_grid)

figure_1_a <- ggplot()+
  geom_sf(basemap,
          mapping=aes(), fill=NA)+
  geom_sf(dplyr::slice_sample(my_occ_dat_sf, n=10000),
          mapping=aes(color=basisOfRecord), alpha=0.5)+
  scale_color_discrete_qualitative(palette="dynamic", name="Type of Record")+
  theme_void()+
  theme(legend.position="none",
        plot.background=element_rect(fill="white", color="white"))

# Figure One B: Lineplot of records over time by the basis of each record:
# Get count of records per species by occupancy interval
my_occ_counts <- my_occ_dat %>%
  dplyr::group_by(year, basisOfRecord) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::select(year, basisOfRecord, n) %>% ungroup() %>% unique() %>%
  dplyr::filter(between(year, 1970, 2019))

figure_1_b <- ggplot()+
  geom_line(my_occ_counts,
            mapping=aes(x=year, y=n, group=basisOfRecord, color=basisOfRecord),
            size=1)+
  geom_point(my_occ_counts,
            mapping=aes(x=year, y=n, group=basisOfRecord, color=basisOfRecord),
            size=2)+
  scale_color_discrete_qualitative(palette="dynamic", name="Type of Record")+
  scale_y_log10(labels=scales::comma)+
  labs(x="Year Of Occurrence", y="Number of Records")+
  theme_cowplot()+
  theme(plot.background=element_rect(fill="white", color="white"),
        legend.position=c(0.05, 0.9))

# Figure One C: Map of net climatic change from OI_1 to OI_10:
my_grid$PCA1 <- my_data_100_1$my.data$PCA1[,10]-my_data_100_1$my.data$PCA1[,1]
my_grid$PCA2 <- my_data_100_1$my.data$PCA2[,10]-my_data_100_1$my.data$PCA2[,1]
my_grid$temp <- ((my_data_100_1$my.data$temp[,10]*
                    attr(my_data_100_1$my.info$temp, 'scaled:scale')[10]+
                    attr(my_data_100_1$my.info$temp, 'scaled:center')[10])-
                   (my_data_100_1$my.data$temp[,1]*
                      attr(my_data_100_1$my.info$temp, 'scaled:scale')[1]+
                      attr(my_data_100_1$my.info$temp, 'scaled:center')[1]))
my_grid$precip <- ((my_data_100_1$my.data$precip[,10]*
                    attr(my_data_100_1$my.info$precip, 'scaled:scale')[10]+
                    attr(my_data_100_1$my.info$precip, 'scaled:center')[10])-
                   (my_data_100_1$my.data$precip[,1]*
                      attr(my_data_100_1$my.info$precip, 'scaled:scale')[1]+
                      attr(my_data_100_1$my.info$precip, 'scaled:center')[1]))

my_grid <- my_grid %>%
  dplyr::mutate(trend=ifelse(temp>0 & precip>0, "Hotter-Wetter", 
                       ifelse(temp>0 & precip<0, "Hotter-Drier",
                              ifelse(temp<0 & precip>0, "Colder-Wetter",
                                     ifelse(temp<0 & precip<0, "Colder-Drier", "No Change")))))

basemap_int <- basemap %>%
  st_intersection(my_grid)
basemap <- basemap %>%
  st_crop(basemap_int)

figure_1_c <- ggplot()+
  geom_sf(basemap,
          mapping=aes(), fill=NA)+
  geom_sf(my_grid,
          mapping=aes(fill=as.factor(trend)), color=NA,
          alpha=0.8)+
  scale_fill_manual(values=c("Hotter-Wetter"="#9db469",
                             "Hotter-Drier"="#db9d85",
                             "Colder-Wetter"="#87aedf",
                             "Colder-Drier"="#da95cc",
                             "No Change"="grey"),
                    name="Climatic Shift\n(1970-1974 to\n2015-2019)")+
  theme_void()+
  theme(legend.position=c(0.125,0.3),
        plot.background=element_rect(fill="white", color="white"))

# Combine the panels into a single figure
figure_one_top <- cowplot::plot_grid(figure_1_a, figure_1_c, ncol=2, labels=c("(a)", "(b)"))
figure_one <- cowplot::plot_grid(figure_one_top, figure_1_b, nrow=2, labels=c("", "(c)"))
ggsave2("../figures/main/figure_one.png", figure_one, dpi=400, height=8, width=13)

####################################################################################################
# FIGURE TWO                                                                            
####################################################################################################

plot_figure_two(sims_matrix_100, my_traits, scale="100")
plot_figure_two(sims_matrix_200, my_traits, scale="200")

####################################################################################################
# Figure Three                                                                                     #
####################################################################################################

plot_figure_three(occShift_100_all, my_traits, sims_matrix_100, "100")
plot_figure_three(occShift_200_all, my_traits, sims_matrix_200, "200")


####################################################################################################
# Figure Five                                                                                      #
####################################################################################################

plot_figure_five(occShift_100_all, "100")
plot_figure_five(occShift_200_all, "200")


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
