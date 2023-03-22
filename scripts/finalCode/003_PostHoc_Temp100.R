library(tidyverse)
library(rjags)
library(brms)
source("000_Initialization.R")
source("MISC_PrepModelData.R")
source("MISC_ModelSpecification.R")
source("MISC_PlotCode_temp100.R")

# Read in a example raster
my_raster <- raster::raster("../../data/climate/BIO1_2.5min.tif") %>%
  raster::crop(sf::st_bbox(basemap_wgs)) %>%
  raster::projectRaster(crs=crs_1) %>%
  raster::aggregate(fact=5)

my_raster <- raster::raster(resolution=c(100*1000, 100*1000),
                            crs=crs_1,
                            ext=raster::extent(my_raster))

# Read in the species trait data and tree
sp_traits <- read.csv("../../data/taxa/species_traits.csv")
sp_tree <- readRDS("../../output/tree_vcv.rds")

# Make the model-ready data
my_data_100_1 <- make.data(100)

# Read in the grid data
my_grid_100 <- readRDS("../../output/data/grid_100.rds") %>%
  dplyr::filter(GID %in% my_data_100_1$my.info$kept.sites)

# Load in the 100x100 km temperature model results
my_res_100_1 <- readRDS("my_res_100_1_temp.RDS")
MCMCvis::MCMCtrace(my_res_100_1, Rhat=TRUE, filename="../../figures/supplemental/SupplementalFile_S1.pdf")
my_sum_100_1 <- MCMCvis::MCMCsummary(my_res_100_1)
my_sim_mat_100_1 <- as.matrix(my_res_100_1)

response_100 <- plot_response_curves(my_sim_mat_100_1, my_data_100_1, sp_traits)

gridded_occ_dx <- grid_species_specific_occupancy(my_sim_mat_100_1, 
                                                  my_data_100_1, covariate="temp") %>%
  dplyr::left_join(my_grid_100, by="GID") %>%
  sf::st_as_sf()

raster_occ_dx <- list()
raster_occ_dx_inv <- list()
for(ss in 1:my_data_100_1$my.constants$nsp){
  my_grid <- dplyr::filter(gridded_occ_dx, SPID==ss) %>% st_cast("POLYGON")
  raster_occ_dx[[ss]] <- fasterize::fasterize(my_grid, 
                                              raster=my_raster, 
                                              field="Ratio")
  raster_occ_dx_inv[[ss]] <- raster_occ_dx[[ss]]*-1
}

# Grab the occupancy shifts for the northernmost, core, and southernmost
# sites in each species range and also produce occupancy shift maps for all species
northern_dx <- list() # upper quartile
southern_dx <- list() # lower quartile
core_dx <- list() # interquartile range
total_dx <- list()
for(i in 1:length(raster_occ_dx)){
  
  # Print what species is being worked on in the iteration
  message(paste("Processing shift for:", sp_traits$species[i], "..."))
  
  # Grab the species-specific shift data
  my_rast_sf <- as.data.frame(raster_occ_dx[[i]], xy=TRUE) %>%
    na.omit()
  
  if(i==84){
    northern_dx[[i]] <- c(SPID=i,
                          mean=NA,
                          sd=NA,
                          se=NA)
    southern_dx[[i]] <- c(SPID=i,
                          mean=NA,
                          sd=NA,
                          se=NA)
    core_dx[[i]] <- c(SPID=i,
                      mean=mean(my_rast_sf$layer),
                      sd=sd(my_rast_sf$layer),
                      se=sd(my_rast_sf$layer)/
                        sqrt(nrow(my_rast_sf)))
    total_dx[[i]] <- c(SPID=i,
                       mean=mean(my_rast_sf$layer),
                       sd=sd(my_rast_sf$layer),
                       se=sd(my_rast_sf$layer)/
                         sqrt(nrow(my_rast_sf)))
  } else{
    # Calculate and save off the northern, core, and southern shift data
    my_lat_groups <- Hmisc::cut2(my_rast_sf$y, g=3)
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
  
  # Create the occupancy shift map and save to a .pdf file
  # my_occ_dx_map <- ggplot()+
  #   geom_sf(basemap, mapping=aes(), fill=NA)+
  #   geom_tile(my_rast_sf, mapping=aes(x=x, y=y, fill=layer),
  #             alpha=0.7)+
  #   colorspace::scale_fill_continuous_divergingx(palette="Zissou 1", rev=TRUE,
  #                                                na.value=NA, labels=scales::percent,
  #                                                name="Shift in 50-year\nOccupancy Probability")+
  #   ggtitle(sp_traits[i,]$species)+
  #   theme_map()+
  #   theme(legend.position=c(0.05, 0.3),
  #         plot.title=element_text(face="italic"))
  # ggsave2(paste0("../../figures/supplemental/rangeMaps/", sp_traits[i,]$binomial, ".pdf"),
  #         dpi=400, height=6, width=8)
}

# Compile the trend data and save off the trends. ##################################################
# Merge the trend results for each geographic context.
total_dx_df <- do.call(rbind, total_dx) %>% as.data.frame() %>%
  dplyr::arrange(mean) %>%
  dplyr::mutate(ordering=row_number(),
                geo="total") %>%
  dplyr::left_join(sp_traits, by="SPID") %>%
  dplyr::mutate(code=stringr::str_to_upper(paste0(substr(species, 1, 3),
                                                  substr(word(species, 2), 1, 3)))) %>%
  dplyr::mutate(rangeTempBin=ifelse(rangeTemp <= quantile(.$rangeTemp, probs=0.25),
                                    "Coldest",
                                    ifelse(rangeTemp >= quantile(.$rangeTemp, probs=0.75),
                                           "Warmest", "Average"))) %>%
  dplyr::mutate(rangeTempBin=factor(rangeTempBin, levels=c("Coldest", "Average", "Warmest")))

core_dx_df <- do.call(rbind, core_dx) %>% as.data.frame()
core_dx_df <- core_dx_df[match(total_dx_df$SPID, core_dx_df$SPID),] %>%
  dplyr::mutate(ordering=row_number(),
                geo="core") %>%
  dplyr::left_join(sp_traits, by="SPID") %>%
  dplyr::mutate(code=stringr::str_to_upper(paste0(substr(species, 1, 3),
                                                  substr(word(species, 2), 1, 3)))) %>%
  dplyr::mutate(rangeTempBin=ifelse(rangeTemp <= quantile(.$rangeTemp, probs=0.25),
                                    "Coldest",
                                    ifelse(rangeTemp >= quantile(.$rangeTemp, probs=0.75),
                                           "Warmest", "Average"))) %>%
  dplyr::mutate(rangeTempBin=factor(rangeTempBin, levels=c("Coldest", "Average", "Warmest")))

northern_dx_df <- do.call(rbind, northern_dx) %>% as.data.frame()
northern_dx_df <- northern_dx_df[match(total_dx_df$SPID, northern_dx_df$SPID),] %>%
  dplyr::mutate(ordering=row_number(),
                geo="north") %>%
  dplyr::left_join(sp_traits, by="SPID") %>%
  dplyr::mutate(code=stringr::str_to_upper(paste0(substr(species, 1, 3),
                                                  substr(word(species, 2), 1, 3)))) %>%
  dplyr::mutate(rangeTempBin=ifelse(rangeTemp <= quantile(.$rangeTemp, probs=0.25),
                                    "Coldest",
                                    ifelse(rangeTemp >= quantile(.$rangeTemp, probs=0.75),
                                           "Warmest", "Average"))) %>%
  dplyr::mutate(rangeTempBin=factor(rangeTempBin, levels=c("Coldest", "Average", "Warmest")))

southern_dx_df <- do.call(rbind, southern_dx) %>% as.data.frame()
southern_dx_df <- southern_dx_df[match(total_dx_df$SPID, southern_dx_df$SPID),] %>%
  dplyr::mutate(ordering=row_number(),
                geo="south") %>%
  dplyr::left_join(sp_traits, by="SPID") %>%
  dplyr::mutate(code=stringr::str_to_upper(paste0(substr(species, 1, 3),
                                                  substr(word(species, 2), 1, 3)))) %>%
  dplyr::mutate(rangeTempBin=ifelse(rangeTemp <= quantile(.$rangeTemp, probs=0.25),
                                    "Coldest",
                                    ifelse(rangeTemp >= quantile(.$rangeTemp, probs=0.75),
                                           "Warmest", "Average"))) %>%
  dplyr::mutate(rangeTempBin=factor(rangeTempBin, levels=c("Coldest", "Average", "Warmest")))


total_dx_df$diapauseStage <- factor(total_dx_df$diapauseStage,
                                    levels=c("Egg", "Larva", "Pupa", "Adult", ""),
                                    ordered=TRUE)
beepr::beep()

# Calculate the number of species where the occupancy shift is more positive in the northern
# than in the southern part of their ranges.
master_dx_df <- do.call(rbind, list(core_dx_df, northern_dx_df, southern_dx_df, total_dx_df))
sum_test_df <- left_join(southern_dx_df, northern_dx_df, by="SPID") %>%
  dplyr::mutate(southOverNorth=ifelse(mean.x > mean.y, 1, 0),
                northOverSouth=ifelse(mean.y > mean.x, 1, 0))
print(paste("Number of species more positive in North:", sum(sum_test_df$northOverSouth, na.rm=TRUE), "\n",
            "Number of species more positive in South:", sum(sum_test_df$southOverNorth, na.rm=TRUE)))

FIGURE_THREE_A <- ggplot()+
  geom_hline(yintercept=0, linetype=2)+
  geom_linerange(NULL,
                 mapping=aes(x=seq(1,90,2), 
                             ymin=-Inf, 
                             ymax=total_dx_df$mean[seq(1,90,2)]),
                 alpha=0.05)+
  geom_pointrange(total_dx_df,
                  mapping=aes(x=ordering, 
                              y=mean, 
                              ymax=mean+sd, 
                              ymin=mean-sd,
                              color=rangeTempBin),
                  size=0.3)+
  scale_color_manual(values=c("dodgerblue", "black", "firebrick1"),
                     name="Range Classification")+
  scale_x_continuous(breaks=seq(1, 90, 1), 
                     labels=total_dx_df$code,
                     expand=expansion(add=1),
                     name="Species Code")+
  scale_y_continuous(labels=scales::label_percent(add_plusses=TRUE),
                     name=" ",
                     limits=c(-0.20, 0.25))+
  theme_cowplot()+
  theme(legend.position=c(0.5,0.1),
        legend.direction="horizontal",
        plot.background=element_rect(color=NULL, fill="white"),
        axis.title=element_text(size=18),
        axis.text=element_text(size=14),
        axis.text.x=element_text(size=10, angle=90, vjust=0.45,
                                 hjust=1))

FIGURE_THREE_B <- ggplot()+
  geom_hline(yintercept=0, linetype=2)+
  geom_linerange(NULL,
                 mapping=aes(x=seq(1,90,2), 
                             ymin=-Inf, 
                             ymax=southern_dx_df$mean[seq(1,90,2)]),
                 alpha=0.05)+
  geom_pointrange(southern_dx_df,
                  mapping=aes(x=ordering, 
                              y=mean, 
                              ymax=mean+sd, 
                              ymin=mean-sd,
                              color=rangeTempBin),
                  size=0.3)+
  scale_color_manual(values=c("dodgerblue", "black", "firebrick1"),
                     name="Range Classification")+
  scale_x_continuous(breaks=seq(1, 90, 1), 
                     labels=total_dx_df$code,
                     expand=expansion(add=1),
                     name=" ")+
  scale_y_continuous(labels=scales::label_percent(add_plusses=TRUE),
                     name=" ",
                     limits=c(-0.20, 0.25))+
  theme_cowplot()+
  theme(legend.position="none",
        legend.direction="horizontal",
        plot.background=element_rect(color=NULL, fill="white"),
        axis.title=element_text(size=18),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=14))

FIGURE_THREE_C <- ggplot()+
  geom_hline(yintercept=0, linetype=2)+
  geom_linerange(NULL,
                 mapping=aes(x=seq(1,90,2), 
                             ymin=-Inf, 
                             ymax=core_dx_df$mean[seq(1,90,2)]),
                 alpha=0.05)+
  geom_pointrange(core_dx_df,
                  mapping=aes(x=ordering, 
                              y=mean, 
                              ymax=mean+sd, 
                              ymin=mean-sd,
                              color=rangeTempBin),
                  size=0.3)+
  scale_color_manual(values=c("dodgerblue", "black", "firebrick1"),
                     name="Range Classification")+
  scale_x_continuous(breaks=seq(1, 90, 1), 
                     labels=total_dx_df$code,
                     expand=expansion(add=1),
                     name=" ")+
  scale_y_continuous(labels=scales::label_percent(add_plusses=TRUE),
                     name=" ",
                     limits=c(-0.20, 0.25))+
  theme_cowplot()+
  theme(legend.position="none",
        legend.direction="horizontal",
        plot.background=element_rect(color=NULL, fill="white"),
        axis.title=element_text(size=18),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=14))

FIGURE_THREE_D <- ggplot()+
  geom_hline(yintercept=0, linetype=2)+
  geom_linerange(NULL,
                 mapping=aes(x=seq(1,90,2), 
                             ymin=-Inf, 
                             ymax=northern_dx_df$mean[seq(1,90,2)]),
                 alpha=0.05)+
  geom_pointrange(northern_dx_df,
                  mapping=aes(x=ordering, 
                              y=mean, 
                              ymax=mean+sd, 
                              ymin=mean-sd,
                              color=rangeTempBin),
                  size=0.3)+
  scale_color_manual(values=c("dodgerblue", "black", "firebrick1"),
                     name="Range Classification")+
  scale_x_continuous(breaks=seq(1, 90, 1), 
                     labels=total_dx_df$code,
                     expand=expansion(add=1),
                     name=" ")+
  scale_y_continuous(labels=scales::label_percent(add_plusses=TRUE),
                     name=" ",
                     limits=c(-0.20, 0.25))+
  theme_cowplot()+
  theme(legend.position="none",
        legend.direction="horizontal",
        plot.background=element_rect(color=NULL, fill="white"),
        axis.title=element_text(size=18),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=14))

FIGURE_THREE_LEGEND <- cowplot::get_legend(FIGURE_THREE_A)
FIGURE_THREE_SUPP <- cowplot::plot_grid(FIGURE_THREE_D+labs(title="(a) Northern Sites in Range above 45°N"),
                                        NULL,
                                        FIGURE_THREE_C+labs(title="(b) Mid-latitude Sites in Range above 45°N"), 
                                        NULL,
                                        FIGURE_THREE_B+labs(title="(c) Southern Sites in Range above 45°N"),
                                        nrow=5,
                                        rel_heights=c(1,-0.2,1,-0.2,1))
FIGURE_THREE_MAIN <- cowplot::plot_grid(FIGURE_THREE_SUPP,
                                        FIGURE_THREE_A+
                                          theme(legend.position=c(0.05,0.8),
                                                legend.direction="horizontal")+
                                          labs(title="(d) All Sites in Range above 45°N"), 
                                        nrow=2, rel_heights=c(1.25,0.75))
FIGURE_THREE <- FIGURE_THREE_MAIN+
  draw_label("Shift in Occupancy Probability from 1970s to 2010s", 
             x=0, y=0.5, vjust=1.5, angle=90, size=18)+
  theme(plot.background=element_rect(color=NULL, fill="white"))
ggsave2("../../figures/main/FIGURE_003.png", FIGURE_THREE, dpi=350, 
        height=8, width=12, units="in")

## TRAIT-BASED ANALYSES ############################################################################
# Scale the data and set reference levels for all categorical predictors
total_dx_df_model <- total_dx_df %>%
  dplyr::mutate(rangeSize_z = scale(rangeSize),
                aveWingspan_z = scale(aveWingspan),
                numReportedHostplantFamilies_z = scale(numReportedHostplantFamilies),
                diapauseStage_z = factor(diapauseStage, levels=c("Egg", "Larva", "Pupa", "Adult"), order=TRUE),
                disturbanceAffinity = ifelse(disturbanceAffinity=="Mixed", "Generalist", disturbanceAffinity)) %>%
  dplyr::mutate(disturbanceAffinity_z = factor(disturbanceAffinity, levels=c("Generalist", "Avoidant", "Associated")),
                edgeAffinity_z = factor(edgeAffinity, levels=c("Genralist", "Avoidant", "Associated")),
                canopyAffinity_z = factor(canopyAffinity, levels=c("Generalist", "Mixed", "Open")),
                voltinism_z = factor(voltinism, levels=c("Univoltine", "Multivoltine"), order=TRUE),
                rangeTemp_z = scale(rangeTemp),
                rangePrecip_z = scale(rangePrecip)) %>%
  dplyr::select(species, mean, se, sd, rangeTemp_z, 
                rangePrecip_z, rangeSize_z, aveWingspan_z, 
                numReportedHostplantFamilies_z,
                diapauseStage_z, disturbanceAffinity_z) %>%
  dplyr::filter(complete.cases(.))

sp_tree <- ape::read.tree("../../data/taxa/SupDryad_treepl.tre")
sp_tree$tip.label <- str_replace_all(sp_tree$tip.label, "_", " ")
sp_tree$tip.label <- stringr::word(sp_tree$tip.label, start=3, end=4)

sp_tree <- sp_tree %>%
  ape::keep.tip(total_dx_df$species) %>%
  ape::vcv.phylo()

# INTERCEPT ONLY MODELS (MODEL A)
MODEL_A <- brms::brm(mean~1,
                     data=total_dx_df_model,
                     iter=200000,
                     warmup=100000,
                     thin=50,
                     family=gaussian(),
                     prior=c(prior(normal(0,10), "Intercept")),
                     control=list(max_treedepth=15,
                                  adapt_delta=0.9999),
                     cores=5,
                     save_pars = save_pars(all = TRUE))
saveRDS(MODEL_A, "../../output/modelFiles/100kmTemp_ModelA.rds")

# INTERCEPT + PHYLOGENY MODELS (MODEL B)
MODEL_B <- brms::brm(mean~1+
                       (1|gr(species, cov=sp_tree)),
                     data=total_dx_df_model,
                     data2=list(sp_tree=sp_tree),
                     iter=200000,
                     warmup=100000,
                     thin=50,
                     family=gaussian(),
                     prior=c(prior(normal(0,10), "Intercept")),
                     control=list(max_treedepth=15,
                                  adapt_delta=0.9999),
                     cores=5,
                     save_pars = save_pars(all = TRUE))
saveRDS(MODEL_B, "../../output/modelFiles/100kmTemp_ModelB.rds")

# TEMPERATURE MODEL (MODEL C)
MODEL_C <- brms::brm(mean~1+
                       rangeTemp_z,
                     data=total_dx_df_model,
                     iter=200000,
                     warmup=100000,
                     thin=50,
                     family=gaussian(),
                     prior=c(prior(normal(0,10), "Intercept")),
                     control=list(max_treedepth=15,
                                  adapt_delta=0.9999),
                     cores=5,
                     save_pars = save_pars(all = TRUE))
saveRDS(MODEL_C, "../../output/modelFiles/100kmTemp_ModelC.rds")

# TEMPERATURE + PHYLOGENY MODEL (MODEL D)
MODEL_D <- brms::brm(mean~1+
                       rangeTemp_z+
                       (1|gr(species, cov=sp_tree)),
                     data=total_dx_df_model,
                     data2=list(sp_tree=sp_tree),
                     iter=200000,
                     warmup=100000,
                     thin=50,
                     family=gaussian(),
                     prior=c(prior(normal(0,10), "Intercept")),
                     control=list(max_treedepth=15,
                                  adapt_delta=0.9999),
                     cores=5,
                     save_pars = save_pars(all = TRUE))
saveRDS(MODEL_D, "../../output/modelFiles/100kmTemp_ModelD.rds")

# RANGE SIZE MODEL (MODEL E)
MODEL_E <- brms::brm(mean~1+
                       rangeSize_z,
                     data=total_dx_df_model,
                     iter=200000,
                     warmup=100000,
                     thin=50,
                     family=gaussian(),
                     prior=c(prior(normal(0,10), "Intercept")),
                     control=list(max_treedepth=15,
                                  adapt_delta=0.9999),
                     cores=5,
                     save_pars = save_pars(all = TRUE))
saveRDS(MODEL_E, "../../output/modelFiles/100kmTemp_ModelE.rds")

# RANGE SIZE + PHYLOGENY MODEL (MODEL F)
MODEL_F <- brms::brm(mean~1+
                       rangeSize_z+
                       (1|gr(species, cov=sp_tree)),
                     data=total_dx_df_model,
                     data2=list(sp_tree=sp_tree),
                     iter=200000,
                     warmup=100000,
                     thin=50,
                     family=gaussian(),
                     prior=c(prior(normal(0,10), "Intercept")),
                     control=list(max_treedepth=15,
                                  adapt_delta=0.9999),
                     cores=5,
                     save_pars = save_pars(all = TRUE))
saveRDS(MODEL_F, "../../output/modelFiles/100kmTemp_ModelF.rds")

# WINGSPAN MODEL (MODEL G)
MODEL_G <- brms::brm(mean~1+
                       aveWingspan_z,
                     data=total_dx_df_model,
                     iter=200000,
                     warmup=100000,
                     thin=50,
                     family=gaussian(),
                     prior=c(prior(normal(0,10), "Intercept")),
                     control=list(max_treedepth=15,
                                  adapt_delta=0.9999),
                     cores=5,
                     save_pars = save_pars(all = TRUE))
saveRDS(MODEL_G, "../../output/modelFiles/100kmTemp_ModelG.rds")

# WINGSPAN + PHYLOGENY MODEL (MODEL H)
MODEL_H <- brms::brm(mean~1+
                       aveWingspan_z+
                       (1|gr(species, cov=sp_tree)),
                     data=total_dx_df_model,
                     data2=list(sp_tree=sp_tree),
                     iter=200000,
                     warmup=100000,
                     thin=50,
                     family=gaussian(),
                     prior=c(prior(normal(0,10), "Intercept")),
                     control=list(max_treedepth=15,
                                  adapt_delta=0.9999),
                     cores=5,
                     save_pars = save_pars(all = TRUE))
saveRDS(MODEL_H, "../../output/modelFiles/100kmTemp_ModelH.rds")

# HOSTPLANT MODEL (MODEL I)
MODEL_I <- brms::brm(mean~1+
                       numReportedHostplantFamilies_z,
                     data=total_dx_df_model,
                     iter=200000,
                     warmup=100000,
                     thin=50,
                     family=gaussian(),
                     prior=c(prior(normal(0,10), "Intercept")),
                     control=list(max_treedepth=15,
                                  adapt_delta=0.9999),
                     cores=5,
                     save_pars = save_pars(all = TRUE))
saveRDS(MODEL_I, "../../output/modelFiles/100kmTemp_ModelI.rds")

# HOSTPLANT + PHYLOGENY MODEL (MODEL J)
MODEL_J <- brms::brm(mean~1+
                       numReportedHostplantFamilies_z+
                       (1|gr(species, cov=sp_tree)),
                     data=total_dx_df_model,
                     data2=list(sp_tree=sp_tree),
                     iter=200000,
                     warmup=100000,
                     thin=50,
                     family=gaussian(),
                     prior=c(prior(normal(0,10), "Intercept")),
                     control=list(max_treedepth=15,
                                  adapt_delta=0.9999),
                     cores=5,
                     save_pars = save_pars(all = TRUE))
saveRDS(MODEL_J, "../../output/modelFiles/100kmTemp_ModelJ.rds")

# OVERWINTERING MODEL (MODEL K)
MODEL_K <- brms::brm(mean~1+
                       rangeTemp_z+
                       diapauseStage_z+
                       rangeTemp_z:diapauseStage_z,
                     data=total_dx_df_model,
                     iter=200000,
                     warmup=100000,
                     thin=50,
                     family=gaussian(),
                     prior=c(prior(normal(0,10), "Intercept")),
                     control=list(max_treedepth=15,
                                  adapt_delta=0.9999),
                     cores=5,
                     save_pars = save_pars(all = TRUE))
saveRDS(MODEL_K, "../../output/modelFiles/100kmTemp_ModelK.rds")

# OVERWINTERING + PHYLOGENY MODEL (MODEL L)
MODEL_L <- brms::brm(mean~1+
                       rangeTemp_z+
                       diapauseStage_z+
                       rangeTemp_z:diapauseStage_z+
                     (1|gr(species, cov=sp_tree)),
                     data=total_dx_df_model,
                     data2=list(sp_tree=sp_tree),
                     iter=200000,
                     warmup=100000,
                     thin=50,
                     family=gaussian(),
                     prior=c(prior(normal(0,10), "Intercept")),
                     control=list(max_treedepth=15,
                                  adapt_delta=0.9999),
                     cores=5,
                     save_pars = save_pars(all = TRUE))
saveRDS(MODEL_L, "../../output/modelFiles/100kmTemp_ModelL.rds")

# COMPLEX MODEL (MODEL M)
MODEL_M <- brms::brm(mean~1+
                       rangeTemp_z+
                       rangeSize_z+
                       diapauseStage_z+
                       aveWingspan_z+
                       numReportedHostplantFamilies_z+
                       disturbanceAffinity_z,
                     data=total_dx_df_model,
                     iter=200000,
                     warmup=100000,
                     thin=50,
                     family=gaussian(),
                     prior=c(prior(normal(0,10), "Intercept")),
                     control=list(max_treedepth=15,
                                  adapt_delta=0.9999),
                     cores=5,
                     save_pars = save_pars(all = TRUE))
saveRDS(MODEL_M, "../../output/modelFiles/100kmTemp_ModelM.rds")

# COMPLEX + PHYLOGENY MODEL (MODEL N)
MODEL_N <- brms::brm(mean~1+
                       rangeTemp_z+
                       rangeSize_z+
                       diapauseStage_z+
                       aveWingspan_z+
                       numReportedHostplantFamilies_z+
                       disturbanceAffinity_z+
                       (1|gr(species, cov=sp_tree)),
                     data=total_dx_df_model,
                     data2=list(sp_tree=sp_tree),
                     iter=200000,
                     warmup=100000,
                     thin=50,
                     family=gaussian(),
                     prior=c(prior(normal(0,10), "Intercept")),
                     control=list(max_treedepth=15,
                                  adapt_delta=0.99999),
                     cores=5,
                     save_pars = save_pars(all = TRUE))
saveRDS(MODEL_N, "../../output/modelFiles/100kmTemp_ModelN.rds")

# TEMPERATURE BY WINGSPAN MODEL (MODEL O)
MODEL_O <- brms::brm(mean~1+
                       rangeTemp_z+
                       aveWingspan_z+
                       rangeTemp_z:aveWingspan_z,
                     data=total_dx_df_model,
                     iter=200000,
                     warmup=100000,
                     thin=50,
                     family=gaussian(),
                     prior=c(prior(normal(0,10), "Intercept")),
                     control=list(max_treedepth=15,
                                  adapt_delta=0.9999),
                     cores=5,
                     save_pars = save_pars(all = TRUE))

####################################################################################################
## COMPARE ALL OF THE MODELS TO ONE ANOTHER AND VISUALIZE
####################################################################################################

loo::loo_compare(loo::loo(MODEL_A, moment_match=TRUE),
                 loo::loo(MODEL_B, moment_match=TRUE),
                 loo::loo(MODEL_C, moment_match=TRUE),
                 loo::loo(MODEL_D, moment_match=TRUE),
                 loo::loo(MODEL_E, moment_match=TRUE),
                 loo::loo(MODEL_F, moment_match=TRUE),
                 loo::loo(MODEL_G, moment_match=TRUE),
                 loo::loo(MODEL_H, moment_match=TRUE),
                 loo::loo(MODEL_I, moment_match=TRUE),
                 loo::loo(MODEL_J, moment_match=TRUE),
                 loo::loo(MODEL_K, moment_match=TRUE),
                 loo::loo(MODEL_L, moment_match=TRUE),
                 loo::loo(MODEL_M, moment_match=TRUE),
                 loo::loo(MODEL_N, moment_match=TRUE),
                 loo::loo(MODEL_O, moment_match=TRUE))



# GRAB DRAWS
MODEL_C_DRAWS <- tidybayes::gather_draws(MODEL_C,
                                            b_Intercept, 
                                            b_rangeTemp_z, regex=TRUE)  %>%
  dplyr::mutate(.variable=str_replace(.variable, "b_", "")) %>%
  dplyr::mutate(.variable=factor(.variable, levels=c("rangeTemp_z",
                                                     "Intercept")))

MODEL_G_DRAWS <- tidybayes::gather_draws(MODEL_G,
                                         b_Intercept, 
                                         b_aveWingspan_z, regex=TRUE)  %>%
  dplyr::mutate(.variable=str_replace(.variable, "b_", "")) %>%
  dplyr::mutate(.variable=factor(.variable, levels=c("aveWingspan_z",
                                                     "Intercept")))

conditional_effects(MODEL_O)

# FIGURE FOUR ########################################
MODEL_C_LINES <- total_dx_df_model %>%
  modelr::data_grid(rangeTemp_z=seq(-3, 2, length.out=50)) %>%
  tidybayes::add_predicted_draws(MODEL_C)
MODEL_G_LINES <- total_dx_df_model %>%
  modelr::data_grid(aveWingspan_z=seq(-2, 5, length.out=50)) %>%
  tidybayes::add_predicted_draws(MODEL_G)

MODEL_C_PLOT <- ggplot(MODEL_C_LINES,
                       mapping=aes(x=rangeTemp_z, y=mean))+
  tidybayes::stat_lineribbon(mapping=aes(y=.prediction),
                             .width=c(0.95, 0.8, 0.5))+
  scale_fill_brewer(palette="Greens", name="Credible Interval")+
  geom_point(data=total_dx_df_model,
             mapping=aes(x=rangeTemp_z, y=mean),
             size=2)+
  geom_hline(yintercept=0, linetype=2)+
  scale_y_continuous(labels=scales::percent_format(),
                     name="Mean Occupancy Shift from\nthe 1970s to 2010s",
                     limits=c(-0.12, 0.12))+
  scale_x_continuous(name="Range-wide Average Annual Temperature [°C]",
                     limits=c(-3,2),
                     breaks=seq(-3,2,1),
                     labels=round((seq(-3,2,1)*sd(sp_traits$rangeTemp))+
                                    mean(sp_traits$rangeTemp)), 1)+
  theme_cowplot()+
  theme(plot.background=element_rect(fill="white",
                                     color="white"),
        legend.position=c(0.1, 0.85))

MODEL_G_PLOT <- ggplot(MODEL_G_LINES,
                       mapping=aes(x=aveWingspan_z, y=mean))+
  tidybayes::stat_lineribbon(mapping=aes(y=.prediction),
                             .width=c(0.95, 0.8, 0.5))+
  scale_fill_brewer(palette="Greens")+
  geom_point(data=total_dx_df_model,
             mapping=aes(x=aveWingspan_z, y=mean),
             size=2)+
  geom_hline(yintercept=0, linetype=2)+
  scale_y_continuous(labels=scales::percent_format(),
                     name=" ",
                     limits=c(-0.12, 0.12))+
  scale_x_continuous(name="Average Butterfly Wingspan [cm]",
                     limits=c(-2,5),
                     breaks=seq(-2,5,1),
                     labels=round((seq(-2,5,1)*sd(sp_traits$aveWingspan))+
                                    mean(sp_traits$aveWingspan)), 1)+
  theme_cowplot()+
  theme(plot.background=element_rect(fill="white",
                                     color="white"),
        legend.position="none")

cowplot::plot_grid(MODEL_C_PLOT, MODEL_G_PLOT,
                   ncol=2,
                   labels=c("(a)", "(b)"))
ggsave2("../../figures/main/FIGURE_004.png", dpi=400, height=5, width=10)

# Pagel's Lambda Estimates
hyp <- "sd_species__Intercept^2 / (sd_species__Intercept^2 + sigma^2) = 0"
(hyp <- brms::hypothesis(MODEL_D, hyp, class = NULL))

hyp <- "sd_species__Intercept^2 / (sd_species__Intercept^2 + sigma^2) = 0"
(hyp <- brms::hypothesis(MODEL_H, hyp, class = NULL))


# Other Phylogenetic Tests (outside of the BRMS model)
sp_tree <- ape::read.tree("../../data/taxa/SupDryad_treepl.tre")
sp_tree$tip.label <- str_replace_all(sp_tree$tip.label, "_", " ")
sp_tree$tip.label <- stringr::word(sp_tree$tip.label, start=3, end=4)

sp_tree <- sp_tree %>%
  ape::keep.tip(total_dx_df_model$species)
rownames(total_dx_df_model) <- total_dx_df_model$species
my_tree <- phylobase::phylo4d(sp_tree, total_dx_df_model[,2:9])

phylosignal::phyloSignal(my_tree, methods="all")

phylosim <- phylosignal::phyloSim(tree=my_tree, method="all", nsim=100, reps=99)
plot(phylosim, stacked.methods=TRUE, what="pval")

my_tree <- phylobase::phylo4d(sp_tree, total_dx_df_model[,2:9])
phylosignal::barplot.phylo4d(my_tree, 
                             center=FALSE, 
                             scale=FALSE)
mean.crig <- phylosignal::phyloCorrelogram(my_tree, trait="mean")

mean.lipa <- phylosignal::lipaMoran(my_tree)
mean.lipa.p4d <- phylosignal::lipaMoran(my_tree, as.p4d=TRUE)
phylosignal::barplot.phylo4d(my_tree, 
                             bar.col=(mean.lipa$p.value<0.05)+1, 
                             center=FALSE, 
                             scale=FALSE)



