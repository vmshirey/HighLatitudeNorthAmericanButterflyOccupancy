library(tidyverse)
library(rjags)
source("000_Initialization.R")
source("MISC_PrepModelData.R")
source("MISC_ModelSpecification.R")
source("MISC_PlotCode_precip100.R")

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
my_res_100_1 <- readRDS("my_res_100_1_precip.RDS")
#MCMCvis::MCMCtrace(my_res_100_1, Rhat=TRUE, filename="../../figures/supplemental/SupplementalFile_S3.pdf")
my_sum_100_1 <- MCMCvis::MCMCsummary(my_res_100_1)
my_sim_mat_100_1 <- as.matrix(my_res_100_1)

response_100 <- plot_response_curves(my_sim_mat_100_1, my_data_100_1, sp_traits)

gridded_occ_dx <- grid_species_specific_occupancy(my_sim_mat_100_1, 
                                                  my_data_100_1, covariate="precip") %>%
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
  Sys.sleep(1)
  
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
  # ggsave2(paste0("../../figures/supplemental/rangeMaps/", sp_traits[i,]$binomial, "_100precip.pdf"),
  #         dpi=400, height=6, width=8)
}
core_dx_df <- do.call(rbind, core_dx) %>% as.data.frame() %>%
  dplyr::arrange(mean) %>%
  dplyr::mutate(ordering=row_number(),
                geo="core") %>%
  dplyr::left_join(sp_traits, by="SPID") %>%
  dplyr::mutate(code=stringr::str_to_upper(paste0(substr(species, 1, 3),
                                                  substr(word(species, 2), 1, 3)))) %>%
  dplyr::mutate(rangePrecipBin=ifelse(rangePrecip <= quantile(.$rangePrecip, probs=0.25),
                                      "Driest",
                                      ifelse(rangePrecip >= quantile(.$rangePrecip, probs=0.75),
                                             "Wettest", "Average"))) %>%
  dplyr::mutate(rangePrecipBin=factor(rangePrecipBin, levels=c("Driest", "Average", "Wettest")))
write.csv(core_dx_df, "../../output/trends/core_dx_df_precip100.csv")

northern_dx_df <- do.call(rbind, northern_dx) %>% as.data.frame()
northern_dx_df <- northern_dx_df[match(core_dx_df$SPID, northern_dx_df$SPID),] %>%
  dplyr::mutate(ordering=row_number(),
                geo="north") %>%
  dplyr::left_join(sp_traits, by="SPID") %>%
  dplyr::mutate(code=stringr::str_to_upper(paste0(substr(species, 1, 3),
                                                  substr(word(species, 2), 1, 3)))) %>%
  dplyr::mutate(rangePrecipBin=ifelse(rangePrecip <= quantile(.$rangePrecip, probs=0.25),
                                    "Driest",
                                    ifelse(rangePrecip >= quantile(.$rangePrecip, probs=0.75),
                                           "Wettest", "Average"))) %>%
  dplyr::mutate(rangePrecipBin=factor(rangePrecipBin, levels=c("Driest", "Average", "Wettest")))
write.csv(northern_dx_df, "../../output/trends/northern_dx_df_precip100.csv")

southern_dx_df <- do.call(rbind, southern_dx) %>% as.data.frame()
southern_dx_df <- southern_dx_df[match(core_dx_df$SPID, southern_dx_df$SPID),] %>%
  dplyr::mutate(ordering=row_number(),
                geo="south") %>%
  dplyr::left_join(sp_traits, by="SPID") %>%
  dplyr::mutate(code=stringr::str_to_upper(paste0(substr(species, 1, 3),
                                                  substr(word(species, 2), 1, 3)))) %>%
  dplyr::mutate(rangePrecipBin=ifelse(rangePrecip <= quantile(.$rangePrecip, probs=0.25),
                                    "Driest",
                                    ifelse(rangePrecip >= quantile(.$rangePrecip, probs=0.75),
                                           "Wettest", "Average"))) %>%
  dplyr::mutate(rangePrecipBin=factor(rangePrecipBin, levels=c("Driest", "Average", "Wettest")))
write.csv(southern_dx_df, "../../output/trends/southern_dx_df_precip100.csv")

total_dx_df <- do.call(rbind, total_dx) %>% as.data.frame()
total_dx_df <- total_dx_df[match(core_dx_df$SPID, total_dx_df$SPID),] %>%
  dplyr::mutate(ordering=row_number(),
                geo="south") %>%
  dplyr::left_join(sp_traits, by="SPID") %>%
  dplyr::mutate(code=stringr::str_to_upper(paste0(substr(species, 1, 3),
                                                  substr(word(species, 2), 1, 3)))) %>%
  dplyr::mutate(rangePrecipBin=ifelse(rangePrecip <= quantile(.$rangePrecip, probs=0.25),
                                    "Driest",
                                    ifelse(rangePrecip >= quantile(.$rangePrecip, probs=0.75),
                                           "Wettest", "Average"))) %>%
  dplyr::mutate(rangePrecipBin=factor(rangePrecipBin, levels=c("Driest", "Average", "Wettest")))
write.csv(total_dx_df, "../../output/trends/total_dx_df_precip100.csv")
beepr::beep()

master_dx_df <- do.call(rbind, list(core_dx_df, northern_dx_df, southern_dx_df))

FIGURE_THREE_A <- ggplot()+
  geom_hline(yintercept=0, linetype=2)+
  geom_linerange(NULL,
                 mapping=aes(x=seq(1,90,2), 
                             ymin=-Inf, 
                             ymax=northern_dx_df$mean[seq(1,90,2)]),
                 alpha=0.05,
                 position=position_nudge(x=0.25))+
  geom_linerange(NULL,
                 mapping=aes(x=seq(1,90,2), 
                             ymin=-Inf, 
                             ymax=core_dx_df$mean[seq(1,90,2)]),
                 alpha=0.05)+
  geom_linerange(NULL,
                 mapping=aes(x=seq(1,90,2), 
                             ymin=-Inf, 
                             ymax=southern_dx_df$mean[seq(1,90,2)]),
                 alpha=0.05,
                 position=position_nudge(x=-0.25))+
  # geom_curve(arrows_df,
  #            mapping=aes(x=x1, y=y1, xend=x2, yend=y2, group=id),
  #            arrow=arrow(length=unit(0.02, "npc")))+
  # geom_text(NULL,
  #           mapping=aes(x=15.1, y=0.05), size=3, hjust=0, fontface="italic",
  #           label="Although sharply declining in the southern and\ncore parts of their modeled range, some species\nexhibit positive/stable northern trends.")+
  geom_pointrange(master_dx_df,
                  mapping=aes(x=ordering, 
                              y=mean, 
                              ymax=mean+sd, 
                              ymin=mean-sd,
                              fill=rangePrecipBin,
                              color=rangePrecipBin,
                              shape=as.factor(geo)),
                  alpha=0)+
  geom_pointrange(northern_dx_df,
                  mapping=aes(x=ordering, 
                              y=mean, 
                              ymax=mean+sd, 
                              ymin=mean-sd,
                              fill=rangePrecipBin,
                              color=rangePrecipBin),
                  shape=24, size=0.1,
                  position=position_nudge(x=0.25))+
  geom_pointrange(core_dx_df,
                  mapping=aes(x=ordering, 
                              y=mean, 
                              ymax=mean+sd, 
                              ymin=mean-sd,
                              fill=rangePrecipBin,
                              color=rangePrecipBin),
                  shape=21, size=0.1)+
  geom_pointrange(southern_dx_df,
                  mapping=aes(x=ordering, 
                              y=mean, 
                              ymax=mean+sd, 
                              ymin=mean-sd,
                              fill=rangePrecipBin,
                              color=rangePrecipBin),
                  shape=25, size=0.1,
                  position=position_nudge(x=-0.25))+
  scale_shape_manual(values=c(21, 24, 25), name="Geographic Context", 
                     labels=c("Core", "Northern", "Southern"),
                     guide=guide_legend(title.position="top",
                                        override.aes=list(fill="black", alpha=1)))+
  scale_color_manual(values=c("#644119", "black", "#195464"),
                     name="Range Classification")+
  scale_fill_manual(values=c("#644119", "black", "#195464"),
                     name="Range Classification")+
  scale_x_continuous(breaks=seq(1, 90, 1), 
                     labels=core_dx_df$code,
                     expand=expansion(add=1),
                     name="Species Code")+
  scale_y_continuous(labels=scales::label_percent(add_plusses=TRUE),
                     name="Average Occupancy Gain/Loss\nfrom the 1970s to 2010s")+
  theme_cowplot()+
  theme(legend.position=c(0.05,0.9),
        legend.direction="horizontal",
        legend.box="horizontal",
        plot.background=element_rect(color=NULL, fill="white"),
        axis.title=element_text(size=18),
        axis.text=element_text(size=16),
        axis.text.x=element_text(size=10, angle=90, hjust=0, vjust=0.4))

sum_test_df <- left_join(southern_dx_df, northern_dx_df, by="SPID") %>%
  dplyr::mutate(southOverNorth=ifelse(mean.x > mean.y, 1, 0),
                northOverSouth=ifelse(mean.y > mean.x, 1, 0))

# FIGURE THREE COMPOSITE
FIGURE_THREE <- cowplot::ggdraw()+
  draw_plot(FIGURE_THREE_A)+
  draw_plot(response_100+theme(axis.title=element_text(size=7),
                               axis.text=element_text(size=7),
                               plot.background=element_rect(fill=NA)),
            x=0.7525, y=0.172,
            width=0.2, height=0.3)
ggsave2("../../figures/supplemental/FIGURE_2_precip100.png", FIGURE_THREE, dpi=400, height=6, width=12)

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
                rangePrecip_z = scale(rangePrecip),
                varTemp_z = scale(rangeTempVar),
                varPrecip_z = scale(rangePrecipVar)) %>%
  dplyr::select(species, mean, se, sd, rangeTemp_z, 
                rangePrecip_z, varTemp_z, varPrecip_z,
                rangeSize_z, aveWingspan_z, 
                numReportedHostplantFamilies_z,
                diapauseStage_z, disturbanceAffinity_z) %>%
  dplyr::filter(complete.cases(.))

sp_tree <- ape::read.tree("../../data/taxa/SupDryad_treepl.tre")
sp_tree$tip.label <- str_replace_all(sp_tree$tip.label, "_", " ")
sp_tree$tip.label <- stringr::word(sp_tree$tip.label, start=3, end=4)

sp_tree <- sp_tree %>%
  ape::keep.tip(total_dx_df_model$species) %>%
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
saveRDS(MODEL_A, "../../output/modelFiles/100kmPrecip_ModelA.rds")

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
saveRDS(MODEL_B, "../../output/modelFiles/100kmPrecip_ModelB.rds")

# TEMPERATURE MODEL (MODEL C)
MODEL_C <- brms::brm(mean~1+
                       rangePrecip_z,
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
saveRDS(MODEL_C, "../../output/modelFiles/100kmPrecip_ModelC.rds")

# TEMPERATURE + PHYLOGENY MODEL (MODEL D)
MODEL_D <- brms::brm(mean~1+
                       rangePrecip_z+
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
saveRDS(MODEL_D, "../../output/modelFiles/100kmPrecip_ModelD.rds")

# TEMPERATURE MODEL (MODEL C2)
MODEL_C2 <- brms::brm(mean~1+
                        varPrecip_z,
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
saveRDS(MODEL_C2, "../../output/modelFiles/100kmPrecip_ModelC2.rds")

# TEMPERATURE + PHYLOGENY MODEL (MODEL D2)
MODEL_D2 <- brms::brm(mean~1+
                        varPrecip_z+
                        (1|gr(species, cov=sp_tree)),
                      data=total_dx_df_model,
                      data2=list(sp_tree=sp_tree),
                      iter=200000,
                      warmup=100000,
                      thin=50,
                      family=gaussian(),
                      prior=c(prior(normal(0,10), "Intercept")),
                      control=list(max_treedepth=16,
                                   adapt_delta=0.99999),
                      cores=5,
                      save_pars = save_pars(all = TRUE))
saveRDS(MODEL_D2, "../../output/modelFiles/100kmPrecip_ModelD2.rds")

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
saveRDS(MODEL_E, "../../output/modelFiles/100kmPrecip_ModelE.rds")

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
saveRDS(MODEL_F, "../../output/modelFiles/100kmPrecip_ModelF.rds")

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
saveRDS(MODEL_G, "../../output/modelFiles/100kmPrecip_ModelG.rds")

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
saveRDS(MODEL_H, "../../output/modelFiles/100kmPrecip_ModelH.rds")

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
saveRDS(MODEL_I, "../../output/modelFiles/100kmPrecip_ModelI.rds")

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
saveRDS(MODEL_J, "../../output/modelFiles/100kmPrecip_ModelJ.rds")

# OVERWINTERING MODEL (MODEL K)
MODEL_K <- brms::brm(mean~1+
                       rangePrecip_z+
                       diapauseStage_z+
                       rangePrecip_z:diapauseStage_z,
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
saveRDS(MODEL_K, "../../output/modelFiles/100kmPrecip_ModelK.rds")

# OVERWINTERING + PHYLOGENY MODEL (MODEL L)
MODEL_L <- brms::brm(mean~1+
                       rangePrecip_z+
                       diapauseStage_z+
                       rangePrecip_z:diapauseStage_z+
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
saveRDS(MODEL_L, "../../output/modelFiles/100kmPrecip_ModelL.rds")

# COMPLEX MODEL (MODEL M)
MODEL_M <- brms::brm(mean~1+
                       rangePrecip_z+
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
saveRDS(MODEL_M, "../../output/modelFiles/100kmPrecip_ModelM.rds")

# COMPLEX + PHYLOGENY MODEL (MODEL N)
MODEL_N <- brms::brm(mean~1+
                       rangePrecip_z+
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
saveRDS(MODEL_N, "../../output/modelFiles/100kmPrecip_ModelN.rds")


####################################################################################################
## COMPARE ALL OF THE MODELS TO ONE ANOTHER AND VISUALIZE
####################################################################################################
MODEL_A <- readRDS("../../output/modelFiles/100kmPrecip_ModelA.rds")
MODEL_B <- readRDS("../../output/modelFiles/100kmPrecip_ModelB.rds")
MODEL_C <- readRDS("../../output/modelFiles/100kmPrecip_ModelC.rds")
MODEL_D <- readRDS("../../output/modelFiles/100kmPrecip_ModelD.rds")
MODEL_E <- readRDS("../../output/modelFiles/100kmPrecip_ModelE.rds")
MODEL_F <- readRDS("../../output/modelFiles/100kmPrecip_ModelF.rds")
MODEL_G <- readRDS("../../output/modelFiles/100kmPrecip_ModelG.rds")
MODEL_H <- readRDS("../../output/modelFiles/100kmPrecip_ModelH.rds")
MODEL_I <- readRDS("../../output/modelFiles/100kmPrecip_ModelI.rds")
MODEL_J <- readRDS("../../output/modelFiles/100kmPrecip_ModelJ.rds")
MODEL_K <- readRDS("../../output/modelFiles/100kmPrecip_ModelK.rds")
MODEL_L <- readRDS("../../output/modelFiles/100kmPrecip_ModelL.rds")
MODEL_M <- readRDS("../../output/modelFiles/100kmPrecip_ModelM.rds")
MODEL_N <- readRDS("../../output/modelFiles/100kmPrecip_ModelN.rds")

loo::loo_compare(loo::loo(MODEL_A, moment_match=TRUE),
                 loo::loo(MODEL_B, moment_match=TRUE),
                 loo::loo(MODEL_C, moment_match=TRUE),
                 loo::loo(MODEL_D, moment_match=TRUE),
                 loo::loo(MODEL_C2, moment_match=TRUE),
                 loo::loo(MODEL_D2, moment_match=TRUE),
                 loo::loo(MODEL_E, moment_match=TRUE),
                 loo::loo(MODEL_F, moment_match=TRUE),
                 loo::loo(MODEL_G, moment_match=TRUE),
                 loo::loo(MODEL_H, moment_match=TRUE),
                 loo::loo(MODEL_I, moment_match=TRUE),
                 loo::loo(MODEL_J, moment_match=TRUE),
                 loo::loo(MODEL_K, moment_match=TRUE),
                 loo::loo(MODEL_L, moment_match=TRUE),
                 loo::loo(MODEL_M, moment_match=TRUE),
                 loo::loo(MODEL_N, moment_match=TRUE))

# GRAB DRAWS
MODEL_C_DRAWS <- tidybayes::gather_draws(MODEL_C,
                                         b_Intercept, 
                                         b_rangePrecip_z, regex=TRUE)  %>%
  dplyr::mutate(.variable=str_replace(.variable, "b_", "")) %>%
  dplyr::mutate(.variable=factor(.variable, levels=c("rangeTemp_z",
                                                     "Intercept")))

# FIGURE FOUR ########################################
MODEL_C_LINES <- total_dx_df_model %>%
  modelr::data_grid(rangePrecip_z=seq(-3, 2, length.out=50)) %>%
  tidybayes::add_predicted_draws(MODEL_C)


ggplot(MODEL_C_LINES,
                       mapping=aes(x=rangePrecip_z, y=mean))+
  tidybayes::stat_lineribbon(mapping=aes(y=.prediction),
                             .width=c(0.95, 0.8, 0.5))+
  scale_fill_brewer(palette="Greens", name="Posterior Predictive Interval")+
  geom_point(data=total_dx_df_model,
             mapping=aes(x=rangePrecip_z, y=mean),
             size=2)+
  geom_hline(yintercept=0, linetype=2)+
  scale_y_continuous(labels=scales::percent_format(),
                     name="Mean Occupancy Shift from\nthe 1970s to 2010s",
                     limits=c(-0.12, 0.12))+
  scale_x_continuous(name="Range-wide Average Annual Precipitation [mm]",
                     limits=c(-2,2),
                     breaks=seq(-2,2,0.5),
                     labels=round((seq(-2,2,0.5)*sd(sp_traits$rangePrecip))+
                                    mean(sp_traits$rangePrecip)), 1)+
  theme_cowplot()+
  theme(plot.background=element_rect(fill="white",
                                     color="white"),
        legend.position=c(0.1, 0.85))
ggsave2("../../figures/supplemental/FIGURE_004_precip100.png", dpi=400, height=5, width=10)

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