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
MCMCvis::MCMCtrace(my_res_100_1, Rhat=TRUE, filename="../../figures/supplemental/SupplementalFile_S3.pdf")
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
  my_occ_dx_map <- ggplot()+
    geom_sf(basemap, mapping=aes(), fill=NA)+
    geom_tile(my_rast_sf, mapping=aes(x=x, y=y, fill=layer),
              alpha=0.7)+
    colorspace::scale_fill_continuous_divergingx(palette="Zissou 1", rev=TRUE,
                                                 na.value=NA, labels=scales::percent,
                                                 name="Shift in 50-year\nOccupancy Probability")+
    ggtitle(sp_traits[i,]$species)+
    theme_map()+
    theme(legend.position=c(0.05, 0.3),
          plot.title=element_text(face="italic"))
  ggsave2(paste0("../../figures/supplemental/rangeMaps/", sp_traits[i,]$binomial, "_100precip.pdf"),
          dpi=400, height=6, width=8)
}
core_dx_df <- do.call(rbind, core_dx) %>% as.data.frame() %>%
  dplyr::arrange(mean) %>%
  dplyr::mutate(ordering=row_number(),
                geo="core") %>%
  dplyr::left_join(sp_traits, by="SPID") %>%
  dplyr::mutate(code=stringr::str_to_upper(paste0(substr(species, 1, 3),
                                                  substr(word(species, 2), 1, 3))))
write.csv(core_dx_df, "../../output/trends/core_dx_df_precip100.csv")

northern_dx_df <- do.call(rbind, northern_dx) %>% as.data.frame()
northern_dx_df <- northern_dx_df[match(core_dx_df$SPID, northern_dx_df$SPID),] %>%
  dplyr::mutate(ordering=row_number(),
                geo="north") %>%
  dplyr::left_join(sp_traits, by="SPID") %>%
  dplyr::mutate(code=stringr::str_to_upper(paste0(substr(species, 1, 3),
                                                  substr(word(species, 2), 1, 3))))
write.csv(northern_dx_df, "../../output/trends/northern_dx_df_precip100.csv")

southern_dx_df <- do.call(rbind, southern_dx) %>% as.data.frame()
southern_dx_df <- southern_dx_df[match(core_dx_df$SPID, southern_dx_df$SPID),] %>%
  dplyr::mutate(ordering=row_number(),
                geo="south") %>%
  dplyr::left_join(sp_traits, by="SPID") %>%
  dplyr::mutate(code=stringr::str_to_upper(paste0(substr(species, 1, 3),
                                                  substr(word(species, 2), 1, 3))))
write.csv(southern_dx_df, "../../output/trends/southern_dx_df_precip100.csv")

total_dx_df <- do.call(rbind, total_dx) %>% as.data.frame()
total_dx_df <- total_dx_df[match(core_dx_df$SPID, total_dx_df$SPID),] %>%
  dplyr::mutate(ordering=row_number(),
                geo="south") %>%
  dplyr::left_join(sp_traits, by="SPID") %>%
  dplyr::mutate(code=stringr::str_to_upper(paste0(substr(species, 1, 3),
                                                  substr(word(species, 2), 1, 3))))
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
                              fill=rangePrecip,
                              color=rangePrecip,
                              shape=as.factor(geo)),
                  alpha=0)+
  geom_pointrange(northern_dx_df,
                  mapping=aes(x=ordering, 
                              y=mean, 
                              ymax=mean+sd, 
                              ymin=mean-sd,
                              fill=rangePrecip,
                              color=rangePrecip),
                  shape=24, size=0.1,
                  position=position_nudge(x=0.25))+
  geom_pointrange(core_dx_df,
                  mapping=aes(x=ordering, 
                              y=mean, 
                              ymax=mean+sd, 
                              ymin=mean-sd,
                              fill=rangePrecip,
                              color=rangePrecip),
                  shape=21, size=0.1)+
  geom_pointrange(southern_dx_df,
                  mapping=aes(x=ordering, 
                              y=mean, 
                              ymax=mean+sd, 
                              ymin=mean-sd,
                              fill=rangePrecip,
                              color=rangePrecip),
                  shape=25, size=0.1,
                  position=position_nudge(x=-0.25))+
  scale_shape_manual(values=c(21, 24, 25), name="Geographic Context", 
                     labels=c("Core", "Northern", "Southern"),
                     guide=guide_legend(title.position="top",
                                        override.aes=list(fill="black", alpha=1)))+
  scale_fill_gradientn(colors=c("#644119", "#a16928", "grey45", "#2887a1", "#195464"),
                       name="Ave. Range Precip. [mm]",
                       guide=guide_colorbar(title.position="top",
                                            barwidth=unit(5, "cm")))+
  scale_color_gradientn(colors=c("#644119", "#a16928", "grey45", "#2887a1", "#195464"),
                        name="Ave. Range Precip. [mm]",
                        guide=guide_colorbar(title.position="top",
                                             barwidth=unit(5, "cm")))+
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
core_dx_df <- master_dx_df %>%
  dplyr::filter(geo=="core") %>%
  dplyr::mutate(rangeSize_z = scale(rangeSize),
                rangePrecip_z = scale(rangePrecip),
                aveWingspan_z = scale(aveWingspan),
                numReportedHostplantFamilies_z = scale(numReportedHostplantFamilies),
                diapauseStage_z = factor(diapauseStage, levels=c("Larva", "Egg", "Pupa", "Adult")),
                disturbanceAffinity = ifelse(disturbanceAffinity=="Mixed", "Generalist", disturbanceAffinity)) %>%
  dplyr::mutate(disturbanceAffinity_z = factor(disturbanceAffinity, levels=c("Generalist", "Avoidant", "Associated")),
                edgeAffinity_z = factor(edgeAffinity, levels=c("Geenralist", "Avoidant", "Associated")),
                canopyAffinity_z = factor(canopyAffinity, levels=c("Generalist", "Mixed", "Open")),
                voltinism_z = factor(voltinism, levels=c("Univoltine", "Multivoltine"))) %>%
  dplyr::select(species, mean, se, rangePrecip_z, rangeSize_z, aveWingspan_z, numReportedHostplantFamilies_z,
                diapauseStage_z, disturbanceAffinity_z) %>%
  dplyr::filter(complete.cases(.))

northern_dx_df <- master_dx_df %>%
  dplyr::filter(geo=="north") %>%
  dplyr::mutate(rangeSize_z = scale(rangeSize),
                rangePrecip_z = scale(rangePrecip),
                aveWingspan_z = scale(aveWingspan),
                numReportedHostplantFamilies_z = scale(numReportedHostplantFamilies),
                diapauseStage_z = factor(diapauseStage, levels=c("Larva", "Egg", "Pupa", "Adult")),
                disturbanceAffinity = ifelse(disturbanceAffinity=="Mixed", "Generalist", disturbanceAffinity)) %>%
  dplyr::mutate(disturbanceAffinity_z = factor(disturbanceAffinity, levels=c("Generalist", "Avoidant", "Associated")),
                edgeAffinity_z = factor(edgeAffinity, levels=c("Geenralist", "Avoidant", "Associated")),
                canopyAffinity_z = factor(canopyAffinity, levels=c("Generalist", "Mixed", "Open")),
                voltinism_z = factor(voltinism, levels=c("Univoltine", "Multivoltine"))) %>%
  dplyr::select(species, mean, se, rangePrecip_z, rangeSize_z, aveWingspan_z, numReportedHostplantFamilies_z,
                diapauseStage_z, disturbanceAffinity_z) %>%
  dplyr::filter(complete.cases(.))

southern_dx_df <- master_dx_df %>%
  dplyr::filter(geo=="south") %>%
  dplyr::mutate(rangeSize_z = scale(rangeSize),
                rangePrecip_z = scale(rangePrecip),
                aveWingspan_z = scale(aveWingspan),
                numReportedHostplantFamilies_z = scale(numReportedHostplantFamilies),
                diapauseStage_z = factor(diapauseStage, levels=c("Larva", "Egg", "Pupa", "Adult")),
                disturbanceAffinity = ifelse(disturbanceAffinity=="Mixed", "Generalist", disturbanceAffinity)) %>%
  dplyr::mutate(disturbanceAffinity_z = factor(disturbanceAffinity, levels=c("Generalist", "Avoidant", "Associated")),
                edgeAffinity_z = factor(edgeAffinity, levels=c("Geenralist", "Avoidant", "Associated")),
                canopyAffinity_z = factor(canopyAffinity, levels=c("Generalist", "Mixed", "Open")),
                voltinism_z = factor(voltinism, levels=c("Univoltine", "Multivoltine"))) %>%
  dplyr::select(species, mean, se, rangePrecip_z, rangeSize_z, aveWingspan_z, numReportedHostplantFamilies_z,
                diapauseStage_z, disturbanceAffinity_z) %>%
  dplyr::filter(complete.cases(.))

sp_tree <- ape::read.tree("../../data/taxa/SupDryad_treepl.tre")
sp_tree$tip.label <- str_replace_all(sp_tree$tip.label, "_", " ")
sp_tree$tip.label <- stringr::word(sp_tree$tip.label, start=3, end=4)

sp_tree <- sp_tree %>%
  ape::keep.tip(core_dx_df$species) %>%
  ape::vcv.phylo()

# INTERCEPT ONLY MODELS
library(brms)
my_fit_1_core <- brms::brm(mean~1,
                           data=core_dx_df,
                           iter=200000,
                           warmup=100000,
                           thin=50,
                           family=gaussian(),
                           prior=c(prior(normal(0,10), "Intercept")),
                           control=list(max_treedepth=15,
                                        adapt_delta=0.99),
                           cores=5,
                           save_pars = save_pars(all = TRUE))
my_fit_1_north <- brms::brm(mean~1,
                            data=northern_dx_df,
                            iter=200000,
                            warmup=100000,
                            thin=50,
                            family=gaussian(),
                            prior=c(prior(normal(0,10), "Intercept")),
                            control=list(max_treedepth=15,
                                         adapt_delta=0.99),
                            cores=5,
                            save_pars = save_pars(all = TRUE))
my_fit_1_south <- brms::brm(mean~1,
                            data=southern_dx_df,
                            iter=200000,
                            warmup=100000,
                            thin=50,
                            family=gaussian(),
                            prior=c(prior(normal(0,10), "Intercept")),
                            control=list(max_treedepth=15,
                                         adapt_delta=0.99),
                            cores=5,
                            save_pars = save_pars(all = TRUE))
my_fit_1_list <- list(my_fit_1_core, my_fit_1_south, my_fit_1_north)
saveRDS(my_fit_1_list, "../../output/posthoc_100km_precip_fit1.rds")

# PHYLOGENETIC INTERCEPT MODELS
my_fit_2_core <- brms::brm(mean~1+
                             (1|gr(species, cov=sp_tree)),
                           data=core_dx_df,
                           data2=list(sp_tree=sp_tree),
                           iter=200000,
                           warmup=100000,
                           thin=50,
                           family=gaussian(),
                           prior=c(prior(normal(0,10), class="Intercept")),
                           control=list(max_treedepth=15,
                                        adapt_delta=0.99),
                           cores=5,
                           save_pars = save_pars(all = TRUE))
my_fit_2_north <- brms::brm(mean~1+
                              (1|gr(species, cov=sp_tree)),
                            data=northern_dx_df,
                            data2=list(sp_tree=sp_tree),
                            iter=200000,
                            warmup=100000,
                            thin=50,
                            family=gaussian(),
                            prior=c(prior(normal(0,10), "Intercept")),
                            control=list(max_treedepth=15,
                                         adapt_delta=0.99),
                            cores=5,
                            save_pars = save_pars(all = TRUE))
my_fit_2_south <- brms::brm(mean~1+
                              (1|gr(species, cov=sp_tree)),
                            data=southern_dx_df,
                            data2=list(sp_tree=sp_tree),
                            iter=200000,
                            warmup=100000,
                            thin=50,
                            family=gaussian(),
                            prior=c(prior(normal(0,10), "Intercept")),
                            control=list(max_treedepth=15,
                                         adapt_delta=0.99),
                            cores=5,
                            save_pars = save_pars(all = TRUE))
my_fit_2_list <- list(my_fit_2_core, my_fit_2_south, my_fit_2_north)
saveRDS(my_fit_2_list, "../../output/posthoc_100km_precip_fit2.rds")

# TEMPERATURE ONLY MODELS
my_fit_3_core <- brms::brm(mean~1+
                             rangePrecip_z,
                           data=core_dx_df,
                           iter=200000,
                           warmup=100000,
                           thin=50,
                           family=gaussian(),
                           prior=c(prior(normal(0,10), "b"),
                                   prior(normal(0,10), "Intercept")),
                           control=list(max_treedepth=15,
                                        adapt_delta=0.99),
                           cores=5,
                           save_pars = save_pars(all = TRUE))
my_fit_3_north <- brms::brm(mean~1+
                              rangePrecip_z,
                            data=northern_dx_df,
                            iter=200000,
                            warmup=100000,
                            thin=50,
                            family=gaussian(),
                            prior=c(prior(normal(0,10), "b"),
                                    prior(normal(0,10), "Intercept")),
                            control=list(max_treedepth=15,
                                         adapt_delta=0.99),
                            cores=5,
                            save_pars = save_pars(all = TRUE))
my_fit_3_south <- brms::brm(mean~1+
                              rangePrecip_z,
                            data=southern_dx_df,
                            iter=200000,
                            warmup=100000,
                            thin=50,
                            family=gaussian(),
                            prior=c(prior(normal(0,10), "b"),
                                    prior(normal(0,10), "Intercept")),
                            control=list(max_treedepth=15,
                                         adapt_delta=0.99),
                            cores=5,
                            save_pars = save_pars(all = TRUE))
my_fit_3_list <- list(my_fit_3_core, my_fit_3_south, my_fit_3_north)
saveRDS(my_fit_3_list, "../../output/posthoc_100km_precip_fit3.rds")

# TEMPERATURE AND PHYLOGENY MODEL
my_fit_4_core <- brms::brm(mean~1+
                             rangePrecip_z+
                             (1|gr(species, cov=sp_tree)),
                           data=core_dx_df,
                           data2=list(sp_tree=sp_tree),
                           iter=200000,
                           warmup=100000,
                           thin=50,
                           family=gaussian(),
                           prior=c(prior(normal(0,10), "b"),
                                   prior(normal(0,10), "Intercept")),
                           control=list(max_treedepth=15,
                                        adapt_delta=0.9999),
                           cores=5,
                           save_pars = save_pars(all = TRUE))
my_fit_4_north <- brms::brm(mean~1+
                              rangePrecip_z+
                              (1|gr(species, cov=sp_tree)),
                            data=northern_dx_df,
                            data2=list(sp_tree=sp_tree),
                            iter=200000,
                            warmup=100000,
                            thin=50,
                            family=gaussian(),
                            prior=c(prior(normal(0,10), "b"),
                                    prior(normal(0,10), "Intercept")),
                            control=list(max_treedepth=15,
                                         adapt_delta=0.9999),
                            cores=5,
                            save_pars = save_pars(all = TRUE))
my_fit_4_south <- brms::brm(mean~1+
                              rangePrecip_z+
                              (1|gr(species, cov=sp_tree)),
                            data=southern_dx_df,
                            data2=list(sp_tree=sp_tree),
                            iter=200000,
                            warmup=100000,
                            thin=50,
                            family=gaussian(),
                            prior=c(prior(normal(0,10), "b"),
                                    prior(normal(0,10), "Intercept")),
                            control=list(max_treedepth=15,
                                         adapt_delta=0.999),
                            cores=5,
                            save_pars = save_pars(all = TRUE))
my_fit_4_list <- list(my_fit_4_core, my_fit_4_south, my_fit_4_north)
saveRDS(my_fit_4_list, "../../output/posthoc_100km_precip_fit4.rds")

# TEMPERATURE AND OVERWINTERING STAGE INTERACTION
my_fit_5_core <- brms::brm(mean~1+
                             rangePrecip_z*diapauseStage_z,
                           data=core_dx_df,
                           iter=200000,
                           warmup=100000,
                           thin=50,
                           family=gaussian(),
                           prior=c(prior(normal(0,10), "b"),
                                   prior(normal(0,10), "Intercept")),
                           control=list(max_treedepth=15,
                                        adapt_delta=0.99),
                           cores=5,
                           save_pars = save_pars(all = TRUE))
my_fit_5_north <- brms::brm(mean~1+
                              rangePrecip_z*diapauseStage_z,
                            data=northern_dx_df,
                            iter=200000,
                            warmup=100000,
                            thin=50,
                            family=gaussian(),
                            prior=c(prior(normal(0,10), "b"),
                                    prior(normal(0,10), "Intercept")),
                            control=list(max_treedepth=15,
                                         adapt_delta=0.99),
                            cores=5,
                            save_pars = save_pars(all = TRUE))
my_fit_5_south <- brms::brm(mean~1+
                              rangePrecip_z*diapauseStage_z,
                            data=southern_dx_df,
                            iter=200000,
                            warmup=100000,
                            thin=50,
                            family=gaussian(),
                            prior=c(prior(normal(0,10), "b"),
                                    prior(normal(0,10), "Intercept")),
                            control=list(max_treedepth=15,
                                         adapt_delta=0.99),
                            cores=5,
                            save_pars = save_pars(all = TRUE))
my_fit_5_list <- list(my_fit_5_core, my_fit_5_south, my_fit_5_north)
saveRDS(my_fit_5_list, "../../output/posthoc_100km_precip_fit5.rds")

# TEMPERATURE AND OVERWINTERING STAGE INTERACTION WITH PHYLOGENY
my_fit_6_core <- brms::brm(mean~1+
                             rangePrecip_z*diapauseStage_z+
                             (1|gr(species, cov=sp_tree)),
                           data=core_dx_df,
                           data2=list(sp_tree=sp_tree),
                           iter=200000,
                           warmup=100000,
                           thin=50,
                           family=gaussian(),
                           prior=c(prior(normal(0,10), "b"),
                                   prior(normal(0,10), "Intercept")),
                           control=list(max_treedepth=15,
                                        adapt_delta=0.99),
                           cores=5,
                           save_pars = save_pars(all = TRUE))
my_fit_6_north <- brms::brm(mean~1+
                              rangePrecip_z*diapauseStage_z+
                              (1|gr(species, cov=sp_tree)),
                            data=northern_dx_df,
                            data2=list(sp_tree=sp_tree),
                            iter=200000,
                            warmup=100000,
                            thin=50,
                            family=gaussian(),
                            prior=c(prior(normal(0,10), "b"),
                                    prior(normal(0,10), "Intercept")),
                            control=list(max_treedepth=15,
                                         adapt_delta=0.99),
                            cores=5,
                            save_pars = save_pars(all = TRUE))
my_fit_6_south <- brms::brm(mean~1+
                              rangePrecip_z*diapauseStage_z+
                              (1|gr(species, cov=sp_tree)),
                            data=southern_dx_df,
                            data2=list(sp_tree=sp_tree),
                            iter=200000,
                            warmup=100000,
                            thin=50,
                            family=gaussian(),
                            prior=c(prior(normal(0,10), "b"),
                                    prior(normal(0,10), "Intercept")),
                            control=list(max_treedepth=15,
                                         adapt_delta=0.99),
                            cores=5,
                            save_pars = save_pars(all = TRUE))
my_fit_6_list <- list(my_fit_6_core, my_fit_6_south, my_fit_6_north)
saveRDS(my_fit_6_list, "../../output/posthoc_100km_precip_fit6.rds")

# TRAIT ONLY MODELS (ALL TRAITS NO INTERACTIONS)
my_fit_7_core <- brms::brm(mean~1+
                             rangePrecip_z+
                             rangeSize_z+
                             aveWingspan_z+
                             numReportedHostplantFamilies_z+
                             diapauseStage_z+
                             disturbanceAffinity_z,
                           data=core_dx_df,
                           iter=200000,
                           warmup=100000,
                           thin=50,
                           family=gaussian(),
                           prior=c(prior(normal(0,10), "b"),
                                   prior(normal(0,10), "Intercept")),
                           control=list(max_treedepth=15,
                                        adapt_delta=0.99),
                           cores=5,
                           save_pars = save_pars(all = TRUE))
my_fit_7_north <- brms::brm(mean~1+
                              rangePrecip_z+
                              rangeSize_z+
                              aveWingspan_z+
                              numReportedHostplantFamilies_z+
                              diapauseStage_z+
                              disturbanceAffinity_z,
                            data=northern_dx_df,
                            iter=200000,
                            warmup=100000,
                            thin=50,
                            family=gaussian(),
                            prior=c(prior(normal(0,10), "b"),
                                    prior(normal(0,10), "Intercept")),
                            control=list(max_treedepth=15,
                                         adapt_delta=0.99),
                            cores=5,
                            save_pars = save_pars(all = TRUE))
my_fit_7_south <- brms::brm(mean~1+
                              rangePrecip_z+
                              rangeSize_z+
                              aveWingspan_z+
                              numReportedHostplantFamilies_z+
                              diapauseStage_z+
                              disturbanceAffinity_z,
                            data=southern_dx_df,
                            iter=200000,
                            warmup=100000,
                            thin=50,
                            family=gaussian(),
                            prior=c(prior(normal(0,10), "b"),
                                    prior(normal(0,10), "Intercept")),
                            control=list(max_treedepth=15,
                                         adapt_delta=0.99),
                            cores=5,
                            save_pars = save_pars(all = TRUE))
my_fit_7_list <- list(my_fit_7_core, my_fit_7_south, my_fit_7_north)
saveRDS(my_fit_7_list, "../../output/posthoc_100km_precip_fit7.rds")

# TRAIT AND PHYLOGENY MODELS
my_fit_8_core <- brms::brm(mean~1+
                             rangePrecip_z+
                             rangeSize_z+
                             aveWingspan_z+
                             numReportedHostplantFamilies_z+
                             diapauseStage_z+
                             disturbanceAffinity_z+
                             (1|gr(species, cov=sp_tree)),
                           data=core_dx_df,
                           data2=list(sp_tree=sp_tree),
                           iter=200000,
                           warmup=100000,
                           thin=50,
                           family=gaussian(),
                           prior=c(prior(normal(0,10), "b"),
                                   prior(normal(0,10), "Intercept")),
                           control=list(max_treedepth=15,
                                        adapt_delta=0.9999),
                           cores=5,
                           save_pars = save_pars(all = TRUE))
my_fit_8_north <- brms::brm(mean~1+
                              rangePrecip_z+
                              rangeSize_z+
                              aveWingspan_z+
                              numReportedHostplantFamilies_z+
                              diapauseStage_z+
                              disturbanceAffinity_z+
                              (1|gr(species, cov=sp_tree)),
                            data=northern_dx_df,
                            data2=list(sp_tree=sp_tree),
                            iter=200000,
                            warmup=100000,
                            thin=50,
                            family=gaussian(),
                            prior=c(prior(normal(0,10), "b"),
                                    prior(normal(0,10), "Intercept")),
                            control=list(max_treedepth=15,
                                         adapt_delta=0.9999),
                            cores=5,
                            save_pars = save_pars(all = TRUE))
my_fit_8_south <- brms::brm(mean~1+
                              rangePrecip_z+
                              rangeSize_z+
                              aveWingspan_z+
                              numReportedHostplantFamilies_z+
                              diapauseStage_z+
                              disturbanceAffinity_z+
                              (1|gr(species, cov=sp_tree)),
                            data=southern_dx_df,
                            data2=list(sp_tree=sp_tree),
                            iter=200000,
                            warmup=100000,
                            thin=50,
                            family=gaussian(),
                            prior=c(prior(normal(0,10), "b"),
                                    prior(normal(0,10), "Intercept")),
                            control=list(max_treedepth=15,
                                         adapt_delta=0.999),
                            cores=5,
                            save_pars = save_pars(all = TRUE))
my_fit_8_list <- list(my_fit_8_core, my_fit_8_south, my_fit_8_north)
saveRDS(my_fit_8_list, "../../output/posthoc_100km_precip_fit8.rds")

####################################################################################################
## COMPARE ALL OF THE MODELS TO ONE ANOTHER AND VISUALIZE
####################################################################################################
my_fit_1_list <- readRDS("../../output/posthoc_100km_precip_fit1.rds")
my_fit_2_list <- readRDS("../../output/posthoc_100km_precip_fit2.rds")
my_fit_3_list <- readRDS("../../output/posthoc_100km_precip_fit3.rds")
my_fit_4_list <- readRDS("../../output/posthoc_100km_precip_fit4.rds")
my_fit_5_list <- readRDS("../../output/posthoc_100km_precip_fit5.rds")
my_fit_6_list <- readRDS("../../output/posthoc_100km_precip_fit6.rds")
my_fit_7_list <- readRDS("../../output/posthoc_100km_precip_fit7.rds")
my_fit_8_list <- readRDS("../../output/posthoc_100km_precip_fit8.rds")

loo::loo_compare(loo::loo(my_fit_1_core, moment_match=TRUE), 
                 loo::loo(my_fit_2_core, moment_match=TRUE), 
                 loo::loo(my_fit_3_core, moment_match=TRUE), 
                 loo::loo(my_fit_4_core, moment_match=TRUE),
                 loo::loo(my_fit_5_core, moment_match=TRUE),
                 loo::loo(my_fit_6_core, moment_match=TRUE),
                 loo::loo(my_fit_7_core, moment_match=TRUE),
                 loo::loo(my_fit_8_core, moment_match=TRUE))

loo::loo_compare(loo::loo(my_fit_1_north, moment_match=TRUE), 
                 loo::loo(my_fit_2_north, moment_match=TRUE), 
                 loo::loo(my_fit_3_north, moment_match=TRUE), 
                 loo::loo(my_fit_4_north, moment_match=TRUE),
                 loo::loo(my_fit_5_north, moment_match=TRUE),
                 loo::loo(my_fit_6_north, moment_match=TRUE),
                 loo::loo(my_fit_7_north, moment_match=TRUE),
                 loo::loo(my_fit_8_north, moment_match=TRUE))

loo::loo_compare(loo::loo(my_fit_1_south, moment_match=TRUE), 
                 loo::loo(my_fit_2_south, moment_match=TRUE), 
                 loo::loo(my_fit_3_south, moment_match=TRUE), 
                 loo::loo(my_fit_4_south, moment_match=TRUE),
                 loo::loo(my_fit_5_south, moment_match=TRUE),
                 loo::loo(my_fit_6_south, moment_match=TRUE),
                 loo::loo(my_fit_7_south, moment_match=TRUE),
                 loo::loo(my_fit_8_south, moment_match=TRUE))


# Grab the draws from the top candidate models
sp_traits <- sp_traits %>%
  dplyr::select(family, species) %>%
  arrange(family) %>%
  group_by(family) %>%
  arrange(species, .by_group=TRUE) %>%
  ungroup() %>%
  dplyr::mutate(globalOrdering=row_number())


# Pull the trait parameter estimates
my_draws_4_south <- tidybayes::gather_draws(my_fit_4_south,
                                            b_Intercept, b_rangePrecip_z, regex=TRUE) %>%
  dplyr::mutate(.variable=str_replace(.variable, "b_", "")) %>%
  dplyr::mutate(.variable=factor(.variable, levels=c("rangePrecip_z",
                                                     "Intercept")))

my_draws_4_core <- tidybayes::gather_draws(my_fit_4_core,
                                           b_Intercept, b_rangePrecip_z, regex=TRUE)  %>%
  dplyr::mutate(.variable=str_replace(.variable, "b_", "")) %>%
  dplyr::mutate(.variable=factor(.variable, levels=c("rangePrecip_z",
                                                     "Intercept")))

my_draws_4_north <- tidybayes::gather_draws(my_fit_4_north,
                                            b_Intercept, b_rangePrecip_z, regex=TRUE)  %>%
  dplyr::mutate(.variable=str_replace(.variable, "b_", "")) %>%
  dplyr::mutate(.variable=factor(.variable, levels=c("rangePrecip_z",
                                                     "Intercept")))

# FIGURE FOUR ########################################
ggplot()+
  tidybayes::stat_interval(my_draws_4_north,
                           mapping=aes(x=.value, y=.variable),
                           position=position_nudge(y=+0.2))+
  scale_color_brewer(palette="Purples", name="Cred. Int.")+
  ggnewscale::new_scale_color()+
  tidybayes::stat_interval(my_draws_4_core,
                           mapping=aes(x=.value, y=.variable))+
  scale_color_brewer(palette="Greys", name="Cred. Int.")+
  ggnewscale::new_scale_color()+
  tidybayes::stat_interval(my_draws_4_south,
                           mapping=aes(x=.value, y=.variable),
                           position=position_nudge(y=-0.2))+
  scale_color_brewer(palette="Purples", name="Cred. Int.")+
  geom_vline(xintercept=0, linetype=2)+
  labs(x="Parameter Estimate", y="Parameter")+
  scale_y_discrete(labels=c("Range-wide\nPrecip.", "Intercept"))+
  theme_cowplot()+
  theme(plot.background=element_rect(fill="white", color="white"),
        axis.text.y=element_text(angle=90, hjust=0.5),
        legend.position="top")
ggsave2("../../figures/supplemental/FIGURE_004_precip100.png", dpi=400, height=3, width=6)

hyp <- "sd_species__Intercept^2 / (sd_species__Intercept^2 + sigma^2) = 0"
(hyp <- brms::hypothesis(my_fit_4_core, hyp, class = NULL))

hyp <- "sd_species__Intercept^2 / (sd_species__Intercept^2 + sigma^2) = 0"
(hyp <- brms::hypothesis(my_fit_4_north, hyp, class = NULL))

hyp <- "sd_species__Intercept^2 / (sd_species__Intercept^2 + sigma^2) = 0"
(hyp <- brms::hypothesis(my_fit_4_south, hyp, class = NULL))

ggsave2("../../figures/main/FIGURE_5_temp100.png", FIGURE_FOUR_TREE, dpi=400,
        height=8, width=8)


