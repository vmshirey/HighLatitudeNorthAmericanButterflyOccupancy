####################################################################################################
# FUNCTION: compute_occ_shift -
# Computes the occupancy interval shift at a collection of sites given the samples derived
# from occupancy-detection modeling.
####################################################################################################
compute_occ_shift <- function(my_data, my_sim_matrix, my_traits, site_type="core",
                              scale="100"){
  
  # Load in the grids depending on scale
  if(scale=="100"){
    my_grid <- readRDS("../output/data/grid_100.rds")
  } else if(scale=="200"){
    my_grid <- readRDS("../output/data/grid_200.rds")
  }
  
  my_occ_list <- list()
  
  # For each species...
  for(sp in 1:nrow(my_traits)){
    
    # Grab the relevant sites and partition into bins
    my_sites <- my_data$my.info$range.list[sp,]
    my_sites_filter <- my_grid %>% 
      dplyr::filter(GID %in% my_sites) %>%
      sf::st_centroid() %>%
      dplyr::mutate(lon=sf::st_coordinates(.)[,1],
                    lat=sf::st_coordinates(.)[,2]) %>%
      dplyr::mutate(latClass=cut(lat, 
                                 breaks=c(quantile(lat, probs=c(0, 0.25, 0.5, 0.75, 1))),
                                 labels=c("Southern", "Core", "Core", "Northern")))
    
    # Filter the sites according to specifications
    if(site_type=="all"){
      
      my_sites <- my_sites
      
    } else if(site_type=="core"){
      
      my_sites <- dplyr::filter(my_sites_filter, latClass=="Core") %>%
        pull(GID)
      
    } else if(site_type=="south"){
      
      my_sites <- dplyr::filter(my_sites_filter, latClass=="Southern") %>%
        pull(GID)
      
    } else if(site_type=="north"){
      
      my_sites <- dplyr::filter(my_sites_filter, latClass=="Northern") %>%
        pull(GID) 
      
    } else{
      message("Invalid site type, must be all, core, south, or north.")
      return(0)
    }
    
    
  
}

####################################################################################################
# FUNCTION: plot_figure_two -
# Creates figure two given a scale of inference.
####################################################################################################
plot_figure_two <- function(sims_matrix, my_traits, 
                            scale="100", makeCommOnly=FALSE){
  # Pull the preciperature parameter estimates per species
  my_temp_samp <- sims_matrix[,grepl("psi.beta.temp", colnames(sims_matrix))]
  
  # Extract summary statistics
  my_temp_mean <- apply(my_temp_samp, 2, mean)
  my_temp_lower <- apply(my_temp_samp, 2, quantile, probs=c(0.025))
  my_temp_upper <- apply(my_temp_samp, 2, quantile, probs=c(0.975))
  
  # Extract community statistics
  my_temp_mean_comm <- mean(my_temp_samp)
  my_temp_int_comm <- quantile(my_temp_samp, probs=c(0.025,0.975))
  
  # Merge the temperature parameter estimates with the species trait data
  my_temp <- data.frame(mean=my_temp_mean, 
                        lower=my_temp_lower, 
                        upper=my_temp_upper) %>%
    cbind(my_traits) %>%
    arrange(mean) %>%
    dplyr::mutate(ordering=row_number()) %>%
    dplyr::mutate(trend=ifelse(sign(lower)==-1 & sign(upper)==1,0,
                               ifelse(sign(lower)==-1 & sign(upper)==-1,-1,1))) %>%
    dplyr::mutate(temperatureClass=cut(ave_temp2, 4))
  
  if(makeCommOnly==TRUE){
    tempEffectsPlotCommOnly <- ggplot()+
      geom_rect(mapping=aes(xmin=my_temp_int_comm[1], 
                            xmax=my_temp_int_comm[2], 
                            ymin=-Inf, ymax=Inf), fill="grey92")+
      geom_vline(xintercept=my_temp_mean_comm, linetype=1, color="grey72")+
      geom_vline(xintercept=0, linetype=2, color="black")+
      geom_pointrange(my_temp, 
                      mapping=aes(y=ordering,
                                  x=mean, 
                                  xmin=lower, 
                                  xmax=upper, 
                                  group=SPID), alpha=0, size=0.9,
                      color=NA)+
      scale_y_continuous(breaks=c(1:nrow(my_traits)), labels=my_temp$speciesCode)+
      labs(y="Species Code", x="Effect of Rising Minimum Temperature (°C) on\nOccupancy Probability")+
      theme_cowplot()+
      theme(legend.position="top",
            axis.text.y=element_text(size=9),
            legend.box="vertical",
            legend.key.width=unit(1, 'cm'),
            legend.background=element_rect(fill=alpha("white", 0.9)),
            plot.background=element_rect(fill="white"))
  }
  
  # Plot the temperature parameter estimates by species
  tempEffectsPlot <- ggplot()+
    geom_rect(mapping=aes(xmin=my_temp_int_comm[1], 
                          xmax=my_temp_int_comm[2], 
                          ymin=-Inf, ymax=Inf), fill="grey92")+
    geom_vline(xintercept=my_temp_mean_comm, linetype=1, color="grey72")+
    geom_vline(xintercept=0, linetype=2, color="black")+
    geom_pointrange(my_temp, 
                    mapping=aes(y=ordering,
                                x=mean, 
                                xmin=lower, 
                                xmax=upper, 
                                group=SPID,
                                color=ave_temp2), alpha=0.9, size=0.9)+
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
  
  # Pull the precipitation parameter estimates per species
  my_precip_samp <- sims_matrix[,grepl("psi.beta.precip", colnames(sims_matrix))]
  my_precip_mean <- apply(my_precip_samp, 2, mean)
  my_precip_lower <- apply(my_precip_samp, 2, quantile, probs=c(0.025))
  my_precip_upper <- apply(my_precip_samp, 2, quantile, probs=c(0.975))
  
  # Extract community statistics
  my_precip_mean_comm <- mean(my_precip_samp)
  my_precip_int_comm <- quantile(my_precip_samp, probs=c(0.025,0.975))
  
  # Merge the precipitation parameter estimates with the species trait data
  my_precip <- data.frame(mean=my_precip_mean, 
                          lower=my_precip_lower, 
                          upper=my_precip_upper) %>%
    cbind(my_traits) %>%
    arrange(mean) %>%
    dplyr::mutate(ordering=row_number()) %>%
    dplyr::mutate(trend=ifelse(sign(lower)==-1 & sign(upper)==1,0,
                               ifelse(sign(lower)==-1 & sign(upper)==-1,-1,1))) %>%
    dplyr::mutate(preciperatureClass=cut(ave_precip2, 4))
  
  if(makeCommOnly==TRUE){
    precipEffectsPlotCommOnly <- ggplot()+
      geom_rect(mapping=aes(xmin=my_precip_int_comm[1], 
                            xmax=my_precip_int_comm[2], 
                            ymin=-Inf, ymax=Inf), fill="grey92")+
      geom_vline(xintercept=my_precip_mean_comm, linetype=1, color="grey72")+
      geom_vline(xintercept=0, linetype=2, color="black")+
      geom_pointrange(my_precip, 
                      mapping=aes(y=ordering,
                                  x=mean, 
                                  xmin=lower, 
                                  xmax=upper, 
                                  group=SPID), size=0.9, alpha=0,
                      color=NA)+
      scale_y_continuous(breaks=c(1:nrow(my_traits)), labels=my_precip$speciesCode)+
      labs(y="Species Code", x="Estimated Effect of Rising Precipitation (cm) on\nOccupancy Probability")+
      theme_cowplot()+
      theme(legend.position="top",
            axis.text.y=element_text(size=9),
            legend.box="vertical",
            legend.key.width=unit(1, 'cm'),
            legend.background=element_rect(fill=alpha("white", 0.9)),
            plot.background=element_rect(fill="white"))
  }
  
  # Plot the precipitation parameter estimates by species
  precipEffectsPlot <- ggplot()+
    geom_rect(mapping=aes(xmin=my_precip_int_comm[1], 
                          xmax=my_precip_int_comm[2], 
                          ymin=-Inf, ymax=Inf), fill="grey92")+
    geom_vline(xintercept=my_precip_mean_comm, linetype=1, color="grey72")+
    geom_vline(xintercept=0, linetype=2, color="black")+
    geom_pointrange(my_precip, 
                    mapping=aes(y=ordering,
                                x=mean, 
                                xmin=lower, 
                                xmax=upper, 
                                group=SPID,
                                color=ave_precip2), size=0.9, alpha=0.9)+
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
  
  figure_two_main <- cowplot::plot_grid(tempEffectsPlot, precipEffectsPlot, ncol=2,
                                        labels=c("(a)", "(b)"))
  
  if(makeCommOnly==TRUE){
    figure_two_mainCommOnly <- cowplot::plot_grid(tempEffectsPlotCommOnly, 
                                                  precipEffectsPlotCommOnly, ncol=2,
                                                  labels=c("(a)", "(b)"))
    ggsave2(paste0("../figures/FIGURE_2_", scale, "_CommOnly.png"), 
            figure_two_mainCommOnly, dpi=400, height=11, width=10)
  }
  ggsave2(paste0("../figures/FIGURE_2_", scale, ".png"), 
          figure_two_main, dpi=400, height=11, width=10)
}


####################################################################################################
# FUNCTION: plot_figure_three -
# Creates figure three given a scale of inference.
####################################################################################################
plot_figure_three <- function(sims_matrix_100, my_data_100, sims_matrix_200, my_data_200, my_traits,
                              scale="100"){
  
}

####################################################################################################
# FUNCTION: plot_figure_four -
# Creates figure four given a scale of inference.
####################################################################################################
plot_figure_four <- function(occ_dx, my_data,
                             scale="100", SPIDs=c(1)){
  
  if(scale=="50"){
    my_grid <- readRDS("../output/data/grid_50.rds")
  } else if(scale=="100"){
    my_grid <- readRDS("../output/data/grid_100.rds")
  } else{
    my_grid <- readRDS("../output/data/grid_200.rds")
  }
  
  plot_list <- list()
  for(i in 1:length(SPIDs)){
    
    my_range_sites <- my_data$my.info$range.list[SPIDs[i],]
    my_kept_sites <- my_data$my.info$kept.sites
    
    my_valid_sites <- intersect(my_range_sites, my_kept_sites)
    
    my_grid <- dplyr::filter(my_grid, GID %in% my_valid_sites)
    
    ave_site_tempRange <- quantile(my_data$my.data$temp[my_valid_sites, 1])[c(1,2)]
    my_cold_sites <- which(between(my_data$my.data$temp[my_valid_sites, 1],
                                   ave_site_tempRange[1], ave_site_tempRange[2]), 
                           arr.ind=TRUE)
    
    ave_site_tempRange <- quantile(my_data$my.data$temp[my_valid_sites, 1])[c(4,5)]
    my_warm_sites <- which(between(my_data$my.data$temp[my_valid_sites, 1],
                                   ave_site_tempRange[1], ave_site_tempRange[2]), 
                           arr.ind=TRUE)
    
    my_grid <- my_grid %>%
      dplyr::mutate(siteType=ifelse(GID %in% my_cold_sites, "Colder", 
                                    ifelse(GID %in% my_warm_sites, "Warmer", "Average")))
    my_grid$siteType <- factor(my_grid$siteType,
                               levels=c("Colder", "Average", "Warmer"))
    
    plot_list[[i]] <- ggplotGrob(ggplot()+
                                   geom_sf(my_grid, mapping=aes(fill=siteType), color=NA)+
                                   geom_sf(basemap, mapping=aes(), fill=NA)+
                                   scale_fill_manual(values=c("#006da2", "grey", "#9f4d48"), name="Site Type")+
                                   theme_void()+
                                   theme(legend.position=c(0.27, 0.4))) 
    
  }
  return(plot_list)
}

####################################################################################################
# FUNCTION: plot_figure_five -
# Creates figure five given a scale of inference.
####################################################################################################
plot_figure_five <- function(occ_dx, scale="100"){
  my_tree <- readRDS("../output/tree_topology.rds")
  A <- ape::vcv.phylo(my_tree)
  
  library(tidybayes)
  
  temp_occ_dx_model <- occ_dx[[1]] %>%
    dplyr::filter(Voltinism!="", !is.na(Voltinism),
                  DisturbanceAffinitySimple!="", !is.na(DisturbanceAffinitySimple))%>%
    dplyr::mutate(phylo=species)
  
  temp_occ_dx_model$Voltinism <- factor(temp_occ_dx_model$Voltinism,
                                        levels=c("Univoltine", "Bivoltine",
                                                 "Multivoltine", "Biennial"))
  temp_occ_dx_model$DisturbanceAffinitySimple <- factor(temp_occ_dx_model$DisturbanceAffinitySimple,
                                                        levels=c("Generalist", "Affinity",
                                                                 "Avoidant"))
  
  trait_model <- brms::brm(meanOcc_lowTemp~
                             DisturbanceAffinity+
                             Voltinism+
                             NumHostplantFamilies+
                             AveWingspan+
                             (1|gr(phylo, cov=A)),
                           data=temp_occ_dx_model,
                           prior = c(
                             prior(normal(0, 10), "b"),
                             prior(normal(0, 50), "Intercept"),
                             prior(student_t(3, 0, 20), "sd"),
                             prior(student_t(3, 0, 20), "sigma")
                           ),
                           family=gaussian(),
                           data2=list(A=A),
                           iter=51000,
                           warmup=1000,
                           thin=100,
                           chains=4,
                           control=list(adapt_delta=0.99))
  
  hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0"
  (hyp <- hypothesis(trait_model, hyp, class = NULL))
  
  hyp_samp <- hyp$samples %>% as.data.frame()
  
  samples_reff <- trait_model %>%
    tidybayes::gather_draws(c(r_phylo[species, term]))
  samples_feff <- trait_model %>%
    tidybayes::gather_draws(c(b_DisturbanceAffinityStrongAffinity, b_DisturbanceAffinityWeakAffinity,
                              b_DisturbanceAffinityWeakAvoidant, b_DisturbanceAffinityStrongAvoidant,
                              b_VoltinismBivoltine, b_VoltinismMultivoltine, b_VoltinismBiennial,
                              b_NumHostplantFamilies, b_AveWingspan, b_Intercept))
  
  figure_5_a_inset <- ggplot()+
    stat_halfeye(hyp_samp, mapping=aes(x=H1), fill="#fb6a4a", alpha=0.5)+
    labs(x=expression(lambda), y="Posterior Density")+
    theme_cowplot()+
    theme(text=element_text(size=9),
          axis.text=element_text(size=9),
          legend.position="none",
          legend.box="vertical",
          legend.key.width=unit(1, 'cm'),
          plot.background=element_rect(fill="white"))
  
  figure_5_a <- ggplot()+
    geom_vline(xintercept=0, linetype=2)+
    stat_interval(samples_reff, mapping=aes(y=.variable, x=.value), .width = c(.5, .89, .95),
                  alpha=0.8)+
    scale_color_brewer(palette="Greys")+
    ggnewscale::new_scale_color()+
    stat_interval(samples_feff, mapping=aes(y=.variable, x=.value), .width = c(.5, .89, .95),
                  alpha=0.8)+
    scale_color_brewer(palette="Blues")+
    scale_y_discrete(limits=rev(c("b_Intercept", "b_VoltinismBivoltine", "b_VoltinismMultivoltine",
                                  "b_VoltinismBiennial", "b_DisturbanceAffinityStrongAvoidant",
                                  "b_DisturbanceAffinityWeakAvoidant", "b_DisturbanceAffinityWeakAffinity",
                                  "b_DisturbanceAffinityStrongAffinity", "b_NumHostplantFamilies", "b_AveWingspan",
                                  "r_phylo")),
                     labels=rev(c("Intercept", "Bivoltine", "Multivoltine", "Biennial",
                                  "Strongly Disturbance Avoidant", "Weakly Disturbance Avoidant",
                                  "Weak Disturbance Affinity", "Strong Disturbance Affinity",
                                  "Number of Hostplant Families", "Average Wingspan",
                                  "Phylogeny")))+
    labs(x="Parameter Estimate", y="")+
    theme_cowplot()+
    theme(axis.text=element_text(size=14),
          legend.position="none",
          legend.box="vertical",
          legend.key.width=unit(1, 'cm'),
          plot.background=element_rect(fill="white"))
  
  samples_reff <- samples_reff %>%
    dplyr::mutate(specificEpithet=sub(".*[.]", "", species)) %>%
    dplyr::mutate(speciesCode=toupper(paste0(substr(species, 1, 3),
                                             substr(specificEpithet, 1, 3)))) %>%
    inner_join(my_traits, by="speciesCode") %>%
    dplyr::mutate(WithinFamilyABCOrdering=factor(WithinFamilyABCOrdering,
                                                 levels=seq(1,70,1))) %>%
    group_by(WithinFamilyABCOrdering)
  
  figure_5_b <- ggplot()+
    geom_hline(yintercept=0, linetype=2)+
    stat_interval(dplyr::filter(samples_reff, Family=="Hesperiidae"), 
                  mapping=aes(x=WithinFamilyABCOrdering, y=.value), 
                  .width = c(.5, .89, .95),
                  alpha=0.8)+
    scale_color_discrete(type=brewer.reds(3))+
    ggnewscale::new_scale_color()+
    stat_interval(dplyr::filter(samples_reff, Family=="Lycaenidae"), 
                  mapping=aes(x=WithinFamilyABCOrdering, y=.value), 
                  .width = c(.5, .89, .95),
                  alpha=0.8)+
    scale_color_discrete(type=brewer.blues(3))+
    ggnewscale::new_scale_color()+
    stat_interval(dplyr::filter(samples_reff, Family=="Nymphalidae"), 
                  mapping=aes(x=WithinFamilyABCOrdering, y=.value), 
                  .width = c(.5, .89, .95),
                  alpha=0.8)+
    scale_color_discrete(type=brewer.greens(3))+
    ggnewscale::new_scale_color()+
    stat_interval(dplyr::filter(samples_reff, Family=="Papilionidae"), 
                  mapping=aes(x=WithinFamilyABCOrdering, y=.value), 
                  .width = c(.5, .89, .95),
                  alpha=0.8)+
    scale_color_discrete(type=brewer.purples(3))+
    ggnewscale::new_scale_color()+
    stat_interval(dplyr::filter(samples_reff, Family=="Pieridae"), 
                  mapping=aes(x=WithinFamilyABCOrdering, y=.value), 
                  .width = c(.5, .89, .95),
                  alpha=0.8)+
    scale_color_discrete(type=brewer.oranges(3))+
    scale_x_discrete(limits=seq(1,70,1), 
                     labels=my_traits[order(my_traits$Family, my_traits$binomial),]$speciesCode)+
    labs(x="Species Code", y="Parameter Estimate")+
    theme_cowplot()+
    theme(axis.text=element_text(size=14),
          axis.text.x=element_text(size=9, angle = 90, vjust = 0.5, hjust=1),
          legend.position="none",
          legend.box="vertical",
          legend.key.width=unit(1, 'cm'),
          plot.background=element_rect(fill="white"))
  
  figure_5_top <- cowplot::plot_grid(NULL, figure_5_a, NULL, ncol=3, rel_widths=c(0.1,0.8,0.1),
                                     labels=c("", "(a)", ""))
  figure_5 <- cowplot::plot_grid(figure_5_top, figure_5_b, nrow=2, labels=c("", "(b)"))
  
  ggsave2(paste0("../figures/FIGURE_5_", scale, ".png"), 
          figure_5, dpi=400, height=10, width=10)
}