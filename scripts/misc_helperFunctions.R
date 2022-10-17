####################################################################################################
# FUNCTION: compute_occ_shift -
# Computes the occupancy interval shift at a collection of sites given the samples derived
# from occupancy-detection modeling.
####################################################################################################
compute_occ_shift <- function(my_data, my_sim_matrix, 
                              my_traits, site_type="core", scale="100"){
  
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
                                 labels=c("Southern", "Core1", "Core2", "Northern"),
                                 include.lowest=TRUE))
    
    # Filter the sites according to specifications
    if(site_type=="all"){
      
      my_sites <- my_sites
      
    } else if(site_type=="core"){
      
      my_sites <- dplyr::filter(my_sites_filter, latClass %in% c("Core1", "Core2")) %>%
        pull(GID)
      
    } else if(site_type=="south"){
      
      my_sites <- dplyr::filter(my_sites_filter, latClass %in% c("Southern")) %>%
        pull(GID)
      
    } else if(site_type=="north"){
      
      my_sites <- dplyr::filter(my_sites_filter, latClass %in% c("Northern")) %>%
        pull(GID) 
      
    } else{
      message("Invalid site type, must be all, core, south, or north.")
      return(0)
    }
    
    meanPCA1_Occ1 <- mean(my_data$my.data$PCA1[my_sites, 1], na.rm=TRUE)
    meanPCA2_Occ1 <- mean(my_data$my.data$PCA2[my_sites, 1], na.rm=TRUE)
    
    meanPCA1_Occ10 <- mean(my_data$my.data$PCA1[my_sites, 10], na.rm=TRUE)
    meanPCA2_Occ10 <- mean(my_data$my.data$PCA2[my_sites, 10], na.rm=TRUE)
    
    meanProb_Occ1 <- plogis(my_sim_matrix[,"mu.psi.0"]+
                              my_sim_matrix[,"psi.area"]*mean(my_data$my.data$gridarea)+
                              my_sim_matrix[,sprintf("psi.beta.pca1[%d]", sp)]*meanPCA1_Occ1+
                              my_sim_matrix[,sprintf("psi.beta.pca2[%d]", sp)]*meanPCA2_Occ1)
    
    meanProb_Occ10 <- plogis(my_sim_matrix[,"mu.psi.0"]+
                               my_sim_matrix[,"psi.area"]*mean(my_data$my.data$gridarea)+
                               my_sim_matrix[,sprintf("psi.beta.pca1[%d]", sp)]*meanPCA1_Occ10+
                               my_sim_matrix[,sprintf("psi.beta.pca2[%d]", sp)]*meanPCA2_Occ10)
    
    meanProb_dx <- meanProb_Occ10-meanProb_Occ1
    
    my_occ_list[[sp]] <- c(sp,
                           mean(meanProb_dx),
                           quantile(meanProb_dx, c(0.05, 0.975))[1],
                           quantile(meanProb_dx, c(0.05, 0.975))[2],
                           sd(meanProb_dx))
  }
  my_occ_df <- do.call(rbind, my_occ_list) %>% as.data.frame()
  colnames(my_occ_df) <- c("SPID", "mean", "lower", "upper", "sd")
  
  return(my_occ_df)
}

####################################################################################################
# FUNCTION: plot_figure_three -
# Creates figure two given a scale of inference.
####################################################################################################
plot_figure_three <- function(sims_matrix, my_traits, 
                              scale="100", makeCommOnly=FALSE){
  # Pull the preciperature parameter estimates per species
  my_pca1_samp <- sims_matrix[,grepl("psi.beta.pca1", colnames(sims_matrix))]
  
  # Extract summary statistics
  my_pca1_mean <- apply(my_pca1_samp, 2, mean)
  my_pca1_lower <- apply(my_pca1_samp, 2, quantile, probs=c(0.025))
  my_pca1_upper <- apply(my_pca1_samp, 2, quantile, probs=c(0.975))
  
  # Extract community statistics
  my_pca1_mean_comm <- mean(my_pca1_samp)
  my_pca1_int_comm <- quantile(my_pca1_samp, probs=c(0.025,0.975))
  
  # Merge the temperature parameter estimates with the species trait data
  my_pca1 <- data.frame(mean=my_pca1_mean, 
                        lower=my_pca1_lower, 
                        upper=my_pca1_upper) %>%
    cbind(my_traits) %>%
    arrange(mean) %>%
    dplyr::mutate(ordering=row_number()) %>%
    dplyr::mutate(trend=ifelse(sign(lower)==-1 & sign(upper)==1,0,
                               ifelse(sign(lower)==-1 & sign(upper)==-1,-1,1))) %>%
    dplyr::mutate(temperatureClass=cut(ave_temp2, 4, include.lowest=TRUE),
                  precipitationClass=cut(ave_precip2, 4, include.lowest=TRUE)) %>%
    dplyr::mutate(trend=ifelse(temperatureClass=="[-12,-5.3]" & precipitationClass=="[195,426]", "Cold/Dry", 
                               ifelse(temperatureClass=="[-12,-5.3]" & precipitationClass=="(426,656]", "Cold/Average",
                                      ifelse(temperatureClass=="[-12,-5.3]" & precipitationClass=="(886,1.12e+03]", "Cold/Wet",
                                             ifelse(temperatureClass=="(8.03,14.7]" & precipitationClass=="[195,426]", "Hot/Dry", 
                                                    ifelse(temperatureClass=="(8.03,14.7]" & precipitationClass=="(426,656]", "Hot/Average",
                                                           ifelse(temperatureClass=="(8.03,14.7]" & precipitationClass=="(886,1.12e+03]", "Hot/Wet", "Average/Average")))))))
  
  # Plot the temperature parameter estimates by species
  pca1_effects_plot <- ggplot()+
    geom_rect(mapping=aes(ymin=my_pca1_int_comm[1], 
                          ymax=my_pca1_int_comm[2], 
                          xmin=-Inf, xmax=Inf), fill="grey92")+
    geom_hline(yintercept=my_pca1_mean_comm, linetype=1, color="grey72")+
    geom_hline(yintercept=0, linetype=2, color="black")+
    geom_pointrange(my_pca1, 
                    mapping=aes(x=ordering,
                                y=mean, 
                                ymin=lower, 
                                ymax=upper, 
                                group=SPID,
                                color=temperatureClass, fill=temperatureClass,
                                shape=precipitationClass), alpha=1, size=0.9)+
    scale_x_continuous(breaks=c(1:nrow(my_traits)), labels=my_pca1$speciesCode)+
    scale_color_manual(values=c("[-12,-5.3]"="#089392",
                                "(-5.3,1.37]"="#9dcd84",
                                "(1.37,8.03]"="#eab672",
                                "(8.03,14.7]"="#cf5a7f"),
                       labels=c("Coldest", "Colder", "Warmer", "Warmest"),
                       name="Temperature Class",
                       guide = guide_legend(
                         direction = "horizontal",
                         title.position = "top"
                       ))+
    scale_fill_manual(values=c("[-12,-5.3]"="#089392",
                               "(-5.3,1.37]"="#9dcd84",
                               "(1.37,8.03]"="#eab672",
                               "(8.03,14.7]"="#cf5a7f"),
                      labels=c("Coldest", "Colder", "Warmer", "Warmest"),
                      name="Temperature Class",
                      guide = guide_legend(
                        direction = "horizontal",
                        title.position = "top"
                      ))+
    scale_shape_manual(values=c("[195,426]"=25,
                                "(426,656]"=22,
                                "(656,886]"=21,
                                "(886,1.12e+03]"=24),
                       labels=c("Driest", "Drier", "Wetter", "Wettest"),
                       name="Precipitation Class",
                       guide = guide_legend(
                         direction = "horizontal",
                         title.position = "top",
                         override.aes = list(fill = "black")
                       ))+
    labs(x="", y="Effect of Climatic PCA1 Score\non Occupancy Probability")+
    theme_cowplot()+
    theme(legend.position=c(0.025, 0.85),
          legend.background=element_blank(),
          axis.text.x=element_text(size=9, angle=90, vjust = 0.5, hjust=1),
          legend.box="vertical",
          legend.key.width=unit(1, 'cm'),
          plot.background=element_rect(fill="white"))
  
  # Pull the precipitation parameter estimates per species
  my_pca2_samp <- sims_matrix[,grepl("psi.beta.pca2", colnames(sims_matrix))]
  my_pca2_mean <- apply(my_pca2_samp, 2, mean)
  my_pca2_lower <- apply(my_pca2_samp, 2, quantile, probs=c(0.025))
  my_pca2_upper <- apply(my_pca2_samp, 2, quantile, probs=c(0.975))
  
  # Extract community statistics
  my_pca2_mean_comm <- mean(my_pca2_samp)
  my_pca2_int_comm <- quantile(my_pca2_samp, probs=c(0.025,0.975))
  
  # Merge the precipitation parameter estimates with the species trait data
  my_pca2 <- data.frame(mean=my_pca2_mean, 
                        lower=my_pca2_lower, 
                        upper=my_pca2_upper) %>%
    cbind(my_traits) %>%
    arrange(mean) %>%
    dplyr::mutate(ordering=row_number()) %>%
    dplyr::mutate(trend=ifelse(sign(lower)==-1 & sign(upper)==1,0,
                               ifelse(sign(lower)==-1 & sign(upper)==-1,-1,1))) %>%
    dplyr::mutate(temperatureClass=cut(ave_temp2, 4, include.lowest=TRUE),
                  precipitationClass=cut(ave_precip2, 4, include.lowest=TRUE))
  
  # Plot the precipitation parameter estimates by species
  pca2_effects_plot <- ggplot()+
    geom_rect(mapping=aes(ymin=my_pca2_int_comm[1], 
                          ymax=my_pca2_int_comm[2], 
                          xmin=-Inf, xmax=Inf), fill="grey92")+
    geom_hline(yintercept=my_pca2_mean_comm, linetype=1, color="grey72")+
    geom_hline(yintercept=0, linetype=2, color="black")+
    geom_pointrange(my_pca2, 
                    mapping=aes(x=ordering,
                                y=mean, 
                                ymin=lower, 
                                ymax=upper, 
                                group=SPID,
                                color=temperatureClass, fill=temperatureClass,
                                shape=precipitationClass), alpha=1, size=0.9)+
    scale_x_continuous(breaks=c(1:nrow(my_traits)), labels=my_pca2$speciesCode)+
    scale_color_manual(values=c("[-12,-5.3]"="#089392",
                                "(-5.3,1.37]"="#9dcd84",
                                "(1.37,8.03]"="#eab672",
                                "(8.03,14.7]"="#cf5a7f"),
                       labels=c("Coldest", "Colder", "Warmer", "Warmest"),
                       name="Temperature Class",
                       guide = guide_legend(
                         direction = "horizontal",
                         title.position = "top"
                       ))+
    scale_fill_manual(values=c("[-12,-5.3]"="#089392",
                                "(-5.3,1.37]"="#9dcd84",
                                "(1.37,8.03]"="#eab672",
                                "(8.03,14.7]"="#cf5a7f"),
                       labels=c("Coldest", "Colder", "Warmer", "Warmest"),
                       name="Temperature Class",
                       guide = guide_legend(
                         direction = "horizontal",
                         title.position = "top"
                       ))+
    scale_shape_manual(values=c("[195,426]"=25,
                                "(426,656]"=22,
                                "(656,886]"=21,
                                "(886,1.12e+03]"=24),
                       labels=c("Driest", "Drier", "Wetter", "Wettest"),
                       name="Precipitation Class",
                       guide = guide_legend(
                         direction = "horizontal",
                         title.position = "top"
                       ))+
    labs(x="Species Code", y="Effect of Climatic PCA2 Score\non Occupancy Probability")+
    theme_cowplot()+
    theme(legend.position="none",
          axis.text.x=element_text(size=9, angle=90, vjust = 0.5, hjust=1),
          plot.background=element_rect(fill="white"))
  
  # Place each panel together into the main figure and save
  figure_two <- cowplot::plot_grid(pca1_effects_plot, pca2_effects_plot, nrow=2,
                                        labels=c("(a)", "(b)"))
  ggsave2(paste0("../figures/main/FIGURE_3_", scale, ".png"), 
          figure_two, dpi=400, height=10, width=16)
}


####################################################################################################
# FUNCTION: plot_figure_two -
# Creates figure two given a scale of inference.
####################################################################################################
plot_figure_two <- function(occ_dx_core, occ_dx_north, occ_dx_south, 
                            my_traits, sims_matrix,
                              scale="100"){
  
  occ_dx_core$N_trend <- occ_dx_north[order(match(occ_dx_north[[1]], occ_dx_core$SPID)),] %>%
    pull(trend)
  occ_dx_core$S_trend <- occ_dx_south[order(match(occ_dx_south[[1]], occ_dx_core$SPID)),] %>%
    pull(trend)
    
  figure_3 <- ggplot()+
    geom_hline(yintercept=0, linetype=2)+
    geom_pointrange(occ_dx_core,
                    mapping=aes(x=order, y=mean, ymin=lower, ymax=upper, group=SPID,
                                color=as.factor(trend)),
                    size=1)+
    scale_color_discrete_diverging(palette="Berlin", rev=TRUE,
                                   labels=c("Declining", "Overlap Zero", "Increasing"),
                                   name="Trend")+
    scale_x_continuous(breaks=seq(1:nrow(my_traits)), labels=occ_dx$speciesCode,
                       name="Species Code")+
    scale_y_continuous(labels=scales::percent, 
                       name="Mean In-range Occupancy Shift\n(1970-1974 to 2015-2019)")+
    theme_cowplot()+
    theme(axis.text.x=element_text(size=9, angle=90, vjust=0.5, hjust=1),
          legend.position=c(0.025, 0.95),
          legend.direction="horizontal",
          plot.background=element_rect(fill="white"))
  
  ggsave2(paste0("../figures/main/FIGURE_2_", scale, ".png"), figure_3, dpi=400, height=6, width=12)
  
}

####################################################################################################
# FUNCTION: plot_figure_five -
# Creates figure five given a scale of inference.
####################################################################################################
plot_figure_five <- function(occ_dx, scale="100"){
  my_tree <- readRDS("../output/tree_topology.rds")
  A <- ape::vcv.phylo(my_tree)
  
  library(tidybayes)
  
  temp_occ_dx_model <- occ_dx %>%
    dplyr::filter(Voltinism!="", !is.na(Voltinism),
                  DisturbanceAffinitySimple!="", !is.na(DisturbanceAffinitySimple))%>%
    dplyr::mutate(phylo=species)
  
  temp_occ_dx_model$Voltinism <- factor(temp_occ_dx_model$Voltinism,
                                        levels=c("Univoltine", "Bivoltine",
                                                 "Multivoltine", "Biennial"))
  temp_occ_dx_model$DisturbanceAffinitySimple <- factor(temp_occ_dx_model$DisturbanceAffinitySimple,
                                                        levels=c("Generalist", "Affinity",
                                                                 "Avoidant"))
  
  
  trait_model <- brms::brm(mean|mi(sd)~DisturbanceAffinity+
                             Voltinism+
                             scale(NumHostplantFamilies)+
                             scale(AveWingspan)+
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
                              b_scaleNumHostplantFamilies, b_scaleAveWingspan, b_Intercept))
  
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
                                  "b_DisturbanceAffinityStrongAffinity", "b_scaleNumHostplantFamilies",
                                  "b_scaleAveWingspan","r_phylo")),
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