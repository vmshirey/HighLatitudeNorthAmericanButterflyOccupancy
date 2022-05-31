# Create a function to visualize figure one in the manuscript
visFigureOne <- function(){
  
  sf_use_s2(FALSE)
  
  basemap <- st_read("../data/shapefile/ne_10m_land.shp") %>%
    st_crop(xmin=-180, xmax=-50, ymin=45, ymax=90)
  occ_dat <- readRDS("../output/occ_sf.rds")
  
  occ_sum <- occ_dat %>%
    st_drop_geometry() %>%
    dplyr::group_by(year, basisOfRecord) %>%
    dplyr::mutate(count=n()) %>%
    dplyr::select(year, basisOfRecord, count) %>%
    dplyr::ungroup() %>% unique() %>%
    dplyr::filter(between(year, 1970, 2019))
  
  occ_pie <- data.frame(group=c("Community Events", "Singleton Events"),
                        perc=c(1-0.3373132, 0.3373132)) %>%
    dplyr::mutate(csum = rev(cumsum(rev(perc))), 
                  pos = perc/2 + lead(csum, 1),
                  pos = if_else(is.na(pos), perc/2, pos))
  
  grid_100 <- readRDS("../output/data/grid_100.rds")
  temp_ras <- readRDS("../output/data/temp100.rds") %>%
    left_join(grid_100) %>%
    dplyr::mutate(temp2019=rowMeans(dplyr::select(temp_ras, matches("201")), na.rm = TRUE),
                  temp1970=rowMeans(dplyr::select(temp_ras, matches("197")), na.rm = TRUE)) %>%
    dplyr::mutate(tempDiff=temp2019-temp1970) %>% st_as_sf()
  
  FIGURE_ONE_INSET_ONE <- ggplot()+
    geom_line(occ_sum, mapping=aes(x=year, y=count, color=basisOfRecord))+
    geom_point(occ_sum, mapping=aes(x=year, y=count, color=basisOfRecord))+
    geom_text_repel(occ_sum %>% group_by(year, basisOfRecord) %>% dplyr::filter(year==2019), 
                    mapping=aes(x=year+1, y=count, color=basisOfRecord, label=basisOfRecord),
                    fontface="bold", direction="y", hjust=0, xlim = c(2020.8, NA),
                    segment.size=1, segment.linetype="dotted", box.padding=0.4,
                    segment.curvature = -0.1, segment.ncp = 3, segment.angle = 20)+
    coord_cartesian(clip = "off")+
    scale_color_manual(values=c("#1fa07a", "#d95f02"), labels=c("Observation", "Specimen"))+
    scale_x_continuous(limits=c(1970,2050), breaks=seq(1970, 2020, by=10), expand=c(0,0))+
    scale_y_log10()+
    labs(x="Year of Observation/Collection", y="Frequency of Records")+
    theme_cowplot()+
    theme(legend.position="none", 
          plot.background=element_rect(fill=alpha("white", 0.4)), 
          panel.background=element_rect(fill=alpha("white", 0.4)))
  
  FIGURE_ONE_INSET_TWO <- ggplot()+
    geom_bar(occ_pie, mapping=aes(x="", y=perc, fill=group), 
             stat="identity", width=1, color="white", alpha=0.4)+
    geom_label_repel(occ_pie, mapping=aes(x="", y=pos, color=group,
                                          label=paste0(group, ": ", round(perc*100, 2), "%")),
                     size=4.5, nudge_x=0.5, show.legend=FALSE)+
    scale_fill_manual(values=c("#9694bf", "#e72a8a"))+
    scale_color_manual(values=c("#9694bf", "#e72a8a"))+
    coord_polar("y", start=0)+
    theme_void()+
    theme(legend.position="none")
  
  FIGURE_ONE_INSET_THREE <- ggplot()+
    geom_sf(temp_ras, mapping=aes(fill=tempDiff), color=NA)+
    scale_fill_continuous_divergingx(palette="RdYlBu", rev=TRUE, 
                                     name="Change in Average Minimum Temperature\nfrom 1970s to 2010s (Celsius)",
                                     guide=guide_colorbar(title.position = "top"),
                                     limits=c(-2,4), mid=0, breaks=c(-2,0,2,4))+
    theme_map()+
    theme(legend.position = c(0.05,0.05),
          legend.direction = "horizontal",
          legend.box = "horizontal",
          legend.key.width= unit(1, 'cm'))
  
  FIGURE_ONE_MAP <- ggplot()+
    geom_sf(basemap, mapping=aes(), fill="#BBBBBB", color=NA)+
    geom_sf(dplyr::slice_sample(occ_dat, prop=0.5), 
            mapping=aes(color=basisOfRecord), shape=20, alpha=0.5,
            size=0.001)+
    scale_color_manual(values=c("#1fa07a", "#d95f02"), labels=c("Observation", "Specimen"),
                       name="Type of Record")+
    theme_map()+
    theme(legend.position="none")
  
  FIGURE_ONE <- ggdraw()+
    draw_plot(FIGURE_ONE_INSET_THREE)+
    draw_plot(FIGURE_ONE_INSET_ONE, x=0.5, y=0.67, width=0.5, height=0.3)+
    draw_plot(FIGURE_ONE_INSET_TWO, x=0.075, y=0.2, width=0.3, height=0.3)+
    draw_plot(FIGURE_ONE_MAP, x=0.8, y=0.73, width=0.2, height=0.25)+
    draw_plot_label(c("(a)", "(b)", "(c)"),
                    x=c(0.15,0.47,0.1),
                    y=c(0.8,0.95,0.45),
                    size=16, fontface="bold")
  
  # Save the final figure off
  ggsave2("../output/main/final/FIGURE_ONE.png", FIGURE_ONE, dpi=350, width=14, height=10)
  
}

# Create a function to visualize figure two in the manuscript
visFigureTwo <- function(my_res, my_data, my_traits){
  
  # Create species codes for plotting
  my_traits <- my_traits %>%
    dplyr::mutate(specificEpithet=word(species, 2)) %>%
    dplyr::mutate(speciesCode=toupper(paste0(substr(species, 1, 3),
                                             substr(specificEpithet, 1, 3)))) %>%
    dplyr::mutate(tempClass=cut(ave_temp2, 4)) %>%
    dplyr::mutate(tempClass=ifelse(tempClass==levels(tempClass)[1], "Coldest Ranges",
                                   ifelse(tempClass==levels(tempClass)[4], "Warmest Ranges", "All Other Ranges")))
  
  # Decompress the R2JAGS output into a single object for summaries and plotting
  message("Reconfiguring parallel output as single JAGS run...")
  my_res$chains <- mcmc.list(as.mcmc(my_res[[1]]$BUGSoutput$sims.matrix),
                             as.mcmc(my_res[[2]]$BUGSoutput$sims.matrix),
                             as.mcmc(my_res[[3]]$BUGSoutput$sims.matrix))
  my_res$summary <- MCMCvis::MCMCsummary(my_res$chains)
  my_res$sims.matrix <- as.matrix(my_res$chains)
  my_sims_mat <- my_res$sims.matrix
  
  # Pull the model summary
  my_sum <- my_res$summary %>% as.data.frame() %>%
    rownames_to_column("Parameter")
  message(".")
  
  # Pull all the samples
  my_samp <- my_res$sims.matrix
  my_samp_beta_temp <- data.frame(samples=my_samp[,"psi.temp.env"])
  my_samp_beta_temp_trait <- data.frame(samples=my_samp[,"psi.temp.trait"])
  message(".")
  
  # Pull the temp effects
  my_temp <- dplyr::filter(my_sum, grepl("psi.beta.temp", Parameter)) %>%
    cbind(my_traits) %>%
    arrange(mean) %>%
    dplyr::mutate(ordering=row_number()) 
  message(".")
  
  # Visualize the temp effects as a forest plot
  message("Visualizing posterior densities for effect subcomponents...")
  beta_env_temp <- ggplot()+
    geom_vline(xintercept=0, linetype=2, color="white")+
    geom_density(my_samp_beta_temp, mapping=aes(x=samples), color="white", fill="white", alpha=0.2)+
    annotate("text", x=Inf, y=Inf, parse=TRUE,
             label="beta[1]", color="white", fontface="bold", size=6, hjust=1, vjust=1)+
    scale_x_continuous(limits=c(-0.6,0.6))+
    labs(x="Effect Estimate", y="Posterior Density")+
    theme_cowplot()+
    theme(text = element_text(color="white"),
          axis.line = element_line(color="white"),
          axis.ticks = element_line(color="white"),
          axis.text.y = element_text(color="white"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.35, size=9, color="white"),
          plot.background = element_rect(fill="#202020"),
          legend.position = c(0.05,0.9),
          legend.direction = "horizontal",
          legend.box = "horizontal",
          legend.key.width= unit(1, 'cm'))
  
  beta_trait_temp <- ggplot()+
    geom_vline(xintercept=0, linetype=2, color="white")+
    geom_density(my_samp_beta_temp_trait, mapping=aes(x=samples), color="white", fill="white", alpha=0.2)+
    annotate("text", x=Inf, y=Inf, parse=TRUE,
             label="beta[2]", color="white", fontface="bold", size=6, hjust=1, vjust=1)+
    scale_x_continuous(limits=c(-0.6,0.6))+
    labs(x="Effect Estimate", y="")+
    theme_cowplot()+
    theme(text = element_text(color="white"),
          axis.line = element_line(color="white"),
          axis.ticks = element_line(color="white"),
          axis.text.y = element_text(color="white"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.35, size=9, color="white"),
          plot.background = element_rect(fill="#202020"),
          legend.position = c(0.05,0.9),
          legend.direction = "horizontal",
          legend.box = "horizontal",
          legend.key.width= unit(1, 'cm'))
  message(".")
  
  # RESPONSE CURVES FOR ALL SPECIES
  message("Plotting response curves for all species...")
  spp_lines <- list()
  for(ss in 1:nrow(my_traits)){
    get.y.val <- function(dd){
      chains <- plogis(
        my_sims_mat[,sprintf("sp.psi.0[%d]", ss)]+
          my_sims_mat[,"p.area"]*mean(my_data$my.constants$gridarea)+
          my_sims_mat[,sprintf("psi.beta.temp[%d]", ss)]*dd+
          (my_sims_mat[,"psi.precip.env"]+
             my_sims_mat[,"psi.precip.trait"]*mean(my_data$my.constants$precip)))
      c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
    }
    x.vals <- seq(from=-20, to=20, length=100)
    
    y.vals <- lapply(x.vals, get.y.val)
    y.vals2_spp <- do.call(rbind, y.vals)
    
    y.vals2_spp <- y.vals2_spp %>% as.data.frame() %>%
      dplyr::mutate(spp=ss, x=x.vals)
    spp_lines[[ss]] <- y.vals2_spp
  }
  y.vals2_spp <- do.call(rbind, spp_lines)
  y.vals2_spp <- y.vals2_spp %>% left_join(my_traits, by=c("spp"="SPID"))
  
  get.y.val <- function(dd){
    chains <- plogis(mean(my_sims_mat[,grepl("sp.psi.0", colnames(my_sims_mat))])+
                       my_sims_mat[,"p.area"]*mean(my_data$my.constants$gridarea)+
                       (my_sims_mat[,"psi.temp.env"]+
                          my_sims_mat[,"psi.temp.trait"])*dd+
                       (my_sims_mat[,"psi.precip.env"]+
                          my_sims_mat[,"psi.precip.trait"]*mean(my_data$my.constants$precip)))
    
    c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
  }
  x.vals <- seq(from=-20, to=20, length=100)
  y.vals <- lapply(x.vals, get.y.val)
  y.vals2_com <- do.call(rbind, y.vals) %>% as.data.frame()
  y.vals2_com$x <- x.vals
  
  allResponseCurves <- ggplot()+
    geom_line(y.vals2_spp, mapping=aes(x=x, y=mean, group=spp), 
              alpha=0.05, color="white")+
    labs(x="", y="Occupancy Probability")+
    geom_line(y.vals2_com, mapping=aes(x=x, y=mean), color="white")+
    geom_line(y.vals2_com, mapping=aes(x=x, y=`2.5%`), color="white", linetype=2)+
    geom_line(y.vals2_com, mapping=aes(x=x, y=`97.5%`), color="white", linetype=2)+
    scale_y_continuous(labels=scales::percent)+
    theme_cowplot()+
    theme(text = element_text(color="white"),
          axis.line = element_line(color="white"),
          axis.ticks = element_line(color="white"),
          axis.text.y = element_text(color="white"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.35, size=9, color="white"),
          plot.background = element_rect(fill="#202020"),
          legend.position = c(0.05,0.9),
          legend.direction = "horizontal",
          legend.box = "horizontal",
          legend.key.width= unit(1, 'cm'))
  message(".")
  
  # RESPONSE CURVES FOR COLDEST ADAPTED SPECIES  
  message("Plotting response curves for coldest-adapted species...")
  spp_lines <- list()
  coldest_species <- dplyr::filter(my_traits, tempClass=="Coldest Ranges")
  for(ss in 1:nrow(coldest_species)){
    get.y.val <- function(dd){
      chains <- plogis(
        my_sims_mat[,sprintf("sp.psi.0[%d]", coldest_species$SPID[ss])]+
          my_sims_mat[,"p.area"]*mean(my_data$my.constants$gridarea)+
          my_sims_mat[,sprintf("psi.beta.temp[%d]", coldest_species$SPID[ss])]*dd+
          my_sims_mat[,sprintf("psi.beta.precip[%d]", coldest_species$SPID[ss])]*
          mean(my_data$my.constants$precip))
      c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
    }
    x.vals <- seq(from=-20, to=20, length=100)
    
    y.vals <- lapply(x.vals, get.y.val)
    y.vals2_spp <- do.call(rbind, y.vals)
    
    y.vals2_spp <- y.vals2_spp %>% as.data.frame() %>%
      dplyr::mutate(spp=ss, x=x.vals)
    spp_lines[[ss]] <- y.vals2_spp
  }
  y.vals2_spp <- do.call(rbind, spp_lines)
  y.vals2_spp <- y.vals2_spp %>% left_join(my_traits, by=c("spp"="SPID"))
  
  get.y.val <- function(dd){
    chains <- plogis(mean(my_sims_mat[,sprintf("sp.psi.0[%d]", coldest_species$SPID)])+
                       my_sims_mat[,"p.area"]*mean(my_data$my.constants$gridarea)+
                       mean(my_sims_mat[,sprintf("psi.beta.temp[%d]", coldest_species$SPID)])*dd+
                       mean(my_sims_mat[,sprintf("psi.beta.precip[%d]", coldest_species$SPID)])*
                       mean(my_data$my.constants$precip))
    
    c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
  }
  x.vals <- seq(from=-20, to=20, length=100)
  y.vals <- lapply(x.vals, get.y.val)
  y.vals2_com <- do.call(rbind, y.vals) %>% as.data.frame()
  y.vals2_com$x <- x.vals
  
  coldResponseCurves <- ggplot()+
    geom_line(y.vals2_spp, mapping=aes(x=x, y=mean, group=spp), 
              alpha=0.05, color="#04d9ff")+
    labs(x="", y="")+
    geom_line(y.vals2_com, mapping=aes(x=x, y=mean), color="#04d9ff")+
    geom_line(y.vals2_com, mapping=aes(x=x, y=`2.5%`), color="#04d9ff", linetype=2)+
    geom_line(y.vals2_com, mapping=aes(x=x, y=`97.5%`), color="#04d9ff", linetype=2)+
    scale_y_continuous(labels=scales::percent)+
    theme_cowplot()+
    theme(text = element_text(color="white"),
          axis.line = element_line(color="white"),
          axis.ticks = element_line(color="white"),
          axis.text.y = element_text(color="white"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.35, size=9, color="white"),
          plot.background = element_rect(fill="#202020"),
          legend.position = c(0.05,0.9),
          legend.direction = "horizontal",
          legend.box = "horizontal",
          legend.key.width= unit(1, 'cm'))
  message(".")
  
  # RESPONSE CURVES FOR WARMEST ADAPTED SPECIES
  message("Plotting response curves for warmest-adapted species...")
  spp_lines <- list()
  warmest_species <- dplyr::filter(my_traits, tempClass=="Warmest Ranges")
  for(ss in 1:nrow(warmest_species)){
    get.y.val <- function(dd){
      chains <- plogis(
        my_sims_mat[,sprintf("sp.psi.0[%d]", warmest_species$SPID[ss])]+
          my_sims_mat[,"p.area"]*mean(my_data$my.constants$gridarea)+
          my_sims_mat[,sprintf("psi.beta.temp[%d]", warmest_species$SPID[ss])]*dd+
          mean(my_sims_mat[,sprintf("psi.beta.precip[%d]", warmest_species$SPID[ss])])*
          mean(my_data$my.constants$precip))
      c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
    }
    x.vals <- seq(from=-20, to=20, length=100)
    
    y.vals <- lapply(x.vals, get.y.val)
    y.vals2_spp <- do.call(rbind, y.vals)
    
    y.vals2_spp <- y.vals2_spp %>% as.data.frame() %>%
      dplyr::mutate(spp=ss, x=x.vals)
    spp_lines[[ss]] <- y.vals2_spp
  }
  y.vals2_spp <- do.call(rbind, spp_lines)
  y.vals2_spp <- y.vals2_spp %>% left_join(my_traits, by=c("spp"="SPID"))
  
  get.y.val <- function(dd){
    chains <- plogis(mean(my_sims_mat[,sprintf("sp.psi.0[%d]", warmest_species$SPID)])+
                       my_sims_mat[,"p.area"]*mean(my_data$my.constants$gridarea)+
                       mean(my_sims_mat[,sprintf("psi.beta.temp[%d]", warmest_species$SPID)])*dd+
                       mean(my_sims_mat[,sprintf("psi.beta.precip[%d]", warmest_species$SPID)])*
                       mean(my_data$my.constants$precip))
    
    c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
  }
  x.vals <- seq(from=-20, to=20, length=100)
  y.vals <- lapply(x.vals, get.y.val)
  y.vals2_com <- do.call(rbind, y.vals) %>% as.data.frame()
  y.vals2_com$x <- x.vals
  
  warmResponseCurves <- ggplot()+
    geom_line(y.vals2_spp, mapping=aes(x=x, y=mean, group=spp), 
              alpha=0.05, color="#ff9f37")+
    labs(x="", y="")+
    geom_line(y.vals2_com, mapping=aes(x=x, y=mean), color="#ff9f37")+
    geom_line(y.vals2_com, mapping=aes(x=x, y=`2.5%`), color="#ff9f37", linetype=2)+
    geom_line(y.vals2_com, mapping=aes(x=x, y=`97.5%`), color="#ff9f37", linetype=2)+
    scale_y_continuous(labels=scales::percent)+
    theme_cowplot()+
    theme(text = element_text(color="white"),
          axis.line = element_line(color="white"),
          axis.ticks = element_line(color="white"),
          axis.text.y = element_text(color="white"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.35, size=9, color="white"),
          plot.background = element_rect(fill="#202020"),
          legend.position = c(0.05,0.9),
          legend.direction = "horizontal",
          legend.box = "horizontal",
          legend.key.width= unit(1, 'cm'))
  message(".")
  
  # CONSTRUCT THE MAIN FOREST PLOT
  message("Plotting the main forest plot of effect estimates...")
  main_temp <- ggplot()+
    geom_hline(yintercept=0, linetype=2, color="white")+
    geom_pointrange(my_temp, 
                    mapping=aes(x=ordering,
                                y=mean, 
                                ymin=`2.5%`, 
                                ymax=`97.5%`, 
                                group=SPID, color=ave_temp2, shape=tempClass))+
    annotate("text", x=17, y=1.75, parse=TRUE, hjust=0, vjust=1,
               label="Effect == beta[1] + beta[2]*TemperatureTrait", 
             size=6, color="white")+
    scale_x_continuous(breaks=c(1:nrow(my_traits)), labels=my_temp$speciesCode)+
    scale_y_continuous(limits=c(-2,2))+
    scale_shape_manual(values=c(1, 15, 19), guide=guide_legend(override.aes=list(fill="white", color="white")),
                       name="Range Wide Temperature Class:")+
    scale_color_continuous_divergingx(palette="RdYlBu", rev=TRUE,
                                      guide=guide_colorbar(title.position = "top"),
                                      name="Average Annual Range-wide Temperature (Celsius)",
                                      mid=mean(my_traits$ave_temp2), breaks=c(-15, -10, -5, 0, 5, 10, 15),
                                      limits=c(-15,15))+
    labs(x="Species Code", y="Estimated Effect of Minimum Temperature on Occupancy Probability")+
    theme_cowplot()+
    theme(text = element_text(color="white"),
          axis.line = element_line(color="white"),
          axis.ticks = element_line(color="white"),
          axis.text.y = element_text(color="white"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.35, size=9, color="white"),
          plot.background = element_rect(fill="#202020"),
          legend.position = c(0.05,0.9),
          legend.direction = "horizontal",
          legend.box = "horizontal",
          legend.key.width= unit(1, 'cm'))
  message(".")
  
  FIGURE_TWO <- ggdraw()+
    draw_plot(main_temp)+
    draw_plot(beta_env_temp, x=0.05, y=0.65, height=0.2, width=0.2)+
    draw_plot(beta_trait_temp, x=0.25, y=0.65, height=0.2, width=0.2)+
    draw_plot(allResponseCurves, x=0.4, y=0.1, height=0.3, width=0.2)+
    draw_plot(coldResponseCurves, x=0.6, y=0.1, height=0.3, width=0.2)+
    draw_plot(warmResponseCurves, x=0.8, y=0.1, height=0.3, width=0.2)+
    draw_plot_label(label="(a)", x=0, y=1, size=16, color="#FFFFFF")+
    draw_plot_label(label=c("(b)", "(c)"), 
                    x=c(0.08, 0.28), 
                    y=c(0.85, 0.85), size=16, color="#FFFFFF")+
    draw_plot_label(label=c("(d)", "(e)", "(f)"), 
                    x=c(0.42, 0.62, 0.82),
                    y=c(0.425, 0.425, 0.425), size=16, color="#FFFFFF")
  
  
  ggsave2("../output/main/final/FIGURE_TWO.png", FIGURE_TWO, dpi=350, height=10, width=16)
}

# Create a function to visualize figure three in the manuscript
visFigureThree <- function(my_res, my_data, my_traits){
  
  SPIDS <- my_traits$SPID
  
  # Decompress the R2JAGS output into a single object for summaries and plotting
  message("Reconfiguring parallel output as single JAGS run...")
  my_res$chains <- mcmc.list(as.mcmc(my_res[[1]]$BUGSoutput$sims.matrix),
                             as.mcmc(my_res[[2]]$BUGSoutput$sims.matrix),
                             as.mcmc(my_res[[3]]$BUGSoutput$sims.matrix))
  my_res$summary <- MCMCvis::MCMCsummary(my_res$chains)
  my_res$sims.matrix <- as.matrix(my_res$chains)
  my_sims_mat <- my_res$sims.matrix

  message("Calculating species/grid/occupancy interval probabilities...")
  aggregate_occupancy <- list()
  for(sp in 1:length(SPIDS)){
    
    SPID <- SPIDS[sp]
    
    my_sites <- my_data$my.info$range.list[SPID,] %>% na.omit()
    nsites[sp] <- length(my_sites)
    
    if(length(my_sites)==0){
      sp=sp+1
      SPID <- SPIDS[sp]
      my_sites <- my_data$my.info$range.list[SPID,] %>% na.omit()
    }
    
    for(era in 1:my_data$my.constants$nyr){
      for(site in 1:length(my_sites)){
        
        chains <-  plogis(
          my_sims_mat[,sprintf("sp.psi.0[%d]", SPID)]+
            my_sims_mat[,"p.area"]*my_data$my.constants$gridarea[my_sites[site]]+
            my_sims_mat[,sprintf("psi.beta.temp[%d]", SPID)]*my_data$my.constants$temp[site,era]+
            my_sims_mat[,sprintf("psi.beta.precip[%d]", SPID)]*my_data$my.constants$precip[site,era]+
            my_sims_mat[,sprintf("psi.site[%d]", site)]+
            my_sims_mat[,sprintf("psi.yr[%d]", era)])
        
        aggregate_occupancy <- aggregate_occupancy %>%
          append(list(data.frame(SPID=SPID, era=era, site=my_sites[site],
                                 occ=mean(chains), lwr_bound=quantile(chains, probs=c(0.025,0.975))[1],
                                 upr_bound=quantile(chains, probs=c(0.025,0.975))[2])))
      }
    }
  }
  aggregate_occupancy <- do.call(rbind, aggregate_occupancy)
  overall_occupancy <- aggregate_occupancy %>%
    group_by(era, SPID) %>%
    dplyr::mutate(mean_occ=mean(occ), lwr=mean(lwr_bound), upr=mean(upr_bound)) %>%
    ungroup() %>%
    dplyr::select(SPID, era, mean_occ, lwr, upr) %>% unique()
  
  ggplot()
    
}

visSupplementalFigureTwo <- function(my_res, my_data, my_traits){
  # Create species codes for plotting
  my_traits <- my_traits %>%
    dplyr::mutate(specificEpithet=word(species, 2)) %>%
    dplyr::mutate(speciesCode=toupper(paste0(substr(species, 1, 3),
                                             substr(specificEpithet, 1, 3)))) %>%
    dplyr::mutate(precipClass=cut(ave_precip2, 4)) %>%
    dplyr::mutate(precipClass=ifelse(precipClass==levels(precipClass)[1], "Driest Ranges",
                                   ifelse(precipClass==levels(precipClass)[4], "Wettest Ranges", "All Other Ranges")))
  
  # Decompress the R2JAGS output into a single object for summaries and plotting
  message("Reconfiguring parallel output as single JAGS run...")
  my_res$chains <- mcmc.list(as.mcmc(my_res[[1]]$BUGSoutput$sims.matrix),
                             as.mcmc(my_res[[2]]$BUGSoutput$sims.matrix),
                             as.mcmc(my_res[[3]]$BUGSoutput$sims.matrix))
  my_res$summary <- MCMCvis::MCMCsummary(my_res$chains)
  my_res$sims.matrix <- as.matrix(my_res$chains)
  my_sims_mat <- my_res$sims.matrix
  
  # Pull the model summary
  my_sum <- my_res$summary %>% as.data.frame() %>%
    rownames_to_column("Parameter")
  message(".")
  
  # Pull all the samples
  my_samp <- my_res$sims.matrix
  my_samp_beta_precip <- data.frame(samples=my_samp[,"psi.precip.env"])
  my_samp_beta_precip_trait <- data.frame(samples=my_samp[,"psi.precip.trait"])
  message(".")
  
  # Pull the precip effects
  my_precip <- dplyr::filter(my_sum, grepl("psi.beta.precip", Parameter)) %>%
    cbind(my_traits) %>%
    arrange(mean) %>%
    dplyr::mutate(ordering=row_number()) 
  message(".")
  
  # Visualize the precip effects as a forest plot
  message("Visualizing posterior densities for effect subcomponents...")
  beta_env_precip <- ggplot()+
    geom_vline(xintercept=0, linetype=2, color="white")+
    geom_density(my_samp_beta_precip, mapping=aes(x=samples), color="white", fill="white", alpha=0.2)+
    annotate("text", x=Inf, y=Inf, parse=TRUE,
             label="beta[1]", color="white", fontface="bold", size=6, hjust=1, vjust=1)+
    scale_x_continuous(limits=c(-0.1,0.1))+
    labs(x="Effect Estimate", y="Posterior Density")+
    theme_cowplot()+
    theme(text = element_text(color="white"),
          axis.line = element_line(color="white"),
          axis.ticks = element_line(color="white"),
          axis.text.y = element_text(color="white"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.35, size=9, color="white"),
          plot.background = element_rect(fill="#202020"),
          legend.position = c(0.05,0.9),
          legend.direction = "horizontal",
          legend.box = "horizontal",
          legend.key.width= unit(1, 'cm'))
  
  beta_trait_precip <- ggplot()+
    geom_vline(xintercept=0, linetype=2, color="white")+
    geom_density(my_samp_beta_precip_trait, mapping=aes(x=samples), color="white", fill="white", alpha=0.2)+
    annotate("text", x=Inf, y=Inf, parse=TRUE,
             label="beta[2]", color="white", fontface="bold", size=6, hjust=1, vjust=1)+
    scale_x_continuous(limits=c(-0.1,0.1))+
    labs(x="Effect Estimate", y="")+
    theme_cowplot()+
    theme(text = element_text(color="white"),
          axis.line = element_line(color="white"),
          axis.ticks = element_line(color="white"),
          axis.text.y = element_text(color="white"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.35, size=9, color="white"),
          plot.background = element_rect(fill="#202020"),
          legend.position = c(0.05,0.9),
          legend.direction = "horizontal",
          legend.box = "horizontal",
          legend.key.width= unit(1, 'cm'))
  message(".")
  
  # RESPONSE CURVES FOR ALL SPECIES
  message("Plotting response curves for all species...")
  spp_lines <- list()
  for(ss in 1:nrow(my_traits)){
    get.y.val <- function(dd){
      chains <- plogis(
        my_sims_mat[,sprintf("sp.psi.0[%d]", ss)]+
          my_sims_mat[,"p.area"]*mean(my_data$my.constants$gridarea)+
          my_sims_mat[,sprintf("psi.beta.precip[%d]", ss)]*dd+
          (my_sims_mat[,"psi.temp.env"]+
             my_sims_mat[,"psi.temp.trait"]*mean(my_data$my.constants$temp)))
      c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
    }
    x.vals <- seq(from=0, to=300, length=300)
    
    y.vals <- lapply(x.vals, get.y.val)
    y.vals2_spp <- do.call(rbind, y.vals)
    
    y.vals2_spp <- y.vals2_spp %>% as.data.frame() %>%
      dplyr::mutate(spp=ss, x=x.vals)
    spp_lines[[ss]] <- y.vals2_spp
  }
  y.vals2_spp <- do.call(rbind, spp_lines)
  y.vals2_spp <- y.vals2_spp %>% left_join(my_traits, by=c("spp"="SPID"))
  
  get.y.val <- function(dd){
    chains <- plogis(mean(my_sims_mat[,grepl("sp.psi.0", colnames(my_sims_mat))])+
                       my_sims_mat[,"p.area"]*mean(my_data$my.constants$gridarea)+
                       (my_sims_mat[,"psi.precip.env"]+
                          my_sims_mat[,"psi.precip.trait"])*dd+
                       (my_sims_mat[,"psi.temp.env"]+
                          my_sims_mat[,"psi.temp.trait"]*mean(my_data$my.constants$temp)))
    
    c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
  }
  x.vals <- seq(from=0, to=300, length=300)
  y.vals <- lapply(x.vals, get.y.val)
  y.vals2_com <- do.call(rbind, y.vals) %>% as.data.frame()
  y.vals2_com$x <- x.vals
  
  allResponseCurves <- ggplot()+
    geom_line(y.vals2_spp, mapping=aes(x=x, y=mean, group=spp), 
              alpha=0.05, color="white")+
    labs(x="", y="Occupancy Probability")+
    geom_line(y.vals2_com, mapping=aes(x=x, y=mean), color="white")+
    geom_line(y.vals2_com, mapping=aes(x=x, y=`2.5%`), color="white", linetype=2)+
    geom_line(y.vals2_com, mapping=aes(x=x, y=`97.5%`), color="white", linetype=2)+
    scale_y_continuous(labels=scales::percent)+
    theme_cowplot()+
    theme(text = element_text(color="white"),
          axis.line = element_line(color="white"),
          axis.ticks = element_line(color="white"),
          axis.text.y = element_text(color="white"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.35, size=9, color="white"),
          plot.background = element_rect(fill="#202020"),
          legend.position = c(0.05,0.9),
          legend.direction = "horizontal",
          legend.box = "horizontal",
          legend.key.width= unit(1, 'cm'))
  message(".")
  
  # RESPONSE CURVES FOR driest ADAPTED SPECIES  
  message("Plotting response curves for driest-adapted species...")
  spp_lines <- list()
  driest_species <- dplyr::filter(my_traits, precipClass=="Driest Ranges")
  for(ss in 1:nrow(driest_species)){
    get.y.val <- function(dd){
      chains <- plogis(
        my_sims_mat[,sprintf("sp.psi.0[%d]", driest_species$SPID[ss])]+
          my_sims_mat[,"p.area"]*mean(my_data$my.constants$gridarea)+
          my_sims_mat[,sprintf("psi.beta.precip[%d]", driest_species$SPID[ss])]*dd+
          my_sims_mat[,sprintf("psi.beta.temp[%d]", driest_species$SPID[ss])]*
          mean(my_data$my.constants$temp))
      c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
    }
    x.vals <- seq(from=0, to=300, length=300)
    
    y.vals <- lapply(x.vals, get.y.val)
    y.vals2_spp <- do.call(rbind, y.vals)
    
    y.vals2_spp <- y.vals2_spp %>% as.data.frame() %>%
      dplyr::mutate(spp=ss, x=x.vals)
    spp_lines[[ss]] <- y.vals2_spp
  }
  y.vals2_spp <- do.call(rbind, spp_lines)
  y.vals2_spp <- y.vals2_spp %>% left_join(my_traits, by=c("spp"="SPID"))
  
  get.y.val <- function(dd){
    chains <- plogis(mean(my_sims_mat[,sprintf("sp.psi.0[%d]", driest_species$SPID)])+
                       my_sims_mat[,"p.area"]*mean(my_data$my.constants$gridarea)+
                       mean(my_sims_mat[,sprintf("psi.beta.precip[%d]", driest_species$SPID)])*dd+
                       mean(my_sims_mat[,sprintf("psi.beta.temp[%d]", driest_species$SPID)])*
                       mean(my_data$my.constants$temp))
    
    c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
  }
  x.vals <- seq(from=0, to=300, length=300)
  y.vals <- lapply(x.vals, get.y.val)
  y.vals2_com <- do.call(rbind, y.vals) %>% as.data.frame()
  y.vals2_com$x <- x.vals
  
  dryResponseCurves <- ggplot()+
    geom_line(y.vals2_spp, mapping=aes(x=x, y=mean, group=spp), 
              alpha=0.05, color="#edf8b1")+
    labs(x="", y="")+
    geom_line(y.vals2_com, mapping=aes(x=x, y=mean), color="#edf8b1")+
    geom_line(y.vals2_com, mapping=aes(x=x, y=`2.5%`), color="#edf8b1", linetype=2)+
    geom_line(y.vals2_com, mapping=aes(x=x, y=`97.5%`), color="#edf8b1", linetype=2)+
    scale_y_continuous(labels=scales::percent)+
    theme_cowplot()+
    theme(text = element_text(color="white"),
          axis.line = element_line(color="white"),
          axis.ticks = element_line(color="white"),
          axis.text.y = element_text(color="white"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.35, size=9, color="white"),
          plot.background = element_rect(fill="#202020"),
          legend.position = c(0.05,0.9),
          legend.direction = "horizontal",
          legend.box = "horizontal",
          legend.key.width= unit(1, 'cm'))
  message(".")
  
  # RESPONSE CURVES FOR WETTEST ADAPTED SPECIES
  message("Plotting response curves for wetest-adapted species...")
  spp_lines <- list()
  wettest_species <- dplyr::filter(my_traits, precipClass=="Wettest Ranges")
  for(ss in 1:nrow(wettest_species)){
    get.y.val <- function(dd){
      chains <- plogis(
        my_sims_mat[,sprintf("sp.psi.0[%d]", wettest_species$SPID[ss])]+
          my_sims_mat[,"p.area"]*mean(my_data$my.constants$gridarea)+
          my_sims_mat[,sprintf("psi.beta.precip[%d]", wettest_species$SPID[ss])]*dd+
          mean(my_sims_mat[,sprintf("psi.beta.temp[%d]", wettest_species$SPID[ss])])*
          mean(my_data$my.constants$temp))
      c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
    }
    x.vals <- seq(from=0, to=300, length=300)
    
    y.vals <- lapply(x.vals, get.y.val)
    y.vals2_spp <- do.call(rbind, y.vals)
    
    y.vals2_spp <- y.vals2_spp %>% as.data.frame() %>%
      dplyr::mutate(spp=ss, x=x.vals)
    spp_lines[[ss]] <- y.vals2_spp
  }
  y.vals2_spp <- do.call(rbind, spp_lines)
  y.vals2_spp <- y.vals2_spp %>% left_join(my_traits, by=c("spp"="SPID"))
  
  get.y.val <- function(dd){
    chains <- plogis(mean(my_sims_mat[,sprintf("sp.psi.0[%d]", wettest_species$SPID)])+
                       my_sims_mat[,"p.area"]*mean(my_data$my.constants$gridarea)+
                       mean(my_sims_mat[,sprintf("psi.beta.precip[%d]", wettest_species$SPID)])*dd+
                       mean(my_sims_mat[,sprintf("psi.beta.temp[%d]", wettest_species$SPID)])*
                       mean(my_data$my.constants$temp))
    
    c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
  }
  x.vals <- seq(from=0, to=300, length=300)
  y.vals <- lapply(x.vals, get.y.val)
  y.vals2_com <- do.call(rbind, y.vals) %>% as.data.frame()
  y.vals2_com$x <- x.vals
  
  wetResponseCurves <- ggplot()+
    geom_line(y.vals2_spp, mapping=aes(x=x, y=mean, group=spp), 
              alpha=0.05, color="#1d91c0")+
    labs(x="", y="")+
    geom_line(y.vals2_com, mapping=aes(x=x, y=mean), color="#1d91c0")+
    geom_line(y.vals2_com, mapping=aes(x=x, y=`2.5%`), color="#1d91c0", linetype=2)+
    geom_line(y.vals2_com, mapping=aes(x=x, y=`97.5%`), color="#1d91c0", linetype=2)+
    scale_y_continuous(labels=scales::percent)+
    theme_cowplot()+
    theme(text = element_text(color="white"),
          axis.line = element_line(color="white"),
          axis.ticks = element_line(color="white"),
          axis.text.y = element_text(color="white"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.35, size=9, color="white"),
          plot.background = element_rect(fill="#202020"),
          legend.position = c(0.05,0.9),
          legend.direction = "horizontal",
          legend.box = "horizontal",
          legend.key.width= unit(1, 'cm'))
  message(".")
  
  # CONSTRUCT THE MAIN FOREST PLOT
  message("Plotting the main forest plot of effect estimates...")
  main_precip <- ggplot()+
    geom_hline(yintercept=0, linetype=2, color="white")+
    geom_pointrange(my_precip, 
                    mapping=aes(x=ordering,
                                y=mean, 
                                ymin=`2.5%`, 
                                ymax=`97.5%`, 
                                group=SPID, color=ave_precip2, shape=precipClass))+
    annotate("text", x=17, y=1.75, parse=TRUE, hjust=0, vjust=1,
             label="Effect == beta[1] + beta[2]*PrecipitationTrait", 
             size=6, color="white")+
    scale_x_continuous(breaks=c(1:nrow(my_traits)), labels=my_precip$speciesCode)+
    scale_y_continuous(limits=c(-2,2))+
    scale_shape_manual(values=c(1, 15, 19), guide=guide_legend(override.aes=list(fill="white", color="white")),
                       name="Range Wide preciperature Class:")+
    scale_color_distiller(palette="YlGnBu", direction=1,
                                      guide=guide_colorbar(title.position = "top"),
                                      name="Average Annual Range-wide Precipitation (mm)")+
    labs(x="Species Code", y="Estimated Effect of Precipitation on Occupancy Probability")+
    theme_cowplot()+
    theme(text = element_text(color="white"),
          axis.line = element_line(color="white"),
          axis.ticks = element_line(color="white"),
          axis.text.y = element_text(color="white"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.35, size=9, color="white"),
          plot.background = element_rect(fill="#202020"),
          legend.position = c(0.05,0.9),
          legend.direction = "horizontal",
          legend.box = "horizontal",
          legend.key.width= unit(1, 'cm'))
  message(".")
  
  FIGURE_SUPP_TWO <- ggdraw()+
    draw_plot(main_precip)+
    draw_plot(beta_env_precip, x=0.05, y=0.65, height=0.2, width=0.2)+
    draw_plot(beta_trait_precip, x=0.25, y=0.65, height=0.2, width=0.2)+
    draw_plot(allResponseCurves, x=0.4, y=0.1, height=0.3, width=0.2)+
    draw_plot(dryResponseCurves, x=0.6, y=0.1, height=0.3, width=0.2)+
    draw_plot(wetResponseCurves, x=0.8, y=0.1, height=0.3, width=0.2)+
    draw_plot_label(label="(a)", x=0, y=1, size=16, color="#FFFFFF")+
    draw_plot_label(label=c("(b)", "(c)"), 
                    x=c(0.08, 0.28), 
                    y=c(0.85, 0.85), size=16, color="#FFFFFF")+
    draw_plot_label(label=c("(d)", "(e)", "(f)"), 
                    x=c(0.42, 0.62, 0.82),
                    y=c(0.425, 0.425, 0.425), size=16, color="#FFFFFF")
  
  
  ggsave2("../output/main/supp/FIGURE_S2.png", FIGURE_SUPP_TWO, dpi=400, height=10, width=16)
  ggsave2("../output/main/supp/FIGURE_S2.pdf", FIGURE_SUPP_TWO, dpi=400, height=10, width=16)
  
}



