# FUNCTION: species_specific_occupancy - computes the occupancy estimate for
# each species in the analysis for all time periods. 
species_specific_occupancy <- function(sims.mat, modelData){
  
  sp_occ_site_trends <- list()
  
  my_grid_200 <- readRDS("../../output/data/grid_200.rds") %>%
    dplyr::filter(GID %in% modelData$my.info$kept.sites)
  
  ranges <- readRDS("../../output/data/grid_200_range.rds") %>%
    as.data.frame() %>%
    dplyr::mutate(GID=row_number()) %>%
    as.data.frame() %>%
    left_join(my_grid_200, by="GID") %>%
    dplyr::filter(GID %in% modelData$my.info$kept.sites) %>%
    sf::st_drop_geometry() 
  
  for(ss in 1:modelData$my.constants$nsp){
    
    relevant_sites <- ranges[,"GID"][ranges[,ss]]
    relevant_sites <- relevant_sites[relevant_sites %in% modelData$my.info$kept.sites]
    
    occChains <- list()
    
    if(covariate=="temp"){
      for(site in 1:length(relevant_sites)){
        
        occChains[[site]] <- c(ss, relevant_sites[site], mean(c(plogis(sims.mat[,"mu.psi.0"]+
                                                                         sims.mat[,sprintf("psi.sp[%d]", ss)]+
                                                                         sims.mat[,"psi.area"]*modelData$my.data$gridarea[site]+
                                                                         sims.mat[,sprintf("psi.beta.temp[%d]", ss)]*modelData$my.data$temp[site,10]+
                                                                         sims.mat[,sprintf("psi.beta.temp[%d]", ss)]*modelData$my.data$temp[site,10]^2),
                                                                plogis(sims.mat[,"mu.psi.0"]+
                                                                         sims.mat[,sprintf("psi.sp[%d]", ss)]+
                                                                         sims.mat[,"psi.area"]*modelData$my.data$gridarea[site]+
                                                                         sims.mat[,sprintf("psi.beta.temp[%d]", ss)]*modelData$my.data$temp[site,9]+
                                                                         sims.mat[,sprintf("psi.beta.temp[%d]", ss)]*modelData$my.data$temp[site,9]^2)))-
                                 mean(c(plogis(sims.mat[,"mu.psi.0"]+
                                                 sims.mat[,sprintf("psi.sp[%d]", ss)]+
                                                 sims.mat[,"psi.area"]*modelData$my.data$gridarea[site]+
                                                 sims.mat[,sprintf("psi.beta.temp[%d]", ss)]*modelData$my.data$temp[site,1]+
                                                 sims.mat[,sprintf("psi.beta.temp[%d]", ss)]*modelData$my.data$temp[site,1]^2),
                                        plogis(sims.mat[,"mu.psi.0"]+
                                                 sims.mat[,sprintf("psi.sp[%d]", ss)]+
                                                 sims.mat[,"psi.area"]*modelData$my.data$gridarea[site]+
                                                 sims.mat[,sprintf("psi.beta.temp[%d]", ss)]*modelData$my.data$temp[site,2]+
                                                 sims.mat[,sprintf("psi.beta.temp[%d]", ss)]*modelData$my.data$temp[site,2]^2))))
        
      }
      sp_occ_site_trends[[ss]] <- do.call(rbind, occChains) 
    } else if(covariate=="precip"){
      for(site in 1:length(relevant_sites)){
        
        occChains[[site]] <- c(ss, relevant_sites[site], mean(c(plogis(sims.mat[,"mu.psi.0"]+
                                                                         sims.mat[,sprintf("psi.sp[%d]", ss)]+
                                                                         sims.mat[,"psi.area"]*modelData$my.data$gridarea[site]+
                                                                         sims.mat[,sprintf("psi.beta.precip[%d]", ss)]*modelData$my.data$precip[site,10]),
                                                                plogis(sims.mat[,"mu.psi.0"]+
                                                                         sims.mat[,sprintf("psi.sp[%d]", ss)]+
                                                                         sims.mat[,"psi.area"]*modelData$my.data$gridarea[site]+
                                                                         sims.mat[,sprintf("psi.beta.precip[%d]", ss)]*modelData$my.data$precip[site,9])))-
                                 mean(c(plogis(sims.mat[,"mu.psi.0"]+
                                                 sims.mat[,sprintf("psi.sp[%d]", ss)]+
                                                 sims.mat[,"psi.area"]*modelData$my.data$gridarea[site]+
                                                 sims.mat[,sprintf("psi.beta.precip[%d]", ss)]*modelData$my.data$precip[site,1]),
                                        plogis(sims.mat[,"mu.psi.0"]+
                                                 sims.mat[,sprintf("psi.sp[%d]", ss)]+
                                                 sims.mat[,"psi.area"]*modelData$my.data$gridarea[site]+
                                                 sims.mat[,sprintf("psi.beta.precip[%d]", ss)]*modelData$my.data$precip[site,2]))))
        
      }
      sp_occ_site_trends[[ss]] <- do.call(rbind, occChains) 
    }
  }
  sp_trends <- do.call(rbind, sp_occ_site_trends) %>% as.data.frame()
  colnames(sp_trends) <- c("SPID", "GID", "Ratio")
  
  y_vals_occ <- sp_trends %>%
    dplyr::group_by(SPID) %>%
    dplyr::mutate(nsite=n()) %>%
    dplyr::mutate(mean=mean(Ratio),
                  lower=quantile(Ratio, probs=0.25),
                  upper=quantile(Ratio, probs=0.75)) %>%
    dplyr::select(SPID, mean, lower, upper) %>%
    dplyr::ungroup() %>%
    unique() %>%
    dplyr::arrange(mean) %>%
    dplyr::mutate(ordering=row_number()) %>%
    left_join(sp_traits) %>%
    dplyr::mutate(code=stringr::str_to_upper(paste0(substr(species, 1, 3),
                                                    substr(word(species, 2), 1, 3))))
  
  return(list(y_vals_occ, trend_plot))
  
}

# FUNCTION: plot_response_curves - plots the species-specific and overall
# butterfly response to a given environmental predictor for values that are
# relevant over the species range (absolute min/max for each species).
plot_response_curves <- function(sims.mat, modelData, sp_traits){
  
  # Calculate the overall butterfly response to average minimum temperature
  x_temp <- seq(-2, 6, length.out=100)
  getYVal_temp_com <- function(dd, sims.mat){
    chains <- plogis(sims.mat[,"mu.psi.0"]+
                       rowMeans(sims.mat[,grepl("psi.sp", colnames(sims.mat))])+
                       rowMeans(sims.mat[,grepl("psi.beta.precip", colnames(sims.mat))])*dd)
    
    y_vals <- c(mean=mean(chains), quantile(chains, probs=c(0.05,0.95)))
    names(y_vals) <- c("mean", "lower", "upper")
    return(y_vals)
  } 
  y_vals_com <- lapply(x_temp, getYVal_temp_com, sims.mat)
  y_vals_com <- do.call(rbind, y_vals_com) %>% as.data.frame()
  
  # Calculate the overall cold-adapted butterfly response to average minimum temperature
  # for the coldest quarter of species
  cold_indices <- which(sp_traits$rangePrecip <= quantile(sp_traits$rangePrecip, probs=0.25))
  getYVal_temp_cold <- function(dd, sims.mat){
    chains <- plogis(sims.mat[,"mu.psi.0"]+
                       rowMeans(sims.mat[,grepl("psi.sp", colnames(sims.mat))][,cold_indices])+
                       rowMeans(sims.mat[,grepl("psi.beta.precip", colnames(sims.mat))][,cold_indices])*dd)
    
    y_vals <- c(mean=mean(chains), quantile(chains, probs=c(0.05,0.95)))
    names(y_vals) <- c("mean", "lower", "upper")
    return(y_vals)
  } 
  y_vals_cold <- lapply(x_temp, getYVal_temp_cold, sims.mat)
  y_vals_cold <- do.call(rbind, y_vals_cold) %>% as.data.frame()
  
  # Calculate the overall warm-adapted butterfly response to average minimum temperature
  # for the warmest quarter of species
  warm_indices <- which(sp_traits$rangePrecip >= quantile(sp_traits$rangePrecip, probs=0.75))
  getYVal_temp_warm <- function(dd, sims.mat){
    chains <- plogis(sims.mat[,"mu.psi.0"]+
                       rowMeans(sims.mat[,grepl("psi.sp", colnames(sims.mat))][,warm_indices])+
                       rowMeans(sims.mat[,grepl("psi.beta.precip", colnames(sims.mat))][,warm_indices])*dd)
    
    y_vals <- c(mean=mean(chains), quantile(chains, probs=c(0.05,0.95)))
    names(y_vals) <- c("mean", "lower", "upper")
    return(y_vals)
  } 
  y_vals_warm <- lapply(x_temp, getYVal_temp_warm, sims.mat)
  y_vals_warm <- do.call(rbind, y_vals_warm) %>% as.data.frame()
  
  # Calculate the overall butterfly response again, but keep the chains intact
  getYVal_temp_com_chains <- function(dd, sims.mat){
    chains <- plogis(sims.mat[,"mu.psi.0"]+
                       rowMeans(sims.mat[,grepl("psi.sp", colnames(sims.mat))])+
                       rowMeans(sims.mat[,grepl("psi.beta.precip", colnames(sims.mat))])*dd)
    
    return(data.frame(dd, chains))
  } 
  y_vals_com_chains <- lapply(x_temp, getYVal_temp_com_chains, sims.mat)
  y_vals_com_chains <- do.call(rbind, y_vals_com_chains) %>% as.data.frame()
  
  # Calculate the species-specific responses to average minimum temperature.
  # Grab the minimum and maximum temperature experienced by the species within
  # its range.
  spp_temp_range <- list()
  for(ss in 1:modelData$my.constants$nsp){
    spp_temp_range[[ss]] <- c(min=min(modelData$my.data$precip[modelData$my.info$range.list[ss,]], na.rm=TRUE),
                              max=max(modelData$my.data$precip[modelData$my.info$range.list[ss,]], na.rm=TRUE))
  }
  
  spp_lines <- list()
  for(ss in 1:modelData$my.constants$nsp){
    get.y.val <- function(dd){
      chains <- plogis(sims.mat[,"mu.psi.0"]+
                         sims.mat[,sprintf("psi.sp[%d]", ss)]+
                         sims.mat[,sprintf("psi.beta.precip[%d]", ss)]*dd)
      
      c(mean=mean(chains), quantile(chains, probs=c(0.05,0.095)))
    }
    y_vals <- lapply(x_temp, get.y.val)
    y_vals2_spp <- do.call(rbind, y_vals)
    
    y_vals2_spp <- y_vals2_spp %>% as.data.frame() %>%
      dplyr::mutate(spp=ss, x=x_temp) %>%
      dplyr::mutate(mean=ifelse(between(x, 
                                        spp_temp_range[[ss]][1], 
                                        spp_temp_range[[ss]][2]), 
                                mean, 
                                NA)) %>%
      dplyr::filter(!is.na(mean))
    spp_lines[[ss]] <- y_vals2_spp
  }
  y_vals2_spp <- do.call(rbind, spp_lines) %>%
    left_join(sp_traits, by=c("spp"="SPID"))
  
  response_curves <- ggplot()+
    geom_line(y_vals2_spp,
              mapping=aes(x=x, y=mean, 
                          group=spp),
              alpha=0.1)+
    # OVERALL AVERAGE RESPONSE FROM ALL BUTTERFLIES
    geom_line(y_vals_com,
              mapping=aes(x=x_temp, y=mean),
              size=1)+
    geom_ribbon(y_vals_com,
                mapping=aes(x=x_temp, ymin=lower, ymax=upper),
                alpha=0.4)+
    # COLD-ADAPTED SPECIES RESPONSES
    geom_line(y_vals_cold,
              mapping=aes(x=x_temp, y=mean),
              size=1, color="#a16928")+
    geom_ribbon(y_vals_cold,
                mapping=aes(x=x_temp, ymin=lower, ymax=upper),
                alpha=0.4, color=NA, fill="#a16928")+
    # WARM-ADAPTED SPECIES RESPONSES
    geom_line(y_vals_warm,
              mapping=aes(x=x_temp, y=mean),
              size=1, color="#2887a1")+
    geom_ribbon(y_vals_warm,
                mapping=aes(x=x_temp, ymin=lower, ymax=upper),
                alpha=0.4, color=NA, fill="#2887a1")+
    scale_y_continuous(limits=c(0,1), labels=scales::percent)+
    scale_x_continuous(limits=c(-1.223907,6), breaks=c(-1.223907, 0, 2, 4, 6),
                       labels=round((c(-1.223907, 0, 2, 4, 6)*sd(my_data_200_1$master.data$precip, na.rm=TRUE))+
                                      mean(my_data_200_1$master.data$precip, na.rm=TRUE)))+
    labs(x="Average Precipitation (mm)",
         y="Occupancy Probability")+
    theme_cowplot()+
    theme(plot.background=element_rect(color=NULL, fill="white"),
          axis.title=element_text(size=18),
          axis.text=element_text(size=16),
          legend.position="none")
  ggsave2("../../figures/supplemental/precip_responseCurve_200.png", dpi=400, height=8, width=8)
  return(response_curves)
  
}

# FUNCTION: grid_species_specific_occupancy - computes the difference in occupancy
# between two occupancy intervals (specified by parameters). This can be used to
# calculate slope fields to assess the directionality of occupancy shift over time.
grid_species_specific_occupancy <- function(sims.mat, modelData, covariate="temp"){
  
  sp_occ_site_trends <- list()
  
  my_grid_200 <- readRDS("../../output/data/grid_200.rds") %>%
    dplyr::filter(GID %in% modelData$my.info$kept.sites)
  
  ranges <- readRDS("../../output/data/grid_200_range.rds") %>%
    as.data.frame() %>%
    dplyr::mutate(GID=row_number()) %>%
    as.data.frame() %>%
    left_join(my_grid_200, by="GID") %>%
    dplyr::filter(GID %in% modelData$my.info$kept.sites) %>%
    sf::st_drop_geometry() 
  
  for(ss in 1:modelData$my.constants$nsp){
    
    relevant_sites <- ranges[,"GID"][ranges[,ss]]
    relevant_sites <- relevant_sites[relevant_sites %in% modelData$my.info$kept.sites]
    
    occChains <- list()
    if(covariate=="temp"){
      for(site in 1:length(relevant_sites)){
        occChains[[site]] <- c(ss, relevant_sites[site], mean(mean(c(plogis(sims.mat[,"mu.psi.0"]+
                                                                              sims.mat[,sprintf("psi.sp[%d]", ss)]+
                                                                              sims.mat[,"psi.area"]*modelData$my.data$gridarea[site]+
                                                                              sims.mat[,sprintf("psi.beta.temp[%d]", ss)]*modelData$my.data$temp[site,10]+
                                                                              sims.mat[,sprintf("psi.beta.temp2[%d]", ss)]*modelData$my.data$temp[site,10]^2),
                                                                     plogis(sims.mat[,"mu.psi.0"]+
                                                                              sims.mat[,sprintf("psi.sp[%d]", ss)]+
                                                                              sims.mat[,"psi.area"]*modelData$my.data$gridarea[site]+
                                                                              sims.mat[,sprintf("psi.beta.temp[%d]", ss)]*modelData$my.data$temp[site,9]+
                                                                              sims.mat[,sprintf("psi.beta.temp2[%d]", ss)]*modelData$my.data$temp[site,9]^2)))-
                                                                mean(c(plogis(sims.mat[,"mu.psi.0"]+
                                                                                sims.mat[,sprintf("psi.sp[%d]", ss)]+
                                                                                sims.mat[,"psi.area"]*modelData$my.data$gridarea[site]+
                                                                                sims.mat[,sprintf("psi.beta.temp[%d]", ss)]*modelData$my.data$temp[site,1]+
                                                                                sims.mat[,sprintf("psi.beta.temp2[%d]", ss)]*modelData$my.data$temp[site,1]^2),
                                                                       plogis(sims.mat[,"mu.psi.0"]+
                                                                                sims.mat[,sprintf("psi.sp[%d]", ss)]+
                                                                                sims.mat[,"psi.area"]*modelData$my.data$gridarea[site]+
                                                                                sims.mat[,sprintf("psi.beta.temp[%d]", ss)]*modelData$my.data$temp[site,2]+
                                                                                sims.mat[,sprintf("psi.beta.temp2[%d]", ss)]*modelData$my.data$temp[site,2]^2)))))
        
      }
      sp_occ_site_trends[[ss]] <- do.call(rbind, occChains)  
    } else if(covariate=="precip"){
      for(site in 1:length(relevant_sites)){
        occChains[[site]] <- c(ss, relevant_sites[site], mean(mean(c(plogis(sims.mat[,"mu.psi.0"]+
                                                                              sims.mat[,sprintf("psi.sp[%d]", ss)]+
                                                                              sims.mat[,"psi.area"]*modelData$my.data$gridarea[site]+
                                                                              sims.mat[,sprintf("psi.beta.precip[%d]", ss)]*modelData$my.data$precip[site,10]),
                                                                     plogis(sims.mat[,"mu.psi.0"]+
                                                                              sims.mat[,sprintf("psi.sp[%d]", ss)]+
                                                                              sims.mat[,"psi.area"]*modelData$my.data$gridarea[site]+
                                                                              sims.mat[,sprintf("psi.beta.precip[%d]", ss)]*modelData$my.data$precip[site,9])))-
                                                                mean(c(plogis(sims.mat[,"mu.psi.0"]+
                                                                                sims.mat[,sprintf("psi.sp[%d]", ss)]+
                                                                                sims.mat[,"psi.area"]*modelData$my.data$gridarea[site]+
                                                                                sims.mat[,sprintf("psi.beta.precip[%d]", ss)]*modelData$my.data$precip[site,1]),
                                                                       plogis(sims.mat[,"mu.psi.0"]+
                                                                                sims.mat[,sprintf("psi.sp[%d]", ss)]+
                                                                                sims.mat[,"psi.area"]*modelData$my.data$gridarea[site]+
                                                                                sims.mat[,sprintf("psi.beta.precip[%d]", ss)]*modelData$my.data$precip[site,2])))))
        
      }
      sp_occ_site_trends[[ss]] <- do.call(rbind, occChains)  
    }
    
  }
  sp_trends <- do.call(rbind, sp_occ_site_trends) %>% as.data.frame()
  colnames(sp_trends) <- c("SPID", "GID", "Ratio")
  return(sp_trends)
}
