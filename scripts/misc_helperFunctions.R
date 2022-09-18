####################################################################################################
# FUNCTION: compute_occ_shift -
# Computes the occupancy interval shift at a collection of sites given the samples derived
# from occupancy-detection modeling. Returns results at both the 100 x 100 and 200 x 200
# kilometer scales.
####################################################################################################

compute_occ_shift <- function(my_data_100_1, sim_matrix_100, 
                              my_data_200_1, sim_matrix_200,
                              my_traits, site_type="mean"){
  # Main panel figure (overall species occupancy shifts from OI 1 to OI 10)
  sp_occ_diff_100 <- list()
  sp_occ_trend_100 <- list()
  
  sp_occ_diff_200 <- list()
  sp_occ_trend_200 <- list()
  
  for(i in 1:nrow(my_traits)){
    
    sp_sites <- my_data_200_1$my.info$range.list[i,] %>% na.omit()
    
    # AVERAGE SITE IN 1970-1974
    if(site_type=="mean"){
      
      ave_site_tempRange <- quantile(my_data_200_1$my.data$temp[sp_sites, 1])[c(2,4)]
      ave_site_tempInd <- which(between(my_data_200_1$my.data$temp[sp_sites, 1],
                                        ave_site_tempRange[1], ave_site_tempRange[2]), 
                                arr.ind=TRUE)
      
      ave_site_precipRange <- quantile(my_data_200_1$my.data$precip[sp_sites, 1])[c(2,4)]
      ave_site_precipInd <- which(between(my_data_200_1$my.data$precip[sp_sites, 1],
                                          ave_site_precipRange[1], ave_site_precipRange[2]), 
                                  arr.ind=TRUE)
      
      ave_site_Ind <- intersect(ave_site_tempInd, ave_site_precipInd)
      
      ave_site_size <- my_data_200_1$my.data$gridarea[ave_site_Ind]
      
      site_yr_effects_OI_1 <-  sims_matrix_200[,sprintf("psi.site[%d,1]", ave_site_Ind)]
      site_yr_effects_OI_10 <- sims_matrix_200[,sprintf("psi.site[%d,10]", ave_site_Ind)]
      
      ave_site_temp_1 <- my_data_200_1$my.data$temp[ave_site_Ind, 1]
      ave_site_temp_2 <- my_data_200_1$my.data$temp[ave_site_Ind, 2]
      ave_site_temp_3 <- my_data_200_1$my.data$temp[ave_site_Ind, 3]
      ave_site_temp_4 <- my_data_200_1$my.data$temp[ave_site_Ind, 4]
      ave_site_temp_5 <- my_data_200_1$my.data$temp[ave_site_Ind, 5]
      ave_site_temp_6 <- my_data_200_1$my.data$temp[ave_site_Ind, 6]
      ave_site_temp_7 <- my_data_200_1$my.data$temp[ave_site_Ind, 7]
      ave_site_temp_8 <- my_data_200_1$my.data$temp[ave_site_Ind, 8]
      ave_site_temp_9 <- my_data_200_1$my.data$temp[ave_site_Ind, 9]
      ave_site_temp_10 <- my_data_200_1$my.data$temp[ave_site_Ind, 10]
      
      ave_site_precip_1 <- my_data_200_1$my.data$precip[ave_site_Ind, 1]
      ave_site_precip_2 <- my_data_200_1$my.data$precip[ave_site_Ind, 2]
      ave_site_precip_3 <- my_data_200_1$my.data$precip[ave_site_Ind, 3]
      ave_site_precip_4 <- my_data_200_1$my.data$precip[ave_site_Ind, 4]
      ave_site_precip_5 <- my_data_200_1$my.data$precip[ave_site_Ind, 5]
      ave_site_precip_6 <- my_data_200_1$my.data$precip[ave_site_Ind, 6]
      ave_site_precip_7 <- my_data_200_1$my.data$precip[ave_site_Ind, 7]
      ave_site_precip_8 <- my_data_200_1$my.data$precip[ave_site_Ind, 8]
      ave_site_precip_9 <- my_data_200_1$my.data$precip[ave_site_Ind, 9]
      ave_site_precip_10 <- my_data_200_1$my.data$precip[ave_site_Ind, 10]
      
    } else if(site_type=="low_temp_mean_precip"){
      
      ave_site_tempRange <- quantile(my_data_200_1$my.data$temp[sp_sites, 1])[c(1,2)]
      ave_site_tempInd <- which(between(my_data_200_1$my.data$temp[sp_sites, 1],
                                        ave_site_tempRange[1], ave_site_tempRange[2]), 
                                arr.ind=TRUE)
      
      ave_site_precipRange <- quantile(my_data_200_1$my.data$precip[sp_sites, 1])[c(2,4)]
      ave_site_precipInd <- which(between(my_data_200_1$my.data$precip[sp_sites, 1],
                                          ave_site_precipRange[1], ave_site_precipRange[2]), 
                                  arr.ind=TRUE)
      
      if(length(intersect(ave_site_tempInd, ave_site_precipInd))==0){
        ave_site_temp_1 <- my_data_200_1$my.data$temp[ave_site_tempInd, 1]
        ave_site_temp_2 <- my_data_200_1$my.data$temp[ave_site_tempInd, 2]
        ave_site_temp_3 <- my_data_200_1$my.data$temp[ave_site_tempInd, 3]
        ave_site_temp_4 <- my_data_200_1$my.data$temp[ave_site_tempInd, 4]
        ave_site_temp_5 <- my_data_200_1$my.data$temp[ave_site_tempInd, 5]
        ave_site_temp_6 <- my_data_200_1$my.data$temp[ave_site_tempInd, 6]
        ave_site_temp_7 <- my_data_200_1$my.data$temp[ave_site_tempInd, 7]
        ave_site_temp_8 <- my_data_200_1$my.data$temp[ave_site_tempInd, 8]
        ave_site_temp_9 <- my_data_200_1$my.data$temp[ave_site_tempInd, 9]
        ave_site_temp_10 <- my_data_200_1$my.data$temp[ave_site_tempInd, 10]
        
        ave_site_precip_1 <- my_data_200_1$my.data$precip[sp_sites, 1]
        ave_site_precip_2 <- my_data_200_1$my.data$precip[sp_sites, 2]
        ave_site_precip_3 <- my_data_200_1$my.data$precip[sp_sites, 3]
        ave_site_precip_4 <- my_data_200_1$my.data$precip[sp_sites, 4]
        ave_site_precip_5 <- my_data_200_1$my.data$precip[sp_sites, 5]
        ave_site_precip_6 <- my_data_200_1$my.data$precip[sp_sites, 6]
        ave_site_precip_7 <- my_data_200_1$my.data$precip[sp_sites, 7]
        ave_site_precip_8 <- my_data_200_1$my.data$precip[sp_sites, 8]
        ave_site_precip_9 <- my_data_200_1$my.data$precip[sp_sites, 9]
        ave_site_precip_10 <- my_data_200_1$my.data$precip[sp_sites, 10]
        
        ave_site_size <- my_data_200_1$my.data$gridarea[ave_site_tempInd]
        
        site_yr_effects_OI_1 <-  sims_matrix_200[,sprintf("psi.site[%d,1]", ave_site_tempInd)]
        site_yr_effects_OI_10 <- sims_matrix_200[,sprintf("psi.site[%d,10]", ave_site_tempInd)]
        
        ave_site_Ind <- ave_site_tempInd
        
      } else{
        ave_site_Ind <- intersect(ave_site_tempInd, ave_site_precipInd)
        
        ave_site_temp_1 <- my_data_200_1$my.data$temp[ave_site_Ind, 1]
        ave_site_temp_2 <- my_data_200_1$my.data$temp[ave_site_Ind, 2]
        ave_site_temp_3 <- my_data_200_1$my.data$temp[ave_site_Ind, 3]
        ave_site_temp_4 <- my_data_200_1$my.data$temp[ave_site_Ind, 4]
        ave_site_temp_5 <- my_data_200_1$my.data$temp[ave_site_Ind, 5]
        ave_site_temp_6 <- my_data_200_1$my.data$temp[ave_site_Ind, 6]
        ave_site_temp_7 <- my_data_200_1$my.data$temp[ave_site_Ind, 7]
        ave_site_temp_8 <- my_data_200_1$my.data$temp[ave_site_Ind, 8]
        ave_site_temp_9 <- my_data_200_1$my.data$temp[ave_site_Ind, 9]
        ave_site_temp_10 <- my_data_200_1$my.data$temp[ave_site_Ind, 10]
        
        ave_site_precip_1 <- my_data_200_1$my.data$precip[ave_site_Ind, 1]
        ave_site_precip_2 <- my_data_200_1$my.data$precip[ave_site_Ind, 2]
        ave_site_precip_3 <- my_data_200_1$my.data$precip[ave_site_Ind, 3]
        ave_site_precip_4 <- my_data_200_1$my.data$precip[ave_site_Ind, 4]
        ave_site_precip_5 <- my_data_200_1$my.data$precip[ave_site_Ind, 5]
        ave_site_precip_6 <- my_data_200_1$my.data$precip[ave_site_Ind, 6]
        ave_site_precip_7 <- my_data_200_1$my.data$precip[ave_site_Ind, 7]
        ave_site_precip_8 <- my_data_200_1$my.data$precip[ave_site_Ind, 8]
        ave_site_precip_9 <- my_data_200_1$my.data$precip[ave_site_Ind, 9]
        ave_site_precip_10 <- my_data_200_1$my.data$precip[ave_site_Ind, 10]
        
        ave_site_size <- my_data_200_1$my.data$gridarea[ave_site_Ind]
        
        site_yr_effects_OI_1 <-  sims_matrix_200[,sprintf("psi.site[%d,1]", ave_site_Ind)]
        site_yr_effects_OI_10 <- sims_matrix_200[,sprintf("psi.site[%d,10]", ave_site_Ind)]
        
      }
      
    } else if(site_type=="mean_temp_low_precip"){
      
      ave_site_tempRange <- quantile(my_data_200_1$my.data$temp[sp_sites, 1])[c(2,4)]
      ave_site_tempInd <- which(between(my_data_200_1$my.data$temp[sp_sites, 1],
                                        ave_site_tempRange[1], ave_site_tempRange[2]), 
                                arr.ind=TRUE)
      
      ave_site_precipRange <- quantile(my_data_200_1$my.data$precip[sp_sites, 1])[c(1,2)]
      ave_site_precipInd <- which(between(my_data_200_1$my.data$precip[sp_sites, 1],
                                          ave_site_precipRange[1], ave_site_precipRange[2]), 
                                  arr.ind=TRUE)
      
      if(length(intersect(ave_site_tempInd, ave_site_precipInd))==0){
        ave_site_temp_1 <- my_data_200_1$my.data$temp[sp_sites, 1]
        ave_site_temp_2 <- my_data_200_1$my.data$temp[sp_sites, 2]
        ave_site_temp_3 <- my_data_200_1$my.data$temp[sp_sites, 3]
        ave_site_temp_4 <- my_data_200_1$my.data$temp[sp_sites, 4]
        ave_site_temp_5 <- my_data_200_1$my.data$temp[sp_sites, 5]
        ave_site_temp_6 <- my_data_200_1$my.data$temp[sp_sites, 6]
        ave_site_temp_7 <- my_data_200_1$my.data$temp[sp_sites, 7]
        ave_site_temp_8 <- my_data_200_1$my.data$temp[sp_sites, 8]
        ave_site_temp_9 <- my_data_200_1$my.data$temp[sp_sites, 9]
        ave_site_temp_10 <- my_data_200_1$my.data$temp[sp_sites, 10]
        
        ave_site_precip_1 <- my_data_200_1$my.data$precip[ave_site_precipInd, 1]
        ave_site_precip_2 <- my_data_200_1$my.data$precip[ave_site_precipInd, 2]
        ave_site_precip_3 <- my_data_200_1$my.data$precip[ave_site_precipInd, 3]
        ave_site_precip_4 <- my_data_200_1$my.data$precip[ave_site_precipInd, 4]
        ave_site_precip_5 <- my_data_200_1$my.data$precip[ave_site_precipInd, 5]
        ave_site_precip_6 <- my_data_200_1$my.data$precip[ave_site_precipInd, 6]
        ave_site_precip_7 <- my_data_200_1$my.data$precip[ave_site_precipInd, 7]
        ave_site_precip_8 <- my_data_200_1$my.data$precip[ave_site_precipInd, 8]
        ave_site_precip_9 <- my_data_200_1$my.data$precip[ave_site_precipInd, 9]
        ave_site_precip_10 <- my_data_200_1$my.data$precip[ave_site_precipInd, 10]
        
        ave_site_size <- my_data_200_1$my.data$gridarea[ave_site_precipInd]
        
        site_yr_effects_OI_1 <-  sims_matrix_200[,sprintf("psi.site[%d,1]", ave_site_precipInd)]
        site_yr_effects_OI_10 <- sims_matrix_200[,sprintf("psi.site[%d,10]", ave_site_precipInd)]
        
        ave_site_Ind <- ave_site_precipInd
        
      } else{
        ave_site_Ind <- intersect(ave_site_tempInd, ave_site_precipInd)
        
        ave_site_temp_1 <- my_data_200_1$my.data$temp[ave_site_Ind, 1]
        ave_site_temp_2 <- my_data_200_1$my.data$temp[ave_site_Ind, 2]
        ave_site_temp_3 <- my_data_200_1$my.data$temp[ave_site_Ind, 3]
        ave_site_temp_4 <- my_data_200_1$my.data$temp[ave_site_Ind, 4]
        ave_site_temp_5 <- my_data_200_1$my.data$temp[ave_site_Ind, 5]
        ave_site_temp_6 <- my_data_200_1$my.data$temp[ave_site_Ind, 6]
        ave_site_temp_7 <- my_data_200_1$my.data$temp[ave_site_Ind, 7]
        ave_site_temp_8 <- my_data_200_1$my.data$temp[ave_site_Ind, 8]
        ave_site_temp_9 <- my_data_200_1$my.data$temp[ave_site_Ind, 9]
        ave_site_temp_10 <- my_data_200_1$my.data$temp[ave_site_Ind, 10]
        
        ave_site_precip_1 <- my_data_200_1$my.data$precip[ave_site_Ind, 1]
        ave_site_precip_2 <- my_data_200_1$my.data$precip[ave_site_Ind, 2]
        ave_site_precip_3 <- my_data_200_1$my.data$precip[ave_site_Ind, 3]
        ave_site_precip_4 <- my_data_200_1$my.data$precip[ave_site_Ind, 4]
        ave_site_precip_5 <- my_data_200_1$my.data$precip[ave_site_Ind, 5]
        ave_site_precip_6 <- my_data_200_1$my.data$precip[ave_site_Ind, 6]
        ave_site_precip_7 <- my_data_200_1$my.data$precip[ave_site_Ind, 7]
        ave_site_precip_8 <- my_data_200_1$my.data$precip[ave_site_Ind, 8]
        ave_site_precip_9 <- my_data_200_1$my.data$precip[ave_site_Ind, 9]
        ave_site_precip_10 <- my_data_200_1$my.data$precip[ave_site_Ind, 10]
        
        ave_site_size <- my_data_200_1$my.data$gridarea[ave_site_Ind]
        
        site_yr_effects_OI_1 <-  sims_matrix_200[,sprintf("psi.site[%d,1]", ave_site_Ind)]
        site_yr_effects_OI_10 <- sims_matrix_200[,sprintf("psi.site[%d,10]", ave_site_Ind)]
        
      }
      
    } else if(site_type=="high_temp_mean_precip"){
      
      ave_site_tempRange <- quantile(my_data_200_1$my.data$temp[sp_sites, 1])[c(4,5)]
      ave_site_tempInd <- which(between(my_data_200_1$my.data$temp[sp_sites, 1],
                                        ave_site_tempRange[1], ave_site_tempRange[2]), 
                                arr.ind=TRUE)
      
      ave_site_precipRange <- quantile(my_data_200_1$my.data$precip[sp_sites, 1])[c(2,4)]
      ave_site_precipInd <- which(between(my_data_200_1$my.data$precip[sp_sites, 1],
                                          ave_site_precipRange[1], ave_site_precipRange[2]), 
                                  arr.ind=TRUE)
      
      if(length(intersect(ave_site_tempInd, ave_site_precipInd))==0){
        ave_site_temp_1 <- my_data_200_1$my.data$temp[ave_site_tempInd, 1]
        ave_site_temp_2 <- my_data_200_1$my.data$temp[ave_site_tempInd, 2]
        ave_site_temp_3 <- my_data_200_1$my.data$temp[ave_site_tempInd, 3]
        ave_site_temp_4 <- my_data_200_1$my.data$temp[ave_site_tempInd, 4]
        ave_site_temp_5 <- my_data_200_1$my.data$temp[ave_site_tempInd, 5]
        ave_site_temp_6 <- my_data_200_1$my.data$temp[ave_site_tempInd, 6]
        ave_site_temp_7 <- my_data_200_1$my.data$temp[ave_site_tempInd, 7]
        ave_site_temp_8 <- my_data_200_1$my.data$temp[ave_site_tempInd, 8]
        ave_site_temp_9 <- my_data_200_1$my.data$temp[ave_site_tempInd, 9]
        ave_site_temp_10 <- my_data_200_1$my.data$temp[ave_site_tempInd, 10]
        
        ave_site_precip_1 <- my_data_200_1$my.data$precip[sp_sites, 1]
        ave_site_precip_2 <- my_data_200_1$my.data$precip[sp_sites, 2]
        ave_site_precip_3 <- my_data_200_1$my.data$precip[sp_sites, 3]
        ave_site_precip_4 <- my_data_200_1$my.data$precip[sp_sites, 4]
        ave_site_precip_5 <- my_data_200_1$my.data$precip[sp_sites, 5]
        ave_site_precip_6 <- my_data_200_1$my.data$precip[sp_sites, 6]
        ave_site_precip_7 <- my_data_200_1$my.data$precip[sp_sites, 7]
        ave_site_precip_8 <- my_data_200_1$my.data$precip[sp_sites, 8]
        ave_site_precip_9 <- my_data_200_1$my.data$precip[sp_sites, 9]
        ave_site_precip_10 <- my_data_200_1$my.data$precip[sp_sites, 10]
        
        ave_site_size <- my_data_200_1$my.data$gridarea[ave_site_tempInd]
        
        site_yr_effects_OI_1 <-  sims_matrix_200[,sprintf("psi.site[%d,1]", ave_site_tempInd)]
        site_yr_effects_OI_10 <- sims_matrix_200[,sprintf("psi.site[%d,10]", ave_site_tempInd)]
        
        
        ave_site_Ind <- ave_site_tempInd
        
      } else{
        ave_site_Ind <- intersect(ave_site_tempInd, ave_site_precipInd)
        
        ave_site_temp_1 <- my_data_200_1$my.data$temp[ave_site_Ind, 1]
        ave_site_temp_2 <- my_data_200_1$my.data$temp[ave_site_Ind, 2]
        ave_site_temp_3 <- my_data_200_1$my.data$temp[ave_site_Ind, 3]
        ave_site_temp_4 <- my_data_200_1$my.data$temp[ave_site_Ind, 4]
        ave_site_temp_5 <- my_data_200_1$my.data$temp[ave_site_Ind, 5]
        ave_site_temp_6 <- my_data_200_1$my.data$temp[ave_site_Ind, 6]
        ave_site_temp_7 <- my_data_200_1$my.data$temp[ave_site_Ind, 7]
        ave_site_temp_8 <- my_data_200_1$my.data$temp[ave_site_Ind, 8]
        ave_site_temp_9 <- my_data_200_1$my.data$temp[ave_site_Ind, 9]
        ave_site_temp_10 <- my_data_200_1$my.data$temp[ave_site_Ind, 10]
        
        ave_site_precip_1 <- my_data_200_1$my.data$precip[ave_site_Ind, 1]
        ave_site_precip_2 <- my_data_200_1$my.data$precip[ave_site_Ind, 2]
        ave_site_precip_3 <- my_data_200_1$my.data$precip[ave_site_Ind, 3]
        ave_site_precip_4 <- my_data_200_1$my.data$precip[ave_site_Ind, 4]
        ave_site_precip_5 <- my_data_200_1$my.data$precip[ave_site_Ind, 5]
        ave_site_precip_6 <- my_data_200_1$my.data$precip[ave_site_Ind, 6]
        ave_site_precip_7 <- my_data_200_1$my.data$precip[ave_site_Ind, 7]
        ave_site_precip_8 <- my_data_200_1$my.data$precip[ave_site_Ind, 8]
        ave_site_precip_9 <- my_data_200_1$my.data$precip[ave_site_Ind, 9]
        ave_site_precip_10 <- my_data_200_1$my.data$precip[ave_site_Ind, 10]
        
        ave_site_size <- my_data_200_1$my.data$gridarea[ave_site_Ind]
        
        site_yr_effects_OI_1 <-  sims_matrix_200[,sprintf("psi.site[%d,1]", ave_site_Ind)]
        site_yr_effects_OI_10 <- sims_matrix_200[,sprintf("psi.site[%d,10]", ave_site_Ind)]
        
      }
      
    } else if(site_type=="mean_temp_high_precip"){
      
      ave_site_tempRange <- quantile(my_data_200_1$my.data$temp[sp_sites, 1])[c(2,4)]
      ave_site_tempInd <- which(between(my_data_200_1$my.data$temp[sp_sites, 1],
                                        ave_site_tempRange[1], ave_site_tempRange[2]), 
                                arr.ind=TRUE)
      
      ave_site_precipRange <- quantile(my_data_200_1$my.data$precip[sp_sites, 1])[c(4,5)]
      ave_site_precipInd <- which(between(my_data_200_1$my.data$precip[sp_sites, 1],
                                          ave_site_precipRange[1], ave_site_precipRange[2]), 
                                  arr.ind=TRUE)
      
      if(length(intersect(ave_site_tempInd, ave_site_precipInd))==0){
        ave_site_temp_1 <- my_data_200_1$my.data$temp[sp_sites, 1]
        ave_site_temp_2 <- my_data_200_1$my.data$temp[sp_sites, 2]
        ave_site_temp_3 <- my_data_200_1$my.data$temp[sp_sites, 3]
        ave_site_temp_4 <- my_data_200_1$my.data$temp[sp_sites, 4]
        ave_site_temp_5 <- my_data_200_1$my.data$temp[sp_sites, 5]
        ave_site_temp_6 <- my_data_200_1$my.data$temp[sp_sites, 6]
        ave_site_temp_7 <- my_data_200_1$my.data$temp[sp_sites, 7]
        ave_site_temp_8 <- my_data_200_1$my.data$temp[sp_sites, 8]
        ave_site_temp_9 <- my_data_200_1$my.data$temp[sp_sites, 9]
        ave_site_temp_10 <- my_data_200_1$my.data$temp[sp_sites, 10]
        
        ave_site_precip_1 <- my_data_200_1$my.data$precip[ave_site_precipInd, 1]
        ave_site_precip_2 <- my_data_200_1$my.data$precip[ave_site_precipInd, 2]
        ave_site_precip_3 <- my_data_200_1$my.data$precip[ave_site_precipInd, 3]
        ave_site_precip_4 <- my_data_200_1$my.data$precip[ave_site_precipInd, 4]
        ave_site_precip_5 <- my_data_200_1$my.data$precip[ave_site_precipInd, 5]
        ave_site_precip_6 <- my_data_200_1$my.data$precip[ave_site_precipInd, 6]
        ave_site_precip_7 <- my_data_200_1$my.data$precip[ave_site_precipInd, 7]
        ave_site_precip_8 <- my_data_200_1$my.data$precip[ave_site_precipInd, 8]
        ave_site_precip_9 <- my_data_200_1$my.data$precip[ave_site_precipInd, 9]
        ave_site_precip_10 <- my_data_200_1$my.data$precip[ave_site_precipInd, 10]
        
        ave_site_size <- my_data_200_1$my.data$gridarea[ave_site_precipInd]
        
        site_yr_effects_OI_1 <-  sims_matrix_200[,sprintf("psi.site[%d,1]", ave_site_precipInd)]
        site_yr_effects_OI_10 <- sims_matrix_200[,sprintf("psi.site[%d,10]", ave_site_precipInd)]
        
        
        ave_site_Ind <- ave_site_precipInd
        
      } else{
        ave_site_Ind <- intersect(ave_site_tempInd, ave_site_precipInd)
        
        ave_site_temp_1 <- my_data_200_1$my.data$temp[ave_site_Ind, 1]
        ave_site_temp_2 <- my_data_200_1$my.data$temp[ave_site_Ind, 2]
        ave_site_temp_3 <- my_data_200_1$my.data$temp[ave_site_Ind, 3]
        ave_site_temp_4 <- my_data_200_1$my.data$temp[ave_site_Ind, 4]
        ave_site_temp_5 <- my_data_200_1$my.data$temp[ave_site_Ind, 5]
        ave_site_temp_6 <- my_data_200_1$my.data$temp[ave_site_Ind, 6]
        ave_site_temp_7 <- my_data_200_1$my.data$temp[ave_site_Ind, 7]
        ave_site_temp_8 <- my_data_200_1$my.data$temp[ave_site_Ind, 8]
        ave_site_temp_9 <- my_data_200_1$my.data$temp[ave_site_Ind, 9]
        ave_site_temp_10 <- my_data_200_1$my.data$temp[ave_site_Ind, 10]
        
        ave_site_precip_1 <- my_data_200_1$my.data$precip[ave_site_Ind, 1]
        ave_site_precip_2 <- my_data_200_1$my.data$precip[ave_site_Ind, 2]
        ave_site_precip_3 <- my_data_200_1$my.data$precip[ave_site_Ind, 3]
        ave_site_precip_4 <- my_data_200_1$my.data$precip[ave_site_Ind, 4]
        ave_site_precip_5 <- my_data_200_1$my.data$precip[ave_site_Ind, 5]
        ave_site_precip_6 <- my_data_200_1$my.data$precip[ave_site_Ind, 6]
        ave_site_precip_7 <- my_data_200_1$my.data$precip[ave_site_Ind, 7]
        ave_site_precip_8 <- my_data_200_1$my.data$precip[ave_site_Ind, 8]
        ave_site_precip_9 <- my_data_200_1$my.data$precip[ave_site_Ind, 9]
        ave_site_precip_10 <- my_data_200_1$my.data$precip[ave_site_Ind, 10]
        
        ave_site_size <- my_data_200_1$my.data$gridarea[ave_site_Ind]
        
        site_yr_effects_OI_1 <-  sims_matrix_200[,sprintf("psi.site[%d,1]", ave_site_Ind)]
        site_yr_effects_OI_10 <- sims_matrix_200[,sprintf("psi.site[%d,10]", ave_site_Ind)]
        
      }
    }
    
    sp_occ_OI_1 <- plogis(sims_matrix_200[,"mu.psi.0"]+
                            sims_matrix_200[,"psi.area"]*mean(ave_site_size)+ 
                            sims_matrix_200[,sprintf("psi.beta.temp[%d]", i)]*mean(ave_site_temp_1)+
                            sims_matrix_200[,sprintf("psi.beta.precip[%d]", i)]*mean(ave_site_precip_1)+
                            site_yr_effects_OI_1+
                            sims_matrix_200[,sprintf("psi.sp[%d]", i)])
    
    sp_occ_OI_10 <- plogis(sims_matrix_200[,"mu.psi.0"]+
                             sims_matrix_200[,"psi.area"]*mean(ave_site_size)+ 
                             sims_matrix_200[,sprintf("psi.beta.temp[%d]", i)]*mean(ave_site_temp_10)+
                             sims_matrix_200[,sprintf("psi.beta.precip[%d]", i)]*mean(ave_site_precip_10)+
                             site_yr_effects_OI_10)
    
    sp_occ_OI_diff <- sp_occ_OI_10 - sp_occ_OI_1
    
    
    sp_occ_diff_200[[i]] <- c(mean(sp_occ_OI_diff), 
                              quantile(sp_occ_OI_diff, probs=0.025, na.rm=TRUE), 
                              quantile(sp_occ_OI_diff, probs=0.975, na.rm=TRUE))
    
    ave_site_Ind_200 <- ave_site_Ind
    
    # 100 x 100 Kilometer Analysis
    sp_sites <- my_data_100_1$my.info$range.list[i,] %>% na.omit()
    
    # AVERAGE SITE IN 1970-1974
    if(site_type=="mean"){
      
      ave_site_tempRange <- quantile(my_data_100_1$my.data$temp[sp_sites, 1])[c(2,4)]
      ave_site_tempInd <- which(between(my_data_100_1$my.data$temp[sp_sites, 1],
                                        ave_site_tempRange[1], ave_site_tempRange[2]), 
                                arr.ind=TRUE)
      
      ave_site_precipRange <- quantile(my_data_100_1$my.data$precip[sp_sites, 1])[c(2,4)]
      ave_site_precipInd <- which(between(my_data_100_1$my.data$precip[sp_sites, 1],
                                          ave_site_precipRange[1], ave_site_precipRange[2]), 
                                  arr.ind=TRUE)
      
      ave_site_Ind <- intersect(ave_site_tempInd, ave_site_precipInd)
      
      ave_site_size <- my_data_100_1$my.data$gridarea[ave_site_Ind]
      
      site_yr_effects_OI_1 <-  sims_matrix_100[,sprintf("psi.site[%d,1]", ave_site_Ind)]
      site_yr_effects_OI_10 <- sims_matrix_100[,sprintf("psi.site[%d,10]", ave_site_Ind)]
      
      ave_site_temp_1 <- my_data_100_1$my.data$temp[ave_site_Ind, 1]
      ave_site_temp_2 <- my_data_100_1$my.data$temp[ave_site_Ind, 2]
      ave_site_temp_3 <- my_data_100_1$my.data$temp[ave_site_Ind, 3]
      ave_site_temp_4 <- my_data_100_1$my.data$temp[ave_site_Ind, 4]
      ave_site_temp_5 <- my_data_100_1$my.data$temp[ave_site_Ind, 5]
      ave_site_temp_6 <- my_data_100_1$my.data$temp[ave_site_Ind, 6]
      ave_site_temp_7 <- my_data_100_1$my.data$temp[ave_site_Ind, 7]
      ave_site_temp_8 <- my_data_100_1$my.data$temp[ave_site_Ind, 8]
      ave_site_temp_9 <- my_data_100_1$my.data$temp[ave_site_Ind, 9]
      ave_site_temp_10 <- my_data_100_1$my.data$temp[ave_site_Ind, 10]
      
      ave_site_precip_1 <- my_data_100_1$my.data$precip[ave_site_Ind, 1]
      ave_site_precip_2 <- my_data_100_1$my.data$precip[ave_site_Ind, 2]
      ave_site_precip_3 <- my_data_100_1$my.data$precip[ave_site_Ind, 3]
      ave_site_precip_4 <- my_data_100_1$my.data$precip[ave_site_Ind, 4]
      ave_site_precip_5 <- my_data_100_1$my.data$precip[ave_site_Ind, 5]
      ave_site_precip_6 <- my_data_100_1$my.data$precip[ave_site_Ind, 6]
      ave_site_precip_7 <- my_data_100_1$my.data$precip[ave_site_Ind, 7]
      ave_site_precip_8 <- my_data_100_1$my.data$precip[ave_site_Ind, 8]
      ave_site_precip_9 <- my_data_100_1$my.data$precip[ave_site_Ind, 9]
      ave_site_precip_10 <- my_data_100_1$my.data$precip[ave_site_Ind, 10]
      
    } else if(site_type=="low_temp_mean_precip"){
      
      ave_site_tempRange <- quantile(my_data_100_1$my.data$temp[sp_sites, 1])[c(1,2)]
      ave_site_tempInd <- which(between(my_data_100_1$my.data$temp[sp_sites, 1],
                                        ave_site_tempRange[1], ave_site_tempRange[2]), 
                                arr.ind=TRUE)
      
      ave_site_precipRange <- quantile(my_data_100_1$my.data$precip[sp_sites, 1])[c(2,4)]
      ave_site_precipInd <- which(between(my_data_100_1$my.data$precip[sp_sites, 1],
                                          ave_site_precipRange[1], ave_site_precipRange[2]), 
                                  arr.ind=TRUE)
      
      if(length(intersect(ave_site_tempInd, ave_site_precipInd))==0){
        ave_site_temp_1 <- my_data_100_1$my.data$temp[ave_site_tempInd, 1]
        ave_site_temp_2 <- my_data_100_1$my.data$temp[ave_site_tempInd, 2]
        ave_site_temp_3 <- my_data_100_1$my.data$temp[ave_site_tempInd, 3]
        ave_site_temp_4 <- my_data_100_1$my.data$temp[ave_site_tempInd, 4]
        ave_site_temp_5 <- my_data_100_1$my.data$temp[ave_site_tempInd, 5]
        ave_site_temp_6 <- my_data_100_1$my.data$temp[ave_site_tempInd, 6]
        ave_site_temp_7 <- my_data_100_1$my.data$temp[ave_site_tempInd, 7]
        ave_site_temp_8 <- my_data_100_1$my.data$temp[ave_site_tempInd, 8]
        ave_site_temp_9 <- my_data_100_1$my.data$temp[ave_site_tempInd, 9]
        ave_site_temp_10 <- my_data_100_1$my.data$temp[ave_site_tempInd, 10]
        
        ave_site_precip_1 <- my_data_100_1$my.data$precip[sp_sites, 1]
        ave_site_precip_2 <- my_data_100_1$my.data$precip[sp_sites, 2]
        ave_site_precip_3 <- my_data_100_1$my.data$precip[sp_sites, 3]
        ave_site_precip_4 <- my_data_100_1$my.data$precip[sp_sites, 4]
        ave_site_precip_5 <- my_data_100_1$my.data$precip[sp_sites, 5]
        ave_site_precip_6 <- my_data_100_1$my.data$precip[sp_sites, 6]
        ave_site_precip_7 <- my_data_100_1$my.data$precip[sp_sites, 7]
        ave_site_precip_8 <- my_data_100_1$my.data$precip[sp_sites, 8]
        ave_site_precip_9 <- my_data_100_1$my.data$precip[sp_sites, 9]
        ave_site_precip_10 <- my_data_100_1$my.data$precip[sp_sites, 10]
        
        ave_site_size <- my_data_100_1$my.data$gridarea[ave_site_tempInd]
        
        site_yr_effects_OI_1 <-  sims_matrix_100[,sprintf("psi.site[%d,1]", ave_site_tempInd)]
        site_yr_effects_OI_10 <- sims_matrix_100[,sprintf("psi.site[%d,10]", ave_site_tempInd)]
        
        ave_site_Ind <- ave_site_tempInd
        
      } else{
        ave_site_Ind <- intersect(ave_site_tempInd, ave_site_precipInd)
        
        ave_site_temp_1 <- my_data_100_1$my.data$temp[ave_site_Ind, 1]
        ave_site_temp_2 <- my_data_100_1$my.data$temp[ave_site_Ind, 2]
        ave_site_temp_3 <- my_data_100_1$my.data$temp[ave_site_Ind, 3]
        ave_site_temp_4 <- my_data_100_1$my.data$temp[ave_site_Ind, 4]
        ave_site_temp_5 <- my_data_100_1$my.data$temp[ave_site_Ind, 5]
        ave_site_temp_6 <- my_data_100_1$my.data$temp[ave_site_Ind, 6]
        ave_site_temp_7 <- my_data_100_1$my.data$temp[ave_site_Ind, 7]
        ave_site_temp_8 <- my_data_100_1$my.data$temp[ave_site_Ind, 8]
        ave_site_temp_9 <- my_data_100_1$my.data$temp[ave_site_Ind, 9]
        ave_site_temp_10 <- my_data_100_1$my.data$temp[ave_site_Ind, 10]
        
        ave_site_precip_1 <- my_data_100_1$my.data$precip[ave_site_Ind, 1]
        ave_site_precip_2 <- my_data_100_1$my.data$precip[ave_site_Ind, 2]
        ave_site_precip_3 <- my_data_100_1$my.data$precip[ave_site_Ind, 3]
        ave_site_precip_4 <- my_data_100_1$my.data$precip[ave_site_Ind, 4]
        ave_site_precip_5 <- my_data_100_1$my.data$precip[ave_site_Ind, 5]
        ave_site_precip_6 <- my_data_100_1$my.data$precip[ave_site_Ind, 6]
        ave_site_precip_7 <- my_data_100_1$my.data$precip[ave_site_Ind, 7]
        ave_site_precip_8 <- my_data_100_1$my.data$precip[ave_site_Ind, 8]
        ave_site_precip_9 <- my_data_100_1$my.data$precip[ave_site_Ind, 9]
        ave_site_precip_10 <- my_data_100_1$my.data$precip[ave_site_Ind, 10]
        
        ave_site_size <- my_data_100_1$my.data$gridarea[ave_site_Ind]
        
        site_yr_effects_OI_1 <-  sims_matrix_100[,sprintf("psi.site[%d,1]", ave_site_Ind)]
        site_yr_effects_OI_10 <- sims_matrix_100[,sprintf("psi.site[%d,10]", ave_site_Ind)]
        
      }
      
    } else if(site_type=="mean_temp_low_precip"){
      
      ave_site_tempRange <- quantile(my_data_100_1$my.data$temp[sp_sites, 1])[c(2,4)]
      ave_site_tempInd <- which(between(my_data_100_1$my.data$temp[sp_sites, 1],
                                        ave_site_tempRange[1], ave_site_tempRange[2]), 
                                arr.ind=TRUE)
      
      ave_site_precipRange <- quantile(my_data_100_1$my.data$precip[sp_sites, 1])[c(1,2)]
      ave_site_precipInd <- which(between(my_data_100_1$my.data$precip[sp_sites, 1],
                                          ave_site_precipRange[1], ave_site_precipRange[2]), 
                                  arr.ind=TRUE)
      
      if(length(intersect(ave_site_tempInd, ave_site_precipInd))==0){
        ave_site_temp_1 <- my_data_100_1$my.data$temp[sp_sites, 1]
        ave_site_temp_2 <- my_data_100_1$my.data$temp[sp_sites, 2]
        ave_site_temp_3 <- my_data_100_1$my.data$temp[sp_sites, 3]
        ave_site_temp_4 <- my_data_100_1$my.data$temp[sp_sites, 4]
        ave_site_temp_5 <- my_data_100_1$my.data$temp[sp_sites, 5]
        ave_site_temp_6 <- my_data_100_1$my.data$temp[sp_sites, 6]
        ave_site_temp_7 <- my_data_100_1$my.data$temp[sp_sites, 7]
        ave_site_temp_8 <- my_data_100_1$my.data$temp[sp_sites, 8]
        ave_site_temp_9 <- my_data_100_1$my.data$temp[sp_sites, 9]
        ave_site_temp_10 <- my_data_100_1$my.data$temp[sp_sites, 10]
        
        ave_site_precip_1 <- my_data_100_1$my.data$precip[ave_site_precipInd, 1]
        ave_site_precip_2 <- my_data_100_1$my.data$precip[ave_site_precipInd, 2]
        ave_site_precip_3 <- my_data_100_1$my.data$precip[ave_site_precipInd, 3]
        ave_site_precip_4 <- my_data_100_1$my.data$precip[ave_site_precipInd, 4]
        ave_site_precip_5 <- my_data_100_1$my.data$precip[ave_site_precipInd, 5]
        ave_site_precip_6 <- my_data_100_1$my.data$precip[ave_site_precipInd, 6]
        ave_site_precip_7 <- my_data_100_1$my.data$precip[ave_site_precipInd, 7]
        ave_site_precip_8 <- my_data_100_1$my.data$precip[ave_site_precipInd, 8]
        ave_site_precip_9 <- my_data_100_1$my.data$precip[ave_site_precipInd, 9]
        ave_site_precip_10 <- my_data_100_1$my.data$precip[ave_site_precipInd, 10]
        
        ave_site_size <- my_data_100_1$my.data$gridarea[ave_site_precipInd]
        
        site_yr_effects_OI_1 <-  sims_matrix_100[,sprintf("psi.site[%d,1]", ave_site_precipInd)]
        site_yr_effects_OI_10 <- sims_matrix_100[,sprintf("psi.site[%d,10]", ave_site_precipInd)]
        
        ave_site_Ind <- ave_site_precipInd
        
      } else{
        ave_site_Ind <- intersect(ave_site_tempInd, ave_site_precipInd)
        
        ave_site_temp_1 <- my_data_100_1$my.data$temp[ave_site_Ind, 1]
        ave_site_temp_2 <- my_data_100_1$my.data$temp[ave_site_Ind, 2]
        ave_site_temp_3 <- my_data_100_1$my.data$temp[ave_site_Ind, 3]
        ave_site_temp_4 <- my_data_100_1$my.data$temp[ave_site_Ind, 4]
        ave_site_temp_5 <- my_data_100_1$my.data$temp[ave_site_Ind, 5]
        ave_site_temp_6 <- my_data_100_1$my.data$temp[ave_site_Ind, 6]
        ave_site_temp_7 <- my_data_100_1$my.data$temp[ave_site_Ind, 7]
        ave_site_temp_8 <- my_data_100_1$my.data$temp[ave_site_Ind, 8]
        ave_site_temp_9 <- my_data_100_1$my.data$temp[ave_site_Ind, 9]
        ave_site_temp_10 <- my_data_100_1$my.data$temp[ave_site_Ind, 10]
        
        ave_site_precip_1 <- my_data_100_1$my.data$precip[ave_site_Ind, 1]
        ave_site_precip_2 <- my_data_100_1$my.data$precip[ave_site_Ind, 2]
        ave_site_precip_3 <- my_data_100_1$my.data$precip[ave_site_Ind, 3]
        ave_site_precip_4 <- my_data_100_1$my.data$precip[ave_site_Ind, 4]
        ave_site_precip_5 <- my_data_100_1$my.data$precip[ave_site_Ind, 5]
        ave_site_precip_6 <- my_data_100_1$my.data$precip[ave_site_Ind, 6]
        ave_site_precip_7 <- my_data_100_1$my.data$precip[ave_site_Ind, 7]
        ave_site_precip_8 <- my_data_100_1$my.data$precip[ave_site_Ind, 8]
        ave_site_precip_9 <- my_data_100_1$my.data$precip[ave_site_Ind, 9]
        ave_site_precip_10 <- my_data_100_1$my.data$precip[ave_site_Ind, 10]
        
        ave_site_size <- my_data_100_1$my.data$gridarea[ave_site_Ind]
        
        site_yr_effects_OI_1 <-  sims_matrix_100[,sprintf("psi.site[%d,1]", ave_site_Ind)]
        site_yr_effects_OI_10 <- sims_matrix_100[,sprintf("psi.site[%d,10]", ave_site_Ind)]
        
      }
      
    } else if(site_type=="high_temp_mean_precip"){
      
      ave_site_tempRange <- quantile(my_data_100_1$my.data$temp[sp_sites, 1])[c(4,5)]
      ave_site_tempInd <- which(between(my_data_100_1$my.data$temp[sp_sites, 1],
                                        ave_site_tempRange[1], ave_site_tempRange[2]), 
                                arr.ind=TRUE)
      
      ave_site_precipRange <- quantile(my_data_100_1$my.data$precip[sp_sites, 1])[c(2,4)]
      ave_site_precipInd <- which(between(my_data_100_1$my.data$precip[sp_sites, 1],
                                          ave_site_precipRange[1], ave_site_precipRange[2]), 
                                  arr.ind=TRUE)
      
      if(length(intersect(ave_site_tempInd, ave_site_precipInd))==0){
        ave_site_temp_1 <- my_data_100_1$my.data$temp[ave_site_tempInd, 1]
        ave_site_temp_2 <- my_data_100_1$my.data$temp[ave_site_tempInd, 2]
        ave_site_temp_3 <- my_data_100_1$my.data$temp[ave_site_tempInd, 3]
        ave_site_temp_4 <- my_data_100_1$my.data$temp[ave_site_tempInd, 4]
        ave_site_temp_5 <- my_data_100_1$my.data$temp[ave_site_tempInd, 5]
        ave_site_temp_6 <- my_data_100_1$my.data$temp[ave_site_tempInd, 6]
        ave_site_temp_7 <- my_data_100_1$my.data$temp[ave_site_tempInd, 7]
        ave_site_temp_8 <- my_data_100_1$my.data$temp[ave_site_tempInd, 8]
        ave_site_temp_9 <- my_data_100_1$my.data$temp[ave_site_tempInd, 9]
        ave_site_temp_10 <- my_data_100_1$my.data$temp[ave_site_tempInd, 10]
        
        ave_site_precip_1 <- my_data_100_1$my.data$precip[sp_sites, 1]
        ave_site_precip_2 <- my_data_100_1$my.data$precip[sp_sites, 2]
        ave_site_precip_3 <- my_data_100_1$my.data$precip[sp_sites, 3]
        ave_site_precip_4 <- my_data_100_1$my.data$precip[sp_sites, 4]
        ave_site_precip_5 <- my_data_100_1$my.data$precip[sp_sites, 5]
        ave_site_precip_6 <- my_data_100_1$my.data$precip[sp_sites, 6]
        ave_site_precip_7 <- my_data_100_1$my.data$precip[sp_sites, 7]
        ave_site_precip_8 <- my_data_100_1$my.data$precip[sp_sites, 8]
        ave_site_precip_9 <- my_data_100_1$my.data$precip[sp_sites, 9]
        ave_site_precip_10 <- my_data_100_1$my.data$precip[sp_sites, 10]
        
        ave_site_size <- my_data_100_1$my.data$gridarea[ave_site_tempInd]
        
        site_yr_effects_OI_1 <-  sims_matrix_100[,sprintf("psi.site[%d,1]", ave_site_tempInd)]
        site_yr_effects_OI_10 <- sims_matrix_100[,sprintf("psi.site[%d,10]", ave_site_tempInd)]
        
        
        ave_site_Ind <- ave_site_tempInd
        
      } else{
        ave_site_Ind <- intersect(ave_site_tempInd, ave_site_precipInd)
        
        ave_site_temp_1 <- my_data_100_1$my.data$temp[ave_site_Ind, 1]
        ave_site_temp_2 <- my_data_100_1$my.data$temp[ave_site_Ind, 2]
        ave_site_temp_3 <- my_data_100_1$my.data$temp[ave_site_Ind, 3]
        ave_site_temp_4 <- my_data_100_1$my.data$temp[ave_site_Ind, 4]
        ave_site_temp_5 <- my_data_100_1$my.data$temp[ave_site_Ind, 5]
        ave_site_temp_6 <- my_data_100_1$my.data$temp[ave_site_Ind, 6]
        ave_site_temp_7 <- my_data_100_1$my.data$temp[ave_site_Ind, 7]
        ave_site_temp_8 <- my_data_100_1$my.data$temp[ave_site_Ind, 8]
        ave_site_temp_9 <- my_data_100_1$my.data$temp[ave_site_Ind, 9]
        ave_site_temp_10 <- my_data_100_1$my.data$temp[ave_site_Ind, 10]
        
        ave_site_precip_1 <- my_data_100_1$my.data$precip[ave_site_Ind, 1]
        ave_site_precip_2 <- my_data_100_1$my.data$precip[ave_site_Ind, 2]
        ave_site_precip_3 <- my_data_100_1$my.data$precip[ave_site_Ind, 3]
        ave_site_precip_4 <- my_data_100_1$my.data$precip[ave_site_Ind, 4]
        ave_site_precip_5 <- my_data_100_1$my.data$precip[ave_site_Ind, 5]
        ave_site_precip_6 <- my_data_100_1$my.data$precip[ave_site_Ind, 6]
        ave_site_precip_7 <- my_data_100_1$my.data$precip[ave_site_Ind, 7]
        ave_site_precip_8 <- my_data_100_1$my.data$precip[ave_site_Ind, 8]
        ave_site_precip_9 <- my_data_100_1$my.data$precip[ave_site_Ind, 9]
        ave_site_precip_10 <- my_data_100_1$my.data$precip[ave_site_Ind, 10]
        
        ave_site_size <- my_data_100_1$my.data$gridarea[ave_site_Ind]
        
        site_yr_effects_OI_1 <-  sims_matrix_100[,sprintf("psi.site[%d,1]", ave_site_Ind)]
        site_yr_effects_OI_10 <- sims_matrix_100[,sprintf("psi.site[%d,10]", ave_site_Ind)]
        
      }
      
    } else if(site_type=="mean_temp_high_precip"){
      
      ave_site_tempRange <- quantile(my_data_100_1$my.data$temp[sp_sites, 1])[c(2,4)]
      ave_site_tempInd <- which(between(my_data_100_1$my.data$temp[sp_sites, 1],
                                        ave_site_tempRange[1], ave_site_tempRange[2]), 
                                arr.ind=TRUE)
      
      ave_site_precipRange <- quantile(my_data_100_1$my.data$precip[sp_sites, 1])[c(4,5)]
      ave_site_precipInd <- which(between(my_data_100_1$my.data$precip[sp_sites, 1],
                                          ave_site_precipRange[1], ave_site_precipRange[2]), 
                                  arr.ind=TRUE)
      
      if(length(intersect(ave_site_tempInd, ave_site_precipInd))==0){
        ave_site_temp_1 <- my_data_100_1$my.data$temp[sp_sites, 1]
        ave_site_temp_2 <- my_data_100_1$my.data$temp[sp_sites, 2]
        ave_site_temp_3 <- my_data_100_1$my.data$temp[sp_sites, 3]
        ave_site_temp_4 <- my_data_100_1$my.data$temp[sp_sites, 4]
        ave_site_temp_5 <- my_data_100_1$my.data$temp[sp_sites, 5]
        ave_site_temp_6 <- my_data_100_1$my.data$temp[sp_sites, 6]
        ave_site_temp_7 <- my_data_100_1$my.data$temp[sp_sites, 7]
        ave_site_temp_8 <- my_data_100_1$my.data$temp[sp_sites, 8]
        ave_site_temp_9 <- my_data_100_1$my.data$temp[sp_sites, 9]
        ave_site_temp_10 <- my_data_100_1$my.data$temp[sp_sites, 10]
        
        ave_site_precip_1 <- my_data_100_1$my.data$precip[ave_site_precipInd, 1]
        ave_site_precip_2 <- my_data_100_1$my.data$precip[ave_site_precipInd, 2]
        ave_site_precip_3 <- my_data_100_1$my.data$precip[ave_site_precipInd, 3]
        ave_site_precip_4 <- my_data_100_1$my.data$precip[ave_site_precipInd, 4]
        ave_site_precip_5 <- my_data_100_1$my.data$precip[ave_site_precipInd, 5]
        ave_site_precip_6 <- my_data_100_1$my.data$precip[ave_site_precipInd, 6]
        ave_site_precip_7 <- my_data_100_1$my.data$precip[ave_site_precipInd, 7]
        ave_site_precip_8 <- my_data_100_1$my.data$precip[ave_site_precipInd, 8]
        ave_site_precip_9 <- my_data_100_1$my.data$precip[ave_site_precipInd, 9]
        ave_site_precip_10 <- my_data_100_1$my.data$precip[ave_site_precipInd, 10]
        
        ave_site_size <- my_data_100_1$my.data$gridarea[ave_site_precipInd]
        
        site_yr_effects_OI_1 <-  sims_matrix_100[,sprintf("psi.site[%d,1]", ave_site_precipInd)]
        site_yr_effects_OI_10 <- sims_matrix_100[,sprintf("psi.site[%d,10]", ave_site_precipInd)]
        
        
        ave_site_Ind <- ave_site_precipInd
        
      } else{
        ave_site_Ind <- intersect(ave_site_tempInd, ave_site_precipInd)
        
        ave_site_temp_1 <- my_data_100_1$my.data$temp[ave_site_Ind, 1]
        ave_site_temp_2 <- my_data_100_1$my.data$temp[ave_site_Ind, 2]
        ave_site_temp_3 <- my_data_100_1$my.data$temp[ave_site_Ind, 3]
        ave_site_temp_4 <- my_data_100_1$my.data$temp[ave_site_Ind, 4]
        ave_site_temp_5 <- my_data_100_1$my.data$temp[ave_site_Ind, 5]
        ave_site_temp_6 <- my_data_100_1$my.data$temp[ave_site_Ind, 6]
        ave_site_temp_7 <- my_data_100_1$my.data$temp[ave_site_Ind, 7]
        ave_site_temp_8 <- my_data_100_1$my.data$temp[ave_site_Ind, 8]
        ave_site_temp_9 <- my_data_100_1$my.data$temp[ave_site_Ind, 9]
        ave_site_temp_10 <- my_data_100_1$my.data$temp[ave_site_Ind, 10]
        
        ave_site_precip_1 <- my_data_100_1$my.data$precip[ave_site_Ind, 1]
        ave_site_precip_2 <- my_data_100_1$my.data$precip[ave_site_Ind, 2]
        ave_site_precip_3 <- my_data_100_1$my.data$precip[ave_site_Ind, 3]
        ave_site_precip_4 <- my_data_100_1$my.data$precip[ave_site_Ind, 4]
        ave_site_precip_5 <- my_data_100_1$my.data$precip[ave_site_Ind, 5]
        ave_site_precip_6 <- my_data_100_1$my.data$precip[ave_site_Ind, 6]
        ave_site_precip_7 <- my_data_100_1$my.data$precip[ave_site_Ind, 7]
        ave_site_precip_8 <- my_data_100_1$my.data$precip[ave_site_Ind, 8]
        ave_site_precip_9 <- my_data_100_1$my.data$precip[ave_site_Ind, 9]
        ave_site_precip_10 <- my_data_100_1$my.data$precip[ave_site_Ind, 10]
        
        ave_site_size <- my_data_100_1$my.data$gridarea[ave_site_Ind]
        
        site_yr_effects_OI_1 <-  sims_matrix_100[,sprintf("psi.site[%d,1]", ave_site_Ind)]
        site_yr_effects_OI_10 <- sims_matrix_100[,sprintf("psi.site[%d,10]", ave_site_Ind)]
        
      }
    }
    
    sp_occ_OI_1 <- plogis(sims_matrix_100[,"mu.psi.0"]+
                            sims_matrix_100[,"psi.area"]*mean(ave_site_size)+ 
                            sims_matrix_100[,sprintf("psi.beta.temp[%d]", i)]*mean(ave_site_temp_1)+
                            sims_matrix_100[,sprintf("psi.beta.precip[%d]", i)]*mean(ave_site_precip_1)+
                            site_yr_effects_OI_1+
                            sims_matrix_100[,sprintf("psi.sp[%d]", i)])
    
    sp_occ_OI_10 <- plogis(sims_matrix_100[,"mu.psi.0"]+
                             sims_matrix_100[,"psi.area"]*mean(ave_site_size)+ 
                             sims_matrix_100[,sprintf("psi.beta.temp[%d]", i)]*mean(ave_site_temp_10)+
                             sims_matrix_100[,sprintf("psi.beta.precip[%d]", i)]*mean(ave_site_precip_10)+
                             site_yr_effects_OI_10)
    
    sp_occ_OI_diff <- sp_occ_OI_10 - sp_occ_OI_1
    
    
    sp_occ_diff_100[[i]] <- c(mean(sp_occ_OI_diff), 
                              quantile(sp_occ_OI_diff, probs=0.025, na.rm=TRUE), 
                              quantile(sp_occ_OI_diff, probs=0.975, na.rm=TRUE))
    
    ave_site_Ind_100 <- ave_site_Ind
  }
  
  # Merge the occupancy shift estimates into a single dataframe for plotting
  sp_occ_dx_100 <- do.call(rbind, sp_occ_diff_100) %>% as.data.frame()
  colnames(sp_occ_dx_100) <- c("mean", "lower", "upper")
  sp_occ_dx_100$species <- my_traits$species
  
  sp_occ_dx_200 <- do.call(rbind, sp_occ_diff_200) %>% as.data.frame()
  colnames(sp_occ_dx_200) <- c("mean", "lower", "upper")
  sp_occ_dx_200$species <- my_traits$species
  
  return(list(sp_occ_dx_100, sp_occ_dx_200, ave_site_Ind_100, ave_site_Ind_200))
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
  mean_dx <- compute_occ_shift(my_data_100_1, sim_matrix_100,
                               my_data_200_1, sim_matrix_200,
                               my_traits, site_type="mean")
  
  low_temp_dx <- compute_occ_shift(my_data_100_1, sim_matrix_100,
                                   my_data_200_1, sim_matrix_200,
                                   my_traits, site_type="low_temp_mean_precip")
  
  high_temp_dx <- compute_occ_shift(my_data_100_1, sim_matrix_100,
                                    my_data_200_1, sim_matrix_200,
                                    my_traits, site_type="high_temp_mean_precip")
  
  low_precip_dx <- compute_occ_shift(my_data_100_1, sim_matrix_100,
                                     my_data_200_1, sim_matrix_200,
                                     my_traits, site_type="mean_temp_low_precip")
  
  high_precip_dx <- compute_occ_shift(my_data_100_1, sim_matrix_100,
                                      my_data_200_1, sim_matrix_200,
                                      my_traits, site_type="mean_temp_high_precip")
  
  if(scale=="100"){
    # AVERAGE SITES IN 1970-1974
    cross_scale_occ <- mean_dx[[1]] %>%
      dplyr::select(species, mean, lower, upper)
    
    # LOW TEMPERATURE SITES IN 1970-1974
    cross_scale_occ_lowTemp <- low_temp_dx[[1]] %>%
      dplyr::select(species, mean, lower, upper)
    
    # HIGH TEMPERATURE SITES IN 1970-1974
    cross_scale_occ_highTemp <- high_temp_dx[[1]] %>%
      dplyr::select(species, mean, lower, upper)
    
    # LOW PRECIPITATION SITES IN 1970-1974
    cross_scale_occ_lowPrecip <- low_precip_dx[[1]] %>%
      dplyr::select(species, mean, lower, upper)
    
    # HIGH PRECIPITATION SITES IN 1970-1974
    cross_scale_occ_highPrecip <- high_precip_dx[[1]] %>%
      dplyr::select(species, mean, lower, upper)
  } else{
    # AVERAGE SITES IN 1970-1974
    cross_scale_occ <- mean_dx[[2]] %>%
      dplyr::select(species, mean, lower, upper)
    
    # LOW TEMPERATURE SITES IN 1970-1974
    cross_scale_occ_lowTemp <- low_temp_dx[[2]] %>%
      dplyr::select(species, mean, lower, upper)
    
    # HIGH TEMPERATURE SITES IN 1970-1974
    cross_scale_occ_highTemp <- high_temp_dx[[2]] %>%
      dplyr::select(species, mean, lower, upper)
    
    # LOW PRECIPITATION SITES IN 1970-1974
    cross_scale_occ_lowPrecip <- low_precip_dx[[2]] %>%
      dplyr::select(species, mean, lower, upper)
    
    # HIGH PRECIPITATION SITES IN 1970-1974
    cross_scale_occ_highPrecip <- high_precip_dx[[2]] %>%
      dplyr::select(species, mean, lower, upper)
  }
  
  
  # TEMPERATURE
  temp_occ_dx <- data.frame(species=cross_scale_occ$species,
                            meanOcc_mean=cross_scale_occ$mean, lowerOcc_mean=cross_scale_occ$lower, upperOcc_mean=cross_scale_occ$upper,
                            meanOcc_lowTemp=cross_scale_occ_lowTemp$mean, lowerOcc_lowTemp=cross_scale_occ_lowTemp$lower, upperOcc_lowTemp=cross_scale_occ_lowTemp$upper,
                            meanOcc_highTemp=cross_scale_occ_highTemp$mean, lowerOcc_highTemp=cross_scale_occ_highTemp$lower, upperOcc_highTemp=cross_scale_occ_highTemp$upper) %>%
    dplyr::arrange(meanOcc_mean) %>%
    dplyr::mutate(order=row_number()) %>%
    inner_join(my_traits, by="species") %>%
    dplyr::mutate(trend=ifelse(sign(meanOcc_mean)==-1, "Decreasing", "Increasing")) %>%
    dplyr::mutate(trend=ifelse(between(0, lowerOcc_mean, upperOcc_mean), "Stable", trend))
  temp_occ_dx$trend <- factor(temp_occ_dx$trend,
                              levels=c("Decreasing", "Stable", "Increasing"))
  
  temp_occ_dx_long <- temp_occ_dx %>%
    dplyr::select(species,  meanOcc_lowTemp, meanOcc_mean, meanOcc_highTemp) %>%
    pivot_longer(-c("species")) %>%
    inner_join(my_traits, by="species")
  temp_occ_dx_long$name <- factor(temp_occ_dx_long$name,
                                  levels=c("meanOcc_lowTemp",
                                           "meanOcc_mean",
                                           "meanOcc_highTemp"))
  
  temp_anova <- aov(value~name, dplyr::filter(temp_occ_dx_long,
                                              name %in% c("meanOcc_lowTemp",
                                                          "meanOcc_mean",
                                                          "meanOcc_highTemp")))
  summary(temp_anova)
  
  figure_3_panel_a <- ggplot()+
    geom_vline(xintercept=0, linetype=2, alpha=0.5)+
    geom_pointrange(temp_occ_dx,
                    mapping=aes(y=order, 
                                x=meanOcc_mean, xmin=lowerOcc_mean, xmax=upperOcc_mean, 
                                group="species", color=as.factor(trend)))+
    scale_color_manual(labels=c("Decreasing", "Stable", "Increasing"),
                       values=c("#f8a29e", "black", "#7fbff5"),
                       guide=guide_legend(direction="Vertical",
                                          title.position="top"),
                       name="Average Site Trend")+
    scale_y_continuous(breaks=c(1:nrow(my_traits)), labels=temp_occ_dx$speciesCode)+
    scale_x_continuous(labels=scales::percent)+
    labs(y="Species Code", x="Occupancy Probability Shift from \n1970-1974 to 2015-2019")+
    theme_cowplot()+
    theme(legend.position=c(0.55, 0.1),
          axis.text.y=element_text(size=9),
          legend.box="vertical",
          legend.key.width=unit(1, 'cm'),
          legend.background=element_rect(fill=alpha("white", 0.9)),
          plot.background=element_rect(fill="white"))
  
  figure_3_panel_b <- ggplot()+
    geom_hline(yintercept=0, linetype=2, alpha=0.25)+
    geom_rect(mapping=aes(xmin=0.75, xmax=1.25, 
                          ymin=min(dplyr::filter(temp_occ_dx_long, name=="meanOcc_lowTemp")$value), 
                          ymax=max(dplyr::filter(temp_occ_dx_long, name=="meanOcc_lowTemp")$value)), 
              color="grey",
              alpha=0.2)+
    geom_rect(mapping=aes(xmin=1.75, xmax=2.25, 
                          ymin=min(dplyr::filter(temp_occ_dx_long, name=="meanOcc_mean")$value), 
                          ymax=max(dplyr::filter(temp_occ_dx_long, name=="meanOcc_mean")$value)), 
              color="grey",
              alpha=0.2)+
    geom_rect(mapping=aes(xmin=2.75, xmax=3.25, 
                          ymin=min(dplyr::filter(temp_occ_dx_long, name=="meanOcc_highTemp")$value), 
                          ymax=max(dplyr::filter(temp_occ_dx_long, name=="meanOcc_highTemp")$value)), 
              color="grey",
              alpha=0.2)+
    geom_line(temp_occ_dx_long,
              mapping=aes(x=name, y=value, 
                          group=species, 
                          color=ave_temp2), alpha=0.7)+
    geom_point(temp_occ_dx_long,
               mapping=aes(x=name, y=value, group=species,
                           color=ave_temp2), alpha=0.7)+
    scale_color_continuous_divergingx(palette="Temps",
                                      mid=mean(my_traits$ave_temp2),
                                      guide=guide_colorbar(
                                        direction = "horizontal",
                                        title.position = "top"
                                      ), name="Average Annual Range-\nwide Temperature (°C)")+
    labs(x="", y="Occupancy Probability Shift")+
    scale_x_discrete(labels=c("Low Temp.", "Average Temp.", "High Temp."))+
    scale_y_continuous(labels=scales::percent)+
    theme_cowplot()+
    theme(axis.text=element_text(size=14),
          legend.position="none",
          legend.box="vertical",
          legend.key.width=unit(1, 'cm'),
          plot.background=element_rect(fill="white"))
  
  # Precipitation 
  precip_occ_dx <- data.frame(species=cross_scale_occ$species,
                              meanOcc_mean=cross_scale_occ$mean, lowerOcc_mean=cross_scale_occ$lower, upperOcc_mean=cross_scale_occ$upper,
                              meanOcc_lowPrecip=cross_scale_occ_lowPrecip$mean, lowerOcc_lowPrecip=cross_scale_occ_lowPrecip$lower, upperOcc_lowPrecip=cross_scale_occ_lowPrecip$upper,
                              meanOcc_highPrecip=cross_scale_occ_highPrecip$mean, lowerOcc_highPrecip=cross_scale_occ_highPrecip$lower, upperOcc_highPrecip=cross_scale_occ_highPrecip$upper) %>%
    dplyr::arrange(meanOcc_mean) %>%
    dplyr::mutate(order=row_number())
  
  precip_occ_dx_long <- precip_occ_dx %>%
    dplyr::select(species,  meanOcc_lowPrecip, meanOcc_mean, meanOcc_highPrecip) %>%
    pivot_longer(-c("species")) %>%
    inner_join(my_traits, by="species")
  precip_occ_dx_long$name <- factor(precip_occ_dx_long$name,
                                    levels=c("meanOcc_lowPrecip",
                                             "meanOcc_mean",
                                             "meanOcc_highPrecip"))
  
  precip_anova <- aov(value~name, dplyr::filter(precip_occ_dx_long,
                                                name %in% c("meanOcc_lowPrecip",
                                                            "meanOcc_mean",
                                                            "meanOcc_highPrecip")))
  summary(precip_anova)
  
  figure_3_panel_c <- ggplot()+
    geom_hline(yintercept=0, linetype=2, alpha=0.25)+
    geom_rect(mapping=aes(xmin=0.75, xmax=1.25, 
                          ymin=min(dplyr::filter(precip_occ_dx_long, name=="meanOcc_lowPrecip")$value), 
                          ymax=max(dplyr::filter(precip_occ_dx_long, name=="meanOcc_lowPrecip")$value)), 
              color="grey",
              alpha=0.2)+
    geom_rect(mapping=aes(xmin=1.75, xmax=2.25, 
                          ymin=min(dplyr::filter(precip_occ_dx_long, name=="meanOcc_mean")$value), 
                          ymax=max(dplyr::filter(precip_occ_dx_long, name=="meanOcc_mean")$value)), 
              color="grey",
              alpha=0.2)+
    geom_rect(mapping=aes(xmin=2.75, xmax=3.25, 
                          ymin=min(dplyr::filter(precip_occ_dx_long, name=="meanOcc_highPrecip")$value), 
                          ymax=max(dplyr::filter(precip_occ_dx_long, name=="meanOcc_highPrecip")$value)), 
              color="grey",
              alpha=0.2)+
    geom_line(precip_occ_dx_long,
              mapping=aes(x=name, y=value, 
                          group=species, 
                          color=ave_precip2), alpha=0.7)+
    geom_point(precip_occ_dx_long,
               mapping=aes(x=name, y=value, group=species,
                           color=ave_precip2), alpha=0.7)+
    scale_color_continuous_divergingx(palette="Earth",
                                      mid=mean(my_traits$ave_precip2),
                                      guide=guide_colorbar(
                                        direction = "horizontal",
                                        title.position = "top"
                                      ), name="Average Annual Range-\nwide Precipitation (cm)")+
    labs(x="In-range Site Type", y="Occupancy Probability Shift")+
    scale_x_discrete(labels=c("Low Prec.", "Average Prec.", "High Prec."))+
    scale_y_continuous(labels=scales::percent)+
    theme_cowplot()+
    theme(axis.text=element_text(size=14),
          legend.position="none",
          legend.box="vertical",
          legend.key.width=unit(1, 'cm'),
          plot.background=element_rect(fill="white"))
  
  figure_3_right <- cowplot::plot_grid(figure_3_panel_b, figure_3_panel_c,
                                       labels=c("(b)", "(c)"), nrow=2)
  figure_3 <- cowplot::plot_grid(figure_3_panel_a, figure_3_right, ncol=2,
                                 labels=c("(a)", ""), rel_widths=c(1,1))
  
  ggsave(paste0("../figures/FIGURE_3_", scale, ".png"), 
         figure_3, dpi=400, height=10, width=10)
  
  return(list(temp_occ_dx, precip_occ_dx))
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