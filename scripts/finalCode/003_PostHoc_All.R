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
my_res_100_1 <- readRDS("my_res_100_1_all.RDS")
MCMCvis::MCMCtrace(my_res_100_1, Rhat=TRUE, filename="../../figures/supplemental/SupplementalFile_S5.pdf")
my_sum_100_1 <- MCMCvis::MCMCsummary(my_res_100_1)
my_sim_mat_100_1 <- as.matrix(my_res_100_1)

covariate_terms <- colnames(my_sim_mat_100_1)[grepl("psi.beta", colnames(my_sim_mat_100_1))]
my_tolerances <- data.frame(term=covariate_terms,
                            tol=NA)
for(x in 1:nrow(my_tolerances)){
  
  focalx <- my_tolerances$term[x]
  SPID <- grep("(?<=\\[)([0-9]{1,2})(?=])", focalx, perl=TRUE)
  focaly <- my_tolerances$term[grepl(paste0("(?<=\\[)(", SPID, ")(?=])"), my_tolerances$term,
                                     perl=TRUE)]
  focaly <- focaly[focaly %!in% focalx]
  focaly <- paste0("`", focaly, "`")
  
  mx <- lm(as.formula(paste0("`", focalx, "`", "~",
                             paste(focaly, collapse=" + "))),
           data=as.data.frame(my_sim_mat_100_1))
  my_tolerances[x, "tol"] <- 1/(1-summary(mx)$r.squared)
}

ggplot()+
  geom_histogram(my_tolerances,
                 mapping=aes(x=tol),
                 binwidth=0.05)+
  geom_vline(xintercept=5, linetype=2)+
  theme_cowplot()+
  theme(plot.background=element_rect(fill="white", color="white"))
