library(tidyverse)
library(rjags)
source("000_Initialization.R")
source("MISC_PrepModelData.R")
source("MISC_ModelSpecification.R")

# Read in the species trait data
sp_traits <- read.csv("../../data/taxa/species_traits.csv")

# Prepare the model-ready data
my_data_100_1 <- make.data(100)
my_data_200_1 <- make.data(200)

# Temperature analysis
my_res_100_1 <- run.R2jags.model(d=my_data_100_1, ni=1.5e5, nb=5e4, nt=100, nc=4, precip=FALSE)
saveRDS(my_res_100_1, "my_res_100_1_temp.RDS")
my_res_200_1 <- run.R2jags.model(d=my_data_200_1, ni=1.5e5, nb=5e4, nt=100, nc=4, precip=FALSE)
saveRDS(my_res_200_1, "my_res_200_1_temp.RDS")

# Precipitation Analysis
my_res_100_1 <- run.R2jags.model(d=my_data_100_1, ni=1.5e5, nb=5e4, nt=100, nc=4, precip=TRUE)
saveRDS(my_res_100_1, "my_res_100_1_precip.RDS")
my_res_200_1 <- run.R2jags.model(d=my_data_200_1, ni=1.5e5, nb=5e4, nt=100, nc=4, precip=TRUE)
saveRDS(my_res_200_1, "my_res_200_1_precip.RDS")

# Combined Analysis
my_res_100_1 <- run.R2jags.model(d=my_data_100_1, ni=1.5e5, nb=5e4, nt=100, nc=4, precip="all")
saveRDS(my_res_100_1, "my_res_100_1_all.RDS")
my_res_200_1 <- run.R2jags.model(d=my_data_200_1, ni=1.5e5, nb=5e4, nt=100, nc=4, precip="all")
saveRDS(my_res_200_1, "my_res_200_1_all.RDS")







