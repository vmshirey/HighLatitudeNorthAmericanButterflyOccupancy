library(tidyverse); library(rjags); library(jagsUI); library(parallel)
source("002_PrepData.R")
source("003_ModelSpec.R")

# Load the pruned tree
my_tree <- readRDS("../output/tree_vcv.rds")

# 200 x 200 km analysis
my_data_200_1 <- make.data(my_tree, 200, 1)

my_res_200_1 <- run.R2jags.model(d=my_data_200_1, ni=5.1e4, nb=1e3, nt=100, nc=3, phy=FALSE)
saveRDS(my_res_200_1, "../output/samples/my_res_200_1_LONGRUN_Intercept.RDS")

# my_res_200_1_phy <- run.R2jags.model(d=my_data_200_1, ni=8.1e4, nb=1e3, nt=200, nc=4, phy=TRUE)
# saveRDS(my_res_200_1_phy, "../output/samples/my_res_200_1_phy_LONGRUN_MAINE.RDS")

# 100 x 100 km analysis
my_data_100_1 <- make.data(my_tree, 100, 1)

my_res_100_1 <- run.R2jags.model(d=my_data_100_1, ni=5.1e4, nb=1e3, nt=100, nc=3, phy=FALSE)
saveRDS(my_res_100_1, "../output/samples/my_res_100_1_LONGRUN_Intercept.RDS")

# my_res_100_1_phy <- run.R2jags.model(d=my_data_100_1, ni=8.1e4, nb=1e3, nt=200, nc=4, phy=TRUE)
# saveRDS(my_res_100_1_phy, "../output/samples/my_res_100_1_phy_LONGRUN_MAINE.RDS")