memory.limit(50000)

library(tidyverse); library(rjags); library(jagsUI)
source("002_PrepData.R")
source("003_ModelSpec.R")

# Load the pruned tree
my_tree <- readRDS("output/tree_vcv.rds")

# 200 x 200 km analysis
my_data_200_1 <- make.data(my_tree, 200, 1)
my_res_200_1 <- run.R2jags.model(my_data_200_1)
saveRDS(my_res_200_1, "output/samples/res_200_1.rds")

my_data_200_2 <- make.data(my_tree, 200, 2)
my_res_200_2 <- run.R2jags.model(my_data_200_2)
saveRDS(my_res_200_2, "output/samples/res_200_2.rds")

# 100 x 100 km analysis
my_data_100_1 <- make.data(my_tree, 100, 1)
my_res_100_1 <- run.R2jags.model(my_data_100_1)
saveRDS(my_res_100_1, "output/samples/res_100_1.rds")

my_data_100_2 <- make.data(my_tree, 100, 2)
my_res_100_2 <- run.R2jags.model(my_data_100_2)
saveRDS(my_res_100_2, "output/samples/res_100_2.rds")

# 50 x 50 km analysis
my_data_50_1 <- make.data(my_tree, 50, 1)
my_res_50_1 <- run.R2jags.model(my_data_50_1)
saveRDS(my_res_50_1, "output/samples/res_50_1.rds")

my_data_50_2 <- make.data(my_tree, 50, 2)
my_res_50_2 <- run.R2jags.model(my_data_50_2)
saveRDS(my_res_50_2, "output/samples/res_50_2.rds")

