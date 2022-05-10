memory.limit(50000)

library(tidyverse); library(rjags); library(jagsUI)
source("scripts/002_PrepData.R")
source("scripts/003_ModelSpec.R")

# Load the pruned tree
my_tree <- readRDS("output/tree_vcv.rds")

# 100 x 100 km analysis
my_data_100_1 <- make.data(my_tree, 100, 1)
my_res_100_1 <- run.R2jags.model(my_data_100_1)
saveRDS(my_res_100_1, "output/samples/res_100_1.rds")

my_data_100_2 <- make.data(my_tree, 100, 2)
my_res_100_2 <- run.R2jags.model(my_data_100_2)
saveRDS(my_res_100_2, "output/samples/res_100_2.rds")



