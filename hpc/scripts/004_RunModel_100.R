library(tidyverse); library(rjags); library(jagsUI); library(parallel)
source("scripts/002_PrepData.R")
source("scripts/003_ModelSpec.R")

# Grab the core IDs
MC_CORES <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
options(mc.cores = MC_CORES)

# Load the pruned tree
my_tree <- readRDS("output/tree_vcv.rds")

# 100 x 100 km analysis
my_data_100_1 <- make.data(my_tree, 100, 1)

my_res_100_1 <- mclapply(list(my_data_100_1,my_data_100_1,my_data_100_1), run.R2jags.model, mc.cores=MC_CORES)

saveRDS(my_res_100_1, "output/samples/res_100_1.rds")
