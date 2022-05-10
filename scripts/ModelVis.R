# Load libraries
library(tidyverse); library(cowplot); library(MCMCvis)

# 100 by 100 km analysis
my_res_100_1 <- readRDS("../output/samples/hpc/res_100_1.rds")
my_res_100_2 <- readRDS("../output/samples/hpc/res_100_2.rds")

# 200 by 200 km analysis
my_res_200_1 <- readRDS("../output/samples/hpc/res_200_1.rds")
my_res_200_2 <- readRDS("../output/samples/hpc/res_200_2.rds")



