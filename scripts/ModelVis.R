# Load libraries
source("000_init.R"); library(colorspace)
source("002_PrepData.R")
source("misc_plotCode.R")

# Load species trait data
my_traits <- read.csv("../data/taxa/species_range_clim.csv")

# Load data that went into the model
my_data_100_1 <- make.data(scale=100, imputeThres=1) 
my_data_100_2 <- make.data(scale=100, imputeThres=2) 

my_data_200_1 <- make.data(scale=200, imputeThres=1) 

# 100 by 100 km analysis
my_res_100_1 <- readRDS("../output/samples/res_100_1.rds")
my_res_100_2 <- readRDS("../output/samples/res_100_2.rds")

# 200 by 200 km analysis
my_res_200_1 <- readRDS("../output/samples/hpc_05_2022/res_200_1.rds")
my_res_200_2 <- readRDS("../output/samples/res_200_2.rds")

# Plot effect lines as a forest plot
visFigureTwo(my_res_100_2, my_data_100_2, my_traits)
visSupplementalFigureTwo(my_res_100_2, my_data_100_2, my_traits)


