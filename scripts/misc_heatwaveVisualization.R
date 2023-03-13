# Load the required libraries
library(tidyverse)
library(data.table)
library(cowplot)
library(colorspace)

# Read in the number of days above 32-degrees and the heatwave data
heatwaves <- fread("../../data/climate/HeatWaves.csv")