test <- readRDS("../output/data/grid_50_occur.rds") %>%
  sf::st_drop_geometry()
saveRDS(test, "../output/data/grid_50_occur.rds")