# Load libraries
library(tidyverse)
library(cowplot)
library(ggdag)

# Create a dag outlining our analysis (and other predictors not measured)
occupancy_dag <- ggdag::dagify(y~b+d+e,
                               b~d+e,
                               d~e,
                               exposure="d",
                               outcome="y",
                               labels=c("y"="Occupancy Probability",
                                        "b"="Hostplant Presence",
                                        "d"="Temperature",
                                        "e"="Precipitation"))

ggdag_adjustment_set(occupancy_dag, text=FALSE, use_labels="label",
                     shadow=TRUE, stylized=TRUE)+
  theme_dag_blank()+
  theme(plot.background=element_rect(color="white", fill="white"),
        legend.position="none")