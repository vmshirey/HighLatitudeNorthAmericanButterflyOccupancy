# Load the required libraries
library(tidyverse)
library(data.table)
library(cowplot)
library(colorspace)

# Read in the number of days above 32-degrees and the heatwave data
heatwaves5 <- fread("../../data/climate/WarmSpellsAbove5C.csv")
warmdays5 <- fread("../../data/climate/WarmDaysAbove5C.csv")

heatwaves10 <- fread("../../data/climate/WarmSpellsAbove10C.csv")
warmdays10 <- fread("../../data/climate/WarmDaysAbove10C.csv")

heatwaves15 <- fread("../../data/climate/WarmSpellsAbove15C.csv")
warmdays15 <- fread("../../data/climate/WarmDaysAbove15C.csv")

# Create a new dataframe which includes one mean record of either
# heatwave or warm day incidence by latitude (averaged across longitudes)
# and year.
heatwaves10_long <- heatwaves10 %>%
  pivot_longer(cols=3:48) %>%
  dplyr::mutate(groupLoc=paste(Lat, Lon))

warmdays10_long <- warmdays10 %>%
  pivot_longer(cols=3:48) %>%
  dplyr::mutate(groupLoc=paste(Lat, Lon))

heatwaves10_sum <- heatwaves10 %>%
  pivot_longer(cols=3:48) %>%
  group_by(name) %>%
  dplyr::mutate(aveVal=mean(value)) %>%
  dplyr::mutate(name=as.numeric(name)) %>%
  dplyr::select(name, aveVal) %>%
  ungroup() %>%
  unique()

warmdays10_sum <- warmdays10 %>%
  pivot_longer(cols=3:48) %>%
  group_by(name) %>%
  dplyr::mutate(aveVal=mean(value)) %>%
  dplyr::mutate(name=as.numeric(name)) %>%
  dplyr::select(name, aveVal) %>%
  ungroup() %>%
  unique()

warmdays10_sum_lat <- warmdays10 %>%
  pivot_longer(cols=3:48) %>%
  group_by(name, Lat) %>%
  dplyr::mutate(aveVal=mean(value)) %>%
  dplyr::mutate(name=as.numeric(name)) %>%
  dplyr::select(name, Lat, aveVal) %>%
  ungroup() %>%
  unique()

ggplot()+
  # # 5-Celsius
  # geom_smooth(warmdays5_sum,
  #             mapping=aes(x=name,
  #                         y=aveVal),
  #             color="#80cdc1",
  #             fill="#80cdc1")+
  # geom_line(warmdays5_sum,
  #           mapping=aes(x=name,
  #                       y=zoo::rollmean(aveVal, 5, na.pad=TRUE)),
  #           size=1)+
  # geom_point(warmdays5_sum,
  #            mapping=aes(x=name,
  #                        y=aveVal),
  #            color="#80cdc1",
  #            size=3)+
  # 10-Celsius
  geom_smooth(warmdays10_sum,
              mapping=aes(x=name,
                          y=aveVal),
              color="#f3b300",
              fill="#f3b300",
              method="lm")+
  geom_line(warmdays10_sum,
            mapping=aes(x=name,
                        y=zoo::rollmean(aveVal, 5, na.pad=TRUE)),
            size=1)+
  geom_point(warmdays10_sum,
             mapping=aes(x=name,
                         y=aveVal),
             color="#f3b300",
             size=3)+
  # # 15-Celsius
  # geom_smooth(warmdays15_sum,
  #             mapping=aes(x=name,
  #                         y=aveVal),
  #             color="#ec6626",
  #             fill="#ec6626")+
  # geom_line(warmdays15_sum,
  #           mapping=aes(x=name,
  #                       y=zoo::rollmean(aveVal, 5, na.pad=TRUE)),
  #           size=1)+
  # geom_point(warmdays15_sum,
  #            mapping=aes(x=name,
  #                        y=aveVal),
  #            color="#ec6626",
  #            size=3)+
  xlab("Year")+
  ylab("Ave. Abnormal Warm Day Incidence")+
  theme_cowplot()+
  theme(axis.text.x=element_text(angle=45, hjust=1),
        plot.background=element_rect(fill="white"))
ggsave2("../../figures/supplemental/warmday_incidence.png", dpi=400, height=5, width=6)

ggplot()+
  # # 5-Celsius
  # geom_smooth(warmdays5_sum,
  #             mapping=aes(x=name,
  #                         y=aveVal),
  #             color="#80cdc1",
  #             fill="#80cdc1")+
  # geom_line(warmdays5_sum,
  #           mapping=aes(x=name,
  #                       y=zoo::rollmean(aveVal, 5, na.pad=TRUE)),
  #           size=1)+
  # geom_point(warmdays5_sum,
#            mapping=aes(x=name,
#                        y=aveVal),
#            color="#80cdc1",
#            size=3)+
# 10-Celsius
geom_smooth(heatwaves10_sum,
            mapping=aes(x=name,
                        y=aveVal),
            color="#f3b300",
            fill="#f3b300",
            method="lm")+
  geom_line(heatwaves10_sum,
            mapping=aes(x=name,
                        y=zoo::rollmean(aveVal, 5, na.pad=TRUE)),
            size=1)+
  geom_point(heatwaves10_sum,
             mapping=aes(x=name,
                         y=aveVal),
             color="#f3b300",
             size=3)+
  # # 15-Celsius
  # geom_smooth(warmdays15_sum,
  #             mapping=aes(x=name,
  #                         y=aveVal),
  #             color="#ec6626",
  #             fill="#ec6626")+
  # geom_line(warmdays15_sum,
  #           mapping=aes(x=name,
  #                       y=zoo::rollmean(aveVal, 5, na.pad=TRUE)),
  #           size=1)+
  # geom_point(warmdays15_sum,
#            mapping=aes(x=name,
#                        y=aveVal),
#            color="#ec6626",
#            size=3)+
xlab("Year")+
  ylab("Ave. Abnormal Warm Day Incidence")+
  theme_cowplot()+
  theme(axis.text.x=element_text(angle=45, hjust=1),
        plot.background=element_rect(fill="white"))
ggsave2("../../figures/supplemental/heatwave_incidence.png", dpi=400, height=6, width=6)
