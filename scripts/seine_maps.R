library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(tidyverse)
library(cowplot)
library(ggpubr)
library(brms)
library(rgdal)

theme_set(theme_bw())

# set palette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


## cod map -----------------------------------------------

# load data sets
bays <- read.csv("./data/bay_lat_long.csv", row.names = 1)

# remove Cooks (repeat) and Kujulik (not resampled) and Ugak (not sampled in 2023)
drop <- bays$Bay %in% c("Cooks", "Kujulik", "Ugak")
bays <- bays[!drop,]

ak <- ne_countries(scale = "large", returnclass = "sf", continent="north america")

# use this version unless the high-res version is entered!
# ak <- ne_countries(scale = "medium", returnclass = "sf", continent="north america")
world <- ne_countries(scale='medium', returnclass = "sf")

box <- data.frame(long = c(-163, -163, -151, -151, -163), lat = c(54.5, 59.2, 59.2, 54.5, 54.5))


map.plot <- ggplot(ak) +
  geom_sf(fill="darkgoldenrod3", color=NA) + 
  coord_sf(xlim = c(-163, -151), ylim = c(54.5, 59.5), expand = FALSE) +
  geom_point(data = bays, aes(-lon, lat), fill = cb[4], size=3, shape=21) +
  theme(axis.title = element_blank()) +
  scale_fill_manual(values=cb[4]) +
  scale_x_continuous(breaks = c(-160, -156, -152)) +
  scale_y_continuous(breaks = c(55, 56, 57, 58, 59))

map.plot


ggsave("./figs/cod_map.png", width = 4, height = 3, units = 'in')

## pollock map -----------------------------------------------

# load data sets
bays <- read.csv("./data/bay_lat_long.csv", row.names = 1)

# remove Cooks (repeat) and Kujulik (not resampled) and Ugak (not sampled in 2023)
keep <- c("Agripina", "Anton Larson Bay", "Balboa", "Cook Bay", "Mitrofania", "Port Wrangell") 

bays <- bays %>%
  filter(Bay %in% keep)

ak <- ne_countries(scale = "large", returnclass = "sf", continent="north america")


map.plot <- ggplot(ak) +
  geom_sf(fill="darkgoldenrod3", color=NA) + 
  coord_sf(xlim = c(-163, -151), ylim = c(54.5, 59.5), expand = FALSE) +
  geom_point(data = bays, aes(-lon, lat), fill = cb[4], size=3, shape=21) +
  theme(axis.title = element_blank()) +
  scale_fill_manual(values=cb[4]) +
  scale_x_continuous(breaks = c(-160, -156, -152)) +
  scale_y_continuous(breaks = c(55, 56, 57, 58, 59))

map.plot


ggsave("./figs/pollock_map.png", width = 4, height = 3, units = 'in')
