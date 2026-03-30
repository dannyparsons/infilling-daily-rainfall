library(here)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

source(here("src", "helper_funs.R"))

# Information http://www.meteotemplate.com/template/plugins/climateClassification/koppen.php
# Data from http://hanschen.org/koppen

koppen <- read.table(here("data", "koppen_1901-2010.tsv"), header = TRUE)

koppen_zim <- koppen %>%
  filter(longitude >= 25 & longitude <= 34 & latitude  >= -23 & latitude  <= -15)

zim <- ne_countries(scale = "medium", country = "Zimbabwe", returnclass = "sf")

ggplot() +
  geom_tile(data = koppen_zim,
            aes(x = longitude, y = latitude, fill = p1901_2010)) +
  geom_sf(data = zim, fill = NA, color = "black", linewidth = 0.7) +
  coord_sf(xlim = c(25, 34), ylim = c(-23, -15), expand = FALSE) +
  scale_fill_discrete(drop = TRUE) +
  theme_minimal() +
  labs(fill = "Koppen class")

zimbabwe_stations <- readr::read_csv(here("data", "zimbabwe_stations.csv"))

zimbabwe_stations$station <- gsub("_", " ", zimbabwe_stations$station)

grd_pts <- closest_point(points = koppen %>% dplyr::select(longitude, latitude), 
                         target = zimbabwe_stations %>% dplyr::select(lon, lat))

grd_pts$station <- zimbabwe_stations$station

grd_pts <- left_join(grd_pts, koppen, by = c("longitude", "latitude"))
