library(sf)
library(ggplot2)
library(ggrepel)
library(ggspatial)
library(rnaturalearth)
library(dplyr)
library(readr)
library(here)

source(here("src", "helper_funs.R"))

zimbabwe_stations <- readr::read_csv(here("data", "zimbabwe_stations.csv"))

zimbabwe_stations$station <- gsub("_", " ", zimbabwe_stations$station)

world <- ne_countries(scale = "medium", returnclass = "sf")

zimbabwe <- world %>% filter(admin == "Zimbabwe")

neighbors <- world[st_touches(zimbabwe, world, sparse = FALSE), ]
neighbors <- rbind(zimbabwe, neighbors)

neighbor_labels <- data.frame(
  admin = c("Zimbabwe", "Botswana", "Mozambique", "South Africa", "Zambia", "Namibia"),
  lon   = c(30, 25.5, 34, 30, 27.5, 24.5),
  lat   = c(-19, -20, -19, -23, -16, -17.7)
)

neighbor_labels_sf <- st_as_sf(
  neighbor_labels,
  coords = c("lon", "lat"),
  crs = 4326
)

stations_sf <- st_as_sf(
  zimbabwe_stations,
  coords = c("lon", "lat"),
  crs = 4326
)

ggplot() +
  geom_sf(data = neighbors,
          fill = "grey97",
          color = "grey60",
          linewidth = 0.4) +
  geom_sf(data = zimbabwe, fill = "antiquewhite",
          color = "black", linewidth = 0.6) +
  geom_sf(data = stations_sf, color = "#0B1F5C", size = 2.8) +
  geom_label_repel(data = zimbabwe_stations,
    aes(x = lon, y = lat, label = station),
    size = 5, box.padding = 0.35, point.padding = 0.3, 
    colour = "#0B1F5C",
    segment.color = "grey40") +
  geom_sf_text(
    data = neighbor_labels_sf,
    aes(label = admin),
    size = 3.5,
  ) +
  coord_sf(xlim = c(24.5, 34.5), ylim = c(-23, -15.5)) +
  annotation_scale(location = "bl", width_hint = 0.28) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  base_theme() +
  labs(x = "Longitude", y = "Latitude")

ggsave(here("results", "Fig1.jpeg"),
       width = 12, height = 6)


# Station information -----------------------------------------------------

zimbabwe_bc <- read_rds(here("data", "BC_data", "zimbabwe_agera5_bc.RDS"))
zimbabwe_bc$station <- gsub("_", " ", zimbabwe_bc$station)

zimbabwe_station_info <- zimbabwe_bc %>%
  group_by(station) %>%
  summarise(`Complete Days (%)` = mean(!is.na(rain)) * 100)

table1 <- left_join(zimbabwe_stations, zimbabwe_station_info)
names(table1)[1:3] <- c("Station", "Latitude", "Longitude")

table1 %>%
  mutate(Latitude = sprintf("%.2f", Latitude),
         Longitude = sprintf("%.2f", Longitude),
         `Complete Days (%)` = sprintf("%.1f", `Complete Days (%)`)) %>%
  write.csv(here("results", "Table1.csv"), row.names = FALSE)
