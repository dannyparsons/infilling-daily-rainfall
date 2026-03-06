library(here)
library(dplyr)

source(here("src", "helper_funs.R"))

# Information http://www.meteotemplate.com/template/plugins/climateClassification/koppen.php
# Data from http://hanschen.org/koppen

koppen <- read.table(here("data", "koppen_1901-2010.tsv"), header = TRUE)

zimbabwe_stations <- readr::read_csv(here("data", "zimbabwe_stations.csv"))

zimbabwe_stations$station <- gsub("_", " ", zimbabwe_stations$station)

grd_pts <- closest_point(points = koppen %>% dplyr::select(longitude, latitude), 
                         target = zimbabwe_stations %>% dplyr::select(lon, lat))

grd_pts$station <- zimbabwe_stations$station

grd_pts <- left_join(grd_pts, koppen, by = c("longitude", "latitude"))
