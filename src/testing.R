library(here)
library(readr)
library(dplyr)
library(fitdistrplus)

set.seed(12)

x <- c(0, 0, 5, 6, 2, NA, 2, 7, 0.8, 0.2, 0.9, 6, 30, 50, NA, 3, 0, 0, NA)
y <- c(0, 5, 7, 7, 1, 8,  2, 8, 1.8, 6.5, 4.1, 9, 12, 20, 3,  8, 2, 8, 25)

source(here("src", "methods.R"))

thresh(x, y)
loci(x, y)

rwanda <- readr::read_csv(here("data", "rwanda_daily.csv"))
kigali <- rwanda %>%
  filter(station_name == "KIGALI AERO") %>%
  mutate(rain_sim = precip + rnorm(n(), 4, 4),
         rain_sim = ifelse(rain_sim < 0, 0, rain_sim))

thresh_calc(kigali$precip, kigali$rain_sim)
