library(here)
library(dplyr)
library(ggplot2)
library(lubridate)
library(readr)

source(here("src", "helper_funs.R"))

# Read in the data
zim_five_stations <- read_csv(here("data", "zim_five_stations.csv"))

# Set up for QC cod
zim_five_stations$station <- factor(zim_five_stations$station)
zim_five_stations$date <- as.Date(zim_five_stations$date)

zim_five_stations <- zim_five_stations %>%
  dplyr::filter(date >= as.Date("1979-01-01") & date <= as.Date("2023-06-30"))

zim_five_stations <- zim_five_stations %>% 
  mutate(month = factor(lubridate::month(date)),
         year = lubridate::year(date))

# Rainfall amounts =======================================================================
# Yearly Rainfall
# rainy days (n_rain) and total rainfall (t_rain)
zim_five_stations_year <- zim_five_stations %>%
  group_by(station, year) %>%
  summarise(n_rain = sum(rain > 0.85),
            t_rain = sum(rain))

# Monthly-yearly rainfall - rain days and total rainfall
zim_five_stations_month <- zim_five_stations %>%
  group_by(station, year, month) %>%
  summarise(n_rain = sum(rain > 0.85),
            t_rain = sum(rain))

# Rain amounts
ggplot(zim_five_stations %>% filter(rain > 0), aes(x = month, y = rain)) +
  geom_boxplot() +
  facet_wrap(~station)

# Number rain days
ggplot(zim_five_stations_month, aes(x = month, y = n_rain)) +
  geom_boxplot() +
  facet_wrap(~station)

# Inventory plot =======================================================================
# There are 46 missing values.
# This includes a long period of NAs in December 1992 in Buffalo Range. Otherwise some random missing values throughout for Mt. Darwin (6 values), Buffalo Range, and Plumtree (3 values)
ggplot(zim_five_stations, aes(x = date, y = station, fill = !is.na(rain))) +
  geom_tile() +
  geom_hline(yintercept = seq(0.5, by = 1, length.out = length(unique(zim_five_stations$station)) + 1))

# Fill Date gaps =======================================================================
dates_list <- list()
for(s in unique(zim_five_stations$station)) {
  
  dates <- seq(min((zim_five_stations %>% filter(station == s))$date), 
               max((zim_five_stations %>% filter(station == s))$date),
               by = 1)
  dd <- data.frame(station = s, date = dates)
  dates_list[[length(dates_list) + 1]] <- dd
}
date_df <- bind_rows(dates_list)

nr <- nrow(zim_five_stations)
zim_five_stations <- full_join(date_df, zim_five_stations, by = c("station", "date"))
print(paste("Filled", nrow(zim_five_stations) - nr, "rows"))

# No date gaps found

# Large or negative values check =======================================================================

large_check <- zim_five_stations %>% 
  filter(rain < 0 | rain > 200)
if(nrow(large_check) > 0) View(large_check)

# Three instances, all within realistic range.

# 1. Buffalo Range February 2019
# Has a big values on the next day also, 68.
# Nearby stations also have high rainfall (Chisumbanje has 24.7mm that day; Masvingo has 20.1mm that day; 25.9mm the next).
# Cyclone Dineo hit Mozambique over that day: "Widespread flooding took place in Zimbabwe, with Mutare, Chiredzi, and Beitbridge particularly hard-hit.". Chiredzi is ~18km from Buffalo Range.

# 2. Chisumbanje, January 2005
# In the middle of a 4 day rain spell
# Nearby stations also have high rainfall (Buffalo Range 42.9mm; Masvingo 69.0mm)
# Potentially related to Storm Chedza - "heavy rainfall over Mozambique occurred between 11 and 13 January with some stations recording about 80% of their total January 2015 rainfall as resulting from this event."

# 3. Masvingo, March 2003
# In the middle of a big 6 day rain spell
# No large values in Buffalo Range or Chisumbanje around then but Mt. Darwin has 62.7mm rainfall.
# Other data sources show that Rupike and Zaka had rainfall >100mm that day, within 100km of Masvingo.
# Cyclone Japhet hit over that time, which seems to have been pretty devastating and was reported to be in the South of Zimbabwe.

# Consecutive non-zero values check =======================================================================
consec_check <- zim_five_stations %>% 
  group_by(station) %>%
  mutate(same = rep(rle(as.numeric(rain))$lengths, rle(as.numeric(rain))$lengths)) %>%
  filter(rain > 1.5 & same >= 2)
if(nrow(consec_check) > 0) View(consec_check)

# No instances more than 2 days in length.

# Consecutive rain days check =======================================================================
raindays_check <- zim_five_stations %>%
  group_by(station) %>%
  mutate(raindays = cumsum(rain > 0) - cummax(((rain > 0) == 0) * cumsum(rain > 0))) %>%
  filter(raindays > 15)
if(nrow(raindays_check) > 0) View(raindays_check)

# Consecutive rainy days > 15 days and <= 20 days all in the middle of the rainy season.
# No incorrect variable entry detected.

# Dry months check - strict =======================================================================

drymonths_check <- zim_five_stations %>%
  filter(!month %in% 4:10) %>%
  group_by(station, year, month) %>%
  summarise(t_rain = sum(rain)) %>%
  filter(t_rain == 0)
if(nrow(drymonths_check) > 0) View(drymonths_check)

# Some instances of 0mm total monthly rainfall during the rainy season.
# These values are suspected missing data incorrectly entered as 0s.
# These values are replaced with missing values unless there is strong evidence that the 0 values are correct.

# Major drought years in 1991-92, 1994-95, 2002-03, 2015-16
# from this: https://www.mdpi.com/2071-1050/12/3/752#sustainability-12-00752-f003
# and supported here: https://www.weatherzw.org.zw/news/drought-occurrence-in-zimbabwe/
# Buffalo Range, Chisumbanje, and Masvingo are not *too* vulnerable.

# 1. Buffalo Range, 1992
# No rain in February.
# Trace rainfall only from 24th Jan - 12th March.
# Drought was reported / Severe El Nino year [https://www.weatherzw.org.zw/news/drought-occurrence-in-zimbabwe/]
# From this image it looks like a correct drought period https://www.mdpi.com/2071-1050/12/3/752#sustainability-12-00752-f003
display_daily(zim_five_stations %>% filter(station == "Buffalo_Range" & year == 1992),
              Stations = "Buffalo_Range", Years = 1992, Variables = "rain")
# Values will be left as 0s.

# 2. Buffalo Range, 2004 and 2021: Both had no rain in March

# 2004: Rainfall at end of February and beginning of April.
# 2004: Not a major drought year.
# 2004: Set March as NAs
zim_five_stations <- zim_five_stations %>%
  dplyr::mutate(rain = ifelse(station == "Buffalo_Range" & year == 2004 & month == 3,
                              NA,
                              rain))

# 2021: Rainfall at end of February
# 2021: Not a major drought year.
# 2021: Set March as NAs.
zim_five_stations <- zim_five_stations %>%
  dplyr::mutate(rain = ifelse(station == "Buffalo_Range" & year == 2021 & month == 3,
                              NA,
                              rain))

# 3. Chisumbanje, 1982, 2005, 2008, 2011: No rain in March

# 1982: Rainfall in February.
# 1982: Not a major drought year.
# 1982: Set March as NAs.
display_daily(zim_five_stations %>% filter(station == "Chisumbanje" & year == 1982),
              Stations = "Chisumbanje", Years = 1982, Variables = "rain")

zim_five_stations <- zim_five_stations %>%
  dplyr::mutate(rain = ifelse(station == "Chisumbanje" & year == 1982 & month == 3,
                              NA,
                              rain))

# 2005: Rainfall in February.
# 2005: Not a major drought year.
# 2005: Set March as NAs.
display_daily(zim_five_stations %>% filter(station == "Chisumbanje" & year == 2005),
              Stations = "Chisumbanje", Years = 2005, Variables = "rain")

zim_five_stations <- zim_five_stations %>%
  dplyr::mutate(rain = ifelse(station == "Chisumbanje" & year == 2005 & month == 3,
                              NA,
                              rain))

# 2008: Little rainfall in February.
# 2008: Not a major drought year.
# 2008: Set March as NAs.
display_daily(zim_five_stations %>% filter(station == "Chisumbanje" & year == 2008),
              Stations = "Chisumbanje", Years = 2008, Variables = "rain")

zim_five_stations <- zim_five_stations %>%
  dplyr::mutate(rain = ifelse(station == "Chisumbanje" & year == 2008 & month == 3,
                              NA,
                              rain))

# 2011: Little rainfall in February.
# 2011: Not a major drought year.
# 2011: Set March as NAs.
display_daily(zim_five_stations %>% filter(station == "Chisumbanje" & year == 2011),
              Stations = "Chisumbanje", Years = 2011, Variables = "rain")

zim_five_stations <- zim_five_stations %>%
  dplyr::mutate(rain = ifelse(station == "Chisumbanje" & year == 2011 & month == 3,
                              NA,
                              rain))

# 4. Chisumbanje, 2002-2003; 2009 March & 2009-2010: 0s reported for entire rainy season October - March

# 2002-2003: Is drought year, but no reported rainfall in entire rainy season is likely incorrectly entered missing values.
# 2002-2003: Set all as NA in 2002-2003 rainy season October - March
display_daily(zim_five_stations %>% filter(station == "Chisumbanje" & year %in% 2002:2003),
              Stations = "Chisumbanje", Years = 2002:2003, Variables = "rain")

zim_five_stations <- zim_five_stations %>%
  dplyr::mutate(rain = ifelse(station == "Chisumbanje" & year == 2002 & month %in% 10:12, NA, rain)) %>%
  dplyr::mutate(rain = ifelse(station == "Chisumbanje" & year == 2003 & month %in% 1:3, NA, rain))
  
# 2009 March & 2009-2010: Not a major drought year.
# 2009 March & 2009-2010: No reported in March 2009 and no rainfall in entire 2009-2010 rainy season.
# 2009 March & 2002-2003: Set all as NA from March 2009 - March 2010 as all values in this period are suspected incorrectly entered missing values.
display_daily(zim_five_stations %>% filter(station == "Chisumbanje" & year %in% 2009:2010),
              Stations = "Chisumbanje", Years = 2009:2010, Variables = "rain")

zim_five_stations <- zim_five_stations %>%
  dplyr::mutate(rain = ifelse(station == "Chisumbanje" & year == 2009 & month %in% 3:12, NA, rain)) %>%
  dplyr::mutate(rain = ifelse(station == "Chisumbanje" & year == 2010 & month %in% 1:3, NA, rain))

# 5. Chisumbanje: 2015 December
# 2015: Is a drought year, but some rainfall in November and January.
# 2015: Set December as NAs.

display_daily(zim_five_stations %>% filter(station == "Chisumbanje" & year %in% 2015:2016),
              Stations = "Chisumbanje", Years = 2015:2016, Variables = "rain")

zim_five_stations <- zim_five_stations %>%
  dplyr::mutate(rain = ifelse(station == "Chisumbanje" & year == 2015 & month == 12, NA, rain))

# 6. Masvingo: 2004 March
# 2004: Not a major drought year.
# 2004: Rainfall at end of February and beginning of April.
# 2004: Set March as NAs.
display_daily(zim_five_stations %>% filter(station == "Masvingo" & year == 2004),
              Stations = "Masvingo", Years = 2004, Variables = "rain")

zim_five_stations <- zim_five_stations %>%
  dplyr::mutate(rain = ifelse(station == "Masvingo" & year == 2004 & month == 3, NA, rain))

# 7. Mt Darwin: 1994 November
# 1994: Not a major drought year.
# 1994: Rainfall in October and December.
# 1994: Set November as NAs
display_daily(zim_five_stations %>% filter(station == "Mt_Darwin" & year == 1994),
              Stations = "Mt_Darwin", Years = 1994, Variables = "rain")

zim_five_stations <- zim_five_stations %>%
  dplyr::mutate(rain = ifelse(station == "Mt_Darwin" & year == 1994 & month == 11, NA, rain))

# 8. Mt Darwin: 2020 March
# 2020: Not a major drought year.
# 2020: Rainfall in February.
# 2020: Set March as NAs
display_daily(zim_five_stations %>% filter(station == "Mt_Darwin" & year == 2020),
              Stations = "Mt_Darwin", Years = 2020, Variables = "rain")

zim_five_stations <- zim_five_stations %>%
  dplyr::mutate(rain = ifelse(station == "Mt_Darwin" & year == 2020 & month == 3, NA, rain))

# 9. Plumtree: 2019 March
# 2019: Not a major drought year.
# 2019: Rainfall in February.
# 2019: Set March as NAs
display_daily(zim_five_stations %>% filter(station == "Plumtree" & year == 2019),
              Stations = "Plumtree", Years = 2019, Variables = "rain")

zim_five_stations <- zim_five_stations %>%
  dplyr::mutate(rain = ifelse(station == "Plumtree" & year == 2019 & month == 3, NA, rain))

saveRDS(zim_five_stations, here("data", "zim_five_stations_1979_qc.RDS")) 
