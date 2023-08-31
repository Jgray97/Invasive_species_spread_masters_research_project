# IMPORTING AND PREPARING EBIRD DATA FOR CALCULATION OF NATIVE SPECIES
# ENCOUNTER RATE BETWEEN 2009 and 2010 IN TEXAS
# Author = John Gray
# Email = greyjohn15@gmail.com
# Last Edit = 31/08/2023

### Loading necessary packages ----

library(auk)
library(sf)
library(rnaturalearth)
library(dplyr)
library(lubridate)
library(gridExtra)
library(tidyverse)
library(dggridR)

# resolve namespace conflicts
select <- dplyr::select

### Setting up data directory ----

# set ebd path (probably no need to overwrite if it says there's already a path)
auk::auk_set_ebd_path("/Users/john/Desktop")

# set up data directory
dir.create("data", showWarnings = FALSE)

ebd <- auk_ebd("ebd_US-TX_199901_202304_relJun-2023.txt",  # change to desired file
               file_sampling = "ebd_sampling_relJun-2023.txt") # change to desired file


### Filtering imported data for analysis ----

# define features observations
ebd_filters <- ebd %>%
  auk_state("US-TX") %>%
  auk_date(date = c("2009-01-01", "2010-12-31")) %>%
  auk_protocol(protocol = c("Stationary", "Traveling")) %>% # Standard eBird protocol for data analysis
  auk_complete()

ebd_filters

### Creating output files ----

# output files
data_dir <- "data"
if (!dir.exists(data_dir)) {
  dir.create(data_dir)
}
f_ebd <- file.path(data_dir, "ebd_all_species_TX_2009-2010.txt")
f_sampling <- file.path(data_dir, "ebd_checklists_all_species_TX_2009-2010.txt")

# only run if the files don't already exist
if (!file.exists(f_ebd)) {
  auk_filter(ebd_filters, file = f_ebd, file_sampling = f_sampling)
}

ebd_zf <- auk_zerofill(f_ebd, f_sampling, collapse = TRUE)

### Data prep for analysis ----

# Steps here follow good practice for eBird data analysis

# function to convert time observation to hours since midnight
time_to_decimal <- function(x) {
  x <- hms(x, quiet = TRUE)
  hour(x) + minute(x) / 60 + second(x) / 3600
}  

# clean up variables
ebd_zf <- ebd_zf %>% 
  mutate(
    # convert X to NA
    observation_count = if_else(observation_count == "X", 
                                NA_character_, observation_count),
    observation_count = as.integer(observation_count),
    # effort_distance_km to 0 for non-travelling counts
    effort_distance_km = if_else(protocol_type != "Traveling", 
                                 0, effort_distance_km),
    # convert time to decimal hours since midnight
    time_observations_started = time_to_decimal(time_observations_started),
    # split date into year and day of year
    year = year(observation_date),
    day_of_year = yday(observation_date)
  )

# additional filtering
ebd_zf_filtered <- ebd_zf %>% 
  filter(
    # effort filters
    duration_minutes <= 5 * 60,
    effort_distance_km <= 5,
    # 10 or fewer observers
    number_observers <= 10)

# remove redundant variables
ebird_TX_all_species <- ebd_zf_filtered %>% # adjust name according to target species and country
  select(checklist_id, observer_id, sampling_event_identifier,
         scientific_name,
         observation_count, species_observed, 
         state_code, locality_id, latitude, longitude,
         protocol_type, all_species_reported,
         observation_date, year, day_of_year,
         time_observations_started, 
         duration_minutes, effort_distance_km,
         number_observers)

### Adding in hexagonal grids to data ----

dggs <- dgconstruct(res = 8, projection = "ISEA", metric = TRUE, resround = 'nearest')

ebird_TX_all_species$cell <- dgGEO_to_SEQNUM(dggs, ebird_TX_all_species$longitude, ebird_TX_all_species$latitude)$seqnum

cell_stats <- ebird_TX_all_species %>%
  group_by(year, cell) %>%
  summarise(n_checklists = n(), mean_detection = mean(species_observed))

### Producing output csv files ----

write_csv(ebird_TX_all_species, "data/ebd_TX_all_species_2009-2010_zf.csv", na = "")

write_csv(cell_stats, "data/ebd_TX_all_species_2009-2010_cell_stats.csv", na = "")
