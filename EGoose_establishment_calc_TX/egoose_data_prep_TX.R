# IMPORTING AND PREPARING EBIRD DATA FOR ANALYSIS OF THE SPREAD OF ESTABLISHMENT
# OF AN INVASIVE SPECIES
# Author = John Gray
# Email = greyjohn15@gmail.com
# Last Edit = 10/05/2023

# This is a neater version of 'Myna_NZ_spread_analysis' - if there are any errors
# refer to original script

### Installing and loading necessary packages ----

# install necessary packages
install.packages("remotes")
remotes::install_github("mstrimas/ebppackages")

# load packages
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

# set up data directory
dir.create("data", showWarnings = FALSE)

ebd <- auk_ebd("ebd_US-TX_199901_202304_relJun-2023.txt",  # change to desired file
               file_sampling = "ebd_sampling_relJun-2023.txt") # change to desired file


### Filtering imported data for analysis ----

# define features observations
ebd_filters <- ebd %>%
  auk_species("Egyptian Goose") %>% # Change this to species you are looking at 
  auk_state("US-TX") %>%
  auk_protocol(protocol = c("Stationary", "Traveling")) %>% # Standard eBird protocol for data analysis
  auk_complete()

ebd_filters

### Creating output files ----

# output files
data_dir <- "data"
if (!dir.exists(data_dir)) {
  dir.create(data_dir)
}
f_ebd <- file.path(data_dir, "ebd_egoose_TX.txt")
f_sampling <- file.path(data_dir, "ebd_checklists_egoose_TX.txt")

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
ebird_TX_egoose <- ebd_zf_filtered %>% # adjust name according to target species and country
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

ebird_TX_egoose$cell <- dgGEO_to_SEQNUM(dggs, ebird_TX_egoose$longitude, ebird_TX_egoose$latitude)$seqnum

cell_stats <- ebird_TX_egoose %>%
  group_by(year, cell) %>%
  summarise(n_checklists = n(), mean_detection = mean(species_observed))

### Producing output csv files ----

write_csv(ebird_TX_egoose, "data/ebd_TX_egoose_zf.csv", na = "")

write_csv(cell_stats, "data/ebd_TX_egoose_cell_stats.csv", na = "")

