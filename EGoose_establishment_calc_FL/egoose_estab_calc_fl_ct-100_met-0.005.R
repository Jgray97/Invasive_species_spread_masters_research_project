# CALCULATING ESTABLISHMENT SCORE VALUES FOR THE EGYPTIAN GOOSE ACROSS FLORIDA
# USING CHECKLIST THRESHOLD = 100 AND MINIMUM ENCOUNTER RATE THRESHOLD = 0.005
# Author = John Gray
# Email = greyjohn15@gmail.com
# Last edit = 29/08/2023


# loading necessary packages ----
library(rgbif)
library(tidyr)
library(ggplot2)
library(dggridR)
library(dplyr)
library(sp)
library(gganimate)
library(transformr)
library(terra)
library(tidyterra)

### Importing dataset and defining necessary variables from data prep script ----

# set ebd path (probably no need to overwrite if it says there's already a path)
auk::auk_set_ebd_path("D:/John_Gray_research_project/ebd_NZ_relMar-2023")

# import target species and country dataset
# change name of object accordingly

ebird_FL_egoose <- read.csv("data/ebd_FL_egoose_zf.csv")
cell_stats <- read.csv("data/ebd_FL_egoose_cell_stats.csv")

dggs <- dgconstruct(res = 8, projection = "ISEA", metric = TRUE, 
                    resround = 'nearest')

### Creating function to find year of establishment in a given cell ----

# Function takes the target hexagonal cell, the threshold for minimum number of
# checklists to be considered, and the maximum drop in encounter rate (which 
# defines the criterion for establishment) as inputs and returns the year of 
# establishment in that cell

establishment_year_calc <- function(x, checklist_threshold, min_encounter_threshold) {
  
  ## subset myna dataset to target hexagon
  subject_cell <- subset(ebird_FL_egoose, ebird_FL_egoose$cell == x)
  
  ## get yearly detection stats for hexagon
  subject_cell_stats <- subject_cell %>%
    group_by(year) %>%
    summarise(n_checklists = n(), mean_detection = mean(species_observed))
  
  ## find quartiles of the number of checklists for the years after which 
  ## n_checklists is always greater than threshold value 
  
  # first ignore 2023 as these data are incomplete
  
  subject_cell_stats <- subset(subject_cell_stats, 
                               subject_cell_stats$year != 2023)
  
  # then get rid of rows if there's any year after in which n_checklists is
  # lower than threshold value
  
  subject_cell_stats$min_checklist <- numeric(length = nrow(subject_cell_stats))
  
  for (i in 1:nrow(subject_cell_stats)) {
    subject_cell_stats$min_checklist[i] <- min(subject_cell_stats$n_checklists[i:nrow(subject_cell_stats)])
  }
  
  big_cell_stats <- subset(subject_cell_stats, 
                           subject_cell_stats$min_checklist >= checklist_threshold)

  
  # create a get out option if n_checklists < checklist threshold for all years
  
  if (nrow(big_cell_stats) > 0) {
    
    # Find quartiles for n_checklists
    quartile_A = ceiling(summary(big_cell_stats$n_checklists)[2])
    quartile_B = ceiling(summary(big_cell_stats$n_checklists)[3])
    quartile_C = ceiling(summary(big_cell_stats$n_checklists)[5])
    
    ## generate a dataframe of encounter rates for each year which contains the
    ## mean number of species from a sample (with replacement) equal in size
    ## to the quartiles - this process repeated 100 times for each year 
    # increase from 100 to get more accuracy but this can drastically increase
    # power reqirements
    subject_cell_encounters <- subject_cell %>%
      group_by(year) %>%
      sample_n(size = quartile_A, replace = TRUE) %>%
      summarise(encounter_rate = mean(species_observed))
    
    encounters <- data.frame(matrix(NA, 
                                    nrow = length(subject_cell_encounters$encounter_rate), 
                                    ncol=301))
    
    for (i in 2:101) {
      encounters[i] <- (subject_cell %>%
                          group_by(year) %>%
                          sample_n(size = quartile_A, replace = TRUE) %>%
                          summarise(encounter_rate = mean(species_observed)))$encounter_rate
    }
    
    for (i in 102:201) {
      encounters[i] <- (subject_cell %>%
                          group_by(year) %>%
                          sample_n(size = quartile_B, replace = TRUE) %>%
                          summarise(encounter_rate = mean(species_observed)))$encounter_rate
    }
    
    for (i in 202:301) {
      encounters[i] <- (subject_cell %>%
                          group_by(year) %>%
                          sample_n(size = quartile_C, replace = TRUE) %>%
                          summarise(encounter_rate = mean(species_observed)))$encounter_rate
    }
    
    # create new dataframe with 3 columns for encounter rate according to 10th, 
    # 50th, and 90th percentile of encounter rates
    
    encounters[1] <- subject_cell_encounters$year
    
    colnames(encounters)[1] <- "year"
    
    encounters_switched_quartile_A <- as.data.frame(t(encounters[2:101]))
    
    quantile_1 <- apply(encounters_switched_quartile_A, 2,
                        function(x) quantile(x, probs= c(0.1, 0.5, 0.9)))
    
    quantile_1 <- as.data.frame(t(quantile_1))
    
    encounters_switched_quartile_B <- as.data.frame(t(encounters[102:201]))
    
    quantile_2 <- apply(encounters_switched_quartile_B, 2,
                        function(x) quantile(x, probs= c(0.1, 0.5, 0.9)))
    
    quantile_2 <- as.data.frame(t(quantile_2))
    
    encounters_switched_quartile_C <- as.data.frame(t(encounters[202:301]))
    
    quantile_3 <- apply(encounters_switched_quartile_C, 2,
                        function(x) quantile(x, probs= c(0.1, 0.5, 0.9)))
    
    quantile_3 <- as.data.frame(t(quantile_3))
    
    encounter_prob <- data.frame(encounters$year, quantile_1, quantile_2, quantile_3)
    
    encounter_prob <- subset(encounter_prob,
                             encounter_prob$encounters.year != 2023)
    
    encounter_prob$n_checklists <- subject_cell_stats$n_checklists
    encounter_prob$min_checklists <- subject_cell_stats$min_checklist
    
    colnames(encounter_prob) <- c("year", "s_sample_percentile_10", 
                                  "s_sample_median", "s_sample_percentile_90", 
                                  "m_sample_percentile_10", "m_sample_median", 
                                  "m_sample_percentile_90", "l_sample_percentile_10",
                                  "l_sample_median", "l_sample_percentile_90",
                                  "n_checklists", "min_checklists")
    
    ## CLEAN DATA
    
    # Only keep years where every henceforth year has a number of checklists above threshold
    
    encounter_prob <- subset(encounter_prob, 
                             encounter_prob$min_checklist >= checklist_threshold)
    
    ## DEFINING NON-BINARY, CONVOLUTION-ESQUE FUNCTION FOR ESTABLISHMENT
    
    establishment <- numeric(length = nrow(encounter_prob))
    
    for (i in 1: nrow(encounter_prob)) {
      if (i > (nrow(encounter_prob)-3)) {
        establishment[i] <- 0
      } else {
        if (encounter_prob$m_sample_median[i] <= min_encounter_threshold |
            encounter_prob$m_sample_median[i+1] <= min_encounter_threshold |
            encounter_prob$m_sample_median[i+2] <= min_encounter_threshold |
            encounter_prob$m_sample_median[i+3] <= min_encounter_threshold) {
          establishment[i] <- 0
        } else {
            establishment[i] <- (((3/4) * encounter_prob$m_sample_median[i]) +
              ((3/16) * encounter_prob$m_sample_median[i+1]) +
                  ((3/64) * encounter_prob$m_sample_median[i+2]) +
                  ((3/256) * encounter_prob$m_sample_median[i+3]))
          }
        } 
      }
    
    establishment_with_years <- data.frame(encounter_prob$year, establishment)
    
    colnames(establishment_with_years) <- c("year", "establishment_score")
    
    modern_years <- data.frame( years = c(2000:2019), 
                                dummy_establishment = numeric(length = 20))
    
    for (i in 2000:2019) {
      if (i %in% establishment_with_years$year) {
        modern_years$dummy_establishment[i-1999] <- establishment_with_years$establishment_score[which(establishment_with_years$year == i)] 
      } else {
        modern_years$dummy_establishment[i-1999] <- 0
      }
    }
    
    final_establishment <- modern_years$dummy_establishment
    
    return(final_establishment)
  } else {
    return(numeric(20))
  }
}

### Apply function to every hexagonal cell in targetted country ebird data ----

FL_cells <- unique(ebird_FL_egoose$cell)

FL_establishments <- sapply(FL_cells, establishment_year_calc,
                            checklist_threshold = 100, 
                            min_encounter_threshold = 0.005)

colnames(FL_establishments) <- FL_cells

year <- c(2000:2019)

FL_establishments <- cbind(year, FL_establishments)

# Create additional similar function that finds number of checklists for each
# year

checklist_finder <- function(x, checklist_threshold){
  ## subset egoose dataset to target hexagon
  subject_cell <- subset(ebird_FL_egoose, ebird_FL_egoose$cell == x)
  
  ## get yearly detection stats for hexagon
  subject_cell_stats <- subject_cell %>%
    group_by(year) %>%
    summarise(n_checklists = n())
  
  subject_cell_stats <- subset(subject_cell_stats, 
                               subject_cell_stats$year != 2023)
  
  subject_cell_stats$min_checklist <- numeric(length = nrow(subject_cell_stats))
  
  for (i in 1:nrow(subject_cell_stats)) {
    subject_cell_stats$min_checklist[i] <- min(subject_cell_stats$n_checklists[i:nrow(subject_cell_stats)])
  }
  
  subject_cell_stats <- subset(subject_cell_stats, 
                           subject_cell_stats$min_checklist >= checklist_threshold)
  
  ## keep only rows in 2000-2019 year bracket
  subject_cell_stats <- subset(subject_cell_stats, 
                               subject_cell_stats$year > 1999 & 
                                 subject_cell_stats$year < 2020)
  
  yearly_min_checklists <- numeric(20)
  
  for (i in 2000:2019) {
    if (i %in% subject_cell_stats$year) {
      yearly_min_checklists[i-1999] <- subject_cell_stats$min_checklist[which(subject_cell_stats$year == i)] 
    } else {
      yearly_min_checklists[i-1999] <- 0
    }
  }
  return(yearly_min_checklists)
}

# Apply checklist function to eBird data for Egyptian Goose
FL_min_checklist_numbers <- sapply(FL_cells, checklist_finder, checklist_threshold = 100)

colnames(FL_min_checklist_numbers) <- FL_cells

FL_min_checklist_numbers <- cbind(year, FL_min_checklist_numbers)

# create additional similar function that finds encounter rate for each year

encounter_rate_calc <- function(x, checklist_threshold) {
  
  ## subset egoose dataset to target hexagon
  subject_cell <- subset(ebird_FL_egoose, ebird_FL_egoose$cell == x)
  
  ## get yearly detection stats for hexagon
  subject_cell_stats <- subject_cell %>%
    group_by(year) %>%
    summarise(n_checklists = n(), mean_detection = mean(species_observed))
  
  ## find quartiles of the number of checklists for the years after which 
  ## n_checklists is always greater than threshold value 
  
  # first ignore 2023 as these data are incomplete
  
  subject_cell_stats <- subset(subject_cell_stats, 
                               subject_cell_stats$year != 2023)
  
  # then get rid of rows if there's any year after in which n_checklists is
  # lower than threshold value
  
  subject_cell_stats$min_checklist <- numeric(length = nrow(subject_cell_stats))
  
  for (i in 1:nrow(subject_cell_stats)) {
    subject_cell_stats$min_checklist[i] <- min(subject_cell_stats$n_checklists[i:nrow(subject_cell_stats)])
  }
  
  big_cell_stats <- subset(subject_cell_stats, 
                           subject_cell_stats$min_checklist >= checklist_threshold)
  
  
  # create a get out option if n_checklists < checklist threshold for all years
  
  if (nrow(big_cell_stats) > 0) {
    
    # Find quartiles for n_checklists
    quartile_A = ceiling(summary(big_cell_stats$n_checklists)[2])
    quartile_B = ceiling(summary(big_cell_stats$n_checklists)[3])
    quartile_C = ceiling(summary(big_cell_stats$n_checklists)[5])
    
    ## generate a dataframe of encounter rates for each year which contains the
    ## mean number of species from a sample (with replacement) equal in size
    ## to the quartiles - this process repeated 100 times for each year 
    # increase from 100 to get more accuracy but this can drastically increase
    # power reqirements
    subject_cell_encounters <- subject_cell %>%
      group_by(year) %>%
      sample_n(size = quartile_A, replace = TRUE) %>%
      summarise(encounter_rate = mean(species_observed))
    
    encounters <- data.frame(matrix(NA, 
                                    nrow = length(subject_cell_encounters$encounter_rate), 
                                    ncol=301))
    
    for (i in 2:101) {
      encounters[i] <- (subject_cell %>%
                          group_by(year) %>%
                          sample_n(size = quartile_A, replace = TRUE) %>%
                          summarise(encounter_rate = mean(species_observed)))$encounter_rate
    }
    
    for (i in 102:201) {
      encounters[i] <- (subject_cell %>%
                          group_by(year) %>%
                          sample_n(size = quartile_B, replace = TRUE) %>%
                          summarise(encounter_rate = mean(species_observed)))$encounter_rate
    }
    
    for (i in 202:301) {
      encounters[i] <- (subject_cell %>%
                          group_by(year) %>%
                          sample_n(size = quartile_C, replace = TRUE) %>%
                          summarise(encounter_rate = mean(species_observed)))$encounter_rate
    }
    
    # create new dataframe with 3 columns for encounter rate according to 10th, 
    # 50th, and 90th percentile of encounter rates
    
    encounters[1] <- subject_cell_encounters$year
    
    colnames(encounters)[1] <- "year"
    
    encounters_switched_quartile_A <- as.data.frame(t(encounters[2:101]))
    
    quantile_1 <- apply(encounters_switched_quartile_A, 2,
                        function(x) quantile(x, probs= c(0.1, 0.5, 0.9)))
    
    quantile_1 <- as.data.frame(t(quantile_1))
    
    encounters_switched_quartile_B <- as.data.frame(t(encounters[102:201]))
    
    quantile_2 <- apply(encounters_switched_quartile_B, 2,
                        function(x) quantile(x, probs= c(0.1, 0.5, 0.9)))
    
    quantile_2 <- as.data.frame(t(quantile_2))
    
    encounters_switched_quartile_C <- as.data.frame(t(encounters[202:301]))
    
    quantile_3 <- apply(encounters_switched_quartile_C, 2,
                        function(x) quantile(x, probs= c(0.1, 0.5, 0.9)))
    
    quantile_3 <- as.data.frame(t(quantile_3))
    
    encounter_prob <- data.frame(encounters$year, quantile_1, quantile_2, quantile_3)
    
    encounter_prob <- subset(encounter_prob,
                             encounter_prob$encounters.year != 2023)
    
    encounter_prob$n_checklists <- subject_cell_stats$n_checklists
    encounter_prob$min_checklists <- subject_cell_stats$min_checklist
    
    colnames(encounter_prob) <- c("year", "s_sample_percentile_10", 
                                  "s_sample_median", "s_sample_percentile_90", 
                                  "m_sample_percentile_10", "m_sample_median", 
                                  "m_sample_percentile_90", "l_sample_percentile_10",
                                  "l_sample_median", "l_sample_percentile_90",
                                  "n_checklists", "min_checklists")
    
    ## CLEAN DATA
    
    # Only keep years where every henceforth year has a number of checklists above threshold
    
    encounter_prob <- subset(encounter_prob, 
                             encounter_prob$min_checklist >= checklist_threshold)
    
    encounter_rate <- data.frame(encounter_prob$year, encounter_prob$m_sample_median)
    
    colnames(encounter_rate) <- c("year", "encounter_rate")
    
    modern_years <- data.frame( years = c(2000:2019), 
                                dummy_rate = numeric(length = 20))
    
    for (i in 2000:2019) {
      if (i %in% encounter_rate$year) {
        modern_years$dummy_rate[i-1999] <- encounter_rate$encounter_rate[which(encounter_rate$year == i)] 
      } else {
        modern_years$dummy_rate[i-1999] <- NA
      }
    }
    
    final_encounter_rate <- modern_years$dummy_rate
    
    return(final_encounter_rate)
  } else {
    return(rep(NA,20))
  }
}

## Apply encounter rate function to eBird data

FL_encounter_rates <- sapply(FL_cells, encounter_rate_calc, 
                             checklist_threshold = 100)

colnames(FL_encounter_rates) <- FL_cells

FL_encounter_rates <- cbind(year, FL_encounter_rates)

# pivot longer the dataframes for establishments, checklists and encounter rates
# for better plotting

FL_establishments <- as.data.frame(FL_establishments)

FL_establishments <- pivot_longer(FL_establishments, !year, 
                                  names_to = "grid_cell",
                                  values_to = "establishment_score")

FL_min_checklist_numbers <- as.data.frame(FL_min_checklist_numbers)

FL_min_checklist_numbers <- pivot_longer(FL_min_checklist_numbers, !year,
                                     names_to = "grid_cell",
                                     values_to = "min_checklists")

FL_encounter_rates <- as.data.frame(FL_encounter_rates)

FL_encounter_rates <- pivot_longer(FL_encounter_rates, !year, 
                                  names_to = "grid_cell",
                                  values_to = "encounter_rate")

# Group together the three dataframes for a final "establishment_tracker" df

establishment_tracker <- cbind(FL_establishments, 
                               FL_min_checklist_numbers$min_checklists,
                               FL_encounter_rates$encounter_rate)

# Rename columns + make grid cell numeric

colnames(establishment_tracker)[4] <- "min_checklists"

colnames(establishment_tracker)[5] <- "encounter_rate"

establishment_tracker$grid_cell <- as.numeric(establishment_tracker$grid_cell)

# adding column to display if enough checklists are there for a given cell
# MAKE SURE TO CHANGE IF STATEMENT CRITERION IF CHECKLIST CRITERION IS CHANGED

checklist_threshold <- 100

establishment_tracker$enough_checklists <- numeric(nrow(establishment_tracker))

for (i in 1:nrow(establishment_tracker)) {
  if (establishment_tracker$min_checklists[i] >= checklist_threshold) {
    establishment_tracker$enough_checklists[i] <- 1
  } else {
    establishment_tracker$enough_checklists[i] <- 0
  }
}

# save dataframe as csv

write.csv(establishment_tracker, "data/FINAL-VERSION-FL_egoose_establishment_ct-100_met-0.005.csv", 
          row.names = FALSE)

### Plotting invasive species establishment ----

# Load in the establishment tracker df - this can be commented out if continuing
# from previous point
establishment_tracker <- read.csv("data/FINAL-VERSION-FL_egoose_establishment_ct-100_met-0.005.csv")

# Create dggridR grid - this isn't necessary if running the whole script all in
# one go

dggs <- dgconstruct(res = 8, projection = "ISEA", metric = TRUE, resround = 'nearest')

year <- c(2000:2019)

# Apply dggridR grid to establishment tracker df

cellcenters <- dgSEQNUM_to_GEO(dggs,establishment_tracker$grid_cell)

cellcenters_df <- data.frame(cellcenters)

establishment_tracker <- cbind(establishment_tracker, cellcenters_df)

FL_egoose_grid <- dgcellstogrid(dggs, establishment_tracker$grid_cell)

FL_egoose_grid <- sp::merge(FL_egoose_grid, establishment_tracker, by.x="seqnum", by.y="grid_cell")

# Load in a map of the world for plotting purposes

countries <- map_data("world")

# alter df so that any points close to international dateline don't look weird
# when plotting

wrapped_grid = st_wrap_dateline(FL_egoose_grid, options = c("WRAPDATELINE=YES","DATELINEOFFSET=180"), quiet = TRUE)

# subset grid dataframe by year to enabling plotting by year

for (i in year) {
  subsetted_egeese <- FL_egoose_grid[FL_egoose_grid$year == i,]
  wrapped_subsetted_egeese <- st_wrap_dateline(subsetted_egeese, options = c("WRAPDATELINE=YES","DATELINEOFFSET=180"), quiet = TRUE)
  
  #store the subsetted dataframes as objects with years in names
  assign(paste0("FL_egoose_grid_", i), subsetted_egeese)
  assign(paste0("wrapped_grid_", i), wrapped_subsetted_egeese)

}

# Manually creating a plot of establishment for each year

#### FOR CREATION OF GIF OF ESTABLISHMENT THROUGH THE YEARS, SEE 'spread animation creator' 
# python file

# Year 2000 plot - ignore for analysis as we consider from 2002

ggplot() +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2000, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2000, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-88,-79), ylim = c(23,32)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in FL_ct-100_met-0.005",
        subtitle = "2000")

ggsave(filename = "plots_fl_egoose_final/ct-100_met-0.005/2000.png", width = 10, height= 10)

# Year 2001 plot - ignore for analysis as we consider from 2002

ggplot() +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2001, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2001, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-88,-79), ylim = c(23,32)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in FL_ct-100_met-0.005",
        subtitle = "2001")

ggsave(filename = "plots_fl_egoose_final/ct-100_met-0.005/2001.png", width = 10, height= 10)

# Year 2002 plot

ggplot() +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2002, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2002, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-88,-79), ylim = c(23,32)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in FL_ct-100_met-0.005",
        subtitle = "2002")

ggsave(filename = "plots_fl_egoose_final/ct-100_met-0.005/2002.png", width = 10, height= 10)

# Year 2003 plot

ggplot() +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2003, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2003, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-88,-79), ylim = c(23,32)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in FL_ct-100_met-0.005",
        subtitle = "2003")

ggsave(filename = "plots_fl_egoose_final/ct-100_met-0.005/2003.png", width = 10, height= 10)

# Year 2004 plot

ggplot() +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2004, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2004, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-88,-79), ylim = c(23,32)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in FL_ct-100_met-0.005",
        subtitle = "2004")

ggsave(filename = "plots_fl_egoose_final/ct-100_met-0.005/2004.png", width = 10, height= 10)

# Year 2005 plot

ggplot() +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2005, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2005, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-88,-79), ylim = c(23,32)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in FL_ct-100_met-0.005",
        subtitle = "2005")

ggsave(filename = "plots_fl_egoose_final/ct-100_met-0.005/2005.png", width = 10, height= 10)

# Year 2006 plot

ggplot() +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2006, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2006, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-88,-79), ylim = c(23,32)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in FL_ct-100_met-0.005",
        subtitle = "2006")

ggsave(filename = "plots_fl_egoose_final/ct-100_met-0.005/2006.png", width = 10, height= 10)

# Year 2007 plot

ggplot() +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2007, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2007, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-88,-79), ylim = c(23,32)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in FL_ct-100_met-0.005",
        subtitle = "2007")

ggsave(filename = "plots_fl_egoose_final/ct-100_met-0.005/2007.png", width = 10, height= 10)

# Year 2008 plot

ggplot() +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2008, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2008, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-88,-79), ylim = c(23,32)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in FL_ct-100_met-0.005",
        subtitle = "2008")

ggsave(filename = "plots_fl_egoose_final/ct-100_met-0.005/2008.png", width = 10, height= 10)

# Year 2009 plot

ggplot() +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2009, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2009, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-88,-79), ylim = c(23,32)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in FL_ct-100_met-0.005",
        subtitle = "2009")

ggsave(filename = "plots_fl_egoose_final/ct-100_met-0.005/2009.png", width = 10, height= 10)

# Year 2010 plot

ggplot() +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2010, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2010, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-88,-79), ylim = c(23,32)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in FL_ct-100_met-0.005",
        subtitle = "2010")

ggsave(filename = "plots_fl_egoose_final/ct-100_met-0.005/2010.png", width = 10, height= 10)

# Year 2011 plot

ggplot() +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2011, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2011, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-88,-79), ylim = c(23,32)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in FL_ct-100_met-0.005",
        subtitle = "2011")

ggsave(filename = "plots_fl_egoose_final/ct-100_met-0.005/2011.png", width = 10, height= 10)

# Year 2012 plot

ggplot() +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2012, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2012, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-88,-79), ylim = c(23,32)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in FL_ct-100_met-0.005",
        subtitle = "2012")

ggsave(filename = "plots_fl_egoose_final/ct-100_met-0.005/2012.png", width = 10, height= 10)

# Year 2013 plot

ggplot() +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2013, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2013, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-88,-79), ylim = c(23,32)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in FL_ct-100_met-0.005",
        subtitle = "2013")

ggsave(filename = "plots_fl_egoose_final/ct-100_met-0.005/2013.png", width = 10, height= 10)

# Year 2014 plot

ggplot() +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2014, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2014, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-88,-79), ylim = c(23,32)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in FL_ct-100_met-0.005",
        subtitle = "2014")

ggsave(filename = "plots_fl_egoose_final/ct-100_met-0.005/2014.png", width = 10, height= 10)

# Year 2015 plot

ggplot() +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2015, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2015, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-88,-79), ylim = c(23,32)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in FL_ct-100_met-0.005",
        subtitle = "2015")

ggsave(filename = "plots_fl_egoose_final/ct-100_met-0.005/2015.png", width = 10, height= 10)

# Year 2016 plot

ggplot() +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2016, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2016, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) + 
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-88,-79), ylim = c(23,32)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in FL_ct-100_met-0.005",
        subtitle = "2016")

ggsave(filename = "plots_fl_egoose_final/ct-100_met-0.005/2016.png", width = 10, height= 10)

# Year 2017 plot

ggplot() +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2017, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2017, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-88,-79), ylim = c(23,32)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in FL_ct-100_met-0.005",
        subtitle = "2017")

ggsave(filename = "plots_fl_egoose_final/ct-100_met-0.005/2017.png", width = 10, height= 10)

# Year 2018 plot

ggplot() +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2018, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2018, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-88,-79), ylim = c(23,32)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in FL_ct-100_met-0.005",
        subtitle = "2018")

ggsave(filename = "plots_fl_egoose_final/ct-100_met-0.005/2018.png", width = 10, height= 10)

# Year 2019 plot

ggplot() +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2019, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2019, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-88,-79), ylim = c(23,32)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in FL_ct-100_met-0.005",
        subtitle = "2019")

ggsave(filename = "plots_fl_egoose_final/ct-100_met-0.005/2019.png", width = 10, height= 10)


