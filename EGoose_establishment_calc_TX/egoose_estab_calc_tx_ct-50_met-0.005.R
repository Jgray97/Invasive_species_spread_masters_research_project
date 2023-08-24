# PRODUCING RESULTS ON THE SPREAD OF ESTABLISHMENT OF INVASIVE SPECIES FROM
# CLEANED EBIRD DATA
# Author = John Gray
# Email = greyjohn15@gmail.com
# Last edit = 24/05/2023


# Installing and loading necessary packages ----
install.packages("gganimate")
install.packages("transformr")
install.packages("rgbif")
install.packages("gganimate")
install.packages("transformr")
install.packages("tidyterra")

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

ebird_TX_egoose <- read.csv("data/ebd_TX_egoose_zf.csv")
cell_stats <- read.csv("data/ebd_TX_egoose_cell_stats.csv")

dggs <- dgconstruct(res = 8, projection = "ISEA", metric = TRUE, 
                    resround = 'nearest')

### Creating function to find year of establishment in a given cell ----

# Function takes the target hexagonal cell, the threshold for minimum number of
# checklists to be considered, and the maximum drop in encounter rate (which 
# defines the criterion for establishment) as inputs and returns the year of 
# establishment in that cell

establishment_year_calc <- function(x, checklist_threshold, min_encounter_threshold) {
  
  ## subset myna dataset to target hexagon
  subject_cell <- subset(ebird_TX_egoose, ebird_TX_egoose$cell == x)
  
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

TX_cells <- unique(ebird_TX_egoose$cell)

TX_establishments <- sapply(TX_cells, establishment_year_calc,
                            checklist_threshold = 50, 
                            min_encounter_threshold = 0.005)

colnames(TX_establishments) <- TX_cells

year <- c(2000:2019)

TX_establishments <- cbind(year, TX_establishments)

# Create additional similar dataframe that finds number of checklists for each
# year

checklist_finder <- function(x, checklist_threshold){
  ## subset myna dataset to target hexagon
  subject_cell <- subset(ebird_TX_egoose, ebird_TX_egoose$cell == x)
  
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

TX_min_checklist_numbers <- sapply(TX_cells, checklist_finder, checklist_threshold = 50)

colnames(TX_min_checklist_numbers) <- TX_cells

TX_min_checklist_numbers <- cbind(year, TX_min_checklist_numbers)

# create additional similar dataframe that finds encounter rate for each year

encounter_rate_calc <- function(x, checklist_threshold) {
  
  ## subset myna dataset to target hexagon
  subject_cell <- subset(ebird_TX_egoose, ebird_TX_egoose$cell == x)
  
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

TX_encounter_rates <- sapply(TX_cells, encounter_rate_calc, 
                             checklist_threshold = 50)

colnames(TX_encounter_rates) <- TX_cells

TX_encounter_rates <- cbind(year, TX_encounter_rates)

# pivot longer the dataframes for better plotting

TX_establishments <- as.data.frame(TX_establishments)

TX_establishments <- pivot_longer(TX_establishments, !year, 
                                  names_to = "grid_cell",
                                  values_to = "establishment_score")

TX_min_checklist_numbers <- as.data.frame(TX_min_checklist_numbers)

TX_min_checklist_numbers <- pivot_longer(TX_min_checklist_numbers, !year,
                                     names_to = "grid_cell",
                                     values_to = "min_checklists")

TX_encounter_rates <- as.data.frame(TX_encounter_rates)

TX_encounter_rates <- pivot_longer(TX_encounter_rates, !year, 
                                  names_to = "grid_cell",
                                  values_to = "encounter_rate")

establishment_tracker <- cbind(TX_establishments, 
                               TX_min_checklist_numbers$min_checklists,
                               TX_encounter_rates$encounter_rate)

colnames(establishment_tracker)[4] <- "min_checklists"

colnames(establishment_tracker)[5] <- "encounter_rate"

establishment_tracker$grid_cell <- as.numeric(establishment_tracker$grid_cell)

# adding column to display if enough checklists are there for a given cell
# MAKE SURE TO CHANGE IF STATEMENT CRITERION IF CHECKLIST CRITERION IS CHANGED

checklist_threshold <- 50

establishment_tracker$enough_checklists <- numeric(nrow(establishment_tracker))

for (i in 1:nrow(establishment_tracker)) {
  if (establishment_tracker$min_checklists[i] >= checklist_threshold) {
    establishment_tracker$enough_checklists[i] <- 1
  } else {
    establishment_tracker$enough_checklists[i] <- 0
  }
}

# save dataframe as csv

write.csv(establishment_tracker, "data/TX_egoose_establishment_ct-50_met-0.005.csv", 
          row.names = FALSE)


# read in csv file if previous step is already done
# comment below line out if running for the 1st time
# establishment_tracker <- read.csv("data/FL_egoose_establishment_ct-50_met-0.005.csv")

### Plotting invasive species establishment ----

# dggridR stuff

# dggs <- dgconstruct(res = 8, projection = "ISEA", metric = TRUE, resround = 'nearest') <- Necessary if loaded in from establishment tracker

year <- c(2000:2019)

establishment_tracker <- read.csv("data/TX_egoose_establishment_ct-50_met-0.005.csv")

dggs <- dgconstruct(res = 8, projection = "ISEA", metric = TRUE, 
                    resround = 'nearest')

cellcenters <- dgSEQNUM_to_GEO(dggs,establishment_tracker$grid_cell)

cellcenters_df <- data.frame(cellcenters)

establishment_tracker <- cbind(establishment_tracker, cellcenters_df)

write.csv(establishment_tracker, "data/TX_ct-50_met-0.005_establishment-tracker.csv")
TX_egoose_grid <- dgcellstogrid(dggs, establishment_tracker$grid_cell)

write.csv(TX_egoose_grid, "data/TX_egoose_grid.csv")

TX_egoose_grid <- sp::merge(TX_egoose_grid, establishment_tracker, by.x="seqnum", by.y="grid_cell")

save(TX_egoose_grid, file = "data/TX_egoose_grid.Rdata")

countries <- map_data("world")

USA <- map_data("state")

wrapped_grid = st_wrap_dateline(TX_egoose_grid, options = c("WRAPDATELINE=YES","DATELINEOFFSET=180"), quiet = TRUE)

# subset grid dataframe by year to enabling plotting by year

for (i in year) {
  subsetted_egeese <- TX_egoose_grid[TX_egoose_grid$year == i,]
  wrapped_subsetted_egeese <- st_wrap_dateline(subsetted_egeese, options = c("WRAPDATELINE=YES","DATELINEOFFSET=180"), quiet = TRUE)
  
  #store the subsetted dataframes as objects with years in names
  assign(paste0("TX_egoose_grid_", i), subsetted_egeese)
  assign(paste0("wrapped_grid_", i), wrapped_subsetted_egeese)

}

# Manually creating every graph (I couldn't work out a quicker way to do this)

#### FOR GIF OF ESTABLISHMENT THROUGH THE YEARS, SEE 'greating gif map' 
# python file

ggplot() +
  geom_polygon(data=USA, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2000, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2000, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-107,-92), ylim = c(25,37)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in TX_ct-50_met-0.005",
        subtitle = "2000")

ggsave(filename = "plots_TX_egoose/2000.png", width = 10, height= 10)

ggplot() +
  geom_polygon(data=USA, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2001, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2001, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-107,-92), ylim = c(25,37)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in TX_ct-50_met-0.005",
        subtitle = "2001")

ggsave(filename = "plots_TX_egoose/2001.png", width = 10, height= 10)

ggplot() +
  geom_polygon(data=USA, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2002, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2002, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-107,-92), ylim = c(25,37)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in TX_ct-50_met-0.005",
        subtitle = "2002")

ggsave(filename = "plots_TX_egoose/2002.png", width = 10, height= 10)

ggplot() +
  geom_polygon(data=USA, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2003, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2003, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-107,-92), ylim = c(25,37)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in TX_ct-50_met-0.005",
        subtitle = "2003")

ggsave(filename = "plots_TX_egoose/2003.png", width = 10, height= 10)

ggplot() +
  geom_polygon(data=USA, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2004, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2004, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-107,-92), ylim = c(25,37)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in TX_ct-50_met-0.005",
        subtitle = "2004")

ggsave(filename = "plots_TX_egoose/2004.png", width = 10, height= 10)

ggplot() +
  geom_polygon(data=USA, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2005, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2005, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-107,-92), ylim = c(25,37)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in TX_ct-50_met-0.005",
        subtitle = "2005")

ggsave(filename = "plots_TX_egoose/2005.png", width = 10, height= 10)

ggplot() +
  geom_polygon(data=USA, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2006, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2006, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-107,-92), ylim = c(25,37)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in TX_ct-50_met-0.005",
        subtitle = "2006")

ggsave(filename = "plots_TX_egoose/2006.png", width = 10, height= 10)

ggplot() +
  geom_polygon(data=USA, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2007, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2007, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-107,-92), ylim = c(25,37)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in TX_ct-50_met-0.005",
        subtitle = "2007")

ggsave(filename = "plots_TX_egoose/2007.png", width = 10, height= 10)

ggplot() +
  geom_polygon(data=USA, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2008, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2008, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-107,-92), ylim = c(25,37)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in TX_ct-50_met-0.005",
        subtitle = "2008")

ggsave(filename = "plots_TX_egoose/2008.png", width = 10, height= 10)

ggplot() +
  geom_polygon(data=USA, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2009, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2009, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-107,-92), ylim = c(25,37)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in TX_ct-50_met-0.005",
        subtitle = "2009")

ggsave(filename = "plots_TX_egoose/2009.png", width = 10, height= 10)

ggplot() +
  geom_polygon(data=USA, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2010, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2010, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-107,-92), ylim = c(25,37)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in TX_ct-50_met-0.005",
        subtitle = "2010")

ggsave(filename = "plots_TX_egoose/2010.png", width = 10, height= 10)

ggplot() +
  geom_polygon(data=USA, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2011, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2011, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-107,-92), ylim = c(25,37)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in TX_ct-50_met-0.005",
        subtitle = "2011")

ggsave(filename = "plots_TX_egoose/2011.png", width = 10, height= 10)

ggplot() +
  geom_polygon(data=USA, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2012, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2012, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-107,-92), ylim = c(25,37)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in TX_ct-50_met-0.005",
        subtitle = "2012")

ggsave(filename = "plots_TX_egoose/2012.png", width = 10, height= 10)

ggplot() +
  geom_polygon(data=USA, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2013, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2013, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-107,-92), ylim = c(25,37)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in TX_ct-50_met-0.005",
        subtitle = "2013")

ggsave(filename = "plots_TX_egoose/2013.png", width = 10, height= 10)

ggplot() +
  geom_polygon(data=USA, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2014, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2014, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-107,-92), ylim = c(25,37)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in TX_ct-50_met-0.005",
        subtitle = "2014")

ggsave(filename = "plots_TX_egoose/2014.png", width = 10, height= 10)

ggplot() +
  geom_polygon(data=USA, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2015, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2015, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-107,-92), ylim = c(25,37)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in TX_ct-50_met-0.005",
        subtitle = "2015")

ggsave(filename = "plots_TX_egoose/2015.png", width = 10, height= 10)

ggplot() +
  geom_polygon(data=USA, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2016, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2016, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) + 
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-107,-92), ylim = c(25,37)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in TX_ct-50_met-0.005",
        subtitle = "2016")

ggsave(filename = "plots_TX_egoose/2016.png", width = 10, height= 10)

ggplot() +
  geom_polygon(data=USA, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2017, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2017, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-107,-92), ylim = c(25,37)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in TX_ct-50_met-0.005",
        subtitle = "2017")

ggsave(filename = "plots_TX_egoose/2017.png", width = 10, height= 10)

ggplot() +
  geom_polygon(data=USA, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2018, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2018, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-107,-92), ylim = c(25,37)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in TX_ct-50_met-0.005",
        subtitle = "2018")

ggsave(filename = "plots_TX_egoose/2018.png", width = 10, height= 10)

ggplot() +
  geom_polygon(data=USA, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_sf(data = wrapped_grid_2019, aes(fill = enough_checklists), color=alpha("white", 0.4)) +
  geom_point(data = wrapped_grid_2019, aes( x = lon_deg, y = lat_deg, 
                                            size = establishment_score)) +
  scale_size_continuous(limits = c(0,0.3), 
                        breaks = c(0.1,0.2,0.3),
                        range = c(-0.5,3)) +
  scale_fill_gradient(low=NA, high = alpha("Red",0.5)) +
  coord_sf(xlim=c(-107,-92), ylim = c(25,37)) + 
  labs( x = "Longitude", y = "Latitude", 
        title = "Establishment spread of the Egyptian Goose in TX_ct-50_met-0.005",
        subtitle = "2019")

ggsave(filename = "plots_TX_egoose/2019.png", width = 10, height= 10)
