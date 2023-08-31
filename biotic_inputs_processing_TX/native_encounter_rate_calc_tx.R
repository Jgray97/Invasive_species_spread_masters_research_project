# PRODUCING A TABLE FOR OVERALL ENCOUNTER RATE VALUES ACROSS GRID CELLS IN 
# TEXAS BETWEEN 2002 AND 2015
# Author = John Gray
# Email = greyjohn15@gmail.com
# Last edit = 31/08/2023

# Load necessary packages ----

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

# Calculation of native species encounter rates for 2002-2005 ----

# Load in zerofill datafeame where each entry corresponds to a particular
# checklist and whether one of all observed species in Texas between 2002 and 
# 2005 was observed during that checklist
ebird_TX_all_species <- read.csv("data/ebd_TX_all_species_2002-2005_zf.csv")

# Filter zf dataframe to just entries where the species was observed
just_positives <- filter(ebird_TX_all_species, species_observed == TRUE)

# Create dataframe which lists all checklists and the grid cell they belong to
checklist_list <- ebird_TX_all_species %>%
  group_by(checklist_id) %>%
  summarise(cell = first(cell))

# Create dataframe which lists number of checklists per grid cell
grid_cell_checklists <- checklist_list %>%
  group_by(cell) %>%
  summarise(n_checklists = n())

# calculate number of total checklists between 2002 and 2005
n_checklists <- length(unique(ebird_TX_all_species$checklist_id))

# create a list of all species observed between 2002 and 2005
species <- unique(ebird_TX_all_species$scientific_name)

# create list of all cells in Texas
texas_cells <- unique(ebird_TX_all_species$cell)

# create placeholder df to be populated with encounter rate values
encounter_rates <- data.frame("species" = character(length = (length(species) * length(texas_cells))), 
                              "cell" = numeric(length = (length(species) * length(texas_cells))),
                              "encounter_rate" = numeric(length = (length(species) * length(texas_cells))))

# populate dataframe species column with list of species for each cell
for (i in 1:length(species)) {
  encounter_rates$species[(length(texas_cells)*(i-1)+1):(length(texas_cells)*i)] <- species[i]
}

# populate dataframe cell column with a cell for each species
for (i in 1:length(texas_cells)) {
  encounter_rates$cell[c(length(texas_cells)*(0:(length(species)-1))+i)] <- texas_cells[i]
}

# calculate encounter rate for each species in each cell
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$encounter_rate[i] <- (nrow(filter(just_positives, 
                                                    (cell == encounter_rates$cell[i]) & 
                                                      (scientific_name == encounter_rates$species[i])))) /
    (nrow(filter(checklist_list, cell == encounter_rates$cell[i])))
}

# add a "weight" value for each entry according to the number of checklists
# in each cell
encounter_rates$weight <- 0
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$weight[i] <- grid_cell_checklists$n_checklists[grid_cell_checklists$cell == encounter_rates$cell[i]]
}

# save encounter rate df for 2002-2005 in a csv
write.csv(encounter_rates, "data/TX_encounter_rates_2002-2005.csv")

# remove objects to speed up processing
rm(ebird_TX_all_species)
rm(encounter_rates)
rm(just_positives)

# Calculation of native species encounter rates for 2006-2008 ----

# Load in zerofill datafeame where each entry corresponds to a particular
# checklist and whether one of all observed species in Texas between 2006 and 
# 2008 was observed during that checklist
ebird_TX_all_species <- read.csv("data/ebd_TX_all_species_2006-2008_zf.csv")

# Filter zf dataframe to just entries where the species was observed
just_positives <- filter(ebird_TX_all_species, species_observed == TRUE)

# Create dataframe which lists all checklists and the grid cell they belong to
checklist_list <- ebird_TX_all_species %>%
  group_by(checklist_id) %>%
  summarise(cell = first(cell))

# Create dataframe which lists number of checklists per grid cell
grid_cell_checklists <- checklist_list %>%
  group_by(cell) %>%
  summarise(n_checklists = n())

# calculate number of total checklists between 2006 and 2008
n_checklists <- length(unique(ebird_TX_all_species$checklist_id))

# create a list of all species observed between 2006 and 2008
species <- unique(ebird_TX_all_species$scientific_name)

# create list of all cells in Texas
texas_cells <- unique(ebird_TX_all_species$cell)

# create placeholder df to be populated with encounter rate values
encounter_rates <- data.frame("species" = character(length = (length(species) * length(texas_cells))), 
                              "cell" = numeric(length = (length(species) * length(texas_cells))),
                              "encounter_rate" = numeric(length = (length(species) * length(texas_cells))))

# populate dataframe species column with list of species for each cell
for (i in 1:length(species)) {
  encounter_rates$species[(length(texas_cells)*(i-1)+1):(length(texas_cells)*i)] <- species[i]
}

# populate dataframe cell column with a cell for each species
for (i in 1:length(texas_cells)) {
  encounter_rates$cell[c(length(texas_cells)*(0:(length(species)-1))+i)] <- texas_cells[i]
}

# calculate encounter rate for each species in each cell
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$encounter_rate[i] <- (nrow(filter(just_positives, 
                                                    (cell == encounter_rates$cell[i]) & 
                                                      (scientific_name == encounter_rates$species[i])))) /
    (nrow(filter(checklist_list, cell == encounter_rates$cell[i])))
}

# add a "weight" value for each entry according to the number of checklists
# in each cell
encounter_rates$weight <- 0
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$weight[i] <- grid_cell_checklists$n_checklists[grid_cell_checklists$cell == encounter_rates$cell[i]]
}

# save encounter rate df for 2006-2008 in a csv
write.csv(encounter_rates, "data/TX_encounter_rates_2006-2008.csv")

# remove objects to speed up processing
rm(ebird_TX_all_species)
rm(encounter_rates)
rm(just_positives)

# Calculation of native species encounter rates for 2009-2010 ----

# Load in zerofill datafeame where each entry corresponds to a particular
# checklist and whether one of all observed species in Texas between 2009 and 
# 2010 was observed during that checklist
ebird_TX_all_species <- read.csv("data/ebd_TX_all_species_2009-2010_zf.csv")

# Filter zf dataframe to just entries where the species was observed
just_positives <- filter(ebird_TX_all_species, species_observed == TRUE)

# Create dataframe which lists all checklists and the grid cell they belong to
checklist_list <- ebird_TX_all_species %>%
  group_by(checklist_id) %>%
  summarise(cell = first(cell))

# Create dataframe which lists number of checklists per grid cell
grid_cell_checklists <- checklist_list %>%
  group_by(cell) %>%
  summarise(n_checklists = n())

# calculate number of total checklists between 2009 and 2010
n_checklists <- length(unique(ebird_TX_all_species$checklist_id))

# create a list of all species observed between 2009 and 2010
species <- unique(ebird_TX_all_species$scientific_name)

# create list of all cells in Texas
texas_cells <- unique(ebird_TX_all_species$cell)

# create placeholder df to be populated with encounter rate values
encounter_rates <- data.frame("species" = character(length = (length(species) * length(texas_cells))), 
                              "cell" = numeric(length = (length(species) * length(texas_cells))),
                              "encounter_rate" = numeric(length = (length(species) * length(texas_cells))))

# populate dataframe species column with list of species for each cell
for (i in 1:length(species)) {
  encounter_rates$species[(length(texas_cells)*(i-1)+1):(length(texas_cells)*i)] <- species[i]
}

# populate dataframe cell column with a cell for each species
for (i in 1:length(texas_cells)) {
  encounter_rates$cell[c(length(texas_cells)*(0:(length(species)-1))+i)] <- texas_cells[i]
}

# calculate encounter rate for each species in each cell
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$encounter_rate[i] <- (nrow(filter(just_positives, 
                                                    (cell == encounter_rates$cell[i]) & 
                                                      (scientific_name == encounter_rates$species[i])))) /
    (nrow(filter(checklist_list, cell == encounter_rates$cell[i])))
}

# add a "weight" value for each entry according to the number of checklists
# in each cell
encounter_rates$weight <- 0
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$weight[i] <- grid_cell_checklists$n_checklists[grid_cell_checklists$cell == encounter_rates$cell[i]]
}

# save encounter rate df for 2009-2010 in a csv
write.csv(encounter_rates, "data/TX_encounter_rates_2009-2010.csv")

# remove objects to speed up processing
rm(ebird_TX_all_species)
rm(encounter_rates)
rm(just_positives)

# Calculation of native species encounter rates in 2011 ----

# Load in zerofill datafeame where each entry corresponds to a particular
# checklist and whether one of all observed species in Texas in 2011 was 
# observed during that checklist
ebird_TX_all_species <- read.csv("data/ebd_TX_all_species_2011_zf.csv")

# Filter zf dataframe to just entries where the species was observed
just_positives <- filter(ebird_TX_all_species, species_observed == TRUE)

# Create dataframe which lists all checklists and the grid cell they belong to
checklist_list <- ebird_TX_all_species %>%
  group_by(checklist_id) %>%
  summarise(cell = first(cell))

# Create dataframe which lists number of checklists per grid cell
grid_cell_checklists <- checklist_list %>%
  group_by(cell) %>%
  summarise(n_checklists = n())

# calculate number of total checklists in 2011
n_checklists <- length(unique(ebird_TX_all_species$checklist_id))

# create a list of all species observed in 2011
species <- unique(ebird_TX_all_species$scientific_name)

# create list of all cells in Texas
texas_cells <- unique(ebird_TX_all_species$cell)

# create placeholder df to be populated with encounter rate values
encounter_rates <- data.frame("species" = character(length = (length(species) * length(texas_cells))), 
                              "cell" = numeric(length = (length(species) * length(texas_cells))),
                              "encounter_rate" = numeric(length = (length(species) * length(texas_cells))))

# populate dataframe species column with list of species for each cell
for (i in 1:length(species)) {
  encounter_rates$species[(length(texas_cells)*(i-1)+1):(length(texas_cells)*i)] <- species[i]
}

# populate dataframe cell column with a cell for each species
for (i in 1:length(texas_cells)) {
  encounter_rates$cell[c(length(texas_cells)*(0:(length(species)-1))+i)] <- texas_cells[i]
}

# calculate encounter rate for each species in each cell
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$encounter_rate[i] <- (nrow(filter(just_positives, 
                                                    (cell == encounter_rates$cell[i]) & 
                                                      (scientific_name == encounter_rates$species[i])))) /
    (nrow(filter(checklist_list, cell == encounter_rates$cell[i])))
}

# add a "weight" value for each entry according to the number of checklists
# in each cell
encounter_rates$weight <- 0
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$weight[i] <- grid_cell_checklists$n_checklists[grid_cell_checklists$cell == encounter_rates$cell[i]]
}

# save encounter rate df for 2011 in a csv
write.csv(encounter_rates, "data/TX_encounter_rates_2011.csv")

# remove objects to speed up processing
rm(ebird_TX_all_species)
rm(encounter_rates)
rm(just_positives)

# Calculation of native species encounter rates in 2012 ----

# Load in zerofill datafeame where each entry corresponds to a particular
# checklist and whether one of all observed species in Texas in 2012 was 
# observed during that checklist
ebird_TX_all_species <- read.csv("data/ebd_TX_all_species_2012_zf.csv")

# Filter zf dataframe to just entries where the species was observed
just_positives <- filter(ebird_TX_all_species, species_observed == TRUE)

# Create dataframe which lists all checklists and the grid cell they belong to
checklist_list <- ebird_TX_all_species %>%
  group_by(checklist_id) %>%
  summarise(cell = first(cell))

# Create dataframe which lists number of checklists per grid cell
grid_cell_checklists <- checklist_list %>%
  group_by(cell) %>%
  summarise(n_checklists = n())

# calculate number of total checklists in 2012
n_checklists <- length(unique(ebird_TX_all_species$checklist_id))

# create a list of all species observed in 2012
species <- unique(ebird_TX_all_species$scientific_name)

# create list of all cells in Texas
texas_cells <- unique(ebird_TX_all_species$cell)

# create placeholder df to be populated with encounter rate values
encounter_rates <- data.frame("species" = character(length = (length(species) * length(texas_cells))), 
                              "cell" = numeric(length = (length(species) * length(texas_cells))),
                              "encounter_rate" = numeric(length = (length(species) * length(texas_cells))))

# populate dataframe species column with list of species for each cell
for (i in 1:length(species)) {
  encounter_rates$species[(length(texas_cells)*(i-1)+1):(length(texas_cells)*i)] <- species[i]
}

# populate dataframe cell column with a cell for each species
for (i in 1:length(texas_cells)) {
  encounter_rates$cell[c(length(texas_cells)*(0:(length(species)-1))+i)] <- texas_cells[i]
}

# calculate encounter rate for each species in each cell
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$encounter_rate[i] <- (nrow(filter(just_positives, 
                                                    (cell == encounter_rates$cell[i]) & 
                                                      (scientific_name == encounter_rates$species[i])))) /
    (nrow(filter(checklist_list, cell == encounter_rates$cell[i])))
}

# add a "weight" value for each entry according to the number of checklists
# in each cell
encounter_rates$weight <- 0
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$weight[i] <- grid_cell_checklists$n_checklists[grid_cell_checklists$cell == encounter_rates$cell[i]]
}

# save encounter rate df for 2012 in a csv
write.csv(encounter_rates, "data/TX_encounter_rates_2012.csv")

# remove objects to speed up processing
rm(ebird_TX_all_species)
rm(encounter_rates)
rm(just_positives)

# Calculation of native species encounter rates in 2013 ----

# Load in zerofill datafeame where each entry corresponds to a particular
# checklist and whether one of all observed species in Texas in 2013 was 
# observed during that checklist
ebird_TX_all_species <- read.csv("data/ebd_TX_all_species_2013_zf.csv")

# Filter zf dataframe to just entries where the species was observed
just_positives <- filter(ebird_TX_all_species, species_observed == TRUE)

# Create dataframe which lists all checklists and the grid cell they belong to
checklist_list <- ebird_TX_all_species %>%
  group_by(checklist_id) %>%
  summarise(cell = first(cell))

# Create dataframe which lists number of checklists per grid cell
grid_cell_checklists <- checklist_list %>%
  group_by(cell) %>%
  summarise(n_checklists = n())

# calculate number of total checklists in 2013
n_checklists <- length(unique(ebird_TX_all_species$checklist_id))

# create a list of all species observed in 2013
species <- unique(ebird_TX_all_species$scientific_name)

# create list of all cells in Texas
texas_cells <- unique(ebird_TX_all_species$cell)

# create placeholder df to be populated with encounter rate values
encounter_rates <- data.frame("species" = character(length = (length(species) * length(texas_cells))), 
                              "cell" = numeric(length = (length(species) * length(texas_cells))),
                              "encounter_rate" = numeric(length = (length(species) * length(texas_cells))))

# populate dataframe species column with list of species for each cell
for (i in 1:length(species)) {
  encounter_rates$species[(length(texas_cells)*(i-1)+1):(length(texas_cells)*i)] <- species[i]
}

# populate dataframe cell column with a cell for each species
for (i in 1:length(texas_cells)) {
  encounter_rates$cell[c(length(texas_cells)*(0:(length(species)-1))+i)] <- texas_cells[i]
}

# calculate encounter rate for each species in each cell
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$encounter_rate[i] <- (nrow(filter(just_positives, 
                                                    (cell == encounter_rates$cell[i]) & 
                                                      (scientific_name == encounter_rates$species[i])))) /
    (nrow(filter(checklist_list, cell == encounter_rates$cell[i])))
}

# add a "weight" value for each entry according to the number of checklists
# in each cell
encounter_rates$weight <- 0
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$weight[i] <- grid_cell_checklists$n_checklists[grid_cell_checklists$cell == encounter_rates$cell[i]]
}

# save encounter rate df for 2013 in a csv
write.csv(encounter_rates, "data/TX_encounter_rates_2013.csv")

# remove objects to speed up processing
rm(ebird_TX_all_species)
rm(encounter_rates)
rm(just_positives)

# Calculation of native species encounter rates in 2014 ----

# Load in zerofill datafeame where each entry corresponds to a particular
# checklist and whether one of all observed species in Texas in 2014 was 
# observed during that checklist
ebird_TX_all_species <- read.csv("data/ebd_TX_all_species_2014_zf.csv")

# Filter zf dataframe to just entries where the species was observed
just_positives <- filter(ebird_TX_all_species, species_observed == TRUE)

# Create dataframe which lists all checklists and the grid cell they belong to
checklist_list <- ebird_TX_all_species %>%
  group_by(checklist_id) %>%
  summarise(cell = first(cell))

# Create dataframe which lists number of checklists per grid cell
grid_cell_checklists <- checklist_list %>%
  group_by(cell) %>%
  summarise(n_checklists = n())

# calculate number of total checklists in 2014
n_checklists <- length(unique(ebird_TX_all_species$checklist_id))

# create a list of all species observed in 2014
species <- unique(ebird_TX_all_species$scientific_name)

# create list of all cells in Texas
texas_cells <- unique(ebird_TX_all_species$cell)

# create placeholder df to be populated with encounter rate values
encounter_rates <- data.frame("species" = character(length = (length(species) * length(texas_cells))), 
                              "cell" = numeric(length = (length(species) * length(texas_cells))),
                              "encounter_rate" = numeric(length = (length(species) * length(texas_cells))))

# populate dataframe species column with list of species for each cell
for (i in 1:length(species)) {
  encounter_rates$species[(length(texas_cells)*(i-1)+1):(length(texas_cells)*i)] <- species[i]
}

# populate dataframe cell column with a cell for each species
for (i in 1:length(texas_cells)) {
  encounter_rates$cell[c(length(texas_cells)*(0:(length(species)-1))+i)] <- texas_cells[i]
}

# calculate encounter rate for each species in each cell
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$encounter_rate[i] <- (nrow(filter(just_positives, 
                                                    (cell == encounter_rates$cell[i]) & 
                                                      (scientific_name == encounter_rates$species[i])))) /
    (nrow(filter(checklist_list, cell == encounter_rates$cell[i])))
}

# add a "weight" value for each entry according to the number of checklists
# in each cell
encounter_rates$weight <- 0
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$weight[i] <- grid_cell_checklists$n_checklists[grid_cell_checklists$cell == encounter_rates$cell[i]]
}

# save encounter rate df for 2014 in a csv
write.csv(encounter_rates, "data/TX_encounter_rates_2014.csv")

# remove objects to speed up processing
rm(ebird_TX_all_species)
rm(encounter_rates)
rm(just_positives)

# Calculation of native species encounter rates in 2015 ----

# Load in zerofill datafeame where each entry corresponds to a particular
# checklist and whether one of all observed species in Texas in 2015 was 
# observed during that checklist
ebird_TX_all_species <- read.csv("data/ebd_TX_all_species_2015_zf.csv")

# Filter zf dataframe to just entries where the species was observed
just_positives <- filter(ebird_TX_all_species, species_observed == TRUE)

# Create dataframe which lists all checklists and the grid cell they belong to
checklist_list <- ebird_TX_all_species %>%
  group_by(checklist_id) %>%
  summarise(cell = first(cell))

# Create dataframe which lists number of checklists per grid cell
grid_cell_checklists <- checklist_list %>%
  group_by(cell) %>%
  summarise(n_checklists = n())

# calculate number of total checklists in 2015
n_checklists <- length(unique(ebird_TX_all_species$checklist_id))

# create a list of all species observedin 2015
species <- unique(ebird_TX_all_species$scientific_name)

# create list of all cells in Texas
texas_cells <- unique(ebird_TX_all_species$cell)

# create placeholder df to be populated with encounter rate values
encounter_rates <- data.frame("species" = character(length = (length(species) * length(texas_cells))), 
                              "cell" = numeric(length = (length(species) * length(texas_cells))),
                              "encounter_rate" = numeric(length = (length(species) * length(texas_cells))))

# populate dataframe species column with list of species for each cell
for (i in 1:length(species)) {
  encounter_rates$species[(length(texas_cells)*(i-1)+1):(length(texas_cells)*i)] <- species[i]
}

# populate dataframe cell column with a cell for each species
for (i in 1:length(texas_cells)) {
  encounter_rates$cell[c(length(texas_cells)*(0:(length(species)-1))+i)] <- texas_cells[i]
}

# calculate encounter rate for each species in each cell
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$encounter_rate[i] <- (nrow(filter(just_positives, 
                                                    (cell == encounter_rates$cell[i]) & 
                                                      (scientific_name == encounter_rates$species[i])))) /
    (nrow(filter(checklist_list, cell == encounter_rates$cell[i])))
}

# add a "weight" value for each entry according to the number of checklists
# in each cell
encounter_rates$weight <- 0
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$weight[i] <- grid_cell_checklists$n_checklists[grid_cell_checklists$cell == encounter_rates$cell[i]]
}

# save encounter rate df for 2015 in a csv
write.csv(encounter_rates, "data/TX_encounter_rates_2015.csv")