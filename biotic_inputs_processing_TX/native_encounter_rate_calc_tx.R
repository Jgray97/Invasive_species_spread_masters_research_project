##### THE 2013 WAY IS THE QUICKEST WAY ADOPT THIS APPROACH GOING FORWARD


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

# 2000-2005 ----

ebird_TX_all_species <- read.csv("data/ebd_TX_all_species_2000-2005_zf.csv")

just_positives <- filter(ebird_TX_all_species, species_observed == TRUE)

checklist_list <- ebird_TX_all_species %>%
  group_by(checklist_id) %>%
  summarise(cell = first(cell))

grid_cell_checklists <- checklist_list %>%
  group_by(cell) %>%
  summarise(n_checklists = n())

n_checklists <- length(unique(ebird_TX_all_species$checklist_id))

species <- unique(ebird_TX_all_species$scientific_name)
texas_cells <- unique(ebird_TX_all_species$cell)

encounter_rates <- data.frame("species" = character(length = (length(species) * length(texas_cells))), 
                              "cell" = numeric(length = (length(species) * length(texas_cells))),
                              "encounter_rate" = numeric(length = (length(species) * length(texas_cells))))

for (i in 1:length(species)) {
  encounter_rates$species[(length(texas_cells)*(i-1)+1):(length(texas_cells)*i)] <- species[i]
}

for (i in 1:length(texas_cells)) {
  encounter_rates$cell[c(length(texas_cells)*(0:(length(species)-1))+i)] <- texas_cells[i]
}

for (i in 1:nrow(encounter_rates)) {
  encounter_rates$encounter_rate[i] <- (nrow(filter(just_positives, 
                                                    (cell == encounter_rates$cell[i]) & 
                                                      (scientific_name == encounter_rates$species[i])))) /
    (nrow(filter(checklist_list, cell == encounter_rates$cell[i])))
}

encounter_rates$weight <- 0
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$weight[i] <- grid_cell_checklists$n_checklists[grid_cell_checklists$cell == encounter_rates$cell[i]]
}

write.csv(encounter_rates, "data/TX_encounter_rates_2000-2005.csv")

rm(ebird_TX_all_species)
rm(encounter_rates)
rm(just_positives)

# 2006 - 2008 ----

ebird_TX_all_species <- read.csv("data/ebd_TX_all_species_2006-2008_zf.csv")

just_positives <- filter(ebird_TX_all_species, species_observed == TRUE)

checklist_list <- ebird_TX_all_species %>%
  group_by(checklist_id) %>%
  summarise(cell = first(cell))

grid_cell_checklists <- checklist_list %>%
  group_by(cell) %>%
  summarise(n_checklists = n())

n_checklists <- length(unique(ebird_TX_all_species$checklist_id))

species <- unique(ebird_TX_all_species$scientific_name)
texas_cells <- unique(ebird_TX_all_species$cell)

encounter_rates <- data.frame("species" = character(length = (length(species) * length(texas_cells))), 
                              "cell" = numeric(length = (length(species) * length(texas_cells))),
                              "encounter_rate" = numeric(length = (length(species) * length(texas_cells))))

for (i in 1:length(species)) {
  encounter_rates$species[(length(texas_cells)*(i-1)+1):(length(texas_cells)*i)] <- species[i]
}

for (i in 1:length(texas_cells)) {
  encounter_rates$cell[c(length(texas_cells)*(0:(length(species)-1))+i)] <- texas_cells[i]
}

for (i in 1:nrow(encounter_rates)) {
  encounter_rates$encounter_rate[i] <- (nrow(filter(just_positives, 
                                                    (cell == encounter_rates$cell[i]) & 
                                                      (scientific_name == encounter_rates$species[i])))) /
    (nrow(filter(checklist_list, cell == encounter_rates$cell[i])))
}

encounter_rates$weight <- 0
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$weight[i] <- grid_cell_checklists$n_checklists[grid_cell_checklists$cell == encounter_rates$cell[i]]
}

write.csv(encounter_rates, "data/TX_encounter_rates_2006-2008.csv")

rm(ebird_TX_all_species)
rm(encounter_rates)
rm(just_positives)

# 2009 - 2010 ----

ebird_TX_all_species <- read.csv("data/ebd_TX_all_species_2009-2010_zf.csv")

just_positives <- filter(ebird_TX_all_species, species_observed == TRUE)

checklist_list <- ebird_TX_all_species %>%
  group_by(checklist_id) %>%
  summarise(cell = first(cell))

grid_cell_checklists <- checklist_list %>%
  group_by(cell) %>%
  summarise(n_checklists = n())

n_checklists <- length(unique(ebird_TX_all_species$checklist_id))

species <- unique(ebird_TX_all_species$scientific_name)
texas_cells <- unique(ebird_TX_all_species$cell)

encounter_rates <- data.frame("species" = character(length = (length(species) * length(texas_cells))), 
                              "cell" = numeric(length = (length(species) * length(texas_cells))),
                              "encounter_rate" = numeric(length = (length(species) * length(texas_cells))))

for (i in 1:length(species)) {
  encounter_rates$species[(length(texas_cells)*(i-1)+1):(length(texas_cells)*i)] <- species[i]
}

for (i in 1:length(texas_cells)) {
  encounter_rates$cell[c(length(texas_cells)*(0:(length(species)-1))+i)] <- texas_cells[i]
}

for (i in 1:nrow(encounter_rates)) {
  encounter_rates$encounter_rate[i] <- (nrow(filter(just_positives, 
                                                    (cell == encounter_rates$cell[i]) & 
                                                      (scientific_name == encounter_rates$species[i])))) /
    (nrow(filter(checklist_list, cell == encounter_rates$cell[i])))
}

encounter_rates$weight <- 0
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$weight[i] <- grid_cell_checklists$n_checklists[grid_cell_checklists$cell == encounter_rates$cell[i]]
}

write.csv(encounter_rates, "data/TX_encounter_rates_2009-2010.csv")

rm(ebird_TX_all_species)
rm(encounter_rates)
rm(just_positives)

# 2011 ----

ebird_TX_all_species <- read.csv("data/ebd_TX_all_species_2011_zf.csv")

just_positives <- filter(ebird_TX_all_species, species_observed == TRUE)

checklist_list <- ebird_TX_all_species %>%
  group_by(checklist_id) %>%
  summarise(cell = first(cell))

grid_cell_checklists <- checklist_list %>%
  group_by(cell) %>%
  summarise(n_checklists = n())

n_checklists <- length(unique(ebird_TX_all_species$checklist_id))

species <- unique(ebird_TX_all_species$scientific_name)
texas_cells <- unique(ebird_TX_all_species$cell)

encounter_rates <- data.frame("species" = character(length = (length(species) * length(texas_cells))), 
                              "cell" = numeric(length = (length(species) * length(texas_cells))),
                              "encounter_rate" = numeric(length = (length(species) * length(texas_cells))))

for (i in 1:length(species)) {
  encounter_rates$species[(length(texas_cells)*(i-1)+1):(length(texas_cells)*i)] <- species[i]
}

for (i in 1:length(texas_cells)) {
  encounter_rates$cell[c(length(texas_cells)*(0:(length(species)-1))+i)] <- texas_cells[i]
}

for (i in 1:nrow(encounter_rates)) {
  encounter_rates$encounter_rate[i] <- (nrow(filter(just_positives, 
                                                    (cell == encounter_rates$cell[i]) & 
                                                      (scientific_name == encounter_rates$species[i])))) /
    (nrow(filter(checklist_list, cell == encounter_rates$cell[i])))
}

encounter_rates$weight <- 0
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$weight[i] <- grid_cell_checklists$n_checklists[grid_cell_checklists$cell == encounter_rates$cell[i]]
}

write.csv(encounter_rates, "data/TX_encounter_rates_2011.csv")

rm(ebird_TX_all_species)
rm(encounter_rates)
rm(just_positives)

# 2012 ----

ebird_TX_all_species <- read.csv("data/ebd_TX_all_species_2012_zf.csv")

just_positives <- filter(ebird_TX_all_species, species_observed == TRUE)

checklist_list <- ebird_TX_all_species %>%
  group_by(checklist_id) %>%
  summarise(cell = first(cell))

grid_cell_checklists <- checklist_list %>%
  group_by(cell) %>%
  summarise(n_checklists = n())

n_checklists <- length(unique(ebird_TX_all_species$checklist_id))

species <- unique(ebird_TX_all_species$scientific_name)
texas_cells <- unique(ebird_TX_all_species$cell)

encounter_rates <- data.frame("species" = character(length = (length(species) * length(texas_cells))), 
                              "cell" = numeric(length = (length(species) * length(texas_cells))),
                              "encounter_rate" = numeric(length = (length(species) * length(texas_cells))))

for (i in 1:length(species)) {
  encounter_rates$species[(length(texas_cells)*(i-1)+1):(length(texas_cells)*i)] <- species[i]
}

for (i in 1:length(texas_cells)) {
  encounter_rates$cell[c(length(texas_cells)*(0:(length(species)-1))+i)] <- texas_cells[i]
}

for (i in 1:nrow(encounter_rates)) {
  encounter_rates$encounter_rate[i] <- (nrow(filter(just_positives, 
                                                    (cell == encounter_rates$cell[i]) & 
                                                      (scientific_name == encounter_rates$species[i])))) /
    (nrow(filter(checklist_list, cell == encounter_rates$cell[i])))
}

encounter_rates$weight <- 0
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$weight[i] <- grid_cell_checklists$n_checklists[grid_cell_checklists$cell == encounter_rates$cell[i]]
}

write.csv(encounter_rates, "data/TX_encounter_rates_2012.csv")

rm(ebird_TX_all_species)
rm(encounter_rates)
rm(just_positives)

# 2013----

ebird_TX_all_species <- read.csv("data/ebd_TX_all_species_2013_zf.csv")

just_positives <- filter(ebird_TX_all_species, species_observed == TRUE)

checklist_list <- ebird_TX_all_species %>%
  group_by(checklist_id) %>%
  summarise(cell = first(cell))

grid_cell_checklists <- checklist_list %>%
  group_by(cell) %>%
  summarise(n_checklists = n())

n_checklists <- length(unique(ebird_TX_all_species$checklist_id))

species <- unique(ebird_TX_all_species$scientific_name)
texas_cells <- unique(ebird_TX_all_species$cell)

encounter_rates <- data.frame("species" = character(length = (length(species) * length(texas_cells))), 
                              "cell" = numeric(length = (length(species) * length(texas_cells))),
                              "encounter_rate" = numeric(length = (length(species) * length(texas_cells))))

for (i in 1:length(species)) {
  encounter_rates$species[(length(texas_cells)*(i-1)+1):(length(texas_cells)*i)] <- species[i]
}

for (i in 1:length(texas_cells)) {
  encounter_rates$cell[c(length(texas_cells)*(0:(length(species)-1))+i)] <- texas_cells[i]
}

for (i in 1:nrow(encounter_rates)) {
  encounter_rates$encounter_rate[i] <- (nrow(filter(just_positives, 
                                                    (cell == encounter_rates$cell[i]) & 
                                                      (scientific_name == encounter_rates$species[i])))) /
    (nrow(filter(checklist_list, cell == encounter_rates$cell[i])))
}

encounter_rates$weight <- 0
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$weight[i] <- grid_cell_checklists$n_checklists[grid_cell_checklists$cell == encounter_rates$cell[i]]
}

write.csv(encounter_rates, "data/TX_encounter_rates_2013.csv")

rm(ebird_TX_all_species)
rm(encounter_rates)
rm(just_positives)

# 2014 ----

ebird_TX_all_species <- read.csv("data/ebd_TX_all_species_2014_zf.csv")

just_positives <- filter(ebird_TX_all_species, species_observed == TRUE)

checklist_list <- ebird_TX_all_species %>%
  group_by(checklist_id) %>%
  summarise(cell = first(cell))

grid_cell_checklists <- checklist_list %>%
  group_by(cell) %>%
  summarise(n_checklists = n())

n_checklists <- length(unique(ebird_TX_all_species$checklist_id))

species <- unique(ebird_TX_all_species$scientific_name)
texas_cells <- unique(ebird_TX_all_species$cell)

encounter_rates <- data.frame("species" = character(length = (length(species) * length(texas_cells))), 
                              "cell" = numeric(length = (length(species) * length(texas_cells))),
                              "encounter_rate" = numeric(length = (length(species) * length(texas_cells))))

for (i in 1:length(species)) {
  encounter_rates$species[(length(texas_cells)*(i-1)+1):(length(texas_cells)*i)] <- species[i]
}

for (i in 1:length(texas_cells)) {
  encounter_rates$cell[c(length(texas_cells)*(0:(length(species)-1))+i)] <- texas_cells[i]
}

for (i in 1:nrow(encounter_rates)) {
  encounter_rates$encounter_rate[i] <- (nrow(filter(just_positives, 
                                                    (cell == encounter_rates$cell[i]) & 
                                                      (scientific_name == encounter_rates$species[i])))) /
    (nrow(filter(checklist_list, cell == encounter_rates$cell[i])))
}

encounter_rates$weight <- 0
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$weight[i] <- grid_cell_checklists$n_checklists[grid_cell_checklists$cell == encounter_rates$cell[i]]
}

write.csv(encounter_rates, "data/TX_encounter_rates_2014.csv")

rm(ebird_TX_all_species)
rm(encounter_rates)
rm(just_positives)

# 2015 ----

ebird_TX_all_species <- read.csv("data/ebd_TX_all_species_2015_zf.csv")

just_positives <- filter(ebird_TX_all_species, species_observed == TRUE)

checklist_list <- ebird_TX_all_species %>%
  group_by(checklist_id) %>%
  summarise(cell = first(cell))

grid_cell_checklists <- checklist_list %>%
  group_by(cell) %>%
  summarise(n_checklists = n())

n_checklists <- length(unique(ebird_TX_all_species$checklist_id))

species <- unique(ebird_TX_all_species$scientific_name)
texas_cells <- unique(ebird_TX_all_species$cell)

encounter_rates <- data.frame("species" = character(length = (length(species) * length(texas_cells))), 
                              "cell" = numeric(length = (length(species) * length(texas_cells))),
                              "encounter_rate" = numeric(length = (length(species) * length(texas_cells))))

for (i in 1:length(species)) {
  encounter_rates$species[(length(texas_cells)*(i-1)+1):(length(texas_cells)*i)] <- species[i]
}

for (i in 1:length(texas_cells)) {
  encounter_rates$cell[c(length(texas_cells)*(0:(length(species)-1))+i)] <- texas_cells[i]
}

for (i in 1:nrow(encounter_rates)) {
  encounter_rates$encounter_rate[i] <- (nrow(filter(just_positives, 
                                                    (cell == encounter_rates$cell[i]) & 
                                                      (scientific_name == encounter_rates$species[i])))) /
    (nrow(filter(checklist_list, cell == encounter_rates$cell[i])))
}

encounter_rates$weight <- 0
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$weight[i] <- grid_cell_checklists$n_checklists[grid_cell_checklists$cell == encounter_rates$cell[i]]
}

write.csv(encounter_rates, "data/TX_encounter_rates_2015.csv")

rm(ebird_TX_all_species)
rm(encounter_rates)
rm(just_positives)
