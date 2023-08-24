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

# import 1st dataset

ebird_FL_all_species_2000_2004 <- read.csv("data/ebd_FL_all_species_2000-2004_zf.csv")

# get rid of duplicate entries - there don't seem to be any so this step was removed

species <- unique(ebird_FL_all_species_2000_2004$scientific_name)

dggs <- dgconstruct(res = 8, projection = "ISEA", metric = TRUE, 
                    resround = 'nearest')

florida_cells <- unique(ebird_FL_all_species_2000_2004$cell)

encounter_rates <- data.frame("species" = character(length = (length(species) * length(florida_cells))), 
                              "cell" = numeric(length = (length(species) * length(florida_cells))),
                              "encounter_rate" = numeric(length = (length(species) * length(florida_cells))))

# fill cells and species names for df

for (i in 1:length(species)) {
  encounter_rates$species[(length(florida_cells)*(i-1)+1):(length(florida_cells)*i)] <- species[i]
}

for (i in 1:length(florida_cells)) {
  encounter_rates$cell[c(length(florida_cells)*(0:(length(species)-1))+i)] <- florida_cells[i]
}


# calculate encounter rate

for (i in 1:length(species)) {
  for (j in 1:length(florida_cells)) {
    subject_cell_bird <- subset(ebird_FL_all_species_2000_2004, (ebird_FL_all_species_2000_2004$cell == florida_cells[j]) & (ebird_FL_all_species_2000_2004$scientific_name == species[i]))
    
    encounter_rates$encounter_rate[which((encounter_rates$species == species[i]) & (encounter_rates$cell == florida_cells[j]))] <- mean(subject_cell_bird$species_observed)
  }
}

check <- filter(encounter_rates, encounter_rates$cell == 33045)


# add in weight metric for data this should actually be weight per cell

encounter_rates$weight <- 0

for (i in 1:nrow(encounter_rates)) {
  encounter_rates$weight[i] <- length(unique(filter(ebird_FL_all_species_2000_2004, 
                                                    cell == encounter_rates$cell[i])$checklist_id))
}

# save dataframe

write.csv(encounter_rates, "data/FL_encounter_rates_2000_2004.csv")

# delete original dataframe to make room

rm(ebird_FL_all_species_2000_2004)
rm(encounter_rates)

## Repeat process with other datasets - starting with 05-08 ----

ebird_FL_all_species <- read.csv("data/ebd_FL_all_species_2005-2008_zf.csv")
species <- unique(ebird_FL_all_species$scientific_name)
florida_cells <- unique(ebird_FL_all_species$cell)

encounter_rates <- data.frame("species" = character(length = (length(species) * length(florida_cells))), 
                              "cell" = numeric(length = (length(species) * length(florida_cells))),
                              "encounter_rate" = numeric(length = (length(species) * length(florida_cells))))

for (i in 1:length(species)) {
  encounter_rates$species[(length(florida_cells)*(i-1)+1):(length(florida_cells)*i)] <- species[i]
}

for (i in 1:length(florida_cells)) {
  encounter_rates$cell[c(length(florida_cells)*(0:(length(species)-1))+i)] <- florida_cells[i]
}


for (i in 1:length(species)) {
  for (j in 1:length(florida_cells)) {
    subject_cell_bird <- subset(ebird_FL_all_species, (ebird_FL_all_species$cell == florida_cells[j]) & (ebird_FL_all_species$scientific_name == species[i]))
    
    encounter_rates$encounter_rate[which((encounter_rates$species == species[i]) & (encounter_rates$cell == florida_cells[j]))] <- mean(subject_cell_bird$species_observed)
  }
}

encounter_rates$weight <- 0
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$weight[i] <- length(unique(filter(ebird_FL_all_species, 
                                                    cell == encounter_rates$cell[i])$checklist_id))
}

write.csv(encounter_rates, "data/FL_encounter_rates_2005_2008.csv")

rm(ebird_FL_all_species)
rm(encounter_rates)

# 2009 - 2010 ----

ebird_FL_all_species <- read.csv("data/ebd_FL_all_species_2009-2010_zf.csv")
species <- unique(ebird_FL_all_species$scientific_name)
florida_cells <- unique(ebird_FL_all_species$cell)

encounter_rates <- data.frame("species" = character(length = (length(species) * length(florida_cells))), 
                              "cell" = numeric(length = (length(species) * length(florida_cells))),
                              "encounter_rate" = numeric(length = (length(species) * length(florida_cells))))

for (i in 1:length(species)) {
  encounter_rates$species[(length(florida_cells)*(i-1)+1):(length(florida_cells)*i)] <- species[i]
}

for (i in 1:length(florida_cells)) {
  encounter_rates$cell[c(length(florida_cells)*(0:(length(species)-1))+i)] <- florida_cells[i]
}


for (i in 1:length(species)) {
  for (j in 1:length(florida_cells)) {
    subject_cell_bird <- subset(ebird_FL_all_species, (ebird_FL_all_species$cell == florida_cells[j]) & (ebird_FL_all_species$scientific_name == species[i]))
    
    encounter_rates$encounter_rate[which((encounter_rates$species == species[i]) & (encounter_rates$cell == florida_cells[j]))] <- mean(subject_cell_bird$species_observed)
  }
}

encounter_rates$weight <- 0
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$weight[i] <- length(unique(filter(ebird_FL_all_species, 
                                                    cell == encounter_rates$cell[i])$checklist_id))
}

write.csv(encounter_rates, "data/FL_encounter_rates_2009_2010.csv")

rm(ebird_FL_all_species)
rm(encounter_rates)

# 2011 ----

ebird_FL_all_species <- read.csv("data/ebd_FL_all_species_2011_zf.csv")
species <- unique(ebird_FL_all_species$scientific_name)
florida_cells <- unique(ebird_FL_all_species$cell)

encounter_rates <- data.frame("species" = character(length = (length(species) * length(florida_cells))), 
                              "cell" = numeric(length = (length(species) * length(florida_cells))),
                              "encounter_rate" = numeric(length = (length(species) * length(florida_cells))))

for (i in 1:length(species)) {
  encounter_rates$species[(length(florida_cells)*(i-1)+1):(length(florida_cells)*i)] <- species[i]
}

for (i in 1:length(florida_cells)) {
  encounter_rates$cell[c(length(florida_cells)*(0:(length(species)-1))+i)] <- florida_cells[i]
}


for (i in 1:length(species)) {
  for (j in 1:length(florida_cells)) {
    subject_cell_bird <- subset(ebird_FL_all_species, (ebird_FL_all_species$cell == florida_cells[j]) & (ebird_FL_all_species$scientific_name == species[i]))
    
    encounter_rates$encounter_rate[which((encounter_rates$species == species[i]) & (encounter_rates$cell == florida_cells[j]))] <- mean(subject_cell_bird$species_observed)
  }
}

encounter_rates$weight <- 0
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$weight[i] <- length(unique(filter(ebird_FL_all_species, 
                                                    cell == encounter_rates$cell[i])$checklist_id))
}

write.csv(encounter_rates, "data/FL_encounter_rates_2011.csv")

rm(ebird_FL_all_species)
rm(encounter_rates)

# 2012 ----

ebird_FL_all_species <- read.csv("data/ebd_FL_all_species_2012_zf.csv")
species <- unique(ebird_FL_all_species$scientific_name)
florida_cells <- unique(ebird_FL_all_species$cell)

encounter_rates <- data.frame("species" = character(length = (length(species) * length(florida_cells))), 
                              "cell" = numeric(length = (length(species) * length(florida_cells))),
                              "encounter_rate" = numeric(length = (length(species) * length(florida_cells))))

for (i in 1:length(species)) {
  encounter_rates$species[(length(florida_cells)*(i-1)+1):(length(florida_cells)*i)] <- species[i]
}

for (i in 1:length(florida_cells)) {
  encounter_rates$cell[c(length(florida_cells)*(0:(length(species)-1))+i)] <- florida_cells[i]
}


for (i in 1:length(species)) {
  for (j in 1:length(florida_cells)) {
    subject_cell_bird <- subset(ebird_FL_all_species, (ebird_FL_all_species$cell == florida_cells[j]) & (ebird_FL_all_species$scientific_name == species[i]))
    
    encounter_rates$encounter_rate[which((encounter_rates$species == species[i]) & (encounter_rates$cell == florida_cells[j]))] <- mean(subject_cell_bird$species_observed)
  }
}

encounter_rates$weight <- 0
for (i in 18136:nrow(encounter_rates)) {
  encounter_rates$weight[i] <- length(unique(filter(ebird_FL_all_species, 
                                                    cell == encounter_rates$cell[i])$checklist_id))
}

write.csv(encounter_rates, "FL_encounter_rates_2012.csv")

rm(ebird_FL_all_species)
rm(encounter_rates)

# 2013 ----

ebird_FL_all_species <- read.csv("data/ebd_FL_all_species_2013_zf.csv")

just_positives <- filter(ebird_FL_all_species, species_observed == TRUE)

checklist_list <- ebird_FL_all_species %>%
  group_by(checklist_id) %>%
  summarise(cell = first(cell))

grid_cell_checklists <- checklist_list %>%
  group_by(cell) %>%
  summarise(n_checklists = n())

n_checklists <- length(unique(ebird_FL_all_species$checklist_id))

species <- unique(ebird_FL_all_species$scientific_name)
florida_cells <- unique(ebird_FL_all_species$cell)

encounter_rates <- data.frame("species" = character(length = (length(species) * length(florida_cells))), 
                              "cell" = numeric(length = (length(species) * length(florida_cells))),
                              "encounter_rate" = numeric(length = (length(species) * length(florida_cells))))

for (i in 1:length(species)) {
  encounter_rates$species[(length(florida_cells)*(i-1)+1):(length(florida_cells)*i)] <- species[i]
}

for (i in 1:length(florida_cells)) {
  encounter_rates$cell[c(length(florida_cells)*(0:(length(species)-1))+i)] <- florida_cells[i]
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

write.csv(encounter_rates, "data/FL_encounter_rates_2013.csv")

rm(ebird_FL_all_species)
rm(encounter_rates)

# 2014a ----

ebird_FL_all_species <- read.csv("data/ebd_FL_all_species_2014a_zf.csv")
species <- unique(ebird_FL_all_species$scientific_name)
florida_cells <- unique(ebird_FL_all_species$cell)

encounter_rates <- data.frame("species" = character(length = (length(species) * length(florida_cells))), 
                              "cell" = numeric(length = (length(species) * length(florida_cells))),
                              "encounter_rate" = numeric(length = (length(species) * length(florida_cells))))

for (i in 1:length(species)) {
  encounter_rates$species[(length(florida_cells)*(i-1)+1):(length(florida_cells)*i)] <- species[i]
}

for (i in 1:length(florida_cells)) {
  encounter_rates$cell[c(length(florida_cells)*(0:(length(species)-1))+i)] <- florida_cells[i]
}


for (i in 1:length(species)) {
  for (j in 1:length(florida_cells)) {
    subject_cell_bird <- subset(ebird_FL_all_species, (ebird_FL_all_species$cell == florida_cells[j]) & (ebird_FL_all_species$scientific_name == species[i]))
    
    encounter_rates$encounter_rate[which((encounter_rates$species == species[i]) & (encounter_rates$cell == florida_cells[j]))] <- mean(subject_cell_bird$species_observed)
  }
}

encounter_rates$weight <- 0
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$weight[i] <- length(unique(filter(ebird_FL_all_species, 
                                                    cell == encounter_rates$cell[i])$checklist_id))
}

write.csv(encounter_rates, "data/FL_encounter_rates_2014a.csv")

rm(ebird_FL_all_species)
rm(encounter_rates)

# 2014b ----

ebird_FL_all_species <- read.csv("data/ebd_FL_all_species_2014b_zf.csv")
species <- unique(ebird_FL_all_species$scientific_name)
florida_cells <- unique(ebird_FL_all_species$cell)

encounter_rates <- data.frame("species" = character(length = (length(species) * length(florida_cells))), 
                              "cell" = numeric(length = (length(species) * length(florida_cells))),
                              "encounter_rate" = numeric(length = (length(species) * length(florida_cells))))

for (i in 1:length(species)) {
  encounter_rates$species[(length(florida_cells)*(i-1)+1):(length(florida_cells)*i)] <- species[i]
}

for (i in 1:length(florida_cells)) {
  encounter_rates$cell[c(length(florida_cells)*(0:(length(species)-1))+i)] <- florida_cells[i]
}


for (i in 1:length(species)) {
  for (j in 1:length(florida_cells)) {
    subject_cell_bird <- subset(ebird_FL_all_species, (ebird_FL_all_species$cell == florida_cells[j]) & (ebird_FL_all_species$scientific_name == species[i]))
    
    encounter_rates$encounter_rate[which((encounter_rates$species == species[i]) & (encounter_rates$cell == florida_cells[j]))] <- mean(subject_cell_bird$species_observed)
  }
}

encounter_rates$weight <- 0
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$weight[i] <- length(unique(filter(ebird_FL_all_species, 
                                                    cell == encounter_rates$cell[i])$checklist_id))
}

write.csv(encounter_rates, "data/FL_encounter_rates_2014b.csv")

rm(ebird_FL_all_species)
rm(encounter_rates)

# 2015a ----

ebird_FL_all_species <- read.csv("data/ebd_FL_all_species_2015a_zf.csv")
species <- unique(ebird_FL_all_species$scientific_name)
florida_cells <- unique(ebird_FL_all_species$cell)

encounter_rates <- data.frame("species" = character(length = (length(species) * length(florida_cells))), 
                              "cell" = numeric(length = (length(species) * length(florida_cells))),
                              "encounter_rate" = numeric(length = (length(species) * length(florida_cells))))

for (i in 1:length(species)) {
  encounter_rates$species[(length(florida_cells)*(i-1)+1):(length(florida_cells)*i)] <- species[i]
}

for (i in 1:length(florida_cells)) {
  encounter_rates$cell[c(length(florida_cells)*(0:(length(species)-1))+i)] <- florida_cells[i]
}


for (i in 1:length(species)) {
  for (j in 1:length(florida_cells)) {
    subject_cell_bird <- subset(ebird_FL_all_species, (ebird_FL_all_species$cell == florida_cells[j]) & (ebird_FL_all_species$scientific_name == species[i]))
    
    encounter_rates$encounter_rate[which((encounter_rates$species == species[i]) & (encounter_rates$cell == florida_cells[j]))] <- mean(subject_cell_bird$species_observed)
  }
}

encounter_rates$weight <- 0
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$weight[i] <- length(unique(filter(ebird_FL_all_species, 
                                                    cell == encounter_rates$cell[i])$checklist_id))
}

write.csv(encounter_rates, "data/FL_encounter_rates_2015a.csv")

rm(ebird_FL_all_species)
rm(encounter_rates)

# 2015b ----

ebird_FL_all_species <- read.csv("data/ebd_FL_all_species_2015b_zf.csv")
species <- unique(ebird_FL_all_species$scientific_name)
florida_cells <- unique(ebird_FL_all_species$cell)

encounter_rates <- data.frame("species" = character(length = (length(species) * length(florida_cells))), 
                              "cell" = numeric(length = (length(species) * length(florida_cells))),
                              "encounter_rate" = numeric(length = (length(species) * length(florida_cells))))

for (i in 1:length(species)) {
  encounter_rates$species[(length(florida_cells)*(i-1)+1):(length(florida_cells)*i)] <- species[i]
}

for (i in 1:length(florida_cells)) {
  encounter_rates$cell[c(length(florida_cells)*(0:(length(species)-1))+i)] <- florida_cells[i]
}


for (i in 1:length(species)) {
  for (j in 1:length(florida_cells)) {
    subject_cell_bird <- subset(ebird_FL_all_species, (ebird_FL_all_species$cell == florida_cells[j]) & (ebird_FL_all_species$scientific_name == species[i]))
    
    encounter_rates$encounter_rate[which((encounter_rates$species == species[i]) & (encounter_rates$cell == florida_cells[j]))] <- mean(subject_cell_bird$species_observed)
  }
}

encounter_rates$weight <- 0
for (i in 1:nrow(encounter_rates)) {
  encounter_rates$weight[i] <- length(unique(filter(ebird_FL_all_species, 
                                                    cell == encounter_rates$cell[i])$checklist_id))
}

write.csv(encounter_rates, "data/FL_encounter_rates_2015b.csv")

rm(ebird_FL_all_species)
rm(encounter_rates)
