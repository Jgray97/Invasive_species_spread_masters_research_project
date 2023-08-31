# GENERATING FLORIDA BIOTIC INPUT VARIABLE VALUES FOR RF ALGORITHM
# Author = John Gray
# Email = greyjohn15@gmail.com
# Last edit = 31/08/2023

# Loading necessary packages ----

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

## Import encounter rates dataframes ----

encounter_rates_02_04 <- read.csv("data/FL_encounter_rates_2002_2004.csv")
encounter_rates_05_08 <- read.csv("data/FL_encounter_rates_2005_2008.csv")
encounter_rates_09_10 <- read.csv("data/FL_encounter_rates_2009_2010.csv")
encounter_rates_11 <- read.csv("data/FL_encounter_rates_2011.csv")
encounter_rates_12 <- read.csv("data/FL_encounter_rates_2012.csv")
encounter_rates_13 <- read.csv("data/FL_encounter_rates_2013.csv")
encounter_rates_14a <- read.csv("data/FL_encounter_rates_2014a.csv")
encounter_rates_14b <- read.csv("data/FL_encounter_rates_2014b.csv")
encounter_rates_15a <- read.csv("data/FL_encounter_rates_2015a.csv")
encounter_rates_15b <- read.csv("data/FL_encounter_rates_2015b.csv")

# Create an overall encounter rates df ----

# Rowbind all individual dataframes
encounter_rates <- rbind(encounter_rates_02_04, encounter_rates_05_08, 
                         encounter_rates_09_10, encounter_rates_11, 
                         encounter_rates_12, encounter_rates_13,
                         encounter_rates_14a, encounter_rates_14b,
                         encounter_rates_15a, encounter_rates_15b)

# get rid of 1st unnecessary column that arises from importing csvs
encounter_rates <- encounter_rates[,-1]

# create a weighted encounter rate column to account for number of observations
# associated with each encounter rate value
encounter_rates$weighted_er <- encounter_rates$encounter_rate * encounter_rates$weight

# Create a list of all unique species observed between 2002 and 2015
species <- unique(encounter_rates$species)

# create a list of all florida cells in 2002 - 2015 data
cells <- unique(encounter_rates$cell)

# create overall encounter rates dataframe to be populated with values
er_overall <- data.frame("species" = character(length = (length(species) * length(cells))), 
                              "cell" = numeric(length = (length(species) * length(cells))),
                              "encounter_rate" = numeric(length = (length(species) * length(cells))))

# fill cells and species columns for df
for (i in 1:length(species)) {
  er_overall$species[(length(cells)*(i-1)+1):(length(cells)*i)] <- species[i]
}

for (i in 1:length(cells)) {
  er_overall$cell[c(length(cells)*(0:(length(species)-1))+i)] <- cells[i]
}

# calculate final encounter rate for each cell for each species using weighted
# encounter rates from each time period subset
for (i in 1:nrow(er_overall)) {
  subject_cell_bird <- filter(encounter_rates, ((species == er_overall$species[i]) & (cell == er_overall$cell[i])))
  
  if(nrow(subject_cell_bird) != 0) {
    total_weight <- sum(subject_cell_bird$weight)
    total_weighted_er <- sum(subject_cell_bird$weighted_er)
    er_overall$encounter_rate[i] <- total_weighted_er/total_weight
  } else {
    er_overall$encounter_rate[i] <- 0
  }
}

# get rid of all species whose detection rates are lower than 0.005, i.e if a
# species does not meet the β criterion it is not considered to be part of the
# native species assemblage
common_birds <- filter(er_overall, encounter_rate > 0.005)

# save overall encounter rates dataframe in a csv file
write.csv(common_birds, "data/common_bird_encounter_rates_2002_2015.csv")

# get rid of invasive species - these were determined manually using eBird 
# website
FL_invasive_species <- read.csv("data/invasive_species_list.csv")
common_native_birds <- filter(common_birds, !(species %in% FL_invasive_species$species))

# Save encounter rates for just native birds table
write.csv(common_native_birds, "data/common_native_bird_encounter_rates_2002_2015.csv")

# BIOTIC INPUT PARAMETERS CALCULATION ----

# comment below line out if continuing analysis from code above
common_native_birds <- read.csv("data/common_native_bird_encounter_rates_2002_2015.csv")

# Load in establishment tracker csv (as nearby encounter rates are used to 
# calculate biotic input parameters)
establishment_tracker <- read.csv("data/FINAL-VERSION-FL_egoose_establishment_ct-50_met-0.005.csv")

# create biotic input parameter df

rf_biotic_data <- filter(establishment_tracker, establishment_tracker$year > 2002)[,c(1,2,3,6)]

rf_biotic_data$bs_ring_1 <- 0

rf_biotic_data$bs_ring_2 <- 0

# THROUGH MANUAL INVESTIGATION I HAVE FOUND THAT:

# 1st/adjacent ring cell values given central cell ID = n: n+1, n-1, n+81, n+82,
# n-81, n-82

# 2nd/1-cell-separated ring cell values for central cell ID = n: n+2, n-2, n+80, 
# n+83, n-83, n-80, n-162, n-163, n-164, n+162, n+163, n+164

# create list of grid cells in florida for for loop
florida_cells <- unique(rf_biotic_data$grid_cell)

### 1st ring biotic input variable is defined as the sum of establishment score
### in the previous year, multiplied by the proportion of unique species in cell 
### and target cell which are present in both cells for all cells adjacent to
### the target cell

for (i in 1:length(florida_cells)) {
  subject_cell <- filter(common_native_birds, cell == florida_cells[i])
  
  n_plus_1 <- filter(common_native_birds, cell == (florida_cells[i]+1))
  
  n_minus_1 <- filter(common_native_birds, cell == (florida_cells[i]-1))
  
  n_plus_81 <- filter(common_native_birds, cell == (florida_cells[i]+81))
  
  n_plus_82 <- filter(common_native_birds, cell == (florida_cells[i]+82))
  
  n_minus_81 <- filter(common_native_birds, cell == (florida_cells[i]-81))
  
  n_minus_82 <- filter(common_native_birds, cell == (florida_cells[i]-82))
  
  n_plus_1_similarity <- (length(subject_cell$species) + length(n_plus_1$species) - length(unique(rbind(subject_cell,n_plus_1)$species)))/length(unique(rbind(subject_cell,n_plus_1)$species))
  
  n_minus_1_similarity <- (length(subject_cell$species) + length(n_minus_1$species) - length(unique(rbind(subject_cell,n_minus_1)$species)))/length(unique(rbind(subject_cell,n_minus_1)$species))
  
  n_plus_81_similarity <- (length(subject_cell$species) + length(n_plus_81$species) - length(unique(rbind(subject_cell,n_plus_81)$species)))/length(unique(rbind(subject_cell,n_plus_81)$species))
  
  n_plus_82_similarity <- (length(subject_cell$species) + length(n_plus_82$species) - length(unique(rbind(subject_cell,n_plus_82)$species)))/length(unique(rbind(subject_cell,n_plus_82)$species))
  
  n_minus_81_similarity <- (length(subject_cell$species) + length(n_minus_81$species) - length(unique(rbind(subject_cell,n_minus_81)$species)))/length(unique(rbind(subject_cell,n_minus_81)$species))
  
  n_minus_82_similarity <- (length(subject_cell$species) + length(n_minus_82$species) - length(unique(rbind(subject_cell,n_minus_82)$species)))/length(unique(rbind(subject_cell,n_minus_82)$species))
  
  for (j in 2003:2019) {
    rf_biotic_data$bs_ring_1[(rf_biotic_data$year == j) & 
                                     (rf_biotic_data$grid_cell == florida_cells[i])] <- 
      (if((florida_cells[i] + 1) %in% florida_cells) {
        (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                     (establishment_tracker$grid_cell == (florida_cells[i]+1))] * n_plus_1_similarity)
      } else {
        0
      }) + 
      (if((florida_cells[i] - 1) %in% florida_cells) {
      (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                   (establishment_tracker$grid_cell == (florida_cells[i]-1))] * n_minus_1_similarity)
        } else {
          0
        }) +
      (if((florida_cells[i] + 81) %in% florida_cells) {
      (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                   (establishment_tracker$grid_cell == (florida_cells[i]+81))] * n_plus_81_similarity)
        } else {
          0
        }) +
      (if((florida_cells[i] + 82) %in% florida_cells) {
      (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                   (establishment_tracker$grid_cell == (florida_cells[i]+82))] * n_plus_82_similarity)
        } else {
          0
        }) +
      (if((florida_cells[i] - 81) %in% florida_cells) {
      (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                   (establishment_tracker$grid_cell == (florida_cells[i]-81))] * n_minus_81_similarity)
        } else {
          0
        }) +
      (if((florida_cells[i] - 82) %in% florida_cells) {
      (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                   (establishment_tracker$grid_cell == (florida_cells[i]-82))] * n_minus_82_similarity)
      } else {
        0})
      }
}

## 2nd ring biotic input variable is the same except sum is over cells ----
## which are 1 degree separated from the target cell ----

for (i in 1:length(florida_cells)) {
  subject_cell <- filter(common_native_birds, cell == florida_cells[i])
  
  n_plus_2 <- filter(common_native_birds, cell == (florida_cells[i]+2))
  
  n_minus_2 <- filter(common_native_birds, cell == (florida_cells[i]-2))
  
  n_plus_80 <- filter(common_native_birds, cell == (florida_cells[i]+80))
  
  n_plus_83 <- filter(common_native_birds, cell == (florida_cells[i]+83))
  
  n_minus_80 <- filter(common_native_birds, cell == (florida_cells[i]-80))
  
  n_minus_83 <- filter(common_native_birds, cell == (florida_cells[i]-83))
  
  n_plus_162 <- filter(common_native_birds, cell == (florida_cells[i]+162))
  
  n_plus_163 <- filter(common_native_birds, cell == (florida_cells[i]+163))
  
  n_plus_164 <- filter(common_native_birds, cell == (florida_cells[i]+164))
  
  n_minus_162 <- filter(common_native_birds, cell == (florida_cells[i]-162))
  
  n_minus_163 <- filter(common_native_birds, cell == (florida_cells[i]-163))
  
  n_minus_164 <- filter(common_native_birds, cell == (florida_cells[i]-164))
  
  n_plus_2_similarity <- (length(subject_cell$species) + length(n_plus_2$species) - length(unique(rbind(subject_cell,n_plus_2)$species)))/length(unique(rbind(subject_cell,n_plus_2)$species))
  
  n_minus_2_similarity <- (length(subject_cell$species) + length(n_minus_2$species) - length(unique(rbind(subject_cell,n_minus_2)$species)))/length(unique(rbind(subject_cell,n_minus_2)$species))
  
  n_plus_80_similarity <- (length(subject_cell$species) + length(n_plus_80$species) - length(unique(rbind(subject_cell,n_plus_80)$species)))/length(unique(rbind(subject_cell,n_plus_80)$species))
  
  n_plus_83_similarity <- (length(subject_cell$species) + length(n_plus_83$species) - length(unique(rbind(subject_cell,n_plus_83)$species)))/length(unique(rbind(subject_cell,n_plus_83)$species))
  
  n_minus_80_similarity <- (length(subject_cell$species) + length(n_minus_80$species) - length(unique(rbind(subject_cell,n_minus_80)$species)))/length(unique(rbind(subject_cell,n_minus_80)$species))
  
  n_minus_83_similarity <- (length(subject_cell$species) + length(n_minus_83$species) - length(unique(rbind(subject_cell,n_minus_83)$species)))/length(unique(rbind(subject_cell,n_minus_83)$species))
  
  n_plus_162_similarity <- (length(subject_cell$species) + length(n_plus_162$species) - length(unique(rbind(subject_cell,n_plus_162)$species)))/length(unique(rbind(subject_cell,n_plus_162)$species))
  
  n_plus_163_similarity <- (length(subject_cell$species) + length(n_plus_163$species) - length(unique(rbind(subject_cell,n_plus_163)$species)))/length(unique(rbind(subject_cell,n_plus_163)$species))
  
  n_plus_164_similarity <- (length(subject_cell$species) + length(n_plus_164$species) - length(unique(rbind(subject_cell,n_plus_164)$species)))/length(unique(rbind(subject_cell,n_plus_164)$species))
  
  n_minus_162_similarity <- (length(subject_cell$species) + length(n_minus_162$species) - length(unique(rbind(subject_cell,n_minus_162)$species)))/length(unique(rbind(subject_cell,n_minus_162)$species))
  
  n_minus_163_similarity <- (length(subject_cell$species) + length(n_minus_163$species) - length(unique(rbind(subject_cell,n_minus_163)$species)))/length(unique(rbind(subject_cell,n_minus_163)$species))
  
  n_minus_164_similarity <- (length(subject_cell$species) + length(n_minus_164$species) - length(unique(rbind(subject_cell,n_minus_164)$species)))/length(unique(rbind(subject_cell,n_minus_164)$species))
  
  for (j in 2003:2019) {
    rf_biotic_data$bs_ring_2[(rf_biotic_data$year == j) & 
                               (rf_biotic_data$grid_cell == florida_cells[i])] <- 
      (if((florida_cells[i] + 2) %in% florida_cells) {
        (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                     (establishment_tracker$grid_cell == (florida_cells[i]+2))] * n_plus_2_similarity)
      } else {
        0
      }) + 
      (if((florida_cells[i] - 2) %in% florida_cells) {
        (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                     (establishment_tracker$grid_cell == (florida_cells[i]-2))] * n_minus_2_similarity)
      } else {
        0
      }) +
      (if((florida_cells[i] + 80) %in% florida_cells) {
        (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                     (establishment_tracker$grid_cell == (florida_cells[i]+80))] * n_plus_80_similarity)
      } else {
        0
      }) +
      (if((florida_cells[i] + 83) %in% florida_cells) {
        (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                     (establishment_tracker$grid_cell == (florida_cells[i]+83))] * n_plus_83_similarity)
      } else {
        0
      }) +
      (if((florida_cells[i] - 80) %in% florida_cells) {
        (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                     (establishment_tracker$grid_cell == (florida_cells[i]-80))] * n_minus_80_similarity)
      } else {
        0
      }) +
      (if((florida_cells[i] - 83) %in% florida_cells) {
        (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                     (establishment_tracker$grid_cell == (florida_cells[i]-83))] * n_minus_83_similarity)
      } else {
        0
        }) +
      (if((florida_cells[i] + 162) %in% florida_cells) {
        (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                     (establishment_tracker$grid_cell == (florida_cells[i]+162))] * n_plus_162_similarity)
      } else {
        0
      }) + 
      (if((florida_cells[i] + 163) %in% florida_cells) {
        (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                     (establishment_tracker$grid_cell == (florida_cells[i]+163))] * n_plus_163_similarity)
      } else {
        0
      }) +
      (if((florida_cells[i] + 164) %in% florida_cells) {
        (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                     (establishment_tracker$grid_cell == (florida_cells[i]+164))] * n_plus_164_similarity)
      } else {
        0
      }) +
      (if((florida_cells[i] - 162) %in% florida_cells) {
        (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                     (establishment_tracker$grid_cell == (florida_cells[i]-162))] * n_minus_162_similarity)
      } else {
        0
      }) +
      (if((florida_cells[i] - 163) %in% florida_cells) {
        (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                     (establishment_tracker$grid_cell == (florida_cells[i]-163))] * n_minus_163_similarity)
      } else {
        0
      }) +
      (if((florida_cells[i] - 164) %in% florida_cells) {
        (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                     (establishment_tracker$grid_cell == (florida_cells[i]-164))] * n_minus_164_similarity)
      } else {
        0
      })
  }
}

# don't need to worry about na's because they'll get factored out anyway

# save rf biotic inputs for Florida to a csv file
write.csv(rf_biotic_data, "data/rf_biotic_data.csv")


