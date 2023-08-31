# GENERATING TEXAS BIOTIC INPUT VARIABLE VALUES FOR RF ALGORITHM
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

# import encounter rates dataframes ----

encounter_rates_02_05 <- read.csv("data/TX_encounter_rates_2002-2005.csv")
encounter_rates_06_08 <- read.csv("data/TX_encounter_rates_2006-2008.csv")
encounter_rates_09_10 <- read.csv("data/TX_encounter_rates_2009-2010.csv")
encounter_rates_11 <- read.csv("data/TX_encounter_rates_2011.csv")
encounter_rates_12 <- read.csv("data/TX_encounter_rates_2012.csv")
encounter_rates_13 <- read.csv("data/TX_encounter_rates_2013.csv")
encounter_rates_14 <- read.csv("data/TX_encounter_rates_2014.csv")
encounter_rates_15 <- read.csv("data/TX_encounter_rates_2015.csv")

# Create an overall encounter rates df ----

# Rowbind all individual dataframes
encounter_rates <- rbind(encounter_rates_02_05, encounter_rates_06_08, 
                         encounter_rates_09_10, encounter_rates_11, 
                         encounter_rates_12, encounter_rates_13,
                         encounter_rates_14, encounter_rates_15)

# get rid of first unnecessary column that arises from importing csvs
encounter_rates <- encounter_rates[,-1]

# create a weighted encounter rate column to account for number of observations
# associated with each encounter rate value
encounter_rates$weighted_er <- encounter_rates$encounter_rate * encounter_rates$weight

# Create a list of all unique species observed between 2002 and 2015
species <- unique(encounter_rates$species)

# create a list of all texas cells in 2002 - 2015 data
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
write.csv(common_birds, "data/texas_common_bird_encounter_rates_2002_2015.csv")

# get rid of invasive species - these were determined manually using eBird 
# website
TX_invasive_species <- read.csv("data/Texas_invasives.csv")
common_native_birds <- filter(common_birds, !(species %in% TX_invasive_species$Species))

# Save encounter rates for just native birds table
write.csv(common_native_birds, "data/texas_common_native_bird_encounter_rates_2002_2015.csv")

# BIOTIC INPUT PARAMETERS CALCULATION ----

# comment below line out if continuing analysis from code above
common_native_birds <- read.csv("data/texas_common_native_bird_encounter_rates_2002_2015.csv")


# Load in establishment tracker csv (as nearby encounter rates are used to 
# calculate biotic input parameters)
establishment_tracker <- read.csv("data/TX_egoose_establishment_ct-50_met-0.005.csv")

# 2 scores will be defined by sum of population size multiplied by shared 
# species in inner ring
# + the same for the outer ring

# create biotic input parameter df

rf_biotic_data <- filter(establishment_tracker, establishment_tracker$year > 2009)[,c(1,2,3,6)]

rf_biotic_data$bs_ring_1 <- 0

rf_biotic_data$bs_ring_2 <- 0

# THROUGH MANUAL INVESTIGATION I HAVE FOUND THAT:

# 1st/adjacent ring cell values given central cell ID = n: n+1, n-1, n+80, n+81, 
# n-80, n-81

# 2nd/1-cell-separated ring cell values for central cell ID = n: n+2, n-2, n+79, 
# n+82, n-82, n-79, n-160, n-161, n-162, n+160, n+161, n+162

# create list of grid cells in texas for for loop
texas_cells <- unique(rf_biotic_data$grid_cell)

### 1st ring biotic input variable is defined as the sum of establishment score
### in the previous year, multiplied by the proportion of unique species in cell 
### and target cell which are present in both cells for all cells adjacent to
### the target cell

for (i in 1:length(texas_cells)) {
  subject_cell <- filter(common_native_birds, cell == texas_cells[i])
  
  n_plus_1 <- filter(common_native_birds, cell == (texas_cells[i]+1))
  
  n_minus_1 <- filter(common_native_birds, cell == (texas_cells[i]-1))
  
  n_plus_80 <- filter(common_native_birds, cell == (texas_cells[i]+80))
  
  n_plus_81 <- filter(common_native_birds, cell == (texas_cells[i]+81))
  
  n_minus_80 <- filter(common_native_birds, cell == (texas_cells[i]-80))
  
  n_minus_81 <- filter(common_native_birds, cell == (texas_cells[i]-81))
  
  n_plus_1_similarity <- (length(subject_cell$species) + length(n_plus_1$species) - length(unique(rbind(subject_cell,n_plus_1)$species)))/length(unique(rbind(subject_cell,n_plus_1)$species))
  
  n_minus_1_similarity <- (length(subject_cell$species) + length(n_minus_1$species) - length(unique(rbind(subject_cell,n_minus_1)$species)))/length(unique(rbind(subject_cell,n_minus_1)$species))
  
  n_plus_80_similarity <- (length(subject_cell$species) + length(n_plus_80$species) - length(unique(rbind(subject_cell,n_plus_80)$species)))/length(unique(rbind(subject_cell,n_plus_80)$species))
  
  n_plus_81_similarity <- (length(subject_cell$species) + length(n_plus_81$species) - length(unique(rbind(subject_cell,n_plus_81)$species)))/length(unique(rbind(subject_cell,n_plus_81)$species))
  
  n_minus_80_similarity <- (length(subject_cell$species) + length(n_minus_80$species) - length(unique(rbind(subject_cell,n_minus_80)$species)))/length(unique(rbind(subject_cell,n_minus_80)$species))
  
  n_minus_81_similarity <- (length(subject_cell$species) + length(n_minus_81$species) - length(unique(rbind(subject_cell,n_minus_81)$species)))/length(unique(rbind(subject_cell,n_minus_81)$species))
  
  for (j in 2003:2019) {
    rf_biotic_data$bs_ring_1[(rf_biotic_data$year == j) & 
                                     (rf_biotic_data$grid_cell == texas_cells[i])] <- 
      (if((texas_cells[i] + 1) %in% texas_cells) {
        (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                     (establishment_tracker$grid_cell == (texas_cells[i]+1))] * n_plus_1_similarity)
      } else {
        0
      }) + 
      (if((texas_cells[i] - 1) %in% texas_cells) {
      (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                   (establishment_tracker$grid_cell == (texas_cells[i]-1))] * n_minus_1_similarity)
        } else {
          0
        }) +
      (if((texas_cells[i] + 80) %in% texas_cells) {
      (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                   (establishment_tracker$grid_cell == (texas_cells[i]+80))] * n_plus_80_similarity)
        } else {
          0
        }) +
      (if((texas_cells[i] + 81) %in% texas_cells) {
      (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                   (establishment_tracker$grid_cell == (texas_cells[i]+81))] * n_plus_81_similarity)
        } else {
          0
        }) +
      (if((texas_cells[i] - 80) %in% texas_cells) {
      (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                   (establishment_tracker$grid_cell == (texas_cells[i]-80))] * n_minus_80_similarity)
        } else {
          0
        }) +
      (if((texas_cells[i] - 81) %in% texas_cells) {
      (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                   (establishment_tracker$grid_cell == (texas_cells[i]-81))] * n_minus_81_similarity)
      } else {
        0})
      }
}

## 2nd ring biotic input variable is the same except sum is over cells ----
## which are 1 degree separated from the target cell ----

for (i in 1:length(texas_cells)) {
  subject_cell <- filter(common_native_birds, cell == texas_cells[i])
  
  n_plus_2 <- filter(common_native_birds, cell == (texas_cells[i]+2))
  
  n_minus_2 <- filter(common_native_birds, cell == (texas_cells[i]-2))
  
  n_plus_79 <- filter(common_native_birds, cell == (texas_cells[i]+79))
  
  n_plus_82 <- filter(common_native_birds, cell == (texas_cells[i]+82))
  
  n_minus_79 <- filter(common_native_birds, cell == (texas_cells[i]-79))
  
  n_minus_82 <- filter(common_native_birds, cell == (texas_cells[i]-82))
  
  n_plus_160 <- filter(common_native_birds, cell == (texas_cells[i]+160))
  
  n_plus_161 <- filter(common_native_birds, cell == (texas_cells[i]+161))
  
  n_plus_162 <- filter(common_native_birds, cell == (texas_cells[i]+162))
  
  n_minus_160 <- filter(common_native_birds, cell == (texas_cells[i]-160))
  
  n_minus_161 <- filter(common_native_birds, cell == (texas_cells[i]-161))
  
  n_minus_162 <- filter(common_native_birds, cell == (texas_cells[i]-162))
  
  n_plus_2_similarity <- (length(subject_cell$species) + length(n_plus_2$species) - length(unique(rbind(subject_cell,n_plus_2)$species)))/length(unique(rbind(subject_cell,n_plus_2)$species))
  
  n_minus_2_similarity <- (length(subject_cell$species) + length(n_minus_2$species) - length(unique(rbind(subject_cell,n_minus_2)$species)))/length(unique(rbind(subject_cell,n_minus_2)$species))
  
  n_plus_79_similarity <- (length(subject_cell$species) + length(n_plus_79$species) - length(unique(rbind(subject_cell,n_plus_79)$species)))/length(unique(rbind(subject_cell,n_plus_79)$species))
  
  n_plus_82_similarity <- (length(subject_cell$species) + length(n_plus_82$species) - length(unique(rbind(subject_cell,n_plus_82)$species)))/length(unique(rbind(subject_cell,n_plus_82)$species))
  
  n_minus_79_similarity <- (length(subject_cell$species) + length(n_minus_79$species) - length(unique(rbind(subject_cell,n_minus_79)$species)))/length(unique(rbind(subject_cell,n_minus_79)$species))
  
  n_minus_82_similarity <- (length(subject_cell$species) + length(n_minus_82$species) - length(unique(rbind(subject_cell,n_minus_82)$species)))/length(unique(rbind(subject_cell,n_minus_82)$species))
  
  n_plus_160_similarity <- (length(subject_cell$species) + length(n_plus_160$species) - length(unique(rbind(subject_cell,n_plus_160)$species)))/length(unique(rbind(subject_cell,n_plus_160)$species))
  
  n_plus_161_similarity <- (length(subject_cell$species) + length(n_plus_161$species) - length(unique(rbind(subject_cell,n_plus_161)$species)))/length(unique(rbind(subject_cell,n_plus_161)$species))
  
  n_plus_162_similarity <- (length(subject_cell$species) + length(n_plus_162$species) - length(unique(rbind(subject_cell,n_plus_162)$species)))/length(unique(rbind(subject_cell,n_plus_162)$species))
  
  n_minus_160_similarity <- (length(subject_cell$species) + length(n_minus_160$species) - length(unique(rbind(subject_cell,n_minus_160)$species)))/length(unique(rbind(subject_cell,n_minus_160)$species))
  
  n_minus_161_similarity <- (length(subject_cell$species) + length(n_minus_161$species) - length(unique(rbind(subject_cell,n_minus_161)$species)))/length(unique(rbind(subject_cell,n_minus_161)$species))
  
  n_minus_162_similarity <- (length(subject_cell$species) + length(n_minus_162$species) - length(unique(rbind(subject_cell,n_minus_162)$species)))/length(unique(rbind(subject_cell,n_minus_162)$species))
  
  for (j in 2003:2019) {
    rf_biotic_data$bs_ring_2[(rf_biotic_data$year == j) & 
                               (rf_biotic_data$grid_cell == texas_cells[i])] <- 
      (if((texas_cells[i] + 2) %in% texas_cells) {
        (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                     (establishment_tracker$grid_cell == (texas_cells[i]+2))] * n_plus_2_similarity)
      } else {
        0
      }) + 
      (if((texas_cells[i] - 2) %in% texas_cells) {
        (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                     (establishment_tracker$grid_cell == (texas_cells[i]-2))] * n_minus_2_similarity)
      } else {
        0
      }) +
      (if((texas_cells[i] + 79) %in% texas_cells) {
        (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                     (establishment_tracker$grid_cell == (texas_cells[i]+79))] * n_plus_79_similarity)
      } else {
        0
      }) +
      (if((texas_cells[i] + 82) %in% texas_cells) {
        (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                     (establishment_tracker$grid_cell == (texas_cells[i]+82))] * n_plus_82_similarity)
      } else {
        0
      }) +
      (if((texas_cells[i] - 79) %in% texas_cells) {
        (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                     (establishment_tracker$grid_cell == (texas_cells[i]-79))] * n_minus_79_similarity)
      } else {
        0
      }) +
      (if((texas_cells[i] - 82) %in% texas_cells) {
        (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                     (establishment_tracker$grid_cell == (texas_cells[i]-82))] * n_minus_82_similarity)
      } else {
        0
        }) +
      (if((texas_cells[i] + 160) %in% texas_cells) {
        (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                     (establishment_tracker$grid_cell == (texas_cells[i]+160))] * n_plus_160_similarity)
      } else {
        0
      }) + 
      (if((texas_cells[i] + 161) %in% texas_cells) {
        (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                     (establishment_tracker$grid_cell == (texas_cells[i]+161))] * n_plus_161_similarity)
      } else {
        0
      }) +
      (if((texas_cells[i] + 162) %in% texas_cells) {
        (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                     (establishment_tracker$grid_cell == (texas_cells[i]+162))] * n_plus_162_similarity)
      } else {
        0
      }) +
      (if((texas_cells[i] - 160) %in% texas_cells) {
        (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                     (establishment_tracker$grid_cell == (texas_cells[i]-160))] * n_minus_160_similarity)
      } else {
        0
      }) +
      (if((texas_cells[i] - 161) %in% texas_cells) {
        (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                     (establishment_tracker$grid_cell == (texas_cells[i]-161))] * n_minus_161_similarity)
      } else {
        0
      }) +
      (if((texas_cells[i] - 162) %in% texas_cells) {
        (establishment_tracker$establishment_score[(establishment_tracker$year == (j-1)) & 
                                                     (establishment_tracker$grid_cell == (texas_cells[i]-162))] * n_minus_162_similarity)
      } else {
        0
      })
  }
}

# don't need to worry about na's because they'll get factored out anyway

# save rf biotic inputs for Texas to a csv file
write.csv(rf_biotic_data, "data/texas_rf_biotic_data.csv")


