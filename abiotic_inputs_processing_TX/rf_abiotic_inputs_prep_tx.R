## PREPARING DATA FOR ABIOTIC RANDOM FOREST INVASIVE SPECIES ESTABLISHMNET
## ALGORITHM

# Author = John Gray
# Email = johnpatrickgray97@gmail.com
# Last edit = 10/07/23

# Packages ----

library(dggridR)
library(geosphere)
library(ggplot2)
library(dplyr)

# Loading in data ----

establishment_tracker <- read.csv("data/TX_egoose_establishment_ct-50_met-0.005.csv")

lc_vectors <- read.csv("data/texas_lc_feature_vectors.csv")

prcp_vector <- read.csv("data/prcp_feature_vector.csv")

tmax_vector <- read.csv("data/tmax_feature_vector.csv")

# Find centre-point of each cell and work out how far cells are apart ----

dggs <- dgconstruct(res = 8, projection = "ISEA", metric = TRUE, 
                    resround = 'nearest')

cellcenters <- dgSEQNUM_to_GEO(dggs,establishment_tracker$grid_cell)

cellcenters_df <- as.data.frame(cellcenters)

establishment_tracker <- cbind(establishment_tracker, cellcenters_df)

# 1st ring cell values for central cell = n: n+1, n-1, n+80, n+81, n-80, n-81

# 2nd ring cell values for central cell = n: n+2, n-2, n+80, n+83, n-83, n-80,
# n-162, n-163, n-164, n+162, n+160, n+164

## Create new dataframe with feature vectors

rf_abiotic_data <- filter(establishment_tracker, establishment_tracker$year > 2009)[,c(1,2,3,6)]

rf_abiotic_data$tmax_vector <- 0

rf_abiotic_data$prcp_vector <- 0

# maybe change land cover vector in the future to consider all components separately

rf_abiotic_data$land_cover <- 0

rf_abiotic_data$previous_establishment <- 0

#1st ring tmax vector calc

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$tmax_vector[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                       rf_abiotic_data$year[i] > 2010 &
                                       (rf_abiotic_data$grid_cell[i]+1) %in% rf_abiotic_data$grid_cell &
                                       (rf_abiotic_data$grid_cell[i]+1) %in% tmax_vector$cell &
                                       rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
      (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                             (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+1))]) / 
        abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]+1)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])]))
    } else {
      0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
         rf_abiotic_data$year[i] > 2010 &
         (rf_abiotic_data$grid_cell[i]-1) %in% rf_abiotic_data$grid_cell &
         (rf_abiotic_data$grid_cell[i]-1) %in% tmax_vector$cell &
         rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-1))]) / 
      abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]-1)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
         rf_abiotic_data$year[i] > 2010 &
         (rf_abiotic_data$grid_cell[i]-81) %in% rf_abiotic_data$grid_cell &
         (rf_abiotic_data$grid_cell[i]-81) %in% tmax_vector$cell &
         rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-81))]) / 
      abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]-81)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
       rf_abiotic_data$year[i] > 2010 &
       (rf_abiotic_data$grid_cell[i]-80) %in% rf_abiotic_data$grid_cell &
       (rf_abiotic_data$grid_cell[i]-80) %in% tmax_vector$cell &
       rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-80))]) / 
      abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]-80)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
         rf_abiotic_data$year[i] > 2010 &
         (rf_abiotic_data$grid_cell[i]+81) %in% rf_abiotic_data$grid_cell &
         (rf_abiotic_data$grid_cell[i]+81) %in% tmax_vector$cell &
         rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+81))]) / 
      abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]+81)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
         rf_abiotic_data$year[i] > 2010 &
         (rf_abiotic_data$grid_cell[i]+80) %in% rf_abiotic_data$grid_cell &
         (rf_abiotic_data$grid_cell[i]+80) %in% tmax_vector$cell &
         rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+80))]) / 
      abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]+80)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  })
}

#1st ring precipitation vector calc

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$prcp_vector[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                        rf_abiotic_data$year[i] > 2010 &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% rf_abiotic_data$grid_cell &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% prcp_vector$cell &
                                        rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+1))]) / 
      abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]+1)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-1) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-1) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-1))]) / 
      abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]-1)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-81) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-81))]) / 
      abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]-81)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-80) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-80))]) / 
      abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]-80)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+81) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+81))]) / 
      abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]+81)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+80) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+80))]) / 
      abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]+80)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  })
}

# 1st ring lc_type vector calc

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$land_cover[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                        rf_abiotic_data$year[i] > 2010 &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% rf_abiotic_data$grid_cell &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% lc_vectors$cell &
                                        rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+1))]) / 
      (abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+1)]) - 
             (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
        abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+1)]) - 
              (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
        abs((lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+1)]) - 
              (lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
        abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+1)]) - 
              (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
        abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+1)]) - 
              (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
        abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+1)]) - 
              (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
        abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+1)]) - 
              (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-1) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-1) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-1))]) / 
      (abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-1)]) - 
             (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-1)]) - 
               (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-1)]) - 
               (lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-1)]) - 
               (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-1)]) - 
               (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-1)]) - 
               (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-1)]) - 
               (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-81) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-81))]) / 
      (abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-81)]) - 
             (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-81)]) - 
               (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-81)]) - 
               (lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-81)]) - 
               (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-81)]) - 
               (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-81)]) - 
               (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-81)]) - 
               (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-80) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-80))]) / 
      (abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-80)]) - 
             (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-80)]) - 
               (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-80)]) - 
               (lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-80)]) - 
               (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-80)]) - 
               (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-80)]) - 
               (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-80)]) - 
               (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+81) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+81))]) / 
      (abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+81)]) - 
             (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+81)]) - 
               (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+81)]) - 
               (lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+81)]) - 
               (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+81)]) - 
               (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+81)]) - 
               (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+81)]) - 
               (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+80) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+80))]) / 
      (abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+80)]) - 
             (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+80)]) - 
               (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+80)]) - 
               (lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+80)]) - 
               (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+80)]) - 
               (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+80)]) - 
               (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+80)]) - 
               (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  })
}

# previous establishment value

for (i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$previous_establishment[i] <- if (rf_abiotic_data$enough_checklists[i] == 1 &
                                                   rf_abiotic_data$year[i] > 2010) {
    rf_abiotic_data$establishment_score[which((rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]) & 
                                                                                           (rf_abiotic_data$year == (rf_abiotic_data$year[i] - 1)))]
  } else {
    0
  }
}

colnames(rf_abiotic_data)[5:7] <- c("tmax_ring_1", "prcp_ring_1", "lc_ring_1")

rf_abiotic_data$tmax_ring_2 <- 0

rf_abiotic_data$prcp_ring_2 <- 0

rf_abiotic_data$lc_ring_2 <- 0

# create tmax_ring_2 variable

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$tmax_ring_2[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                        rf_abiotic_data$year[i] > 2010 &
                                        (rf_abiotic_data$grid_cell[i]+2) %in% rf_abiotic_data$grid_cell &
                                        (rf_abiotic_data$grid_cell[i]+2) %in% tmax_vector$cell &
                                        rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+2)))]) / 
      abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]+2)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-2) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-2) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-2)))]) / 
      abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]-2)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-82) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-82)))]) / 
      abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]-82)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-79) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-79) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-79)))]) / 
      abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]-79)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+79) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+79) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+79)))]) / 
      abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]+79)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+82) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+82))]) / 
      abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]+82)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-162) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-162))]) / 
      abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]-162)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-160) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-160) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-160))]) / 
      abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]-160)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-161) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-161) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-161))]) / 
      abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]-161)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+162) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+162))]) / 
      abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]+162)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+160) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+160) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+160))]) / 
      abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]+160)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+161) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+161) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+161))]) / 
      abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]+161)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) 
}

# repeat for 2nd ring of prcp

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$prcp_ring_2[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                        rf_abiotic_data$year[i] > 2010 &
                                        (rf_abiotic_data$grid_cell[i]+2) %in% rf_abiotic_data$grid_cell &
                                        (rf_abiotic_data$grid_cell[i]+2) %in% prcp_vector$cell &
                                        rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+2)))]) / 
      abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]+2)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-2) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-2) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-2)))]) / 
      abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]-2)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-82) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-82)))]) / 
      abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]-82)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-79) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-79) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-79)))]) / 
      abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]-79)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+79) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+79) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+79)))]) / 
      abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]+79)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+82) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+82))]) / 
      abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]+82)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-162) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-162))]) / 
      abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]-162)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-160) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-160) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-160))]) / 
      abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]-160)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-161) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-161) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-161))]) / 
      abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]-161)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+162) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+162))]) / 
      abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]+162)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+160) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+160) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+160))]) / 
      abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]+160)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+161) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+161) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+161))]) / 
      abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]+161)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])]))
  } else {
    0
  }) 
}

# now for lc ring 2 features 

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$lc_ring_2[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                        rf_abiotic_data$year[i] > 2010 &
                                        (rf_abiotic_data$grid_cell[i]+2) %in% rf_abiotic_data$grid_cell &
                                        (rf_abiotic_data$grid_cell[i]+2) %in% lc_vectors$cell &
                                        rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+2)))]) / 
      (abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+2)]) - 
             (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+2)]) - 
               (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+2)]) - 
               (lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+2)]) - 
               (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+2)]) - 
               (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+2)]) - 
               (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+2)]) - 
               (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-2) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-2) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-2)))]) / 
      (abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-2)]) - 
             (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-2)]) - 
               (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-2)]) - 
               (lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-2)]) - 
               (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-2)]) - 
               (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-2)]) - 
               (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-2)]) - 
               (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-82) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-82)))]) / 
      (abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-82)]) - 
             (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-82)]) - 
               (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-82)]) - 
               (lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-82)]) - 
               (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-82)]) - 
               (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-82)]) - 
               (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+-82)]) - 
               (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-79) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-79) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-79)))]) / 
      (abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-79)]) - 
             (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-79)]) - 
               (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-79)]) - 
               (lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-79)]) - 
               (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-79)]) - 
               (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-79)]) - 
               (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-79)]) - 
               (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+79) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+79) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+79)))]) / 
      (abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+79)]) - 
             (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+79)]) - 
               (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+79)]) - 
               (lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+79)]) - 
               (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+79)]) - 
               (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+79)]) - 
               (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+79)]) - 
               (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+82) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+82))]) / 
      (abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+82)]) - 
             (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+82)]) - 
               (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+82)]) - 
               (lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+82)]) - 
               (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+82)]) - 
               (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+82)]) - 
               (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+82)]) - 
               (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-162) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-162))]) / 
      (abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-162)]) - 
             (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-162)]) - 
               (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-162)]) - 
               (lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-162)]) - 
               (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-162)]) - 
               (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-162)]) - 
               (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-162)]) - 
               (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-160) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-160) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-160))]) / 
      (abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-160)]) - 
             (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-160)]) - 
               (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-160)]) - 
               (lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-160)]) - 
               (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-160)]) - 
               (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-160)]) - 
               (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-160)]) - 
               (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-161) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-161) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-161))]) / 
      (abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-161)]) - 
             (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-161)]) - 
               (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-161)]) - 
               (lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-161)]) - 
               (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-161)]) - 
               (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-161)]) - 
               (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-161)]) - 
               (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+162) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+162))]) / 
      (abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+162)]) - 
             (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+162)]) - 
               (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+162)]) - 
               (lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+162)]) - 
               (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+162)]) - 
               (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+162)]) - 
               (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+162)]) - 
               (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+160) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+160) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+160))]) / 
      (abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+160)]) - 
             (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+160)]) - 
               (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+160)]) - 
               (lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+160)]) - 
               (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+160)]) - 
               (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+160)]) - 
               (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+160)]) - 
               (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+161) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+161) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+161))]) / 
      (abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+161)]) - 
             (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+161)]) - 
               (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+161)]) - 
               (lc_vectors$mean_grass_etc[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+161)]) - 
               (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+161)]) - 
               (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) +
         abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+161)]) - 
               (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])) + 
         abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+161)]) - 
               (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) 
}

# now save data

write.csv(rf_abiotic_data, "data/rf_model_abiotic_data_texas.csv")
