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

establishment_tracker <- read.csv("data/FINAL-VERSION-FL_egoose_establishment_ct-50_met-0.005.csv")

lc_vectors <- read.csv("data/lc_feature_vectors.csv")

prcp_vector <- read.csv("data/prcp_feature_vector.csv")

tmax_vector <- read.csv("data/tmax_feature_vector.csv")

# Find centre-point of each cell and work out how far cells are apart ----

dggs <- dgconstruct(res = 8, projection = "ISEA", metric = TRUE, 
                    resround = 'nearest')

cellcenters <- dgSEQNUM_to_GEO(dggs,establishment_tracker$grid_cell)

cellcenters_df <- as.data.frame(cellcenters)

establishment_tracker <- cbind(establishment_tracker, cellcenters_df)

# calculate haversine distance for all cells from 33453 ----

establishment_tracker$haversine <- 0

for (i in 1:nrow(establishment_tracker)) {
  establishment_tracker$haversine[i] <-distHaversine(c(establishment_tracker$lon_deg[i],
                                                       establishment_tracker$lat_deg[i]),
                                                     c(-80.58184,25.63718))
}

ggplot(establishment_tracker, aes(x=haversine)) + geom_histogram(bins = 50)


establishment_tracker_2000 <- filter(establishment_tracker, 
                                     establishment_tracker$year == 2000)

# should have 6 which are about the same, then 12 which are about the same

et_2000_arranged <- establishment_tracker_2000 %>%
  arrange(haversine)

ggplot(establishment_tracker, aes(x = lon_deg, y = lat_deg)) + 
  geom_point() +
  scale_y_continuous( breaks = seq(23, 32, by = 0.5)) +
  scale_x_continuous( breaks = seq(-88, -79, by = 1))

# 1st ring cell values for central cell = n: n+1, n-1, n+81, n+82, n-81, n-82

# 2nd ring cell values for central cell = n: n+2, n-2, n+80, n+83, n-83, n-80,
# n-162, n-163, n-164, n+162, n+163, n+164

## Create new dataframe with feature vectors

rf_abiotic_data <- filter(establishment_tracker, establishment_tracker$year > 2002)[,c(1,2,3,6)]

rf_abiotic_data$tmax_vector <- 0

rf_abiotic_data$prcp_vector <- 0

rf_abiotic_data$water_cover <- 0

rf_abiotic_data$forest_cover <- 0

rf_abiotic_data$grass_cover <- 0

rf_abiotic_data$wetland_cover <- 0

rf_abiotic_data$farming_cover <- 0

rf_abiotic_data$urban_cover <- 0

rf_abiotic_data$barren_cover <- 0

rf_abiotic_data$previous_establishment <- 0

#1st ring tmax vector calc

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$tmax_vector[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                       rf_abiotic_data$year[i] > 2003 &
                                       (rf_abiotic_data$grid_cell[i]+1) %in% rf_abiotic_data$grid_cell &
                                       (rf_abiotic_data$grid_cell[i]+1) %in% tmax_vector$cell &
                                       rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
      (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                             (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+1))]) / 
        exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]+1)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
    } else {
      0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
         rf_abiotic_data$year[i] > 2003 &
         (rf_abiotic_data$grid_cell[i]-1) %in% rf_abiotic_data$grid_cell &
         (rf_abiotic_data$grid_cell[i]-1) %in% tmax_vector$cell &
         rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-1))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]-1)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
         rf_abiotic_data$year[i] > 2003 &
         (rf_abiotic_data$grid_cell[i]-81) %in% rf_abiotic_data$grid_cell &
         (rf_abiotic_data$grid_cell[i]-81) %in% tmax_vector$cell &
         rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-81))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]-81)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
       rf_abiotic_data$year[i] > 2003 &
       (rf_abiotic_data$grid_cell[i]-82) %in% rf_abiotic_data$grid_cell &
       (rf_abiotic_data$grid_cell[i]-82) %in% tmax_vector$cell &
       rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-82))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]-82)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
         rf_abiotic_data$year[i] > 2003 &
         (rf_abiotic_data$grid_cell[i]+81) %in% rf_abiotic_data$grid_cell &
         (rf_abiotic_data$grid_cell[i]+81) %in% tmax_vector$cell &
         rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+81))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]+81)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
         rf_abiotic_data$year[i] > 2003 &
         (rf_abiotic_data$grid_cell[i]+82) %in% rf_abiotic_data$grid_cell &
         (rf_abiotic_data$grid_cell[i]+82) %in% tmax_vector$cell &
         rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+82))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]+82)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  })
}

#1st ring precipitation vector calc

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$prcp_vector[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                        rf_abiotic_data$year[i] > 2003 &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% rf_abiotic_data$grid_cell &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% prcp_vector$cell &
                                        rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+1))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]+1)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-1) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-1) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-1))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]-1)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-81) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-81))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]-81)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-82) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-82))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]-82)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+81) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+81))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]+81)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+82) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+82))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]+82)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  })
}

#1st ring water cover vector calc

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$water_cover[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                        rf_abiotic_data$year[i] > 2003 &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% rf_abiotic_data$grid_cell &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% lc_vectors$cell &
                                        rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+1))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+1)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-1) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-1) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-1))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-1)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-81) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-81))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-81)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-82) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-82))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-82)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+81) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+81))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+81)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+82) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+82))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+82)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  })
}

#1st ring forest cover vector calc

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$forest_cover[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                        rf_abiotic_data$year[i] > 2003 &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% rf_abiotic_data$grid_cell &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% lc_vectors$cell &
                                        rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+1))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+1)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-1) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-1) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-1))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-1)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-81) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-81))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-81)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-82) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-82))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-82)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+81) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+81))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+81)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+82) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+82))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+82)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  })
}

#1st ring grass cover vector calc

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$grass_cover[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                        rf_abiotic_data$year[i] > 2003 &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% rf_abiotic_data$grid_cell &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% lc_vectors$cell &
                                        rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+1))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+1)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-1) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-1) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-1))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-1)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-81) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-81))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-81)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-82) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-82))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-82)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+81) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+81))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+81)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+82) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+82))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+82)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  })
}

#1st ring wetland cover vector calc

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$wetland_cover[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                        rf_abiotic_data$year[i] > 2003 &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% rf_abiotic_data$grid_cell &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% lc_vectors$cell &
                                        rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+1))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+1)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-1) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-1) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-1))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-1)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-81) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-81))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-81)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-82) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-82))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-82)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+81) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+81))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+81)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+82) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+82))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+82)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  })
}

#1st ring farming cover vector calc

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$farming_cover[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                        rf_abiotic_data$year[i] > 2003 &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% rf_abiotic_data$grid_cell &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% lc_vectors$cell &
                                        rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+1))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+1)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-1) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-1) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-1))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-1)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-81) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-81))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-81)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-82) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-82))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-82)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+81) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+81))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+81)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+82) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+82))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+82)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  })
}

#1st ring urban cover vector calc

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$urban_cover[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                        rf_abiotic_data$year[i] > 2003 &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% rf_abiotic_data$grid_cell &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% lc_vectors$cell &
                                        rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+1))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+1)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-1) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-1) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-1))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-1)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-81) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-81))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-81)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-82) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-82))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-82)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+81) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+81))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+81)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+82) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+82))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+82)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  })
}

#1st ring barren cover vector calc

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$barren_cover[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                        rf_abiotic_data$year[i] > 2003 &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% rf_abiotic_data$grid_cell &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% lc_vectors$cell &
                                        rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+1))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+1)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-1) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-1) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-1))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-1)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-81) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-81))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-81)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-82) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-82))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-82)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+81) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+81))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+81)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+82) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+82))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+82)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  })
}



# previous establishment value

for (i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$previous_establishment[i] <- if (rf_abiotic_data$enough_checklists[i] == 1 &
                                                   rf_abiotic_data$year[i] > 2003) {
    rf_abiotic_data$establishment_score[which((rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]) & 
                                                                                           (rf_abiotic_data$year == (rf_abiotic_data$year[i] - 1)))]
  } else {
    0
  }
}

## Start with 2nd ring ----

colnames(rf_abiotic_data)[5:13] <- c("tmax_ring_1", "prcp_ring_1", "water_ring_1",
                                    "forest_ring_1", "grass_ring_1", 
                                    "wetland_ring_1", "farming_ring_1", 
                                    "urban_ring_1", "barren_ring_1")

rf_abiotic_data$tmax_ring_2 <- 0

rf_abiotic_data$prcp_ring_2 <- 0

rf_abiotic_data$water_ring_2 <- 0

rf_abiotic_data$forest_ring_2 <- 0

rf_abiotic_data$grass_ring_2 <- 0

rf_abiotic_data$wetland_ring_2 <- 0

rf_abiotic_data$farming_ring_2 <- 0

rf_abiotic_data$urban_ring_2 <- 0

rf_abiotic_data$barren_ring_2 <- 0

# create tmax_ring_2 variable

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$tmax_ring_2[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                        rf_abiotic_data$year[i] > 2003 &
                                        (rf_abiotic_data$grid_cell[i]+2) %in% rf_abiotic_data$grid_cell &
                                        (rf_abiotic_data$grid_cell[i]+2) %in% tmax_vector$cell &
                                        rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+2)))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]+2)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-2) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-2) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-2)))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]-2)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-83) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-83) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-83)))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]-83)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-80) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-80)))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]-80)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+80) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+80)))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]+80)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+83) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+83) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+83))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]+83)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-162) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-162))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]-162)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-163) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-163) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-163))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]-163)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-164) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-164) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-164))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]-164)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+162) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+162))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]+162)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+163) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+163) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+163))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]+163)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+164) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+164) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+164))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]+164)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) 
}

# repeat for 2nd ring of prcp

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$prcp_ring_2[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                        rf_abiotic_data$year[i] > 2003 &
                                        (rf_abiotic_data$grid_cell[i]+2) %in% rf_abiotic_data$grid_cell &
                                        (rf_abiotic_data$grid_cell[i]+2) %in% prcp_vector$cell &
                                        rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+2)))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]+2)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-2) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-2) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-2)))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]-2)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-83) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-83) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-83)))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]-83)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-80) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-80)))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]-80)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+80) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+80)))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]+80)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+83) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+83) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+83))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]+83)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-162) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-162))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]-162)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-163) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-163) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-163))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]-163)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-164) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-164) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-164))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]-164)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+162) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+162))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]+162)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+163) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+163) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+163))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]+163)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+164) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+164) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+164))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]+164)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) 
}

# water_ring_2

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$water_ring_2[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                        rf_abiotic_data$year[i] > 2003 &
                                        (rf_abiotic_data$grid_cell[i]+2) %in% rf_abiotic_data$grid_cell &
                                        (rf_abiotic_data$grid_cell[i]+2) %in% lc_vectors$cell &
                                        rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+2)))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+2)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-2) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-2) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-2)))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-2)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-83) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-83) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-83)))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-83)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-80) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-80)))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-80)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+80) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+80)))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+80)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+83) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+83) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+83))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+83)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-162) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-162))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-162)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-163) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-163) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-163))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-163)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-164) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-164) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-164))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-164)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+162) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+162))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+162)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+163) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+163) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+163))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+163)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+164) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+164) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+164))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+164)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) 
}

# forest_ring_2

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$forest_ring_2[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                         rf_abiotic_data$year[i] > 2003 &
                                         (rf_abiotic_data$grid_cell[i]+2) %in% rf_abiotic_data$grid_cell &
                                         (rf_abiotic_data$grid_cell[i]+2) %in% lc_vectors$cell &
                                         rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+2)))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+2)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-2) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-2) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-2)))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-2)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-83) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-83) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-83)))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-83)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-80) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-80)))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-80)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+80) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+80)))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+80)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+83) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+83) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+83))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+83)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-162) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-162))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-162)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-163) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-163) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-163))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-163)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-164) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-164) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-164))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-164)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+162) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+162))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+162)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+163) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+163) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+163))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+163)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+164) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+164) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+164))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+164)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) 
}

# grass_ring_2

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$grass_ring_2[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                         rf_abiotic_data$year[i] > 2003 &
                                         (rf_abiotic_data$grid_cell[i]+2) %in% rf_abiotic_data$grid_cell &
                                         (rf_abiotic_data$grid_cell[i]+2) %in% lc_vectors$cell &
                                         rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+2)))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+2)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-2) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-2) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-2)))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-2)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-83) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-83) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-83)))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-83)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-80) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-80)))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-80)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+80) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+80)))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+80)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+83) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+83) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+83))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+83)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-162) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-162))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-162)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-163) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-163) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-163))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-163)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-164) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-164) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-164))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-164)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+162) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+162))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+162)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+163) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+163) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+163))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+163)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+164) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+164) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+164))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+164)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) 
}

# wetland_ring_2

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$wetland_ring_2[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                         rf_abiotic_data$year[i] > 2003 &
                                         (rf_abiotic_data$grid_cell[i]+2) %in% rf_abiotic_data$grid_cell &
                                         (rf_abiotic_data$grid_cell[i]+2) %in% lc_vectors$cell &
                                         rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+2)))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+2)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-2) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-2) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-2)))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-2)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-83) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-83) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-83)))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-83)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-80) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-80)))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-80)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+80) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+80)))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+80)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+83) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+83) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+83))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+83)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-162) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-162))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-162)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-163) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-163) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-163))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-163)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-164) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-164) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-164))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-164)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+162) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+162))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+162)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+163) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+163) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+163))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+163)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+164) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+164) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+164))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+164)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) 
}

# farming_ring_2

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$farming_ring_2[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                         rf_abiotic_data$year[i] > 2003 &
                                         (rf_abiotic_data$grid_cell[i]+2) %in% rf_abiotic_data$grid_cell &
                                         (rf_abiotic_data$grid_cell[i]+2) %in% lc_vectors$cell &
                                         rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+2)))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+2)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-2) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-2) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-2)))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-2)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-83) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-83) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-83)))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-83)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-80) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-80)))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-80)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+80) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+80)))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+80)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+83) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+83) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+83))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+83)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-162) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-162))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-162)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-163) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-163) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-163))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-163)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-164) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-164) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-164))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-164)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+162) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+162))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+162)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+163) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+163) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+163))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+163)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+164) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+164) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+164))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+164)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) 
}

# urban_ring_2

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$urban_ring_2[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                         rf_abiotic_data$year[i] > 2003 &
                                         (rf_abiotic_data$grid_cell[i]+2) %in% rf_abiotic_data$grid_cell &
                                         (rf_abiotic_data$grid_cell[i]+2) %in% lc_vectors$cell &
                                         rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+2)))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+2)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-2) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-2) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-2)))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-2)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-83) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-83) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-83)))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-83)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-80) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-80)))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-80)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+80) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+80)))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+80)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+83) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+83) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+83))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+83)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-162) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-162))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-162)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-163) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-163) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-163))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-163)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-164) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-164) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-164))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-164)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+162) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+162))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+162)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+163) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+163) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+163))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+163)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+164) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+164) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+164))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+164)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) 
}

# barren_ring_2

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$barren_ring_2[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                         rf_abiotic_data$year[i] > 2003 &
                                         (rf_abiotic_data$grid_cell[i]+2) %in% rf_abiotic_data$grid_cell &
                                         (rf_abiotic_data$grid_cell[i]+2) %in% lc_vectors$cell &
                                         rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+2)))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+2)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-2) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-2) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-2)))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-2)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-83) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-83) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-83)))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-83)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-80) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-80)))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-80)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+80) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+80)))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+80)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+83) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+83) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+83))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+83)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-162) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-162))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-162)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-163) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-163) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-163))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-163)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]-164) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-164) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-164))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-164)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+162) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+162))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+162)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+163) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+163) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+163))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+163)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2003 &
           (rf_abiotic_data$grid_cell[i]+164) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+164) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+164))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+164)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) 
}

# now save data

write.csv(rf_abiotic_data, "data/11.08.23-rf_model_abiotic_data.csv")
