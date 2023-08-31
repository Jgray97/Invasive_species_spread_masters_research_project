## PREPARING TEXAS ABIOTIC INPUT VARIABLES FOR RF INVASIVE SPECIES 
## ESTABLISHMENT ALGORITHM 

# Author = John Gray
# Email = greyjohn15@gmail.com
# Last edit = 30/08/23

# Load in packages ----

library(dggridR)
library(geosphere)
library(ggplot2)
library(dplyr)

# Loading in data ----

# Establishment tracker data
establishment_tracker <- read.csv("data/TX_egoose_establishment_ct-50_met-0.005.csv")

# Land cover type data
lc_vectors <- read.csv("data/texas_lc_feature_vectors.csv")

# precipitation data
prcp_vector <- read.csv("data/texas_prcp_feature_vector.csv")

# temperature data
tmax_vector <- read.csv("data/texas_tmax_feature_vector.csv")

# Find centre-point of each cell and work out how far cells are apart ----

# create hexagon grid
dggs <- dgconstruct(res = 8, projection = "ISEA", metric = TRUE, 
                    resround = 'nearest')

# determine cell centre points
cellcenters <- dgSEQNUM_to_GEO(dggs,establishment_tracker$grid_cell)

# turn cell centre points into df
cellcenters_df <- as.data.frame(cellcenters)

# add centre points to establishment tracker df
establishment_tracker <- cbind(establishment_tracker, cellcenters_df)

# AFTER SOME MANUAL INVESTIGATION I HAVE FOUND THAT

# 1st ring cell values for central cell = n: n+1, n-1, n+80, n+81, n-80, n-81

# 2nd ring cell values for central cell = n: n+2, n-2, n+79, n+82, n-82, n-79,
# n-160, n-161, n-162, n+160, n+161, n+162

## Create new dataframe with establishment tracker + 1st ring abiotic rf ----
## abiotic vectors ----

rf_abiotic_data <- filter(establishment_tracker, establishment_tracker$year > 2009)[,c(1,2,3,6)]

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

### For loop that calculates 1st ring temperature related rf input values
# Value is defined by the sum of establishment score in previous years / 
# exponent of difference in temperature compared with adjacent cell for every 
# cell in 1st ring

# We ignore cells where there aren't enough checklists and which don't belong to a
# year in which there is any establishment and only calculate dissimilarity if
# both target and adjacent cell are present in the rf_abiotic data dataframe

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$tmax_vector[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                       rf_abiotic_data$year[i] > 2010 &
                                       (rf_abiotic_data$grid_cell[i]+1) %in% rf_abiotic_data$grid_cell &
                                       (rf_abiotic_data$grid_cell[i]+1) %in% tmax_vector$cell &
                                       rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
      (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                             (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+1))]) / 
        exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]+1)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
    } else {
      0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
         rf_abiotic_data$year[i] > 2010 &
         (rf_abiotic_data$grid_cell[i]-1) %in% rf_abiotic_data$grid_cell &
         (rf_abiotic_data$grid_cell[i]-1) %in% tmax_vector$cell &
         rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-1))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]-1)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
         rf_abiotic_data$year[i] > 2010 &
         (rf_abiotic_data$grid_cell[i]-81) %in% rf_abiotic_data$grid_cell &
         (rf_abiotic_data$grid_cell[i]-81) %in% tmax_vector$cell &
         rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-81))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]-81)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
       rf_abiotic_data$year[i] > 2010 &
       (rf_abiotic_data$grid_cell[i]-80) %in% rf_abiotic_data$grid_cell &
       (rf_abiotic_data$grid_cell[i]-80) %in% tmax_vector$cell &
       rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-80))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]-80)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
         rf_abiotic_data$year[i] > 2010 &
         (rf_abiotic_data$grid_cell[i]+81) %in% rf_abiotic_data$grid_cell &
         (rf_abiotic_data$grid_cell[i]+81) %in% tmax_vector$cell &
         rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+81))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]+81)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
         rf_abiotic_data$year[i] > 2010 &
         (rf_abiotic_data$grid_cell[i]+80) %in% rf_abiotic_data$grid_cell &
         (rf_abiotic_data$grid_cell[i]+80) %in% tmax_vector$cell &
         rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+80))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]+80)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  })
}

### For loop that calculates 1st ring precipitation related rf input values
# Defined same as for temperature but with precipitation instead
# Same cleaning/checks/assumptions made

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$prcp_vector[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                        rf_abiotic_data$year[i] > 2010 &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% rf_abiotic_data$grid_cell &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% prcp_vector$cell &
                                        rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+1))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]+1)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-1) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-1) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-1))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]-1)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-81) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-81))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]-81)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-80) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-80))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]-80)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+81) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+81))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]+81)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+80) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+80))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]+80)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  })
}

# Land cover input variables, starting with water cover are all calculated in 
# the same way as temperature and precipitation, with the same assumptions/checks

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$water_cover[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                        rf_abiotic_data$year[i] > 2010 &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% rf_abiotic_data$grid_cell &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% lc_vectors$cell &
                                        rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+1))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+1)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-1) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-1) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-1))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-1)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-81) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-81))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-81)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-80) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-80))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-80)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+81) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+81))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+81)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+80) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+80))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+80)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  })
}

#1st ring forest cover input variable calc

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$forest_cover[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                        rf_abiotic_data$year[i] > 2010 &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% rf_abiotic_data$grid_cell &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% lc_vectors$cell &
                                        rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+1))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+1)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-1) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-1) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-1))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-1)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-81) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-81))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-81)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-80) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-80))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-80)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+81) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+81))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+81)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+80) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+80))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+80)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  })
}

#1st ring grass cover input variable calc

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$grass_cover[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                        rf_abiotic_data$year[i] > 2010 &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% rf_abiotic_data$grid_cell &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% lc_vectors$cell &
                                        rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+1))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+1)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-1) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-1) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-1))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-1)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-81) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-81))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-81)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-80) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-80))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-80)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+81) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+81))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+81)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+80) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+80))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+80)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  })
}

#1st ring wetland cover input variable calc

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$wetland_cover[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                        rf_abiotic_data$year[i] > 2010 &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% rf_abiotic_data$grid_cell &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% lc_vectors$cell &
                                        rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+1))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+1)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-1) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-1) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-1))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-1)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-81) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-81))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-81)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-80) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-80))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-80)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+81) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+81))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+81)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+80) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+80))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+80)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  })
}

#1st ring farming cover input variable calc

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$farming_cover[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                        rf_abiotic_data$year[i] > 2010 &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% rf_abiotic_data$grid_cell &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% lc_vectors$cell &
                                        rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+1))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+1)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-1) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-1) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-1))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-1)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-81) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-81))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-81)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-80) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-80))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-80)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+81) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+81))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+81)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+80) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+80))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+80)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  })
}

#1st ring urban cover input variable calc

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$urban_cover[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                        rf_abiotic_data$year[i] > 2010 &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% rf_abiotic_data$grid_cell &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% lc_vectors$cell &
                                        rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+1))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+1)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-1) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-1) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-1))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-1)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-81) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-81))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-81)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-80) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-80))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-80)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+81) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+81))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+81)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+80) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+80))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+80)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  })
}

#1st ring barren cover input variable calc

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$barren_cover[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                        rf_abiotic_data$year[i] > 2010 &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% rf_abiotic_data$grid_cell &
                                        (rf_abiotic_data$grid_cell[i]+1) %in% lc_vectors$cell &
                                        rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+1))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+1)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-1) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-1) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-1))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-1)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-81) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-81))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-81)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-80) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-80))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-80)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+81) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+81) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+81))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+81)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+80) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+80) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+80))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+80)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  })
}

# Add in an rf input variable for the establishment score in the target cell in
# the previous year

for (i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$previous_establishment[i] <- if (rf_abiotic_data$enough_checklists[i] == 1 &
                                                   rf_abiotic_data$year[i] > 2010) {
    rf_abiotic_data$establishment_score[which((rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]) & 
                                                                                           (rf_abiotic_data$year == (rf_abiotic_data$year[i] - 1)))]
  } else {
    0
  }
}

## Start with 2nd ring ----

# First off, rename columns for all 1st ring rf abiotic input variables

colnames(rf_abiotic_data)[5:13] <- c("tmax_ring_1", "prcp_ring_1", "water_ring_1",
                                    "forest_ring_1", "grass_ring_1", 
                                    "wetland_ring_1", "farming_ring_1",
                                    "urban_ring_1", "barren_ring_1")

# now create new columns for 2nd ring abiotic input variables

rf_abiotic_data$tmax_ring_2 <- 0

rf_abiotic_data$prcp_ring_2 <- 0

rf_abiotic_data$water_ring_2 <- 0

rf_abiotic_data$forest_ring_2 <- 0

rf_abiotic_data$grass_ring_2 <- 0

rf_abiotic_data$wetland_ring_2 <- 0

rf_abiotic_data$farming_ring_2 <- 0

rf_abiotic_data$urban_ring_2 <- 0

rf_abiotic_data$barren_ring_2 <- 0

# Starting with temperature related input, the 2nd ring input variables are 
# calculated in exactly the same way, except that the cells included in the sum 
# are those which are separated from the target cell by a single cell

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$tmax_ring_2[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                        rf_abiotic_data$year[i] > 2010 &
                                        (rf_abiotic_data$grid_cell[i]+2) %in% rf_abiotic_data$grid_cell &
                                        (rf_abiotic_data$grid_cell[i]+2) %in% tmax_vector$cell &
                                        rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+2)))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]+2)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-2) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-2) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-2)))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]-2)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-82) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-82)))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]-82)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-79) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-79) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-79)))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]-79)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+79) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+79) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+79)))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]+79)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+82) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+82))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]+82)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-162) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-162))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]-162)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-160) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-160) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-160))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]-160)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-161) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-161) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-161))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]-161)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+162) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+162))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]+162)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+160) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+160) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+160))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]+160)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+161) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+161) %in% tmax_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% tmax_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+161))]) / 
      exp(abs((tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i]+161)]) - (tmax_vector$mean_tmax[which(tmax_vector$cell == rf_abiotic_data$grid_cell[i])])))
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
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]+2)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-2) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-2) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-2)))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]-2)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-82) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-82)))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]-82)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-79) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-79) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-79)))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]-79)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+79) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+79) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+79)))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]+79)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+82) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+82))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]+82)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-162) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-162))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]-162)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-160) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-160) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-160))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]-160)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-161) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-161) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-161))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]-161)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+162) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+162))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]+162)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+160) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+160) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+160))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]+160)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+161) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+161) %in% prcp_vector$cell &
           rf_abiotic_data$grid_cell[i] %in% prcp_vector$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+161))]) / 
      exp(abs((prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i]+161)]) - (prcp_vector$mean_prcp[which(prcp_vector$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) 
}

# water ring 2

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$water_ring_2[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                        rf_abiotic_data$year[i] > 2010 &
                                        (rf_abiotic_data$grid_cell[i]+2) %in% rf_abiotic_data$grid_cell &
                                        (rf_abiotic_data$grid_cell[i]+2) %in% lc_vectors$cell &
                                        rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+2)))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+2)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-2) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-2) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-2)))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-2)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-82) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-82)))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-82)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-79) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-79) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-79)))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-79)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+79) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+79) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+79)))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+79)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+82) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+82))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+82)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-162) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-162))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-162)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-160) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-160) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-160))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-160)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-161) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-161) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-161))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-161)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+162) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+162))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+162)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+160) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+160) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+160))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+160)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+161) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+161) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+161))]) / 
      exp(abs((lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+161)]) - (lc_vectors$mean_water[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) 
}

# forest ring 2

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$forest_ring_2[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                         rf_abiotic_data$year[i] > 2010 &
                                         (rf_abiotic_data$grid_cell[i]+2) %in% rf_abiotic_data$grid_cell &
                                         (rf_abiotic_data$grid_cell[i]+2) %in% lc_vectors$cell &
                                         rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+2)))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+2)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-2) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-2) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-2)))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-2)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-82) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-82)))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-82)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-79) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-79) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-79)))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-79)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+79) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+79) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+79)))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+79)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+82) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+82))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+82)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-162) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-162))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-162)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-160) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-160) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-160))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-160)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-161) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-161) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-161))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-161)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+162) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+162))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+162)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+160) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+160) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+160))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+160)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+161) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+161) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+161))]) / 
      exp(abs((lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+161)]) - (lc_vectors$mean_forest[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) 
}

# grass ring 2

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$grass_ring_2[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                         rf_abiotic_data$year[i] > 2010 &
                                         (rf_abiotic_data$grid_cell[i]+2) %in% rf_abiotic_data$grid_cell &
                                         (rf_abiotic_data$grid_cell[i]+2) %in% lc_vectors$cell &
                                         rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+2)))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+2)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-2) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-2) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-2)))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-2)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-82) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-82)))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-82)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-79) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-79) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-79)))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-79)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+79) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+79) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+79)))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+79)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+82) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+82))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+82)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-162) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-162))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-162)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-160) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-160) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-160))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-160)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-161) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-161) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-161))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-161)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+162) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+162))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+162)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+160) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+160) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+160))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+160)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+161) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+161) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+161))]) / 
      exp(abs((lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+161)]) - (lc_vectors$mean_grass[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) 
}

# wetland ring 2

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$wetland_ring_2[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                         rf_abiotic_data$year[i] > 2010 &
                                         (rf_abiotic_data$grid_cell[i]+2) %in% rf_abiotic_data$grid_cell &
                                         (rf_abiotic_data$grid_cell[i]+2) %in% lc_vectors$cell &
                                         rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+2)))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+2)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-2) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-2) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-2)))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-2)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-82) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-82)))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-82)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-79) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-79) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-79)))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-79)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+79) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+79) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+79)))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+79)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+82) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+82))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+82)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-162) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-162))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-162)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-160) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-160) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-160))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-160)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-161) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-161) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-161))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-161)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+162) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+162))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+162)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+160) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+160) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+160))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+160)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+161) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+161) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+161))]) / 
      exp(abs((lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+161)]) - (lc_vectors$mean_wetland[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) 
}

# farming ring 2

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$farming_ring_2[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                         rf_abiotic_data$year[i] > 2010 &
                                         (rf_abiotic_data$grid_cell[i]+2) %in% rf_abiotic_data$grid_cell &
                                         (rf_abiotic_data$grid_cell[i]+2) %in% lc_vectors$cell &
                                         rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+2)))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+2)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-2) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-2) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-2)))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-2)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-82) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-82)))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-82)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-79) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-79) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-79)))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-79)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+79) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+79) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+79)))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+79)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+82) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+82))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+82)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-162) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-162))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-162)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-160) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-160) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-160))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-160)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-161) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-161) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-161))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-161)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+162) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+162))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+162)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+160) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+160) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+160))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+160)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+161) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+161) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+161))]) / 
      exp(abs((lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+161)]) - (lc_vectors$mean_farming[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) 
}

# urban ring 2

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$urban_ring_2[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                         rf_abiotic_data$year[i] > 2010 &
                                         (rf_abiotic_data$grid_cell[i]+2) %in% rf_abiotic_data$grid_cell &
                                         (rf_abiotic_data$grid_cell[i]+2) %in% lc_vectors$cell &
                                         rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+2)))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+2)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-2) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-2) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-2)))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-2)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-82) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-82)))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-82)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-79) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-79) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-79)))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-79)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+79) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+79) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+79)))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+79)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+82) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+82))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+82)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-162) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-162))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-162)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-160) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-160) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-160))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-160)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-161) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-161) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-161))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-161)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+162) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+162))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+162)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+160) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+160) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+160))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+160)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+161) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+161) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+161))]) / 
      exp(abs((lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+161)]) - (lc_vectors$mean_urban[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) 
}

# barren ring 2

for(i in 1:nrow(rf_abiotic_data)) {
  rf_abiotic_data$barren_ring_2[i] <- (if(rf_abiotic_data$enough_checklists[i] == 1 & 
                                         rf_abiotic_data$year[i] > 2010 &
                                         (rf_abiotic_data$grid_cell[i]+2) %in% rf_abiotic_data$grid_cell &
                                         (rf_abiotic_data$grid_cell[i]+2) %in% lc_vectors$cell &
                                         rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+2)))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+2)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-2) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-2) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-2)))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-2)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-82) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-82)))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-82)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-79) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-79) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]-79)))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-79)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+79) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+79) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == (rf_abiotic_data$grid_cell[i]+79)))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+79)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+82) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+82) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+82))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+82)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-162) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-162))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-162)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-160) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-160) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-160))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-160)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]-161) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]-161) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]-161))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]-161)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+162) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+162) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+162))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+162)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+160) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+160) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+160))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+160)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) + (if(rf_abiotic_data$enough_checklists[i] == 1 & 
           rf_abiotic_data$year[i] > 2010 &
           (rf_abiotic_data$grid_cell[i]+161) %in% rf_abiotic_data$grid_cell &
           (rf_abiotic_data$grid_cell[i]+161) %in% lc_vectors$cell &
           rf_abiotic_data$grid_cell[i] %in% lc_vectors$cell) {
    (rf_abiotic_data$establishment_score[which((rf_abiotic_data$year == (rf_abiotic_data$year[i]-1)) & 
                                                 (rf_abiotic_data$grid_cell == rf_abiotic_data$grid_cell[i]+161))]) / 
      exp(abs((lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i]+161)]) - (lc_vectors$mean_barren[which(lc_vectors$cell == rf_abiotic_data$grid_cell[i])])))
  } else {
    0
  }) 
}

# save resultant abiotic rf input variable dataframe in a csv file

write.csv(rf_abiotic_data, "data/rf_model_abiotic_data_texas.csv")
