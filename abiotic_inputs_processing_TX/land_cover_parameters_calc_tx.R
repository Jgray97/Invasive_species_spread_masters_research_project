# CALCULATING LAND COVER RF PARAMETER VALUES FOR GRID CELLS ACROSS TEXAS
# Author = John Gray
# Email = greyjohn15@gmail.com
# Last edit = 30/08/23

# Load in packages + setwd() ----

library(terra)
library(dggridR)
library(ggplot2)
library(tidyverse)
setwd("D:/John_Gray_research_project/ebd_NZ_relMar-2023")

# Load in first case raster (just a test to make sure everything's functioning
# as it should)

land_cover_2001 <- rast('land_cover_data_texas/MCD12Q1.061_LC_Type2_doy2001001_aid0001.tif')

# create hexagon grid

dggs <- dgconstruct(res = 8, projection = "ISEA", metric = TRUE, resround = 'nearest')

# assign cell values to each point

lc_2001_df <- as.data.frame(land_cover_2001, xy = TRUE)

lc_2001_df$cell <- dgGEO_to_SEQNUM(dggs, lc_2001_df$x, lc_2001_df$y)$seqnum

# filter to just area of interest

establishment_tracker <- read.csv("data/TX_egoose_establishment_ct-50_met-0.005.csv")

lc_2001_texas <- filter(lc_2001_df, lc_2001_df$cell %in% establishment_tracker$grid_cell)

# change col names

colnames(lc_2001_texas) <- c("lon", "lat", "land_cover_type", "cell")

# plot to check it looks ok

ggplot() +
  geom_tile(data = lc_2001_texas, aes(x = lon, y = lat, fill = land_cover_type))

# Now calculate average land cover proportion for each cell type in texas in
# each year

# set up a list for each years texas cell land cover dataframes to populate

data_list <- list()

for (i in 2001:2021) {
  
  # create filename object to iteratively load in each year's rasters
  file_name <- paste0("land_cover_data/MCD12Q1.061_LC_Type2_doy", i, "001_aid0001.tif")
  
  # load in raster object for the land cover type data
  lc <- rast(file_name)
  
  # turn raster into a dataframe so that cell values can be assigned
  lc_df <- as.data.frame(lc, xy = TRUE)
  
  # assign grid cell value to each entry in dataframe
  lc_df$cell <- dgGEO_to_SEQNUM(dggs, lc_df$x, lc_df$y)$seqnum
  
  # filter dataframe to just include entries in Texas
  lc_texas <- filter(lc_df, lc_df$cell %in% establishment_tracker$grid_cell)
  
  # rename dataframe column titles
  colnames(lc_texas) <- c("lon", "lat", "land_cover_type", "cell")
  
  # add in value column and then pivot wider to create columns for each land
  # land cover type with a 1 value in that column if raster point belongs to
  # that land cover type
  lc_texas$value <- 1
  
  lc_texas_wider <- lc_texas %>%
    pivot_wider(names_from = land_cover_type, values_from = value)
  
  # change na values to 0 for columns which don't correspond to the correct land
  # cover type
  lc_texas_wider[is.na(lc_texas_wider)]<-0
  
  # Create columns for grouped land cover types, if function accounts for the
  # fact that a specific land cover type might not be present in the data
  lc_texas_wider$water <- if("0" %in% colnames(lc_texas_wider)) {
    lc_texas_wider$"0"
  } else {
    0
  }
  
  lc_texas_wider$forest <- if ("1" %in% colnames(lc_texas_wider)) {
    lc_texas_wider$"1"
  } else {
    0
    } + if ("2" %in% colnames(lc_texas_wider)) {
      lc_texas_wider$"2"
    } else {
      0
      } + if( "3" %in% colnames(lc_texas_wider)) {
        lc_texas_wider$"3"
      } else {
        0
        } + if( "4" %in% colnames(lc_texas_wider)) {
          lc_texas_wider$"4"
        } else {
          0
          } + if("5" %in% colnames(lc_texas_wider)) {
            lc_texas_wider$"5"
          } else {
            0
          }
  
  lc_texas_wider$grass_etc <- if ("6" %in% colnames(lc_texas_wider)) {
    lc_texas_wider$"6"
    } else {
    0
    } + if ("7" %in% colnames(lc_texas_wider)) {
      lc_texas_wider$"7"
      } else {
      0
      } + if( "8" %in% colnames(lc_texas_wider)) {
        lc_texas_wider$"8"
        } else {
        0
        } + if( "9" %in% colnames(lc_texas_wider)) {
          lc_texas_wider$"9"
          } else {
          0
          } + if("10" %in% colnames(lc_texas_wider)) {
            lc_texas_wider$"10"
            } else {
            0
            }
  
  lc_texas_wider$wetland <- if ("11" %in% colnames(lc_texas_wider)) {
    lc_texas_wider$"11"
  } else {
    0
  }
  
  
  lc_texas_wider$farming <- if ("12" %in% colnames(lc_texas_wider)) {
    lc_texas_wider$"12"
  } else {
    0
  } + if ("14" %in% colnames(lc_texas_wider)) {
    lc_texas_wider$"14"
  } else {
    0
  }
  
  lc_texas_wider$urban <- lc_texas_wider$"13"
  
  lc_texas_wider$barren <- if ("15" %in% colnames(lc_texas_wider)) {
    lc_texas_wider$"15"
  } else {
    0
  }
  
  # use group_function to calculate proportion of land cover type coverage in each grid cell
  lc_cells <- lc_texas_wider %>%
    group_by(cell) %>%
    summarise(water_cover = mean(water), forest_cover = mean(forest), 
              grass_etc_cover = mean(grass_etc), wetland_cover = mean(wetland), 
              farming_cover = mean(farming), urban_cover = mean(urban), 
              barren_cover = mean(barren))
  
  # Create column according to year of data
  lc_cells$year <- i
  
  # fill datalist with each year's data
  data_list[[i-2000]] <- lc_cells
}

# Turn each entry in data list into its own individual dataframe object
lc_cells_2001 <- data_list[[1]]
lc_cells_2002 <- data_list[[2]]
lc_cells_2003 <- data_list[[3]]
lc_cells_2004 <- data_list[[4]]
lc_cells_2005 <- data_list[[5]]
lc_cells_2006 <- data_list[[6]]
lc_cells_2007 <- data_list[[7]]
lc_cells_2008 <- data_list[[8]]
lc_cells_2009 <- data_list[[9]]
lc_cells_2010 <- data_list[[10]]
lc_cells_2011 <- data_list[[11]]
lc_cells_2012 <- data_list[[12]]
lc_cells_2013 <- data_list[[13]]
lc_cells_2014 <- data_list[[14]]
lc_cells_2015 <- data_list[[15]]
lc_cells_2016 <- data_list[[16]]
lc_cells_2017 <- data_list[[17]]
lc_cells_2018 <- data_list[[18]]
lc_cells_2019 <- data_list[[19]]
lc_cells_2020 <- data_list[[20]]
lc_cells_2021 <- data_list[[21]]

# rowbind objects for each year to produce 1 big dataframe
lc_cells <- rbind(lc_cells_2001, lc_cells_2002, 
                    lc_cells_2003, lc_cells_2004, lc_cells_2005,
                    lc_cells_2006, lc_cells_2007, lc_cells_2008,
                    lc_cells_2009, lc_cells_2010, lc_cells_2011,
                    lc_cells_2012, lc_cells_2013, lc_cells_2014,
                    lc_cells_2015, lc_cells_2016, lc_cells_2017,
                    lc_cells_2018, lc_cells_2019, lc_cells_2020,
                    lc_cells_2021)

# Calculate average coverage of each land cover type for each cell across across
# all years
lc_cells_overall <- lc_cells %>%
  group_by(cell) %>%
  summarise(mean_water = mean(water_cover), mean_forest = mean(forest_cover), 
            mean_grass_etc = mean(grass_etc_cover), 
            mean_wetland = mean(wetland_cover), mean_farming = mean(farming_cover),
            mean_urban = mean(urban_cover), mean_barren = mean(barren_cover))

# save resultant land cover type coverage dataframe for each cell in Texas
# as a csv file
write.csv(lc_cells_overall, "data/lc_feature_vectors.csv")
