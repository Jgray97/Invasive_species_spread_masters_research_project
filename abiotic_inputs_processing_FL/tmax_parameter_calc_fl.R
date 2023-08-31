# CALCULATING TEMPERATURE RF PARAMETER VALUES FOR GRID CELLS ACROSS FLORIDA
# Author = John Gray
# Email = johnpatrickgray97@gmail.com
# Last edit = 08/07/23

# Load in packages + setwd() ----

library(terra)
library(dggridR)
library(ggplot2)

setwd("D:/John_Gray_research_project/ebd_NZ_relMar-2023")

# Create hexagon grid

dggs <- dgconstruct(res = 8, projection = "ISEA", metric = TRUE, resround = 'nearest')

# Load in temperature raster for the year 2000

tmax_2000 <- rast('tmax/daymet_v4_tmax_annavg_na_2000.tif')

# Load in another spatial object which I know is in the right coordinate 
# reference system and project temperature raster onto its projection
monk_map <- terra::vect("Monk_Parakeet")

tmax_2000 <- project(tmax_2000, crs(monk_map))

# turn temperature raster object into a dataframe to allow assignment of grid 
# cell values
tmax_2000_df <- as.data.frame(tmax_2000, xy = TRUE)

# Assign a grid cell value to each entry in dataframe
tmax_2000_df$cell <- dgGEO_to_SEQNUM(dggs, tmax_2000_df$x, tmax_2000_df$y)$seqnum

# Load in establishment tracker to provide list of Florida grid cells
establishment_tracker <- read.csv("data/FINAL-VERSION-FL_egoose_establishment_ct-50_met-0.005.csv")

# Filter temp dataframe to just include Florida cells
tmax_2000_florida <- filter(tmax_2000_df, tmax_2000_df$cell %in% establishment_tracker$grid_cell)

# change col names
colnames(tmax_2000_florida) <- c("lon", "lat", "tmax", "cell")

# calculate average temperature value for each cell in Florida
tmax_cells_2000 <- tmax_2000_florida %>%
  group_by(cell) %>%
  summarise(mean_tmax = mean(tmax))

# Create a new column in the dataframe which lists the appropriate year for the
# entry
tmax_cells_2000$year <- 2000

# create empty datalist to be filled with each year's dataframe
data_list <- list()

# For loop populates list with cell-grouped temperature dataframes for each year
for (i in 2001:2022) {
  
  # create object referring to downloaded tmax data for each year depending on i
  file_name <- paste0("tmax/daymet_v4_tmax_annavg_na_", i, ".tif")
  
  #Generate raster object for temperature levels across North America
  tmax <- rast(file_name)
  
  # Project raster object onto the right coordinate reference system
  tmax_proj <- project(tmax, crs(monk_map))
  
  # turn raster object into a dataframe so that grid cell values can be assigned
  tmax_df <- as.data.frame(tmax_proj, xy = TRUE)
  
  # Assign grid cell to each point in dataframe
  tmax_df$cell <- dgGEO_to_SEQNUM(dggs, tmax_df$x, tmax_df$y)$seqnum
  
  # Filter to just include entries belonging to cells in Florida
  tmax_florida <- filter(tmax_df, tmax_df$cell %in% establishment_tracker$grid_cell)
  
  # Update column names
  colnames(tmax_florida) <- c("lon", "lat", "tmax", "cell")
  
  # Calculate average values for each Florida cell and store in a dataframe
  tmax_cells <- tmax_florida %>%
    group_by(cell) %>%
    summarise(mean_tmax = mean(tmax))
  
  # add column which tells you the year for a particular entry
  tmax_cells$year <- i
  
  # Add the resultant dataframe to the overarching datalist
  data_list[[i-2000]] <- tmax_cells 
}

# rename datalist items to individual dataframes
tmax_cells_2001 <- data_list[[1]]
tmax_cells_2002 <- data_list[[2]]
tmax_cells_2003 <- data_list[[3]]
tmax_cells_2004 <- data_list[[4]]
tmax_cells_2005 <- data_list[[5]]
tmax_cells_2006 <- data_list[[6]]
tmax_cells_2007 <- data_list[[7]]
tmax_cells_2008 <- data_list[[8]]
tmax_cells_2009 <- data_list[[9]]
tmax_cells_2010 <- data_list[[10]]
tmax_cells_2011 <- data_list[[11]]
tmax_cells_2012 <- data_list[[12]]
tmax_cells_2013 <- data_list[[13]]
tmax_cells_2014 <- data_list[[14]]
tmax_cells_2015 <- data_list[[15]]
tmax_cells_2016 <- data_list[[16]]
tmax_cells_2017 <- data_list[[17]]
tmax_cells_2018 <- data_list[[18]]
tmax_cells_2019 <- data_list[[19]]
tmax_cells_2020 <- data_list[[20]]
tmax_cells_2021 <- data_list[[21]]
tmax_cells_2022 <- data_list[[22]]

# create one big dataframe from all the individual yearly dataframes
tmax_cells <- rbind(tmax_cells_2000, tmax_cells_2001, tmax_cells_2002, 
                    tmax_cells_2003, tmax_cells_2004, tmax_cells_2005,
                    tmax_cells_2006, tmax_cells_2007, tmax_cells_2008,
                    tmax_cells_2009, tmax_cells_2010, tmax_cells_2011,
                    tmax_cells_2012, tmax_cells_2013, tmax_cells_2014,
                    tmax_cells_2015, tmax_cells_2016, tmax_cells_2017,
                    tmax_cells_2018, tmax_cells_2019, tmax_cells_2020,
                    tmax_cells_2021, tmax_cells_2022)

# create final dataframe where temperature in each cell corresponds to average 
# value for that cell from year 2000 to 2022
tmax_cells_overall <- tmax_cells %>%
  group_by(cell) %>%
  summarise(mean_tmax = mean(mean_tmax))

# save resultant dataframe in csv file
write.csv(tmax_cells_overall, "data/tmax_feature_vector.csv")
