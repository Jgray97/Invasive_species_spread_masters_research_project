### CALCULATING TMAX FOR ALL GRID CELLS

# Author = John Gray
# Email = johnpatrickgray97@gmail.com
# Last edit = 08/07/23

# Load in packages + setwd() ----

library(terra)
library(dggridR)
library(ggplot2)

setwd("D:/John_Gray_research_project/ebd_NZ_relMar-2023")


# Load in grid and temp data

dggs <- dgconstruct(res = 8, projection = "ISEA", metric = TRUE, resround = 'nearest')

tmax_2000 <- rast('tmax/daymet_v4_tmax_annavg_na_2000.tif')

# sort out temp data coord reference system

monk_map <- terra::vect("Monk_Parakeet")

tmax_2000 <- project(tmax_2000, crs(monk_map))

# assign grid cell to each point

tmax_2000_df <- as.data.frame(tmax_2000, xy = TRUE)

tmax_2000_df$cell <- dgGEO_to_SEQNUM(dggs, tmax_2000_df$x, tmax_2000_df$y)$seqnum

# filter temp dataframe to just florida cells

establishment_tracker <- read.csv("data/FINAL-FINAL-VERSION-FL_egoose_establishment_ct-50_met-0.005.csv")

tmax_2000_florida <- filter(tmax_2000_df, tmax_2000_df$cell %in% establishment_tracker$grid_cell)

# change col names

colnames(tmax_2000_florida) <- c("lon", "lat", "tmax", "cell")

# check that it looks ok

ggplot() +
  geom_tile(data = tmax_2000_florida, aes(x = lon, y = lat, fill = tmax)) +
  xlim(-88,-79) + 
  ylim(23,32) +
  scale_fill_gradient(limits = c(20, 31), low = alpha("Blue"), 
                      high = alpha("Red"))

# investigate what's happening for 0 values

min(tmax_2000_florida$tmax)

nrow(filter(tmax_2000_florida, is.na(tmax_2000_florida$tmax)))

summary(tmax_2000_florida$tmax)

# calculate cell values

tmax_cells_2000 <- tmax_2000_florida %>%
  group_by(cell) %>%
  summarise(mean_tmax = mean(tmax))

tmax_cells_2000$year <- 2000

# load in rasters for all years

tmax_2001 <- rast('tmax/daymet_v4_tmax_annavg_na_2001.tif')
tmax_2002 <- rast('tmax/daymet_v4_tmax_annavg_na_2002.tif')
tmax_2003 <- rast('tmax/daymet_v4_tmax_annavg_na_2003.tif')
tmax_2004 <- rast('tmax/daymet_v4_tmax_annavg_na_2004.tif')
tmax_2005 <- rast('tmax/daymet_v4_tmax_annavg_na_2005.tif')
tmax_2006 <- rast('tmax/daymet_v4_tmax_annavg_na_2006.tif')
tmax_2007 <- rast('tmax/daymet_v4_tmax_annavg_na_2007.tif')
tmax_2008 <- rast('tmax/daymet_v4_tmax_annavg_na_2008.tif')
tmax_2009 <- rast('tmax/daymet_v4_tmax_annavg_na_2009.tif')
tmax_2010 <- rast('tmax/daymet_v4_tmax_annavg_na_2010.tif')
tmax_2011 <- rast('tmax/daymet_v4_tmax_annavg_na_2011.tif')
tmax_2012 <- rast('tmax/daymet_v4_tmax_annavg_na_2012.tif')
tmax_2013 <- rast('tmax/daymet_v4_tmax_annavg_na_2013.tif')
tmax_2014 <- rast('tmax/daymet_v4_tmax_annavg_na_2014.tif')
tmax_2015 <- rast('tmax/daymet_v4_tmax_annavg_na_2015.tif')
tmax_2016 <- rast('tmax/daymet_v4_tmax_annavg_na_2016.tif')
tmax_2017 <- rast('tmax/daymet_v4_tmax_annavg_na_2017.tif')
tmax_2018 <- rast('tmax/daymet_v4_tmax_annavg_na_2018.tif')
tmax_2019 <- rast('tmax/daymet_v4_tmax_annavg_na_2019.tif')
tmax_2020 <- rast('tmax/daymet_v4_tmax_annavg_na_2020.tif')
tmax_2021 <- rast('tmax/daymet_v4_tmax_annavg_na_2021.tif')
tmax_2022 <- rast('tmax/daymet_v4_tmax_annavg_na_2022.tif')

data_list <- list()

# create for loop to calculate tmax_cells for each year

for (i in 2001:2022) {
  
  file_name <- paste0("tmax/daymet_v4_tmax_annavg_na_", i, ".tif")
  
  tmax <- rast(file_name)
  
  tmax_proj <- project(tmax, crs(monk_map))
  
  tmax_df <- as.data.frame(tmax_proj, xy = TRUE)
  
  tmax_df$cell <- dgGEO_to_SEQNUM(dggs, tmax_df$x, tmax_df$y)$seqnum
  
  tmax_florida <- filter(tmax_df, tmax_df$cell %in% establishment_tracker$grid_cell)
  
  colnames(tmax_florida) <- c("lon", "lat", "tmax", "cell")
  
  tmax_cells <- tmax_florida %>%
    group_by(cell) %>%
    summarise(mean_tmax = mean(tmax))
  
  tmax_cells$year <- i
  
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

tmax_cells <- rbind(tmax_cells_2000, tmax_cells_2001, tmax_cells_2002, 
                    tmax_cells_2003, tmax_cells_2004, tmax_cells_2005,
                    tmax_cells_2006, tmax_cells_2007, tmax_cells_2008,
                    tmax_cells_2009, tmax_cells_2010, tmax_cells_2011,
                    tmax_cells_2012, tmax_cells_2013, tmax_cells_2014,
                    tmax_cells_2015, tmax_cells_2016, tmax_cells_2017,
                    tmax_cells_2018, tmax_cells_2019, tmax_cells_2020,
                    tmax_cells_2021, tmax_cells_2022)

tmax_cells_overall <- tmax_cells %>%
  group_by(cell) %>%
  summarise(mean_tmax = mean(mean_tmax))

write.csv(tmax_cells_overall, "data/tmax_feature_vector.csv")
