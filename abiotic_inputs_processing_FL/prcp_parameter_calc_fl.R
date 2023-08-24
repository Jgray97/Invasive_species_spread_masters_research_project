### CALCULATING PRCP FOR ALL GRID CELLS

# Author = John Gray
# Email = johnpatrickgray97@gmail.com
# Last edit = 08/07/23

# Load in packages + setwd() ----

library(terra)
library(dggridR)
library(ggplot2)

setwd("D:/John_Gray_research_project/ebd_NZ_relMar-2023")

dggs <- dgconstruct(res = 8, projection = "ISEA", metric = TRUE, resround = 'nearest')

monk_map <- terra::vect("Monk_Parakeet")

establishment_tracker <- read.csv("data/FINAL-FINAL-VERSION-FL_egoose_establishment_ct-50_met-0.005.csv")

data_list_prcp <- list()

for (i in 2000:2022) {
  
  file_name <- paste0("prcp/daymet_v4_prcp_annttl_na_", i, ".tif")
  
  prcp <- rast(file_name)
  
  prcp_proj <- project(prcp, crs(monk_map))
  
  prcp_df <- as.data.frame(prcp_proj, xy = TRUE)
  
  prcp_df$cell <- dgGEO_to_SEQNUM(dggs, prcp_df$x, prcp_df$y)$seqnum
  
  prcp_florida <- filter(prcp_df, prcp_df$cell %in% establishment_tracker$grid_cell)
  
  colnames(prcp_florida) <- c("lon", "lat", "prcp", "cell")
  
  prcp_cells <- prcp_florida %>%
    group_by(cell) %>%
    summarise(mean_prcp = mean(prcp))
  
  prcp_cells$year <- i
  
  data_list_prcp[[i-1999]] <- prcp_cells 
}

## CONTINUE AFTER LUNCH FROM HERE

# rename datalist items to individual dataframes

prcp_cells_2000 <- data_list_prcp[[1]]
prcp_cells_2001 <- data_list_prcp[[2]]
prcp_cells_2002 <- data_list_prcp[[3]]
prcp_cells_2003 <- data_list_prcp[[4]]
prcp_cells_2004 <- data_list_prcp[[5]]
prcp_cells_2005 <- data_list_prcp[[6]]
prcp_cells_2006 <- data_list_prcp[[7]]
prcp_cells_2007 <- data_list_prcp[[8]]
prcp_cells_2008 <- data_list_prcp[[9]]
prcp_cells_2009 <- data_list_prcp[[10]]
prcp_cells_2010 <- data_list_prcp[[11]]
prcp_cells_2011 <- data_list_prcp[[12]]
prcp_cells_2012 <- data_list_prcp[[13]]
prcp_cells_2013 <- data_list_prcp[[14]]
prcp_cells_2014 <- data_list_prcp[[15]]
prcp_cells_2015 <- data_list_prcp[[16]]
prcp_cells_2016 <- data_list_prcp[[17]]
prcp_cells_2017 <- data_list_prcp[[18]]
prcp_cells_2018 <- data_list_prcp[[19]]
prcp_cells_2019 <- data_list_prcp[[20]]
prcp_cells_2020 <- data_list_prcp[[21]]
prcp_cells_2021 <- data_list_prcp[[22]]
prcp_cells_2022 <- data_list_prcp[[23]]

prcp_cells <- rbind(prcp_cells_2000, prcp_cells_2001, prcp_cells_2002, 
                    prcp_cells_2003, prcp_cells_2004, prcp_cells_2005,
                    prcp_cells_2006, prcp_cells_2007, prcp_cells_2008,
                    prcp_cells_2009, prcp_cells_2010, prcp_cells_2011,
                    prcp_cells_2012, prcp_cells_2013, prcp_cells_2014,
                    prcp_cells_2015, prcp_cells_2016, prcp_cells_2017,
                    prcp_cells_2018, prcp_cells_2019, prcp_cells_2020,
                    prcp_cells_2021, prcp_cells_2022)

prcp_cells_overall <- prcp_cells %>%
  group_by(cell) %>%
  summarise(mean_prcp = mean(mean_prcp))

write.csv(prcp_cells_overall, "data/prcp_feature_vector.csv")


