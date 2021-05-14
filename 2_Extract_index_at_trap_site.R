#This script reads all HLS index files, 
#extracts index values at trap sites, 
### (for this, the user neefs to change the bands needed for index calcualtion based on whether landsat or sentinel data is being read)
#and saves all data to a time series data frame

rm(list = ls(all = TRUE))


#libraries
library(chron)
library(tidyverse)
library(sf)
library(raster)

###--------------------- user input
index_path <- "./data_Crunch/SatData/HLS/2_indices"

trap_shp <- "./data/trapping_data/traps_2019_2020.shp"
#if the shapefile with the trap location has a column for "trap_site", "site_id" etc, this needs to be changed in line 35

#buffer around trap site for sat data summary. Enter a value in m
buffer_size <- 1000

#name for output csv file. This file will contain the average NDVI/NMDI values within a the given buffer buffer for each trap site and trap date
out_name <- "./data_Crunch/SatData/HLS/HLS_time_series_1000m_cloudmasked.csv"


### --------------------- end of user input





#set tif image source
tif_files <- list.files(path = index_path, pattern = "tif$", recursive = TRUE, full.names = T)

#read shapefile
traps_sf <- st_read(trap_shp) %>%
  dplyr::select(c("site_id")) %>%
  st_difference()

traps <- st_transform(traps_sf, crs =  "+proj=utm +zone=14 +datum=WGS84 +units=m +no_defs")

unique_sites <- unique(traps_sf$site_id)

#create empty dataframe as template
out_dat <- data.frame()


i <- 4
j <- 3

for (i in 1:length(tif_files)){
  #looping through each hdf file
  myfile <- tif_files[i]
  print(myfile)
  mysensor <- str_sub(myfile,-40, -38)
  #mydate <- str_sub(myfile,-29, -23) 
  myyear <- as.numeric(str_sub(myfile,-29, -26))
  mydoy <- as.numeric(str_sub(myfile,-25, -23))
  
  orig <- c(01, 01, myyear)
  tile_date_vec <- unlist(month.day.year(mydoy, orig))
  tile_date <- paste0(tile_date_vec[3],  "-", tile_date_vec[1], "-", tile_date_vec[2])
  print(tile_date)

  myraster <- stack(myfile)
  
  ndvi <- myraster[[1]]
  ndwi_ga <- myraster[[2]]
  ndwi_mc <- myraster[[3]]
  ndwi_xu <- myraster[[4]]
  
  #plot(ndvi)
  #plot(traps, colour= 'red', add=T)
  
  indices <- c("ndvi", "ndwi_ga", "ndwi_mc", "ndwi_xu", "ndvi_1000", "ndwi_ga_1000", "ndwi_mc_1000", "ndwi_xu_1000")
  ndvi_traps <- raster::extract(ndvi, traps)
  ndwi_ga_traps <- raster::extract(ndwi_ga, traps)
  ndwi_mc_traps <- raster::extract(ndwi_mc, traps)
  ndwi_xu_traps <- raster::extract(ndwi_xu, traps)
  ndvi_buffer <- raster::extract(ndvi, traps, buffer = buffer_size, fun=mean)
  ndwi_ga_buffer <- raster::extract(ndwi_ga, traps, buffer = buffer_size, fun=mean)
  ndwi_mc_buffer <- raster::extract(ndwi_mc, traps, buffer = buffer_size, fun=mean)
  ndwi_xu_buffer <- raster::extract(ndwi_xu, traps, buffer = buffer_size, fun=mean)
  
  temp_dat <- data.frame(rbind(ndvi_traps, ndwi_ga_traps, ndwi_mc_traps, ndwi_xu_traps, 
                               ndvi_buffer, ndwi_ga_buffer, ndwi_mc_buffer, ndwi_xu_buffer))
  temp_dat <- temp_dat %>% 
    mutate(obs_date = tile_date, 
           index = indices, 
           sensor = mysensor)

  out_dat <- rbind(out_dat, temp_dat)
}



colnames(out_dat) <- c(unique_sites, 
                       "obs_date", "index", "sensor")

out_long <- pivot_longer(out_dat, cols = starts_with("T"), 
                         names_to = "trap", 
                         values_to = "value")
 
print(out_long)
write_csv(out_long, out_name)
