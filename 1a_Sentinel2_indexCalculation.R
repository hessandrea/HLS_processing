#This script loads sentinel imagery from the harmonized landsat-sentinel dataset, 
#extracts individual bands from hdf filtes
#screens out clouds via the qa filter, 
#saves cloud free tif images

rm(list=ls())

start.time <- Sys.time()


library(rgdal)
library(gdalUtils)
library(raster)
library(R.utils)
library(binaryLogic)
library(tidyverse)
library(sf)
library(chron)


### ----------------- user input 
#the path to a study area shapefile
aoi_sf <- st_read("./data_crunch/Norman_shapes/NormanDissolved.shp")

#full path name of one .hdf subset, simply for the purpose of making a test map to check for correct overlaying
sds <- get_subdatasets('./data/SatData/HLS/L30/2019/14SPD/HLS.L30.T14SPD.2019004.v1.4.hdf')

#the path where all S30 .hdf files are stored
sentinel_path <-"./data/SatData/HLS/S30"

#the path where all index ndvi/ndwi will be stored
out_folder <- "./data_Crunch/SatData/HLS/2_indices/S30"

#the path to a temp folder for 
temp_folder <- "./data_Crunch/SatData/HLS/temp"

### ------------------- end of user input



####
gdal_translate(sds[1], dst_dataset = paste0(temp_folder, "/projection_band.tif"))
rast <- raster("./data_Crunch/SatData/HLS/temp/projectband.tif")

aoi_proj <- st_transform(aoi_sf, crs = crs(rast))

plot(rast)
plot(aoi_proj, add=T)



#loop through all .hdf files
s_files <- list.files(path = sentinel_path, pattern = "hdf$", recursive = TRUE, full.names = T)
s_files



# Sentinel index creation -------------------------------------------------
i<-4
for (i in 1:length(s_files)){
  #looping through each hdf file
  myfile <- s_files[i]
  filename <- str_sub(myfile,-31, -4)
  print(filename)
  
  jul <- as.integer(str_sub(myfile,-12, -10))
  year <- as.integer(str_sub(myfile,-16, -13))
  orig <- c(01, 01, year)
  tile_date <- unlist(month.day.year(jul, orig))
  tile_date
  
  sds <- get_subdatasets(myfile)
  
  #subsetting QA band
  s_qa <- sds[14]
  s_qa_name <- paste0(temp_folder, "/qa_band_crop.tif")
  gdal_translate(s_qa, dst_dataset = s_qa_name)
  s_qa_rast <- raster(s_qa_name)
  s_qa_rast_crop <- crop(s_qa_rast, aoi_proj)
  
  ### creating cloud mask
  #QA bits (bit layer is unint8, fill value 255)
  qa_matrix <- as.matrix(s_qa_rast_crop)
  
  #converting integer into bit values
  first_k_bits <- function(int, k=8, reverse=T) {
    ## Bit 0 is the least significant (read bit words right to left)
    integer_vector <- as.integer(intToBits(int))[1:k]
    if(reverse) integer_vector <- rev(integer_vector)
    return(paste(as.character(integer_vector), collapse=""))
  }
  #creating data frame with bit values
  df <- data.frame(bits=sapply(s_qa_rast_crop[], function(x) first_k_bits(x, k=8, reverse=T)))
  df$bits[s_qa_rast_crop[] == 255] <- NA  
  
  #adding 0/1 columns for each cloud class
  df$cirrus <- substring(df$bits, 8, 8)   #cirrus: bit number 0 (0 no /1 yes). reading right to left, puts it at position 8
  df$clouds <- substring(df$bits, 7, 7)   #clouds: bit number 1 (0 no /1 yes)
  df$adjacent <- substring(df$bits, 6, 6) #adjacent: bit number 2 (0 no/1 yes)
  df <- df %>%
    mutate(screen_value = case_when(
      cirrus == "1" | clouds == "1"  | adjacent == "1" ~ 1
    ))
  #creating an overall cloud screen column
  df$screen_value[is.na(df$screen_value)] <- 0
  
  #### to check how many pixels of each value are present in each file
  # table(df$cloud) 
  # table(df$adjacent)
  # table(df$cirrus)
  # table(df$screen_value)
  
  #create screening matrix and raster
  dims <- dim(s_qa_rast_crop) 
  mask <- raster(ncol=dims[2], nrow=dims[1], 
                 xmn=extent(s_qa_rast_crop)[1], xmx=extent(s_qa_rast_crop)[2],
                 ymn=extent(s_qa_rast_crop)[3], ymx=extent(s_qa_rast_crop)[4], 
                 crs=crs(s_qa_rast_crop))
  values(mask) <- df$screen_value
  #as.raster(mask)

  
  cirrus <- as.integer(df$cirrus)
  cirrus_mask <- raster(ncol=dims[2], nrow=dims[1], 
                 xmn=extent(s_qa_rast_crop)[1], xmx=extent(s_qa_rast_crop)[2],
                 ymn=extent(s_qa_rast_crop)[3], ymx=extent(s_qa_rast_crop)[4], 
                 crs=crs(s_qa_rast_crop))
  values(cirrus_mask) <- cirrus

  
  clouds <- as.integer(df$clouds)
  clouds_mask <- raster(ncol=dims[2], nrow=dims[1], 
                        xmn=extent(s_qa_rast_crop)[1], xmx=extent(s_qa_rast_crop)[2],
                        ymn=extent(s_qa_rast_crop)[3], ymx=extent(s_qa_rast_crop)[4], 
                        crs=crs(s_qa_rast_crop))
  values(clouds_mask) <- clouds

  ### subset relevant bands
  s_blue <- sds[2]
  s_blue_name <- paste0(temp_folder, "/blue_band_crop.tif")
  gdal_translate(s_blue, dst_dataset = s_blue_name)
  s_blue_rast <- raster(s_blue_name)
  s_blue_rast_crop <- crop(s_blue_rast, aoi_proj)
  
  s_green <- sds[3]
  s_green_name <- paste0(temp_folder, "/green_band_crop.tif")
  gdal_translate(s_green, dst_dataset = s_green_name)
  s_green_rast <- raster(s_green_name)
  s_green_rast_crop <- crop(s_green_rast, aoi_proj)
  #plot(s_green_rast_crop)
  
  s_red <- sds[4]
  s_red_name <- paste0(temp_folder, "/red_band_crop.tif")
  gdal_translate(s_red, dst_dataset = s_red_name)
  s_red_rast <- raster(s_red_name)
  s_red_rast_crop <- crop(s_red_rast, aoi_proj)
  #plot(s_red_rast_crop)
  
  s_nir <- sds[8]  
  s_nir_name <- paste0(temp_folder, "/temp/nir_band_crop.tif")
  gdal_translate(s_nir, dst_dataset = s_nir_name)
  s_nir_rast <- raster(s_nir_name)
  s_nir_rast_crop <- crop(s_nir_rast, aoi_proj)
  #plot(s_nir_rast_crop)
  
  s_swir <- sds[12] #band 11 in sentinel2 is layer 12 in sds
  s_swir_name <- paste0(temp_folder, "/swir_band_crop.tif")
  gdal_translate(s_swir, dst_dataset = s_swir_name)
  s_swir_rast <- raster(s_swir_name)
  s_swir_rast_crop <- crop(s_swir_rast, aoi_proj)
  #plot(s_swir_rast_crop)
  
  ### QA plotting
  #plotting
  #plot(s_qa_rast_crop)
  #plot(cirrus_mask)
  #plot(clouds_mask)
  #plot(mask)
  
  #rgb_stack <- stack (s_red_rast_crop, s_green_rast_crop, s_blue_rast_crop)
  #plotRGB(rgb_stack, r=1, g=2, b=3, stretch = "hist")
  
  #false_clor_stack <- stack(s_nir_rast_crop, s_red_rast_crop, s_green_rast_crop)
  #plotRGB(false_clor_stack, r=1, g=2, b=3, stretch = "hist")
  
  
  
  #mask bands
  green_masked <- mask(s_green_rast_crop, mask, inverse=FALSE, 
                      maskvalue=1, updatevalue=NA)
  red_masked <- mask(s_red_rast_crop, mask, inverse=FALSE, 
                         maskvalue=1, updatevalue=NA)
  nir_masked <- mask(s_nir_rast_crop, mask, inverse=FALSE, 
                     maskvalue=1, updatevalue=NA)
  swir_masked <- mask(s_swir_rast_crop, mask, inverse=FALSE, 
                     maskvalue=1, updatevalue=NA)
  
  
  ## calculate indices  
  ndvi <- (nir_masked-red_masked)/(nir_masked+red_masked)    #for sentinal 2: (B8-B4)/(B8+B4)
  ndwi_ga <- (nir_masked - swir_masked)/(nir_masked + swir_masked)
  ndwi_mc <- (green_masked - nir_masked)/(green_masked + nir_masked)
  ndwi_xu <- (green_masked - swir_masked)/(green_masked + swir_masked)


  #stack indices and save
  index_stack <- stack(ndvi, ndwi_ga, ndwi_mc, ndwi_xu) 
  outname <- paste0(out_folder, "/", filename, "indices.crop.tif")
  writeRaster(index_stack, filename=outname)
}


end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
