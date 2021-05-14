# Extracting environmental time series from HLS data
Collection of scripts to process harmonized Landsat 8 & Sentinel 2 data to study mosquito habitats

1) calculate spectral indices (NDVI, NDMI) from a multiple observations of Sentinel 2 imagery
2) extract index values at mosquito trap locations, given my a shapefile with trap locatoins, and convert time series data into a .csv file
3) interpolate missing observations from the spectral index time series data via a robustified linear regression model
