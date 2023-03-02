# Extracting fire danger indices from https://cds.climate.copernicus.eu/cdsapp#!/dataset/cems-fire-historical?tab=overview within fire patches
# This script extracts the maximum value for each index across the fire patch over the fire's lifetime

# Clear data
rm(list = ls())

# Load libraries 
library(sf)
library(raster)
library(lubridate)
library(dplyr)

# Read in shapefiles
setwd("~/BTO projects/Polesia wildfires/fire shapefiles/fire patches")
shp <- st_read("fire_patches.shp")

# Calculate months 3, 6 and 12 months before grid cell dates and match dates to raster filenames by removing the dash
## Name of date column ('mindate') differs from grid cell dateset, indicates the start date of the fire
shp$month3 <- as.Date(shp$mindate) %m-% months(3)
shp$month6 <- as.Date(shp$mindate)  %m-% months(6)
shp$month12 <- as.Date(shp$mindate) %m-% months(12) 

# Create vector of unique lagged dates and remove dash to match raster filenames
lag.dates <- unique(gsub("-", "", c(shp$month3, shp$month6, shp$month12)))

# Create vector of fire dates
dates <- list()
for (i in 1:length(shp$mindate)) {
  dates[[i]] <- gsub("-", "", seq(as.Date(shp$mindate[i]), as.Date(shp$maxdate[i]), by="days")) 
}

## Create vector of fire length (i.e., number of days burned)  
length.vec <- 0
for (i in 1:length(dates)) {
  length.vec[i] <- length(dates[[i]]) 
}

# DROUGHT CODE ##########################

# Name the main folder where the input variables are stored
setwd("~/BTO projects/Polesia wildfires/input variables/dataset-cems-fire-historical-DC")

# Create empty lists
drought.mx <- list() ## Fire date
drought.3.mx <- list() ## 3 months
drought.6.mx <- list() ## 6 months
drought.12.mx <- list() ## 12 months

# Create empty list to import rasters 
filenames <- vector("list", length(dates))

# Loop through each date reading in rasters and calculating max (or mean) drought code over lifetime of fire
for (i in 1:length(dates)) {
  for (t in 1:length.vec[i]) {
    filenames[[i]][t] <-  list.files(path=getwd(),pattern=dates[[i]][t], full.names=TRUE, recursive=TRUE) 
    s1 <- max(stack(filenames[[i]]))
    crs(s1) <- "+proj=longlat +ellps=WGS84 +lon_wrap=180 +datum=WGS84 +no_defs" # Set coordinate reference system
    s1 <- projectRaster(s1, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # Reproject
    drought.mx[[i]] <- as.data.frame(extract(s1, shp[i,], fun = max)) 
  }  
} 

# And lagged DC (3, 6 and 12 months before)
for (i in 1:length(unique(shp$month3))) {
  s2 <- raster(list.files(path=getwd(), pattern = unique(gsub("-", "", shp$month3))[i], full.names=TRUE, recursive=TRUE)) # Find raster file 
  crs(s2) <- "+proj=longlat +ellps=WGS84 +lon_wrap=180 +datum=WGS84 +no_defs" # Set coordinate reference system
  s2 <- projectRaster(s1, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # Reproject
  drought.3.mx[[i]] <- as.data.frame(extract(s2, shp[shp$month3 == unique(shp$month3)[i],], # Extract drought code within each fire patch area
      fun = max)) # Maximum drought code
}
# Repeat for 6 and 12 months before fire
for (i in 1:length(unique(shp$month6))) {
  s3 <- raster(list.files(path=getwd(), pattern = unique(gsub("-", "", shp$month6))[i], full.names=TRUE, recursive=TRUE)) # Find raster file 
  crs(s3) <- "+proj=longlat +ellps=WGS84 +lon_wrap=180 +datum=WGS84 +no_defs" # Set coordinate reference system
  s3 <- projectRaster(s3, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # Reproject
  drought.6.mx[[i]] <- as.data.frame(extract(s3, shp[shp$month6 == unique(shp$month6)[i],], fun = max))
}
for (i in 1:length(unique(shp$month12))) {
  s4 <- raster(list.files(path=getwd(), pattern = unique(gsub("-", "", shp$month12))[i], full.names=TRUE, recursive=TRUE)) # Find raster file 
  crs(s4) <- "+proj=longlat +ellps=WGS84 +lon_wrap=180 +datum=WGS84 +no_defs" # Set coordinate reference system
  s4 <- projectRaster(s4, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # Reproject
  drought.12.mx[[i]] <- as.data.frame(extract(s4, shp[shp$month12 == unique(shp$month12)[i],], fun = max))
}

# Merge dataframes in list    
drought.mx.df <- do.call("rbind", drought.mx)
drought.3.mx.df <- do.call("rbind", drought.3.mx)
drought.6.mx.df <- do.call("rbind", drought.6.mx)
drought.12.mx.df <- do.call("rbind", drought.12.mx)

# Change column names
colnames(drought.mx.df) <- "dc"
colnames(drought.3.mx.df) <- "dc_3"
colnames(drought.6.mx.df) <- "dc_6"
colnames(drought.12.mx.df) <- "dc_12"

# KBDI ####

# Name the main folder where the input variables are stored
setwd("~/BTO projects/Polesia wildfires/input variables/dataset-cems-fire-historical-KBDI")

# Create empty lists
kbdi.mx <- list() 
kbdi.3.mx <- list()
kbdi.6.mx <- list()
kbdi.12.mx <- list() 

# Re-create emtpty list to import rasters 
filenames <- vector("list", length(dates))

# Loop through each date reading in all rasters and calculating the maximum KBDI
for (i in 1:length(dates)) {
  for (t in 1:length.vec[i]) {
    filenames[[i]][t] <-  list.files(path=getwd(),pattern=dates[[i]][t], full.names=TRUE, recursive=TRUE) 
    s5 <- max(stack(filenames[[i]]))
    crs(s5) <- "+proj=longlat +ellps=WGS84 +lon_wrap=180 +datum=WGS84 +no_defs" # Set coordinate reference system
    s5 <- projectRaster(s5, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # Reproject
    s5 <- max(stack(filenames[[i]]))
    kbdi.mx[[i]] <- as.data.frame(extract(s5, shp[i,], fun = max)) 
  }  
}

# And lagged KBDI (3, 6 and 12 months before)
for (i in 1:length(unique(shp$month3))) {
  s6 <- raster(list.files(path=getwd(), pattern = unique(gsub("-", "", shp$month3))[i], full.names=TRUE, recursive=TRUE)) # Find raster file 
  crs(s6) <- "+proj=longlat +ellps=WGS84 +lon_wrap=180 +datum=WGS84 +no_defs" # Set coordinate reference system
  s6 <- projectRaster(s6, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # Reproject
  kbdi.3.mx[[i]] <- as.data.frame(extract(s6, shp[shp$month3 == unique(shp$month3)[i],], fun = max))
}
for (i in 1:length(unique(shp$month6))) {
  s7 <- raster(list.files(path=getwd(), pattern = unique(gsub("-", "", shp$month6))[i], full.names=TRUE, recursive=TRUE)) # Find raster file 
  crs(s7) <- "+proj=longlat +ellps=WGS84 +lon_wrap=180 +datum=WGS84 +no_defs" # Set coordinate reference system
  s7 <- projectRaster(s7, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # Reproject
  kbdi.6.mx[[i]] <- as.data.frame(extract(s7, shp[shp$month6 == unique(shp$month6)[i],], fun = max))
}
for (i in 1:length(unique(shp$month12))) {
  s8 <- raster(list.files(path=getwd(), pattern = unique(gsub("-", "", shp$month12))[i], full.names=TRUE, recursive=TRUE)) # Find raster file 
  crs(s8) <- "+proj=longlat +ellps=WGS84 +lon_wrap=180 +datum=WGS84 +no_defs" # Set coordinate reference system
  s8 <- projectRaster(s8, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # Reproject
  kbdi.12.mx[[i]] <- as.data.frame(extract(s8, shp[shp$month12 == unique(shp$month12)[i],], fun = mean))
}

# Merge dataframes in list    
kbdi.mx.df <- do.call("rbind", kbdi.mx)
kbdi.3.mx.df <- do.call("rbind", kbdi.3.mx)
kbdi.6.mx.df <- do.call("rbind", kbdi.6.mx)
kbdi.12.mx.df <- do.call("rbind", kbdi.12.mx)

# Change column names
colnames(kbdi.mx.df) <- "kbdi"
colnames(kbdi.3.mx.df) <- "kbdi_3"
colnames(kbdi.6.mx.df) <- "kbdi_6"
colnames(kbdi.12.mx.df) <- "kbdi_12"

# FINE FUEL MOISTURE CONTENT ##########################

# Name the main folder where the input variables are stored
setwd("~/BTO projects/Polesia wildfires/input variables/dataset-cems-fire-historical-FFMC")

# Create empty list to fill with FFMC data
ffmc.lst <- list()
filenames <- vector("list", length(dates))

# Loop through each date reading in all rasters and calculating the max FFMC
for (i in 1:length(dates)) {
  for (t in 1:length.vec[i]) {
    filenames[[i]][t] <-  list.files(path=getwd(),pattern=dates[[i]][t],full.names=TRUE,recursive=TRUE) 
    s9 <- max(stack(filenames[[i]]))
    crs(s9) <- "+proj=longlat +ellps=WGS84 +lon_wrap=180 +datum=WGS84 +no_defs" # Set coordinate reference system
    s9 <- projectRaster(s9, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # Reproject
    s9 <- max(stack(filenames[[i]]))
    ffmc.lst[[i]] <- as.data.frame(extract(s9, shp[i,], fun = max)) 
  }  
}

# Merge dataframes in list    
ffmc.df <- do.call("rbind", ffmc.lst)

# Change column name
colnames(ffmc.df)[1] <- "ffmc"

# FWI ##########################

# Name the main folder where the input variables are stored
setwd("~/BTO projects/Polesia wildfires/input variables/dataset-cems-fire-historical-FWI")

# Create empty list 
filenames <- vector("list", length(dates))
fwi.lst <- list()

# Loop through each date reading in all rasters and calculating the maximum over the lifetime of the fire
for (i in 1:length(dates)) {
  for (t in 1:length.vec[i]) {
    filenames[[i]][t] <-  list.files(path=getwd(),pattern=dates[[i]][t],full.names=TRUE,recursive=TRUE) 
    s10 <- max(stack(filenames[[i]]))
    crs(s10) <- "+proj=longlat +ellps=WGS84 +lon_wrap=180 +datum=WGS84 +no_defs" # Set coordinate reference system
    s10 <- projectRaster(s10, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # Reproject
    s10 <- max(stack(filenames[[i]]))
    fwi.lst[[i]] <- as.data.frame(extract(s10, shp[i,], fun = max)) 
  }  
}

# Merge dataframes in list    
fwi.df <- do.call("rbind", fwi.lst)

# Change column name
colnames(fwi.df)[1] <- "fwi"

#### Merge all dataframes ####
dfs <- sapply(.GlobalEnv, is.data.frame) 
df <- do.call(cbind, mget(names(dfs)[dfs]))
# Rename ID column
colnames(df)[4] <- c("z")
df <- select(df, -contains("shp.")) # Remove dates, coordinates and geometry
# Export as csv file
setwd("~/BTO projects/Polesia wildfires/input variables") # Change working directory to where you are storing your output files
write.csv(df, "FWI.csv", row.names = F) ## Outputs are dataframe of fire weather indices and fire patch IDs
