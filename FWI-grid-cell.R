# Extracting fire danger indices from https://cds.climate.copernicus.eu/cdsapp#!/dataset/cems-fire-historical?tab=overview within grid cells

# Clear data
rm(list = ls())

# Load libraries 
library(sf)
library(raster)
library(lubridate)
library(dplyr)

# Read in shapefiles
setwd("~/BTO projects/Polesia wildfires/fire shapefiles/pixels")
pix <- st_read("all_pixels.shp")

# Calculate months 3, 6 and 12 months before grid cell dates and match dates to raster filenames by removing the dash
pix$month3 <- gsub("-", "", as.Date(pix$dates) %m-% months(3))
pix$month6 <- gsub("-", "", as.Date(pix$dates)  %m-% months(6))
pix$month12 <- gsub("-", "", as.Date(pix$dates) %m-% months(12))

# Match grid cell dates to raster filenames by removing the dash
pix$dates <- gsub("-", "", pix$dates)

# Create vector of unique dates
month3 <- unique(pix$month3)
month6 <- unique(pix$month6)
month12 <- unique(pix$month12)
dates <- unique(pix$dates)

# DROUGHT CODE ##########################

# Name the main folder where the input variables are stored
dc.wd <- "~/BTO projects/Polesia wildfires/input variables/dataset-cems-fire-historical-DC"

# Create empty lists
dc.lst <- list() 
dc.3.lst <- list() ## 3 months before grid cell date
dc.6.lst <- list() ## 6 months
dc.12.lst <- list() ## 12 months

# Loop through each date reading in rasters calculating the mean value within grid cells
for (i in 1:length(dates)) {
    dc.lst[[i]] <- as.data.frame(extract(raster(list.files(path = dc.wd, pattern = dates[i], full.names = T,recursive = T)), 
                        pix[pix$dates == dates[i],], fun = mean)) 
    # And lagged drought code (3, 6 and 12 months before)
    dc.3.lst[[i]] <- as.data.frame(extract(raster(list.files(path = dc.wd, pattern = month3[i], full.names = T,recursive = T)), 
                        pix[pix$month3 == month3[i],], fun = mean))
    dc.6.lst[[i]] <- as.data.frame(extract(raster(list.files(path = dc.wd, pattern = month6[i],full.names = T,recursive = T)), 
                        pix[pix$month6 == month6[i],],  fun = mean))
    dc.12.lst[[i]] <- as.data.frame(extract(raster(list.files(path = dc.wd, pattern = month12[i], full.name = T,recursive  = T)), 
                        pix[pix$month12 == month12[i],], fun = mean))
}

# Merge dataframes in list    
dc.df <- do.call("rbind", dc.lst)
dc.3.df <- do.call("rbind", dc.3.lst)
dc.6.df <- do.call("rbind", dc.6.lst)
dc.12.df <- do.call("rbind", dc.12.lst)

# Change column names
colnames(dc.df) <- "dc"
colnames(dc.3.df) <- "3_mnth_dc"
colnames(dc.6.df) <- "6_mnth_dc"
colnames(dc.12.df) <- "12_mnth_dc"

# KBDI ####

# Name the main folder where the input variables are stored
kbdi.wd <- "~/BTO projects/Polesia wildfires/input variables/dataset-cems-fire-historical-KBDI"

# Create empty lists
kbdi.lst <- list() 
kbdi.3.lst <- list() 
kbdi.6.lst <- list() 
kbdi.12.lst <- list() 

# Loop through each date reading in rasters calculating the mean value within grid cells
for (i in 1:length(dates)) {
  kbdi.lst[[i]] <- as.data.frame(extract(raster(list.files(path = kbdi.wd, pattern = dates[i], full.names = T,recursive = T)), 
                          pix[pix$dates == dates[i],], fun = mean)) 
  # And lagged drought code (3, 6 and 12 months before)
  kbdi.3.lst[[i]] <- as.data.frame(extract(raster(list.files(path = kbdi.wd, pattern = month3[i], full.names = T,recursive = T)), 
                          pix[pix$month3 == month3[i],], fun = mean))
  kbdi.6.lst[[i]] <- as.data.frame(extract(raster(list.files(path = kbdi.wd, pattern = month6[i], full.names = T,recursive = T)), 
                          pix[pix$month6 == month6[i],], fun = mean))
  kbdi.12.lst[[i]] <- as.data.frame(extract(raster(list.files(path = kbdi.wd, pattern = month12[i], full.names = T,recursive = T)), 
                          pix[pix$month12 == month12[i],], fun = mean))
}

# Merge dataframes in list    
kbdi.df <- do.call("rbind", kbdi.lst)
kbdi.3.df <- do.call("rbind", kbdi.3.lst)
kbdi.6.df <- do.call("rbind", kbdi.6.lst)
kbdi.12.df <- do.call("rbind", kbdi.12.lst)

# Change column names
colnames(kbdi.df) <- "kbdi"
colnames(kbdi.3.df) <- "3_mnth_kbdi"
colnames(kbdi.6.df) <- "6_mnth_kbdi"
colnames(kbdi.12.df) <- "12_mnth_kbdi"

# FINE FUEL MOISTURE CONTENT ##########################

# Name the main folder where the input variables are stored
ffmc.wd <- "~/BTO projects/Polesia wildfires/input variables/dataset-cems-fire-historical-FFMC"

# Create empty list to fill with FFMC data
ffmc.lst <- list()

# Loop through each date reading in rasters calculating the mean value within grid cells
for (i in 1:length(dates)) {
  ffmc.lst[[i]] <- as.data.frame(extract(raster(list.files(path = ffmc.wd, pattern = dates[i], full.names = T,recursive = T)), 
                                        pix[pix$dates == dates[i],], fun = mean)) 
} ## No lagged data for moisture content

# Merge dataframes in list    
ffmc.df <- do.call("rbind", ffmc.lst)
# Change column name
colnames(ffmc.df) <- "ffmc"

# FWI ##########################

# Name the main folder where the input variables are stored
fwi.wd <- "~/BTO projects/Polesia wildfires/input variables/dataset-cems-fire-historical-FWI"

# Create empty list 
fwi.lst <- list()

# Loop through each date reading in all rasters and calculating the mean value within grid cells
for (i in 1:length(dates)) {
  fwi.lst[[i]] <- as.data.frame(extract(raster(list.files(path = fwi.wd, pattern = dates[i], full.names = T,recursive = T)), 
                                         pix[pix$dates == dates[i],], fun = mean)) 
} ## No lagged data for fire weather index

# Merge dataframes in list    
fwi.df <- do.call("rbind", fwi.lst)

# Change column name
colnames(fwi.df) <- "fwi"

#### Merge all dataframes ####
dfs <- sapply(.GlobalEnv, is.data.frame) 
df <- do.call(cbind, mget(names(dfs)[dfs]))
# Rename columns
colnames(df)[11:12] <- c("z", "pix")
check <- select(df, -contains("pix.")) # Remove dates and geometry
# Export as csv file
setwd("~/BTO projects/Polesia wildfires/input variables") # Change working directory to where you are storing your output files
write.csv(df, "FWI.csv", row.names = F) ## Outputs are dataframe of fire weather indices and grid cell IDs
