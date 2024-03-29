# Extracting distance to nearest road and road density for grid cells and summarising by fire patch, 
# using data from https://www.globio.info/download-grip-dataset

# Clear data
rm(list = ls())

# Load libraries 
library(sf)
library(raster)
library(dplyr)

# Read in shapefiles
setwd("~/BTO projects/Polesia wildfires/fire shapefiles/pixels")
shp <- st_read("all_pixels.shp") ## Voxel dataset
## Can be done instead for ignition points
## setwd("~/BTO projects/Polesia wildfires/fire shapefiles/fire shapefiles")
## shp <- st_read("fire_patches.shp")
## shp <- st_as_sf(shp, coords = c("I_LON", "I_LAT")) # Select ignition coordinates

# Read in road shapefiles from Asia and Europe (continental shapefiles rather than global to speed up processing) and global road density raster
roads.europe <- st_read("~/BTO projects/Polesia wildfires/input variables/GRIP4_Region4_vector_shp/GRIP4_region4.shp")
roads.asia <- st_read("~/BTO projects/Polesia wildfires/input variables/GRIP4_Region5_vector_shp/GRIP4_region5.shp")
dens <- raster("~/BTO projects/Polesia wildfires/input variables/GRIP_density/grip4_total_dens_m_km2.asc")

# Combine both road shapefiles into one
roads <- rbind(roads.asia, roads.europe)

# Find nearest road to grid cell (or ignition point)
n.ptch <- st_nearest_feature(shp, roads)
# Create empty vector
distances <- 0

# Calculate distances to nearest road (in m)
for (i in 1:nrow(shp)) {
  distances[i] <- st_distance(shp[i,], roads[n.ptch[i],])
}

# Add distances to dataset and convert to km
shp$dist_road <- df$distances/1000 

# Extract road density 
shp$dens_road <- extract(dens, shp, fun = mean, weights = T)

# Drop geometry and export as csv
df <- st_drop_geometry(shp)
setwd("~/BTO projects/Polesia wildfires/input files") # Change working directory to where you are storing output files
write.csv(df, "grid-cells-roads.csv", row.names = F) # Output file is a dataframe with grid cell IDs, distance to nearest road and road density

# Can summarise for each fire patch by calculating minimum distance if using voxel dataset
fires <- df[!grepl("^CP", df$z),] # Select burnt grid cells
fires <- fires %>% 
  group_by(z) %>% 
  summarise(dist_road  = min(dist_road), 
            dens_road = mean(dens_road))

# Export as csv
write.csv(fires, "fire-patches-roads.csv", row.names = F) # Output file is same as above but grouped by fires with only fire patch ID
