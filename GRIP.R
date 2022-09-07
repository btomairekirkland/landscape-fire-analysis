# Extracting distance to nearest road and road density for grid cells and summarising by fire patch, 
# using data from https://www.globio.info/download-grip-dataset

# Clear data
rm(list = ls())

# Load libraries 
library(sf)
library(raster)
library(dplyr)

# Read in fire shapefiles
setwd("~/BTO projects/Polesia wildfires/fire shapefiles/pixels")
pix <- st_read("all_pixels.shp") 

# Read in road shapefiels from Asia and Europe (continental shapefiles rather than global to speed up processing) and global road density raster
roads.europe <- st_read("~/BTO projects/Polesia wildfires/input variables/GRIP4_Region4_vector_shp/GRIP4_region4.shp")
roads.asia <- st_read("~/BTO projects/Polesia wildfires/input variables/GRIP4_Region5_vector_shp/GRIP4_region5.shp")
dens <- raster("~/BTO projects/Polesia wildfires/input variables/GRIP_density/grip4_total_dens_m_km2.asc")

# Combine both road shapefiles into one
roads <- rbind(roads.asia, roads.europe)

# Find nearest road to grid cell
n.ptch <- st_nearest_feature(pix, roads)
# Create empty vector
distances <- 0

# Calculate distances between grid cells and nearest road (in m)
for (i in 1:nrow(pix)) {
  distances[i] <- st_distance(pix[i,], roads[n.ptch[i],])
}

# Add distances to dataset and convert to km
pix$dist_road <- df$distances/1000 

# Extract road density 
pix$dens_road <- extract(dens, pix, fun = mean, weights = T)

# Drop geometry and export as csv
df <- st_drop_geometry(pix)
write.csv(df, "grid-cells-roads.csv", row.names = F)

# Summarise for each fire patch calculating minimum distance
fires <- df[!grepl("^CP", df$z),] # Select burnt grid cells
fires <- fires %>% 
  group_by(z) %>% 
  summarise(dist_road  = min(dist_road), 
            dens_road = mean(dens_road))

# Export as csv
write.csv(fires, "fire-patches-roads.csv", row.names = F)
