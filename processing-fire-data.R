# Processing FRY data obtained from https://www.nature.com/articles/sdata2018132 and creating non-fire control points

# Clear data
rm(list = ls())

# Load libraries 
library(sf)
library(raster)
library(lubridate)
library(ggplot2)
library(dplyr)

# First read in data, including fire data spreadsheets and shapefiles and land cover data
# Set working directory
setwd("~/BTO projects/Polesia wildfires")
# Read in wildfire trait data obtained from https://www.nature.com/articles/sdata2018132 
fires <- read.csv("fires-clipped.csv")
# Rename patch ID column
colnames(fires)[which(names(fires) == "ptch_id_Oom")] <- "z"
# Read in land cover and extract land cover type at centre point of pixel
lc <- raster("land_cover.tif") ## Land cover map obtained from https://github.com/tpfd/Polesia-Landcover and 

# Obtain study period dates 
start <- min(as.Date(fires$mindt))
end <- max(as.Date(fires$maxdt))

# Name the main folder where the individual shapefiles for fires are stored
dirname <- "~/BTO projects/Polesia wildfires/fire shapefiles/polygons_6D_fixed"
# List files within sub directories
filenames <- list.files(path = dirname, pattern = ".shp", recursive = T)
# Create empty list to fill with polygon data 
list <- list(data.frame(matrix(NA, ncol = 7)))
# Read in polygon data information and save in list
for( i in 1:length(filenames) ){
  list[[i]] <- st_read(paste0(dirname, "/", filenames[i]))
}
# Merge list into single sf object
sf.ob <- do.call(rbind, list)
# Select columns
sf.ob <- subset(sf.ob, select = c(z, mindate, maxdate, geometry))

# Select large fires (>1km2)
large.fires <- fires[fires$area > 1000000,]

# Filter out smaller fires in shapefiles
fires.sp <- sf.ob[sf.ob$z %in% large.fires$z,]

# Extract polygons from multipolygons
fires.sp <- st_collection_extract(fires.sp, "POLYGON")

# Export shapefile of fire patches
setwd("~/BTO projects/Polesia wildfires/fire shapefiles/fire patches")
st_write(large.fires, "fire_patches.shp", append = F)

# Create and export fire ignition points 
ign <- st_as_sf(x = fires, 
                coords = c("i.centre.x", "i.centre.y"),
                crs = NA)
ign <- st_buffer(ign, 0.0022457331/2, endCapStyle = "SQUARE", nQuadSegs = 1) # Add 250m2 buffer
st_crs(ign) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0" # Add crs info 
setwd("~/BTO projects/Polesia wildfires/fire shapefiles/ignition points") # Export
st_write(ign, "ignition_points.shp", append = F) 

# Grid fire data
grid <- st_make_grid(large.fires, cellsize = 0.0022457331) ## 250m2
grid.sf <- st_as_sf(grid)
st_geometry(grid.sf) <- "geometry" # Rneame geometry column
grid.fires <- st_join(x = st_centroid(grid.sf), y = fires, left = F) # Select grid cells whose centre point overlaps 
                                                                  # with a fire
st_crs(grid.fires) <- NA # Remove CRS 
grid.fires <- st_buffer(grid.fires, 0.0022457331/2, endCapStyle = "SQUARE", nQuadSegs = 1) # Add 250m2 buffer
                                                                                           # around centre point
st_crs(grid.fires) <- st_crs(grid.sf) # Specify CRS

# Create control points (non-fire observations)
# Extract land cover type at centre point of control points
grid.sf$lc <- extract(lc, st_centroid(grid.sf))
# Remove points with unclassified (0 & 10) and water (8)
grid.sub <- grid.sf[!(is.na(grid.sf$lc)) & grid.sf$lc != 0 & 
                      grid.sf$lc != 8 & grid.sf$lc != 10,]
# Evenly sample control points across land cover types
grid.eq <- grid.sub %>%
  group_by(lc) %>%
  sample_n(250000, replace = T) 

# Randomly assign date to point between the start date of your study and the end date
grid.eq$dates <- sample(seq(start, end, by = 1), nrow(grid.eq), replace = T)

# Remove duplicates
grid.eq <- grid.eq[!duplicated(grid.eq),]

# Add ID variable to control points
grid.eq$x <- row.names(grid.eq)

# Is the grid cell burning on this date? 
# Get small fires that have no polygon (n cells < 4)
small.fires <- fires[!(fires$z %in% sf.ob$z),]
# Convert dataframe to spatial object
small.fires <- st_as_sf(small.fires, coords = c("CENTRE.LON", "CENTRE.LAT"), crs = NA)
# Add buffer to small fires for each cell that makes up the fire
small.fires <- st_buffer(small.fires, (0.0022457331/2)*small.fires$n.cell, endCapStyle = "SQUARE", nQuadSegs = 1) 
st_crs(small.fires) <- st_crs(grid.eq)

# Combine with other polygons of large fires, ensuring they have the same columns
small.fires <- subset(small.fires, select = c(z, mindt, maxdt, geometry))
colnames(small.fires)[2:3] <- colnames(sf.ob)[2:3]
comb <- rbind(sf.ob, small.fires)

# Vector of dates of control points
dates <- unique(grid.eq$dates)

# Create empty list 
int <- list()
# Intersect points with fires to find points that fall within a burnt area within a year of the fire
for (i in 1:length(dates)) {
  fires.sub <- comb[comb$mindate <= dates[i] & as.Date(comb$maxdate) %m+% months(12) >= dates[i] ,]
  pts <- grid.eq[grid.eq$dates == dates[i],]
  int[[i]] <- pts[st_intersects(pts, fires.sub) %>% lengths > 0,]
}

# Merge dataframes in list    
int.df <- do.call("rbind", int)

# Remove grid cells that overlap in time and space with fires to select control points
ctrl.pts <- grid.eq[!(grid.eq$x %in% int.df$x),]
plot(ctrl.pts[1,], max.plot = 1) # Visually inspect control points

# Explore distribution of control points by land cover type and resample to select 1.5x as many fire observations 
table(ctrl.pts$lc)
ctrl.pts <- ctrl.pts %>%
  group_by(lc) %>%
  sample_n(173415, replace = F)

# Check temporal distribution
temp.trend <- ctrl.pts %>% 
  st_drop_geometry() %>% 
  mutate(month = substr(dates, 1, 7),
         n = 1) %>% 
  group_by(month) %>% 
  summarise(sum = sum(n))

ggplot(temp.trend, aes(x = month, y = sum)) +
  geom_point() +
  ylab("Number of control points") +
  xlab("Month") + 
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(margin = unit(c(0, 2.5, 0, 0), "mm"))) +
  scale_x_discrete(expand = c(.01,.01)) 

# Combine with gridded pixels
grid.fires <- subset(grid.fires, select = c("z", "mindate", "lc", "geometry")) # Select columns
colnames(grid.fires)[2] <- "dates" # Rename start date of fire to match date column in control point dataset
# Add ID to control points
ctrl.pts$z <- paste0("CP", rownames(ctrl.pts))
# Arrange columns in control point dataset
ctrl.pts <- ctrl.pts[,c(5,3,2,1)]
# Convert date column to character
ctrl.pts$dates <- as.character(ctrl.pts$dates)

# Merge datasets
fin.pix <- rbind(grid.fires, ctrl.pts) 
# Add new ID to each grid cell (both fire and non-fire)
fin.pix$pix <- paste0("PX", rownames(fin.pix))

# Create variable for weights to be used in model based on the proportion of land cover types sampled
# What is the distribution of land cover types across all sample grid cells after removing duplicated grid cells?
non.dup <- fin.pix
non.dup[,5:6] <- st_coordinates(st_centroid(non.dup))
non.dup <- non.dup[!duplicated(non.dup[,5:6]),]
weights <- as.data.frame(table(grid.sf$lc)/table(non.dup$lc)); colnames(weights) <- c("lc", "weights")
fin.pix <- merge(fin.pix, weights, by = "lc")

# Add ID to each pixel
fin.pix$pix <- paste0("PX", rownames(fin.pix))

# Reorder columns and remove land cover column
fin.pix <- fin.pix[,c(2,6,3,4,5)]

# Export dated grid cell dataset
setwd("~/BTO projects/Polesia wildfires/fire shapefiles/pixels")
st_write(fin.pix, "all_pixels.shp", append = F)

# Extract centroid coordinates and temporal information of grid cells and save as .csv to merge with covariate data prior to analysis
info <- as.data.frame(st_coordinates(st_centroid(fin.pix)))
info <- cbind(info, st_drop_geometry(fin.pix))
info$mnth <- substr(fin.pix$dates, 6, 7)
info$year <- substr(fin.pix$dates, 1, 4)
info$n_year <- info$year - min(info$year) + 1
info$day <- yday(as.Date(fin.pix$dates))
setwd("~/BTO projects/Polesia wildfires/input files") # Change working directory to where you are storing output files
write.csv(info, "grid-cell-info.csv", row.names = F) 
