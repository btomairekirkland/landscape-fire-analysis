# Clear data
rm(list = ls())

# Load libraries
library(rgee)
library(reticulate)
library(lubridate)
library(dplyr)
library(tidyr)

ee_Initialize() ##  Initialize GEE language

# Read in grid cells from GEE assets 
pixels <- ee$FeatureCollection(paste0(ee_get_assethome(), '/all_pixels'))

# View property names 
pixels$first()$propertyNames()$getInfo() ## Each grid cell has a 'z' and a 'pix' ID value 
## All 'pix' values are unique, but grid cells within a fire patch have same 'z' value
# Vector of dates 
dts <- unique(pixels$aggregate_array('dates')$getInfo())

# Extract covariates within grid cells using datasets from the GEE catalogue

#### NDVI ##########################

ndvi.terra <- ee$ImageCollection("MODIS/061/MOD13Q1")$
  select('NDVI')
ndvi.aqua <- ee$ImageCollection("MODIS/061/MYD13Q1")$
  select('NDVI')

# Merge feature collections
ndvi <- ndvi.terra$merge(ndvi.aqua)

# Get start dates of NDVI composite image
ndvi.str.date <- as.Date(gsub("_", "-", substring(ndvi$aggregate_array('system:index')$getInfo() ,3)))
ndvi.end.date <- as.Date(ndvi.str.date) %m-% days(1)

# Create empty vector to fill with new dates
new.dts <- dts
# Match NDVI date to fire date so that composite date does not include fire or post fire period
for (i in 1:length(ndvi.str.date)) {
  new.dts[as.Date(dts) >= ndvi.str.date[i+1] & as.Date(dts) <= ndvi.end.date[i+2]] <- 
    as.character(ndvi.str.date[i])
}

# Create empty list to fill with NDVI values
ndvi.lst <- list()

# Extract NDVI values
for (i in 1:length(new.dts)) {
  ndvi.lst[[i]] <- 
    pixels$filter(ee$Filter$eq('dates', dts[i]))$map(function (feature) { 
      # Map over each pixel with given date as computationally intensive
      ee$Feature(feature$geometry(), ndvi$filter(ee$Filter$date(new.dts[i]))$toBands()$rename('NDVI')$
                   reduceRegion(
        # Filter NDVI images for 16 to 8 days before fire starts to get closest NDVI image that spans period
        # before fire
        reducer = ee$Reducer$mean(),
        geometry = feature$geometry(),
        scale = 250
      ))$
        set('z', feature$get('z'))$
        set('pix', feature$get('pix'))
    })
}

# Merge features and export
ndvi.df <- ee$FeatureCollection(ndvi.lst)$flatten()
ndvi.exp <- ee$batch$Export$table(ndvi.df, 'NDVI-250', list(driveFolder = 'pixels', # Save to Google drive
                                                        selectors = c('z', 'pix', 'NDVI'))) # Select columns to export
ndvi.exp$start()
ee_monitoring()

#### CLIMATE ##########################

clim <- ee$ImageCollection('ECMWF/ERA5/DAILY')

# Temperature ~~~~~~~~~~~~~~~~~~~~~~~~~
mx.tmp <- clim$select('maximum_2m_air_temperature')
mean.tmp <- clim$select('mean_2m_air_temperature')

## Create empty lists
mx.tmp.lst.1 <- list()
mx.tmp.lst.2 <- list()
mean.tmp.lst.1 <- list()
mean.tmp.lst.2 <- list()

# Extract temperature values on the day the fire started
# Maximum daily temperature
# Split grid cells in half and do half first as computationally intensive
for (i in 1:(length(dts)/2)) {
  mx.tmp.lst.1[[i]] <- pixels$filter(ee$Filter$eq('dates', dts[i]))$map(function (feature) { 
    ee$Feature(feature$geometry(), mx.tmp$filter(ee$Filter$date(dts[i]))$toBands()$
                 rename("max_temp")$
                 reduceRegion(
                   reducer = ee$Reducer$mean(),
                   geometry = feature$geometry(),
                   scale = 250
                 ))$
      set('z', feature$get('z'))$
      set('pix', feature$get('pix'))
  })
} 

# Repeat for second half
for (i in ((length(dts)/2)+1):length(dts)) {
  mx.tmp.lst.2[[i-(length(dts)/2)]] <- pixels$filter(ee$Filter$eq('dates', dts[i]))$map(function (feature) { 
    ee$Feature(feature$geometry(), mx.tmp$filter(ee$Filter$date(dts[i]))$toBands()$
                 rename("max_temp")$
                 reduceRegion(
                   reducer = ee$Reducer$mean(),
                   geometry = feature$geometry(),
                   scale = 250
                 ))$
      set('z', feature$get('z'))$
      set('pix', feature$get('pix'))
  })
} 

# Repeat for mean temperature
for (i in 1:(length(dts)/2)) {
  mean.tmp.lst.1[[i]] <- pixels$filter(ee$Filter$eq('dates', dts[i]))$map(function (feature) { 
    ee$Feature(feature$geometry(), mean.tmp$filter(ee$Filter$date(dts[i]))$toBands()$
                 rename("mean_temp")$
                 reduceRegion(
                   reducer = ee$Reducer$mean(),
                   geometry = feature$geometry(),
                   scale = 250
                 ))$
      set('z', feature$get('z'))$
      set('pix', feature$get('pix'))
  })
} 

for (i in ((length(dts)/2)+1):length(dts)) {
  mean.tmp.lst.2[[i-(length(dts)/2)]] <- pixels$filter(ee$Filter$eq('dates', dts[i]))$map(function (feature) { 
    ee$Feature(feature$geometry(), mean.tmp$filter(ee$Filter$date(dts[i]))$toBands()$
                 rename("mean_temp")$
                 reduceRegion(
                   reducer = ee$Reducer$mean(),
                   geometry = feature$geometry(),
                   scale = 250
                 ))$
      set('z', feature$get('z'))$
      set('pix', feature$get('pix'))
  })
} 

# Merge list of feature collections
mx.tmp.df.1 <- ee$FeatureCollection(mx.tmp.lst.1)$flatten()
mx.tmp.df.2 <- ee$FeatureCollection(mx.tmp.lst.2)$flatten()
mean.tmp.df.1 <- ee$FeatureCollection(mean.tmp.lst.1)$flatten()
mean.tmp.df.2 <- ee$FeatureCollection(mean.tmp.lst.2)$flatten()

# Export 
mx.tmp.exp.1 <- ee$batch$Export$table(mx.tmp.df.1, 'max-temp-1', list(driveFolder = 'pixels', 
                                                                      selectors = c("z", "pix", "max_temp"))) 
mx.tmp.exp.2 <- ee$batch$Export$table(mx.tmp.df.2, 'max-temp-2', list(driveFolder = 'pixels', 
                                                                      selectors = c("z", "pix", "max_temp"))) 
mean.tmp.exp.1 <- ee$batch$Export$table(mean.tmp.df.1, 'mean-temp-1', list(driveFolder = 'pixels', 
                                                                           selectors = c("z", "pix", "mean_temp")))
mean.tmp.exp.2 <- ee$batch$Export$table(mean.tmp.df.2, 'mean-temp-2', list(driveFolder = 'pixels', 
                                                                           selectors = c("z", "pix", "mean_temp")))

mx.tmp.exp.1$start()
mx.tmp.exp.2$start()
mean.tmp.exp.1$start()
mean.tmp.exp.2$start()

# Precipitation ~~~~~~~~~~~~~~~~~~~~~~~

rain <- clim$select('total_precipitation')

# Create empty lists
rain.lst.1 <- list()
rain.lst.2 <- list()

# Extract total rainfall on the day of the fire
# Again, do this in two halves
for (i in 1:(length(dts)/2)) {
  rain.lst.1[[i]] <- pixels$filter(ee$Filter$eq('dates', dts[i]))$map(function (feature) { 
    ee$Feature(feature$geometry(), rain$filter(ee$Filter$date(dts[i]))$toBands()$
                 rename("rain")$
                 reduceRegion(
                   reducer = ee$Reducer$mean(),
                   geometry = feature$geometry(),
                   scale = 250
                 ))$
      set('z', feature$get('z'))$
      set('pix', feature$get('pix'))
  })
} 

# Merge list of feature collections and export
rain.df.1 <- ee$FeatureCollection(rain.lst.1)$flatten()
rain.exp.1 <- ee$batch$Export$table(rain.df.1, 'rainfall-1', list(driveFolder = 'pixels', selectors = c("z", "pix", "rain"))) 
rain.exp.1$start()

# Repeat for second half of grid cells
for (i in ((length(dts)/2)+1):length(dts)) {
  rain.lst.2[[i-(length(dts)/2)]] <- pixels$filter(ee$Filter$eq('dates', dts[i]))$map(function (feature) { 
    # Map over each pixel with given date
    ee$Feature(feature$geometry(), rain$filter(ee$Filter$date(dts[i]))$toBands()$
                 rename("rain")$
                 reduceRegion(
                   reducer = ee$Reducer$mean(),
                   geometry = feature$geometry(),
                   scale = 250
                 ))$
      set('z', feature$get('z'))$
      set('pix', feature$get('pix'))
  })
} 

rain.df.2 <- ee$FeatureCollection(rain.lst.2)$flatten()
rain.exp.2 <- ee$batch$Export$table(rain.df.2, 'rainfall-2', list(driveFolder = 'pixels', selectors = c("z", "pix", "rain"))) 
rain.exp.2$start()

#### SOIL ##########################

soil <- ee$ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')$
  select('soil')

# Create property in grid cell dataset that describes the month for which to collect data 
## That month if it is late in the month, or previous month if it is early in the month
pixels <- pixels$map(function (feature) {
  day = ee$Number$parse(ee$String(feature$get('dates'))$slice(8,10)) # Create property with day ofthe month
  feature$set('day', day)
})

late <- pixels$filter(ee$Filter$gt('day', 15))$map(function (feature) {
  ym = ee$String(ee$Date$parse('YYYY-MM-dd', feature$get('dates'))$format('YYYY-MM'))
  feature$set('ym', ym)
})

early <- pixels$filter(ee$Filter$lte('day', 15))$map(function (feature) {
  ym = ee$String(ee$Date$parse('YYYY-MM-dd', feature$get('dates'))$advance(-1, 'month')$format('YYYY-MM'))
  feature$set('ym', ym)
})

# Merge into feature collection
pixels <- ee$FeatureCollection(ee$List(c(early,late)))$flatten()

# Vector of month of or preceding fire, depending on whether the fire is at the start or end of the month, and year
month <- unique(pixels$aggregate_array('ym')$getInfo())

# Create empty list to fill with soil moisture date
soil.lst <- list()

# Extract soil moisture content
for (i in 1:length(month)) {
  soil.lst[[i]] <- 
   pixels$filter(ee$Filter$eq('ym', month[i]))$map(function (feature) { 
      # Map over each pixel with given date
      ee$Feature(feature$geometry(), soil$filter(ee$Filter$eq('system:index', gsub("-", "", month[i])))$toBands()$
                   rename("soil_moist")$
                   reduceRegion(
                     reducer = ee$Reducer$mean(),
                     geometry = feature$geometry(),
                     scale = 250
                   ))$
        set('z', feature$get('z'))$
        set('pix', feature$get('pix'))
    })
}

# Merge list of feature collections and export
soil.df <- ee$FeatureCollection(soil.lst)$flatten()
soil.exp <- ee$batch$Export$table(soil.df, 'soil-moisture', list(driveFolder = 'pixels', 
                                                                            selectors = c('z', 'pix', 'soil_moist'))) 
soil.exp$start()
ee_monitoring()

#### SNOW COVER ##########################

snow.aqua <- ee$ImageCollection("MODIS/006/MYD10A1")$
  select('NDSI_Snow_Cover')
snow.terra <- ee$ImageCollection("MODIS/006/MOD10A1")$
  select('NDSI_Snow_Cover')

# Merge feature collections
snow <- snow.terra$merge(snow.aqua)

# Create empty list to fill with maximum, minimum, and mean snow cover during the winter period
snow.mean <- list()
snow.mx <- list()
snow.min <- list()

# Calculate snow date
# First add month variable to dataset to identify winter and non-winter pixels
pixels <- pixels$map(function (feature) {
  mnth = ee$Number$parse(ee$String(feature$get('dates'))$slice(5,7)) 
  sd = ee$String(ee$Date$parse('YYYY-MM-dd', feature$get('dates'))$advance(1, 'day')$format('YYYY-MM-dd'))
  feature$set('month', mnth)$
    set('snow_dt', sd) # Add snow date as day after fire date because the last date when filtering images is not included
})

# Get vector for start and end date to calculate snow cover
months <- as.numeric(substr(dts, 6, 7)) # First extract month - is the fire in December?
year <- as.numeric(substr(dts, 1, 4)) # Extract year to paste into snow date 
snow.date <- 0 # Create vector of dates up to which to calculate snow cover (1st of March if outside of winter)
str.date <- 0 # Vector of start dates to calculate snow cover

for (i in 1:length(dts)) {
  snow.date[i] <- ifelse(months[i] != 12 & months[i] != 1 & months[i] != 2, paste0(year[i],"-03-01"), 
                         as.character(as.Date(dts[i]) + 1))
  str.date[i] <- ifelse(months[i] == 12, paste0(year[i],"-12-01"), paste0(year[i]-1,"-12-01"))
}

# Simplify vectors retaining unique dates to speed up processing
str.date <- str.date[!duplicated(snow.date)]
snow.date <- snow.date[!duplicated(snow.date)]

# Change snow date for fires outside of winter to last day of winter
non.winter <- pixels$filter(ee$Filter$rangeContains('month', 3, 11))$
  map(function (feature) {
  feature$set('snow_dt', ee$String(feature$get('year'))$cat('-03-01'))
})
# Select winter pixels
winter <- pixels$filter(ee$Filter$inList('month', ee$List(c(1,2,12))))

# Merge into feature collection
pixels <- ee$FeatureCollection(ee$List(c(winter,non.winter)))$flatten()

for (i in 1) {
  snow.mean[[i]] <- pixels$filter(ee$Filter$eq('snow_dt', snow.date[i]))$map(function (feature) { 
  ee$Feature(feature$geometry(), snow$filter(ee$Filter$date(str.date[i], snow.date[i]))$mean()
             # Merge images to calculate mean snow cover over the winter period 
             $rename('mean_snow')$
               reduceRegion(
                 reducer = ee$Reducer$mean(),
                 geometry = feature$geometry(),
                 scale = 250
               ))$
      set('z', feature$get('z'))$
      set('pix', feature$get('pix'))
    }) 
  
  snow.mx[[i]] <- pixels$filter(ee$Filter$eq('snow_dt', snow.date[i]))$map(function (feature) { 
    ee$Feature(feature$geometry(), snow$filter(ee$Filter$date(str.date[i], snow.date[i]))$max()
               # Repeat for maximum snow cover
               $rename('max_snow')$
                 reduceRegion(
                   reducer = ee$Reducer$mean(),
                   geometry = feature$geometry(),
                   scale = 250
                 ))$
      set('z', feature$get('z'))$
      set('pix', feature$get('pix'))
    })

  snow.min[[i]] <- pixels$filter(ee$Filter$eq('snow_dt', snow.date[i]))$map(function (feature) { 
    ee$Feature(feature$geometry(), snow$filter(ee$Filter$date(str.date[i], snow.date[i]))$min()
               # Repeat for minimum snow cover
               $rename('min_snow')$
                 reduceRegion(
                   reducer = ee$Reducer$mean(),
                   geometry = feature$geometry(),
                   scale = 250
                 ))$
      set('z', feature$get('z'))$
      set('pix', feature$get('pix'))
    })
}

# Merge list of feature collections
snow.mx.df <- ee$FeatureCollection(snow.mx)$flatten()
snow.min.df <- ee$FeatureCollection(snow.min)$flatten()
snow.mean.df <- ee$FeatureCollection(snow.mean)$flatten()

# Export
snow.mean.exp <- ee$batch$Export$table(snow.mean.df, 'mean-snow-cover', list(driveFolder = 'pixels', 
                                                                             selectors = c('z', 'pix', 'mean_snow'))) 
snow.mx.exp <- ee$batch$Export$table(snow.mx.df, 'max-snow-cover', list(driveFolder = 'pixels', 
                                                                        selectors = c('z', 'pix', 'max_snow'))) 
snow.min.exp <- ee$batch$Export$table(snow.min.df, 'min-snow-cover', list(driveFolder = 'pixels', 
                                                                          selectors = c('z', 'pix', 'min_snow'))) 
snow.mean.exp$start()
snow.mx.exp$start()
snow.min.exp$start()
ee_monitoring()

#### GLOBAL FRICTION SURFACE ##########################

## Static dataset so no need to match dates

fric <- ee$Image('Oxford/MAP/friction_surface_2019')$
  select('friction')

# Mean travel time within pixels
fric.df <- pixels$map(function (feature) { 
  # Map over each pixel with given date
  ee$Feature(feature$geometry(), fric$
               reduceRegion(
                 reducer = ee$Reducer$mean(),
                 geometry = feature$geometry()
               ))$set('z', feature$get('z'))$
                  set('pix', feature$get('pix'))
})

# Export 
fric.exp <- ee$batch$Export$table(fric.df, 'friction', list(driveFolder = 'pixels', selectors = c('z', 'pix', 'friction'))) 
fric.exp$start()
ee_monitoring()

#### POPULATION DENSITY ##########################

pop <- ee$ImageCollection("CIESIN/GPWv411/GPW_Population_Density")$
  select('population_density')

# Extract image ID from collection
id <- pop$aggregate_array('system:index')$getInfo()

# Add year to pixels properties
pixels <- pixels$map(function (feature) {
  num = ee$Number$parse(ee$String(feature$get('dates'))$slice(0,4)) # Selects between the first and 4 elements of a string 
  # and converts to a number
  feature$set('year', num);
})

# Available years are every 5 years between 2000-2020
pop.yrs <- seq(2000,2020, 5)

# Divide pixels to match the 5 yearly intervals
first <- pixels$filter(ee$Filter$rangeContains('year', 2001, 2002))
second <- pixels$filter(ee$Filter$rangeContains('year', 2003, 2007))
third <- pixels$filter(ee$Filter$rangeContains('year', 2008, 2012))
fourth <- pixels$filter(ee$Filter$rangeContains('year', 2013, 2017))
fifth <- pixels$filter(ee$Filter$rangeContains('year', 2018, 2019))

# Combine in list
pix.lst <- list(first, second, third, fourth, fifth)

# Create empty list 
pop.lst <- list()

# Extract population density within grid cells
for (i in 1:length(pop.yrs)) {
  pop.lst[[i]] <- pix.lst[[i]]$map(function (feature) { 
  ee$Feature(feature$geometry(), pop$filter(ee$Filter$eq('system:index', id[i]))$toBands()
             $rename('pop_dens')$
               reduceRegion(
                 reducer = ee$Reducer$mean(),
                 geometry = feature$geometry(),
                 scale = 250
               ))$
    set('z', feature$get('z'))$
    set('pix', feature$get('pix'))
  })
}

# Flatten feature collections and export
pop.df <- ee$FeatureCollection(pop.lst)$flatten()
pop.exp <- ee$batch$Export$table(pop.df, 'pop', list(driveFolder = 'pixels', selectors = c('z', 'pix', 'pop_dens'))) 
pop.exp$start()
ee_monitoring()
