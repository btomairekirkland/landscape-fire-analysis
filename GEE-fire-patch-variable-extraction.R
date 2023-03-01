# Extract covariates within fire patches using datasets from the GEE catalogue
# To use the ouput of this file for subsequent analysis you will need to download the files from your Google Drive

# Clear data
rm(list = ls())

# Load libraries
library(rgee)
library(lubridate)
library(date.table)

ee_Initialize() ##  Initialize GEE language

# Read in fire shapefiles from GEE assets
fires.sp <- ee$FeatureCollection(paste0(ee_get_assethome(), '/fire_patches'))

# Get info from fires 
fires.sp$first()$propertyNames()$getInfo()
fire.dts <- unique(fires.sp$aggregate_array('mindate')$getInfo()) ## Unique dates of fire
strt.dts <- fires.sp$aggregate_array('mindate')$getInfo() ## Date of start of fire
end.dts <- as.character(as.Date(fires.sp$aggregate_array('maxdate')$getInfo()) %m+% days(1)) ## Day after fire
id <- fires.sp$aggregate_array('z')$getInfo() ## Fire ID

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

# Order by date
ndvi.str.date <- ndvi.str.date[order(ndvi.str.date)]
ndvi.end.date <- ndvi.str.date[order(ndvi.end.date)]

# Find closest proceeding NDVI image was not taken during or after the fire
# Create vector of unique end dates for fires
fire.tab = data.table(as.Date(fire.dts)) # What was the last day of the fire? 
colnames(fire.tab) <- "date" # Change column name
ndvi.dt = data.table(date = as.Date(ndvi.end.date), key = 'date') # What is the last date of the NDVI composite image?

# For each fire date, join the closest value end date for the NDVI image while rolling to infinity
indx <- ndvi.dt[fire.tab, roll = Inf, which = T]

# 'indx' is the location in the vector of the NDVI end dates at are closest to the end of the fire date
# Extract the start date of the previous NDVI image 
new.dts <- as.character(ndvi.str.date[indx-1])

# How close are these images to the fires?
diff <- difftime(as.Date(ndvi.dt), as.Date(fire.dts), units = "days")
min(diff);max(diff);mean(diff)

## Create empty list to fill with NDVI values
ndvi.lst <- list()

# Extract NDVI values
for (i in 1:length(new.dts)) {
  ndvi.lst[[i]] <- 
    ndvi$filter(ee$Filter$date(new.dts[i]))$toBands()$ 
    # Filter NDVI images for 16 to 8 days before fire starts to get closest NDVI image that spans period before fire
    reduceRegions( 
      collection = fires.sp$filter(ee$Filter$eq('mindate', fire.dts[i]))$select('z'),
      reducer = ee$Reducer$mean()$combine(
        ee$Reducer$max(), sharedInputs = TRUE)$combine(
          ee$Reducer$min(), sharedInputs = TRUE),
      scale = 250
    ) 
}

# Merge list of feature collections
ndvi.df <- ee$FeatureCollection(ndvi.lst)$flatten()

# Export 
ndvi.exp <- ee$batch$Export$table(ndvi.df, '250-NDVI', list(driveFolder =  'fire_variables')) 
## Default format is csv
ndvi.exp$start()

#### CLIMATE ##########################

clim <- ee$ImageCollection('ECMWF/ERA5/DAILY')

# Temperature ~~~~~~~~~~~~~~~~~~~~~~~~~
mx.tmp <- clim$select('maximum_2m_air_temperature')
mean.tmp <- clim$select('mean_2m_air_temperature')

## Create empty lists
mx.tmp.lst <- list()
mean.tmp.lst <- list()

# Extract temperature values over lifetime of fire
# Loop through each fire individually
for (i in 1:length(strt.dts)) {
  # Maximum temperature
  mx.tmp.lst[[i]] <- 
    mx.tmp$filter(ee$Filter$date(strt.dts[i], end.dts[i]))$max()$ 
    reduceRegions(
      collection = fires.sp$filter(ee$Filter$eq('z', id[i]))$select('z'),
      reducer = ee$Reducer$max(),
      scale = 250)
  # Repeat for mean temperature
  mean.tmp.lst[[i]] <-
    mean.tmp$filter(ee$Filter$date(strt.dts[i], end.dts[i]))$mean()$ 
    reduceRegions(
      collection = fires.sp$filter(ee$Filter$eq('z', id[i]))$select('z'),
      reducer = ee$Reducer$mean(),
      scale = 250)
} 

# Merge list of feature collections
mx.tmp.df <- ee$FeatureCollection(mx.tmp.lst)$flatten()
mean.tmp.df <- ee$FeatureCollection(mean.tmp.lst)$flatten()

# Export 
mx.tmp.exp <- ee$batch$Export$table(mx.tmp.df, 'max-temp', list(driveFolder = 'patches')) 
mean.tmp.exp <- ee$batch$Export$table(mean.tmp.df, 'mean-temp', list(driveFolder = 'patches'))

mx.tmp.exp$start()
mean.tmp.exp$start()

# Precipitation ~~~~~~~~~~~~~~~~~~~~~~~

rain <- clim$select('total_precipitation')

# Create empty lists
rain.mean <- list()
rain.min <- list()
rain.end <- list()
rain.max <- list()

# Extract rainfall during the lifetime of the fire and the day after 
for (i in 1:length(strt.dts)) {
  # Maximum rainfall
  rain.max[[i]] <- 
    rain$filter(ee$Filter$date(strt.dts[i], end.dts[i]))$max()$ 
    reduceRegions(
      collection = fires.sp$filter(ee$Filter$eq('z', id[i]))$select('z'),
      reducer = ee$Reducer$max(),
      scale = 250)
  # Repeat for minimum rainfall
  rain.min[[i]] <- 
    rain$filter(ee$Filter$date(strt.dts[i], end.dts[i]))$min()$ 
    reduceRegions(
      collection = fires.sp$filter(ee$Filter$eq('z', id[i]))$select('z'),
      reducer = ee$Reducer$min(),
      scale = 250)
  # Repeat for mean rainfall
  rain.mean[[i]] <- 
    rain$filter(ee$Filter$date(strt.dts[i], end.dts[i]))$mean()$ 
    reduceRegions(
      collection = fires.sp$filter(ee$Filter$eq('z', id[i]))$select('z'),
      reducer = ee$Reducer$mean(),
      scale = 250)
  # Total rainfall the day after fire, which might supress fires?
  rain.end[[i]] <- 
    rain$filter(ee$Filter$date(end.dts[i]))$toBands()$ 
    reduceRegions(
      collection = fires.sp$filter(ee$Filter$eq('z', id[i]))$select('z'),
      reducer = ee$Reducer$mean(),
      scale = 250)
} 

# Merge list of feature collections
rain.min.df <- ee$FeatureCollection(rain.min)$flatten()
rain.mean.df <- ee$FeatureCollection(rain.mean)$flatten()
rain.end.df <- ee$FeatureCollection(rain.end)$flatten()
rain.max.df <- ee$FeatureCollection(rain.max)$flatten()

# Export 
rain.min.exp <- ee$batch$Export$table(rain.min.df, 'min-rainfall', list(driveFolder = 'patches')) 
rain.mean.exp <- ee$batch$Export$table(rain.mean.df, 'mean-rainfall', list(driveFolder = 'patches')) 
rain.end.exp <- ee$batch$Export$table(rain.end.df, 'end-rainfall', list(driveFolder = 'patches')) 
rain.max.exp <- ee$batch$Export$table(rain.max.df, 'max-rainfall', list(driveFolder = 'patches')) 

rain.max.exp$start()
rain.min.exp$start()
rain.end.exp$start()
rain.mean.exp$start()

#### SOIL ##########################

soil <- ee$ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')$
  select('soil')

# Create property in fire dataset that describes the month for which to collect data 
## That month if it is late in the month, or previous month if it is early in the month
fires.sp <- fires.sp$map(function (feature) {
  day = ee$Number$parse(ee$String(feature$get('mindate'))$slice(8,10)) # Create property with day of the month based on fire start date
  feature$set('day', day)
})

late <- fires.sp$filter(ee$Filter$gt('day', 15))$map(function (feature) {
  ym = ee$String(ee$Date$parse('YYYY-MM-dd', feature$get('mindate'))$format('YYYY-MM'))
  feature$set('ym', ym)
})

early <- fires.sp$filter(ee$Filter$lte('day', 15))$map(function (feature) {
  ym = ee$String(ee$Date$parse('YYYY-MM-dd', feature$get('mindate'))$advance(-1, 'month')$format('YYYY-MM'))
  feature$set('ym', ym)
})

# Merge into feature collection
fires.sp <- ee$FeatureCollection(ee$List(c(early,late)))$flatten()

# Vector of month of or preceding fire, depending on whether the fire is at the start or end of the month, and year
soil.month <- unique(fires.sp$aggregate_array('ym')$getInfo())

# Create empty lists
soil.lst <- list()

# Extract soil moisture content
for (i in 1:length(soil.month)) {
  soil.lst[[i]] <- 
    soil$filter(ee$Filter$eq('system:index', gsub("-", "", soil.month[i])))$toBands()$ 
    reduceRegions(
      collection = fires.sp$filter(ee$Filter$eq('ym', soil.month[i]))$select('z'),
      reducer = ee$Reducer$mean()$combine(
        ee$Reducer$max(), sharedInputs = TRUE)$combine(
          ee$Reducer$min(), sharedInputs = TRUE),
      scale = 250)
} 

# Merge list of feature collections
soil.df <- ee$FeatureCollection(soil.lst)$flatten()

# Export 
soil.exp <- ee$batch$Export$table(soil.df, 'soil-moisture', list(driveFolder = 'patches')) 
soil.exp$start()

#### GLOBAL FRICTION SURFACE ##########################

fric <- ee$Image('Oxford/MAP/friction_surface_2019')$
  select('friction')

# Mean and minimum travel time across fire patch
fric.df <- fric$reduceRegions(
  collection = fires.sp$select('z'),
   reducer = ee$Reducer$mean()$combine(
    ee$Reducer$min(), sharedInputs = TRUE),
  scale = 250)

# Export 
fric.exp <- ee$batch$Export$table(fric.df, 'friction', list(driveFolder = 'patches')) 
fric.exp$start()

#### POPULATION DENSITY ##########################

pop <- ee$ImageCollection("CIESIN/GPWv411/GPW_Population_Density")$
  select('population_density')

# Extract image ID from collection
pop.id <- pop$aggregate_array('system:index')$getInfo()

# Available years are every 5 years between 2000-2020
pop.yrs <- seq(2000,2020, 5)

# Divide fires to match the 5 yearly intervals
first <- fires.sp$filter(ee$Filter$rangeContains('year', 2001, 2002))
second <- fires.sp$filter(ee$Filter$rangeContains('year', 2003, 2007))
third <- fires.sp$filter(ee$Filter$rangeContains('year', 2008, 2012))
fourth <- fires.sp$filter(ee$Filter$rangeContains('year', 2013, 2017))
fifth <- fires.sp$filter(ee$Filter$rangeContains('year', 2018, 2019))

# Combine in list
fires.lst <- list(first, second, third, fourth, fifth)

# Create empty list 
pop.lst <- list()

# Extract population density within grid cells
for (i in 1:length(pop.yrs)) {
  pop.lst[[i]] <- pop$filter(ee$Filter$eq('system:index', pop.id[i]))$toBands()$
                 reduceRegions(
                   reducer = ee$Reducer$mean(),
                   collection = fires.lst[[i]]$select('z'), 
                   scale = 250)
}

# Flatten feature collections and export
pop.df <- ee$FeatureCollection(pop.lst)$flatten()
pop.exp <- ee$batch$Export$table(pop.df, 'pop', list(driveFolder = 'patches')) 
pop.exp$start()
