# Extract covariates within grid cells using datasets from the GEE catalogue
# To use the ouput of this file for subsequent analysis you will need to download the files from your Google Drive

# Clear data
rm(list = ls())

# Load libraries
library(rgee)
library(lubridate)

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

# Create vector to fill with new dates
new.dts <- fire.dts
# Correspond NDVI date to fire date
for (i in 1:length(ndvi.str.date)) {
  new.dts[as.Date(fire.dts) >= ndvi.str.date[i+1] & as.Date(fire.dts) <= ndvi.end.date[i+2]] <- 
    as.character(ndvi.str.date[i])
}

## Create empty list to fill with NDVI values
ndvi.mean.lst <- list()
ndvi.mx.lst <- list()
ndvi.min.lst <- list()

# Extract NDVI values
for (i in 1:length(new.dts)) {
  # Mean NDVI
  ndvi.mean.lst[[i]] <- 
    ndvi$filter(ee$Filter$date(new.dts[i]))$toBands()$ 
    # Filter NDVI images for 16 to 8 days before fire starts to get closest NDVI image that spans period
    # before fire
    reduceRegions( 
      collection = fires.sp$filter(ee$Filter$eq('mindate', fire.dts[i]))$select('z'),
      reducer = ee$Reducer$mean()) 
  # Repeat for maximum NDVI
  ndvi.mx.lst[[i]] <- 
    ndvi$filter(ee$Filter$date(new.dts[i]))$toBands()$ 
    reduceRegions( 
      collection = fires.sp$filter(ee$Filter$eq('mindate', fire.dts[i]))$select('z'),
      reducer = ee$Reducer$max(),
      scale = 250) 
  # Repeat for minimum NDVI
  ndvi.min.lst[[i]] <- 
    ndvi$filter(ee$Filter$date(new.dts[i]))$toBands()$ 
    reduceRegions( 
      collection = fires.sp$filter(ee$Filter$eq('mindate', fire.dts[i]))$select('z'),
      reducer = ee$Reducer$min(),
      scale = 250) 
}

# Merge list of feature collections
ndvi.mean.df <- ee$FeatureCollection(ndvi.mean.lst)$flatten()
ndvi.min.df <- ee$FeatureCollection(ndvi.min.lst)$flatten()
ndvi.mx.df <- ee$FeatureCollection(ndvi.mx.lst)$flatten()

# Export 
ndvi.mean.exp <- ee$batch$Export$table(ndvi.mean.df, 'fires-mean-250-NDVI', list(driveFolder =  'fire_variables')) 
ndvi.min.exp <- ee$batch$Export$table(ndvi.min.df, 'fires-min-250-NDVI', list(driveFolder =  'fire_variables')) 
ndvi.mx.exp <- ee$batch$Export$table(ndvi.mx.df, 'fires-max-250-NDVI', list(driveFolder =  'fire_variables')) 
## Default format is csv

ndvi.mean.exp$start()
ndvi.min.exp$start()
ndvi.mx.exp$start()

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
soil.mean <- list()
soil.min <- list()
soil.mx <- list()

# Extract soil moisture content
for (i in 1:length(soil.month)) {
  # Mean soil moisture content
  soil.mean[[i]] <- 
    soil$filter(ee$Filter$eq('system:index', gsub("-", "", soil.month[i])))$toBands()$ 
    reduceRegions(
      collection = fires.sp$filter(ee$Filter$eq('ym', soil.month[i]))$select('z'),
      reducer = ee$Reducer$mean()) 
  # Repeat for minimum soil moisture content
  soil.min[[i]] <- 
    soil$filter(ee$Filter$eq('system:index', gsub("-", "", soil.month[i])))$toBands()$ 
    reduceRegions(
      collection = fires.sp$filter(ee$Filter$eq('ym', soil.month[i]))$select('z'),
      reducer = ee$Reducer$min(),
      scale = 250) 
  # Repeat for maximum soil moisture content
  soil.mx[[i]] <- 
    soil$filter(ee$Filter$eq('system:index', gsub("-", "", soil.month[i])))$toBands()$ 
    reduceRegions(
      collection = fires.sp$filter(ee$Filter$eq('ym', soil.month[i]))$select('z'),
      reducer = ee$Reducer$max(),
      scale = 250)
} 

# Merge list of feature collections
soil.mean.df <- ee$FeatureCollection(soil.mean)$flatten()
soil.min.df <- ee$FeatureCollection(soil.min)$flatten()
soil.mx.df <- ee$FeatureCollection(soil.mx)$flatten()

# Export 
soil.mean.exp <- ee$batch$Export$table(soil.mean.df, 'mean-soil-moisture', list(driveFolder = 'patches')) 
soil.min.exp <- ee$batch$Export$table(soil.min.df, 'min-soil-moisture', list(driveFolder = 'patches')) 
soil.mx.exp <- ee$batch$Export$table(soil.mx.df, 'max-soil-moisture', list(driveFolder = 'patches')) 

soil.mean.exp$start()
soil.min.exp$start()
soil.mx.exp$start()

#### GLOBAL FRICTION SURFACE ##########################

fric <- ee$Image('Oxford/MAP/friction_surface_2019')$
  select('friction')

# Mean travel time across fire patch
fric.mean <- fric$reduceRegions(
  collection = fires.sp$select('z'),
  reducer = ee$Reducer$mean(),
  scale = 250)
# Minimum travel time from fire patch
fric.min <- fric$reduceRegions(
  collection = fires.sp$select('z'),
  reducer = ee$Reducer$min(),
  scale = 250)

# Export 
fric.mean.exp <- ee$batch$Export$table(fric.mean, 'mean-friction', list(driveFolder = 'patches')) 
fric.min.exp <- ee$batch$Export$table(fric.min, 'min-friction', list(driveFolder = 'patches')) 

fric.mean.exp$start()
fric.min.exp$start()

#### POPULATION DENSITY ##########################

pop <- ee$ImageCollection("CIESIN/GPWv411/GPW_Population_Density")$
  select('population_density')

# Extract image ID from collection
pop.id <- pop$aggregate_array('system:index')$getInfo()

# Add year to fire properties
fires.sp <- fires.sp$map(function (feature) {
  num = ee$Number$parse(ee$String(feature$get('mindate'))$slice(0,4)) # Selects between the first and 4 elements of a string 
  # and converts to a number
  feature$set('year', num);
})

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

