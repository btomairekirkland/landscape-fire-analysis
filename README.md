# LANDSCAPE FIRE ANALYSIS WORKFLOW IN GOOGLE EARTH ENGINE AND R

We used a combination of Google Earth Engine and R to perform large-scale extraction, manipulation, and analysis of data on wildfires in Polesia. We use a mixture of data sources, some of which are available through the Earth Engine catalogue. 

The code is written predominantly in R version 4.0.5, with the following libraries required on top of a standard R install:

* rgee
* sf
* raster
* dplyr
* lubridate
* ggplot2
* mgcv

Install rgee: https://r-spatial.github.io/rgee/ 

## Overview of methods

The fire datasets used in this work were provided by Wentao Chan and Florent Mouillot from Le Centre National de la Recherche Scientifique, France. They are part of the global FRY dataset, based on the European Space Agency's FireCCI products, and include polygon shapefiles for fires larger than 4x250m2 pixels and a fire patch trait dataset for all fires (Laurent et al., 2018). We focus on larger fire patches over 1km2 occurring in Polesia since 2001, using a cut-off period for 6 days to group burnt pixels. We used a land cover map of Polesia provided by Mark de Jong of the Canadian Forest Service and Thomas Dowling from the United Nations Environment Programme (UNEP). Using Google Earth Engine (Gorelick et al., 2017; Tamiminia et al., 2020), we simplified the classification from 14 to the following ten land cover types:


Complex class  | Simplified class | Class number
:------------- | :--------------- | :-----------
Pine, birch, wideleafed coniferous forest | Coniferous forest | 1
Spruce forest | Coniferous forest | 1
Birch forest | Coniferous forest | 1
Oak, deciduous forests, small leaved deciduous forest | Deciduous forest | 2
Alder forests | Deciduous forest | 2
Deciduous indigenous swamp forests | Deciduous forest | 2
Meadows | Meadows | 3
Agriculture | Agriculture | 4
Raised bog | Raised bog | 5
Fen and transition mire | Fen and transition mire | 6
Scrub, forest cuttings and clearings, cleared ground outside of urban areas | Scrub | 7
Water| Water | 8
Urban, cleared ground in urban areas, buildings, tarmac | Urban | 9
Unclassified | Unclassified | 10


Example code for reclassifying land cover maps can be found in ‘GEE-land-cover’. 

Using the ‘sf’ 1.0-7 package in R version 4.0.5 (Pebesma et al., 2020), we converted the FRY dataset to a 250x250m grid for each day of the study period (with each cell for each day hereafter referred to as a ‘voxel’, following the terminology of Preisler et al. (2004) and Vilar et al. (2010). We assigned burnt voxels the minimum burn date of the fire in which they were encompassed. The grid in our study area consisted of 4,254,289 cells, and there were 6,878 daily observations for each cell during a 19 year study period, resulting in 2.93e10 voxels. To reduce the dataset to a manageable (i.e., computationally feasible) size, while retaining all voxels corresponding to a fire event, we sampled from a small proportion of the non-fire observations (i.e., control points). We sampled control points so that the centre of the voxel was evenly distributed across land cover types to avoid over-sampling the dominant land cover type, using the ‘raster’ 3.5-15 package in R (Hijmans et al., 2020). We excluded voxels as potential control points if a fire occurred in that cell within a year prior to the voxel date, as the chances of a fire re-occurring in that same voxel will be low, until the fuel load has recovered. 

We gathered data, at the level of both the fire patch and voxel, on the following variables: 

1. land cover type 
2. Normalized Difference Vegetation Index (NDVI) 
3. temperature
4. accumulated rainfall 
5. soil moisture content 
6. travel speed 
7. population density 
8. distance to nearest road and road density

and four fire danger metrics of the European Forest Fire Information System, specifically the: 

9. Drought Code (DC)
10. Keetch-Byram drought index (KBDI) 
11. Fire Weather Index (FWI)
12. Fine Fuel Moisture Code (FFMC)

We extracted histograms for each land cover type occurring within voxels using the Earth Engine Code Editor and the aforementioned Polesia land cover map. These histograms can then be converted into proportions and aggregated for fire patches using the fire patch ID. We used publicly available remote sensing imagery and other ready-to-use products from the Google Earth Engine platform to gather data on variables (2) to (7) using ‘rgee’ 1.1.3, a binding package for calling THE Google Earth Engine API from within R (Ayba et al., 2022). We combined Moderate Resolution Imaging Spectroradiometer (MODIS) NDVI datasets available at 15-day intervals (8-day intervals when combined). We extracted NDVI data using the closest image in time that predated the fire or voxel date, ensuring the composite time interval of that image did not overlap with a fire. We collated maximum daily temperature and accumulated daily rainfall from ERA5, the latest climate reanalysis produced by the European Centre for Medium-Range Weather Forecasts (ECMWF) (Hersbach et al., 2020), which combines model data with observations from across the world into a globally complete and consistent dataset of different climatic parameters. Information on soil moisture content is not available at a daily time-scale from as far back as 2001, so we used mean monthly soil moisture data, obtained from the TerraClimate dataset (Abatzoglou et al., 2018). For fires or voxels dated to the beginning of the month (before the 15th), soil moisture content from the previous month was extracted, otherwise data from the same month were used. 

We used a global friction surface to extract travel speeds (Weiss et al., 2018, 2020) and the Gridded Population of the World Version 4.11 to extract population density (CIESIN, 2017). The latter is available at five-year intervals from 2000. We extracted population densities using data from the year closest in time to the date of each fire or voxel. We obtained the vector datasets from the Global Roads Inventory Project, available at https://www.globio.info/, to calculate the shortest straight-line distance, in metres, to the nearest road from the voxel and the road density (Meijer et al., 2018). This can then be summarised for fire patches using the fire patch ID or repeated for ignition points. Shapefiles for ignition points can be created from the FRY fire patch dataset through the code in 'processing-fire-data.R'. We did this using the ‘sf’ package 1.0-7 in R (Pebesma et al., 2020).

We extracted daily fire danger metrics, available from the Copernicus Climate Data Store at  https://cds.climate.copernicus.eu/, which are based on the same ECMWF ERA5 reanalysis as the temperature and rainfall data, using the R package ‘raster’ 3.5-15 (Hijmans et al., 2020). The DC is an indicator of the moisture content in deep compact organic layers. The KBDI indicates the amount of precipitation necessary to return the deep duff and upper soil layers to saturated conditions. Both provide a good measure of prolonged drying conditions over the previous few months. The FWI indicates potential fire intensity while the FFMC is an indicator of the moisture content in litter and other cured fine fuels (e.g., needles, mosses, small twigs). 

We overlaid our fire data and control points with all other maps, and for covariates (2) to (12) extracted mean values within voxels. For each fire patch, we extracted minimum or maximum values for covariates (2) to (5) over the lifetime of the fire, assuming that fires are more likely a product of extreme environmental conditions rather than average tendencies (Amiro et al., 2004). In order to allow processing over the very large gridded dataset, we mapped functions over each individual voxel. Occasionally, we split the dataset into halves and processed each half separately (e.g., when extracting data on temperature and rainfall) so as not to exceed Google Earth Engine’s user memory limit. 

## User guide

Within ‘processing-fire-data.R’, you will need to change the specified working directory to where you are storing your fire trait dataset and shapefiles, as well as your own land cover map. The amount of files and folders you need to rename and set up should be fairly minimal, if using fire data from the FRY database. When sampling non-fire observations across land cover types, you may need to change the values used to remove cells centred around water or unclassified land cover types, depending on the values used in your land cover map. You may also need to change the number of cells you want to sample depending on how many fire observations you have and the ratio of fire to non-fire observations you want to obtain. Once you have processed the fire data, all other scripts build upon the outputs from this file. 

With all of the above fire data processed, all you need to do is make sure that the name of the shapefiles is the same as the name of the variable used in the script e.g., ‘all_pixels.shp’ or ‘fire_patches.shp’. 

You can then run ‘GEE-land-cover’ in the Earth Engine Code Editor to extract land cover types from your land cover map, which will need to be uploaded to your assets and imported to your script with the name ‘lc’. This script also includes code to reclassify and view land cover maps, which is currently commented out, so that running the script will simply extract data from the raw land cover map you provide. You can also run ‘GEE-grid-cell-variable-extraction.R’ and 'GEE-fire-patch-variable-extraction.R' to extract data from the Google Earth Engine catalogue for voxels and fire patches, respectively. ‘GRIP.R’ and ‘FWI-grid-cell.R’ require a vector and raster dataset on roads to be downloaded from https://www.globio.info/download-grip-dataset and a global fire danger indices to be downloaded from https://cds.climate.copernicus.eu/cdsapp#!/dataset/cems-fire-historical?tab=overview. Therefore you will need to change the working directories in these scripts to the relevant folders where the data are being stored, as well as where you would like to export the output to. For the former, you can select specific regions of interest to speed up performance. For the latter, you will need to select the years for which you would like to download data - you may need to download in batches. We chose four fire danger indices, the Drought Code, Keetch-Byram Drought Index, Fine Fuel Moisture Code, and Fire Weather Index. Advanced users may wish to modify the above scripts to extract data for additional or alternative variables. 

Using the datasets provided, you can run a simple generalised addative model following the ‘example-analysis.R’ script. You may need to change your values for k. 

## Acknowledgements

Thanks to Wentao Chen and Florent Mouillot of Le Centre National de la Recherche Scientifique. Their input was invaluable in the acquiring and preparation of the fire data. Thanks also to Mark de Jong of the Canadian Forest Service and Thomas Dowling from UNEP for creating the land cover map for Polesia. Thanks to Mark for his invaluable knowledge of fire ecology and his continuous input and feedback. Thanks to Megan Critchley and Susana Baena from UNEP for their help with managing and handling the remote sensing data. Finally, we are grateful to the Endangered Landscapes Programme for funding this work.

## References

Abatzoglou, J. T., Dobrowski, S. Z., Parks, S. A. & Hegewisch, K. C. (2018) TerraClimate, a high-resolution global dataset of monthly climate and climatic water balance from 1958–2015. Sci. Data 5, 170191.

Amiro, B. D. et al. (2004) Fire weather index system components for large fires in the Canadian boreal forest. Int. J. Wildland Fire 13, 391.

Ayba, C. et al. (2022) rgee: R Bindings for Calling the ‘Earth Engine’ API. 

Center for International Earth Science Information Network (2017) Gridded Population of the World, Version 4 (GPWv4): Population Density, Revision 10. 

Gorelick, N. et al. (2017) Google Earth Engine: Planetary-scale geospatial analysis for everyone. Remote Sens. Environ. 202, 18–27.

Hersbach, H. et al. (2020) The ERA5 global reanalysis. Q. J. R. Meteorol. Soc. 146, 1999–2049.

Hijmans, R. J. et al. (2020) raster: Geographic Data Analysis and Modeling. 

Laurent, P. et al. (2018) FRY, a global database of fire patch functional traits derived from space-borne burned area products. Sci. Data 5, 180132.

Meijer, J. R., Huijbregts, M. A. J., Schotten, K. C. G. J. & Schipper, A. M. (2018) Global patterns of current and future road infrastructure. Environ. Res. Lett. 13, 064006.

Pebesma, E. et al. (2020) sf: Simple Features for R.

Preisler, H., Brillinger, D. R., Burgan, R. E. & Benoit, J. (2004) Probability based models for estimation of wildfire risk. Int. J. Wildland Fire 132 133-142 13, 133–142.

Tamiminia, H. et al. (2020) Google Earth Engine for geo-big data applications: A meta-analysis and systematic review. ISPRS J. Photogramm. Remote Sens. 164, 152–170.

Vilar, L., Woolford, D. G., Martell, D. L. & Martín, M. P. (2010) A model for predicting human-caused wildfire occurrence in the region of Madrid, Spain. Int. J. Wildland Fire 19, 325–337.

Weiss, D. J. et al. (2018) A global map of travel time to cities to assess inequalities in accessibility in 2015. Nature (2018) doi:10.1038/nature25181.

Weiss, D. J. et al. (2020) Global maps of travel time to healthcare facilities. Nat. Med. 26, 1835–1838.
