**Land Surface Temperature (LST) Estimation from Landsat 8/9 – R Implementation**  
This repository contains an R-based tool for estimating Land Surface Temperature (LST) from Landsat 8/9 Collection 2 Level-1 imagery, using the Practical Single-Channel Algorithm proposed by Wang et al. (2019).  
It supports both daytime and nighttime image processing and is tailored for environmental and urban climate applications worldwide.    

**Features**  
🔹 Calculates LST from Landsat 8/9 Band 10 (Thermal Infrared Sensor - TIRS)   
🔹 Radiometric and reflectance correction using metadata (MTL & XML)  
🔹 Emissivity estimation using NDVI and reflectance bands based on Li & Jiang (2018)  
🔹 Atmospheric water vapor estimation using meteorological data  
🔹 Full support for daytime and nighttime acquisitions  
🔹 Outputs LST in degrees Celsius as raster files  
🔹 Recommended: use satellite and meteorological data from dates as close as possible to reduce atmospheric variation    

**Input Requirements**  
Landsat 8/9 Collection 2 Level-1 imagery  
Bands: 2–7 (OLI), 10–11 (TIRS)  
Metadata: .MTL.txt or .xml file  
Meteorological data must include: air temperature (°C), relative humidity (%), and atmospheric pressure (mbar)  
Should be from the same day, location, and the hour immediately after the satellite acquisition time (e.g., if the image was acquired at 13:15 PM, use meteorological data from 14:00 PM).    

**How it Works**  
**Preprocessing**  
Reads raster bands and metadata  
Vegetation Indices & Emissivity  
Computes NDVI  
Calculates Land Surface Emissivity (LSE) using coefficients from Li & Jiang (2018):  
Li, Z., & Jiang, J. (2018). A new method for estimating broadband emissivity of land surfaces for land surface temperature retrieval. International Journal of Remote Sensing, 39(13), 4341–4360.    
**Meteorological Adjustment**  
Computes atmospheric water vapor using the empirical formula by Buck (1981):  
Buck, A.L. (1981). New equations for computing vapor pressure and enhancement factor. Journal of Applied Meteorology and Climatology, 20(12), 1527–1532.    
**LST Estimation**  
Applies the Practical Single-Channel Algorithm from Wang et al. (2019):  
Wang, F., Qin, Z., Song, C., Tu, L., Karnieli, A., & Zhao, S. (2019). An improved single-channel algorithm for land surface temperature retrieval from Landsat 8 thermal infrared sensor data. Remote Sensing, 11(5), 522.
Export  
Saves daytime and nighttime LST outputs as GeoTIFF files (LST_day.tif and LST_night.tif)    

**Required R Packages**  
library(terra)  
library(lubridate)  
library(readr)  
library(XML)  
library(magrittr)  
library(stringr)  
library(tibble)  
library(dplyr)    

**Edit the paths at the beginning of the script:**  
dir_day <- "path/to/daytime/scene/"  
dir_night <- "path/to/nighttime/scene/"  
scene_id_day <- "LC08_L1TP_..."  
scene_id_night <- "LC09_L1GT_..."  
output <- "path/to/output/folder/"  

**Edit the meteorological data at the beginning of the script:**  
**Meteorological variables for daytime manually defined by the user**  
t_air_day <-  t_air: Air temperature (°C)  
rel_hum_day <- rel_hum: Relative humidity (decimal, e.g., 0.62 for 62%)  
p_atm_day <- p_atm: Atmospheric pressure (mbar or hPa)    

**Meteorological variables for nighttime manually defined by the user:**  
t_air_night <- t_air: Air temperature (°C)  
rel_hum_night <- rel_hum: Relative humidity (decimal, e.g., 0.62 for 62%)  
p_atm_night <- p_atm: Atmospheric pressure (mbar or hPa)    

**Output**  
LST_day.tif: Raster of daytime land surface temperature (in °C)  
LST_night.tif: Raster of nighttime land surface temperature (in °C)    

**References**  
Wang, F. et al. (2019). An improved single-channel algorithm for land surface temperature retrieval from Landsat 8 thermal infrared sensor data. Remote Sensing, 11(5), 522. https://doi.org/10.3390/rs11050522  
Li, Z., & Jiang, J. (2018). A new method for estimating broadband emissivity of land surfaces. Int. J. Remote Sensing, 39(13), 4341–4360. https://doi.org/10.1080/01431161.2017.1420936  
Buck, A.L. (1981). New equations for computing vapor pressure and enhancement factor. J. Appl. Meteorol., 20(12), 1527–1532.    

**Author**  
Rodrigo Nehara  
r.nehara@usp.br  
GitHub: @RNehara  
