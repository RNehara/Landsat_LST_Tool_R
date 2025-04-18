# Land Surface Temperature (LST) Estimation from Landsat 8/9 – R Implementation  
This repository contains an R-based tool for estimating Land Surface Temperature (LST) from Landsat 8/9 Collection 2 Level-1 imagery, using the Practical Single-Channel Algorithm proposed by Wang et al. (2019).  
It supports both daytime and nighttime image processing and is tailored for environmental and urban climate applications worldwide.    

**Features**  
🔹 Calculates LST from Landsat 8/9 Band 10 (Thermal Infrared Sensor - TIRS)   
🔹 Radiometric and reflectance correction using metadata (TXT & XML)  
🔹 Emissivity estimation using NDVI and reflectance bands based on Li & Jiang (2018)  
🔹 Atmospheric water vapor estimation using meteorological data  
🔹 Full support for daytime and nighttime acquisitions  
🔹 Outputs LST in degrees Celsius as raster files  
🔹 Recommended: use satellite and meteorological data from dates as close as possible to reduce atmospheric variation    

**Input Requirements**  
Landsat 8/9 Collection 2 Level-1 imagery (day time and nighttime imagery)  
Bands: 2–7 (OLI), 10–11 (TIRS)  
Metadata: ANG and MTL.txt and .xml file  
Meteorological data from weather station must include: air temperature (°C), relative humidity (%), and atmospheric pressure (mbar)  
Should be from the same day, location, and the hour immediately after the satellite acquisition time (e.g., if the image was acquired at 13:15 PM, use meteorological data from 14:00 PM).    

**How it Works**  
**Preprocessing**  
Reads raster bands and metadata  
Vegetation Indices & Emissivity  
Computes NDVI  
Calculates Land Surface Emissivity (LSE) using coefficients from Li & Jiang (2018):  
**Meteorological Adjustment**  
Computes atmospheric water vapor using the empirical formula by Buck (1981):  
**LST Estimation**  
Applies the Practical Single-Channel Algorithm from Wang et al. (2019):  

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
Wang, M., Zhang, Z., Hu, T., & Liu, X. (2019). A practical single‐channel algorithm for land surface temperature retrieval: application to landsat series data. Journal of Geophysical Research: Atmospheres, 124(1), 299-316. 10.1029/2018JD029330  
Li, S., & Jiang, G. M. (2018). Land surface temperature retrieval from Landsat-8 data with the generalized split-window algorithm. IEEE Access, 6, 18149-18162. 10.1109/ACCESS.2018.2818741  
Buck, A. L. (1981). New equations for computing vapor pressure and enhancement factor. Journal of Applied Meteorology and Climatology, 20(12), 1527-1532. 10.1175/1520-0450(1981)020<1527:NEFCVP>2.0.CO;2  
Rech, B., Moreira, R. N., Mello, T. A. G., Klouček, T., & Komárek, J. (2024). Assessment of daytime and nighttime surface urban heat islands across local climate zones–A case study in Florianópolis, Brazil. Urban Climate, 55, 101954. https://doi.org/10.1016/j.uclim.2024.101954    

**Author**  
Rodrigo Nehara  
r.nehara@usp.br  
GitHub: @RNehara  
