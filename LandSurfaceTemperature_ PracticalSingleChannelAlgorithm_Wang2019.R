###########################################################################
# ======================================================================
#                                        Rodrigo Nehara
#                                        r.nehara@usp.br
#                                        https://github.com/RNehara
#                                        04-16-2025 (mm/dd/yyyy)
# Land surface temperature estimation for daytime and nighttime
# Landsat 8/9 Collection 2 Level-1 image
# Running R 4.5.x
# ======================================================================

# ------------------------------------------------------------

###########################################################################
# Required packages
library(lubridate)   # for working with dates
library(readr)       # for reading csv files
library(terra)       # raster manipulation
library(XML)         # XML manipulation
library(magrittr)    # pipe operator
library(stringr)     # string manipulation
library(tibble)      # data frames
library(dplyr)       # data manipulation

# ------------------------------------------------------------
# edit here
# ------------------------------------------------------------

    # This script calculates land surface temperature according to the Practical Single-Channel Algorithm proposed by Wang et al. (2019).
    # This script requires Landsat 8/9 Collection 2 Level-1 bands 2, 3, 4, 5, 6, 7, 10 and 11 (daytime and nighttime).
    # Keep only the necessary bands and all .txt and .XML files in the directories.
    # Do not rename the bands and metadata file in the folder.
    # It is recommended that the satellites imagery come from dates as close as possible to minimize atmospheric inconsistencies.
    # Meteorological data from weather station must include: air temperature (°C), relative humidity (%), and atmospheric pressure (mbar)
    # Should be from the same day, location, and the hour immediately after the satellite acquisition time (e.g., if the image was acquired at 13:15 PM, use meteorological data from 14:00 PM).


    # dir: the directory (folder paths) where the images are stored.
    # scene_id: the id of image (ex.: LC08_L1TP_219076_20190614_20200828_02_T1).
    # output: the directory (folder paths) where the images should be exported.

# Define paths for daytime
dir_day <- "C:/Users/nehar/OneDrive/Documentos/Artigo_Clima_Urbano/Landsat_Image/daytime/"
scene_id_day <- "LC09_L1TP_219076_20240619_20240620_02_T1"

# Define paths for nighttime
dir_night <- "C:/Users/nehar/OneDrive/Documentos/Artigo_Clima_Urbano/Landsat_Image/Nightime/"
scene_id_night <- "LC09_L1GT_101168_20230607_20230607_02_T2"

# Meteorological variables for daytime manually defined by the user
t_air_day <- 24.7 # t_air: Air temperature (°C)
rel_hum_day <- 0.423 # rel_hum: Relative humidity (decimal, e.g., 0.62 for 62%)
p_atm_day <- 938.8 # p_atm: Atmospheric pressure (mbar or hPa)

# Meteorological variables for nighttime manually defined by the user
t_air_night <- 13.9 # t_air: Air temperature (°C)
rel_hum_night <- 0.946 # rel_hum: Relative humidity (decimal, e.g., 0.62 for 62%)
p_atm_night <- 937.3 # p_atm: Atmospheric pressure (mbar or hPa)

output <- "C:/Users/nehar/OneDrive/Documentos/Artigo_Clima_Urbano/Landsat_Image/tesste/"

# -------------- edition ends here ---------------------------

cat('Calculating land surface temperature for daytime\n')

cat('Processing daytime images\n')
# Process day images
img_list_day <- list.files(path = dir_day, pattern = '.*TIF$', full.names = TRUE)
image_day <- rast(img_list_day)
names(image_day) <- lapply(names(image_day), function(name) str_sub(name, 42))
metadata_name_day <- list.files(path = dir_day, pattern = '.*xml$', full.names = TRUE)
metadata_day <- XML::xmlToDataFrame(metadata_name_day)
consolidate_rows <- function(df) {
  consolidated <- lapply(seq_along(df), function(i) {
    non_na_values <- na.omit(df[, i])
    if (length(non_na_values) > 0) {
      return(non_na_values[1])
    } else {
      return(NA)
    }
  })
  result_df <- as.data.frame(t(consolidated), stringsAsFactors = FALSE)
  colnames(result_df) <- colnames(df)
  return(result_df)
}
consolidated_metadata_day <- consolidate_rows(metadata_day)
metadata_day <- consolidated_metadata_day

# OLI bands processing
oli_bands <- names(image_day) %>% setdiff(., c('B10', 'B11'))
for (oli in oli_bands) {
  add_ref <- metadata_day %>%  select(starts_with("REFLECTANCE_ADD_BAND_")) %>% select(2:7) %>% sapply(as.numeric) 
  mult_ref <- metadata_day %>% select(starts_with("REFLECTANCE_MULT_BAND_")) %>% select(2:7) %>% sapply(as.numeric) 
  sun_elev <- (metadata_day$SUN_ELEVATION %>% na.omit() %>% as.numeric()) * pi / 180
  image_day$new_ref <- (image_day[[oli]] * mult_ref + add_ref) / sin(sun_elev)
  names(image_day)[length(names(image_day))] <- paste0(oli, '_REF')
}

# Thermal bands processing
tirs_bands <- c('B10', 'B11')
for (tirs in tirs_bands) {
  add_rad <- metadata_day %>% select(starts_with("RADIANCE_ADD_BAND")) %>% select(10:11) %>% sapply(as.numeric)
  mult_rad <- metadata_day %>% select(starts_with("RADIANCE_MULT_BAND")) %>% select(10:11) %>% sapply(as.numeric)
  k1 <- metadata_day %>% select(starts_with("K1_CONSTANT_BAND_")) %>% sapply(as.numeric)
  k2 <- metadata_day %>% select(starts_with("K2_CONSTANT_BAND_")) %>% sapply(as.numeric)
  image_day$new_RAD <- image_day[[tirs]] * mult_rad + add_rad
  names(image_day)[length(names(image_day))] <- paste0(tirs, '_RAD')
  image_day$new_BT <- k2 / log(k1 / image_day[[paste0(tirs, '_RAD')]] + 1)
  names(image_day)[length(names(image_day))] <- paste0(tirs, '_BT')
}
cat('It is done\n')

cat('Calculating land surface emissivity\n')
# Emissivity calculation for day images
# Reference NDVI's
ndvi_s <- 0.2
ndvi_v <- 0.5

# Reference emissivity for B10 (Li & Jiang, 2018)
e_s <- 0.971
e_v <- 0.982

# Surface geometric factor
f <- 0.55

# Coefficients for the LSE model (Li & Jiang, 2018)
coeffs <- c(0.980, -0.140, 0.170, -0.036, -0.083, 0.158, -0.149)

for (i in 1:length(coeffs)) {
  assign(paste0('a', i), coeffs[i])
}
# Bands of reflectance (OLI) and radiance (TIRS)
band_names <- c('B2_REF', 'B3_REF', 'B4_REF', 'B5_REF', 'B6_REF', 'B7_REF', 'B10_RAD')
var_names <- c('B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B10')
for (b in 1:length(band_names)) {
  assign(var_names[b], image_day[[band_names[b]]])
}

# NDVI
ndvi <- (B5 - B4) / (B5 + B4)

# Reclassify NDVI for Pv calculation

# Matrix of reclassification
# NDVI below NDVI of bare soil becomes equal to NDVI of bare soil
# NDVI above NDVI of vegetation becomes equal to NDVI of vegetation

ndvi_reclass_matrix <- matrix(c(-10, ndvi_s, ndvi_s, ndvi_v, 10, ndvi_v), ncol = 3, byrow = TRUE)
ndvi_reclass <- classify(ndvi, rcl = ndvi_reclass_matrix, include.lowest = TRUE)

# Fractional vegetation cover
Pv <- ((ndvi_reclass - ndvi_s) / (ndvi_v - ndvi_s)) ^ 2

# Term accounting for cavity effect
de <- (1 - e_s) * e_v * f * (1 - Pv)
de <- mask(de, ndvi >= ndvi_v, maskvalues = 1, updatevalue = 0.005)

# CASE 1: NDVI < NDVI OF BARE SOIL
e_mask_1 <- ndvi < ndvi_s
emiss_1 <- (
  a1 +
    a2 * mask(B2, e_mask_1, maskvalues = 0) +
    a3 * mask(B3, e_mask_1, maskvalues = 0) +
    a4 * mask(B4, e_mask_1, maskvalues = 0) +
    a5 * mask(B5, e_mask_1, maskvalues = 0) +
    a6 * mask(B6, e_mask_1, maskvalues = 0) +
    a7 * mask(B7, e_mask_1, maskvalues = 0)
)

# CASE 2: NDVI OF BARE SOIL <= NDVI <= NDVI OF VEGETATION
e_mask_2 <- ndvi >= ndvi_s & ndvi <= ndvi_v
emiss_2 <- mask(Pv, e_mask_2, maskvalues = 0) * (e_v - e_s) + e_s + de

# CASE 3: NDVI > NDVI OF VEGETATION
e_mask_3 <- ndvi > ndvi_v
emiss_3 <- e_v + mask(de, e_mask_3, maskvalues = 0)

# FINAL EMISSIVITY
emiss <- app(c(emiss_1, emiss_2, emiss_3), fun = max, na.rm = TRUE)
image_day$emiss <- emiss
cat('It is done\n')

cat('Calculating land surface temperature\n')

# Atmospheric water-vapor content (g/cm²) (Buck, 1981)
awv <- (0.098 * rel_hum_day * 6.1121 * (1 + 7.2e-4 + (3.2e-6 + 5.9e-10 * t_air_day^2) * p_atm_day) * exp(((18.729 - t_air_day / 227.3) * t_air_day) / (t_air_day + 257.87)))
cat('Atmospheric water vapor (g/cm²):', awv)

## TEMPERATURE RETRIEVAL ##

# Coefficients
lst_coeffs <- matrix(c(
  -0.28009, 1.257429, 0.275109, -1.32876, -0.1696,
  0.999069, 0.033453, 0.015232,
  -0.60336, 1.613485, -4.98989, 2.772703, -1.04271,
  1.739598, -0.54978, 0.129006,
  2.280539, 0.918191, -38.3363, 13.82581, -1.75455,
  5.003919, -1.62832, 0.196687,
  -0.4107, 1.493577, 0.278271, -1.22502, -0.31067,
  1.022016, -0.01969, 0.036001
), byrow = FALSE, nrow = 8)

# Selection of coefficients according to AWV ranges
if (awv <= 2) {
  column <- 1
} else if (awv > 2 & awv <= 4) {
  column <- 2
} else if (awv > 4 & awv <= 7) {
  column <- 3
} else {
  column <- 4
}
for (i in 1:nrow(lst_coeffs)) {
  assign(paste0('a', i - 1), lst_coeffs[i, column])
}
emiss <- image_day$emiss
B10 <- image_day$B10_RAD

# Coefficients (Wang et al., 2019)
c1 <- 1.19104e+08
c2 <- 1.43877e+04
lambda <- 10.904 # effective wavelength for B10

# Blackbody radiance from land surface
BT <- a0 + a1 * awv + (a2 + a3 * awv + a4 * awv^2) / emiss + (a5 + a6 * awv + a7 * awv^2) * B10 / emiss

# Land surface temperature (K)
LST <- (c2 / lambda) / log((c1 / (lambda^5 * BT)) + 1)
LST_day <- LST - 273.15
plot(LST_day)
image_day$LST <- LST - 273.15
cat('It is done\n')

cat('Processing nighttime images\n')

# Process night images
img_list_night <- list.files(path = dir_night, pattern = '.*TIF$', full.names = TRUE)
image_night <- rast(img_list_night)
names(image_night) <- lapply(names(image_night), function(name) str_sub(name, 42))
metadata_name_night <- list.files(path = dir_night, pattern = '.*xml$', full.names = TRUE)
metadata_night <- XML::xmlToDataFrame(metadata_name_night)
consolidate_rows <- function(df) {
  consolidated <- lapply(seq_along(df), function(i) {
    non_na_values <- na.omit(df[, i])
    if (length(non_na_values) > 0) {
      return(non_na_values[1])
    } else {
      return(NA)
    }
  })
  result_df <- as.data.frame(t(consolidated), stringsAsFactors = FALSE)
  colnames(result_df) <- colnames(df)
  return(result_df)
}
consolidated_metadata_night <- consolidate_rows(metadata_night)
metadata_night <- consolidated_metadata_night

# Thermal bands processing
tirs_bands <- c('B10', 'B11')
for (tirs in tirs_bands) {
  add_rad <- metadata_night %>% select(starts_with('RADIANCE_ADD_BAND')) %>% select(10:11) %>% sapply(as.numeric)
  mult_rad <- metadata_night %>% select(starts_with('RADIANCE_MULT_BAND')) %>% select(10:11) %>% sapply(as.numeric)
  k1 <- metadata_night %>% select(starts_with('K1_CONSTANT_BAND_')) %>% sapply(as.numeric)
  k2 <- metadata_night %>% select(starts_with('K2_CONSTANT_BAND_')) %>% sapply(as.numeric)
  image_night$new_RAD <- image_night[[tirs]] * mult_rad + add_rad
  names(image_night)[length(names(image_night))] <- paste0(tirs, '_RAD')
  image_night$new_BT <- k2 / log(k1 / image_night[[paste0(tirs, '_RAD')]] + 1)
  names(image_night)[length(names(image_night))] <- paste0(tirs, '_BT')
}

# Bands of radiance (TIRS)
band_names <- c('B10_RAD')
var_names <- c('B10')

for (b in 1:length(band_names)) {
  assign(var_names[b], image_night[[band_names[b]]])
}

# Resampling LSE band
emiss <- resample(emiss, B10, method = "bilinear")
cat('It is done\n')

cat('Calculating land surface temperature\n')

# Atmospheric water-vapor content (g/cm²) (Buck, 1981)
awv <- (0.098 * rel_hum_night * 6.1121 * (1 + 7.2e-4 + (3.2e-6 + 5.9e-10 * t_air_night^2) * p_atm_night) * exp(((18.729 - t_air_night / 227.3) * t_air_night) / (t_air_night + 257.87)))
cat('/nAtmospheric water vapor (g/cm²):', awv)

# Coefficients
lst_coeffs <- matrix(c(
  -0.28009, 1.257429, 0.275109, -1.32876, -0.1696,
  0.999069, 0.033453, 0.015232,
  -0.60336, 1.613485, -4.98989, 2.772703, -1.04271,
  1.739598, -0.54978, 0.129006,
  2.280539, 0.918191, -38.3363, 13.82581, -1.75455,
  5.003919, -1.62832, 0.196687,
  -0.4107, 1.493577, 0.278271, -1.22502, -0.31067,
  1.022016, -0.01969, 0.036001
), byrow = FALSE, nrow = 8)

# Selection of coefficients according to AWV ranges
if (awv <= 2) {
  column <- 1
} else if (awv > 2 & awv <= 4) {
  column <- 2
} else if (awv > 4 & awv <= 7) {
  column <- 3
} else {
  column <- 4
}
for (i in 1:nrow(lst_coeffs)) {
  assign(paste0('a', i - 1), lst_coeffs[i, column])
}

B10 <- image_night$B10_RAD

# Coefficients (Wang et al., 2019)
c1 <- 1.19104e+08
c2 <- 1.43877e+04
lambda <- 10.904 # effective wavelength for B10

# Blackbody radiance from land surface
BT <- a0 + a1 * awv + (a2 + a3 * awv + a4 * awv^2) / emiss + (a5 + a6 * awv + a7 * awv^2) * B10 / emiss

# Land surface temperature (K)
LST <- (c2 / lambda) / log((c1 / (lambda^5 * BT)) + 1)
image_night$LST <- LST
LST_night <- LST - 273.15
plot(LST_night)
image_night$LST <- LST - 273.15
cat('It is done\n')

cat('Exporting LST rasters (ºC)\n')
writeRaster(LST_day, filename = file.path(output, "LST_day.tif"), overwrite = TRUE)
writeRaster(LST_night, filename = file.path(output, "LST_night.tif"), overwrite = TRUE)
cat('It is done\n')



# References
# Wang, M., Zhang, Z., Hu, T., & Liu, X. (2019). 
# A practical single‐channel algorithm for land surface temperature retrieval: application to landsat series data.
# Journal of Geophysical Research: Atmospheres, 124(1), 299-316. 10.1029/2018JD029330
                             
# Li, S., & Jiang, G. M. (2018).
# Land surface temperature retrieval from Landsat-8 data with the generalized split-window algorithm.
# IEEE Access, 6, 18149-18162. 10.1109/ACCESS.2018.2818741
                             
# Buck, A. L. (1981).
# New equations for computing vapor pressure and enhancement factor.
# Journal of Applied Meteorology and Climatology, 20(12), 1527-1532. 10.1175/1520-0450(1981)020<1527:NEFCVP>2.0.CO;2

# Rech, B., Moreira, R. N., Mello, T. A. G., Klouček, T., & Komárek, J. (2024).
# Assessment of daytime and nighttime surface urban heat islands across local climate zones–A case study in Florianópolis, Brazil.
# Urban Climate, 55, 101954. https://doi.org/10.1016/j.uclim.2024.101954
