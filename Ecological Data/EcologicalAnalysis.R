
# Load libraries (some I need to remove)
library(ncdf4)
library(tidyverse)
library(reshape2)
library(future.apply)
library(data.table)
library(DT)
library(terra)
library(plotly)
library(stars)
library(purrr)

# Set active directory to current file location
parent_dir <- setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Temporarily disabled (HUGE FILE FOR SOME REASON)

# process_sst <- function() {
#   # File path
#   file <- "C:/Users/lucas/Desktop/ASU Capstone/Data/BioChemistry/SST_MED_SST_L4_REP_OBSERVATIONS_subset.nc"
#   
#   # Get time and depth info
#   ncData <- nc_open(file)
#   names(ncData$var)
#   analysed_sst_dt <- stars::read_ncdf(file, var = "analysed_sst", proxy = FALSE)
#   analysed_sst_dt <- st_as_stars(analysed_sst_dt)
#   analysed_sst_dt <- as.data.frame(analysed_sst_dt, xy=TRUE)
#   analysed_sst_dt <- analysed_sst_dt %>% na.omit() #%>% mutate(o2 = str_sub(o2, 1, 8))
#   analysed_sst_dt$analysed_sst <- as.numeric(analysed_sst_dt$analysed_sst)
#   
#   
#   analysis_error_dt <- stars::read_ncdf(file, var = "analysis_error", proxy = FALSE)
#   analysis_error_dt <- st_as_stars(analysis_error_dt)
#   analysis_error_dt <- as.data.frame(analysis_error_dt, xy=TRUE)
#   analysis_error_dt <-analysis_error_dt %>% na.omit() %>% mutate(nppv = str_sub(analysis_error, 1, 8))
#   analysis_error_dt$nppv <- as.numeric(analysis_error_dt$analysis_error)
#   
#   nc_close(ncData)
#   gc()
#   # If needed: summarize by time
#   noCoord_O2 <- o2_dt %>% select(time, o2, depth)
#   o2_summary <- o2_dt %>% group_by(time, depth) %>% summarise(mean_DO = mean(o2))
# }


# Function to process Dissolved Oxygen and Net Primary Production
process_bio <- function() {
  
  # File path
  file <- paste0(parent_dir,"/BioChemistry/cmems_mod_med_bgc-bio_anfc_4.2km_P1M-m.nc")
  
  ncData <- nc_open(file)
  o2_dt <- stars::read_ncdf(file, var = "o2", proxy = FALSE)
  o2_dt <- o2_dt[, , , 1:6] # Keep only 13 meters depth of sea
  o2_dt <- st_as_stars(o2_dt)
  o2_dt <- as.data.frame(o2_dt, xy=TRUE)
  o2_dt <<- o2_dt %>% 
    na.omit() %>% 
    mutate(o2 = str_sub(o2, 1, 8)) %>% #strip off unit
    mutate(o2 = as.numeric(o2)) %>%
    rename("Dissolved Oxygen [mmol/m^3]" = o2)
  
  nppv_dt <- stars::read_ncdf(file, var = "nppv", proxy = FALSE)
  nppv_dt <- nppv_dt[, , , 1:6] # Keep only 13 meters depth of sea
  nppv_dt <- st_as_stars(nppv_dt)
  nppv_dt <- as.data.frame(nppv_dt, xy=TRUE) 
  nppv_dt <<- nppv_dt %>% 
    na.omit() %>% 
    mutate(nppv = str_sub(nppv, 1, 8)) %>% #strip off unit
    mutate(nppv = as.numeric(nppv)) %>%
    rename("Net Primary Production [mg/d/m^3]" = nppv)
  
  nc_close(ncData)
  gc()
}

#Function to process Dissolved Inorganic Carbon, Ocean pH and Alkalinity
process_ph_talk_dissic <- function() {
  # File path
  file <- paste0(parent_dir,"/BioChemistry/cmems_mod_med_bgc-car_anfc_4.2km_P1M-m.nc")
  
  # Get time and depth info
  ncData <- nc_open(file)
  names(ncData$var)
  dissic_dt <- stars::read_ncdf(file, var = "dissic", proxy = FALSE)
  dissic_dt <- dissic_dt[, , , 1:6] # Keep only 13 meters depth of sea
  dissic_dt <- st_as_stars(dissic_dt)
  dissic_dt <- as.data.frame(dissic_dt, xy=TRUE)
  dissic_dt <<- dissic_dt %>% 
    na.omit() %>%
    mutate(dissic = str_sub(dissic, 1, 8)) %>%#strip off unit
    mutate(dissic = as.numeric(dissic)) %>% 
    rename("Dissolved Inorganic Carbon [mol/m^3]" = dissic)
  
  
  
  pH_dt <- stars::read_ncdf(file, var = "ph", proxy = FALSE)
  pH_dt <- pH_dt[, , , 1:6]# Keep only 13 meters depth of sea
  pH_dt <- st_as_stars(pH_dt)
  pH_dt <- as.data.frame(pH_dt, xy=TRUE) 
  pH_dt <<- pH_dt %>% 
    na.omit() %>% 
    mutate(ph = str_sub(ph, 1, 8)) %>% #strip off unit
    mutate(ph = as.numeric(ph)) %>%
    rename("Ocean pH [1]" = ph)
  
  talk_dt <- stars::read_ncdf(file, var = "talk", proxy = FALSE)
  talk_dt <- talk_dt[, , , 1:6]# Keep only 13 meters depth of sea
  talk_dt <- st_as_stars(talk_dt)
  talk_dt <- as.data.frame(talk_dt, xy=TRUE) 
  talk_dt <<- talk_dt %>% 
    na.omit() %>% 
    mutate(talk = str_sub(talk, 1, 8)) %>% #strip off unit
    mutate(talk = as.numeric(talk)) %>%
    rename("Alkalinity [mol/m^3]" = talk)
  
  nc_close(ncData)
  gc() #clears unused memory
  
}

# Function for processing Surface partial pressure of CO2 and Surface CO2 flux
process_co2 <- function() {
  # File path
  file <- paste0(parent_dir,"/BioChemistry/cmems_mod_med_bgc-co2_anfc_4.2km_P1M-m.nc")
  
  # Get time and depth info
  ncData <- nc_open(file)
  names(ncData$var)
  fgco2_dt <- stars::read_ncdf(file, var = "fgco2", proxy = FALSE)
  fgco2_dt <- st_as_stars(fgco2_dt)
  fgco2_dt <- as.data.frame(fgco2_dt, xy=TRUE)
  fgco2_dt <<- fgco2_dt %>% 
    na.omit() %>% 
    mutate(fgco2 = as.numeric(fgco2)) %>%
    rename("Surface partial pressure of CO2 [Pa]" = fgco2)
  
  
  spco2_dt <- stars::read_ncdf(file, var = "spco2", proxy = FALSE)
  spco2_dt <- st_as_stars(spco2_dt)
  spco2_dt <- as.data.frame(spco2_dt, xy=TRUE) 
  spco2_dt <<- spco2_dt %>% 
    na.omit() %>% 
    mutate(spco2 = as.numeric(spco2)) %>% 
    rename("Surface CO2 flux [kg m-2 s-1]" = spco2)
  
  nc_close(ncData)
  gc()
  
}

# Function for processing Ammonium, Nitrate, Phosphate and Silicate
process_car <- function() {
  # File path
  file <- paste0(parent_dir,"/BioChemistry/cmems_mod_med_bgc-nut_anfc_4.2km_P1M-m.nc")
  
  # Get time and depth info
  ncData <- nc_open(file)
  names(ncData$var)
  nh4_dt <- stars::read_ncdf(file, var = "nh4", proxy = FALSE)
  nh4_dt <- nh4_dt[, , , 1:6]# Keep only 13 meters depth of sea
  nh4_dt <- st_as_stars(nh4_dt)
  nh4_dt <- as.data.frame(nh4_dt, xy=TRUE) 
  nh4_dt <<- nh4_dt %>% 
    na.omit() %>% 
    mutate(nh4 = str_sub(nh4, 1, 11)) %>% #strip off unit
    mutate(nh4 = as.numeric(nh4)) %>% 
    rename("Ammonium [mmol/m^3]" = nh4)
  
  
  no3_dt <- stars::read_ncdf(file, var = "no3", proxy = FALSE)
  no3_dt <- no3_dt[, , , 1:6]# Keep only 13 meters depth of sea
  no3_dt <- st_as_stars(no3_dt)
  no3_dt <- as.data.frame(no3_dt, xy=TRUE)
  no3_dt <<- no3_dt %>% 
    na.omit() %>% 
    mutate(no3 = str_sub(no3, 1, 12)) %>% #strip off unit
    mutate(no3 = as.numeric(no3)) %>% 
    rename("Nitrate [mmol/m^3]" = no3)
  
  po4_dt <- stars::read_ncdf(file, var = "po4", proxy = FALSE)
  po4_dt <- po4_dt[, , , 1:6] # Keep only 13 meters depth of sea
  po4_dt <- st_as_stars(po4_dt)
  po4_dt <- as.data.frame(po4_dt, xy=TRUE)
  po4_dt <<- po4_dt %>%
    na.omit() %>% 
    mutate(po4 = str_sub(po4, 1, 12)) %>% #strip off unit
    mutate(po4 = as.numeric(po4)) %>%
    rename("Phosphate [mmol/m^3]" = po4)
  
  si_dt <- stars::read_ncdf(file, var = "si", proxy = FALSE)
  si_dt <- si_dt[, , , 1:6]# Keep only 13 meters depth of sea
  si_dt <- st_as_stars(si_dt)
  si_dt <- as.data.frame(si_dt, xy=TRUE)
  si_dt <<- si_dt %>% 
    na.omit() %>% 
    mutate(si = str_sub(si, 1, 8)) %>% #strip off unit
    mutate(si = as.numeric(si)) %>%
    rename("Silicate [mmol/m^3]" = si)
  nc_close(ncData)
  gc()
  
}

#process_sst()

process_bio()

process_ph_talk_dissic()

process_co2()

process_car()


#CO2 measurements arent included yet as they dont have a depth related.

dt_list <- list(dissic_dt, nh4_dt, po4_dt, si_dt, talk_dt, pH_dt, o2_dt, nppv_dt, no3_dt)

BIG_FRAME <- reduce(dt_list, inner_join, by=c("latitude", "longitude", "depth", "time"))
