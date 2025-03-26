import copernicusmarine as cm
import os
import pandas as pd

# This will download biogeochemsitry data and sea surface temperature data

#Set login information
# sign up for an account and put in this information to download the data
USERNAME = ""
PASSWORD = ""

# Authenticate (this creates a configuration file with your credentials)
cm.login(USERNAME, PASSWORD)

# Set dataset to download
dataset_ids = [
    "cmems_mod_med_bgc-bio_anfc_4.2km_P1M-m",
    "cmems_mod_med_bgc-car_anfc_4.2km_P1M-m",
    "cmems_mod_med_bgc-co2_anfc_4.2km_P1M-m",
    "cmems_mod_med_bgc-nut_anfc_4.2km_P1M-m",
    "cmems_mod_med_bgc-optics_anfc_4.2km_P1M-m",
    "cmems_mod_med_bgc-pft_anfc_4.2km_P1M-m",
    "omi_health_chl_medsea_oceancolour_area_averaged_mean",
    "cmems_SST_MED_SST_L4_REP_OBSERVATIONS_010_021"
]

#Download each data set with filters
for bioG in dataset_ids:    
   ds = cm.open_dataset(
       dataset_id=bioG,
       username=USERNAME,
       password=PASSWORD,
       start_datetime="2000-01-01",
       end_datetime="2025-02-20",
       minimum_longitude=23,
       maximum_longitude=37,
       minimum_latitude=30,
       maximum_latitude=37
   )
   
   ds.to_netcdf(f"{bioG}.nc")

print("Datasets have been downloaded")
    
