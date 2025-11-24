##########################################
# A SCRIPT TO LOAD NCSS GFS ANALYSIS AND
# FOR BASIC MAP ANALYSIS
##########################################

import warnings
warnings.filterwarnings("ignore")

# IMPORTS
import time as comp_time
import sys
from siphon.catalog import TDSCatalog
from siphon.ncss import NCSS
from datetime import datetime, timedelta
from xarray.backends import NetCDF4DataStore
import xarray as xr
import numpy as np # Added for consistency in calling np.arange if needed outside this function

def gfs_forecast(center_lat=37.86, center_lon=-98.61, box_size=20, forecast_hours=[6, 12, 18, 24, 36, 48, 60, 72]):
    st = comp_time.time()
    
    print('ACCESSING GFS DATA...')
    
    # define dataset URL
    url = 'https://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/latest.xml'
    
    try:
        cat = TDSCatalog(url)
    except Exception as e:
        sys.exit(f"NCSS URL FAILED -- GFS data may not be available at this time. Error: {e}")
    
    latest_ds = list(cat.datasets.values())[0]
    ncss = NCSS(latest_ds.access_urls['NetcdfSubset'])
    
    # Find analysis time
    start_time = ncss.metadata.time_span['begin']
    base_time = datetime.strptime(start_time, '%Y-%m-%dT%H:%M:%SZ')
    
    # Define the datetime objects for all desired forecast hours
    query_times = [base_time + timedelta(hours=fh) for fh in forecast_hours]

    query = ncss.query()

    # --- FIX APPLIED HERE: Removed the invalid .subset_nearest() method ---
    # Request a time range from the earliest to the latest forecast hour
    query.time_range(query_times[0], query_times[-1])
    
    query.accept('netcdf4')
    
    # Subset variables
    query.variables(
        'MSLP_Eta_model_reduction_msl',
        'Pressure_surface',
        'Geopotential_height_isobaric',
        'Temperature_isobaric',
        'Relative_humidity_isobaric',
        'Temperature_height_above_ground',
        'Relative_humidity_height_above_ground',
        'u-component_of_wind_height_above_ground', 'v-component_of_wind_height_above_ground',
        'u-component_of_wind_isobaric', 'v-component_of_wind_isobaric',
        "Categorical_Rain_surface","Categorical_Freezing_Rain_surface",
        "Categorical_Ice_Pellets_surface", "Categorical_Snow_surface", 
        "Composite_reflectivity_entire_atmosphere"
    ).add_lonlat()
        
    # Subset by lat-lon domain (using your current wide domain)
    query.lonlat_box(north=58.415, west=-139.970, east=-57.269, south=16.209)
    
    # Fetch data (single request for all times)
    ncss_data = ncss.get_data(query)

    # Open the full 4D dataset (Time, Level, Lat, Lon)
    raw_data_full = xr.open_dataset(NetCDF4DataStore(ncss_data)).metpy.parse_cf()
        
    # Iterate through the full dataset to extract individual forecast hours using Xarray's nearest selection
    all_forecasts = {}
    
    print('Data received. Processing individual forecast hours...')

    for fh in forecast_hours:
        time_to_select = base_time + timedelta(hours=fh) 
        print(f"    GETIING GFS FH: {fh}")
        
        # Use Xarray's .sel with method='nearest' to handle any mismatch
        try:
            forecast_slice = raw_data_full.sel(time=time_to_select, method='nearest')
        except:
            try:
                forecast_slice = raw_data_full.sel(time1=time_to_select, method='nearest')
            except:
                forecast_slice = raw_data_full.sel(time2=time_to_select, method='nearest')
                pass
            pass
        
        all_forecasts[fh] = forecast_slice
        print(f'Processed GFS forecast for +{fh} hours: {time_to_select.strftime("%Y-%m-%d %H:%M:%SZ")}')

    elapsed_time = comp_time.time() - st
    print('ALL FORECASTS LOADED. Time elapsed:', comp_time.strftime("%H:%M:%S", comp_time.gmtime(elapsed_time)))
    
    return all_forecasts