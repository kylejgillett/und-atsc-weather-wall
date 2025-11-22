##########################################
# A SCRIPT TO LOAD NCSS GFS ANALYSIS AND
# FOR BASIC MAP ANALYSIS
#
#
#
#
# (C) Kyle J Gillett 2025
##########################################

import warnings
warnings.filterwarnings("ignore")

# IMPORTS
import time as comp_time
import sys
import warnings
from siphon.catalog import TDSCatalog
from siphon.ncss import NCSS
from datetime import datetime, timedelta
from xarray.backends import NetCDF4DataStore
import xarray as xr




def gfs_forecast(center_lat=37.86, center_lon=-98.61, box_size=20, forecast_hours=[6, 12, 18, 24, 36, 48, 60, 72]):
    st = comp_time.time()
    
    print('ACCESSING GFS DATA...')
    
    # define dataset URL
    url = 'https://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/latest.xml'
    
    try:
        cat = TDSCatalog(url)
    except:
        sys.exit("NCSS URL FAILED -- GFS data may not be available at this time.")
    
    latest_ds = list(cat.datasets.values())[0]
    ncss = NCSS(latest_ds.access_urls['NetcdfSubset'])
    
    # Find analysis time
    start_time = ncss.metadata.time_span['begin']
    base_time = datetime.strptime(start_time, '%Y-%m-%dT%H:%M:%SZ')
    
    all_forecasts = {}
    
    for fh in forecast_hours:
        query_time = base_time + timedelta(hours=fh)
        query = ncss.query()
        query.time(query_time).accept('netcdf4')
        
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
            'Convective_available_potential_energy_surface',
            'Storm_relative_helicity_height_above_ground_layer',
            'U-Component_Storm_Motion_height_above_ground_layer',
            'V-Component_Storm_Motion_height_above_ground_layer'
        ).add_lonlat()
        
        # Subset by lat-lon domain
        query.lonlat_box(north=58.415, west=-139.970, east=-57.269, south=16.209)
        
        # Fetch data
        ncss_data = ncss.get_data(query)

        raw_data = xr.open_dataset(NetCDF4DataStore(ncss_data)).metpy.parse_cf()
        
        all_forecasts[fh] = raw_data
        print(f'Loaded GFS forecast for +{fh} hours: {query_time}')
    
    elapsed_time = comp_time.time() - st
    print('ALL FORECASTS LOADED. Time elapsed:', comp_time.strftime("%H:%M:%S", comp_time.gmtime(elapsed_time)))
    
    return all_forecasts




