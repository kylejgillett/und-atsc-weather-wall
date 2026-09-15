##########################################
# A SCRIPT TO LOAD NCSS GFS FORECAST DATA
# FOR BASIC MAP ANALYSIS
# KYLE J GILLETT, UNV. NORTH DAKOTA, 2026
##########################################

import warnings
warnings.filterwarnings("ignore")

# IMPORTS
import gc
import time as comp_time
import sys
from siphon.catalog import TDSCatalog
from siphon.ncss import NCSS
from datetime import datetime, timedelta
from xarray.backends import NetCDF4DataStore
import xarray as xr


def gfs_forecast(center_lat=37.86, center_lon=-98.61, box_size=20, forecast_hours=[12, 24, 36, 48, 60]):
    st = comp_time.time()
    
    print('ACCESSING GFS DATA...')
    
    # Define dataset URL
    url = 'https://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/latest.xml'
    
    try:
        cat = TDSCatalog(url)
    except Exception as e:
        sys.exit(f"NCSS URL FAILED -- GFS data may not be available at this time. Error: {e}")
    
    # Get latest GFS dataset
    latest_ds = list(cat.datasets.values())[0]
    
    try:
        ncss = NCSS(latest_ds.access_urls['NetcdfSubset'])
    except Exception as e:
        sys.exit(f"NCSS DATASET FAILED -- GFS data may not be available at this time. Error: {e}")
    
    # Find analysis time
    start_time = ncss.metadata.time_span['begin']
    base_time = datetime.strptime(start_time, '%Y-%m-%dT%H:%M:%SZ')
    
    print(f'    GFS CYCLE: {base_time.strftime("%Y-%m-%d %H:%M:%SZ")}')
    
    # Retrieve only requested forecast hours
    for fh in forecast_hours:
        
        time_to_select = base_time + timedelta(hours=fh)
        
        print(f"    GETTING GFS FH: {fh}")
        
        # Create NCSS query
        query = ncss.query()
        
        # Request exact forecast time
        query.time(time_to_select)
        
        # Request NetCDF4
        query.accept('netcdf4')
        
        # Subset variables
        query.variables(
            'MSLP_Eta_model_reduction_msl',
            'Geopotential_height_isobaric',
            'Temperature_isobaric',
            'Temperature_height_above_ground',
            'u-component_of_wind_height_above_ground',
            'v-component_of_wind_height_above_ground',
            'u-component_of_wind_isobaric',
            'v-component_of_wind_isobaric',
            'Categorical_Rain_surface',
            'Categorical_Freezing_Rain_surface',
            'Categorical_Ice_Pellets_surface',
            'Categorical_Snow_surface',
            'Composite_reflectivity_entire_atmosphere'
        ).add_lonlat()
        
        # Subset by lat-lon domain
        query.lonlat_box(
            north=center_lat+box_size,
            west=center_lon-box_size,
            east=center_lon+box_size,
            south=center_lat-box_size
        )
        
        # Fetch forecast hour
        try:
            ncss_data = ncss.get_data(query)
        except Exception as e:
            print(f"    GFS FH {fh} FAILED -- {e}")
            continue
        
        # Convert to Xarray dataset
        forecast_data = xr.open_dataset(
            NetCDF4DataStore(ncss_data)
        ).metpy.parse_cf().load()
        
        # Remove singleton time dimensions
        for dim in list(forecast_data.dims):
            if dim.startswith('time') and forecast_data.sizes[dim] == 1:
                forecast_data = forecast_data.squeeze(dim=dim, drop=True)
        
        print(f'    GFS FH {fh} COMPLETE: '
              f'{time_to_select.strftime("%Y-%m-%d %H:%M:%SZ")}')
        
        # Send forecast to plotting script
        yield fh, time_to_select, forecast_data
        forecast_data.close()
        ncss_data.close()
        del forecast_data
        del ncss_data
        gc.collect()
    
    elapsed_time = comp_time.time() - st
    
    print('ALL FORECASTS COMPLETE. Time elapsed:',
        comp_time.strftime("%H:%M:%S", comp_time.gmtime(elapsed_time)))