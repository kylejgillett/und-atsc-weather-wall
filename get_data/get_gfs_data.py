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
from concurrent.futures import ThreadPoolExecutor
from netCDF4 import Dataset
from siphon.catalog import TDSCatalog
from siphon.ncss import NCSS
from datetime import datetime, timedelta
from xarray.backends import NetCDF4DataStore
import xarray as xr


def _fetch_gfs_hour(ncss, fh, time_to_select, center_lat, center_lon, box_size):
    """
    Download one GFS forecast hour from NCSS as raw NetCDF bytes.

    Runs in a background thread so the next forecast hour downloads
    while the current one is being plotted. Only raw bytes are returned
    here -- the NetCDF/HDF5 file is opened on the main thread, since
    HDF5 is not thread-safe.
    """

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

    return ncss.get_data_raw(query)


def gfs_forecast(center_lat=37.86, center_lon=-98.61, box_size=20, forecast_hours=[12, 24, 36, 48, 60], min_pressure=5000):
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
    
    # Retrieve only requested forecast hours.
    # One forecast hour is downloaded ahead in a background thread while
    # the current one is plotted, so at most ~2 hours are held in memory.
    with ThreadPoolExecutor(max_workers=1) as pool:
        
        next_fetch = pool.submit(_fetch_gfs_hour, ncss, forecast_hours[0], base_time + timedelta(hours=forecast_hours[0]),
                                 center_lat, center_lon, box_size)
        
        for i, fh in enumerate(forecast_hours):
            
            time_to_select = base_time + timedelta(hours=fh)
            fetch = next_fetch
            
            # Start downloading the next forecast hour
            if i + 1 < len(forecast_hours):
                next_fh = forecast_hours[i + 1]
                next_fetch = pool.submit(_fetch_gfs_hour, ncss, next_fh, base_time + timedelta(hours=next_fh),
                                         center_lat, center_lon, box_size)
            
            # Wait for this forecast hour
            try:
                raw_data = fetch.result()
            except Exception as e:
                print(f"    GFS FH {fh} FAILED -- {e}")
                continue
            
            # Open the downloaded bytes as NetCDF (in memory)
            ncss_data = Dataset(f"gfs_f{fh:03d}.nc", memory=raw_data)
            
            # Convert to Xarray dataset
            forecast_data = xr.open_dataset(
                NetCDF4DataStore(ncss_data)
            ).metpy.parse_cf().load()
            
            # Remove singleton time dimensions
            for dim in list(forecast_data.dims):
                if dim.startswith('time') and forecast_data.sizes[dim] == 1:
                    forecast_data = forecast_data.squeeze(dim=dim, drop=True)

            # Drop stratospheric levels (1-40 hPa) that are never plotted,
            # keeping 50 hPa so the 2PVU calculation has a level above 100 hPa
            forecast_data = forecast_data.sel(isobaric=slice(min_pressure, None))

            print(f'    GFS FH {fh} COMPLETE: '
                  f'{time_to_select.strftime("%Y-%m-%d %H:%M:%SZ")}')
            
            # Send forecast to plotting script
            yield fh, time_to_select, forecast_data
            forecast_data.close()
            ncss_data.close()
            del forecast_data
            del ncss_data
            del raw_data
            gc.collect()
    
    elapsed_time = comp_time.time() - st
    
    print('ALL FORECASTS COMPLETE. Time elapsed:',
        comp_time.strftime("%H:%M:%S", comp_time.gmtime(elapsed_time)))