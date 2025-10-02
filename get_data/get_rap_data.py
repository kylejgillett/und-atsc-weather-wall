##########################################
# A SCRIPT TO LOAD NCSS RAP ANALYSIS AND
# RAP REANALYSIS DATASETS FOR BASIC
# MAP ANALYSIS
#
#
# To use: import pull_rap_data
# then call pull_rap_data.analysis() or
# pull_rap_data.reanalysis()
#
#
# (C) Kyle J Gillett 2024
##########################################

import warnings
warnings.filterwarnings("ignore")

# IMPORTS
import time as comp_time
import sys
import warnings

from siphon.catalog import TDSCatalog
from siphon.ncss import NCSS
from datetime import datetime

from xarray.backends import NetCDF4DataStore
import xarray as xr



### RAP ANALYSIS ###
#########################################################################################################
def analysis(center_lat=37.86, center_lon=-98.61, box_size=20):
    st = comp_time.time()

    print(f'    ACCESSING RAP DATA')

    # define dataset URL & try to access it to make sure it works
    url = 'http://thredds.ucar.edu/thredds/catalog/grib/NCEP/RAP/CONUS_13km/latest.xml'
    #url = 'https://thredds.ucar.edu/thredds/ncss/grid/grib/NCEP/RAP/CONUS_13km/RR_CONUS_13km_20250315_1800.grib2'
    try:
        cat = TDSCatalog(url)
        source = 'RAP ANALYSIS'
    except:
        sys.exit(
            "NCSS URL FAILED -- THIS HAPPENS WHEN A BAD REQUEST IS MADE. RAP Analysis data may not be available at this time.")
        pass

    # set up TDS query
    latest_ds = list(cat.datasets.values())[0]
    ncss = NCSS(latest_ds.access_urls['NetcdfSubset'])
    query = ncss.query()
    # Find start time
    start_time = ncss.metadata.time_span['begin']
    fcst_date = datetime.strptime(start_time, '%Y-%m-%dT%H:%M:%SZ')
    year1 = fcst_date.strftime('%Y')
    month1 = fcst_date.strftime('%m')
    day1 = fcst_date.strftime('%d')
    hour1 = fcst_date.strftime('%H')
    # Subset data by time
    query.time(fcst_date).accept('netcdf4')
    # Subsets data by variables
    query.variables('MSLP_MAPS_System_Reduction_msl',
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
                    'V-Component_Storm_Motion_height_above_ground_layer').add_lonlat()

    # Subset data by lat-lon domain
    query.lonlat_box(north=58.415, west=-139.970, east=-57.269, south=16.209)

    # Gets data
    ncss_data = ncss.get_data(query)

    raw_data = xr.open_dataset(NetCDF4DataStore(ncss_data)).metpy.parse_cf()


    elapsed_time = comp_time.time() - st
    print('    RAP DATA LOADED.....Time elapsed:', comp_time.strftime("%H:%M:%S", comp_time.gmtime(elapsed_time)))

    return raw_data


#########################################################################################################


### RAP REANALYSIS ###
#########################################################################################################
def reanalysis(year, month, day, hour, center_lat, center_lon, box_size):
    st = comp_time.time()

    print(f'    ACCESSING RAP DATA')

    # create dict of RAP-data urls for the different versions of RAP and RUC from NCEI
    urls = {
        #'RAP_cust': 'https://thredds.ucar.edu/thredds/ncss/grid/grib/NCEP/RAP/CONUS_13km/RR_CONUS_13km_20250315_1800.grib2',
        'RAP_13km': 'https://www.ncei.noaa.gov/thredds/ncss/model-rap130/' + str(year) + str(month) + '/' + str(
            year) + str(month) + str(day) + '/rap_130_' + str(year) + str(month) + str(day) + '_' + str(
            hour) + '00_000.grb2',
        'RAP_13km_old': 'https://www.ncdc.noaa.gov/thredds/ncss/model-rap130-old/' + str(year) + str(month) + '/' + str(
            year) + str(month) + str(day) + '/rap_130_' + str(year) + str(month) + str(day) + '_' + str(
            hour) + '00_000.grb2',

        'RAP_13km_anl': 'https://www.ncei.noaa.gov/thredds/ncss/model-rap130anl/' + str(year) + str(month) + '/' + str(
            year) + str(month) + str(day) + '/rap_130_' + str(year) + str(month) + str(day) + '_' + str(
            hour) + '00_000.grb2',
        'RAP_13km_anl_old': 'https://www.ncdc.noaa.gov/thredds/ncss/model-rap130anl-old/' + str(year) + str(
            month) + '/' + str(year) + str(month) + str(day) + '/rap_130_' + str(year) + str(month) + str(
            day) + '_' + str(hour) + '00_000.grb2',

        'RAP_25km': 'https://www.ncei.noaa.gov/thredds/ncss/model-rap252/' + str(year) + str(month) + '/' + str(
            year) + str(month) + str(day) + '/rap_252_' + str(year) + str(month) + str(day) + '_' + str(
            hour) + '00_000.grb2',
        'RAP_25km_old': 'https://www.ncei.noaa.gov/thredds/ncss/model-rap252-old/' + str(year) + str(month) + '/' + str(
            year) + str(month) + str(day) + '/rap_252_' + str(year) + str(month) + str(day) + '_' + str(
            hour) + '00_000.grb2',

        'RAP_25km_anl': 'https://www.ncei.noaa.gov/thredds/ncss/model-rap252anl/' + str(year) + str(month) + '/' + str(
            year) + str(month) + str(day) + '/rap_252_' + str(year) + str(month) + str(day) + '_' + str(
            hour) + '00_000.grb2',
        'RAP_25km_anl_old': 'https://www.ncei.noaa.gov/thredds/ncss/model-rap252anl-old/' + str(year) + str(
            month) + '/' + str(year) + str(month) + str(day) + '/rap_252_' + str(year) + str(month) + str(
            day) + '_' + str(hour) + '00_000.grb2',

        'RUC_13km': 'https://www.ncei.noaa.gov/thredds/ncss/model-ruc130anl/' + str(year) + str(month) + '/' + str(
            year) + str(month) + str(day) + '/ruc2anl_130_' + str(year) + str(month) + str(day) + '_' + str(
            hour) + '00_000.grb2',
        'RUC_13km_old': 'https://www.ncei.noaa.gov/thredds/ncss/model-ruc130anl-old/' + str(year) + str(
            month) + '/' + str(year) + str(month) + str(day) + '/ruc2anl_130_' + str(year) + str(month) + str(
            day) + '_' + str(hour) + '00_000.grb2',

        'RUC_25km': 'https://www.ncei.noaa.gov/thredds/ncss/model-ruc252anl/' + str(year) + str(month) + '/' + str(
            year) + str(month) + str(day) + '/ruc2anl_252_' + str(year) + str(month) + str(day) + '_' + str(
            hour) + '00_000.grb',
        'RUC_25km_old': 'https://www.ncei.noaa.gov/thredds/ncss/model-ruc252anl/' + str(year) + str(month) + '/' + str(
            year) + str(month) + str(day) + '/ruc2anl_252_' + str(year) + str(month) + str(day) + '_' + str(
            hour) + '00_000.grb'
    }

    # create a simple test for each URL, use the first one that works
    try:
        for url, key in zip(urls.values(), urls.keys()):
            try:
                NCSS(url)
                print(f'> DATASET USED: {key}')
                url_to_use = url
                source = str(key)[0:3]
                break
            except:
                pass
    except:
        pass
    try:
        ncss = NCSS(url_to_use)
    except:
        warnings.filterwarnings("ignore")
        sys.exit(
            'NCSS URL FAILED -- THIS HAPPENS WHEN A BAD REQUEST IS MADE.\n> CHECK TO MAKE SURE YOU ENTERED THE CORRECT DATES.\n> NOTE: DATA IS NOT AVAILIABLE FOR EVERY DATE\n> THIS CATALOG OFTEN EXPERIENCES MISSING DATA/OUTAGES\n> MAKE SURE DATES ARE STRINGS -- MONTH, DAY AND HOUR MUST BE TWO DIGITS (EX: 18, 06, 00)')
        pass

    # set up TDS query
    #query = ncss.query().accept('netcdf4').add_lonlat()
    #query = ncss.query().time(datetime(int(year), int(month), int(day), int(hour))).accept('netcdf4').add_lonlat()
    query = ncss.query().time(datetime(int(year), int(month), int(day), int(hour))).accept('netcdf').add_lonlat()
    # subset data by variable names for RAP & RUC (of course they have to be different)
    # query.variables('MSLP_MAPS_System_Reduction_msl',
    #                 'Pressure_surface',
    #                 'Geopotential_height_isobaric',
    #                 'Temperature_isobaric',
    #                 'Relative_humidity_isobaric',
    #                 'Temperature_height_above_ground',
    #                 'Relative_humidity_height_above_ground',
    #                 'u-component_of_wind_height_above_ground', 'v-component_of_wind_height_above_ground',
    #                 'u-component_of_wind_isobaric', 'v-component_of_wind_isobaric',
    #                 'Convective_available_potential_energy_surface',
    #                 'Storm_relative_helicity_height_above_ground_layer',
    #                 'U-Component_Storm_Motion_height_above_ground_layer',
    #                 'V-Component_Storm_Motion_height_above_ground_layer')
    if source in ['rap', 'RAP']:
            query.variables('MSLP_MAPS_System_Reduction_msl', 'Pressure_surface',
                            'Geopotential_height_isobaric', 'Geopotential_height_surface',
                            'Temperature_isobaric', 'Temperature_height_above_ground',
                            'Relative_humidity_isobaric', 'Relative_humidity_height_above_ground',
                            'Dewpoint_temperature_height_above_ground',
                            'Pseudo-adiabatic_potential_temperature_or_equivalent_potential_temperature_surface',
                            'u-component_of_wind_height_above_ground', 'v-component_of_wind_height_above_ground',
                            'u-component_of_wind_isobaric', 'v-component_of_wind_isobaric',
                            'Convective_available_potential_energy_surface',
                            # 'Convective_Available_Potential_Energy_surface',
                            'Vertical_velocity_pressure_isobaric',
                            'Categorical_Freezing_Rain_surface', 'Categorical_Ice_Pellets_surface',
                            'Categorical_Rain_surface', 'Categorical_Snow_surface',
                            'Composite_reflectivity_entire_atmosphere', 'Lightning_surface')

    else:
        query.variables('MSLP_MAPS_System_Reduction_msl',
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
                        'V-Component_Storm_Motion_height_above_ground_layer').add_lonlat()

    # subset data by requested domain
    # query.lonlat_box(center_lat+box_size, center_lat-box_size,
    #                  center_lon+box_size, center_lon-box_size)

    query.lonlat_box(north=58.415, west=-139.970, east=-57.269, south=16.209)

    # laod the data from TDS
    #raw_data = data.get_data(query)
    ncss_data = ncss.get_data(query)
    raw_data = xr.open_dataset(NetCDF4DataStore(ncss_data)).metpy.parse_cf()


    elapsed_time = comp_time.time() - st
    print('    RAP DATA LOADED.....Time elapsed:', comp_time.strftime("%H:%M:%S", comp_time.gmtime(elapsed_time)))
    return raw_data
#########################################################################################################