##########################################################
#            LATEST METARS RETREVIAL SCRIPT
#  (c) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2025
##########################################################

import warnings
import socket
warnings.filterwarnings("ignore")

import io
import requests
from datetime import datetime, timedelta, timezone
from metpy.calc import reduce_point_density
from awips.dataaccess import DataAccessLayer
from dynamicserialize.dstypes.com.raytheon.uf.common.time import TimeRange
import numpy as np
import metpy.calc as mpcalc
from metpy.units import units
import pandas as pd
import cartopy.crs as ccrs


# seconds to wait on a stalled EDEX connection before giving up
EDEX_TIMEOUT = 60

# if aviationweather.gov returns fewer stations than this, also try EDEX and use whichever has more.
# normally ~5000 stations come back; during the Nov 2026 EDEX outage EDEX returned ~25
MIN_STATIONS = 500

# aviationweather.gov bulk METAR cache: the latest report from every station, updated every minute
AWC_METAR_CACHE = "https://aviationweather.gov/data/cache/metars.cache.csv.gz"

# drop aviationweather.gov reports older than this (hours)
AWC_MAX_AGE_HOURS = 2


def get_cloud_cover(code):
    if 'OVC' in code:
        return 8
    elif 'BKN' in code:
        return 6
    elif 'SCT' in code:
        return 4
    elif 'FEW' in code:
        return 2
    else:
        return 0


### METAR OBSERVATIONS ###
#########################################################################################################
def get_metar_data(date=datetime.utcnow(), reduced_to=100000):

    # aviationweather.gov first (faster, fresher, more stations); fall back to EDEX if it
    # fails or returns only a handful of stations
    try:
        data_df, time_title = _get_awc_metar_data(date)
    except Exception as e:
        print(f'    AVIATIONWEATHER.GOV METARS FAILED -- {e}')
        data_df, time_title = None, ""

    if data_df is None or len(data_df) < MIN_STATIONS:
        awc_count = 0 if data_df is None else len(data_df)
        print(f'    AVIATIONWEATHER.GOV RETURNED {awc_count} STATIONS -- TRYING EDEX METARS')
        try:
            edex_df, edex_time = _get_edex_metar_data(date)
        except Exception as e:
            print(f'    EDEX METARS FAILED -- {e}')
            edex_df, edex_time = None, ""

        # keep whichever source has more stations
        if edex_df is not None and (data_df is None or len(edex_df) > len(data_df)):
            data_df, time_title = edex_df, edex_time

    if data_df is None:
        raise RuntimeError("NO METAR DATA -- aviationweather.gov and EDEX both failed")

    if len(data_df) > 0:
        print(f'    METAR DATA LOADED.....{time_title}')

    proj = ccrs.LambertConformal()
    point_locs = proj.transform_points(ccrs.PlateCarree(), data_df['longitude'],
                                       data_df['latitude'])
    reduced_data = data_df[reduce_point_density(point_locs, reduced_to)]

    return reduced_data, time_title



def _get_edex_metar_data(date):

    # EDEX Request
    edexServer = "edex-cloud.unidata.ucar.edu"
    DataAccessLayer.changeEDEXHost(edexServer)
    request = DataAccessLayer.newDataRequest("obs")
    # define desired parameters
    single_value_params = ["stationName", "longitude", "latitude",
                           "temperature", "dewpoint", "windDir",
                           "windSpeed", "seaLevelPress"]
    multi_value_params = ["skyCover"]
    params = single_value_params + multi_value_params

    # set all parameters on the request
    request.setParameters(*(params))

    # Time range
    lastHourDateTime = date - timedelta(hours=1)
    start = lastHourDateTime.strftime('%Y-%m-%d %H')
    beginRange = datetime.strptime(start + ":00:00", "%Y-%m-%d %H:%M:%S")
    endRange = datetime.strptime(start + ":59:59", "%Y-%m-%d %H:%M:%S")
    timerange = TimeRange(beginRange, endRange)

    # Get response. python-awips opens its EDEX connection with no timeout, so a stalled
    # server would hang the job forever; use a socket timeout for this request only
    previous_timeout = socket.getdefaulttimeout()
    socket.setdefaulttimeout(EDEX_TIMEOUT)
    try:
        response = DataAccessLayer.getGeometryData(request, timerange)
    finally:
        socket.setdefaulttimeout(previous_timeout)

    # define a dictionary and array that will be populated from our for loop below
    obs = dict({params: [] for params in params})
    station_names = []
    time_title = ""
    i = 0

    # cycle through all the data in the response, in reverse order to get the most recent data first
    for ob in reversed(response):
        avail_params = ob.getParameters()
        # print(avail_params)
        # if it has cloud information, we want the last of the 6 entries (most recent)
        if "skyCover" in avail_params:
            if i == 5:
                # store the associated cloud cover int for the skyCover string
                obs['skyCover'].append(get_cloud_cover(ob.getString("skyCover")))
            i = i + 1
        elif "stationName" in avail_params:
            # If we already have a record for this stationName, skip
            if ob.getString('stationName') not in station_names:
                station_names.append(ob.getString('stationName'))
                i = 0
                if time_title == "":
                    time_title = str(ob.getDataTime())
                for param in single_value_params:
                    if param in avail_params:
                        try:
                            obs[param].append(ob.getNumber(param))
                        except TypeError:
                            obs[param].append(ob.getString(param))
                    else:
                        obs[param].append(None)

    data = dict()
    data['stid'] = np.array(obs['stationName'])
    data['latitude'] = np.array(obs['latitude'])
    data['longitude'] = np.array(obs['longitude'])
    data['air_temperature'] = np.array(obs['temperature'], dtype=float) * units.degC
    data['dew_point_temperature'] = np.array(obs['dewpoint'], dtype=float) * units.degC
    data['sea_level_pressure'] = (np.array(obs['seaLevelPress'], dtype=float) / 100) * units.mb

    direction = np.array(obs['windDir'])
    direction[direction == -9999.0] = 'nan'

    METAR_u, METAR_v = mpcalc.wind_components(np.array(obs['windSpeed']) * units('knots'), direction * units.degree)
    data['eastward_wind'], data['northward_wind'] = METAR_u, METAR_v
    data['cloud_coverage'] = np.array(obs['skyCover'])

    data_df = pd.DataFrame.from_dict(data)

    return data_df, time_title



def _get_awc_metar_data(date):
    """
    Latest METAR from every station via the aviationweather.gov bulk cache,
    returned with the same columns and units as the EDEX data.
    """

    response = requests.get(AWC_METAR_CACHE, timeout=60,
                            headers={"User-Agent": "UND-Atmospheric-Sciences-Weather-Wall/1.0"})
    response.raise_for_status()
    metars = pd.read_csv(io.BytesIO(response.content), compression="gzip")

    # keep reports from the last AWC_MAX_AGE_HOURS before `date` (naive UTC, as for EDEX)
    obs_time = pd.to_datetime(metars['observation_time'], utc=True)
    valid_end = pd.Timestamp(date).tz_localize(timezone.utc) + pd.Timedelta(minutes=15)
    keep = (obs_time <= valid_end) & (obs_time >= valid_end - pd.Timedelta(hours=AWC_MAX_AGE_HOURS))
    keep &= metars['latitude'].notna() & metars['longitude'].notna()
    metars = metars[keep].reset_index(drop=True)
    obs_time = obs_time[keep].reset_index(drop=True)

    # total sky cover = the most-covered reported layer (OVC > BKN > SCT > FEW > clear)
    layer_cols = [col for col in metars.columns if col.startswith('sky_cover')]
    cloud_cover = metars[layer_cols].fillna("").astype(str).apply(lambda col: col.map(get_cloud_cover)).max(axis=1)

    data = dict()
    data['stid'] = metars['station_id'].to_numpy()
    data['latitude'] = metars['latitude'].to_numpy(dtype=float)
    data['longitude'] = metars['longitude'].to_numpy(dtype=float)
    data['air_temperature'] = metars['temp_c'].to_numpy(dtype=float) * units.degC
    data['dew_point_temperature'] = metars['dewpoint_c'].to_numpy(dtype=float) * units.degC
    data['sea_level_pressure'] = metars['sea_level_pressure_mb'].to_numpy(dtype=float) * units.mb

    # variable or missing directions come through as NaN (no barb drawn), same as EDEX -9999
    direction = metars['wind_dir_degrees'].to_numpy(dtype=float)
    speed = metars['wind_speed_kt'].to_numpy(dtype=float)
    METAR_u, METAR_v = mpcalc.wind_components(speed * units('knots'), direction * units.degree)
    data['eastward_wind'], data['northward_wind'] = METAR_u, METAR_v
    data['cloud_coverage'] = cloud_cover.to_numpy(dtype=int)

    data_df = pd.DataFrame.from_dict(data)

    # label with the hour most reports fall in, in the same format as the EDEX time
    time_title = ""
    if len(obs_time) > 0:
        time_title = obs_time.median().floor('h').strftime('%Y-%m-%d %H:%M:%S')

    return data_df, time_title