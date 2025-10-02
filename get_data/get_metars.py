##########################################################
#            LATEST METARS RETREVIAL SCRIPT
#  (c) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2025
##########################################################

import warnings
warnings.filterwarnings("ignore")

from datetime import datetime, timedelta
from metpy.calc import reduce_point_density
from awips.dataaccess import DataAccessLayer
from dynamicserialize.dstypes.com.raytheon.uf.common.time import TimeRange
import numpy as np
import metpy.calc as mpcalc
from metpy.units import units
import pandas as pd
import cartopy.crs as ccrs


### METAR OBSERVATIONS ###
#########################################################################################################
def get_metar_data(date=datetime.utcnow(), reduced_to=100000):
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

    # Get response
    response = DataAccessLayer.getGeometryData(request, timerange)

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

    if len(METAR_u) > 0:
        print('    METAR DATA LOADED.....{time_title}')

    data_df = pd.DataFrame.from_dict(data)
    proj = ccrs.LambertConformal()
    point_locs = proj.transform_points(ccrs.PlateCarree(), data_df['longitude'],
                                       data_df['latitude'])
    reduced_data = data_df[reduce_point_density(point_locs, reduced_to)]

    return reduced_data, time_title