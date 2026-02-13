##########################################################
#          AWOS/ASOS TIMESERIES DATA FROM IEM
#  (c) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2025
##########################################################

# imports 
from datetime import datetime, timedelta, timezone
import pandas as pd
import numpy as np
import metpy.calc as mpcalc
from metpy.units import units



# fetch data
station_id = 'GFK'
def get_asos_obs(station_id, hours=48):
    
    # define time params
    end = datetime.now(timezone.utc)
    start = end - timedelta(hours=hours)

    # set up IEM data fetch
    url = (
        "https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py?"
        f"station={station_id}&data=tmpf,dwpf,relh,sknt,gust,vsby,feel,drct,wxcodes,mslp,p01i&"
        f"tz=Etc/UTC&format=csv&latlon=no&"
        f"year1={start.year}&month1={start.month}&day1={start.day}&hour1={start.hour}&"
        f"year2={end.year}&month2={end.month}&day2={end.day}&hour2={end.hour}"
    )

    # fetch data via pandas 
    df = pd.read_csv(url, skiprows=5)

    # clean up df, remove missings, set strings to numerics, drop nan rows, sort by time
    df = df.replace("M", np.nan)
    numeric_cols = ['tmpf', 'dwpf', 'relh', 'sknt', 'gust', 'vsby', 'feel', 'drct', 'mslp', 'p01i']
    df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors='coerce')
    df = df.dropna(subset=['tmpf'])
    df['valid'] = pd.to_datetime(df['valid'])
    df = df.sort_values('valid')

    # compute u and v components for vector plot, set a normalized length 
    df['u'], df['v'] = mpcalc.wind_components(df['sknt'].values*units.kts, df['drct'].values*units.degrees)
    df['u_norm'] = df['u'] / np.sqrt(df['u']**2 + df['v']**2)
    df['v_norm'] = df['v'] / np.sqrt(df['u']**2 + df['v']**2)
    df['y_arrow'] = np.full_like(df['sknt'], 40)

    # return dataframe
    return df