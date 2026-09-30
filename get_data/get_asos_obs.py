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
import time
from io import StringIO
import requests


# fetch data
station_id = 'GFK'
def get_asos_obs(station_id, hours=48, max_retries=4):
    
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


    for attempt in range(1, max_retries + 1):

        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            # read returned CSV into pandas
            df = pd.read_csv(StringIO(response.text), skiprows=5)
            break

        except (requests.exceptions.RequestException, pd.errors.ParserError, pd.errors.EmptyDataError) as e:
            print(f"    {station_id} DATA REQUEST FAILED (attempt {attempt}/{max_retries}).....{e}")

            # If this was the final attempt, send the error
            # back to the calling script.
            if attempt == max_retries:
                raise RuntimeError(f"{station_id} IEM ASOS DATA UNAVAILABLE AFTER {max_retries} ATTEMPTS") from e

            # progressive retry delay:
            # 5 sec, 10 sec, 15 sec...
            wait_time = 5 * attempt

            print("    RETRYING {station_id} IN {wait_time} SECONDS...")

            time.sleep(wait_time)



    # clean up df, remove missings, set strings to numerics, drop nan rows, sort by time
    df = df.replace("M", np.nan)
    numeric_cols = ['tmpf', 'dwpf', 'relh', 'sknt', 'gust', 'vsby', 'feel', 'drct', 'mslp', 'p01i']
    df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors='coerce')
    df = df.dropna(subset=['tmpf'])
    df["valid"] = pd.to_datetime(df["valid"], errors="coerce", utc=True)
    df = df.dropna(subset=["valid"])
    df = df.sort_values('valid')

    # compute u and v components for vector plot, set a normalized length 
    df['u'], df['v'] = mpcalc.wind_components(df['sknt'].values*units.kts, df['drct'].values*units.degrees)
    wind_mag = np.sqrt(df["u"]**2 + df["v"]**2)
    df["u_norm"] = np.where(wind_mag > 0, df["u"] / wind_mag, np.nan)
    df["v_norm"] = np.where(wind_mag > 0, df["v"] / wind_mag, np.nan)
    df["y_arrow"] = np.full(len(df), 40.0)

    # return dataframe
    return df