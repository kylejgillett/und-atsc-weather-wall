##########################################################
#      GFK 7-DAY NWS & GEFS FORECAST DATA GETTER
#  (c) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2026
##########################################################

# imports
import time as comp_time
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


# KGFK location and data endpoints
GFK_LAT = 47.9493
GFK_LON = -97.1761
USER_AGENT = "UND-ATSC-Weather-Wall (weather data display)"
NWS_API = "https://api.weather.gov"
GEFS_API = "https://ensemble-api.open-meteo.com/v1/ensemble"
GEFS_META = "https://ensemble-api.open-meteo.com/data/ncep_gefs025/static/meta.json"

# NWS gridded weather types grouped into the precip types shown on the meteogram
FZRA_TYPES = {"freezing_rain", "freezing_drizzle"}
SLEET_TYPES = {"sleet", "ice_pellets"}
SNOW_TYPES = {"snow", "snow_showers"}
RAIN_TYPES = {"rain", "rain_showers", "drizzle"}
THUNDER_TYPES = {"thunderstorms", "hail"}

# non-precip NWS weather types, shortened to METAR-style codes for the WX row
WX_CODES = {
    "fog": "FG",
    "freezing_fog": "FZFG",
    "ice_fog": "ICE FG",
    "haze": "HZ",
    "smoke": "FU",
    "blowing_snow": "BLSN",
    "blowing_dust": "BLDU",
    "blowing_sand": "BLSA",
    "dust": "DU",
    "volcanic_ash": "VA",
    "ice_crystals": "IC",
    "frost": "FROST",
    "freezing_spray": "FZSPRAY",
}



# shared request session with retries for flaky connections and rate limits
def build_session(accept=None):

    retry = Retry(total=4, connect=4, read=4, backoff_factor=1.0, status_forcelist=(429, 500, 502, 503, 504),
                  allowed_methods=("GET",), respect_retry_after_header=True)
    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT})
    if accept is not None:
        session.headers.update({"Accept": accept})
    session.mount("https://", HTTPAdapter(max_retries=retry))

    return session



# split NWS gridpoint layers ("2026-11-10T05:00:00+00:00/PT6H") into (start, end, value)
def parse_grid_layer(layer):

    intervals = []
    for item in layer.get("values", []):
        start, duration = item["validTime"].split("/")
        start = pd.Timestamp(start).tz_convert("UTC")
        intervals.append((start, start + pd.Timedelta(duration), item["value"]))

    return intervals



# fill an instantaneous gridpoint layer (temp, dewpoint, PoP...) onto an hourly index
def grid_layer_to_hourly(layer, hours):

    series = pd.Series(np.nan, index=hours)
    for start, end, value in parse_grid_layer(layer):
        if value is not None:
            series[(hours >= start) & (hours < end)] = value

    return series.to_numpy()



# place accumulation layers (QPF, snow, ice) at the end hour of each interval so a
# cumulative sum steps up when each accumulation period is complete
def grid_accum_to_hourly(layer, hours):

    series = pd.Series(0.0, index=hours)
    for start, end, value in parse_grid_layer(layer):
        end_hour = end - pd.Timedelta(hours=1)
        if value is None or end <= hours[0] or end_hour not in series.index:
            continue
        # keep only the future part of an interval already in progress
        fraction = (end - max(start, hours[0])) / (end - start)
        series[end_hour] += value * fraction

    return series.to_numpy()



# reduce an hour's NWS weather list to one precip type label
def classify_weather(conditions):

    types = {c.get("weather") for c in conditions if c.get("weather")}
    if types & FZRA_TYPES:
        return "Freezing Rain"
    if types & SLEET_TYPES or (types & SNOW_TYPES and types & RAIN_TYPES):
        return "Wintry Mix"
    if types & SNOW_TYPES:
        return "Snow"
    if types & THUNDER_TYPES:
        return "Thunderstorms"
    if types & RAIN_TYPES:
        return "Rain"

    return ""



# list the non-precip weather in an hour as codes ("FG BLSN"), in WX_CODES order
def obstruction_codes(conditions):

    types = {c.get("weather") for c in conditions if c.get("weather")}
    return " ".join(code for wx, code in WX_CODES.items() if wx in types)



def c_to_f(values):
    return np.asarray(values, dtype=float) * 9 / 5 + 32



###########################
# NWS GRIDPOINT FORECAST
###########################
def get_nws_forecast(lat=GFK_LAT, lon=GFK_LON, hours=168):

    print("    ACCESSING NWS GRIDPOINT FORECAST...")
    session = build_session(accept="application/geo+json")

    # resolve the forecast office grid for the point, then pull the raw gridpoint
    # data (one request holds temp, dewpoint, wind chill, PoP, QPF, snow, and weather)
    point = session.get(f"{NWS_API}/points/{lat:.4f},{lon:.4f}", timeout=20)
    point.raise_for_status()
    grid_url = point.json()["properties"]["forecastGridData"]

    response = session.get(grid_url, timeout=30)
    response.raise_for_status()
    grid = response.json()["properties"]
    issued = pd.Timestamp(grid["updateTime"]).tz_convert("UTC")

    # hourly index from the current hour to the end of the temperature forecast
    temp_end = parse_grid_layer(grid["temperature"])[-1][1]
    start = pd.Timestamp(datetime.now(timezone.utc)).floor("h")
    end = min(temp_end, start + pd.Timedelta(hours=hours))
    valid = pd.date_range(start, end, freq="h", inclusive="left")

    # instantaneous fields (gridpoint temps are in °C)
    df = pd.DataFrame({"valid_time": valid})
    df["temperature_f"] = c_to_f(grid_layer_to_hourly(grid["temperature"], valid))
    df["dewpoint_f"] = c_to_f(grid_layer_to_hourly(grid["dewpoint"], valid))
    df["wind_chill_f"] = c_to_f(grid_layer_to_hourly(grid["windChill"], valid))
    df["pop_percent"] = np.nan_to_num(grid_layer_to_hourly(grid["probabilityOfPrecipitation"], valid))
    df["sky_cover_percent"] = grid_layer_to_hourly(grid["skyCover"], valid)

    # wind (km/h -> kt)
    df["wind_speed_kt"] = grid_layer_to_hourly(grid["windSpeed"], valid) / 1.852
    df["wind_gust_kt"] = grid_layer_to_hourly(grid["windGust"], valid) / 1.852
    df["wind_dir_deg"] = grid_layer_to_hourly(grid["windDirection"], valid)

    # accumulations (mm -> in)
    df["qpf_in"] = grid_accum_to_hourly(grid["quantitativePrecipitation"], valid) / 25.4
    df["snow_in"] = grid_accum_to_hourly(grid["snowfallAmount"], valid) / 25.4
    df["ice_in"] = grid_accum_to_hourly(grid["iceAccumulation"], valid) / 25.4

    # precip type and non-precip weather (fog, haze, blowing snow...) from the weather grid
    precip_type = pd.Series("", index=valid)
    obstruction = pd.Series("", index=valid)
    for w_start, w_end, conditions in parse_grid_layer(grid["weather"]):
        in_window = (valid >= w_start) & (valid < w_end)
        precip_type[in_window] = classify_weather(conditions or [])
        obstruction[in_window] = obstruction_codes(conditions or [])
    df["precip_type"] = precip_type.to_numpy()
    df["obstruction"] = obstruction.to_numpy()

    # official daytime highs and overnight lows (period start/end + value)
    extremes = []
    for kind, layer in (("H", grid["maxTemperature"]), ("L", grid["minTemperature"])):
        for p_start, p_end, value in parse_grid_layer(layer):
            if value is not None and p_end > start and p_start < end:
                extremes.append({"kind": kind, "start": p_start, "end": p_end, "temperature_f": float(c_to_f(value))})
    extremes = pd.DataFrame(extremes, columns=["kind", "start", "end", "temperature_f"])

    df = df.dropna(subset=["temperature_f"]).reset_index(drop=True)
    if df.empty:
        raise RuntimeError("NWS GRIDPOINT FORECAST RETURNED NO TEMPERATURES")

    print(f"    NWS FORECAST COMPLETE.....{len(df)} HOURS, ISSUED {issued:%Y-%m-%d %H:%MZ}")

    return df, extremes, issued



###########################
# GEFS ENSEMBLE MEMBERS
###########################
def get_gefs_plumes(lat=GFK_LAT, lon=GFK_LON, days=8):

    st = comp_time.time()
    print("    ACCESSING GEFS 0.25° ENSEMBLE MEMBERS...")
    session = build_session()

    # Open-Meteo serves all 31 GEFS 0.25° members (control + 30 perturbations)
    # for a point in a single request, instead of hundreds of GRIB subsets
    params = {"latitude": lat, "longitude": lon, "models": "gfs025", "hourly": "temperature_2m,precipitation",
              "temperature_unit": "fahrenheit", "precipitation_unit": "inch", "timezone": "UTC", "forecast_days": days}
    response = session.get(GEFS_API, params=params, timeout=(10, 45))
    response.raise_for_status()
    data = response.json()["hourly"]

    # build [time x member] dataframes; the un-suffixed variable is the control run
    valid = pd.to_datetime(data["time"], utc=True)
    plumes = {}
    for var, name in (("temperature_2m", "temperature_f"), ("precipitation", "precipitation_in")):
        members = {"c00": data[var]}
        members.update({f"p{key[-2:]}": values for key, values in data.items() if key.startswith(f"{var}_member")})
        plumes[name] = pd.DataFrame(members, index=valid, dtype=float).rename_axis("valid_time")

    # model run time (not in the forecast response; read from the dataset metadata)
    try:
        meta = session.get(GEFS_META, timeout=15)
        meta.raise_for_status()
        init = pd.Timestamp(meta.json()["last_run_initialisation_time"], unit="s", tz="UTC")
    except (requests.RequestException, KeyError, ValueError) as e:
        print(f"    GEFS INIT TIME UNAVAILABLE.....{e}")
        init = None

    n_members = plumes["temperature_f"].shape[1]
    if n_members < 2:
        raise RuntimeError("GEFS ENSEMBLE MEMBERS UNAVAILABLE")

    elapsed_time = comp_time.time() - st
    init_str = f"{init:%Y-%m-%d %HZ}" if init is not None else "UNKNOWN"
    print(f"    GEFS COMPLETE.....{n_members} MEMBERS, INIT {init_str} "
          f"({comp_time.strftime('%H:%M:%S', comp_time.gmtime(elapsed_time))})")

    return plumes, init
