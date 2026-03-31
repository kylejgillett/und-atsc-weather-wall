##########################################################
#             NWS HEADLINES LOADING SCRIPT
#  (c) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2026
##########################################################

import pandas as pd
import requests
import geopandas as gpd
from shapely.geometry import shape

url = "https://api.weather.gov/alerts/active"

PHENOM_DESCRIPTIONS = {
    'WI': 'Winter Weather',
    'HW': 'High Wind',
    'BZ': 'Blizzard',
    'FW': 'Flood',
    'WW': 'Winter Weather',
    'FZ': 'Freezing',
    'SC': 'Special Marine',
    'GL': 'Gale',
    'DU': 'Dust Storm',
    'HT': 'Heat',
    'UP': 'Unknown Precip',
    'WS': 'Winter Storm',
    'SR': 'Severe Thunderstorm',
    'FA': 'Flood',
    'FL': 'Flood',
    'CW': 'Coastal Flood',
    'LO': 'Low Water',
    'RP': 'Rip Current',
    'IS': 'Ice Storm',
    'LS': 'Lake Shore',
    'SU': 'Small Craft',
    'LW': 'Lake Wind',
    'TR': 'Tornado',
}

SIG_DESCRIPTIONS = {
    'W': 'Warning',
    'A': 'Advisory',
    'Y': 'Watch',
    'S': 'Statement',
}


def _wwa_event_name(row):
    phenom = str(row.get('phenom', '')).strip()
    sig = str(row.get('sig', '')).strip().upper()
    phenom_text = PHENOM_DESCRIPTIONS.get(phenom, phenom)
    sig_text = SIG_DESCRIPTIONS.get(sig, '')
    if phenom_text and sig_text:
        return f"{phenom_text} {sig_text}"
    if phenom_text:
        return phenom_text
    if sig_text:
        return sig_text
    return 'Unknown'


def get_sbw():

    print(f"    LOADING SBW")

    r = requests.get(url, headers={"User-Agent": "my-weather-app"})
    data = r.json()

    features = data["features"]

    gdf = gpd.GeoDataFrame(
        [
            {
                "event": f["properties"]["event"],
                "severity": f["properties"]["severity"],
                "headline": f["properties"]["headline"],
                "geometry": shape(f["geometry"]) if f["geometry"] else None,
            }
            for f in features if f["geometry"] is not None
        ],
        geometry="geometry",
        crs="EPSG:4326"
    )

    return gdf


def get_wwa():
    
    xmin, ymin, xmax, ymax = -127, 23, -65, 51

    print(f"    LOADING WWA")

    wwa_url = (
        "https://mapservices.weather.noaa.gov/eventdriven/rest/services/"
        "WWA/watch_warn_adv/FeatureServer/1/query?"
        f"geometry={xmin},{ymin},{xmax},{ymax}"
        "&geometryType=esriGeometryEnvelope"
        "&inSR=4326"
        "&spatialRel=esriSpatialRelIntersects"
        "&outFields=*"
        "&resultRecordCount=5000"
        "&f=geojson"
    )

    gdf = gpd.read_file(wwa_url)

    # WWA data uses short phenom/sig codes; set a readable event value for coloring.
    # Overwrite raw numeric codes so the legend picks semantic color names.
    gdf["event"] = gdf.apply(_wwa_event_name, axis=1)

    return gdf