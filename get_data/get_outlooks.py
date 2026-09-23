##########################################################
#              NWS / NOAA OUTLOOK GETTER
#  UND Atmospheric Sciences Weather Wall
#  KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2026
##########################################################

import warnings
warnings.filterwarnings("ignore")

# IMPORTS
import time as comp_time
import io
from datetime import datetime, timedelta, timezone
from numbers import Number

import geopandas as geopandas
import pandas as pd
import requests


# REQUEST SESSION
SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "UND-Atmospheric-Sciences-Weather-Wall/1.0"})

# NOAA / NWS ARCGIS SERVICES
SPC_CONVECTIVE = "https://mapservices.weather.noaa.gov/vector/rest/services/outlooks/SPC_wx_outlks/MapServer"
SPC_FIRE = "https://mapservices.weather.noaa.gov/vector/rest/services/fire_weather/SPC_firewx/MapServer"
WPC_ERO  = "https://mapservices.weather.noaa.gov/vector/rest/services/hazards/wpc_precip_hazards/MapServer"
CPC_610  = "https://mapservices.weather.noaa.gov/vector/rest/services/outlooks/cpc_6_10_day_outlk/MapServer"
WPC_WSO  = "https://mapservices.weather.noaa.gov/experimental/rest/services/wpc_winter_storm_outlook/MapServer"
NHC_TROPICAL = "https://mapservices.weather.noaa.gov/tropical/rest/services/tropical/NHC_tropical_weather/MapServer"
NHC_TROPICAL_SUMMARY = "https://mapservices.weather.noaa.gov/tropical/rest/services/tropical/NHC_tropical_weather_summary/MapServer"
USDM_CURRENT = "https://droughtmonitor.unl.edu/data/json/usdm_current.json"



#########################################################################################################
# INTERNAL ARCGIS / TIME HELPERS
#########################################################################################################
def _get_arcgis_layer(service_url, layer_id, where="1=1"):
    query_url = f"{service_url}/{layer_id}/query"
    params = {"where": where, "outFields": "*", "returnGeometry": "true", "outSR": "4326",
        "geometryPrecision": "5", "f": "geojson",}

    response = SESSION.get(query_url, params=params, timeout=30)
    response.raise_for_status()
    data = geopandas.read_file(io.BytesIO(response.content))
    if data is None or data.empty:
        return None
    if data.crs is None:
        data = data.set_crs("EPSG:4326")
    else:
        data = data.to_crs("EPSG:4326")
    return data


def _parse_datetime(value):
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except Exception:
        pass

    if isinstance(value, datetime):
        if value.tzinfo is None:
            return value.replace(tzinfo=timezone.utc)
        return value.astimezone(timezone.utc)

    if isinstance(value, pd.Timestamp):
        if value.tzinfo is None:
            value = value.tz_localize("UTC")
        else:
            value = value.tz_convert("UTC")
        return value.to_pydatetime()

    if isinstance(value, Number):
        try:
            unit = "ms" if abs(float(value)) > 1.0e11 else "s"
            parsed = pd.to_datetime(value, unit=unit, utc=True, errors="coerce")
            if not pd.isna(parsed):
                return parsed.to_pydatetime()
        except Exception:
            pass

    text = str(value).strip()

    for fmt in ["%Y%m%d%H%M", "%Y%m%d%H", "%Y%m%d"]:
        try:
            return datetime.strptime(text, fmt).replace(tzinfo=timezone.utc)
        except Exception:
            pass

    try:
        parsed = pd.to_datetime(text, utc=True, errors="coerce")
        if not pd.isna(parsed):
            return parsed.to_pydatetime()
    except Exception:
        pass

    return None


def _get_time(data, fields, mode="first"):
    if data is None or data.empty:
        return None

    values = []

    for field in fields:
        if field not in data.columns:
            continue

        for value in data[field].dropna().values:
            parsed = _parse_datetime(value)
            if parsed is not None:
                values.append(parsed)

        if values:
            break

    if not values:
        return None

    if mode == "min":
        return min(values)
    elif mode == "max":
        return max(values)
    else:
        return values[0]


def _set_time_attrs(data, issue_fields=None, start_fields=None, end_fields=None):
    if data is None:
        return data

    issue_fields = issue_fields or []
    start_fields = start_fields or []
    end_fields = end_fields or []

    data.attrs["issue_time"] = _get_time(data, issue_fields, mode="max")
    data.attrs["valid_start"] = _get_time(data, start_fields, mode="min")
    data.attrs["valid_end"] = _get_time(data, end_fields, mode="max")

    return data


def _data_is_stale(data, max_age_hours, fields):
    data_time = _get_time(data, fields, mode="max")
    if data_time is None:
        return False
    now_utc = datetime.now(timezone.utc)
    age = now_utc - data_time

    return age > timedelta(hours=max_age_hours)


#########################################################################################################
### SPC CONVECTIVE OUTLOOK ###
#########################################################################################################
def get_spc_convective_outlook(day=1):
    st = comp_time.time()

    print(f'    ACCESSING SPC DAY {day} CONVECTIVE OUTLOOK')

    layer_ids = {1: 1, 2: 9, 3: 17,}

    try:
        outlook = _get_arcgis_layer(SPC_CONVECTIVE, layer_ids[day])
    except Exception as e:
        print(f'    SPC DAY {day} CONVECTIVE OUTLOOK FAILED -- {e}')
        return None

    if outlook is None:
        print(f'    SPC DAY {day} CONVECTIVE OUTLOOK NOT AVAILABLE')
        return None

    outlook = _set_time_attrs(outlook, issue_fields=["issue", "idp_filedate", "idp_ingestdate"], start_fields=["valid"], end_fields=["expire"],)

    elapsed_time = comp_time.time() - st
    print(f'    SPC DAY {day} CONVECTIVE OUTLOOK LOADED.....Time elapsed:',
          comp_time.strftime("%H:%M:%S", comp_time.gmtime(elapsed_time)))

    return outlook
#########################################################################################################


#########################################################################################################
### SPC FIRE WEATHER OUTLOOK ###
#########################################################################################################
def get_spc_fire_outlook(day=1):
    st = comp_time.time()

    print(f'    ACCESSING SPC DAY {day} FIRE WEATHER OUTLOOK')

    layer_ids = {1: [1, 2], 2: [4, 5],}

    parts = []
    for layer_id, risk_type in zip(layer_ids[day], ["wind_rh", "dry_thunder"]):
        try:
            data = _get_arcgis_layer(SPC_FIRE, layer_id)
        except Exception as e:
            print(f'    SPC DAY {day} FIRE LAYER {layer_id} FAILED -- {e}')
            data = None

        if data is not None:
            data["risk_type"] = risk_type
            parts.append(data)

    if not parts:
        print(f'    SPC DAY {day} FIRE WEATHER OUTLOOK NOT AVAILABLE')
        return None

    outlook = geopandas.GeoDataFrame(
        pd.concat(parts, ignore_index=True),
        geometry="geometry",
        crs=parts[0].crs,
    )

    if _data_is_stale(outlook, 48, ["idp_filedate", "idp_ingestdate"]):
        print(f'    SPC DAY {day} FIRE WEATHER OUTLOOK IS STALE -- skipping')
        return None

    outlook = _set_time_attrs(
        outlook,
        issue_fields=["idp_filedate", "idp_ingestdate"],
        start_fields=["valid"],
        end_fields=["expire"],
    )

    elapsed_time = comp_time.time() - st
    print(f'    SPC DAY {day} FIRE WEATHER OUTLOOK LOADED.....Time elapsed:',
          comp_time.strftime("%H:%M:%S", comp_time.gmtime(elapsed_time)))

    return outlook
#########################################################################################################


#########################################################################################################
### WPC EXCESSIVE RAINFALL OUTLOOK ###
#########################################################################################################
def get_wpc_excessive_rainfall_outlook(day=1):
    st = comp_time.time()

    print(f'    ACCESSING WPC DAY {day} EXCESSIVE RAINFALL OUTLOOK')

    try:
        outlook = _get_arcgis_layer(WPC_ERO, day - 1)
    except Exception as e:
        print(f'    WPC DAY {day} EXCESSIVE RAINFALL OUTLOOK FAILED -- {e}')
        return None

    if outlook is None:
        print(f'    WPC DAY {day} EXCESSIVE RAINFALL OUTLOOK NOT AVAILABLE')
        return None

    if _data_is_stale(outlook, 48, ["issue_time", "idp_filedate", "idp_ingestdate"]):
        print(f'    WPC DAY {day} EXCESSIVE RAINFALL OUTLOOK IS STALE -- skipping')
        return None

    outlook = _set_time_attrs(
        outlook,
        issue_fields=["issue_time", "idp_filedate", "idp_ingestdate"],
        start_fields=["start_time"],
        end_fields=["end_time"])

    elapsed_time = comp_time.time() - st
    print(f'    WPC DAY {day} EXCESSIVE RAINFALL OUTLOOK LOADED.....Time elapsed:',
          comp_time.strftime("%H:%M:%S", comp_time.gmtime(elapsed_time)))

    return outlook
#########################################################################################################


#########################################################################################################
### CPC 6-10 DAY TEMPERATURE OUTLOOK ###
#########################################################################################################
def get_cpc_610_temperature_outlook():
    st = comp_time.time()

    print('    ACCESSING CPC 6-10 DAY TEMPERATURE OUTLOOK')

    try:
        outlook = _get_arcgis_layer(CPC_610, 0)
    except Exception as e:
        print(f'    CPC 6-10 DAY TEMPERATURE OUTLOOK FAILED -- {e}')
        return None

    if outlook is None:
        print('    CPC 6-10 DAY TEMPERATURE OUTLOOK NOT AVAILABLE')
        return None

    if _data_is_stale(outlook, 72, ["fcst_date", "idp_filedate", "idp_ingestdate"]):
        print('    CPC 6-10 DAY TEMPERATURE OUTLOOK IS STALE -- skipping')
        return None

    outlook = _set_time_attrs(
        outlook,
        issue_fields=["fcst_date", "idp_filedate", "idp_ingestdate"],
        start_fields=["start_date"],
        end_fields=["end_date"],
    )

    elapsed_time = comp_time.time() - st
    print('    CPC 6-10 DAY TEMPERATURE OUTLOOK LOADED.....Time elapsed:',
          comp_time.strftime("%H:%M:%S", comp_time.gmtime(elapsed_time)))

    return outlook
#########################################################################################################


#########################################################################################################
### CPC 6-10 DAY PRECIPITATION OUTLOOK ###
#########################################################################################################
def get_cpc_610_precipitation_outlook():
    st = comp_time.time()

    print('    ACCESSING CPC 6-10 DAY PRECIPITATION OUTLOOK')

    try:
        outlook = _get_arcgis_layer(CPC_610, 1)
    except Exception as e:
        print(f'    CPC 6-10 DAY PRECIPITATION OUTLOOK FAILED -- {e}')
        return None

    if outlook is None:
        print('    CPC 6-10 DAY PRECIPITATION OUTLOOK NOT AVAILABLE')
        return None

    if _data_is_stale(outlook, 72, ["fcst_date", "idp_filedate", "idp_ingestdate"]):
        print('    CPC 6-10 DAY PRECIPITATION OUTLOOK IS STALE -- skipping')
        return None

    outlook = _set_time_attrs(
        outlook,
        issue_fields=["fcst_date", "idp_filedate", "idp_ingestdate"],
        start_fields=["start_date"],
        end_fields=["end_date"],
    )

    elapsed_time = comp_time.time() - st
    print('    CPC 6-10 DAY PRECIPITATION OUTLOOK LOADED.....Time elapsed:',
          comp_time.strftime("%H:%M:%S", comp_time.gmtime(elapsed_time)))

    return outlook
#########################################################################################################


#########################################################################################################
### WPC WINTER STORM OUTLOOK ###
#########################################################################################################
def get_wpc_winter_storm_outlook(day=1, hazard="snow"):
    st = comp_time.time()

    hazard = str(hazard).lower().replace(" ", "_")

    print(f'    ACCESSING WPC DAY {day} WINTER STORM OUTLOOK ({hazard.upper()})')

    if day not in [1, 2, 3]:
        print('    WPC WINTER STORM OUTLOOK FAILED -- day must be 1, 2, or 3')
        return None

    if hazard not in ["snow", "freezing_rain", "freezingrain", "ice"]:
        print('    WPC WINTER STORM OUTLOOK FAILED -- hazard must be snow or freezing_rain')
        return None

    if hazard == "snow":
        layer_id = day
        hazard_name = "snow"
    else:
        layer_id = 6 + day
        hazard_name = "freezing_rain"

    try:
        outlook = _get_arcgis_layer(WPC_WSO, layer_id)
    except Exception as e:
        print(f'    WPC DAY {day} WINTER STORM OUTLOOK FAILED -- {e}')
        return None

    if outlook is None:
        print(f'    WPC DAY {day} WINTER STORM OUTLOOK NOT AVAILABLE')
        return None

    if _data_is_stale(outlook, 72, ["issue_time", "idp_filedate", "idp_ingestdate"]):
        print(f'    WPC DAY {day} WINTER STORM OUTLOOK IS STALE -- skipping')
        return None

    outlook["hazard_type"] = hazard_name

    outlook = _set_time_attrs(
        outlook,
        issue_fields=["issue_time", "idp_filedate", "idp_ingestdate"],
        start_fields=["valid_time"],
        end_fields=[],
    )

    elapsed_time = comp_time.time() - st
    print(f'    WPC DAY {day} WINTER STORM OUTLOOK LOADED.....Time elapsed:',
          comp_time.strftime("%H:%M:%S", comp_time.gmtime(elapsed_time)))

    return outlook
#########################################################################################################












#########################################################################################################
### NHC 7-DAY TROPICAL WEATHER OUTLOOK + ACTIVE STORMS ###
#########################################################################################################
def get_nhc_7day_tropical_outlook():
    st = comp_time.time()

    print('    ACCESSING NHC 7-DAY TROPICAL OUTLOOK + ACTIVE STORMS')


    # -------------------------------------------------------------------------
    # 7-DAY GRAPHICAL TROPICAL WEATHER OUTLOOK
    #
    # Use the full NHC tropical service here. This is the service that contains
    # the operational GTWO layers for all basins.
    #
    #   2   current disturbance location
    #   3   potential development region
    #   398 development motion
    # -------------------------------------------------------------------------
    outlook_layers = {
        "points": 2,
        "regions": 3,
        "motion": 398,
    }


    # -------------------------------------------------------------------------
    # ACTIVE TROPICAL CYCLONES
    #
    # Use the summary service here because it combines all currently active
    # storms into single forecast point / track / cone layers.
    #
    #   5   forecast points
    #   6   forecast track
    #   7   forecast cone
    # -------------------------------------------------------------------------
    storm_layers = {
        "storm_points": 5,
        "storm_track": 6,
        "storm_cone": 7,
    }


    outlook = {}


    # -------------------------------------------------------------------------
    # GET 7-DAY OUTLOOK DATA
    # -------------------------------------------------------------------------
    for key, layer_id in outlook_layers.items():

        try:
            data = _get_arcgis_layer(
                NHC_TROPICAL,
                layer_id,
            )

        except Exception as e:
            print(f'    NHC {key.upper()} LAYER FAILED -- {e}')
            data = None


        if data is not None and not data.empty:

            # Keep only Atlantic + Eastern Pacific GTWO features.
            #
            # NHC defines the Eastern Pacific basin as east of 140W.
            # Using geometry here is more reliable than depending on varying
            # basin strings in the ArcGIS attributes.
            bounds = data.geometry.bounds

            data = data[
                bounds["maxx"] >= -140.0
            ].copy()


            if data.empty:
                data = None


        outlook[key] = data


        if data is not None and not data.empty:
            print(f'        {key}: {len(data)} feature(s)')

            if "basin" in data.columns:
                print(
                    f'            basins: '
                    f'{data["basin"].dropna().astype(str).unique().tolist()}'
                )

        else:
            print(f'        {key}: no data')


    # -------------------------------------------------------------------------
    # GET ACTIVE STORM DATA
    # -------------------------------------------------------------------------
    for key, layer_id in storm_layers.items():

        try:
            data = _get_arcgis_layer(
                NHC_TROPICAL_SUMMARY,
                layer_id,
            )

        except Exception as e:
            print(f'    NHC {key.upper()} LAYER FAILED -- {e}')
            data = None


        if data is not None and data.empty:
            data = None


        outlook[key] = data


        if data is not None and not data.empty:
            print(f'        {key}: {len(data)} feature(s)')

        else:
            print(f'        {key}: no data')


    # -------------------------------------------------------------------------
    # CHECK FOR AVAILABLE DATA
    # -------------------------------------------------------------------------
    available = [
        data
        for data in outlook.values()
        if data is not None and not data.empty
    ]


    if not available:
        print('    NHC TROPICAL OUTLOOK / ACTIVE STORMS NOT AVAILABLE')
        return None


    # -------------------------------------------------------------------------
    # GET MOST RECENT UPDATE TIME
    # -------------------------------------------------------------------------
    combined_time = None

    for data in available:

        data_time = _get_time(
            data,
            [
                "idp_filedate",
                "idp_ingestdate",
            ],
            mode="max",
        )


        if data_time is not None:

            if combined_time is None or data_time > combined_time:
                combined_time = data_time


    # -------------------------------------------------------------------------
    # CHECK DATA AGE
    # -------------------------------------------------------------------------
    if combined_time is not None:

        age = datetime.now(timezone.utc) - combined_time

        if age > timedelta(hours=24):
            print('    NHC TROPICAL DATA IS STALE -- skipping')
            return None


    outlook["issue_time"] = combined_time


    elapsed_time = comp_time.time() - st

    print(
        '    NHC TROPICAL OUTLOOK / ACTIVE STORMS LOADED.....Time elapsed:',
        comp_time.strftime(
            "%H:%M:%S",
            comp_time.gmtime(elapsed_time),
        ),
    )


    return outlook
#########################################################################################################


#########################################################################################################
### U.S. DROUGHT MONITOR ###
#########################################################################################################
def get_us_drought_monitor():
    st = comp_time.time()

    print('    ACCESSING U.S. DROUGHT MONITOR')

    # Use the NDMC current GeoJSON directly. This is the same source used by
    # the older plotting workflow and avoids the multipart/date-line behavior
    # seen in some ArcGIS representations of the product.
    try:
        outlook = geopandas.read_file(USDM_CURRENT)
    except Exception as e:
        print(f'    U.S. DROUGHT MONITOR FAILED -- {e}')
        return None

    if outlook is None or outlook.empty:
        print('    U.S. DROUGHT MONITOR NOT AVAILABLE')
        return None

    if outlook.crs is None:
        outlook = outlook.set_crs('EPSG:4326')
    else:
        outlook = outlook.to_crs('EPSG:4326')

    # The NDMC GeoJSON normally contains the DM category field. Retain a
    # fallback to the established feature ordering D0-D4 used by the older
    # script in case the field is ever omitted.
    if 'DM' in outlook.columns:
        outlook['DM'] = pd.to_numeric(
            outlook['DM'].astype(str).str.upper().str.replace('D', '', regex=False),
            errors='coerce',
        )
    elif len(outlook) == 5:
        outlook = outlook.copy()
        outlook['DM'] = range(5)
    else:
        print('    U.S. DROUGHT MONITOR HAS UNRECOGNIZED CATEGORY DATA -- skipping')
        return None

    outlook = outlook[outlook['DM'].isin([0, 1, 2, 3, 4])].copy()

    if outlook.empty:
        print('    U.S. DROUGHT MONITOR HAS NO D0-D4 GEOMETRY -- skipping')
        return None

    # Prefer a valid date supplied in the data itself. If none is present,
    # use the file Last-Modified time as the release date and the standard
    # Tuesday-valid / Thursday-release relationship for the map date.
    valid_time = None
    issue_time = None

    for field in ['MapDate', 'mapdate', 'ValidDate', 'valid_date', 'Date', 'date']:
        if field in outlook.columns:
            parsed = pd.to_datetime(outlook[field], errors='coerce', utc=True).dropna()
            if not parsed.empty:
                valid_time = parsed.max().to_pydatetime()
                break

    try:
        response = SESSION.head(USDM_CURRENT, timeout=20, allow_redirects=True)
        response.raise_for_status()
        last_modified = response.headers.get('Last-Modified')

        if last_modified:
            issue_time = pd.to_datetime(last_modified, utc=True).to_pydatetime()
    except Exception:
        pass

    if valid_time is None and issue_time is not None:
        valid_time = issue_time - timedelta(days=2)

    if issue_time is None and valid_time is not None:
        issue_time = valid_time + timedelta(days=2)

    if issue_time is None:
        issue_time = datetime.now(timezone.utc)

    outlook.attrs['issue_time'] = issue_time
    outlook.attrs['valid_start'] = valid_time
    outlook.attrs['valid_end'] = None

    elapsed_time = comp_time.time() - st
    print('    U.S. DROUGHT MONITOR LOADED.....Time elapsed:',
          comp_time.strftime("%H:%M:%S", comp_time.gmtime(elapsed_time)))

    return outlook
#########################################################################################################
