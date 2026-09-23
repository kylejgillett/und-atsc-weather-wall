##########################################
#
# A SCRIPT TO LOAD HRRR FORECAST DATA
# FOR BASIC SURFACE MAP ANALYSIS
#
# KYLE J GILLETT, UNV. NORTH DAKOTA, 2026
#
##########################################

import warnings
warnings.filterwarnings("ignore")

# IMPORTS

import gc
import os
import sys
import time as comp_time
import tempfile
from datetime import datetime, timedelta, timezone
import requests
import xarray as xr


# NOMADS URLs

FILTER_URL = "https://nomads.ncep.noaa.gov/cgi-bin/filter_hrrr_2d.pl"

DATA_URL = (
    "https://nomads.ncep.noaa.gov/pub/data/nccf/com/hrrr/prod"
)


def _find_latest_hrrr_cycle(max_lookback=8):
    """
    Find the newest HRRR cycle for which forecast hour 18
    is available on NOMADS.
    """

    now = datetime.now(timezone.utc).replace(
        minute=0,
        second=0,
        microsecond=0
    )

    # Start one hour back since the current-hour cycle will
    # generally not be complete yet.
    for lag in range(1, max_lookback + 1):

        cycle = now - timedelta(hours=lag)

        date_string = cycle.strftime("%Y%m%d")
        cycle_hour = cycle.strftime("%H")

        idx_url = (
            f"{DATA_URL}/hrrr.{date_string}/conus/"
            f"hrrr.t{cycle_hour}z.wrfsfcf18.grib2.idx"
        )

        try:

            response = requests.get(
                idx_url,
                timeout=10
            )

            if response.status_code == 200:
                return cycle.replace(tzinfo=None)

        except requests.RequestException:
            continue

    sys.exit(
        "HRRR CYCLE SEARCH FAILED -- "
        "No complete HRRR cycle found during the "
        f"previous {max_lookback} hours."
    )

def _download_hrrr_subset(
    cycle,
    fh,
    center_lat,
    center_lon,
    box_size
):
    """
    Download only the HRRR fields needed for plotting.
    """

    date_string = cycle.strftime("%Y%m%d")
    cycle_hour = cycle.strftime("%H")

    filename = (
        f"hrrr.t{cycle_hour}z."
        f"wrfsfcf{fh:02d}.grib2"
    )

    params = {

        # -------------------------------------------------
        # FILE
        # -------------------------------------------------

        "dir": f"/hrrr.{date_string}/conus",
        "file": filename,

        # -------------------------------------------------
        # VARIABLES
        # -------------------------------------------------

        # Composite simulated reflectivity
        "var_REFC": "on",

        # 2-m temperature / dewpoint
        "var_TMP": "on",
        "var_DPT": "on",

        # 10-m wind
        "var_UGRD": "on",
        "var_VGRD": "on",

        # Total cloud cover
        "var_TCDC": "on",

        # Mean sea-level pressure
        "var_MSLMA": "on",

        # -------------------------------------------------
        # LEVELS
        # -------------------------------------------------

        # REFC + total cloud cover
        "lev_entire_atmosphere": "on",

        # Temperature / dewpoint
        "lev_2_m_above_ground": "on",

        # Wind
        "lev_10_m_above_ground": "on",

        # MSLP
        "lev_mean_sea_level": "on",

        # -------------------------------------------------
        # DOMAIN
        # -------------------------------------------------

        "subregion": "",
        "leftlon": center_lon - box_size,
        "rightlon": center_lon + box_size,
        "toplat": center_lat + box_size,
        "bottomlat": center_lat - box_size,
    }

    try:

        response = requests.get(
            FILTER_URL,
            params=params,
            timeout=120
        )

        response.raise_for_status()

    except requests.RequestException as e:

        raise RuntimeError(
            f"NOMADS DOWNLOAD FAILED -- {e}"
        )

    # NOMADS can occasionally return an HTML error page
    # while still giving an HTTP 200 response.
    if not response.content.startswith(b"GRIB"):

        raise RuntimeError(
            "NOMADS did not return a valid GRIB2 file."
        )

    tmp = tempfile.NamedTemporaryFile(
        suffix=".grib2",
        delete=False
    )

    tmp.write(response.content)
    tmp.close()

    return tmp.name





def _clean_field(field):
    """
    Remove scalar vertical coordinates before fields from
    different HRRR levels are combined into one Dataset.
    """

    field = field.squeeze(drop=True)

    preserve = {
        "latitude",
        "longitude",
        "time",
        "step",
        "valid_time",
    }

    drop_coords = []

    for coord in field.coords:

        if (
            coord not in field.dims
            and coord not in preserve
        ):
            drop_coords.append(coord)

    if drop_coords:
        field = field.drop_vars(
            drop_coords,
            errors="ignore"
        )

    return field

def _open_hrrr_grib(path):
    """
    Open requested HRRR GRIB groups and combine them into
    one simple Xarray Dataset.

    Returned variables
    ------------------
    t2m
        2-m temperature [K]

    d2m
        2-m dewpoint [K]

    u10
        10-m U wind [m/s]

    v10
        10-m V wind [m/s]

    refc
        Composite simulated reflectivity [dBZ]

    cloud_cover
        Total atmospheric cloud cover [%]

    mslp
        Mean sea-level pressure [Pa]
    """

    backend_common = {
        "indexpath": ""
    }


    # =====================================================
    # 2 M TEMPERATURE / DEWPOINT
    # =====================================================

    ds_2m = xr.open_dataset(
        path,
        engine="cfgrib",
        backend_kwargs={
            **backend_common,
            "filter_by_keys": {
                "typeOfLevel": "heightAboveGround",
                "level": 2,
            },
        },
    )


    # =====================================================
    # 10 M WIND
    # =====================================================

    ds_10m = xr.open_dataset(
        path,
        engine="cfgrib",
        backend_kwargs={
            **backend_common,
            "filter_by_keys": {
                "typeOfLevel": "heightAboveGround",
                "level": 10,
            },
        },
    )


    # =====================================================
    # ENTIRE ATMOSPHERE
    #
    # Contains:
    #     REFC = composite reflectivity
    #     TCDC = total cloud cover
    # =====================================================

    ds_atm = xr.open_dataset(
        path,
        engine="cfgrib",
        backend_kwargs={
            **backend_common,
            "filter_by_keys": {
                "typeOfLevel": "atmosphere",
            },
        },
    )


    # =====================================================
    # MEAN SEA LEVEL PRESSURE
    # =====================================================

    ds_msl = xr.open_dataset(
        path,
        engine="cfgrib",
        backend_kwargs={
            **backend_common,
            "filter_by_keys": {
                "typeOfLevel": "meanSea",
            },
        },
    )


    # =====================================================
    # LOAD EVERYTHING INTO MEMORY
    # =====================================================

    ds_2m.load()
    ds_10m.load()
    ds_atm.load()
    ds_msl.load()


    # =====================================================
    # BUILD CLEAN OUTPUT DATASET
    # =====================================================

    data_vars = {}


    # -----------------------------------------------------
    # TEMPERATURE
    # -----------------------------------------------------

    if "t2m" in ds_2m:

        data_vars["t2m"] = _clean_field(
            ds_2m["t2m"]
        )

    else:

        print(
            " WARNING: HRRR 2-m temperature "
            "was not found."
        )


    # -----------------------------------------------------
    # DEWPOINT
    # -----------------------------------------------------

    if "d2m" in ds_2m:

        data_vars["d2m"] = _clean_field(
            ds_2m["d2m"]
        )

    else:

        print(
            " WARNING: HRRR 2-m dewpoint "
            "was not found."
        )


    # -----------------------------------------------------
    # 10-M WIND
    # -----------------------------------------------------

    if "u10" in ds_10m:

        data_vars["u10"] = _clean_field(
            ds_10m["u10"]
        )

    else:

        print(
            " WARNING: HRRR 10-m U wind "
            "was not found."
        )


    if "v10" in ds_10m:

        data_vars["v10"] = _clean_field(
            ds_10m["v10"]
        )

    else:

        print(
            " WARNING: HRRR 10-m V wind "
            "was not found."
        )


    # -----------------------------------------------------
    # COMPOSITE SIMULATED REFLECTIVITY
    # -----------------------------------------------------

    if "refc" in ds_atm:

        data_vars["refc"] = _clean_field(
            ds_atm["refc"]
        )

    else:

        print(
            " WARNING: HRRR composite reflectivity "
            "was not found."
        )


    # -----------------------------------------------------
    # TOTAL CLOUD COVER
    #
    # GRIB TCDC -> cfgrib typically calls this "tcc"
    # Units: %
    # -----------------------------------------------------

    if "tcc" in ds_atm:

        data_vars["cloud_cover"] = _clean_field(
            ds_atm["tcc"]
        )

    else:

        print(
            " WARNING: HRRR total cloud cover "
            "was not found."
        )


    # -----------------------------------------------------
    # MEAN SEA-LEVEL PRESSURE
    #
    # HRRR MSLMA -> cfgrib "mslma"
    # Units: Pa
    # -----------------------------------------------------

    if "mslma" in ds_msl:

        data_vars["mslp"] = _clean_field(
            ds_msl["mslma"]
        )

    else:

        print(
            " WARNING: HRRR MSLMA "
            "was not found."
        )


    # =====================================================
    # FINAL DATASET
    # =====================================================

    forecast_data = xr.Dataset(
        data_vars
    ).load()


    # -----------------------------------------------------
    # METADATA
    # -----------------------------------------------------

    forecast_data.attrs["model"] = "HRRR"

    forecast_data.attrs["temperature_level"] = (
        "2 m AGL"
    )

    forecast_data.attrs["dewpoint_level"] = (
        "2 m AGL"
    )

    forecast_data.attrs["wind_level"] = (
        "10 m AGL"
    )

    forecast_data.attrs["reflectivity"] = (
        "Composite simulated reflectivity"
    )

    forecast_data.attrs["cloud_cover"] = (
        "Total atmospheric cloud cover"
    )

    forecast_data.attrs["mslp"] = (
        "MAPS mean sea-level pressure"
    )


    # =====================================================
    # CLOSE ORIGINAL GRIB DATASETS
    # =====================================================

    ds_2m.close()
    ds_10m.close()
    ds_atm.close()
    ds_msl.close()


    return forecast_data


def hrrr_forecast(
    center_lat=37.86,
    center_lon=-98.61,
    box_size=20,
    forecast_hours=(0,3),#(0, 3, 6, 9, 12, 15, 18),
    request_pause=10,
):


    st = comp_time.time()

    print("ACCESSING HRRR DATA...")

    # -----------------------------------------------------
    # FIND LATEST COMPLETE CYCLE
    # -----------------------------------------------------

    base_time = _find_latest_hrrr_cycle()

    print(
        f' HRRR CYCLE: '
        f'{base_time.strftime("%Y-%m-%d %H:%M:%SZ")}'
    )

    # -----------------------------------------------------
    # RETRIEVE FORECAST HOURS
    # -----------------------------------------------------

    for i, fh in enumerate(forecast_hours):

        valid_time = base_time + timedelta(
            hours=fh
        )

        print(f" GETTING HRRR FH: {fh:02d}")

        grib_path = None
        forecast_data = None

        try:

            # NOMADS asks batch users not to hammer the
            # filtering service with rapid sequential calls.
            if i > 0 and request_pause > 0:
                comp_time.sleep(request_pause)

            grib_path = _download_hrrr_subset(
                cycle=base_time,
                fh=fh,
                center_lat=center_lat,
                center_lon=center_lon,
                box_size=box_size,
            )

            forecast_data = _open_hrrr_grib(
                grib_path
            )

        except Exception as e:

            print(
                f" HRRR FH {fh:02d} FAILED -- {e}"
            )

            if grib_path and os.path.exists(grib_path):
                os.remove(grib_path)

            gc.collect()

            continue

        print(
            f' HRRR FH {fh:02d} COMPLETE: '
            f'{valid_time.strftime("%Y-%m-%d %H:%M:%SZ")}'
        )

        # ---------------------------------------------
        # SEND FORECAST TO PLOTTING SCRIPT
        # ---------------------------------------------

        yield fh, valid_time, forecast_data

        # ---------------------------------------------
        # CLEANUP
        # ---------------------------------------------

        forecast_data.close()

        if (
            grib_path is not None
            and os.path.exists(grib_path)
        ):
            os.remove(grib_path)

        del forecast_data

        gc.collect()

    elapsed_time = comp_time.time() - st

    print(
        "ALL HRRR FORECASTS COMPLETE. Time elapsed:",
        comp_time.strftime(
            "%H:%M:%S",
            comp_time.gmtime(elapsed_time)
        ),
    )