##########################################################
#       GFK NWS / GEFS FORECAST DATA GETTER
#  UNIVERSITY OF NORTH DAKOTA - ATMOSPHERIC SCIENCES
##########################################################

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timedelta, timezone
from pathlib import Path
import tempfile
import time as comp_time

import numpy as np
import pandas as pd
import requests
import cfgrib
import xarray as xr
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


GFK_LAT = 47.9493
GFK_LON = -97.1761
NWS_API = "https://api.weather.gov"
GEFS_FILTER = "https://nomads.ncep.noaa.gov/cgi-bin/filter_gefs_atmos_0p25s.pl"
GEFS_DATA = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/gens/prod"


# ============================================================
# SHARED REQUEST SESSION
# ============================================================

def _session(user_agent: str) -> requests.Session:
    retry = Retry(
        total=3,
        connect=3,
        read=3,
        backoff_factor=0.5,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("GET",),
        respect_retry_after_header=True,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session = requests.Session()
    session.headers.update({"User-Agent": user_agent})
    session.mount("https://", adapter)
    return session


# ============================================================
# NWS HOURLY FORECAST
# ============================================================

def get_nws_hourly_forecast(lat: float = GFK_LAT, lon: float = GFK_LON) -> pd.DataFrame:
    """Return hourly NWS temperature, apparent temperature, dewpoint, and PoP."""
    print("    ACCESSING NWS HOURLY FORECAST...")
    session = _session("UND-ATSC-Weather-Wall (weather data display)")
    session.headers.update({"Accept": "application/geo+json"})

    # Resolve the NWS forecast grid for the requested point.
    points = session.get(f"{NWS_API}/points/{lat:.4f},{lon:.4f}", timeout=20)
    points.raise_for_status()
    hourly_url = points.json()["properties"]["forecastHourly"]
    response = session.get(hourly_url, timeout=30)
    response.raise_for_status()

    # Convert returned hourly periods to a compact UTC dataframe. NWS exposes
    # dewpoint and wind chill as Celsius quantitative values, while the
    # period temperature is in the temperatureUnit field (normally °F).
    rows = response.json()["properties"]["periods"]
    records = []
    for period in rows:
        temperature = period.get("temperature")
        unit = period.get("temperatureUnit", "F")
        if temperature is None:
            continue
        if unit == "C":
            temperature = temperature * 9 / 5 + 32

        def quantitative_f(name):
            value = period.get(name)
            if not isinstance(value, dict):
                return np.nan
            value = value.get("value")
            if value is None:
                return np.nan
            return float(value) * 9 / 5 + 32

        pop = period.get("probabilityOfPrecipitation") or {}
        pop_value = pop.get("value") if isinstance(pop, dict) else None
        short_forecast = period.get("shortForecast", "")
        description = short_forecast.lower()
        if any(word in description for word in ("freezing rain", "ice pellets", "sleet", "wintry mix", "mixed precipitation")):
            precip_type = "Wintry Mix"
        elif "snow" in description or "flurr" in description:
            precip_type = "Snow"
        elif "thunder" in description:
            precip_type = "Thunderstorms"
        elif any(word in description for word in ("rain", "shower", "drizzle")):
            precip_type = "Rain"
        else:
            precip_type = "None"

        records.append({
            "valid_time": pd.Timestamp(period["startTime"]).tz_convert("UTC"),
            "temperature_f": float(temperature),
            "dewpoint_f": quantitative_f("dewpoint"),
            "wind_chill_f": quantitative_f("windChill"),
            "pop_percent": float(pop_value) if pop_value is not None else 0.0,
            "precip_type": precip_type,
        })

    if not records:
        raise RuntimeError("NWS returned no hourly temperatures for GFK.")
    forecast = pd.DataFrame(records).drop_duplicates("valid_time").sort_values("valid_time").reset_index(drop=True)
    print(f"    NWS HOURLY FORECAST COMPLETE: {len(forecast)} PERIODS")
    return forecast


def get_nws_qpf_forecast(lat: float = GFK_LAT, lon: float = GFK_LON) -> pd.DataFrame:
    """Return NWS gridded QPF intervals, converted to inches.

    The NWS hourly forecast does not include quantitative precipitation
    amounts. Those are provided by the point's forecastGridData endpoint.
    Each returned valid_time marks the end of its QPF accumulation interval.
    """
    print("    ACCESSING NWS QUANTITATIVE PRECIPITATION FORECAST...")
    session = _session("UND-ATSC-Weather-Wall (weather data display)")
    session.headers.update({"Accept": "application/geo+json"})
    try:
        point_response = session.get(f"{NWS_API}/points/{lat:.4f},{lon:.4f}", timeout=20)
        point_response.raise_for_status()
        grid_url = point_response.json()["properties"]["forecastGridData"]
        grid_response = session.get(grid_url, timeout=30)
        grid_response.raise_for_status()
        qpf = grid_response.json()["properties"].get("quantitativePrecipitation", {})
    except (requests.RequestException, KeyError, TypeError, ValueError) as exc:
        print(f"    NWS QPF UNAVAILABLE: {exc}")
        return pd.DataFrame(columns=["valid_time", "precipitation_in"])

    unit = str(qpf.get("uom", "")).lower()
    records = []
    for item in qpf.get("values", []):
        try:
            interval_start, duration = item["validTime"].split("/")
            valid_time = pd.Timestamp(interval_start) + pd.Timedelta(duration)
            value = item.get("value")
            if value is None:
                continue
            value = float(value)
            if "mm" in unit:
                value /= 25.4
            elif "inch" in unit or unit.endswith(":in"):
                pass
            else:
                print(f"    NWS QPF UNAVAILABLE: unsupported unit {unit!r}")
                return pd.DataFrame(columns=["valid_time", "precipitation_in"])
            records.append({
                "valid_time": valid_time.tz_convert("UTC"),
                "precipitation_in": value,
            })
        except (KeyError, TypeError, ValueError):
            continue

    if not records:
        print("    NWS QPF UNAVAILABLE: no valid accumulation intervals returned")
        return pd.DataFrame(columns=["valid_time", "precipitation_in"])

    forecast = (
        pd.DataFrame(records)
        .drop_duplicates("valid_time")
        .sort_values("valid_time")
        .reset_index(drop=True)
    )
    print(f"    NWS QPF FORECAST COMPLETE: {len(forecast)} ACCUMULATION INTERVALS")
    return forecast


# ============================================================
# GEFS CYCLE LOOKUP
# ============================================================


def _latest_gefs_cycle(session: requests.Session, lookback_hours: int = 72) -> datetime:
    """Find a GEFS cycle with an available ensemble-mean lead-zero file."""
    now = datetime.now(timezone.utc).replace(minute=0, second=0, microsecond=0)
    # GEFS is initialized every six hours. Begin with the latest cycle
    # expected to have its lead-zero ensemble-mean file available.
    cycle = now.replace(hour=(now.hour // 6) * 6) - timedelta(hours=6)
    lookup_results = []
    for lag in range(0, lookback_hours + 1, 6):
        candidate = cycle - timedelta(hours=lag)
        relative_dir = f"gefs.{candidate:%Y%m%d}/{candidate:%H}/atmos/pgrb2sp25"
        filename = f"geavg.t{candidate:%H}z.pgrb2s.0p25.f000"
        url = f"{GEFS_DATA}/{relative_dir}/{filename}.idx"
        try:
            response = session.get(url, timeout=10)
            if response.status_code == 200:
                print(f"    GEFS CYCLE FOUND: {candidate:%Y-%m-%d %HZ}")
                return candidate
            response_text = response.text.lower()
            if "over rate limit" in response_text:
                raise RuntimeError(
                    "NOMADS is rate-limiting GEFS requests. This build has not "
                    "downloaded any ensemble mean/spread data yet; wait before retrying."
                )
            result = f"{candidate:%m-%d %HZ}: HTTP {response.status_code}"
            if response.history:
                result += f" after redirect to {response.url}"
            location = response.headers.get("Location")
            if location:
                result += f" (Location: {location})"
            elif response.status_code in (301, 302, 303, 307, 308):
                body_excerpt = " ".join(response.text.split())[:180]
                result += f" (final URL: {response.url}"
                if body_excerpt:
                    result += f"; response: {body_excerpt}"
                result += ")"
            lookup_results.append(result)
        except requests.RequestException as exc:
            lookup_results.append(f"{candidate:%m-%d %HZ}: {type(exc).__name__}: {exc}")
            continue
    # Candidates are checked newest first, so the first entries are the most
    # useful when explaining why cycle discovery failed.
    recent_results = "; ".join(lookup_results[:4])
    raise RuntimeError(
        f"Could not find an available GEFS cycle in the previous {lookback_hours} hours. "
        f"Latest checks: {recent_results or 'no responses received'}."
    )


# ============================================================
# GEFS ENSEMBLE MEAN / SPREAD DATA
# ============================================================

def _read_gefs_point(
    session: requests.Session,
    cycle: datetime,
    lead_hour: int,
    product: str,
    lat: float,
    lon: float,
) -> dict:
    # Build one NOMADS request for the ensemble mean or spread at a lead.
    directory = f"/gefs.{cycle:%Y%m%d}/{cycle:%H}/atmos/pgrb2sp25"
    filename = f"{product}.t{cycle:%H}z.pgrb2s.0p25.f{lead_hour:03d}"
    params = {
        "dir": directory,
        "file": filename,
        "var_TMP": "on",
        "var_APCP": "on",
        "var_CRAIN": "on",
        "var_CSNOW": "on",
        "var_CICEP": "on",
        "var_CFRZR": "on",
        "lev_2_m_above_ground": "on",
        "lev_surface": "on",
        "subregion": "on",
        "leftlon": f"{lon - 0.3:.2f}",
        "rightlon": f"{lon + 0.3:.2f}",
        "toplat": f"{lat + 0.3:.2f}",
        "bottomlat": f"{lat - 0.3:.2f}",
    }
    response = session.get(GEFS_FILTER, params=params, timeout=45)
    if "over rate limit" in response.text[:1000].lower():
        raise RuntimeError(
            "NOMADS rate-limited a GEFS mean/spread request. Wait before retrying."
        )
    response.raise_for_status()
    if not response.content.startswith(b"GRIB"):
        raise RuntimeError(f"NOMADS did not return GRIB data for {filename}.")

    with tempfile.TemporaryDirectory(prefix="weatherwall_gefs_") as temp_dir:
        grib_path = Path(temp_dir) / "point.grib2"
        grib_path.write_bytes(response.content)
        datasets = cfgrib.open_datasets(
            grib_path,
            engine="cfgrib",
            backend_kwargs={"indexpath": ""},
        )
        try:
            # cfgrib separates different level types into separate datasets.
            fields = {}
            for dataset in datasets:
                for name in dataset.data_vars:
                    fields.setdefault(name.lower(), (dataset, name))

            def point_value(candidates):
                for candidate in candidates:
                    if candidate not in fields:
                        continue
                    dataset, name = fields[candidate]
                    field = dataset[name]
                    latitude = np.asarray(dataset["latitude"].values)
                    longitude = np.asarray(dataset["longitude"].values)
                    if latitude.ndim == 1 and longitude.ndim == 1:
                        longitude, latitude = np.meshgrid(longitude, latitude)
                    target_lon = lon % 360 if np.nanmax(longitude) > 180 else lon
                    nearest = np.unravel_index(
                        np.nanargmin(
                            (latitude - lat) ** 2
                            + ((longitude - target_lon) * np.cos(np.deg2rad(lat))) ** 2
                        ),
                        latitude.shape,
                    )
                    value = np.asarray(field.values)[nearest]
                    return float(value)
                return np.nan

            # cfgrib commonly names NCEP's 2 m TMP field ``2t``; some
            # backends expose the same field as ``t2m`` or ``t``.
            kelvin = point_value(("2t", "t2m", "t"))
            if not np.isfinite(kelvin):
                available = ", ".join(sorted(fields)) or "none"
                raise RuntimeError(
                    f"No 2 m temperature field found in {filename}; "
                    f"cfgrib returned: {available}."
                )

            # GEFS APCP is an interval accumulation (or its ensemble spread).
            # It is expressed in kg m-2, numerically equivalent to mm liquid
            # water; convert to inches.
            precip_mm = 0.0 if lead_hour == 0 else point_value(("tp", "apcp"))
            if not np.isfinite(precip_mm):
                precip_mm = np.nan
            precip_in = precip_mm / 25.4 if np.isfinite(precip_mm) else np.nan
            rain = point_value(("crain",))
            snow = point_value(("csnow",))
            ice = point_value(("cicep",))
            freezing_rain = point_value(("cfrzr",))
            temperature_f = (
                kelvin * 9 / 5
                if product == "gespr"
                else (kelvin - 273.15) * 9 / 5 + 32
            )
            return {
                "temperature_f": temperature_f,
                "precipitation_in": precip_in,
                "rain": rain,
                "snow": snow,
                "ice_pellets": ice,
                "freezing_rain": freezing_rain,
            }
        finally:
            for dataset in datasets:
                dataset.close()


def get_gefs_temperature_ensemble(
    lat: float = GFK_LAT,
    lon: float = GFK_LON,
    *,
    forecast_hours=range(0, 181, 6),
    max_workers: int = 2,
) -> tuple[pd.DataFrame, datetime]:
    """Return GEFS ensemble mean and spread for temperature and precipitation.

    NCEP's ``geavg`` and ``gespr`` products replace separate requests for all
    31 members. This reduces the normal 31-member, 31-lead retrieval from 961
    GRIB subsets to 62. Spread fields are standard deviations across members.
    """
    start_time = comp_time.time()
    print("    ACCESSING GEFS ENSEMBLE MEAN / SPREAD...")
    session = _session("UND-ATSC-Weather-Wall (weather data display)")
    cycle = _latest_gefs_cycle(session)
    forecast_hours = tuple(int(hour) for hour in forecast_hours)
    products = ("geavg", "gespr")
    work = [(hour, product) for hour in forecast_hours for product in products]
    data_by_lead = {}
    failed = 0
    completed = 0
    failure_samples = []
    print(
        f"    REQUESTING GEFS MEAN + SPREAD AT {len(forecast_hours)} LEADS "
        f"({len(work)} SUBSETS; 2 PER LEAD)..."
    )

    def fetch(item):
        hour, product = item
        return hour, product, _read_gefs_point(session, cycle, hour, product, lat, lon)

    with ThreadPoolExecutor(max_workers=max(1, int(max_workers))) as executor:
        futures = {executor.submit(fetch, item): item for item in work}
        for future in as_completed(futures):
            completed += 1
            try:
                hour, product, values = future.result()
                data_by_lead.setdefault(hour, {})[product] = values
            except Exception as exc:
                failed += 1
                hour, product = futures[future]
                if len(failure_samples) < 5:
                    failure_samples.append(f"{product} FH{hour:03d}: {exc}")

            # Report occasional progress without printing every successful file.
            interval = max(1, len(work) // 10)
            if completed % interval == 0 or completed == len(work):
                print(f"    GEFS SUBSETS COMPLETE: {completed}/{len(work)}")

    records = []
    for hour in forecast_hours:
        products_at_lead = data_by_lead.get(hour, {})
        mean_values = products_at_lead.get("geavg")
        spread_values = products_at_lead.get("gespr")
        if mean_values is None:
            continue

        record = {
            "valid_time": pd.Timestamp(cycle + timedelta(hours=hour)),
            "forecast_hour": hour,
            "temperature_f": mean_values["temperature_f"],
            "temperature_spread_f": spread_values["temperature_f"] if spread_values is not None else np.nan,
            "precipitation_in": mean_values["precipitation_in"],
            "precipitation_spread_in": (
                spread_values["precipitation_in"]
                if spread_values is not None else np.nan
            ),
            "rain": mean_values["rain"],
            "snow": mean_values["snow"],
            "ice_pellets": mean_values["ice_pellets"],
            "freezing_rain": mean_values["freezing_rain"],
        }
        records.append(record)

    if not records:
        for failure in failure_samples:
            print(f"    GEFS SUBSET FAILED: {failure}")
        raise RuntimeError("No GEFS ensemble-mean temperatures could be retrieved.")
    ensemble = pd.DataFrame(records).sort_values("valid_time").reset_index(drop=True)
    elapsed_time = comp_time.time() - start_time
    for failure in failure_samples:
        print(f"    GEFS SUBSET WARNING: {failure}")
    print(
        f"    GEFS MEAN/SPREAD COMPLETE: {len(ensemble)} LEADS, {failed} FAILED "
        f"({comp_time.strftime('%H:%M:%S', comp_time.gmtime(elapsed_time))})"
    )
    return ensemble, cycle
