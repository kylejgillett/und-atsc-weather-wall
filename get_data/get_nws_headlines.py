##########################################################
#             NWS HEADLINES LOADING SCRIPT
#  (c) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2026
##########################################################

from concurrent.futures import ThreadPoolExecutor
from time import perf_counter

import geopandas as gpd
import requests

from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


# ============================================================
# CONFIGURATION
# ============================================================

NWS_WWA_SERVER = ("https://mapservices.weather.noaa.gov/eventdriven/rest/services/WWA/watch_warn_adv/FeatureServer")

SBW_LAYER = 0
WWA_LAYER = 1
CONUS_BBOX = (-127.0, 22.0, -65.0, 51.0)
OUT_FIELDS = "prod_type,cap_id"


# ============================================================
# HELPERS
# ============================================================

def _empty_gdf():
    """Return an empty GeoDataFrame with the expected structure."""

    return gpd.GeoDataFrame({"event": [],"cap_id": [],},geometry=[],crs="EPSG:4326")


def _make_session():
    """
    Create a requests Session with automatic retries.

    This protects the weather-wall process from occasional
    transient NOAA/API failures.
    """

    retry = Retry(
        total=3,
        connect=3,
        read=3,
        backoff_factor=0.4,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("GET",),
        respect_retry_after_header=True,
    )

    adapter = HTTPAdapter(
        max_retries=retry,
        pool_connections=2,
        pool_maxsize=2,
    )

    session = requests.Session()

    session.headers.update(
        {
            "User-Agent": (
                "UND-ATSC-Weather-Wall/"
                "1.0 (University of North Dakota)"
            )
        }
    )

    session.mount("https://", adapter)

    return session


def _fetch_layer(layer, name, bbox=CONUS_BBOX):
    """
    Retrieve one NOAA WWA FeatureServer layer as GeoJSON.

    Parameters
    ----------
    layer : int
        ArcGIS layer number.

    name : str
        Human-readable layer name used for logging.

    bbox : tuple
        (xmin, ymin, xmax, ymax)

    Returns
    -------
    geopandas.GeoDataFrame
    """

    start = perf_counter()

    xmin, ymin, xmax, ymax = bbox

    params = {
        "where": "1=1",

        # Restrict request geographically before NOAA sends data
        "geometry": f"{xmin},{ymin},{xmax},{ymax}",
        "geometryType": "esriGeometryEnvelope",
        "inSR": "4326",
        "outSR": "4326",
        "spatialRel": "esriSpatialRelIntersects",

        # Only retrieve what we need
        "outFields": OUT_FIELDS,

        "returnGeometry": "true",
        "returnZ": "false",
        "returnM": "false",

        # More than sufficient for display graphics while
        # reducing unnecessary coordinate precision.
        "geometryPrecision": "5",

        # Server maximum is currently 4000
        "resultRecordCount": "4000",

        # Native GeoJSON output
        "f": "geojson",
    }

    url = f"{NWS_WWA_SERVER}/{layer}/query"

    try:

        with _make_session() as session:

            response = session.get(
                url,
                params=params,

                # (connection timeout, read timeout)
                timeout=(4, 20),
            )

            response.raise_for_status()

            data = response.json()

    except (
        requests.RequestException,
        ValueError,
    ) as exc:

        print(
            f"    WARNING: {name} request failed: {exc}"
        )

        return _empty_gdf()

    # --------------------------------------------------------
    # ArcGIS error handling
    # --------------------------------------------------------

    if "error" in data:

        print(
            f"    WARNING: NOAA returned an error for "
            f"{name}: {data['error']}"
        )

        return _empty_gdf()

    features = data.get("features", [])

    if not features:

        print(
            f"    {name:<25} 0 found "
            f"({perf_counter() - start:.2f}s)"
        )

        return _empty_gdf()

    # --------------------------------------------------------
    # Convert GeoJSON directly to GeoDataFrame
    # --------------------------------------------------------

    gdf = gpd.GeoDataFrame.from_features(
        features,
        crs="EPSG:4326",
    )

    # NOAA already supplies the complete human-readable
    # alert name:
    #
    # prod_type = "Tornado Warning"
    # prod_type = "Winter Weather Advisory"
    # etc.
    gdf = gdf.rename(
        columns={
            "prod_type": "event",
        }
    )

    # Ensure expected columns exist
    for column in ("event", "cap_id"):

        if column not in gdf.columns:
            gdf[column] = None

    gdf = gdf[
        [
            "event",
            "cap_id",
            "geometry",
        ]
    ]

    # Remove empty geometries now rather than making
    # Matplotlib deal with them later.
    gdf = gdf.loc[
        gdf.geometry.notna()
        & ~gdf.geometry.is_empty
    ].copy()

    print(
        f"    {name:<25} "
        f"{len(gdf):>4} found "
        f"({perf_counter() - start:.2f}s)"
    )

    return gdf


# ============================================================
# PUBLIC FUNCTION
# ============================================================

def get_nws_headlines(bbox=CONUS_BBOX):
    """
    Retrieve all NWS polygons required for the weather-wall map.

    Both NOAA layers are downloaded concurrently.

    Returns
    -------
    sbw : GeoDataFrame
        Short-fused / storm-based warnings.

    wwa : GeoDataFrame
        Watches, warnings, and advisories.
    """

    start = perf_counter()

    print("    LOADING NWS HEADLINES")

    # Both requests are independent network operations,
    # so run them at the same time.
    with ThreadPoolExecutor(max_workers=2) as executor:

        sbw_future = executor.submit(
            _fetch_layer,
            SBW_LAYER,
            "Storm-Based Warnings",
            bbox,
        )

        wwa_future = executor.submit(
            _fetch_layer,
            WWA_LAYER,
            "WWA Headlines",
            bbox,
        )

        sbw = sbw_future.result()
        wwa = wwa_future.result()

    # --------------------------------------------------------
    # Protect against duplicate CAP alerts between layers
    # --------------------------------------------------------

    if (
        not sbw.empty
        and not wwa.empty
        and "cap_id" in sbw.columns
        and "cap_id" in wwa.columns
    ):

        sbw_ids = set(
            sbw["cap_id"].dropna()
        )

        if sbw_ids:

            wwa = wwa.loc[
                ~wwa["cap_id"].isin(sbw_ids)
            ].copy()

    print(
        f"    HEADLINES COMPLETE "
        f"({perf_counter() - start:.2f}s)"
    )

    return sbw, wwa

