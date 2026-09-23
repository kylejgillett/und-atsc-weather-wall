##########################################################
#            NWS / WPC / CPC OUTLOOK PLOTS
#  UND Atmospheric Sciences Weather Wall
##########################################################

print("############\nSCRIPT RUNNING\n############")

import time as comp_time
st = comp_time.time()

import os
import sys
import warnings
from datetime import datetime, timedelta, timezone
from numbers import Number

warnings.filterwarnings("ignore")

import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import requests

# -----------------------------------------------------------------------------
# Project imports
# -----------------------------------------------------------------------------
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(script_dir, ".."))

if project_root not in sys.path:
    sys.path.append(project_root)

from utils.utils import *
from utils.map import map_builder
from utils.figure import figure_builder


# -----------------------------------------------------------------------------
# User settings
# -----------------------------------------------------------------------------
TARGET_DT = datetime.now(timezone.utc)
CONUS_EXTENT = [-119, -74, 23, 50]
OUTPUT_DIR = "staged_figures/conus_spc_outlooks/"

# Which product families to build.
PLOT_WPC_ERO = True
PLOT_WPC_WSO = True
PLOT_CPC = True
PLOT_SPC_FIRE = True
PLOT_US_HAZARDS = True

# WPC ERO currently has Days 1-5.
ERO_DAYS = [1, 2, 3, 4, 5]

# WPC Winter Storm Outlook has Days 1-4 plus a Days 1-4 maximum.
WSO_DAYS = [1, 2, 3, 4]
WSO_INCLUDE_MAX = True

# CPC seasonal services provide 13 overlapping 3-month leads.
# Lead 1 only is a reasonable Weather Wall default; use range(1, 14) for all.
CPC_SEASONAL_LEADS = [1]

# Fire weather can be expanded to range(1, 9). Keeping D1-D2 by default avoids
# generating a large number of extra products.
SPC_FIRE_DAYS = [1, 2]

# Overall transparency multiplier for NOAA polygon fills.
POLYGON_ALPHA = 0.76

# Requests session
SESSION = requests.Session()
SESSION.headers.update({
    "User-Agent": "UND-Atmospheric-Sciences-Weather-Wall/1.0"
})


# -----------------------------------------------------------------------------
# NOAA/NWS ArcGIS services
# -----------------------------------------------------------------------------
WPC_ERO = (
    "https://mapservices.weather.noaa.gov/vector/rest/services/"
    "hazards/wpc_precip_hazards/MapServer"
)

WPC_WSO = (
    "https://mapservices.weather.noaa.gov/experimental/rest/services/"
    "wpc_winter_storm_outlook/MapServer"
)

CPC_610 = (
    "https://mapservices.weather.noaa.gov/vector/rest/services/"
    "outlooks/cpc_6_10_day_outlk/MapServer"
)

CPC_814 = (
    "https://mapservices.weather.noaa.gov/vector/rest/services/"
    "outlooks/cpc_8_14_day_outlk/MapServer"
)

CPC_MONTHLY_TEMP = (
    "https://mapservices.weather.noaa.gov/vector/rest/services/"
    "outlooks/cpc_mthly_temp_outlk/MapServer"
)

CPC_MONTHLY_PRCP = (
    "https://mapservices.weather.noaa.gov/vector/rest/services/"
    "outlooks/cpc_mthly_precip_outlk/MapServer"
)

CPC_SEASONAL_TEMP = (
    "https://mapservices.weather.noaa.gov/vector/rest/services/"
    "outlooks/cpc_sea_temp_outlk/MapServer"
)

CPC_SEASONAL_PRCP = (
    "https://mapservices.weather.noaa.gov/vector/rest/services/"
    "outlooks/cpc_sea_precip_outlk/MapServer"
)

SPC_FIRE = (
    "https://mapservices.weather.noaa.gov/vector/rest/services/"
    "fire_weather/SPC_firewx/MapServer"
)

US_HAZARDS = (
    "https://mapservices.weather.noaa.gov/vector/rest/services/"
    "hazards/cpc_weather_hazards/MapServer"
)


# -----------------------------------------------------------------------------
# Utility functions
# -----------------------------------------------------------------------------
def _rgba(color, alpha_mult=1.0):
    """Convert an ArcGIS [R,G,B,A] color to a Matplotlib RGBA tuple."""
    if not color or len(color) < 3:
        return (0.0, 0.0, 0.0, 0.0)

    alpha = color[3] / 255.0 if len(color) >= 4 else 1.0
    return (
        color[0] / 255.0,
        color[1] / 255.0,
        color[2] / 255.0,
        max(0.0, min(1.0, alpha * alpha_mult)),
    )


def _parse_datetime(value):
    """Parse common ArcGIS/NWS date representations into UTC datetimes."""
    if value is None:
        return None

    try:
        if pd.isna(value):
            return None
    except Exception:
        pass

    # ArcGIS date fields often arrive as epoch milliseconds.
    if isinstance(value, Number):
        value = float(value)
        if value > 1.0e11:
            ts = pd.to_datetime(value, unit="ms", utc=True, errors="coerce")
        elif value > 1.0e9:
            ts = pd.to_datetime(value, unit="s", utc=True, errors="coerce")
        else:
            return None

        if pd.isna(ts):
            return None
        return ts.to_pydatetime()

    text = str(value).strip()
    if not text:
        return None

    # SPC/fire-weather style compact strings.
    compact_formats = {
        8: "%Y%m%d",
        10: "%Y%m%d%H",
        12: "%Y%m%d%H%M",
        14: "%Y%m%d%H%M%S",
    }
    if text.isdigit() and len(text) in compact_formats:
        try:
            return datetime.strptime(text, compact_formats[len(text)]).replace(
                tzinfo=timezone.utc
            )
        except ValueError:
            pass

    ts = pd.to_datetime(text, utc=True, errors="coerce")
    if pd.isna(ts):
        return None
    return ts.to_pydatetime()


def _column(gdf, *candidates):
    """Find a dataframe column case-insensitively."""
    lookup = {str(c).lower(): c for c in gdf.columns}
    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]
    return None


def _first_value(gdf, *candidates):
    col = _column(gdf, *candidates)
    if col is None:
        return None

    values = gdf[col].dropna()
    if len(values) == 0:
        return None
    return values.iloc[0]


def _first_datetime(gdf, *candidates):
    col = _column(gdf, *candidates)
    if col is None:
        return None

    for value in gdf[col].dropna():
        parsed = _parse_datetime(value)
        if parsed is not None:
            return parsed
    return None


def _freshness_datetime(gdf):
    """Prefer an issue/file timestamp rather than a future valid/end time."""
    candidates = (
        "issue_time",
        "issue",
        "fcst_date",
        "forecast_date",
        "idp_filedate",
        "idp_ingestdate",
    )

    for candidate in candidates:
        col = _column(gdf, candidate)
        if col is None:
            continue

        parsed = [_parse_datetime(v) for v in gdf[col].dropna()]
        parsed = [v for v in parsed if v is not None]
        if parsed:
            return max(parsed)

    return None


def _format_dt(dt, include_time=True):
    if dt is None:
        return None
    if include_time:
        return dt.strftime("%HZ %a %b %d, %Y").upper()
    return dt.strftime("%a %b %d, %Y").upper()


def _build_valid_text(gdf):
    """Construct a compact issue/valid line from whichever fields exist."""
    issue_dt = _first_datetime(
        gdf,
        "issue_time",
        "issue",
        "fcst_date",
        "forecast_date",
        "idp_filedate",
    )

    start_dt = _first_datetime(gdf, "start_time", "start_date", "valid")
    end_dt = _first_datetime(gdf, "end_time", "end_date", "expire")

    # WSO carries a human-readable valid_time field.
    raw_valid = _first_value(gdf, "valid_time")

    parts = []
    if issue_dt is not None:
        parts.append(f"Issued: {_format_dt(issue_dt)}")

    if start_dt is not None and end_dt is not None:
        parts.append(f"Valid: {_format_dt(start_dt)} – {_format_dt(end_dt)}")
    elif raw_valid is not None:
        parts.append(f"Valid: {raw_valid}")
    elif start_dt is not None:
        parts.append(f"Valid: {_format_dt(start_dt)}")

    return "  •  ".join(parts) if parts else None


def _arcgis_value(value):
    """Normalize dataframe values to ArcGIS renderer key strings."""
    if value is None:
        return ""

    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass

    if isinstance(value, Number):
        f = float(value)
        if f.is_integer():
            return str(int(f))
    return str(value)


def fetch_arcgis_layer(service_url, layer_id):
    """Download one NOAA ArcGIS feature layer as GeoJSON + renderer metadata."""
    layer_url = f"{service_url.rstrip('/')}/{layer_id}"

    meta = SESSION.get(
        layer_url,
        params={"f": "pjson"},
        timeout=30,
    )
    meta.raise_for_status()
    metadata = meta.json()

    if "error" in metadata:
        raise RuntimeError(metadata["error"])

    query = SESSION.get(
        f"{layer_url}/query",
        params={
            "where": "1=1",
            "outFields": "*",
            "returnGeometry": "true",
            "outSR": 4326,
            "f": "geojson",
        },
        timeout=60,
    )
    query.raise_for_status()
    payload = query.json()

    features = payload.get("features", [])
    if not features:
        return None, metadata

    gdf = gpd.GeoDataFrame.from_features(features, crs="EPSG:4326")
    gdf = gdf[gdf.geometry.notna() & ~gdf.geometry.is_empty].copy()

    if len(gdf) == 0:
        return None, metadata

    return gdf, metadata


def _symbol_style(symbol):
    if not symbol:
        return None

    outline = symbol.get("outline", {}) or {}
    width = outline.get("width", 0.65)

    try:
        width = float(width)
    except Exception:
        width = 0.65

    return {
        "facecolor": _rgba(symbol.get("color"), POLYGON_ALPHA),
        "edgecolor": _rgba(outline.get("color", [0, 0, 0, 255]), 1.0),
        "linewidth": max(0.35, width),
    }


def _plot_unique_renderer(ax, gdf, renderer):
    fields = [
        renderer.get("field1"),
        renderer.get("field2"),
        renderer.get("field3"),
    ]
    fields = [f for f in fields if f]
    delimiter = renderer.get("fieldDelimiter", ",")

    # ArcGIS field names are effectively case-insensitive, while GeoJSON
    # property casing can vary by service.
    actual_fields = []
    for field in fields:
        actual = _column(gdf, field)
        if actual is None:
            raise KeyError(f"Renderer field '{field}' not found in layer")
        actual_fields.append(actual)

    styles = {}
    ordered_items = []

    for item in renderer.get("uniqueValueInfos", []):
        style = _symbol_style(item.get("symbol"))
        if style is None:
            continue

        key = str(item.get("value", ""))
        styles[key] = style
        ordered_items.append((key, item.get("label", key), style))

    default_style = _symbol_style(renderer.get("defaultSymbol"))

    def make_key(row):
        return delimiter.join(_arcgis_value(row[f]) for f in actual_fields)

    keys = gdf.apply(make_key, axis=1)

    # Draw by class instead of one add_geometries call per polygon.
    for key in pd.unique(keys):
        mask = keys == key
        style = styles.get(str(key), default_style)
        if style is None:
            continue

        geoms = list(gdf.loc[mask, "geometry"])
        if not geoms:
            continue

        ax.add_geometries(
            geoms,
            crs=ccrs.PlateCarree(),
            facecolor=style["facecolor"],
            edgecolor=style["edgecolor"],
            linewidth=style["linewidth"],
            zorder=7,
        )

    # Build a deduplicated categorical legend. Skip completely transparent
    # fills (e.g. the WSO <10% background categories).
    legend_colors = []
    legend_labels = []
    seen = set()

    for _, label, style in ordered_items:
        if label in seen:
            continue
        if style["facecolor"][3] <= 0.01:
            continue
        seen.add(label)
        legend_colors.append(style["facecolor"])
        legend_labels.append(label)

    return legend_colors, legend_labels


def _plot_simple_renderer(ax, gdf, renderer):
    style = _symbol_style(renderer.get("symbol"))
    if style is None:
        return [], []

    ax.add_geometries(
        list(gdf.geometry),
        crs=ccrs.PlateCarree(),
        facecolor=style["facecolor"],
        edgecolor=style["edgecolor"],
        linewidth=style["linewidth"],
        zorder=7,
    )
    return [], []


def plot_arcgis_renderer(ax, gdf, metadata):
    renderer = metadata.get("drawingInfo", {}).get("renderer", {})
    renderer_type = renderer.get("type")

    if renderer_type == "uniqueValue":
        return _plot_unique_renderer(ax, gdf, renderer)

    if renderer_type == "simple":
        return _plot_simple_renderer(ax, gdf, renderer)

    raise NotImplementedError(
        f"Unsupported ArcGIS renderer type: {renderer_type!r} "
        f"for layer {metadata.get('name', 'unknown')}"
    )


def _dedupe_legend(colors, labels):
    out_colors = []
    out_labels = []
    seen = set()

    for color, label in zip(colors, labels):
        if label in seen:
            continue
        seen.add(label)
        out_colors.append(color)
        out_labels.append(label)

    return out_colors, out_labels


def build_product(spec):
    """Fetch, freshness-check, plot, and save one configured outlook product."""
    print(f"\n--- {spec['title']} ---")

    fetched = []

    for layer_id in spec["layers"]:
        try:
            gdf, metadata = fetch_arcgis_layer(spec["service"], layer_id)
        except Exception as exc:
            print(f"    layer {layer_id}: unavailable ({exc})")
            continue

        if gdf is None or len(gdf) == 0:
            print(f"    layer {layer_id}: no features")
            continue

        stamp = _freshness_datetime(gdf)
        max_age = timedelta(hours=spec["max_age_hours"])

        if stamp is not None and TARGET_DT - stamp > max_age:
            age_hours = (TARGET_DT - stamp).total_seconds() / 3600.0
            print(
                f"    layer {layer_id}: stale ({age_hours:.0f} h old; "
                f"latest={stamp:%Y-%m-%d %H:%MZ}) -- skipping"
            )
            continue

        fetched.append((gdf, metadata))
        print(
            f"    layer {layer_id}: loaded {len(gdf)} polygon(s)"
            + (f" | latest={stamp:%Y-%m-%d %H:%MZ}" if stamp else "")
        )

    if not fetched:
        print("    NO CURRENT DATA -- product skipped")
        return False

    fig, ax = map_builder(
        extent=CONUS_EXTENT,
        terrain=True,
        counties=True,
        county_alpha=0.9,
    )

    all_colors = []
    all_labels = []

    for gdf, metadata in fetched:
        try:
            colors, labels = plot_arcgis_renderer(ax, gdf, metadata)
            all_colors.extend(colors)
            all_labels.extend(labels)
        except Exception as exc:
            print(
                f"    renderer failure for {metadata.get('name', 'layer')}: {exc}"
            )

    all_colors, all_labels = _dedupe_legend(all_colors, all_labels)

    # Use the first layer for the displayed issue/valid line. Multi-layer
    # products (e.g. SPC fire weather) should share the same valid period.
    primary_gdf = fetched[0][0]
    valid_text = _build_valid_text(primary_gdf)

    issue_dt = _freshness_datetime(primary_gdf) or TARGET_DT

    # If your current build_filename() supports an explicit order= argument,
    # add it here. This remains compatible with the master-branch helper.
    output_filename = build_filename(
        OUTPUT_DIR,
        "outlook",
        issue_dt,
        variant=spec["variant"],
    )

    figure_builder(
        fig,
        ax,
        title=spec["title"],
        subtitle=spec["subtitle"],
        valid=valid_text,
        category_colors=all_colors if all_colors else None,
        category_labels=all_labels if all_labels else None,
        category_title=spec.get("category_title"),
        footer_left=spec["footer_left"],
        save_path=output_filename,
    )

    plt.close(fig)
    print(f"    saved: {output_filename}")
    return True


# -----------------------------------------------------------------------------
# Product configuration
# -----------------------------------------------------------------------------
products = []

# ---- WPC Excessive Rainfall Outlook -----------------------------------------
if PLOT_WPC_ERO:
    # NOAA layers 0-4 correspond to ERO Days 1-5.
    for day in ERO_DAYS:
        products.append({
            "title": f"Day {day} Excessive Rainfall Outlook",
            "subtitle": "NOAA Weather Prediction Center",
            "service": WPC_ERO,
            "layers": [day - 1],
            "variant": f"ero-d{day}",
            "max_age_hours": 48,
            "category_title": "Excessive Rainfall Risk",
            "footer_left": (
                "NOAA Weather Prediction Center Excessive Rainfall Outlook • "
                "https://www.wpc.ncep.noaa.gov/#page=ero"
            ),
        })


# ---- WPC Winter Storm Outlook -----------------------------------------------
if PLOT_WPC_WSO:
    # Snowfall: layer 1-4 = D1-D4, layer 5 = D1-D4 max.
    for day in WSO_DAYS:
        products.append({
            "title": f"Day {day} Winter Storm Outlook — Snow",
            "subtitle": "NOAA Weather Prediction Center • Probability of Exceeding Warning Criteria",
            "service": WPC_WSO,
            "layers": [day],
            "variant": f"wso-snow-d{day}",
            "max_age_hours": 72,
            "category_title": "Probability",
            "footer_left": (
                "NOAA Weather Prediction Center Experimental Winter Storm Outlook • "
                "https://www.wpc.ncep.noaa.gov/wwd/wso/"
            ),
        })

    # Freezing rain: layer 7-10 = D1-D4, layer 11 = D1-D4 max.
    for day in WSO_DAYS:
        products.append({
            "title": f"Day {day} Winter Storm Outlook — Freezing Rain",
            "subtitle": "NOAA Weather Prediction Center • Probability of Exceeding Warning Criteria",
            "service": WPC_WSO,
            "layers": [6 + day],
            "variant": f"wso-ice-d{day}",
            "max_age_hours": 72,
            "category_title": "Probability",
            "footer_left": (
                "NOAA Weather Prediction Center Experimental Winter Storm Outlook • "
                "https://www.wpc.ncep.noaa.gov/wwd/wso/"
            ),
        })

    if WSO_INCLUDE_MAX:
        products.extend([
            {
                "title": "Days 1–4 Winter Storm Outlook — Snow",
                "subtitle": "NOAA Weather Prediction Center • Maximum Probability of Exceeding Warning Criteria",
                "service": WPC_WSO,
                "layers": [5],
                "variant": "wso-snow-max",
                "max_age_hours": 72,
                "category_title": "Probability",
                "footer_left": (
                    "NOAA Weather Prediction Center Experimental Winter Storm Outlook • "
                    "https://www.wpc.ncep.noaa.gov/wwd/wso/"
                ),
            },
            {
                "title": "Days 1–4 Winter Storm Outlook — Freezing Rain",
                "subtitle": "NOAA Weather Prediction Center • Maximum Probability of Exceeding Warning Criteria",
                "service": WPC_WSO,
                "layers": [11],
                "variant": "wso-ice-max",
                "max_age_hours": 72,
                "category_title": "Probability",
                "footer_left": (
                    "NOAA Weather Prediction Center Experimental Winter Storm Outlook • "
                    "https://www.wpc.ncep.noaa.gov/wwd/wso/"
                ),
            },
        ])


# ---- CPC outlooks ------------------------------------------------------------
if PLOT_CPC:
    cpc_specs = [
        (
            "6–10 Day Temperature Outlook",
            CPC_610,
            0,
            "cpc-610-temp",
            "Temperature Probability",
        ),
        (
            "6–10 Day Precipitation Outlook",
            CPC_610,
            1,
            "cpc-610-prcp",
            "Precipitation Probability",
        ),
        (
            "8–14 Day Temperature Outlook",
            CPC_814,
            0,
            "cpc-814-temp",
            "Temperature Probability",
        ),
        (
            "8–14 Day Precipitation Outlook",
            CPC_814,
            1,
            "cpc-814-prcp",
            "Precipitation Probability",
        ),
        (
            "Monthly Temperature Outlook",
            CPC_MONTHLY_TEMP,
            0,
            "cpc-monthly-temp",
            "Temperature Probability",
        ),
        (
            "Monthly Precipitation Outlook",
            CPC_MONTHLY_PRCP,
            0,
            "cpc-monthly-prcp",
            "Precipitation Probability",
        ),
    ]

    for title, service, layer, variant, legend_title in cpc_specs:
        products.append({
            "title": title,
            "subtitle": "NOAA Climate Prediction Center",
            "service": service,
            "layers": [layer],
            "variant": variant,
            "max_age_hours": 24 * 45 if "Monthly" in title else 96,
            "category_title": legend_title,
            "footer_left": (
                "NOAA Climate Prediction Center Outlook • "
                "https://www.cpc.ncep.noaa.gov/"
            ),
        })

    # Seasonal temperature/precipitation: layer 0 = Lead 1 ... layer 12 = Lead 13.
    for lead in CPC_SEASONAL_LEADS:
        if not 1 <= lead <= 13:
            continue

        layer_id = lead - 1
        products.extend([
            {
                "title": f"Seasonal Temperature Outlook — Lead {lead}",
                "subtitle": "NOAA Climate Prediction Center • 3-Month Outlook",
                "service": CPC_SEASONAL_TEMP,
                "layers": [layer_id],
                "variant": f"cpc-seasonal-temp-l{lead:02d}",
                "max_age_hours": 24 * 45,
                "category_title": "Temperature Probability",
                "footer_left": (
                    "NOAA Climate Prediction Center Seasonal Outlook • "
                    "https://www.cpc.ncep.noaa.gov/"
                ),
            },
            {
                "title": f"Seasonal Precipitation Outlook — Lead {lead}",
                "subtitle": "NOAA Climate Prediction Center • 3-Month Outlook",
                "service": CPC_SEASONAL_PRCP,
                "layers": [layer_id],
                "variant": f"cpc-seasonal-prcp-l{lead:02d}",
                "max_age_hours": 24 * 45,
                "category_title": "Precipitation Probability",
                "footer_left": (
                    "NOAA Climate Prediction Center Seasonal Outlook • "
                    "https://www.cpc.ncep.noaa.gov/"
                ),
            },
        ])


# ---- SPC Fire Weather Outlooks ----------------------------------------------
if PLOT_SPC_FIRE:
    # Each day combines the wind/RH categorical layer and dry-thunder layer.
    fire_layers = {
        1: [1, 2],
        2: [4, 5],
        3: [7, 8],
        4: [10, 11],
        5: [13, 14],
        6: [16, 17],
        7: [19, 20],
        8: [22, 23],
    }

    for day in SPC_FIRE_DAYS:
        if day not in fire_layers:
            continue

        products.append({
            "title": f"Day {day} Fire Weather Outlook",
            "subtitle": "NOAA Storm Prediction Center",
            "service": SPC_FIRE,
            "layers": fire_layers[day],
            "variant": f"fire-d{day}",
            "max_age_hours": 72,
            "category_title": "Fire Weather Risk",
            "footer_left": (
                "NOAA Storm Prediction Center Fire Weather Outlook • "
                "https://www.spc.noaa.gov/products/fire_wx/"
            ),
        })


# ---- U.S. Hazards Outlook ----------------------------------------------------
if PLOT_US_HAZARDS:
    # The service separates temperature, precipitation, and wildfire/drought.
    # Overlay those layers into a single map for each outlook period.
    products.extend([
        {
            "title": "Days 3–7 U.S. Hazards Outlook",
            "subtitle": "NOAA Weather Prediction Center",
            "service": US_HAZARDS,
            "layers": [1, 4, 7],
            "variant": "hazards-d3-7",
            # Product is weekday-only, so allow a weekend/holiday buffer.
            "max_age_hours": 120,
            "category_title": "Hazard",
            "footer_left": (
                "NOAA/NWS U.S. Hazards Outlook • "
                "https://www.wpc.ncep.noaa.gov/threats/threats.php"
            ),
        },
        {
            "title": "Days 8–14 U.S. Hazards Outlook",
            "subtitle": "NOAA Climate Prediction Center",
            "service": US_HAZARDS,
            "layers": [3, 6, 8],
            "variant": "hazards-d8-14",
            "max_age_hours": 120,
            "category_title": "Hazard",
            "footer_left": (
                "NOAA/NWS U.S. Hazards Outlook • "
                "https://www.cpc.ncep.noaa.gov/products/predictions/threats/"
            ),
        },
    ])


# -----------------------------------------------------------------------------
# Build all configured products. Missing/stale products are skipped individually.
# -----------------------------------------------------------------------------
built = 0
skipped = 0

for product in products:
    try:
        if build_product(product):
            built += 1
        else:
            skipped += 1
    except Exception as exc:
        skipped += 1
        print(f"    FAILED -- {product['variant']}: {exc}")


elapsed_time = comp_time.time() - st
print(
    "############\n"
    f"SCRIPT FINISHED: built={built}, skipped={skipped}, "
    f"time={comp_time.strftime('%H:%M:%S', comp_time.gmtime(elapsed_time))}\n"
    "############"
)
