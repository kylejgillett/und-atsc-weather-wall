"""
Consistent Cartopy basemaps for the UND ATSC Weather Wall.

Recommended use:
    fig, ax = build_map(extent=extent, style="weather")
    fig, ax = build_map(extent=extent, style="hazards")
    fig, ax = build_map(extent=extent, style="satellite")

Always specify transform=ccrs.PlateCarree() when plotting lat/lon weather data.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Literal, Sequence

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from cartopy.feature import ShapelyFeature
from cartopy.io import img_tiles
from cartopy.io.shapereader import Reader
from metpy.plots import USCOUNTIES


DATA_CRS = ccrs.PlateCarree()
CONUS_EXTENT = (-122.0, -73.0, 21.0, 56.0)

StyleName = Literal["weather", "hazards", "dark", "satellite"]
DetailName = Literal["auto", "conus", "regional", "local"]

UTILS_DIR = Path(__file__).resolve().parent
TILE_CACHE_DIR = UTILS_DIR / ".cartopy_tile_cache"
TERRAIN_CACHE_DIR = UTILS_DIR / ".cartopy_terrain_cache"

ROAD_FILES = {
    "us_primary": UTILS_DIR / "tl_2024_us_primaryroads" / "tl_2024_us_primaryroads.shp",
    "nd": UTILS_DIR / "tl_2024_38_prisecroads" / "tl_2024_38_prisecroads.shp",
    "mn": UTILS_DIR / "tl_2024_27_prisecroads" / "tl_2024_27_prisecroads.shp",
}


# Quiet basemaps: weather data should be the loudest thing on the figure.
MAP_STYLES = {
    # Signature weather-wall look:
    # dark charcoal terrain, blue-black water, silver geography.
    "weather": {
        "figure": "#091119",
        "land": "#1A2229",
        "water": "#071018",
        "state": "#AAB4BC",
        "border": "#C6CDD2",
        "coast": "#C6CDD2",
        "county": "#66727C",
        "road": "#58646D",
        "state_alpha": 0.78,
        "boundary_alpha": 0.88,
        "county_alpha": 0.28,
        "road_alpha": 0.22,
    },

    # Similar identity, but a touch flatter/quieter for categorical hazards.
    "hazards": {
        "figure": "#0B1218",
        "land": "#1D252B",
        "water": "#081017",
        "state": "#9CA6AD",
        "border": "#BCC4C9",
        "coast": "#BCC4C9",
        "county": "#626D75",
        "road": "#59636B",
        "state_alpha": 0.70,
        "boundary_alpha": 0.80,
        "county_alpha": 0.23,
        "road_alpha": 0.18,
    },

    # Higher-contrast version for radar and luminous fields.
    "dark": {
        "figure": "#050B10",
        "land": "#141C22",
        "water": "#040A0F",
        "state": "#BBC4CA",
        "border": "#D9DEE1",
        "coast": "#D9DEE1",
        "county": "#69757E",
        "road": "#606A72",
        "state_alpha": 0.82,
        "boundary_alpha": 0.92,
        "county_alpha": 0.30,
        "road_alpha": 0.24,
    },

    # Satellite keeps the imagery itself dominant while matching boundary style.
    "satellite": {
        "figure": "#050B10",
        "land": None,
        "water": None,
        "state": "#CBD3D8",
        "border": "#E4E8EA",
        "coast": "#E4E8EA",
        "county": "#8C989F",
        "road": "#8A949A",
        "state_alpha": 0.76,
        "boundary_alpha": 0.88,
        "county_alpha": 0.30,
        "road_alpha": 0.24,
    },
}


DETAIL_STYLES = {
    "conus": {
        "scale": "50m",
        "county_scale": "20m",
        "state_lw": 0.80,
        "boundary_lw": 0.90,
        "county_lw": 0.35,
        "road_lw": 0.40,
    },
    "regional": {
        "scale": "50m",
        "county_scale": "20m",
        "state_lw": 0.90,
        "boundary_lw": 1.00,
        "county_lw": 0.42,
        "road_lw": 0.50,
    },
    "local": {
        "scale": "10m",
        "county_scale": "5m",
        "state_lw": 1.05,
        "boundary_lw": 1.10,
        "county_lw": 0.55,
        "road_lw": 0.60,
    },
}


def _validate_extent(extent: Sequence[float]) -> tuple[float, float, float, float]:
    if len(extent) != 4:
        raise ValueError("extent must be [west, east, south, north].")

    west, east, south, north = map(float, extent)

    if east <= west or north <= south:
        raise ValueError("Invalid map extent.")

    return west, east, south, north


def _auto_detail(extent) -> Literal["conus", "regional", "local"]:
    west, east, south, north = extent
    span = max(east - west, north - south)

    if span > 25:
        return "conus"
    if span > 7:
        return "regional"
    return "local"


def _lambert_for_extent(extent):
    west, east, south, north = extent
    center_lon = (west + east) / 2
    center_lat = (south + north) / 2

    lower = max(20.0, min(35.0, south + 4.0))
    upper = min(60.0, max(45.0, north - 4.0))

    if upper <= lower:
        lower, upper = 30.0, 50.0

    return ccrs.LambertConformal(
        central_longitude=center_lon,
        central_latitude=center_lat,
        standard_parallels=(lower, upper),
    )


def _auto_tile_zoom(extent) -> int:
    west, east, south, north = extent
    span = max(east - west, north - south)

    if span > 35:
        return 4
    if span > 18:
        return 5
    if span > 9:
        return 6
    if span > 4.5:
        return 7
    if span > 2:
        return 8
    return 9


def _auto_terrain_zoom(extent) -> int:
    """Use slightly less detail than imagery; enough for visible relief."""
    return max(4, min(8, _auto_tile_zoom(extent) - 1))


def _terrain_tiles():
    """
    Label-free shaded-relief tiles.

    Using a custom tile URL avoids roads/place labels in the terrain layer,
    which keeps the weather-wall map clean.
    """
    TERRAIN_CACHE_DIR.mkdir(parents=True, exist_ok=True)

    return img_tiles.GoogleTiles(
        url=(
            "https://server.arcgisonline.com/ArcGIS/rest/services/"
            "World_Shaded_Relief/MapServer/tile/{z}/{y}/{x}.jpg"
        ),
        cache=TERRAIN_CACHE_DIR,
    )


@lru_cache(maxsize=8)
def _read_geometries(path_string: str):
    path = Path(path_string)
    if not path.exists():
        return tuple()
    return tuple(Reader(str(path)).geometries())


def _add_line_shapefile(
    ax,
    path: Path,
    *,
    color: str,
    linewidth: float,
    alpha: float,
    zorder: float,
):
    geoms = _read_geometries(str(path))
    if not geoms:
        return False

    feature = ShapelyFeature(
        geoms,
        DATA_CRS,
        facecolor="none",
        edgecolor=color,
        linewidth=linewidth,
        alpha=alpha,
    )
    ax.add_feature(feature, zorder=zorder)
    return True


def _add_roads(ax, detail, style_cfg, detail_cfg):
    if detail == "conus":
        return

    kwargs = {
        "color": style_cfg["road"],
        "linewidth": detail_cfg["road_lw"],
        "alpha": style_cfg["road_alpha"],
        "zorder": 9.6,
    }

    if detail == "regional":
        _add_line_shapefile(ax, ROAD_FILES["us_primary"], **kwargs)
        return

    # Current detailed local road files in the repository.
    added = False
    for key in ("nd", "mn"):
        added |= _add_line_shapefile(ax, ROAD_FILES[key], **kwargs)

    if not added:
        _add_line_shapefile(ax, ROAD_FILES["us_primary"], **kwargs)


def _add_flat_background(ax, style_cfg, scale):
    ax.set_facecolor(style_cfg["figure"])

    for feature, color, zorder in (
        (cfeature.LAND, style_cfg["land"], 0.0),
        (cfeature.OCEAN, style_cfg["water"], 0.0),
        (cfeature.LAKES, style_cfg["water"], 0.1),
    ):
        ax.add_feature(
            feature.with_scale(scale),
            facecolor=color,
            edgecolor="none",
            zorder=zorder,
        )


def _add_boundaries(ax, style_cfg, detail_cfg, counties):
    scale = detail_cfg["scale"]

    ax.add_feature(
        cfeature.STATES.with_scale(scale),
        facecolor="none",
        edgecolor=style_cfg["state"],
        linewidth=detail_cfg["state_lw"],
        alpha=style_cfg["state_alpha"],
        zorder=10.0,
    )

    for feature, color, zorder in (
        (cfeature.BORDERS, style_cfg["border"], 10.1),
        (cfeature.COASTLINE, style_cfg["coast"], 10.2),
    ):
        ax.add_feature(
            feature.with_scale(scale),
            facecolor="none",
            edgecolor=color,
            linewidth=detail_cfg["boundary_lw"],
            alpha=style_cfg["boundary_alpha"],
            zorder=zorder,
        )

    if counties:
        ax.add_feature(
            USCOUNTIES.with_scale(detail_cfg["county_scale"]),
            facecolor="none",
            edgecolor=style_cfg["county"],
            linewidth=detail_cfg["county_lw"],
            alpha=style_cfg["county_alpha"],
            zorder=9.0,
        )


def build_map(
    extent: Sequence[float] = CONUS_EXTENT,
    *,
    style: StyleName = "weather",
    detail: DetailName = "auto",
    projection=None,
    counties: bool | None = None,
    roads: bool | None = None,
    satellite_zoom: int | None = None,
    satellite_alpha: float = 0.72,
    terrain: bool | None = None,
    terrain_zoom: int | None = None,
    terrain_alpha: float = 0.42,
    figsize: tuple[float, float] = (20, 10),
    dpi: int = 250,
):
    """
    Create a consistent weather-wall map.

    Parameters
    ----------
    style
        "weather"   : default model/analysis basemap.
        "hazards"   : extra-neutral background for warnings/SPC polygons.
        "dark"      : dark flat basemap.
        "satellite" : cached, slightly muted Google satellite imagery.

    detail
        "auto", "conus", "regional", or "local".

    counties
        None = automatically show only on local maps.

    roads
        None = automatically show only on local satellite maps.

    terrain
        Add a label-free shaded-relief underlay. None enables it automatically
        for the main "weather" style and disables it for the other styles.

    terrain_zoom
        Tile zoom for the relief layer. None selects a conservative level
        automatically from the map extent.

    terrain_alpha
        Relief opacity. The default is intentionally subtle so filled weather
        fields remain visually dominant.

    projection
        None = centered Lambert for flat maps; tile-native Mercator for
        satellite maps.
    """
    extent = _validate_extent(extent)

    if style not in MAP_STYLES:
        raise ValueError(f"style must be one of {tuple(MAP_STYLES)}")

    if detail == "auto":
        detail = _auto_detail(extent)
    elif detail not in DETAIL_STYLES:
        raise ValueError("detail must be 'auto', 'conus', 'regional', or 'local'")

    style_cfg = MAP_STYLES[style]
    detail_cfg = DETAIL_STYLES[detail]

    if counties is None:
        counties = detail == "local"

    if roads is None:
        roads = style == "satellite" and detail == "local"

    if terrain is None:
        terrain = style == "weather"

    satellite = None

    if style == "satellite":
        TILE_CACHE_DIR.mkdir(parents=True, exist_ok=True)

        satellite = img_tiles.GoogleTiles(
            style="satellite",
            cache=TILE_CACHE_DIR,
        )

        if projection is None:
            projection = satellite.crs

    if projection is None:
        projection = _lambert_for_extent(extent)

    fig = plt.figure(
        figsize=figsize,
        dpi=dpi,
        facecolor=style_cfg["figure"],
    )

    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0],projection=projection,)

    # Always state the CRS of the supplied lon/lat extent.
    ax.set_extent(extent, crs=DATA_CRS)
    ax.set_box_aspect(0.6)
    ax.set_facecolor(style_cfg["figure"])

    if style == "satellite":
        if satellite_zoom is None:
            satellite_zoom = _auto_tile_zoom(extent)

        ax.add_image(
            satellite,
            satellite_zoom,
            alpha=satellite_alpha,
            zorder=0,
        )
    else:
        _add_flat_background(
            ax,
            style_cfg,
            detail_cfg["scale"],
        )

        if terrain:
            relief = _terrain_tiles()

            if terrain_zoom is None:
                terrain_zoom = _auto_terrain_zoom(extent)

            ax.add_image(
                relief,
                terrain_zoom,
                alpha=terrain_alpha,
                zorder=0.25,
            )

    _add_boundaries(
        ax,
        style_cfg,
        detail_cfg,
        counties=counties,
    )

    if roads:
        _add_roads(
            ax,
            detail,
            style_cfg,
            detail_cfg,
        )

    # Cleaner than forcing set_box_aspect() + tight_layout() on GeoAxes.
    try:
        ax.spines["geo"].set_visible(False)
    except (KeyError, AttributeError):
        pass

    return fig, ax


__all__ = [
    "build_map",
    "CONUS_EXTENT",
    "DATA_CRS",
    "MAP_STYLES",
]
