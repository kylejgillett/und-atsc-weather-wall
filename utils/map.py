
##########################################################
#               BASIC MAP GENERATION SCRIPT
#  (c) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2026
##########################################################


from functools import lru_cache
from pathlib import Path
import sys
import os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature import ShapelyFeature
from cartopy.io import img_tiles
from cartopy.io.shapereader import Reader
from metpy.plots import USCOUNTIES


# get script dir
script_dir = os.path.dirname(os.path.abspath(__file__))
# get the parent dir 
project_root = os.path.abspath(os.path.join(script_dir, ".."))
if project_root not in sys.path:
    sys.path.append(project_root)
from utils.utils import *
from utils.style import *


# prep cached and saved gis data
UTILS_DIR = Path(__file__).resolve().parent
TERRAIN_CACHE = UTILS_DIR / ".cartopy_terrain_cache"
SATELLITE_CACHE = UTILS_DIR / ".cartopy_satellite_cache"
US_PRIMARY_ROADS = (
    UTILS_DIR
    / "tl_2024_us_primaryroads"
    / "tl_2024_us_primaryroads.shp")
ND_ROADS = (
    UTILS_DIR
    / "tl_2024_38_prisecroads"
    / "tl_2024_38_prisecroads.shp")
MN_ROADS = (
    UTILS_DIR
    / "tl_2024_27_prisecroads"
    / "tl_2024_27_prisecroads.shp")



def terrain_tiles():
    TERRAIN_CACHE.mkdir(parents=True, exist_ok=True)
    tiles = img_tiles.GoogleTiles(url="https://server.arcgisonline.com/ArcGIS/rest/services/World_Shaded_Relief/MapServer/tile/{z}/{y}/{x}.jpg", cache=TERRAIN_CACHE)
    return tiles


@lru_cache(maxsize=8)
def _read_geometries(path_string):
    path = Path(path_string)
    if not path.exists():
        return tuple()
    return tuple(Reader(str(path)).geometries())



def _add_roads(ax, road_type="primary", color=ROAD_COLOR, linewidth=0.45, alpha=0.30):
    if road_type == "primary":
        road_files = [US_PRIMARY_ROADS]

    elif road_type == "local":
        road_files = [ND_ROADS, MN_ROADS]

    for road_file in road_files:
        geometries = _read_geometries(str(road_file))
        if not geometries:
            continue

        roads = ShapelyFeature(geometries, DATA_CRS, facecolor="none", edgecolor=color, linewidth=linewidth, alpha=alpha,)
        ax.add_feature(roads,zorder=9.5)



DATA_CRS = ccrs.PlateCarree()
CONUS_EXTENT = (-122, -73, 21, 56)



def map_builder(
    extent=CONUS_EXTENT,
    projection=None,
    figsize=(15, 10),
    dpi=250,
    terrain=True,
    terrain_alpha=0.55,
    terrain_zoom=4,
    counties=False,
    county_alpha=0.28,
    county_width=0.40,
    roads=False,
    road_type="primary",
    road_color=ROAD_COLOR,
    road_width=0.45,
    road_alpha=0.30,
    map_scale="50m",
    county_scale="20m",
    satellite=False,
    satellite_zoom=10,
    satellite_alpha=0.8,
    state_color=STATE_COLOR,
    border_color=BORDER_COLOR,
    border_factor=1.0):



    # extent pieces
    west, east, south, north = extent

    # build fig
    fig = plt.figure(figsize=figsize, dpi=dpi, facecolor=FIGURE_BG)

    # set proj
    if projection is None:
        projection = ccrs.LambertConformal(
            central_longitude=(west + east) / 2,
            central_latitude=(south + north) / 2,
            standard_parallels=(30, 60))

    # define map axis area on fig
    MAP_RECT = [0.015, 0.065, 0.895, 0.805]

    # init map axis and add basic features
    ax = fig.add_axes(MAP_RECT, projection=projection,)
    ax.set_extent(extent, crs=DATA_CRS)
    ax.set_box_aspect(0.6)
    ax.set_facecolor(WATER_COLOR)
    ax.add_feature(cfeature.LAND.with_scale(map_scale), facecolor=LAND_COLOR, edgecolor="none", zorder=0,)
    ax.add_feature(cfeature.OCEAN.with_scale(map_scale), facecolor=WATER_COLOR, edgecolor="none", zorder=0)
    ax.add_feature(cfeature.LAKES.with_scale(map_scale),facecolor=WATER_COLOR, edgecolor="none", zorder=0.1)

    if terrain:
        ax.add_image( terrain_tiles(),terrain_zoom, alpha=terrain_alpha, zorder=0.25)

    if satellite:
        SATELLITE_CACHE.mkdir(parents=True, exist_ok=True)
        satellite_tiles = img_tiles.GoogleTiles(style="satellite", cache=SATELLITE_CACHE,)
        ax.add_image(satellite_tiles, satellite_zoom, alpha=satellite_alpha, zorder=0.25)

    ax.add_feature(cfeature.STATES.with_scale(map_scale), facecolor="none", edgecolor=state_color, linewidth=1.0*border_factor, alpha=0.80, zorder=10)
    ax.add_feature(cfeature.BORDERS.with_scale(map_scale), facecolor="none", edgecolor=border_color, linewidth=1.1*border_factor, alpha=0.88, zorder=10.1)
    ax.add_feature(cfeature.COASTLINE.with_scale(map_scale), facecolor="none",  edgecolor=border_color, linewidth=1.1*border_factor, alpha=0.88, zorder=10.2)

    if counties:
        ax.add_feature(USCOUNTIES.with_scale(county_scale), facecolor="none", edgecolor=COUNTY_COLOR, linewidth=county_width, alpha=county_alpha, zorder=9)

    if roads:
        _add_roads(ax, road_type=road_type, color=road_color,linewidth=road_width, alpha=road_alpha)

    geo_spine = ax.spines["geo"]

    geo_spine.set_visible(True)
    geo_spine.set_linewidth(2.0)
    geo_spine.set_edgecolor(UND_GREEN)
    geo_spine.set_zorder(20)

    return fig, ax