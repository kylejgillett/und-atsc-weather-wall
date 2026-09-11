


print("############\nSCRIPT RUNNING\n############")
import time as comp_time
st = comp_time.time()

import pandas as pd
import warnings
warnings.filterwarnings("ignore")

import geopandas
from datetime import datetime, timedelta, timezone
import requests
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import calendar
from metpy.plots import USCOUNTIES

import io
import sys
import os


def add_wwa_colors(gdf):
    """Colors for broad NWS watches/warnings/advisories."""

    if gdf.empty:
        return gdf

    gdf = gdf.copy()

    gdf["_plot_color"] = (
        gdf["event"]
        .map(NWS_ALERT_COLORS)
        .fillna(NWS_ALERT_DEFAULT_COLOR)
    )

    return gdf


def add_sbw_colors(gdf):
    """Radar-style colors for storm-based warnings."""

    if gdf.empty:
        return gdf

    gdf = gdf.copy()

    gdf["_plot_color"] = (
        gdf["event"]
        .map(SBW_COLORS)
        .fillna(NWS_ALERT_DEFAULT_COLOR)
    )

    return gdf



# get script dir
script_dir = os.path.dirname(os.path.abspath(__file__))

# get the parent dir 
project_root = os.path.abspath(os.path.join(script_dir, ".."))

if project_root not in sys.path:
    sys.path.append(project_root)

from utils.nws_alert_colors import NWS_ALERT_COLORS, NWS_ALERT_DEFAULT_COLOR, SBW_COLORS



from get_data.get_nws_headlines import get_nws_headlines
sbw, wwa = get_nws_headlines()

wwa = add_wwa_colors(wwa)
sbw = add_sbw_colors(sbw)



print(f"    FOUND {len(sbw)}, ACTIVE STORM BASED WARNINGS")
print(f"    FOUND {len(wwa)}, ACTIVE WATCHES/WARNINGS/ADVISORIES")


def build_map(extent=[-116, -76, 22, 56], projection=ccrs.LambertConformal(), style='black'):

    fig = plt.figure(figsize=(20, 10), dpi=250)
    fig.set_facecolor('#009946')
    ax = plt.axes(projection=projection)

    # apply the map extent (lat/lon bounding box)
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    # axis aspect ratio
    ax.set_box_aspect(0.6)
    # add map features
    
    if style == 'light':
        color = 'gray'
        alpha = 0.7
    else: 
        color = 'black'
        alpha = 0.8
    ax.add_feature(cfeature.STATES, edgecolor='black', alpha=0.5, linestyle='-', linewidth=1, zorder=10)
    ax.add_feature(cfeature.LAND, facecolor=color, alpha=alpha, zorder=1)
    ax.add_feature(cfeature.OCEAN, facecolor=color, alpha=alpha+0.2, zorder=0)
    #ax.add_feature(cfeature.BORDERS, color='white', alpha=1, linestyle='-', linewidth=1, zorder=11)
    ax.add_feature(cfeature.COASTLINE, color='black', alpha=0.5, linestyle='-', linewidth=1, zorder=11)
    # from cartopy.io import img_tiles
    # satellite = img_tiles.GoogleTiles(style='satellite')
    # ax.add_image(satellite, 4)

    # apply tight layout to the figure (keeps things tiddy)
    plt.tight_layout()

    # return the figure axis
    return fig, ax



fig, ax = build_map(style='dark')
ax.add_feature(USCOUNTIES.with_scale('20m'), alpha=0.1, edgecolor='black', linestyle='-', lw=0.5, zorder=12.1)


def prepare_alerts(gdf):

    if gdf.empty:
        return gdf

    gdf = gdf.copy()

    gdf["_plot_color"] = (
        gdf["event"]
        .map(NWS_ALERT_COLORS)
        .fillna(NWS_ALERT_DEFAULT_COLOR)
    )

    return gdf

wwa.plot(
    ax=ax,
    color=wwa["_plot_color"],
    edgecolor="black",
    linewidth=0.4,
    alpha=0.45,
    transform=ccrs.PlateCarree(),
    zorder=12.3,
)

sbw.plot(
    ax=ax,
    color=sbw["_plot_color"],
    edgecolor="black",
    linewidth=0.45,
    alpha=0.65,
    transform=ccrs.PlateCarree(),
    zorder=12.4,
)

#################################
# ADD MAP EXTRAS
#################################
# plot title, add one to the left with model name and data names, add another to the right with time info
plt.figtext(0.08, 1.03, f'   NWS HEADLINES', weight='bold', ha='left', fontsize=20, color='white')
plt.figtext(0.08, 1.00, f'   VALID: {datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M")}z', ha='left', 
            fontsize=18, color='white')
plt.figtext(0.915, 1.04, f' ', ha='left', fontsize=20)

from PIL import Image
img = Image.open('utils/images/und-logo.png')
#                  side-side  up-down  size   size
imgax = fig.add_axes([0.85, 1, 0.06, 0.06], anchor='SE', zorder=3)
imgax.imshow(img)
imgax.axis('off')


plt.savefig(f"staged_figures/nws_headlines/conus_nws_headlines.png", bbox_inches="tight")


elapsed_time = comp_time.time() - st
print(f"############\nSCRIPT FINISHED: time: {comp_time.strftime('%H:%M:%S', comp_time.gmtime(elapsed_time))}\n############")