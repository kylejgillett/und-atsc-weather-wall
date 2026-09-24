##########################################################
#       REGIONAL RAP AND SAT ANALYSIS PLOTS SCRIPT
#  (c) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2026
##########################################################

print("############\nSCRIPT RUNNING\n############")

import time as comp_time
st = comp_time.time()

import warnings
warnings.filterwarnings("ignore")

from metpy.plots import USCOUNTIES
from metpy.units import units
import cartopy.crs   as ccrs
from datetime import datetime, timezone, timedelta
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
from scipy.ndimage import zoom
import numpy as np
import sys
import os

# get script dir
script_dir = os.path.dirname(os.path.abspath(__file__))

# get the parent dir 
project_root = os.path.abspath(os.path.join(script_dir, ".."))

if project_root not in sys.path:
    sys.path.append(project_root)

# import modules from sub dirs
from utils.colormaps import *
from utils.utils import *
from utils.map import map_builder
from utils.figure import figure_builder
from get_data.get_hrrr_data import hrrr_forecast


now_utc = datetime.now(timezone.utc)




import numpy as np
from metpy.units import units
from metpy.plots import StationPlot, StationPlotLayout
from matplotlib.patheffects import withStroke




from matplotlib.colors import LinearSegmentedColormap

cloud_cmap = LinearSegmentedColormap.from_list(
    "hrrr_clouds",
    [
        "#666B6D",  # 20%  - subtle
        "#7E8385",  # 35%
        "#999D9E",  # 50%
        "#B4B7B7",  # 65%
        "#CFD1D0",  # 80%
        "#E1E2DF",  # 90%
        "#EEEDE8",  # 100%
]
)


# GFK CENTERED
center_lat, center_lon = 46.841203, -98.777673
box_size   = 3 # 6.5
west = center_lon  - box_size
east = center_lon  + box_size
south = center_lat - box_size
north = center_lat + box_size


for fh, valid_time, raw_data in hrrr_forecast(center_lat=center_lat, center_lon=center_lon, box_size=8):

    temperature = raw_data["t2m"]
    dewpoint = raw_data["d2m"]
    u = raw_data["u10"]
    v = raw_data["v10"]
    reflectivity = raw_data["refc"]
    reflectivity = reflectivity.where(reflectivity >= 5)
    cloud_cover = raw_data["cloud_cover"]
    mslp = raw_data["mslp"]

    lon = raw_data["longitude"]
    lat = raw_data["latitude"]



    # PLOT DATA
    fig, ax = map_builder(extent=[west, east, south, north], 
                         satellite=True, satellite_zoom=8, state_color='white', border_color='white')
    ax.set_extent((west, east, south, north),crs=ccrs.PlateCarree())
    fig.canvas.draw()

    zoom_factor = 2

    ref_plot = zoom(reflectivity.values,zoom_factor,order=1)
    lat_plot = zoom(lat.values, zoom_factor,order=1)
    lon_plot = zoom(lon.values,zoom_factor,order=1)

    pm = ax.pcolormesh(lon_plot, lat_plot, ref_plot,
                vmin=-32, vmax=95, cmap=rs_expertreflect_cmap, alpha=0.8, shading="auto",
                zorder=1.3, transform=ccrs.PlateCarree())


    cloud_smooth = ndimage.gaussian_filter(cloud_cover.values, sigma=0.8)
    cloud_masked = np.ma.masked_less(cloud_smooth, 20)

    cf = ax.contourf(lon, lat, cloud_masked, levels=[20, 35, 50, 65, 80, 90, 100], cmap=cloud_cmap, alpha=0.3, transform=ccrs.PlateCarree(), zorder=1.15)
        
    cs = ax.contour(lon, lat, mslp/100, np.arange(904, 1054, 4), colors='black',
                    linewidths=2.0, linestyles='-',
                    transform=ccrs.PlateCarree(), zorder=11)
    plt.clabel(cs, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
            rightside_up=True, use_clabeltext=True)



    TEXT_OUTLINE = [withStroke(linewidth=3.0, foreground="#1A1A18")]
    BARB_OUTLINE = [withStroke(linewidth=3.5, foreground="#1A1A18")]
    stride = 28
    station_lon = lon.values[::stride, ::stride].ravel()
    station_lat = lat.values[::stride, ::stride].ravel()
    station_temp = (temperature.values[::stride, ::stride].ravel()* units.kelvin).to("degF")
    station_dewpoint = (dewpoint.values[::stride, ::stride].ravel()* units.kelvin).to("degF")
    station_u = (u.values[::stride, ::stride].ravel()* units("m/s")).to("knots")
    station_v = (v.values[::stride, ::stride].ravel()* units("m/s")).to("knots")
    station_data = {
        "air_temperature": station_temp,
        "dew_point_temperature": station_dewpoint,
        "eastward_wind": station_u,
        "northward_wind": station_v}
    custom_layout = StationPlotLayout()
    custom_layout.add_barb("eastward_wind", "northward_wind", units="knots", path_effects=BARB_OUTLINE,)
    custom_layout.add_value("NW", "air_temperature", fmt=".0f", fontsize=8, weight="bold", color="#FF662F", path_effects=TEXT_OUTLINE,)
    custom_layout.add_value( "SW","dew_point_temperature", fmt=".0f", fontsize=8, weight="bold", color="#8FE388", path_effects=TEXT_OUTLINE)
    stationplot = StationPlot(ax, station_lon, station_lat, clip_on=True, transform=ccrs.PlateCarree(), fontsize=10, zorder=12, alpha=1, color="white", spacing=7)
    custom_layout.plot(stationplot, station_data)


    composite_filename = build_filename("staged_figures/hrrr_forecasts/", f"hrrr_forecast", now_utc, variant=f"{fh:02d}")
    figure_builder(fig, ax,
        title=f"{(valid_time - timedelta(hours=fh)).strftime("%HZ")} HRRR Forecast  •  Surface",
        subtitle=f'Composite Simulated Relfectivity (dBz)  •  Total Cloud Cover (%)  •  MSLP (hPa)  •  Forecast Station Plot',
        valid=f"+ F{fh:03d}hr • VALID {valid_time.strftime('%a %d %b %Y').upper()} - {valid_time.strftime('%HZ')}",
        mappable=pm,
        cbar_title="Reflectivity",
        cbar_units="dBz",
        cbar_ticks=np.arange(-30, 100, 5),
        footer_left=f"HRRR 3km  •  INIT {(valid_time - timedelta(hours=fh)).strftime('%d %b %Y %HZ').upper()}",
        footer_right=' ',
        save_path=composite_filename)

elapsed_time = comp_time.time() - st
print(f"############\nSCRIPT FINISHED: time: {comp_time.strftime("%H:%M:%S", comp_time.gmtime(elapsed_time))}\n############")