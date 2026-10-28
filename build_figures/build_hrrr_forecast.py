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
import glob
from netCDF4 import Dataset
import wrf

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


# WRF Data Directory
wrf_dir = "/data/jordan.i.christian/UND_WRF/WRFV4.6.1/run"
wrf_files = sorted(glob.glob(os.path.join(wrf_dir, "wrfout_d01_*")))

if not wrf_files:
    print("No WRF files found in directory.")
    sys.exit()

# Get initialization time from the first file
first_nc = Dataset(wrf_files[0])
init_time_str = first_nc.START_DATE
init_time = datetime.strptime(init_time_str, '%Y-%m-%d_%H:%M:%S').replace(tzinfo=timezone.utc)


for wrf_file in wrf_files:

    ncfile = Dataset(wrf_file)

    time_dt = wrf.extract_times(ncfile, timeidx=0, do_xtime=False)
    time_str = str(time_dt).split('.')[0]
    valid_time = datetime.strptime(time_str, '%Y-%m-%dT%H:%M:%S').replace(tzinfo=timezone.utc)

    fh = int((valid_time - init_time).total_seconds()/3600)

    temperature = wrf.getvar(ncfile, "T2")
    dewpoint = wrf.getvar(ncfile, "td2", units="K")
    uv10 = wrf.getvar(ncfile, "uvmet10", units="m s-1") 
    u = uv10[0, :]
    v = uv10[1, :]
    refl_3d = wrf.getvar(ncfile, "REFL_10CM")
    reflectivity = refl_3d.max(dim="bottom_top")
    # reflectivity = reflectivity.where(reflectivity >= 5)
    cloud_frac = wrf.getvar(ncfile, "cloudfrac")
    cloud_cover = np.amax(cloud_frac.values, axis=0)*100.0
    mslp = wrf.getvar(ncfile, "slp")

    lon = ncfile["XLONG"][0,:,:]
    lat = ncfile["XLAT"][0,:,:]

    # PLOT DATA
    fig, ax = map_builder(extent=[west, east, south, north], 
                         satellite=True, satellite_zoom=8, state_color='white', border_color='white')
    ax.set_extent((west, east, south, north),crs=ccrs.PlateCarree())
    fig.canvas.draw()

    zoom_factor = 2

    # ref_plot = zoom(reflectivity.values,zoom_factor,order=1)
    # lat_plot = zoom(lat,zoom_factor,order=1)
    # lon_plot = zoom(lon,zoom_factor,order=1)

    refl_masked = np.ma.masked_less(reflectivity.values, 5)
    
    pm = ax.pcolormesh(lon, lat, refl_masked,
                vmin=-32, vmax=95, cmap=rs_expertreflect_cmap, alpha=0.8, shading="auto",
                zorder=1.3, transform=ccrs.PlateCarree())


    cloud_smooth = ndimage.gaussian_filter(cloud_cover, sigma=0.8)
    cloud_masked = np.ma.masked_less(cloud_smooth, 20)

    cf = ax.contourf(lon, lat, cloud_masked, levels=[20, 35, 50, 65, 80, 90, 100], cmap=cloud_cmap, alpha=0.3, transform=ccrs.PlateCarree(), zorder=1.15)
        
    cs = ax.contour(lon, lat, mslp, np.arange(904, 1054, 2), colors='black',
                    linewidths=2.0, linestyles='-',
                    transform=ccrs.PlateCarree(), zorder=11)
    plt.clabel(cs, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
            rightside_up=True, use_clabeltext=True)



    TEXT_OUTLINE = [withStroke(linewidth=3.0, foreground="#1A1A18")]
    BARB_OUTLINE = [withStroke(linewidth=3.5, foreground="#1A1A18")]
    stride = 7
    station_lon = lon[::stride, ::stride].ravel()
    station_lat = lat[::stride, ::stride].ravel()
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


    composite_filename = build_filename("/data/jordan.i.christian/UND_WRF/und-atsc-weather-wall/staged_figures/hrrr_forecasts/", f"hrrr_forecast", now_utc, variant=f"{fh:02d}")
    figure_builder(fig, ax,
        title=f"{(valid_time - timedelta(hours=fh)).strftime('%HZ')} WRF Forecast  •  Surface",
        subtitle=f'Composite Simulated Relfectivity (dBz)  •  Total Cloud Cover (%)  •  MSLP (hPa)  •  Forecast Station Plot',
        valid=f"+ F{fh:03d}hr • VALID {valid_time.strftime('%a %d %b %Y').upper()} - {valid_time.strftime('%HZ')}",
        mappable=pm,
        cbar_title="Reflectivity",
        cbar_units="dBz",
        cbar_ticks=np.arange(-30, 100, 5),
        footer_left=f"WRF 12km  •  INIT {(valid_time - timedelta(hours=fh)).strftime('%d %b %Y %HZ').upper()}",
        footer_right=' ',
        save_path=composite_filename)

elapsed_time = comp_time.time() - st
print(f"############\nSCRIPT FINISHED: time: {comp_time.strftime('%H:%M:%S', comp_time.gmtime(elapsed_time))}\n############")
