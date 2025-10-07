print("############\nSCRIPT RUNNING\n############")

import time as comp_time
st = comp_time.time()

import warnings
warnings.filterwarnings("ignore")

import scipy.ndimage as ndimage
import metpy.calc as mpcalc
from metpy.units import units
import cartopy.crs   as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mcolors
from datetime import datetime, timezone
import numpy as np
import sys
import os
import xarray as xr


# get script dir
script_dir = os.path.dirname(os.path.abspath(__file__))

# get the parent dir 
project_root = os.path.abspath(os.path.join(script_dir, ".."))

if project_root not in sys.path:
    sys.path.append(project_root)

# import modules from sub dirs
from utils.colormaps import *
from get_data.get_metars import get_metar_data
from get_data.get_rap_data import analysis
from get_data.get_goes_from_aws import download_goes_file
from get_data.get_radar_mosaic import get_latest_mosaic
from get_data.get_wpc_bulletin import plot_bulletin




# parse date information for sat data download
utc_date = datetime.now(timezone.utc)
utc_doy = utc_date.timetuple().tm_yday
if utc_doy < 100:
    utc_doy = f'0{utc_doy}'
else:
    utc_doy = str(utc_doy)
utc_now = [utc_date.strftime("%Y"), utc_date.strftime("%m"), utc_date.strftime("%d"), utc_date.strftime("%H"), utc_doy]




# set up rap retrieval 
center_lat = 46.841203
center_lon = -98.777673
box_size   = 6.5
west = center_lon  - box_size
east = center_lon  + box_size
south = center_lat - box_size
north = center_lat + box_size

# pull rap data
raw_data = analysis(center_lat, center_lon, box_size=box_size)

# LATS & LONS
lats = raw_data.variables['lat'][:]
lons = raw_data.variables['lon'][:]

# PRES LEVS
pres_levs = raw_data['isobaric'][:]
pres_levs = pres_levs / 100

# DATE INFO
try:
    data_date = raw_data['time'].values[0]
except:
    data_date = raw_data['time1'].values[0]
    pass
valid_date = f'{data_date}'


# BASIC DATA EXTRACTION
ghgt_iso = ndimage.gaussian_filter(raw_data.variables['Geopotential_height_isobaric'][0], sigma=0.75)
temp_iso = ndimage.gaussian_filter(raw_data.variables['Temperature_isobaric'][0], sigma=0.75) - 273.15
uwnd_iso = ndimage.gaussian_filter(raw_data.variables['u-component_of_wind_isobaric'][0], sigma=0.75) * 1.94384
vwnd_iso = ndimage.gaussian_filter(raw_data.variables['v-component_of_wind_isobaric'][0], sigma=0.75) * 1.94384

pres_sfc = ndimage.gaussian_filter(raw_data.variables['MSLP_MAPS_System_Reduction_msl'][0], sigma=0.75)
temp_sfc = ndimage.gaussian_filter(raw_data.variables['Temperature_height_above_ground'][0], sigma=0.75) - 273.15
uwnd_sfc = ndimage.gaussian_filter(raw_data.variables['u-component_of_wind_height_above_ground'][0], sigma=0.75) * 1.94384
vwnd_sfc = ndimage.gaussian_filter(raw_data.variables['v-component_of_wind_height_above_ground'][0], sigma=0.75) * 1.94384
relh_sfc = ndimage.gaussian_filter(raw_data.variables['Relative_humidity_height_above_ground'][0], sigma=0.75)
dwpt_sfc = mpcalc.dewpoint_from_relative_humidity(temp_sfc*units.degC, relh_sfc*units.percent)





# build mpa function
def build_map(extent=[west, east, south, north], projection=ccrs.LambertConformal(), style='light'):
    fig = plt.figure(figsize=(20, 12))
    fig.set_facecolor('#009946')
    ax = plt.axes(projection=projection)

    ax.set_extent(extent)
    ax.set_box_aspect(0.7)

    if style == 'light':
        color = 'gray'
        alpha = 0.5
    else:
        color = 'black'
        alpha = 0.8

    ax.add_feature(cfeature.STATES, edgecolor='navy', alpha=0.5, linestyle='-', linewidth=1, zorder=10)
    ax.add_feature(cfeature.LAND, facecolor=color, alpha=alpha, zorder=1)
    ax.add_feature(cfeature.OCEAN, facecolor=color, alpha=alpha + 0.2, zorder=0)
    ax.add_feature(cfeature.COASTLINE, color='navy', alpha=0.5, linestyle='-', linewidth=1, zorder=11)

    plt.tight_layout()

    return fig, ax


# get radar mosaic data
radar_data, radar_lat, radar_lon, radar_time = get_latest_mosaic(utc_now[0], utc_now[1], utc_now[2])

# get metar data
metar_obs, metar_time = get_metar_data(reduced_to=50000)

# get satellite data
output_filename = download_goes_file(utc_now[0],utc_now[4],utc_now[3],None)

import time
def wait_for_download(file_path, timeout=240):
    start_time = time.time()
    while time.time() - start_time < timeout:
        if os.path.exists(file_path):
            # check if the file is still being written to
            try:
                # attempt to open the file in read mode
                with open(file_path, 'rb') as f:
                    pass  # if it opens without error, it's done downloading
                return True
            except IOError:
                pass
        time.sleep(5)
    return False

if wait_for_download(output_filename):
    print(f"File '{output_filename}' download complete.")
else:
    print(f"File '{output_filename}' download timed out.")


xrds = xr.open_dataset(output_filename)
radiance = xrds.variables['Rad'][:]

# Set time variables
scan_start = datetime.strptime(xrds.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ')
sat_time = scan_start.strftime('%H%M%S')

data = xrds.metpy.parse_cf('Rad')
data_crs = data.metpy.cartopy_crs
data_xcord = data.x
data_ycord = data.y

globe = ccrs.Globe(semimajor_axis=xrds['goes_imager_projection'].semi_major_axis,
                   semiminor_axis=xrds['goes_imager_projection'].semi_minor_axis)

crs = ccrs.Geostationary(central_longitude=xrds['goes_imager_projection'].longitude_of_projection_origin,
                         satellite_height=xrds['goes_imager_projection'].perspective_point_height,
                         globe=globe)
proj = ccrs.LambertConformal(globe=globe)

xrds.close()













fig, ax = build_map(style='dark')


###################################################################
# METAR STATION PLOTS
###################################################################
from metpy.plots import StationPlot, StationPlotLayout, sky_cover
custom_layout = StationPlotLayout()
custom_layout.add_barb('eastward_wind', 'northward_wind', units='knots')
custom_layout.add_value('NW', 'air_temperature', fmt='.0f', units='degF', color='darkred')
#custom_layout.add_value('NE', 'sea_level_pressure', fmt='.0f', units='mb', color='black')
custom_layout.add_value('SW', 'dew_point_temperature', fmt='.0f', units='degF', color='darkgreen')
custom_layout.add_symbol('C', 'cloud_coverage', sky_cover)
stationplot = StationPlot(ax, metar_obs['longitude'], metar_obs['latitude'], clip_on=True,
                          transform=ccrs.PlateCarree(), fontsize=8, zorder=12, alpha=1)
custom_layout.plot(stationplot, metar_obs)


###################################################################
# RAP MSLP
###################################################################
cs = ax.contour(lons, lats, pres_sfc/100, np.arange(904, 1054, 4), colors='black',
                linewidths=2.0, linestyles='-',
                transform=ccrs.PlateCarree(), zorder=11)
plt.clabel(cs, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
           rightside_up=True, use_clabeltext=True)


###################################################################
# RADAR MOSAIC
###################################################################
pm = ax.pcolormesh(radar_lon+0.05, radar_lat+0.05, radar_data,
              vmin=-15, vmax=95, cmap=rs_expertreflect_cmap, alpha=0.8, zorder=1.3, transform=ccrs.PlateCarree())

###################################################################
# LATEST FRONTS BULLETIN
###################################################################
texts, params, geoms, valid_time = plot_bulletin(ax)




###################################################################
# SATELLITE DATA
###################################################################
ax.imshow(radiance.values, origin='upper', cmap=ir_greys, vmin=50, vmax=130,
           extent=(data_xcord.values.min(), data_xcord.values.max(), data_ycord.values.min(), data_ycord.values.max()),
           regrid_shape=2000,
           aspect='auto',
           interpolation='gaussian',
           transform=crs,
           alpha=0.8, zorder=1)



#################################
# ADD MAP EXTRAS
#################################
# plot title, add one to the left with model name and data names, add another to the right with time info
plt.figtext(0.08, 1.03, f'   RAP Surface Analysis | {valid_date[0:10]} {valid_date[11:-13]}z', weight='bold', ha='left', fontsize=20, color='white')
plt.figtext(0.08, 1.00, f'    RAP MSLP (hPa), {metar_time[11:16]}z METARs, {valid_time}z WPC Fronts, {str(radar_time)[11:16]}z Reflectivity Mosaic, {sat_time[0:2]}:{sat_time[2:4]}z GOES16 Radiance', ha='left', fontsize=18, color='white')
plt.figtext(0.915, 1.04, f' ', ha='left', fontsize=20)
# colorbar for filled contour
cbar = plt.colorbar(pm, aspect=70, fraction=0.02, ax=ax, orientation='horizontal', pad=-0.01, extendrect=True)
cbar.set_label('Reflectivity (dBz)', fontsize=15, color='white')


# add UND logo
from PIL import Image
img = Image.open('utils/images/und-logo.png')
#                  side-side  up-down  size   size
imgax = fig.add_axes([0.83, 1.01, 0.06, 0.06], anchor='SE', zorder=3)
plt.figtext(0.81, 0.995, f'ATMOSPHERIC SCIENCES', ha='left', weight='bold', fontsize=10, color='white')
imgax.imshow(img)
imgax.axis('off')


plt.savefig("staged_figures/left_1_2.png", bbox_inches="tight")


elapsed_time = comp_time.time() - st
print(f"############\nSCRIPT FINISHED: time: {comp_time.strftime("%H:%M:%S", comp_time.gmtime(elapsed_time))}\n############")