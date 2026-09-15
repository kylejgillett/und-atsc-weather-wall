##########################################################
#       REGIONAL RAP AND SAT ANALYSIS PLOTS SCRIPT
#  (c) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2026
##########################################################

print("############\nSCRIPT RUNNING\n############")

import time as comp_time
st = comp_time.time()

import warnings
warnings.filterwarnings("ignore")

from matplotlib.colors import PowerNorm
import scipy.ndimage as ndimage
import metpy.calc as mpcalc
from metpy.plots import USCOUNTIES
from metpy.units import units
import cartopy.crs   as ccrs
import matplotlib.pyplot as plt
from datetime import datetime, timezone
import numpy as np
import sys
import os
import xarray as xr
import goes2go

# get script dir
script_dir = os.path.dirname(os.path.abspath(__file__))

# get the parent dir 
project_root = os.path.abspath(os.path.join(script_dir, ".."))

if project_root not in sys.path:
    sys.path.append(project_root)

# import modules from sub dirs
from utils.colormaps import *
from utils.utils import *
from get_data.get_metars import get_metar_data
from get_data.get_rap_data import analysis
from get_data.get_goes_from_aws import download_goes19_visible, subset_goes_to_map
from get_data.get_radar_mosaic import get_latest_mosaic
from get_data.get_wpc_bulletin import plot_bulletin
from utils.add_nws_headlines import add_nws_headlines
from utils.map import map_builder
from utils.figure import figure_builder


# parse date information for sat data download
utc_date = datetime.now(timezone.utc)
utc_doy = utc_date.timetuple().tm_yday
if utc_doy < 100:
    utc_doy = f'0{utc_doy}'
else:
    utc_doy = str(utc_doy)
utc_now = [utc_date.strftime("%Y"), utc_date.strftime("%m"), utc_date.strftime("%d"), utc_date.strftime("%H"), utc_doy]




# set up rap retrieval 
# GFK CENTERED
center_lat, center_lon = 46.841203, -98.777673

box_size   = 4.5 # 6.5
west = center_lon  - box_size
east = center_lon  + box_size
south = center_lat - box_size
north = center_lat + box_size

# pull rap data
raw_data = analysis(center_lat, center_lon, box_size=12)

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





# get radar mosaic data
radar_data, radar_lat, radar_lon, radar_time = get_latest_mosaic(utc_now[0], utc_now[1], utc_now[2])

# get metar data
try:
    metar_obs, metar_time = get_metar_data(reduced_to=50000)
    filtered_metars = metar_obs[
            (metar_obs['latitude'] >= center_lat - box_size*3) & (metar_obs['latitude'] <= center_lat + box_size*3) &
            (metar_obs['longitude'] >= center_lon - box_size*3) & (metar_obs['longitude'] <= center_lon + box_size*3)]
    filtered_metars['air_temperature'] = (filtered_metars['air_temperature']* 9/5) + 32
    filtered_metars['dew_point_temperature'] = (filtered_metars['dew_point_temperature']* 9/5) + 32
    bad_metar = (filtered_metars['air_temperature'] < -100) | (filtered_metars['dew_point_temperature'] < -100)
    if bad_metar.any():
        filtered_metars = filtered_metars[~bad_metar].reset_index(drop=True)
except:
    pass

# get satellite data
sat_file = download_goes19_visible(utc_now[0], utc_now[4], utc_now[3])
xrds_sat = xr.open_dataset(sat_file)
# sat_crs = xrds_sat.FOV.crs
# sat_x = xrds_sat.FOV.x.values
# sat_y = xrds_sat.FOV.y.values
# sat_extent = (float(np.nanmin(sat_x)), float(np.nanmax(sat_x)),
#               float(np.nanmin(sat_y)), float(np.nanmax(sat_y)))
sat_valid = datetime.fromisoformat(xrds_sat.time_coverage_start.replace("Z", "+00:00"))
sat_type = "GOES-19 Band 02 Visible"
sat_valid_str = sat_valid.strftime("%Y-%m-%d %H:%MZ")
sat_time_str = sat_valid.strftime("%H:%MZ")

visible, sat_crs, sat_extent = subset_goes_to_map(xrds_sat, west, east, south, north, pad_km=250)








# build map
fig, ax = map_builder(extent=[west, east, south, north], 
                      terrain=True, terrain_zoom=8, state_color='white', border_color='white')
ax.set_extent((west, east, south, north),crs=ccrs.PlateCarree())
fig.canvas.draw()

###################################################################
# SATELLITE
###################################################################
# ax.imshow(xrds_sat["CMI"].values, origin="upper", extent=sat_extent, transform=sat_crs,
#           cmap="gray", norm=PowerNorm(gamma=0.70, vmin=0.0, vmax=1.3), interpolation="nearest",
#           regrid_shape=900, alpha=0.90, zorder=1)
# xrds_sat.close()

ax.imshow(visible.values, origin="upper", extent=sat_extent, transform=sat_crs, cmap="gray",
          norm=PowerNorm(gamma=0.55, vmin=0.0, vmax=1.1), interpolation="nearest",
          regrid_shape=1200, alpha=0.85, zorder=1)
xrds_sat.close()


###################################################################
# METAR STATION PLOTS
###################################################################
try:
    from metpy.plots import StationPlot, StationPlotLayout, sky_cover
    from matplotlib.patheffects import withStroke
    TEXT_OUTLINE = [withStroke(linewidth=3.0, foreground="#1A1A18")]
    BARB_OUTLINE = [withStroke(linewidth=3.5, foreground="#1A1A18")]

    custom_layout = StationPlotLayout()
    custom_layout.add_barb('eastward_wind', 'northward_wind', units='knots', path_effects=BARB_OUTLINE)
    custom_layout.add_value('NW', 'air_temperature', fmt='.0f', fontsize=10, weight='bold', color="#FF662F", path_effects=TEXT_OUTLINE)
    custom_layout.add_value('SW', 'dew_point_temperature', fmt='.0f', fontsize=10, weight='bold', color="#8FE388", path_effects=TEXT_OUTLINE)
    custom_layout.add_symbol('C', 'cloud_coverage', sky_cover, path_effects=TEXT_OUTLINE)
    stationplot = StationPlot(ax, filtered_metars['longitude'], filtered_metars['latitude'], clip_on=True,
                            transform=ccrs.PlateCarree(), fontsize=10, zorder=12, alpha=1, color='white')
    custom_layout.plot(stationplot, filtered_metars)
except: 
    pass


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
              vmin=-32, vmax=95, cmap=rs_expertreflect_cmap, alpha=0.8, zorder=1.3, transform=ccrs.PlateCarree())


###################################################################
# LATEST FRONTS BULLETIN
###################################################################
texts, params, geoms, valid_time = plot_bulletin(ax)


# ##################################################################
# ADD NWS HEADLINES
# ##################################################################
# add_nws_headlines(ax, wwa_alpha=0.1, sbw_alpha=0.10, linewidth=0.5, zorder=16, legend=True)




#################################
# ADD MAP EXTRAS
#################################
# # plot title, add one to the left with model name and data names, add another to the right with time info
# plt.figtext(0.08, 1.03, f'   RAP Surface Analysis | {valid_date[0:10]} {valid_date[11:-13]}z', weight='bold', ha='left', fontsize=20, color='white')
# plt.figtext(0.08, 1.00, f'   RAP MSLP (hPa), {metar_time[11:16]}z METARs, {valid_time}z WPC Fronts, {str(radar_time)[11:16]}z Reflectivity Mosaic, {sat_time[0:2]}:{sat_time[2:4]}z GOES16 Radiance', ha='left', fontsize=18, color='white')
# plt.figtext(0.915, 1.04, f' ', ha='left', fontsize=20)
# # # colorbar for filled contour
# # cbar = plt.colorbar(pm, aspect=70, fraction=0.02, ax=ax, orientation='horizontal', pad=-0.01, extendrect=True)
# # cbar.set_label('Reflectivity (dBz)', fontsize=15, color='white')
# plt.figtext(0.915, 1.04, f' ', ha='left', fontsize=20)
# plt.figtext(0.915, -0.01, f' ', ha='left', fontsize=20)
# cax = fig.add_axes([0.91, 0.024, 0.01, 0.95])
# cbar = fig.colorbar(pm, cax=cax, orientation='vertical', ticks=np.arange(-30, 100, 5), extendrect=True)
# cax.text(3, 0.5, 'Reflectivity (dBz)', ha='left',va='center',rotation=270, color='white',fontsize=12,fontweight='bold',transform=cax.transAxes)
# cbar.ax.tick_params(axis='y', labelcolor='white') 
# for t in cbar.ax.get_yticklabels():
#     t.set_fontweight('bold')
#     t.set_fontsize(9)
# cbar.ax.set_facecolor('black')

# # add UND logo
# from PIL import Image
# img = Image.open('utils/images/und-logo.png')
# #                  side-side  up-down  size   size
# imgax = fig.add_axes([0.83, 1.01, 0.06, 0.06], anchor='SE', zorder=3)
# plt.figtext(0.81, 0.995, f'ATMOSPHERIC SCIENCES', ha='left', weight='bold', fontsize=10, color='white')
# imgax.imshow(img)
# imgax.axis('off')


composite_filename = build_filename("staged_figures/regional_surface_analysis/", f"regional_rap_analysis", utc_date)

figure_builder(fig, ax,
    title=f"RAP Analysis • Surface",
    subtitle=f'RAP MSLP (hPa), Surface Observations, WPC Fronts, MRMS Reflectivity Mosaic (dBz), GOES19 Satellite',
    valid=f"Valid • {valid_date[0:10]} {valid_date[11:-13]}z",
    mappable=pm,
    cbar_title="Reflectivity",
    cbar_units="dBz",
    cbar_ticks=np.arange(-30, 100, 5),
    footer_left=f"RAP 13km  •  INIT {valid_date[0:10]} {valid_date[11:-13]}z  •  {metar_time[11:16]}z METARs, {valid_time}z WPC Fronts, {str(radar_time)[11:16]}z Reflectivity Mosaic, {sat_time_str} {sat_type}",
    footer_right=' ',
    save_path=composite_filename)

elapsed_time = comp_time.time() - st
print(f"############\nSCRIPT FINISHED: time: {comp_time.strftime("%H:%M:%S", comp_time.gmtime(elapsed_time))}\n############")
