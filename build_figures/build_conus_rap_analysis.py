##########################################################
#              RAP ANALYSIS FIGURES SCRIPT
#  (c) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2025
##########################################################

print("############\nSCRIPT RUNNING\n############")
import time as comp_time
st = comp_time.time()

import warnings
warnings.filterwarnings("ignore")

import scipy.ndimage as ndimage
import metpy.calc as mpcalc
from metpy.units import units
from metpy.interpolate import interpolate_to_isosurface
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from datetime import datetime, timezone
from matplotlib.colors import Normalize
from matplotlib.colors import PowerNorm
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
from utils.utils import *
from utils.colormaps import *
from get_data.get_metars import get_metar_data
from get_data.get_rap_data import analysis
from get_data.get_goes_from_aws import download_goes19_visible
from get_data.get_radar_mosaic import get_latest_mosaic
from get_data.get_wpc_bulletin import plot_bulletin
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


def c2f(celsius_array):
    return (np.asanyarray(celsius_array) * 1.8) + 32


now_utc = datetime.now(timezone.utc)

#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
# set up rap retrieval 
center_lat = 46.841203
center_lon = -98.777673
box_size   = 50
west = center_lon  - box_size
east = center_lon  + box_size
south = center_lat - box_size
north = center_lat + box_size

projection = ccrs.LambertConformal(
    central_longitude=-95.0,
    central_latitude=25.0,
    standard_parallels=(25.0, 25.0)
)

# pull rap data
raw_data = analysis(box_size=box_size)

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
    try:
        data_date = raw_data[f'time1'].values[0]
    except:
        try:
            data_date = raw_data[f'time2'].values[0]
        except:
            pass
        pass
    pass

valid_date = f'{data_date}'


# BASIC DATA EXTRACTION
sigma = 1.5
ghgt_iso = ndimage.gaussian_filter(raw_data.variables['Geopotential_height_isobaric'][0], sigma=sigma)
temp_iso = ndimage.gaussian_filter(raw_data.variables['Temperature_isobaric'][0], sigma=sigma) - 273.15
uwnd_iso = ndimage.gaussian_filter(raw_data.variables['u-component_of_wind_isobaric'][0], sigma=sigma) * 1.94384
vwnd_iso = ndimage.gaussian_filter(raw_data.variables['v-component_of_wind_isobaric'][0], sigma=sigma) * 1.94384

pres_sfc = ndimage.gaussian_filter(raw_data.variables['MSLP_MAPS_System_Reduction_msl'][0], sigma=sigma)
temp_sfc = ndimage.gaussian_filter(raw_data.variables['Temperature_height_above_ground'][0], sigma=sigma) - 273.15
uwnd_sfc = ndimage.gaussian_filter(raw_data.variables['u-component_of_wind_height_above_ground'][0], sigma=sigma) * 1.94384
vwnd_sfc = ndimage.gaussian_filter(raw_data.variables['v-component_of_wind_height_above_ground'][0], sigma=sigma) * 1.94384
relh_sfc = ndimage.gaussian_filter(raw_data.variables['Relative_humidity_height_above_ground'][0], sigma=sigma)
dwpt_sfc = mpcalc.dewpoint_from_relative_humidity(temp_sfc*units.degC, relh_sfc*units.percent)
cape_ml = ndimage.gaussian_filter(raw_data.variables['Convective_available_potential_energy_pressure_difference_layer'][0,0], sigma=sigma)
cin_ml  = ndimage.gaussian_filter(raw_data.variables['Convective_inhibition_pressure_difference_layer'][0,0], sigma=sigma)

# --- BUNKERS STORM MOTION --------------------------------------------------
bunkers_sm_u = raw_data['U-Component_Storm_Motion_height_above_ground_layer'][0,0]
bunkers_sm_v = raw_data['V-Component_Storm_Motion_height_above_ground_layer'][0,0]

# --- STORM-RELATIVE MEAN WIND (6-9km and 0-2km) -----------------------------
hgt_iso = raw_data['Geopotential_height_isobaric'][0]
mask_6_9km = (hgt_iso >= 6000) & (hgt_iso <= 9000)
u_6_9km = raw_data['u-component_of_wind_isobaric'][0].where(mask_6_9km).mean(dim='isobaric', skipna=True)
v_6_9km = raw_data['v-component_of_wind_isobaric'][0].where(mask_6_9km).mean(dim='isobaric', skipna=True)
u_0_2km = raw_data['u-component_of_wind_height_above_ground'][0].sel(height_above_ground=slice(0,2000)).mean(dim='height_above_ground', skipna=True)
v_0_2km = raw_data['v-component_of_wind_height_above_ground'][0].sel(height_above_ground=slice(0,2000)).mean(dim='height_above_ground', skipna=True)
srw_6_9_u = (u_6_9km - bunkers_sm_u) * 1.94384
srw_6_9_v = (v_6_9km - bunkers_sm_v) * 1.94384
srw_0_2_u = (u_0_2km - bunkers_sm_u) * 1.94384
srw_0_2_v = (v_0_2km - bunkers_sm_v) * 1.94384




#################################
# CALCULATE FRONTOGENESIS
################################# 
# using MetPy, compute frontogensis
fgen = mpcalc.frontogenesis(raw_data['Temperature_isobaric'], 
                            raw_data['u-component_of_wind_isobaric'],
                            raw_data['v-component_of_wind_isobaric'],
                            latitude=raw_data['lon'].values,longitude=raw_data['lon'].values,crs=ccrs.PlateCarree())
# convert units to delta deg C / 100 km / 3 hr
fgen = fgen.metpy.convert_units('delta_degC/km/hour')*3*100
# create `fgen_masked`, a DataArray of fgen values >2, all else are Nan
fgen_masked = fgen.where(fgen > (2 * units('delta_degC/km/hour')))


#################################
# CALCULATE ADVECTION
################################# 
# using MetPy, compute temperature advection
adv = mpcalc.advection(raw_data['Temperature_isobaric'], 
                       raw_data['u-component_of_wind_isobaric'],
                       raw_data['v-component_of_wind_isobaric'],
                       latitude=raw_data['lat'].values,longitude=raw_data['lon'].values,crs=ccrs.PlateCarree())
# convert units to delta deg C / hr
adv = adv.metpy.convert_units('delta_degC/hour')
# apply some smoothing to adv 
adv = ndimage.gaussian_filter(adv, sigma=2, order=0) * units('K/sec')


# # use metpy to compute theta & add it into the `rap-data` DataSet
# raw_data['theta'] = mpcalc.potential_temperature(raw_data['isobaric'],raw_data['Temperature_isobaric'])

# # use metpy to compute PV and add it into the `rap_data` DataSet
# # Let MetPy compute spacing internally from lat/lon coordinates
# raw_data['pv'] = mpcalc.potential_vorticity_baroclinic(raw_data['theta'],
#                                                         raw_data['isobaric'],
#                                                         u=raw_data['u-component_of_wind_isobaric'],
#                                                         v=raw_data['v-component_of_wind_isobaric'],
#                                                         latitude=raw_data['lat'],
#                                                         longitude=raw_data['lon'])

# thta_on_2pvu = interpolate_to_isosurface(raw_data['pv'].values, raw_data['theta'].values,  2*1e-6, bottom_up_search=True)
# u_on_2pvu    = interpolate_to_isosurface(raw_data['pv'][:,0,:,:].values, raw_data['u-component_of_wind_isobaric'][0,:,:,:].values, 2*1e-6, bottom_up_search=True)
# v_on_2pvu    = interpolate_to_isosurface(raw_data['pv'][:,0,:,:].values, raw_data['v-component_of_wind_isobaric'][0,:,:,:].values, 2*1e-6, bottom_up_search=True)

#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################








#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
# get radar mosaic data
radar_data, radar_lat, radar_lon, radar_time = get_latest_mosaic(utc_now[0], utc_now[1], utc_now[2])


try:
    # get metar data
    metar_obs, metar_time = get_metar_data(reduced_to=150000)
    metar_obs['air_temperature'] = (metar_obs['air_temperature']* 9/5) + 32
    metar_obs['dew_point_temperature'] = (metar_obs['dew_point_temperature']* 9/5) + 32
    bad_metar = (metar_obs['air_temperature'] < -100) | (metar_obs['dew_point_temperature'] < -100)
    if bad_metar.any():
        metar_obs = metar_obs[~bad_metar].reset_index(drop=True)
except:
    pass


# get satellite data
sat_file = download_goes19_visible(utc_now[0], utc_now[4], utc_now[3])
xrds_sat = xr.open_dataset(sat_file)
sat_crs = xrds_sat.FOV.crs
sat_x = xrds_sat.FOV.x.values
sat_y = xrds_sat.FOV.y.values
sat_extent = (float(np.nanmin(sat_x)), float(np.nanmax(sat_x)),
              float(np.nanmin(sat_y)), float(np.nanmax(sat_y)))
sat_valid = datetime.fromisoformat(xrds_sat.time_coverage_start.replace("Z", "+00:00"))
sat_type = "GOES-19 Band 02 Visible"
sat_valid_str = sat_valid.strftime("%Y-%m-%d %H:%MZ")
sat_time_str = sat_valid.strftime("%H:%MZ")


#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################






#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
#################################
# BUILD 300 HPA MAP
#################################
fig, ax = map_builder(projection=projection, extent=[-119, -74, 23.5, 53.5], terrain=True, terrain_zoom=6, state_color='navy', border_color='navy')

# slice data
plev300 = np.where(pres_levs == 300)[0][0]
ghgt_300 = ghgt_iso[plev300]
uwnd_300 = uwnd_iso[plev300]
vwnd_300 = vwnd_iso[plev300]
wdsp_300 = np.sqrt(uwnd_300**2 + vwnd_300**2)

# plot 300 hpa heights
contour = ax.contour(lons, lats, ghgt_300, np.arange(0, 12000, 120),
                colors='black', linewidths=3.0, linestyles='-',
                transform=ccrs.PlateCarree(), zorder=11)
plt.clabel(contour, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
           rightside_up=True, use_clabeltext=True)

# plot 300 hpa wind speed
contourf = ax.contourf(lons, lats, wdsp_300, np.arange(50, 160, 5), extend='max',
                 cmap=wdsp_cmap, alpha=0.7, transform=ccrs.PlateCarree(), zorder=4)

# plot 300 hpa wind barbs
every = 20
barbs = ax.barbs(lons.values[0::every, 0::every], lats.values[0::every, 0::every],
                uwnd_300[0::every, 0::every], vwnd_300[0::every, 0::every],
                length=6.5, alpha=0.7, transform=ccrs.PlateCarree(), zorder=12)


# # plot title, add one to the left with model name and data names, add another to the right with time info
# plt.figtext(0.08, 1.03, f'     RAP 300 hPa Analysis | {valid_date[0:10]} {valid_date[11:-13]}z', weight='bold', ha='left', fontsize=20, color='white')
# plt.figtext(0.08, 1.00, f'     Heights (m), Wind (kt)', ha='left', fontsize=18, color='white')
# plt.figtext(0.915, 1.04, f' ', ha='left', fontsize=20)
# plt.figtext(0.915, -0.01, f' ', ha='left', fontsize=20)
# cax = fig.add_axes([0.91, 0.024, 0.01, 0.95])
# cbar = fig.colorbar(contourf, cax=cax, orientation='vertical', ticks=np.arange(50, 160, 5)[::1], extendrect=True)
# cax.text(3, 0.5, f'Wind Speed (kts)',ha='left',va='center',rotation=270, color='white',fontsize=12,fontweight='bold',transform=cax.transAxes)
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


composite_filename = build_filename("staged_figures/conus_rap_analysis/", f"conus_analysis", now_utc, variant='300a')

figure_builder(fig, ax,
    title=f"RAP Analysis • 300hPa",
    subtitle=f"Heights (m), Wind (kt)",
    valid=f"Valid • {valid_date[0:10]} {valid_date[11:-13]}z",
    mappable=contourf,
    cbar_title="Wind Speed",
    cbar_units="kts",
    cbar_ticks=np.arange(50, 160, 5)[::1],
    footer_left=f"RAP 13km • INIT {valid_date[0:10]} {valid_date[11:-13]}z",
    save_path=composite_filename)

print("    FINISHED 300HPA FLOW MAP")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################










# #############################################################################################################################################################################
# #############################################################################################################################################################################
# #############################################################################################################################################################################
# #################################
# # BUILD 300HPA PVA MAP
# #################################
# fig, ax = build_map()

# # use 300hpa slices from above

# # plot 300hpa heights
# contour = ax.contour(lons, lats, ghgt_300, np.arange(0, 12000, 60),
#                 colors='black', linewidths=3.0, linestyles='-',
#                 transform=ccrs.PlateCarree(), zorder=11)
# plt.clabel(contour, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
#            rightside_up=True, use_clabeltext=True)


# # plot 300hpa pv fill
# contourf = ax.contourf(raw_data['lon'], raw_data['lat'], raw_data['pv'][plev300,0,:,:]*1e6, pv_clevs, cmap=pv_cmap,
#                  transform=ccrs.PlateCarree(),extend='both')

# # plot a single dashed contour @ 2PVU
# pv_contour = ax.contour(raw_data['lon'], raw_data['lat'], raw_data['pv'][plev300,0,:,:]*1e6, [2], colors='navy',linestyles='dashed',linewidths=2,
#                  transform=ccrs.PlateCarree())

# # plot 300 hpa wind barbs
# every = 15
# barbs = ax.barbs(lons.values[0::every, 0::every], lats.values[0::every, 0::every],
#                 uwnd_300[0::every, 0::every], vwnd_300[0::every, 0::every],
#                 length=6.5, alpha=0.7, transform=ccrs.PlateCarree(), zorder=11)

# # plot title, add one to the left with model name and data names, add another to the right with time info
# plt.figtext(0.08, 1.03, f'     RAP 300 hPa Analysis | {valid_date[0:10]} {valid_date[11:-13]}z', weight='bold', ha='left', fontsize=20, color='white')
# plt.figtext(0.08, 1.00, f'     Heights (m), Potential Vorticity (PVU), Wind (kt)', ha='left', fontsize=18, color='white')
# plt.figtext(0.915, 1.04, f' ', ha='left', fontsize=20)
# plt.figtext(0.915, -0.01, f' ', ha='left', fontsize=20)
# cax = fig.add_axes([0.91, 0.024, 0.01, 0.95])
# cbar = fig.colorbar(contourf, cax=cax, orientation='vertical', ticks=pv_clevs[::5], extendrect=True)
# cax.text(3, 0.5, r'Potential Vorticity Units (PVU; $\rm{10^{-6}\ K\ kg^{-1}\ m^{2}\ s^{-1}})$' + ' | 2PVU (dashed)', ha='left',va='center',rotation=270, color='white',fontsize=12,fontweight='bold',transform=cax.transAxes)
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

#composite_filename = build_filename("staged_figures/conus_rap_analysis/", f"conus_analysis", data_date.astype('datetime64[us]').item(), variant='300b')
#plt.savefig(composite_filename, bbox_inches="tight")

# print("    FINISHED 300HPA PVA MAP")
# #############################################################################################################################################################################
# #############################################################################################################################################################################
# #############################################################################################################################################################################








#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
#################################
# BUILD 500HPA FLOW MAP
#################################
fig, ax = map_builder(projection=projection, extent=[-119, -74, 23.5, 53.5], terrain=True, state_color='navy', border_color='navy', border_factor=1.5)

# slice data 
plev500 = np.where(pres_levs == 500)[0][0]
ghgt_500 = ghgt_iso[plev500]
uwnd_500 = uwnd_iso[plev500]
vwnd_500 = vwnd_iso[plev500]
wdsp_500 = np.sqrt(uwnd_500**2 + vwnd_500**2)


# plot 500 hpa heights
contour = ax.contour(lons, lats, ghgt_500, np.arange(3000, 7000, 60),
                colors='black', linewidths=3.0, linestyles='-',
                transform=ccrs.PlateCarree(), zorder=11)
plt.clabel(contour, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
           rightside_up=True, use_clabeltext=True)

# plot 500hpa wind speed
contourf = ax.contourf(lons, lats, wdsp_500, np.arange(30, 140, 5), extend='max',
                 cmap=wdsp_cmap, alpha=0.7, transform=ccrs.PlateCarree(), zorder=4)

# plot 500 hpa wind barbs
every = 20
barbs = ax.barbs(lons.values[0::every, 0::every], lats.values[0::every, 0::every],
                uwnd_500[0::every, 0::every], vwnd_500[0::every, 0::every],
                length=6.5, alpha=0.7, transform=ccrs.PlateCarree(), zorder=12)

# # plot title, add one to the left with model name and data names, add another to the right with time info
# plt.figtext(0.08, 1.03, f'     RAP 500 hPa Analysis | {valid_date[0:10]} {valid_date[11:-13]}z', weight='bold', ha='left', fontsize=20, color='white')
# plt.figtext(0.08, 1.00, f'     Heights (m), Wind (kt)', ha='left', fontsize=18, color='white')
# plt.figtext(0.915, 1.04, f' ', ha='left', fontsize=20)
# plt.figtext(0.915, -0.01, f' ', ha='left', fontsize=20)
# cax = fig.add_axes([0.91, 0.024, 0.01, 0.95])
# cbar = fig.colorbar(contourf, cax=cax, orientation='vertical', ticks=np.arange(30, 140, 5), extendrect=True)
# cax.text(3, 0.5, 'Wind Speed (kts)', ha='left',va='center',rotation=270, color='white',fontsize=12,fontweight='bold',transform=cax.transAxes)
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


composite_filename = build_filename("staged_figures/conus_rap_analysis/", f"conus_analysis", now_utc, variant='500a')

figure_builder(fig, ax,
    title=f"RAP Analysis • 500hPa",
    subtitle=f"Heights (m), Wind (kt)",
    valid=f"Valid • {valid_date[0:10]} {valid_date[11:-13]}z",
    mappable=contourf,
    cbar_title="Wind Speed",
    cbar_units="kts",
    cbar_ticks=np.arange(30, 140, 5),
    footer_left=f"RAP 13km • INIT {valid_date[0:10]} {valid_date[11:-13]}z",
    save_path=composite_filename)

print("    FINISHED 500HPA FLOW MAP")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################











#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
#################################
# BUILD 500 HPA REL VORT MAP
#################################
fig, ax = map_builder(projection=projection, extent=[-119, -74, 23.5, 53.5], terrain=True, state_color='navy', border_color='navy', border_factor=1.5)

n_reps = 150

# compute vorticity and vorticity advection
dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)
f = mpcalc.coriolis_parameter(np.deg2rad(lats)).to('1/s')
vor_500 = mpcalc.smooth_n_point(mpcalc.vorticity(uwnd_500*units.kts, vwnd_500*units.kts, dx=dx, dy=dy), 9, n_reps)

# plot 500 hpa heights
contour = ax.contour(lons, lats, ghgt_500, np.arange(3000, 7000, 60),
                colors='black', linewidths=3.0, linestyles='-',
                transform=ccrs.PlateCarree(), zorder=11)
plt.clabel(contour, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
           rightside_up=True, use_clabeltext=True)

# plot relative vorticity fill
norm = mcolors.TwoSlopeNorm(vmin=-30, vcenter=0, vmax=50)
contourf = ax.contourf(lons, lats, vor_500 * 10**5, np.arange(-30, 52, 2), 
                      norm=norm, extend='both', cmap='PuOr_r', zorder=5, alpha=1, transform=ccrs.PlateCarree())

# plot wind barbs
every = 20
barbs = ax.barbs(lons.values[0::every, 0::every], lats.values[0::every, 0::every],
                uwnd_500[0::every, 0::every], vwnd_500[0::every, 0::every],
                length=6.5, alpha=0.7, transform=ccrs.PlateCarree(), zorder=12)

# # plot title, add one to the left with model name and data names, add another to the right with time info
# plt.figtext(0.08, 1.03, f'     RAP 500 hPa Analysis | {valid_date[0:10]} {valid_date[11:-13]}z', weight='bold', ha='left', fontsize=20, color='white')
# plt.figtext(0.08, 1.00, f'     Heights (m), Rel. Vorticity (/sec•10⁵), Wind (kt)', ha='left', fontsize=18, color='white')
# plt.figtext(0.915, 1.04, f' ', ha='left', fontsize=20)
# plt.figtext(0.915, -0.01, f' ', ha='left', fontsize=20)
# cax = fig.add_axes([0.91, 0.024, 0.01, 0.95])
# cbar = fig.colorbar(vort_cf, cax=cax, orientation='vertical', ticks=np.arange(-30, 52, 2), extendrect=True)
# cax.text(3, 0.5, 'Relative Vorticity (/sec•10⁵)', ha='left',va='center',rotation=270, color='white',fontsize=12,fontweight='bold',transform=cax.transAxes)
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

# composite_filename = build_filename("staged_figures/conus_rap_analysis/", f"conus_analysis", data_date.astype('datetime64[us]').item(), variant='500b')
# plt.savefig(composite_filename, bbox_inches="tight")

composite_filename = build_filename("staged_figures/conus_rap_analysis/", f"conus_analysis", now_utc, variant='500b')

figure_builder(fig, ax,
    title=f"RAP Analysis • 500hPa",
    subtitle=f"Heights (m), Rel. Vorticity "+r"($\mathrm{s}^{-1} \times 10^{5}$)" +", Wind (kt)",
    valid=f"Valid • {valid_date[0:10]} {valid_date[11:-13]}z",
    mappable=contourf,
    cbar_title="Relative Vorticity",
    cbar_units=r"$\mathrm{s}^{-1} \times 10^{5}$",
    cbar_ticks=np.arange(-30, 52, 2),
    footer_left=f"RAP 13km • INIT {valid_date[0:10]} {valid_date[11:-13]}z",
    save_path=composite_filename)

print("    FINISHED 500HPA REL VORT MAP")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################








# #############################################################################################################################################################################
# #############################################################################################################################################################################
# #############################################################################################################################################################################
# #################################
# # BUILD 500HPA ABSVORTADV MAP
# #################################
# fig, ax = build_map()

# # use 500hpa data and vort calculations from above 
# relvort_adv = mpcalc.advection(vor_500, uwnd_500, vwnd_500, dx=dx, dy=dy) *1e9

# # plot 500hpa heights
# contour = ax.contour(lons, lats, ghgt_500, np.arange(3000, 7000, 60),
#                 colors='black', linewidths=3.0, linestyles='-',
#                 transform=ccrs.PlateCarree(), zorder=11)
# plt.clabel(contour, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
#            rightside_up=True, use_clabeltext=True)

# # plot 500hpa abs vort adv fill
# vortadv_cf = ax.contourf(lons, lats, ndimage.gaussian_filter(relvort_adv, 4), np.arange(-50, 52, 2),               #np.arange(-6*12**-7, 6*12**-7, 1*10**-9),
#                              extend='both', cmap='bwr', zorder=5, alpha=1, transform=ccrs.PlateCarree())

# # plot 500hpa wind barbs
# every = 20
# barbs = ax.barbs(lons.values[0::every, 0::every], lats.values[0::every, 0::every],
#                 uwnd_500[0::every, 0::every], vwnd_500[0::every, 0::every],
#                 length=6.5, alpha=0.7, transform=ccrs.PlateCarree(), zorder=12)

# # plot title, add one to the left with model name and data names, add another to the right with time info
# plt.figtext(0.08, 1.03, f'     RAP 500 hPa Analysis | {valid_date[0:10]} {valid_date[11:-13]}z', weight='bold', ha='left', fontsize=20, color='white')
# plt.figtext(0.08, 1.00, f'     Heights (m), Rel. Vorticity Adv. (sec⁻²•10⁹), Wind (kt)', ha='left', fontsize=18, color='white')
# plt.figtext(0.915, 1.04, f' ', ha='left', fontsize=20)
# plt.figtext(0.915, -0.01, f' ', ha='left', fontsize=20)
# cax = fig.add_axes([0.91, 0.024, 0.01, 0.95])
# cbar = fig.colorbar(vortadv_cf, cax=cax, orientation='vertical', ticks=np.arange(-50, 52, 2), extendrect=True)
# cax.text(3, 0.5, 'Relative Vorticity Advection (sec⁻²•10⁹)', ha='left',va='center',rotation=270, color='white',fontsize=12,fontweight='bold',transform=cax.transAxes)
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

# composite_filename = build_filename("staged_figures/conus_rap_analysis/", f"conus_analysis", data_date.astype('datetime64[us]').item(), variant='500c')
# plt.savefig(composite_filename, bbox_inches="tight")

# print("    FINISHED 500HPA REL VORT ADV MAP")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################








#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
#################################
# BUILD 700HPA TEMP MAP
#################################
fig, ax = map_builder(projection=projection, extent=[-119, -74, 23.5, 53.5], terrain=True, state_color='navy', border_color='navy', border_factor=1.5)

# slice data
plev700 = np.where(pres_levs == 700)[0][0]
ghgt_700 = ghgt_iso[plev700]
uwnd_700 = uwnd_iso[plev700]
vwnd_700 = vwnd_iso[plev700]
temp_700 = temp_iso[plev700]
wdsp_700 = np.sqrt(uwnd_700**2 + vwnd_700**2)


# plot 700 hpa heights
contour = ax.contour(lons, lats, ghgt_700, np.arange(1800, 4000, 30),
                colors='black', linewidths=3.0, linestyles='-',
                transform=ccrs.PlateCarree(), zorder=11)
plt.clabel(contour, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
           rightside_up=True, use_clabeltext=True)

# plot 0C isotherm
ax.contour(lons, lats, temp_700, levels=[0], linewidths=3, linestyles='--', colors='cyan', transform=ccrs.PlateCarree(), zorder=5)  

tmpcs = ax.contour(lons, lats, temp_700, levels=[6, 8, 10, 12, 14, 16, 18], linewidths=1, linestyles='--', colors='k', transform=ccrs.PlateCarree(), zorder=5)  
plt.clabel(tmpcs, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
           rightside_up=True, use_clabeltext=True)

# plot 700hpa temperature fill
contourf = ax.contourf(lons, lats, temp_700, np.arange(-40, 42, 1), extent='both',
                 cmap=temp_cmap, alpha=1, transform=ccrs.PlateCarree(), zorder=4)

# plot 700hpa wind barbs
every = 15
barbs = ax.barbs(lons.values[0::every, 0::every], lats.values[0::every, 0::every],
                uwnd_700[0::every, 0::every], vwnd_700[0::every, 0::every],
                length=6.5, alpha=0.7, transform=ccrs.PlateCarree(), zorder=11)

# # plot title, add one to the left with model name and data names, add another to the right with time info
# plt.figtext(0.08, 1.03, f'     RAP 700 hPa Analysis | {valid_date[0:10]} {valid_date[11:-13]}z', weight='bold', ha='left', fontsize=20, color='white')
# plt.figtext(0.08, 1.00, f'     Heights (m), Temperature (C), Wind (kts)', ha='left', fontsize=18, color='white')
# plt.figtext(0.915, 1.04, f' ', ha='left', fontsize=20)
# plt.figtext(0.915, -0.01, f' ', ha='left', fontsize=20)
# cax = fig.add_axes([0.91, 0.024, 0.01, 0.95])
# cbar = fig.colorbar(contourf, cax=cax, orientation='vertical', ticks=np.arange(-40, 42, 5), extendrect=True)
# cax.text(3, 0.5, 'Temperature (C)', ha='left',va='center',rotation=270, color='white',fontsize=12,fontweight='bold',transform=cax.transAxes)
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

#composite_filename = build_filename("staged_figures/conus_rap_analysis/", f"conus_analysis", data_date.astype('datetime64[us]').item(), variant='700a')

composite_filename = build_filename("staged_figures/conus_rap_analysis/", f"conus_analysis", now_utc, variant='700a')

figure_builder(fig, ax,
    title=f"RAP Analysis • 700hPa",
    subtitle=f"Heights (m), Temperature (C), Wind (kts)",
    valid=f"Valid • {valid_date[0:10]} {valid_date[11:-13]}z",
    mappable=contourf,
    cbar_title="Temperature",
    cbar_units="°C",
    cbar_ticks=np.arange(-40, 42, 5),
    footer_left=f"RAP 13km • INIT {valid_date[0:10]} {valid_date[11:-13]}z",
    save_path=composite_filename)

print("    FINISHED 700HPA TEMP MAP")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################











#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
#################################
# BUILD 850HPA FLOW MAP
#################################
fig, ax = map_builder(projection=projection, extent=[-119, -74, 23.5, 53.5], terrain=True, state_color='navy', border_color='navy', border_factor=1.5)

# slice data
plev850 = np.where(pres_levs == 850)[0][0]
ghgt_850 = ghgt_iso[plev850]
uwnd_850 = uwnd_iso[plev850]
vwnd_850 = vwnd_iso[plev850]
temp_850 = temp_iso[plev850]
wdsp_850 = np.sqrt(uwnd_850**2 + vwnd_850**2)

# plot 850 hpa heights
contour = ax.contour(lons, lats, ghgt_850, np.arange(0, 1700, 30),
                colors='black', linewidths=3.0, linestyles='-',
                transform=ccrs.PlateCarree(), zorder=11)
plt.clabel(contour, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
           rightside_up=True, use_clabeltext=True)

# plot 850 wind speed fill
contourf = ax.contourf(lons, lats, wdsp_850, np.arange(20, 85, 5), extend='max',
                            cmap=wdsp_cmap, alpha=0.7, transform=ccrs.PlateCarree(), zorder=4)


# plot 850hpa wind barbs
every = 15
barbs = ax.barbs(lons.values[0::every, 0::every], lats.values[0::every, 0::every],
                uwnd_850[0::every, 0::every], vwnd_850[0::every, 0::every],
                length=6.5, alpha=0.7, transform=ccrs.PlateCarree(), zorder=11)

# # plot title, add one to the left with model name and data names, add another to the right with time info
# plt.figtext(0.08, 1.03, f'     RAP 850 hPa Analysis | {valid_date[0:10]} {valid_date[11:-13]}z', weight='bold', ha='left', fontsize=20, color='white')
# plt.figtext(0.08, 1.00, f'     Heights (m), Wind (kts)', ha='left', fontsize=18, color='white')
# plt.figtext(0.915, 1.04, f' ', ha='left', fontsize=20)
# plt.figtext(0.915, -0.01, f' ', ha='left', fontsize=20)
# cax = fig.add_axes([0.91, 0.024, 0.01, 0.95])
# cbar = fig.colorbar(contourf, cax=cax, orientation='vertical', ticks=np.arange(25, 100, 5), extendrect=True)
# cax.text(3, 0.5, 'Wind Speed (kts)', ha='left',va='center',rotation=270, color='white',fontsize=12,fontweight='bold',transform=cax.transAxes)
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

composite_filename = build_filename("staged_figures/conus_rap_analysis/", f"conus_analysis", now_utc, variant='850a')

figure_builder(fig, ax,
    title=f"RAP Analysis • 850hPa",
    subtitle=f"Heights (m), Wind (kt)",
    valid=f"Valid • {valid_date[0:10]} {valid_date[11:-13]}z",
    mappable=contourf,
    cbar_title="Wind Speed",
    cbar_units="kts",
    cbar_ticks=np.arange(25, 100, 5),
    footer_left=f"RAP 13km • INIT {valid_date[0:10]} {valid_date[11:-13]}z",
    save_path=composite_filename)

print("    FINISHED 850HPA FLOW MAP")
#############################################################################################################################################################################
#############################################################################################################################################################################
############################################################################################################################################################################.










#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
#################################
# BUILD 850HPA TEMP MAP
#################################
fig, ax = map_builder(projection=projection, extent=[-119, -74, 23.5, 53.5], terrain=True, state_color='navy', border_color='navy', border_factor=1.5)

# plot 850 hpa heights
contour = ax.contour(lons, lats, ghgt_850, np.arange(0, 1700, 30),
                colors='black', linewidths=3.0, linestyles='-',
                transform=ccrs.PlateCarree(), zorder=11)
plt.clabel(contour, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
           rightside_up=True, use_clabeltext=True)

# plot 0C isotherm
ax.contour(lons, lats, temp_850, levels=[0], linewidths=3, linestyles='--', colors='cyan', transform=ccrs.PlateCarree(), zorder=5)  

# plot 850hpa temperature fill
contourf = ax.contourf(lons, lats, temp_850, np.arange(-40, 42, 1), extent='both',
                 cmap=temp_cmap, alpha=1, transform=ccrs.PlateCarree(), zorder=4)

# plot 850hpa wind barbs
every = 15
barbs = ax.barbs(lons.values[0::every, 0::every], lats.values[0::every, 0::every],
                uwnd_850[0::every, 0::every], vwnd_850[0::every, 0::every],
                length=6.5, alpha=0.7, transform=ccrs.PlateCarree(), zorder=11)

# # plot title, add one to the left with model name and data names, add another to the right with time info
# plt.figtext(0.08, 1.03, f'     RAP 850 hPa Analysis | {valid_date[0:10]} {valid_date[11:-13]}z', weight='bold', ha='left', fontsize=20, color='white')
# plt.figtext(0.08, 1.00, f'     Heights (m), Temperature (C), Wind (kts)', ha='left', fontsize=18, color='white')
# plt.figtext(0.915, 1.04, f' ', ha='left', fontsize=20)
# plt.figtext(0.915, -0.01, f' ', ha='left', fontsize=20)
# cax = fig.add_axes([0.91, 0.024, 0.01, 0.95])
# cbar = fig.colorbar(contourf, cax=cax, orientation='vertical', ticks=np.arange(-40, 42, 5), extendrect=True)
# cax.text(3, 0.5, 'Temperature (C)', ha='left',va='center',rotation=270, color='white',fontsize=12,fontweight='bold',transform=cax.transAxes)
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

# composite_filename = build_filename("staged_figures/conus_rap_analysis/", f"conus_analysis", data_date.astype('datetime64[us]').item(), variant='850b')

composite_filename = build_filename("staged_figures/conus_rap_analysis/", f"conus_analysis", now_utc, variant='850b')

figure_builder(fig, ax,
    title=f"RAP Analysis • 850hPa",
    subtitle=f"Heights (m), Temperature (C), Wind (kts)",
    valid=f"Valid • {valid_date[0:10]} {valid_date[11:-13]}z",
    mappable=contourf,
    cbar_title="Temperature",
    cbar_units="°C",
    cbar_ticks=np.arange(-40, 42, 5),
    footer_left=f"RAP 13km • INIT {valid_date[0:10]} {valid_date[11:-13]}z",
    save_path=composite_filename)

print("    FINISHED 850HPA TEMP MAP")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################










#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
#################################
# BUILD 850HPA TADV MAP
#################################
try:
    fig, ax = map_builder(projection=projection, extent=[-119, -74, 23.5, 53.5], terrain=True, state_color='navy', border_color='navy', border_factor=1.5)

    # slice data
    plev850 = np.where(pres_levs == 850)[0][0]
    ghgt_850 = ghgt_iso[plev850]
    uwnd_850 = uwnd_iso[plev850]
    vwnd_850 = vwnd_iso[plev850]
    temp_850 = temp_iso[plev850]
    wdsp_850 = np.sqrt(uwnd_850**2 + vwnd_850**2)

    # plot 850hpa heights
    contour = ax.contour(lons, lats, ghgt_850, np.arange(0, 1700, 30),
                    colors='black', linewidths=3.0, linestyles='-',
                    transform=ccrs.PlateCarree(), zorder=11)
    plt.clabel(contour, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
            rightside_up=True, use_clabeltext=True)



    n_reps = 80
    # plot 850hpa tadv
    contourf = ax.contourf(raw_data['lon'], raw_data['lat'], 3*(mpcalc.smooth_n_point(adv[0,plev850,:,:], 9, n_reps)),
                    np.arange(-7,7.25,0.25), cmap='bwr', transform=ccrs.PlateCarree(), zorder=4, extend='both')

    ax.contour(lons, lats, temp_850, levels=[0], linewidths=3, linestyles='--', colors='gray', transform=ccrs.PlateCarree(), zorder=5)  

    # plot frontogenesis hatch (> 2 degC / 100km / 3hr)
    fgen_hatch = np.where(fgen_masked[0,plev850,:,:].values >= 2, 1, np.nan)
    ax.contourf(raw_data['lon'], raw_data['lat'], fgen_hatch,
                levels=[0.5, 1.5], colors='none', hatches=['////'], alpha=0,
                transform=ccrs.PlateCarree(), zorder=6)

    # plot 850hpa wind barbs
    every = 15
    barbs = ax.barbs(lons.values[0::every, 0::every], lats.values[0::every, 0::every],
                    uwnd_850[0::every, 0::every], vwnd_850[0::every, 0::every],
                    length=6.5, alpha=0.7, transform=ccrs.PlateCarree(), zorder=11)

    # # plot title, add one to the left with model name and data names, add another to the right with time info
    # plt.figtext(0.08, 1.03, f'     RAP 850 hPa Analysis | {valid_date[0:10]} {valid_date[11:-13]}z', weight='bold', ha='left', fontsize=20, color='white')
    # plt.figtext(0.08, 1.00, f'     Heights (m), 3hr Temperature Adv (C/3hr), Frontogenesis (>2'+u'\xb0'+'C / 100km / 3hr), Wind (kt)', ha='left', fontsize=18, color='white')
    # plt.figtext(0.915, 1.04, f' ', ha='left', fontsize=20)
    # plt.figtext(0.915, -0.01, f' ', ha='left', fontsize=20)
    # cax = fig.add_axes([0.91, 0.024, 0.01, 0.95])
    # cbar = fig.colorbar(tadv_contourf, cax=cax, orientation='vertical', ticks=np.arange(-7,7.25,1), extendrect=True)
    # cax.text(3, 0.5, 'Temperature Advection' + ' ('+u'\xb0'+'C / 3hr)', ha='left',va='center',rotation=270, color='white',fontsize=12,fontweight='bold',transform=cax.transAxes)
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

    # composite_filename = build_filename("staged_figures/conus_rap_analysis/", f"conus_analysis", data_date.astype('datetime64[us]').item(), variant='850c')

    composite_filename = build_filename("staged_figures/conus_rap_analysis/", f"conus_analysis", now_utc, variant='850c')

    figure_builder(fig, ax,
        title=f"RAP Analysis • 850hPa",
        subtitle=r"Heights (m), 3hr Temperature Adv (C/3hr), Frontogenesis ($>2^\circ\mathrm{C}\,(100\,\mathrm{km})^{-1}(3\,\mathrm{hr})^{-1}$), Wind (kt)",
        valid=f"Valid • {valid_date[0:10]} {valid_date[11:-13]}z", 
        mappable=contourf,
        cbar_title="Temperature Advection",
        cbar_units=r"$>2^\circ\mathrm{C}\,(100\,\mathrm{km})^{-1}(3\,\mathrm{hr})^{-1}$",
        cbar_ticks=np.arange(-7,7.25,1),
        footer_left=f"RAP 13km • INIT {valid_date[0:10]} {valid_date[11:-13]}z",
        save_path=composite_filename)

    print("    FINISHED 850HPA TEMP ADV MAP")

except TypeError as e:
    print("    ERROR ENCOUNTERED WHILE BUILDING 850HPA TEMP ADV MAP. SKIPPING....")
    pass
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################










#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
#################################
# SURFACE TEMPERATURE MAP
#################################
fig, ax = map_builder(projection=projection, extent=[-119, -74, 23.5, 53.5], terrain=True, state_color='navy', border_color='navy', border_factor=1.5)


# plot mslp
cs = ax.contour(lons, lats, pres_sfc/100, np.arange(904, 1054, 4), colors='black',
                linewidths=3.0, linestyles='-',
                transform=ccrs.PlateCarree(), zorder=11)
plt.clabel(cs, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
           rightside_up=True, use_clabeltext=True)

ax.contour(lons, lats, c2f(temp_sfc[0,:,:]), levels=[32], linewidths=3, linestyles='--', colors='cyan', transform=ccrs.PlateCarree(), zorder=5)  

contourf = ax.contourf(lons, lats, c2f(temp_sfc[0,:,:]), np.arange(-60, 121, 1), extent='both',
                 cmap=temp_cmap, alpha=1, transform=ccrs.PlateCarree(), zorder=4)

# plot  wind barbs
every = 15
barbs = ax.barbs(lons.values[0::every, 0::every], lats.values[0::every, 0::every],
                uwnd_sfc[0, 0::every, 0::every], vwnd_sfc[0, 0::every, 0::every],
                length=6.5, alpha=0.7, transform=ccrs.PlateCarree(), zorder=11)

# # plot title, add one to the left with model name and data names, add another to the right with time info
# plt.figtext(0.08, 1.03, f'   RAP Surface Analysis | {valid_date[0:10]} {valid_date[11:-13]}z', weight='bold', ha='left', fontsize=20, color='white')
# plt.figtext(0.08, 1.00, f'   MSLP (hPa), Temperature (F), Wind (kts)', ha='left', fontsize=18, color='white')
# plt.figtext(0.915, 1.04, f' ', ha='left', fontsize=20)
# plt.figtext(0.915, -0.01, f' ', ha='left', fontsize=20)
# cax = fig.add_axes([0.91, 0.024, 0.01, 0.95])
# cbar = fig.colorbar(contourf, cax=cax, orientation='vertical', ticks=np.arange(-60, 130, 10), extendrect=True)
# cax.text(3, 0.5, 'Temperature (F)', ha='left',va='center',rotation=270, color='white',fontsize=12,fontweight='bold',transform=cax.transAxes)
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

# composite_filename = build_filename("staged_figures/conus_rap_analysis/", f"conus_analysis", data_date.astype('datetime64[us]').item(), variant='000a')
# plt.savefig(composite_filename, bbox_inches="tight")

composite_filename = build_filename("staged_figures/conus_rap_analysis/", f"conus_analysis", now_utc, variant='000a')

figure_builder(fig, ax,
    title=f"RAP Analysis • Surface",
    subtitle=f"MSLP (hPa), Temperature (F), Wind (kts)",
    valid=f"Valid • {valid_date[0:10]} {valid_date[11:-13]}z",
    mappable=contourf,
    cbar_title="Temperature",
    cbar_units="°F",
    cbar_ticks=np.arange(-60, 130, 10),
    footer_left=f"RAP 13km • INIT {valid_date[0:10]} {valid_date[11:-13]}z",
    save_path=composite_filename)

print("    FINISHED SFC TEMP MAP")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################





#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
#################################
# SURFACE DEWPOINT MAP
#################################
fig, ax = map_builder(projection=projection, extent=[-119, -74, 23.5, 53.5], terrain=True, state_color='navy', border_color='navy', border_factor=1.5)


# plot mslp
cs = ax.contour(lons, lats, pres_sfc/100, np.arange(904, 1054, 4), colors='black',
                linewidths=3.0, linestyles='-',
                transform=ccrs.PlateCarree(), zorder=11)
plt.clabel(cs, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
           rightside_up=True, use_clabeltext=True)

dpcs = ax.contour(lons, lats, c2f(dwpt_sfc[0,:,:]), levels=[45, 50, 55, 60, 65, 70, 75], colors='k',
           linewidths=1, linestyles='--', transform=ccrs.PlateCarree(), zorder=5)  
plt.clabel(dpcs, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
           rightside_up=True, use_clabeltext=True)

levels = np.arange(-40, 85, 5)
norm = Normalize(vmin=-50, vmax=100)
contourf = ax.contourf(lons, lats, c2f(dwpt_sfc[0,:,:]), levels, norm=norm, extend='both',
                 cmap=dwpt_cmap, alpha=1, transform=ccrs.PlateCarree(), zorder=4)

# plot  wind barbs
every = 15
barbs = ax.barbs(lons.values[0::every, 0::every], lats.values[0::every, 0::every],
                uwnd_sfc[0, 0::every, 0::every], vwnd_sfc[0, 0::every, 0::every],
                length=6.5, alpha=0.7, transform=ccrs.PlateCarree(), zorder=11)

# # plot title, add one to the left with model name and data names, add another to the right with time info
# plt.figtext(0.08, 1.03, f'   RAP Surface Analysis | {valid_date[0:10]} {valid_date[11:-13]}z', weight='bold', ha='left', fontsize=20, color='white')
# plt.figtext(0.08, 1.00, f'    MSLP (hPa), Dewpoint (F), Wind (kts)', ha='left', fontsize=18, color='white')
# plt.figtext(0.915, 1.04, f' ', ha='left', fontsize=20)
# plt.figtext(0.915, -0.01, f' ', ha='left', fontsize=20)
# cax = fig.add_axes([0.91, 0.024, 0.01, 0.95])
# cbar = fig.colorbar(contourf, cax=cax, orientation='vertical', ticks=np.arange(-50, 90, 5), extendrect=True)
# cax.text(3, 0.5, 'Dewpoint temperature (F)', ha='left',va='center',rotation=270, color='white',fontsize=12,fontweight='bold',transform=cax.transAxes)
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

# composite_filename = build_filename("staged_figures/conus_rap_analysis/", f"conus_analysis", data_date.astype('datetime64[us]').item(), variant='000b')
# plt.savefig(composite_filename, bbox_inches="tight")

composite_filename = build_filename("staged_figures/conus_rap_analysis/", f"conus_analysis", now_utc, variant='000b')

figure_builder(fig, ax,
    title=f"RAP Analysis • Surface",
    subtitle=f"MSLP (hPa), Dewpoint (F), Wind (kts)",
    valid=f"Valid • {valid_date[0:10]} {valid_date[11:-13]}z",
    mappable=contourf,
    cbar_title="Dewpoint",
    cbar_units="°F",
    cbar_ticks=np.arange(-50, 90, 5),
    footer_left=f"RAP 13km • INIT {valid_date[0:10]} {valid_date[11:-13]}z",
    save_path=composite_filename)

print("    FINISHED SFC DWPT MAP")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################








#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
#################################
# SURFACE CAPE MAP
#################################
fig, ax = map_builder(projection=projection, extent=[-119, -74, 23.5, 53.5], terrain=True, state_color='navy', border_color='navy', border_factor=1.5)


# plot mslp
cs = ax.contour(lons, lats, pres_sfc/100, np.arange(904, 1054, 4), colors='black',
                linewidths=2, linestyles='-',
                transform=ccrs.PlateCarree(), alpha=1, zorder=11)
plt.clabel(cs, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
           rightside_up=True, use_clabeltext=True)



contourf = ax.contourf(lons, lats, cape_ml, levels=np.arange(100, 5000, 100), extend='max',
                 cmap=cape_cmap, transform=ccrs.PlateCarree(), alpha=0.7, zorder=4)




# plot storm-relative mean wind barbs (6-9km and 0-2km)
every = 15

mask = cape_ml > 50
# 6–9 km SRMW
barbs_6_9 = ax.barbs(lons.values[0::every, 0::every], lats.values[0::every, 0::every],
    np.where(mask[0::every, 0::every],srw_6_9_u.values[0::every, 0::every], np.nan),
    np.where(mask[0::every, 0::every],srw_6_9_v.values[0::every, 0::every], np.nan),
    length=5.5, color='darkblue', alpha=0.8, transform=ccrs.PlateCarree(), zorder=14)

# 0–2 km SRMW
barbs_0_2 = ax.barbs(lons.values[0::every, 0::every], lats.values[0::every, 0::every],
    np.where(mask[0::every, 0::every], srw_0_2_u.values[0::every, 0::every], np.nan),
    np.where(mask[0::every, 0::every], srw_0_2_v.values[0::every, 0::every], np.nan),
    length=5.5, color='darkred', alpha=1, transform=ccrs.PlateCarree(), zorder=14)


# plot MLCIN as a hatch
import matplotlib as mpl
mpl.rcParams['hatch.linewidth'] = 0.3
cin_hatch = np.where(cin_ml[:,:] < -40, 1, np.nan)
ax.contourf(lons, lats, cin_hatch, levels=[0.5, 1.5], colors='none', hatches=['////'],
    edgecolors='red', linewidths=0.05, transform=ccrs.PlateCarree(), zorder=13, alpha=0)



# # plot title, add one to the left with model name and data names, add another to the right with time info
# plt.figtext(0.08, 1.03, f'   RAP Surface Analysis | {valid_date[0:10]} {valid_date[11:-13]}z', weight='bold', ha='left', fontsize=20, color='white')
# plt.figtext(0.08, 1.00, f'    MSLP (hPa), MLCAPE (J/kg), MLCIN (<-40 J/kg), 6-9km (gray) & 0-2km (black) SRW Crossover', ha='left', fontsize=18, color='white')
# plt.figtext(0.915, 1.04, f' ', ha='left', fontsize=20)
# plt.figtext(0.915, -0.01, f' ', ha='left', fontsize=20)
# cax = fig.add_axes([0.91, 0.024, 0.01, 0.95])
# cbar = fig.colorbar(contourf, cax=cax, orientation='vertical', extendrect=True)
# # Use a subset of ticks so labels are legible and the bar matches the provided colorbar style
# cbar.set_ticks([100, 1000, 2000, 3000, 4000, 5000, 6000])
# cax.text(3.7, 0.5, '250hPa Mixed Layer Convective Available Potential Energy (J/kg)', ha='left',va='center',rotation=270, color='white',fontsize=12,fontweight='bold',transform=cax.transAxes)
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

# composite_filename = build_filename("staged_figures/conus_rap_analysis/", f"conus_analysis", data_date.astype('datetime64[us]').item(), variant='000c')
composite_filename = build_filename("staged_figures/conus_rap_analysis/", f"conus_analysis", now_utc, variant='000c')

figure_builder(fig, ax,
    title=f"RAP Analysis • Surface",
    subtitle=f"MSLP (hPa), MLCAPE (J/kg), MLCIN (<-40 J/kg), 6-9km (gray) & 0-2km (black) SRW Crossover",
    valid=f"Valid • {valid_date[0:10]} {valid_date[11:-13]}z",
    mappable=contourf,
    cbar_title="250hPa Mixed Layer Convective Available Potential Energy",
    cbar_units="J/kg",
    cbar_ticks=[100, 1000, 2000, 3000, 4000, 5000, 6000],
    footer_left=f"RAP 13km • INIT {valid_date[0:10]} {valid_date[11:-13]}z",
    save_path=composite_filename)

print("    FINISHED SFC CAPE MAP")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################







#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
#################################
# SURFACE OBS MAP
#################################
fig, ax = map_builder(projection=projection, extent=[-118, -74, 24, 52], terrain=True, terrain_zoom=6, state_color='white', border_color='white')


try:
    ###################################################################
    # METAR STATION PLOTS
    ###################################################################
    from metpy.plots import StationPlot, StationPlotLayout, sky_cover
    from matplotlib.patheffects import withStroke
    TEXT_OUTLINE = [withStroke(linewidth=3.0, foreground="#1A1A18")]
    BARB_OUTLINE = [withStroke(linewidth=3.5, foreground="#1A1A18")]

    custom_layout = StationPlotLayout()
    custom_layout.add_barb('eastward_wind', 'northward_wind', units='knots', path_effects=BARB_OUTLINE)
    custom_layout.add_value('NW', 'air_temperature', fmt='.0f', fontsize=5, color="#FF6B35", path_effects=TEXT_OUTLINE)
    custom_layout.add_value('SW', 'dew_point_temperature', fmt='.0f', fontsize=5, color="#8FE388", path_effects=TEXT_OUTLINE)
    custom_layout.add_symbol('C', 'cloud_coverage', sky_cover, path_effects=TEXT_OUTLINE)
    stationplot = StationPlot(ax, metar_obs['longitude'], metar_obs['latitude'], clip_on=True,
                            transform=ccrs.PlateCarree(), fontsize=6, zorder=12, alpha=0.8, color='white')
    custom_layout.plot(stationplot, metar_obs)
except: 
    pass



# plot mslp
cs = ax.contour(lons, lats, pres_sfc/100, np.arange(904, 1054, 4), colors='black',
                linewidths=2.0, linestyles='-',
                transform=ccrs.PlateCarree(), zorder=11)
plt.clabel(cs, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
           rightside_up=True, use_clabeltext=True)

# plot nexrad mosaic
pm = ax.pcolormesh(radar_lon+0.05, radar_lat+0.05, radar_data,
              vmin=-15, vmax=95, cmap=rs_expertreflect_cmap, alpha=0.8, zorder=1.3, transform=ccrs.PlateCarree())

# plot wpc fronts bulletin
texts, params, geoms, valid_time = plot_bulletin(ax)


# plot sat
ax.imshow(xrds_sat["CMI"].values, origin="upper", extent=sat_extent, transform=sat_crs,
        cmap="gray", norm=PowerNorm(gamma=0.55, vmin=0.0, vmax=1.1), interpolation="nearest",
        regrid_shape=700, alpha=0.90, zorder=1)
xrds_sat.close()



# # plot title, add one to the left with model name and data names, add another to the right with time info
# plt.figtext(0.08, 1.03, f'     RAP Surface Analysis | {valid_date[0:10]} {valid_date[11:-13]}z', weight='bold', ha='left', fontsize=20, color='white')
# plt.figtext(0.08, 1.00, f'     RAP MSLP (hPa), {metar_time[11:16]}z METARs, {valid_time}z WPC Fronts, {str(radar_time)[11:16]}z Reflectivity Mosaic, {sat_time[0:2]}:{sat_time[2:4]}z GOES19 Radiance', ha='left', fontsize=18, color='white')
# # plt.figtext(0.915, 1.04, f' ', ha='left', fontsize=20)
# # # colorbar for filled contour
# # cbar = plt.colorbar(pm, aspect=70, fraction=0.02, ax=ax, orientation='horizontal', pad=-0.01, extendrect=True)
# # cbar.set_label('Reflectivity (dBz)',  fontsize=15, color='white', fontweight='bold')
# # cbar.ax.tick_params(labelcolor='white')
# # for t in cbar.ax.get_xticklabels():
# #     t.set_fontweight('bold')
# plt.figtext(0.915, 1.04, f' ', ha='left', fontsize=20)
# plt.figtext(0.915, -0.01, f' ', ha='left', fontsize=20)
# cax = fig.add_axes([0.91, 0.024, 0.01, 0.95])
# cbar = fig.colorbar(pm, cax=cax, orientation='vertical', ticks=np.arange(-15, 95, 5), extendrect=True)
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

# composite_filename = build_filename("staged_figures/conus_rap_analysis/", f"conus_analysis", data_date.astype('datetime64[us]').item(), variant='000d')
composite_filename = build_filename("staged_figures/conus_rap_analysis/", f"conus_analysis", now_utc, variant='000d')

figure_builder(fig, ax,
    title=f"RAP Analysis • Surface",
    subtitle=F'RAP MSLP (hPa), Surface Observations, WPC Fronts, MRMS Reflectivity Mosaic (dBz), GOES19 Band 02 Visible',
    valid=f"Valid • {valid_date[0:10]} {valid_date[11:-13]}z",
    mappable=pm,
    cbar_title="Reflectivity",
    cbar_units="dBz",
    cbar_ticks=np.arange(-15, 95, 5),
    footer_left=f"RAP 13km  •  INIT {valid_date[0:10]} {valid_date[11:-13]}z  •  {metar_time[11:16]}z METARs, {valid_time}z WPC Fronts, {str(radar_time)[11:16]}z Reflectivity Mosaic, {sat_time_str} {sat_type}",
    save_path=composite_filename)

print("    FINISHED SFC ANL MAP")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################



elapsed_time = comp_time.time() - st
print(f"############\nSCRIPT FINISHED: time: {comp_time.strftime('%H:%M:%S', comp_time.gmtime(elapsed_time))}\n############")

