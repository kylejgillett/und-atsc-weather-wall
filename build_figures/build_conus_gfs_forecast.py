##########################################################
#         GFS FORECAST ANALYSIS FIGURES SCRIPT
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
from datetime import datetime, timezone, timedelta
import numpy as np
import sys
import os
import xarray as xr
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm


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
from get_data.get_gfs_data import gfs_forecast
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





#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################

# set up rap retrieval 
box_size   = 50

# pull rap data
#forecast_hours = [int(hr) for hr in np.arange(6, 78, 6)]
datas = gfs_forecast(box_size=box_size, forecast_hours=[21])

for fh, raw_data in datas.items():


    # LATS & LONS
    lats = raw_data.variables['latitude'][:]
    lons = raw_data.variables['longitude'][:]

    # PRES LEVS
    pres_levs = raw_data['isobaric'][:]
    pres_levs = pres_levs / 100

    # DATE INFO
    fh = str(fh).zfill(3)

    # extract date objects and set up valid date title
    run_date = f'{raw_data['reftime'].values}'
    valid_date = f'{raw_data['reftime'].values.astype('datetime64[ms]').astype(datetime) + timedelta(hours=int(fh))}'
    valid_day_name = day_to_abbrev(raw_data['reftime'].values.astype('datetime64[ms]').astype(datetime) + timedelta(hours=int(fh)))
    valid_date_str = f"+ {fh}hr valid {valid_day_name} {valid_date[5:7]}/{valid_date[8:10]} {valid_date[-8:-6]}z"
    print(valid_date)

    # BASIC DATA EXTRACTION
    sigma = 1.5
    ghgt_iso = ndimage.gaussian_filter(raw_data.variables['Geopotential_height_isobaric'], sigma=sigma)
    temp_iso = ndimage.gaussian_filter(raw_data.variables['Temperature_isobaric'], sigma=sigma) - 273.15
    uwnd_iso = ndimage.gaussian_filter(raw_data.variables['u-component_of_wind_isobaric'], sigma=sigma) * 1.94384
    vwnd_iso = ndimage.gaussian_filter(raw_data.variables['v-component_of_wind_isobaric'], sigma=sigma) * 1.94384

    pres_sfc = ndimage.gaussian_filter(raw_data.variables['MSLP_Eta_model_reduction_msl'], sigma=sigma)
    temp_sfc = ndimage.gaussian_filter(raw_data.variables['Temperature_height_above_ground'], sigma=sigma) - 273.15
    uwnd_sfc = ndimage.gaussian_filter(raw_data.variables['u-component_of_wind_height_above_ground'], sigma=sigma) * 1.94384
    vwnd_sfc = ndimage.gaussian_filter(raw_data.variables['v-component_of_wind_height_above_ground'], sigma=sigma) * 1.94384
    relh_sfc = ndimage.gaussian_filter(raw_data.variables['Relative_humidity_height_above_ground'], sigma=sigma)
    dwpt_sfc = mpcalc.dewpoint_from_relative_humidity(temp_sfc*units.degC, relh_sfc*units.percent)
    reft_sfc = ndimage.gaussian_filter(raw_data.variables['Composite_reflectivity_entire_atmosphere'], 0.01)
    rn_sfc = ndimage.gaussian_filter(raw_data.variables['Categorical_Rain_surface'], 0.75)
    sn_sfc = ndimage.gaussian_filter(raw_data.variables['Categorical_Snow_surface'], 0.75)
    zr_sfc = ndimage.gaussian_filter(raw_data.variables['Categorical_Freezing_Rain_surface'], 0.75)
    ip_sfc = ndimage.gaussian_filter(raw_data.variables['Categorical_Ice_Pellets_surface'], 0.75)

    #################################
    # CALCULATE FRONTOGENESIS
    ################################# 
    # using MetPy, compute frontogensis
    fgen = mpcalc.frontogenesis(raw_data['Temperature_isobaric'], 
                                raw_data['u-component_of_wind_isobaric'],
                                raw_data['v-component_of_wind_isobaric'],
                                latitude=raw_data['latitude'].values,
                                longitude=raw_data['longitude'].values,
                                crs=ccrs.PlateCarree())
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
                        latitude=raw_data['latitude'].values,
                        longitude=raw_data['longitude'].values,
                        crs=ccrs.PlateCarree())
    # convert units to delta deg C / hr
    adv = adv.metpy.convert_units('delta_degC/hour')
    # apply some smoothing to adv 
    adv = ndimage.gaussian_filter(adv, sigma=2, order=0) * units('K/sec')


    # use metpy to compute theta & add it into the `rap-data` DataSet
    raw_data['theta'] = mpcalc.potential_temperature(raw_data['isobaric'],raw_data['Temperature_isobaric'])

    # use metpy to compute latitude / longitude grid deltas (dx, dy) for PV calculation
    dx, dy = mpcalc.lat_lon_grid_deltas(raw_data['longitude'].values, raw_data['latitude'].values)

    # use metpy to compute PV and add it into the `rap_data` DataSet
    raw_data['pv'] = mpcalc.potential_vorticity_baroclinic(raw_data['theta'],
                                                            raw_data['isobaric'],
                                                            u=raw_data['u-component_of_wind_isobaric'],
                                                            v=raw_data['v-component_of_wind_isobaric'],
                                                            dx=dx[None, :, :], dy=dy[None, :, :],
                                                            latitude=raw_data['latitude'])

    thta_on_2pvu = interpolate_to_isosurface(raw_data['pv'].values, raw_data['theta'].values,  2*1e-6, bottom_up_search=True)
    u_on_2pvu    = interpolate_to_isosurface(raw_data['pv'][:,:,:].values, raw_data['u-component_of_wind_isobaric'][:,:,:].values, 2*1e-6, bottom_up_search=True)
    v_on_2pvu    = interpolate_to_isosurface(raw_data['pv'][:,:,:].values, raw_data['v-component_of_wind_isobaric'][:,:,:].values, 2*1e-6, bottom_up_search=True)

    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #############################################################################################################################################################################










    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #############################################################################################################################################################################

    # build map function
    def build_map(extent=[-122, -73, 21, 56], projection=ccrs.LambertConformal(), style='light', add_sat=False):
        fig = plt.figure(figsize=(20, 10))
        fig.set_facecolor('#009946')
        ax = plt.axes(projection=projection)

        ax.set_extent(extent)
        ax.set_box_aspect(0.6)

        if style == 'light':
            color = 'gray'
            alpha = 0.5
        else:
            color = 'black'
            alpha = 0.8

        ax.add_feature(cfeature.STATES, edgecolor='navy', alpha=0.4, linestyle='-', linewidth=3, zorder=10)
        ax.add_feature(cfeature.LAND, facecolor=color, alpha=alpha, zorder=0.1)
        ax.add_feature(cfeature.OCEAN, facecolor=color, alpha=alpha + 0.2, zorder=0)
        ax.add_feature(cfeature.COASTLINE, color='navy', alpha=0.4, linestyle='-', linewidth=3, zorder=11)
        if add_sat:
            from cartopy.io import img_tiles
            satellite = img_tiles.GoogleTiles(style='satellite')
            ax.add_image(satellite, 4)

        plt.tight_layout()

        return fig, ax
    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #############################################################################################################################################################################





    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #################################
    # BUILD 300 HPA MAP
    #################################
    fig, ax = build_map(add_sat=True)

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
    every = 10
    barbs = ax.barbs(lons.values[0::every], lats.values[ 0::every],
                    uwnd_300[0::every, 0::every], vwnd_300[0::every, 0::every],
                    length=6.5, alpha=0.7, transform=ccrs.PlateCarree(), zorder=12)

    # plot title, add one to the left with model name and data names, add another to the right with time info
    plt.figtext(0.08, 1.03, f'     {run_date[11:-16]}z GFS Forecast | {valid_date_str}', weight='bold', ha='left', fontsize=20, color='white')
    plt.figtext(0.08, 1.00, f'     Heights (m), Wind (kt)', ha='left', fontsize=18, color='white')
    plt.figtext(0.915, 1.04, f' ', ha='left', fontsize=20)
    plt.figtext(0.915, -0.01, f' ', ha='left', fontsize=20)
    cax = fig.add_axes([0.91, 0.024, 0.01, 0.95])
    cbar = fig.colorbar(contourf, cax=cax, orientation='vertical', ticks=np.arange(50, 160, 5)[::1], extendrect=True)
    cax.text(3, 0.5, f'Wind Speed (kts)',ha='left',va='center',rotation=270, color='white',fontsize=12,fontweight='bold',transform=cax.transAxes)
    cbar.ax.tick_params(axis='y', labelcolor='white') 
    for t in cbar.ax.get_yticklabels():
        t.set_fontweight('bold')
        t.set_fontsize(9)
    cbar.ax.set_facecolor('black')

    # colorbar for filled contour
    # cbar = plt.colorbar(contourf, aspect=70, fraction=0.02, ax=ax, orientation='horizontal', pad=-0.01, extendrect=True)
    # cbar.set_label('Wind Speed (kts)', fontsize=15, color='white', fontweight='bold')
    # cbar.ax.tick_params(labelcolor='white')
    # for t in cbar.ax.get_xticklabels():
    #     t.set_fontweight('bold')

    # add UND logo
    from PIL import Image
    img = Image.open('utils/images/und-logo.png')
    #                  side-side  up-down  size   size
    imgax = fig.add_axes([0.83, 1.01, 0.06, 0.06], anchor='SE', zorder=3)
    plt.figtext(0.81, 0.995, f'ATMOSPHERIC SCIENCES', ha='left', weight='bold', fontsize=10, color='white')
    imgax.imshow(img)
    imgax.axis('off')

    plt.savefig(f"staged_figures/conus_gfs_forecasts/gfs_300_flow_F{fh}.png", bbox_inches="tight")

    print("    FINISHED 300HPA FLOW MAP")
    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #############################################################################################################################################################################










    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #################################
    # BUILD 300HPA PVA MAP
    #################################
    fig, ax = build_map()

    # use 300hpa slices from above

    # plot 300hpa heights
    contour = ax.contour(lons, lats, ghgt_300, np.arange(0, 12000, 60),
                    colors='black', linewidths=3.0, linestyles='-',
                    transform=ccrs.PlateCarree(), zorder=11)
    plt.clabel(contour, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
            rightside_up=True, use_clabeltext=True)


    # plot 300hpa pv fill
    contourf = ax.contourf(raw_data['longitude'], raw_data['latitude'], raw_data['pv'][plev300,:,:]*1e6, pv_clevs, cmap=pv_cmap,
                    transform=ccrs.PlateCarree(),extend='both')

    # plot a single dashed contour @ 2PVU
    pv_contour = ax.contour(raw_data['longitude'], raw_data['latitude'], raw_data['pv'][plev300,:,:]*1e6, [2], colors='navy',linestyles='dashed',linewidths=2,
                    transform=ccrs.PlateCarree())

    # plot 300 hpa wind barbs
    every = 10
    barbs = ax.barbs(lons.values[0::every], lats.values[ 0::every],
                    uwnd_300[0::every, 0::every], vwnd_300[0::every, 0::every],
                    length=6.5, alpha=0.7, transform=ccrs.PlateCarree(), zorder=11)

    # plot title, add one to the left with model name and data names, add another to the right with time info
    plt.figtext(0.08, 1.03, f'     {run_date[11:-16]}z GFS Forecast | {valid_date_str}', weight='bold', ha='left', fontsize=20, color='white')
    plt.figtext(0.08, 1.00, f'     300 hPa Heights (m), Potential Vorticity (PVU), Wind (kt)', ha='left', fontsize=18, color='white')
    plt.figtext(0.915, 1.04, f' ', ha='left', fontsize=20)
    plt.figtext(0.915, -0.01, f' ', ha='left', fontsize=20)
    cax = fig.add_axes([0.91, 0.024, 0.01, 0.95])
    cbar = fig.colorbar(contourf, cax=cax, orientation='vertical', ticks=pv_clevs[::5], extendrect=True)
    cax.text(3, 0.5, r'Potential Vorticity Units (PVU; $\rm{10^{-6}\ K\ kg^{-1}\ m^{2}\ s^{-1}})$' + ' | 2PVU (dashed)', ha='left',va='center',rotation=270, color='white',fontsize=12,fontweight='bold',transform=cax.transAxes)
    cbar.ax.tick_params(axis='y', labelcolor='white') 
    for t in cbar.ax.get_yticklabels():
        t.set_fontweight('bold')
        t.set_fontsize(9)
    cbar.ax.set_facecolor('black')

    # # colorbar for filled contour
    # cbar = plt.colorbar(contourf, aspect=70, fraction=0.02, ax=ax, orientation='horizontal', pad=-0.01, extendrect=True)
    # cbar.set_label(r'Potential Vorticity Units (PVU; $\rm{10^{-6}\ K\ kg^{-1}\ m^{2}\ s^{-1}})$' + ' | 2PVU (dashed)',  fontsize=15, color='white', fontweight='bold')
    # cbar.ax.tick_params(labelcolor='white')
    # for t in cbar.ax.get_xticklabels():
    #     t.set_fontweight('bold')

    # add UND logo
    from PIL import Image
    img = Image.open('utils/images/und-logo.png')
    #                  side-side  up-down  size   size
    imgax = fig.add_axes([0.83, 1.01, 0.06, 0.06], anchor='SE', zorder=3)
    plt.figtext(0.81, 0.995, f'ATMOSPHERIC SCIENCES', ha='left', weight='bold', fontsize=10, color='white')
    imgax.imshow(img)
    imgax.axis('off')

    plt.savefig(f"staged_figures/conus_gfs_forecasts/gfs_300_pv_F{fh}.png", bbox_inches="tight")

    print("    FINISHED 300HPA PVA MAP")
    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #############################################################################################################################################################################








    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #################################
    # BUILD 500HPA FLOW MAP
    #################################
    fig, ax = build_map(add_sat=True)

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
    every = 10
    barbs = ax.barbs(lons.values[0::every], lats.values[ 0::every],
                    uwnd_500[0::every, 0::every], vwnd_500[0::every, 0::every],
                    length=6.5, alpha=0.7, transform=ccrs.PlateCarree(), zorder=12)

    # plot title, add one to the left with model name and data names, add another to the right with time info
    plt.figtext(0.08, 1.03, f'     {run_date[11:-16]}z GFS Forecast | {valid_date_str}', weight='bold', ha='left', fontsize=20, color='white')
    plt.figtext(0.08, 1.00, f'     500 hPa Heights (m), Wind (kt)', ha='left', fontsize=18, color='white')
    plt.figtext(0.915, 1.04, f' ', ha='left', fontsize=20)
    plt.figtext(0.915, -0.01, f' ', ha='left', fontsize=20)
    cax = fig.add_axes([0.91, 0.024, 0.01, 0.95])
    cbar = fig.colorbar(contourf, cax=cax, orientation='vertical', ticks=np.arange(30, 140, 5), extendrect=True)
    cax.text(3, 0.5, 'Wind Speed (kts)', ha='left',va='center',rotation=270, color='white',fontsize=12,fontweight='bold',transform=cax.transAxes)
    cbar.ax.tick_params(axis='y', labelcolor='white') 
    for t in cbar.ax.get_yticklabels():
        t.set_fontweight('bold')
        t.set_fontsize(9)
    cbar.ax.set_facecolor('black')

    # # colorbar for filled contour
    # cbar = plt.colorbar(contourf, aspect=70, fraction=0.02, ax=ax, orientation='horizontal', pad=-0.01, extendrect=True)
    # cbar.set_label('Wind Speed (kts)',  fontsize=15, color='white', fontweight='bold')
    # cbar.ax.tick_params(labelcolor='white')
    # for t in cbar.ax.get_xticklabels():
    #     t.set_fontweight('bold')


    # add UND logo
    from PIL import Image
    img = Image.open('utils/images/und-logo.png')
    #                  side-side  up-down  size   size
    imgax = fig.add_axes([0.83, 1.01, 0.06, 0.06], anchor='SE', zorder=3)
    plt.figtext(0.81, 0.995, f'ATMOSPHERIC SCIENCES', ha='left', weight='bold', fontsize=10, color='white')
    imgax.imshow(img)
    imgax.axis('off')

    plt.savefig(f"staged_figures/conus_gfs_forecasts/gfs_500_flow_F{fh}.png", bbox_inches="tight")

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
    fig, ax = build_map()

    n_reps = 80

    # compute vorticity and vorticity advection
    dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)
    f = mpcalc.coriolis_parameter(np.deg2rad(lats)).to('1/s')
    vor_500 = mpcalc.smooth_n_point(mpcalc.vorticity(uwnd_500*units.kts, vwnd_500*units.kts, dx=dx, dy=dy), 9, n_reps)
    #avor_500 = vor_500 + f

    # plot 500 hpa heights
    contour = ax.contour(lons, lats, ghgt_500, np.arange(3000, 7000, 60),
                    colors='black', linewidths=3.0, linestyles='-',
                    transform=ccrs.PlateCarree(), zorder=11)
    plt.clabel(contour, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
            rightside_up=True, use_clabeltext=True)

    # plot relative vorticity fill
    #avor_500 = np.where((vor_500.m * 10**5 > -5) & (vor_500.m * 10**5 < 15), np.nan, vor_500)
    norm = mcolors.TwoSlopeNorm(vmin=-30, vcenter=0, vmax=50)
    vort_cf = ax.contourf(lons, lats, vor_500 * 10**5, np.arange(-30, 52, 2), 
                        norm=norm, extend='both', cmap='PuOr_r', zorder=5, alpha=1, transform=ccrs.PlateCarree())

    # plot wind barbs
    every = 10
    barbs = ax.barbs(lons.values[0::every], lats.values[ 0::every],
                    uwnd_500[0::every, 0::every], vwnd_500[0::every, 0::every],
                    length=6.5, alpha=0.7, transform=ccrs.PlateCarree(), zorder=12)

    # plot title, add one to the left with model name and data names, add another to the right with time info
    plt.figtext(0.08, 1.03, f'     {run_date[11:-16]}z GFS Forecast | {valid_date_str}', weight='bold', ha='left', fontsize=20, color='white')
    plt.figtext(0.08, 1.00, f'     500 hPa Heights (m), Rel. Vorticity (/sec•10⁵), Wind (kt)', ha='left', fontsize=18, color='white')
    plt.figtext(0.915, 1.04, f' ', ha='left', fontsize=20)
    plt.figtext(0.915, -0.01, f' ', ha='left', fontsize=20)
    cax = fig.add_axes([0.91, 0.024, 0.01, 0.95])
    cbar = fig.colorbar(vort_cf, cax=cax, orientation='vertical', ticks=np.arange(-30, 52, 2), extendrect=True)
    cax.text(3, 0.5, 'Relative Vorticity (/sec•10⁵)', ha='left',va='center',rotation=270, color='white',fontsize=12,fontweight='bold',transform=cax.transAxes)
    cbar.ax.tick_params(axis='y', labelcolor='white') 
    for t in cbar.ax.get_yticklabels():
        t.set_fontweight('bold')
        t.set_fontsize(9)
    cbar.ax.set_facecolor('black')
    # # colorbar for filled contour
    # cbar = plt.colorbar(vort_cf, aspect=70, fraction=0.02, ax=ax, orientation='horizontal', pad=-0.01, extendrect=True)
    # cbar.set_label('Relative Vorticity (/sec•10⁵)',  fontsize=15, color='white', fontweight='bold')
    # cbar.ax.tick_params(labelcolor='white')
    # for t in cbar.ax.get_xticklabels():
    #     t.set_fontweight('bold')



    # add UND logo
    from PIL import Image
    img = Image.open('utils/images/und-logo.png')
    #                  side-side  up-down  size   size
    imgax = fig.add_axes([0.83, 1.01, 0.06, 0.06], anchor='SE', zorder=3)
    plt.figtext(0.81, 0.995, f'ATMOSPHERIC SCIENCES', ha='left', weight='bold', fontsize=10, color='white')
    imgax.imshow(img)
    imgax.axis('off')

    plt.savefig(f"staged_figures/conus_gfs_forecasts/gfs_500_relvort_F{fh}.png", bbox_inches="tight")

    print("    FINISHED 500HPA REL VORT MAP")
    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #############################################################################################################################################################################








    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #################################
    # BUILD 500HPA ABSVORTADV MAP
    #################################
    fig, ax = build_map()

    # use 500hpa data and vort calculations from above 
    relvort_adv = mpcalc.advection(vor_500, uwnd_500, vwnd_500, dx=dx, dy=dy) *1e9
    #absvort_adv = mpcalc.advection(avor_500, uwnd_500, vwnd_500, dx=dx, dy=dy) *1e9

    # plot 500hpa heights
    contour = ax.contour(lons, lats, ghgt_500, np.arange(3000, 7000, 60),
                    colors='black', linewidths=3.0, linestyles='-',
                    transform=ccrs.PlateCarree(), zorder=11)
    plt.clabel(contour, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
            rightside_up=True, use_clabeltext=True)

    # plot 500hpa abs vort adv fill
    vortadv_cf = ax.contourf(lons, lats, relvort_adv, np.arange(-40, 42, 2),               #np.arange(-6*12**-7, 6*12**-7, 1*10**-9),
                                extend='both', cmap='bwr', zorder=5, alpha=1, transform=ccrs.PlateCarree())

    # plot 500hpa wind barbs
    every = 10
    barbs = ax.barbs(lons.values[0::every], lats.values[ 0::every],
                    uwnd_500[0::every, 0::every], vwnd_500[0::every, 0::every],
                    length=6.5, alpha=0.7, transform=ccrs.PlateCarree(), zorder=12)

    # plot title, add one to the left with model name and data names, add another to the right with time info
    plt.figtext(0.08, 1.03, f'     {run_date[11:-16]}z GFS Forecast | {valid_date_str}', weight='bold', ha='left', fontsize=20, color='white')
    plt.figtext(0.08, 1.00, f'     500 hPa Heights (m), Rel. Vorticity Adv. (sec⁻²•10⁹), Wind (kt)', ha='left', fontsize=18, color='white')
    plt.figtext(0.915, 1.04, f' ', ha='left', fontsize=20)
    plt.figtext(0.915, -0.01, f' ', ha='left', fontsize=20)
    cax = fig.add_axes([0.91, 0.024, 0.01, 0.95])
    cbar = fig.colorbar(vortadv_cf, cax=cax, orientation='vertical', ticks=np.arange(-40, 42, 2), extendrect=True)
    cax.text(3, 0.5, 'Relative Vorticity Advection (sec⁻²•10⁹)', ha='left',va='center',rotation=270, color='white',fontsize=12,fontweight='bold',transform=cax.transAxes)
    cbar.ax.tick_params(axis='y', labelcolor='white') 
    for t in cbar.ax.get_yticklabels():
        t.set_fontweight('bold')
        t.set_fontsize(9)
    cbar.ax.set_facecolor('black')
    # # colorbar for filled contour
    # cbar = plt.colorbar(vortadv_cf, aspect=70, fraction=0.02, ax=ax, orientation='horizontal', pad=-0.01, extendrect=True)
    # cbar.set_label('Relative Vorticity Advection (sec⁻²•10⁹)',  fontsize=15, color='white', fontweight='bold')
    # cbar.ax.tick_params(labelcolor='white')
    # for t in cbar.ax.get_xticklabels():
    #     t.set_fontweight('bold')

    # add UND logo
    from PIL import Image
    img = Image.open('utils/images/und-logo.png')
    #                  side-side  up-down  size   size
    imgax = fig.add_axes([0.83, 1.01, 0.06, 0.06], anchor='SE', zorder=3)
    plt.figtext(0.81, 0.995, f'ATMOSPHERIC SCIENCES', ha='left', weight='bold', fontsize=10, color='white')
    imgax.imshow(img)
    imgax.axis('off')

    plt.savefig(f"staged_figures/conus_gfs_forecasts/gfs_500_relvortadv_F{fh}.png", bbox_inches="tight")

    print("    FINISHED 500HPA REL VORT ADV MAP")
    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #############################################################################################################################################################################







    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #################################
    # BUILD 850HPA TADV MAP
    #################################
    fig, ax = build_map()

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



    n_reps = 20
    # plot 850hpa tadv
    tadv_contourf = ax.contourf(raw_data['longitude'], raw_data['latitude'], 3*(mpcalc.smooth_n_point(adv[plev850,:,:], 9, n_reps)),
                    np.arange(-7,7.25,0.25), cmap='bwr', transform=ccrs.PlateCarree(), zorder=4, extend='both')

    ax.contour(lons, lats, temp_850, levels=[0], linewidths=3, linestyles='--', colors='gray', transform=ccrs.PlateCarree(), zorder=5)  

    # plot filled contours pf frontogenesis > 2 delta deg C / hr
    fgen_contourf = ax.contour(raw_data['longitude'], raw_data['latitude'], fgen_masked[plev850,:,:], 
                            np.arange(1, 32, 2), colors='navy', linestyles='-',
                            transform=ccrs.PlateCarree(), zorder=4)

    # plot 850hpa wind barbs
    every = 10
    barbs = ax.barbs(lons.values[0::every], lats.values[ 0::every],
                    uwnd_850[0::every, 0::every], vwnd_850[0::every, 0::every],
                    length=6.5, alpha=0.7, transform=ccrs.PlateCarree(), zorder=11)

    # plot title, add one to the left with model name and data names, add another to the right with time info
    plt.figtext(0.08, 1.03, f'     {run_date[11:-16]}z GFS Forecast | {valid_date_str}', weight='bold', ha='left', fontsize=20, color='white')
    plt.figtext(0.08, 1.00, f'     850 hPa Heights (m), 3hr Temperature Adv (C/3hr), Frontogenesis (>2'+u'\xb0'+'C / 100km / 3hr), Wind (kt)', ha='left', fontsize=18, color='white')
    plt.figtext(0.915, 1.04, f' ', ha='left', fontsize=20)
    plt.figtext(0.915, -0.01, f' ', ha='left', fontsize=20)
    cax = fig.add_axes([0.91, 0.024, 0.01, 0.95])
    cbar = fig.colorbar(tadv_contourf, cax=cax, orientation='vertical', ticks=np.arange(-7,7.25,1), extendrect=True)
    cax.text(3, 0.5, 'Temperature Advection' + ' ('+u'\xb0'+'C / 3hr)', ha='left',va='center',rotation=270, color='white',fontsize=12,fontweight='bold',transform=cax.transAxes)
    cbar.ax.tick_params(axis='y', labelcolor='white') 
    for t in cbar.ax.get_yticklabels():
        t.set_fontweight('bold')
        t.set_fontsize(9)
    cbar.ax.set_facecolor('black')
    # # colorbar for filled contour
    # cbar = plt.colorbar(tadv_contourf, aspect=70, fraction=0.02, ax=ax, orientation='horizontal', pad=-0.01, extendrect=True)
    # cbar.set_label('Temperature Advection' + ' ('+u'\xb0'+'C / 3hr)',  fontsize=15, color='white', fontweight='bold')
    # cbar.ax.tick_params(labelcolor='white')
    # for t in cbar.ax.get_xticklabels():
    #     t.set_fontweight('bold')

    # add UND logo
    from PIL import Image
    img = Image.open('utils/images/und-logo.png')
    #                  side-side  up-down  size   size
    imgax = fig.add_axes([0.83, 1.01, 0.06, 0.06], anchor='SE', zorder=3)
    plt.figtext(0.81, 0.995, f'ATMOSPHERIC SCIENCES', ha='left', weight='bold', fontsize=10, color='white')
    imgax.imshow(img)
    imgax.axis('off')

    plt.savefig(f"staged_figures/conus_gfs_forecasts/gfs_850_tempadv_F{fh}.png", bbox_inches="tight")

    print("    FINISHED 850HPA TEMP ADV MAP")
    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #############################################################################################################################################################################







    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #################################
    # 850 TEMP MAP
    #################################
    fig, ax = build_map()

    # use 850 slices from above


    # plot 850 hpa heights
    contour = ax.contour(lons, lats, ghgt_850, np.arange(0, 1700, 30),
                    colors='black', linewidths=3.0, linestyles='-',
                    transform=ccrs.PlateCarree(), zorder=11)
    plt.clabel(contour, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
            rightside_up=True, use_clabeltext=True)

    ax.contour(lons, lats, temp_850, levels=[0], linewidths=3, linestyles='--', colors='cyan', transform=ccrs.PlateCarree(), zorder=5)  

    # contourf = ax.contourf(lons, lats, wdsp_850, np.arange(25, 100, 5),
    #                  cmap='BuPu', alpha=0.7, transform=ccrs.PlateCarree(), zorder=4)
    contourf = ax.contourf(lons, lats, temp_850, np.arange(-40, 42, 1), extent='both',
                    cmap=temp_cmap, alpha=1, transform=ccrs.PlateCarree(), zorder=4)

    # plot 850hpa wind barbs
    every = 10
    barbs = ax.barbs(lons.values[0::every], lats.values[ 0::every],
                    uwnd_850[0::every, 0::every], vwnd_850[0::every, 0::every],
                    length=6.5, alpha=0.7, transform=ccrs.PlateCarree(), zorder=11)

    # plot title, add one to the left with model name and data names, add another to the right with time info
    plt.figtext(0.08, 1.03, f'     {run_date[11:-16]}z GFS Forecast | {valid_date_str}', weight='bold', ha='left', fontsize=20, color='white')
    plt.figtext(0.08, 1.00, f'     850 hPa Heights (m), Temperature (C), Wind (kts)', ha='left', fontsize=18, color='white')
    plt.figtext(0.915, 1.04, f' ', ha='left', fontsize=20)
    plt.figtext(0.915, -0.01, f' ', ha='left', fontsize=20)
    cax = fig.add_axes([0.91, 0.024, 0.01, 0.95])
    cbar = fig.colorbar(contourf, cax=cax, orientation='vertical', ticks=np.arange(-40, 42, 5), extendrect=True)
    cax.text(3, 0.5, 'Temperature (C)', ha='left',va='center',rotation=270, color='white',fontsize=12,fontweight='bold',transform=cax.transAxes)
    cbar.ax.tick_params(axis='y', labelcolor='white') 
    for t in cbar.ax.get_yticklabels():
        t.set_fontweight('bold')
        t.set_fontsize(9)
    cbar.ax.set_facecolor('black')
    # # colorbar for filled contour
    # cbar = plt.colorbar(contourf, aspect=70, fraction=0.02, ax=ax, orientation='horizontal', pad=-0.01, extendrect=True)
    # cbar.set_label('Temperature (C)',  fontsize=15, color='white', fontweight='bold')
    # cbar.ax.tick_params(labelcolor='white')
    # for t in cbar.ax.get_xticklabels():
    #     t.set_fontweight('bold')

    # add UND logo
    from PIL import Image
    img = Image.open('utils/images/und-logo.png')
    #                  side-side  up-down  size   size
    imgax = fig.add_axes([0.83, 1.01, 0.06, 0.06], anchor='SE', zorder=3)
    plt.figtext(0.81, 0.995, f'ATMOSPHERIC SCIENCES', ha='left', weight='bold', fontsize=10, color='white')
    imgax.imshow(img)
    imgax.axis('off')

    plt.savefig(f"staged_figures/conus_gfs_forecasts/gfs_850_temp_F{fh}.png", bbox_inches="tight")

    print("    FINISHED 850HPA TEMP MAP")
    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #############################################################################################################################################################################









    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #################################
    # SURFACE TEMPERATURE MAP
    #################################
    fig, ax = build_map()


    # plot mslp
    cs = ax.contour(lons, lats, pres_sfc/100, np.arange(904, 1054, 4), colors='black',
                    linewidths=3.0, linestyles='-',
                    transform=ccrs.PlateCarree(), zorder=11)
    plt.clabel(cs, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
            rightside_up=True, use_clabeltext=True)

    ax.contour(lons, lats, temp_sfc[0, :,:], levels=[0], linewidths=3, linestyles='--', colors='cyan', transform=ccrs.PlateCarree(), zorder=5)  

    contourf = ax.contourf(lons, lats, temp_sfc[0, :,:], np.arange(-50, 51, 1), extent='both',
                    cmap=temp_cmap, alpha=1, transform=ccrs.PlateCarree(), zorder=4)

    # plot  wind barbs
    every = 10
    barbs = ax.barbs(lons.values[0::every], lats.values[ 0::every],
                    uwnd_sfc[0, 0::every, 0::every], vwnd_sfc[0, 0::every, 0::every],
                    length=6.5, alpha=0.7, transform=ccrs.PlateCarree(), zorder=11)

    # plot title, add one to the left with model name and data names, add another to the right with time info
    plt.figtext(0.08, 1.03, f'     {run_date[11:-16]}z GFS Forecast | {valid_date_str}', weight='bold', ha='left', fontsize=20, color='white')
    plt.figtext(0.08, 1.00, f'     Surface MSLP (hPa), Temperature (C), Wind (kts)', ha='left', fontsize=18, color='white')
    plt.figtext(0.915, 1.04, f' ', ha='left', fontsize=20)
    plt.figtext(0.915, -0.01, f' ', ha='left', fontsize=20)
    cax = fig.add_axes([0.91, 0.024, 0.01, 0.95])
    cbar = fig.colorbar(contourf, cax=cax, orientation='vertical', ticks=np.arange(-50, 51, 5), extendrect=True)
    cax.text(3, 0.5, 'Temperature (C)', ha='left',va='center',rotation=270, color='white',fontsize=12,fontweight='bold',transform=cax.transAxes)
    cbar.ax.tick_params(axis='y', labelcolor='white') 
    for t in cbar.ax.get_yticklabels():
        t.set_fontweight('bold')
        t.set_fontsize(9)
    cbar.ax.set_facecolor('black')
    # # colorbar for filled contour
    # cbar = plt.colorbar(contourf, aspect=70, fraction=0.02, ax=ax, orientation='horizontal', pad=-0.01, extendrect=True)
    # cbar.set_label('Temperature (C)',  fontsize=15, color='white', fontweight='bold')
    # cbar.ax.tick_params(labelcolor='white')
    # for t in cbar.ax.get_xticklabels():
    #     t.set_fontweight('bold')

    # add UND logo
    from PIL import Image
    img = Image.open('utils/images/und-logo.png')
    #                  side-side  up-down  size   size
    imgax = fig.add_axes([0.83, 1.01, 0.06, 0.06], anchor='SE', zorder=3)
    plt.figtext(0.81, 0.995, f'ATMOSPHERIC SCIENCES', ha='left', weight='bold', fontsize=10, color='white')
    imgax.imshow(img)
    imgax.axis('off')

    plt.savefig(f"staged_figures/conus_gfs_forecasts/gfs_sfc_temp_F{fh}.png", bbox_inches="tight")

    print("    FINISHED SFC TEMP MAP")
    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #############################################################################################################################################################################

    
    
    
    
    
    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #################################
    # SURFACE PTYPE MAP
    #################################
    fig, ax = build_map(add_sat=True)

    # plot mslp
    cs = ax.contour(lons, lats, pres_sfc/100, np.arange(904, 1054, 4), colors='black',
                    linewidths=3.0, linestyles='-',
                    transform=ccrs.PlateCarree(), zorder=11)
    plt.clabel(cs, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
            rightside_up=True, use_clabeltext=True)
    
    ptype = np.zeros_like(rn_sfc)
    ptype = np.where(sn_sfc >= 0.01, 4, ptype)
    ptype = np.where(ip_sfc >= 0.01, 3, ptype)
    ptype = np.where(zr_sfc >= 0.01, 2, ptype)
    ptype = np.where(rn_sfc >= 0.01, 1, ptype)

    # prepare masking, cmaps, levels, norm
    levels = np.arange(-5, 50, 1)
    norm = BoundaryNorm(levels, 256)
    ptype_data_2d = {
        'Rain': (1, np.where(ptype == 1, reft_sfc, np.nan), 'Greens'),
        'FrzRain': (2, np.where(ptype == 2, reft_sfc, np.nan), 'RdPu'),
        'Sleet': (3, np.where(ptype == 3, reft_sfc, np.nan), 'YlOrRd'),
        'Snow': (4, np.where(ptype == 4, reft_sfc, np.nan), 'Blues'),
    }
        
    # store contourf mappable objs in a dict for cbars
    cf_mappables = {} 
    for name, (value, mask, cmap) in ptype_data_2d.items():
        cf = ax.contourf(lons, lats, mask, levels=levels,cmap=cmap,norm=norm,
            transform=ccrs.PlateCarree(),extend='max', zorder=value)
        cf_mappables[name] = cf
    

    cbar_x_start = 0.91 
    cbar_width = 0.01 
    cbar_height = 0.22
    cbar_spacing = 0.02 
    cbar_y_start_top = 0.75

    ptype_order = ['Rain', 'FrzRain', 'Sleet', 'Snow'] 

    for i, name in enumerate(ptype_order):
        cbar_y_start = cbar_y_start_top - (i * (cbar_height + cbar_spacing))
        cax = fig.add_axes([cbar_x_start, cbar_y_start, cbar_width, cbar_height])
        cbar = fig.colorbar(cf_mappables[name],cax=cax,orientation='vertical',ticks=levels[::3],extendrect=True)
        cax.text( 3.0, 0.5, f'{name}',ha='left',va='center',rotation=270,color='white',fontsize=12,fontweight='bold',transform=cax.transAxes)
        cbar.ax.tick_params(axis='y', labelcolor='white') 
        for t in cbar.ax.get_yticklabels():
            t.set_fontweight('bold')
            t.set_fontsize(9)
        cbar.ax.set_facecolor('black')
    

    # plot 500-1000 thickness 
    ghgt_500 = ghgt_iso[np.where(pres_levs == 500)[0][0]]
    ghgt_1000 = ghgt_iso[np.where(pres_levs == 1000)[0][0]]
    thickness_1000_500 = ndimage.gaussian_filter(ghgt_500 - ghgt_1000, sigma=3.0)

    clevs = (np.arange(0, 5400, 60), np.array([5400]), np.arange(5460, 7000, 60))
    colors = ('tab:blue', 'cyan', 'tab:red')
    kw_clabels = {'fontsize': 11, 'inline': True, 'inline_spacing': 5, 'fmt': '%i','rightside_up': True, 'use_clabeltext': True}

    for clevthick, color in zip(clevs, colors):
        if 5400 in clevthick:
            linestyles = 'solid'
            linewidths = 2.0
        else:
            linestyles = 'solid'
            linewidths = 2.0

        cs = ax.contour(lons, lats,thickness_1000_500, levels=clevthick, colors=color,
            linewidths=linewidths, linestyles=linestyles,transform=ccrs.PlateCarree(),zorder=10)
        plt.clabel(cs, **kw_clabels)

    # plot  wind barbs
    every = 10
    barbs = ax.barbs(lons.values[0::every], lats.values[ 0::every],
                    uwnd_sfc[0, 0::every, 0::every], vwnd_sfc[0, 0::every, 0::every],
                    length=6.5, alpha=0.7, transform=ccrs.PlateCarree(), zorder=11)

    # plot title, add one to the left with model name and data names, add another to the right with time info
    plt.figtext(0.08, 1.03, f'     {run_date[11:-16]}z GFS Forecast | {valid_date_str}', weight='bold', ha='left', fontsize=20, color='white')
    plt.figtext(0.08, 1.00, f'     Surface MSLP (hPa), Composite Reflectivity Precip Type (dBZ), 1000-500 hPa Thickness (m), Wind (kts)', ha='left', fontsize=18, color='white')
    plt.figtext(0.915, 1.04, f' ', ha='left', fontsize=20)
    plt.figtext(0.915, -0.01, f' ', ha='left', fontsize=20)

    # add UND logo
    from PIL import Image
    img = Image.open('utils/images/und-logo.png')
    #                  side-side  up-down  size   size
    imgax = fig.add_axes([0.83, 1.01, 0.06, 0.06], anchor='SE', zorder=3)
    plt.figtext(0.81, 0.995, f'ATMOSPHERIC SCIENCES', ha='left', weight='bold', fontsize=10, color='white')
    imgax.imshow(img)
    imgax.axis('off')

    plt.savefig(f"staged_figures/conus_gfs_forecasts/gfs_sfc_ptype_F{fh}.png", bbox_inches="tight")

    print("    FINISHED SFC PTYPE MAP")
    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #############################################################################################################################################################################

    
    
    print(f"    FINISHED F{fh} FIGURES")


elapsed_time = comp_time.time() - st
print(f"############\nSCRIPT FINISHED: time: {comp_time.strftime('%H:%M:%S', comp_time.gmtime(elapsed_time))}\n############")
