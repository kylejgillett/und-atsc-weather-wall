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
import gc
from matplotlib.colors import BoundaryNorm


# get script dir
script_dir = os.path.dirname(os.path.abspath(__file__))

# get the parent dir 
project_root = os.path.abspath(os.path.join(script_dir, ".."))

if project_root not in sys.path:
    sys.path.append(project_root)

# import modules from sub dirs
from utils.colormaps import *
from utils.utils import *
from get_data.get_gfs_data import gfs_forecast
from utils.map import map_builder
from utils.figure import figure_builder

# parse date information for sat data download
utc_date = datetime.now(timezone.utc)




#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################

for fh, forecast_time, raw_data in gfs_forecast(center_lat=37.86, center_lon=-98.61, box_size=45, forecast_hours=[12, 24, 36, 48, 60]):


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
    valid_date_str = f"+ F{fh}hr • Valid • {valid_day_name} {valid_date[5:7]}/{valid_date[8:10]} {valid_date[-8:-6]}z"
    print(valid_date)

    # BASIC DATA EXTRACTION
    # PRES LEVS
    pres_levs = raw_data['isobaric'][:] / 100

    # DEFINE PRESSURE LEVELS
    plev300 = np.where(pres_levs == 300)[0][0]
    plev500 = np.where(pres_levs == 500)[0][0]
    plev850 = np.where(pres_levs == 850)[0][0]
    plev1000 = np.where(pres_levs == 1000)[0][0]

    # BASIC DATA EXTRACTION
    sigma = 1.5
    ghgt_300 = ndimage.gaussian_filter(raw_data['Geopotential_height_isobaric'][plev300].values, sigma=sigma)
    ghgt_500 = ndimage.gaussian_filter(raw_data['Geopotential_height_isobaric'][plev500].values, sigma=sigma)
    ghgt_850 = ndimage.gaussian_filter(raw_data['Geopotential_height_isobaric'][plev850].values, sigma=sigma)
    ghgt_1000 = ndimage.gaussian_filter(raw_data['Geopotential_height_isobaric'][plev1000].values, sigma=sigma)
    temp_850 = ndimage.gaussian_filter(raw_data['Temperature_isobaric'][plev850].values, sigma=sigma) - 273.15
    uwnd_300 = ndimage.gaussian_filter(raw_data['u-component_of_wind_isobaric'][plev300].values, sigma=sigma) * 1.94384
    vwnd_300 = ndimage.gaussian_filter(raw_data['v-component_of_wind_isobaric'][plev300].values, sigma=sigma) * 1.94384
    uwnd_500 = ndimage.gaussian_filter(raw_data['u-component_of_wind_isobaric'][plev500].values, sigma=sigma) * 1.94384
    vwnd_500 = ndimage.gaussian_filter(raw_data['v-component_of_wind_isobaric'][plev500].values, sigma=sigma) * 1.94384
    uwnd_850 = ndimage.gaussian_filter(raw_data['u-component_of_wind_isobaric'][plev850].values, sigma=sigma) * 1.94384
    vwnd_850 = ndimage.gaussian_filter(raw_data['v-component_of_wind_isobaric'][plev850].values, sigma=sigma) * 1.94384

    pres_sfc = ndimage.gaussian_filter(raw_data.variables['MSLP_Eta_model_reduction_msl'], sigma=sigma)
    temp_sfc = ndimage.gaussian_filter(raw_data.variables['Temperature_height_above_ground'], sigma=sigma) - 273.15
    uwnd_sfc = ndimage.gaussian_filter(raw_data.variables['u-component_of_wind_height_above_ground'], sigma=sigma) * 1.94384
    vwnd_sfc = ndimage.gaussian_filter(raw_data.variables['v-component_of_wind_height_above_ground'], sigma=sigma) * 1.94384
    reft_sfc = ndimage.gaussian_filter(raw_data.variables['Composite_reflectivity_entire_atmosphere'], 0.01)
    rn_sfc = ndimage.gaussian_filter(raw_data.variables['Categorical_Rain_surface'], 0.75)
    sn_sfc = ndimage.gaussian_filter(raw_data.variables['Categorical_Snow_surface'], 0.75)
    zr_sfc = ndimage.gaussian_filter(raw_data.variables['Categorical_Freezing_Rain_surface'], 0.75)
    ip_sfc = ndimage.gaussian_filter(raw_data.variables['Categorical_Ice_Pellets_surface'], 0.75)
    thickness_1000_500 = ndimage.gaussian_filter(ghgt_500 - ghgt_1000,sigma=3.0)


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

    #################################
    # CALCULATE 2PVU SURFACE
    #################################
    pvu_level = 2.0e-6
    pres = raw_data['isobaric'].values
    sort_idx = np.argsort(pres)
    pres = pres[sort_idx]
    pv_iso = raw_data['pv'].values[sort_idx, :, :]
    thta_iso = raw_data['theta'].values[sort_idx, :, :]
    uwnd_pvu = raw_data['u-component_of_wind_isobaric'].values[sort_idx, :, :]
    vwnd_pvu = raw_data['v-component_of_wind_isobaric'].values[sort_idx, :, :]
    pvu_mask = (pres >= 10000) & (pres <= 50000)
    pres = pres[pvu_mask]
    pv_iso = pv_iso[pvu_mask, :, :]
    thta_iso = thta_iso[pvu_mask, :, :]
    uwnd_pvu = uwnd_pvu[pvu_mask, :, :]
    vwnd_pvu = vwnd_pvu[pvu_mask, :, :]
    pres_3d = np.broadcast_to(pres[:, None, None], pv_iso.shape)
    valid_2pvu = ((np.nanmin(pv_iso, axis=0) <= pvu_level) &(np.nanmax(pv_iso, axis=0) >= pvu_level))

    # interpolate onto the 2 PVU surface
    thta_on_2pvu = interpolate_to_isosurface(pv_iso, thta_iso, pvu_level, bottom_up_search=False)
    p_on_2pvu = interpolate_to_isosurface(pv_iso, pres_3d, pvu_level, bottom_up_search=False) / 100.0
    u_on_2pvu = interpolate_to_isosurface(pv_iso, uwnd_pvu, pvu_level, bottom_up_search=False) * 1.94384
    v_on_2pvu = interpolate_to_isosurface(pv_iso, vwnd_pvu, pvu_level, bottom_up_search=False) * 1.94384

    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #############################################################################################################################################################################









    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #################################
    # BUILD 300 HPA MAP
    #################################
    fig, ax = map_builder(extent=[-119, -74, 23.5, 53.5], terrain=True, state_color='navy', border_color='navy', border_factor=1.5)

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

    composite_filename = build_filename("staged_figures/conus_gfs_forecasts/", f"gfs_300a", utc_date, variant=fh)

    figure_builder(fig, ax,
        title=f"GFS Forecast • 300hPa",
        subtitle=f"Heights (m), Wind (kt)",
        valid=valid_date_str,
        mappable=contourf,
        cbar_title="Wind Speed",
        cbar_units="kts",
        cbar_ticks=np.arange(50, 160, 5)[::1],
        footer_left=f"GFS 0.25° • INIT {run_date[0:10]} {run_date[11:-13]}z",
        save_path=composite_filename)

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
    fig, ax = map_builder(extent=[-119, -74, 23.5, 53.5], terrain=True, state_color='k', border_color='k', border_factor=1.5)

    n_reps = 20
    p_on_2pvu = mpcalc.smooth_n_point(p_on_2pvu, 9, n_reps)

    # plot pressure on the 2 PVU surface
    contour = ax.contour(lons, lats, p_on_2pvu, np.arange(100, 751, 50),
                    colors='black', linewidths=3.0, linestyles='-',
                    transform=ccrs.PlateCarree(), zorder=11)
    plt.clabel(contour, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
            rightside_up=True, use_clabeltext=True)

    # plot 300hpa pv fill
    contourf = ax.contourf(raw_data['longitude'], raw_data['latitude'], thta_on_2pvu, np.arange(270, 410, 5), cmap=pv_cmap,
                    transform=ccrs.PlateCarree(), extend='both')

    every = 10
    barbs = ax.barbs(lons.values[0::every], lats.values[0::every],
                    u_on_2pvu[0::every, 0::every], v_on_2pvu[0::every, 0::every],
                    length=6.5, alpha=0.7, transform=ccrs.PlateCarree(), zorder=11)

    composite_filename = build_filename("staged_figures/conus_gfs_forecasts/", f"gfs_300b", utc_date, variant=fh)

    figure_builder(fig, ax,
        title=f"GFS Forecast • 2PVU Surface",
        subtitle=f"Pressure (hPa), Potential Temperature (K), Wind (kt)",
        valid=valid_date_str,
        mappable=contourf,
        cbar_title="2PVU Potential Temperature",
        cbar_units=r"K",
        cbar_ticks=np.arange(270, 410, 10),
        footer_left=f"GFS 0.25° • INIT {run_date[0:10]} {run_date[11:-13]}z",
        save_path=composite_filename)
    
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
    fig, ax = map_builder(extent=[-119, -74, 23.5, 53.5], terrain=True, state_color='navy', border_color='navy', border_factor=1.5)

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

    composite_filename = build_filename("staged_figures/conus_gfs_forecasts/", f"gfs_500a", utc_date, variant=fh)

    figure_builder(fig, ax,
        title=f"GFS Forecast • 500hPa",
        subtitle=f"Heights (m), Wind (kt)",
        valid=valid_date_str,
        mappable=contourf,
        cbar_title="Wind Speed",
        cbar_units="kts",
        cbar_ticks=np.arange(30, 140, 5),
        footer_left=f"GFS 0.25° • INIT {run_date[0:10]} {run_date[11:-13]}z",
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
    fig, ax = map_builder(extent=[-119, -74, 23.5, 53.5], terrain=True, state_color='navy', border_color='navy', border_factor=1.5)

    n_reps = 50
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
    every = 10
    barbs = ax.barbs(lons.values[0::every], lats.values[ 0::every],
                    uwnd_500[0::every, 0::every], vwnd_500[0::every, 0::every],
                    length=6.5, alpha=0.7, transform=ccrs.PlateCarree(), zorder=12)


    composite_filename = build_filename("staged_figures/conus_gfs_forecasts/", f"gfs_500b", utc_date, variant=fh)

    figure_builder(fig, ax,
        title=f"GFS Forecast • 500hPa",
        subtitle="Heights (m), Rel. Vorticity "+r"($\mathrm{s}^{-1} \times 10^{5}$)" +", Wind (kt)",
        valid=valid_date_str,
        mappable=contourf,
        cbar_title="Relative Vorticity",
        cbar_units=r"$\mathrm{s}^{-1} \times 10^{5}$",
        cbar_ticks=np.arange(-30, 52, 2),
        footer_left=f"GFS 0.25° • INIT {run_date[0:10]} {run_date[11:-13]}z",
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
    # #absvort_adv = mpcalc.advection(avor_500, uwnd_500, vwnd_500, dx=dx, dy=dy) *1e9

    # # plot 500hpa heights
    # contour = ax.contour(lons, lats, ghgt_500, np.arange(3000, 7000, 60),
    #                 colors='black', linewidths=3.0, linestyles='-',
    #                 transform=ccrs.PlateCarree(), zorder=11)
    # plt.clabel(contour, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
    #         rightside_up=True, use_clabeltext=True)

    # # plot 500hpa abs vort adv fill
    # vortadv_cf = ax.contourf(lons, lats, relvort_adv, np.arange(-40, 42, 2),               #np.arange(-6*12**-7, 6*12**-7, 1*10**-9),
    #                             extend='both', cmap='bwr', zorder=5, alpha=1, transform=ccrs.PlateCarree())

    # # plot 500hpa wind barbs
    # every = 10
    # barbs = ax.barbs(lons.values[0::every], lats.values[ 0::every],
    #                 uwnd_500[0::every, 0::every], vwnd_500[0::every, 0::every],
    #                 length=6.5, alpha=0.7, transform=ccrs.PlateCarree(), zorder=12)

    # # plot title, add one to the left with model name and data names, add another to the right with time info
    # plt.figtext(0.08, 1.03, f'     {run_date[11:-16]}z GFS Forecast | {valid_date_str}', weight='bold', ha='left', fontsize=20, color='white')
    # plt.figtext(0.08, 1.00, f'     500 hPa Heights (m), Rel. Vorticity Adv. (sec⁻²•10⁹), Wind (kt)', ha='left', fontsize=18, color='white')
    # plt.figtext(0.915, 1.04, f' ', ha='left', fontsize=20)
    # plt.figtext(0.915, -0.01, f' ', ha='left', fontsize=20)
    # cax = fig.add_axes([0.91, 0.024, 0.01, 0.95])
    # cbar = fig.colorbar(vortadv_cf, cax=cax, orientation='vertical', ticks=np.arange(-40, 42, 2), extendrect=True)
    # cax.text(3, 0.5, 'Relative Vorticity Advection (sec⁻²•10⁹)', ha='left',va='center',rotation=270, color='white',fontsize=12,fontweight='bold',transform=cax.transAxes)
    # cbar.ax.tick_params(axis='y', labelcolor='white') 
    # for t in cbar.ax.get_yticklabels():
    #     t.set_fontweight('bold')
    #     t.set_fontsize(9)
    # cbar.ax.set_facecolor('black')
    # # # colorbar for filled contour
    # # cbar = plt.colorbar(vortadv_cf, aspect=70, fraction=0.02, ax=ax, orientation='horizontal', pad=-0.01, extendrect=True)
    # # cbar.set_label('Relative Vorticity Advection (sec⁻²•10⁹)',  fontsize=15, color='white', fontweight='bold')
    # # cbar.ax.tick_params(labelcolor='white')
    # # for t in cbar.ax.get_xticklabels():
    # #     t.set_fontweight('bold')

    # # add UND logo
    # from PIL import Image
    # img = Image.open('utils/images/und-logo.png')
    # #                  side-side  up-down  size   size
    # imgax = fig.add_axes([0.83, 1.01, 0.06, 0.06], anchor='SE', zorder=3)
    # plt.figtext(0.81, 0.995, f'ATMOSPHERIC SCIENCES', ha='left', weight='bold', fontsize=10, color='white')
    # imgax.imshow(img)
    # imgax.axis('off')

    # plt.savefig(f"staged_figures/conus_gfs_forecasts/gfs_500_relvortadv_F{fh}.png", bbox_inches="tight")
    #plt.close(fig)
    # print("    FINISHED 500HPA REL VORT ADV MAP")
    # #############################################################################################################################################################################
    # #############################################################################################################################################################################
    # #############################################################################################################################################################################







    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #################################
    # BUILD 850HPA TADV MAP
    #################################
    fig, ax = map_builder(extent=[-119, -74, 23.5, 53.5], terrain=True, state_color='navy', border_color='navy', border_factor=1.5)

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

    composite_filename = build_filename("staged_figures/conus_gfs_forecasts/", f"gfs_850b", utc_date, variant=fh)

    figure_builder(fig, ax,
        title=f"GFS Forecast • 850hPa",
        subtitle=f'Heights (m), 3hr Temperature Adv (C/3hr), Frontogenesis (>2'+u'\xb0'+'C / 100km / 3hr), Wind (kt)',
        valid=valid_date_str,
        mappable=tadv_contourf,
        cbar_title="Temperature Advection",
        cbar_units='('+u'\xb0'+'C / 3hr)',
        cbar_ticks=np.arange(-7,7.25,1),
        footer_left=f"GFS 0.25° • INIT {run_date[0:10]} {run_date[11:-13]}z",
        save_path=composite_filename)
    
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
    fig, ax = map_builder(extent=[-119, -74, 23.5, 53.5], terrain=True, state_color='navy', border_color='navy', border_factor=1.5)

    # plot 850 hpa heights
    contour = ax.contour(lons, lats, ghgt_850, np.arange(0, 1700, 30),
                    colors='black', linewidths=3.0, linestyles='-',
                    transform=ccrs.PlateCarree(), zorder=11)
    plt.clabel(contour, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
            rightside_up=True, use_clabeltext=True)
    ax.contour(lons, lats, temp_850, levels=[0], linewidths=3, linestyles='--', colors='cyan', transform=ccrs.PlateCarree(), zorder=5)  
    contourf = ax.contourf(lons, lats, temp_850, np.arange(-40, 42, 1), extent='both',
                    cmap=temp_cmap, alpha=1, transform=ccrs.PlateCarree(), zorder=4)

    # plot 850hpa wind barbs
    every = 10
    barbs = ax.barbs(lons.values[0::every], lats.values[ 0::every],
                    uwnd_850[0::every, 0::every], vwnd_850[0::every, 0::every],
                    length=6.5, alpha=0.7, transform=ccrs.PlateCarree(), zorder=11)

    composite_filename = build_filename("staged_figures/conus_gfs_forecasts/", f"gfs_850a", utc_date, variant=fh)

    figure_builder(fig, ax,
        title=f"GFS Forecast • 850hPa",
        subtitle=f'Heights (m), Temperature (C), Wind (kts)',
        valid=valid_date_str,
        mappable=contourf,
        cbar_title="Temperature",
        cbar_units='°C',
        cbar_ticks=np.arange(-40, 42, 5),
        footer_left=f"GFS 0.25° • INIT {run_date[0:10]} {run_date[11:-13]}z",
        save_path=composite_filename)

    
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
    fig, ax = map_builder(extent=[-119, -74, 23.5, 53.5], terrain=True, state_color='navy', border_color='navy', border_factor=1.5)

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

    composite_filename = build_filename("staged_figures/conus_gfs_forecasts/", f"gfs_000a", utc_date, variant=fh)

    figure_builder(fig, ax,
        title=f"GFS Forecast • Surface",
        subtitle=f'MSLP (hPa), 2m Temperature (C), 10m Wind (kts)',
        valid=valid_date_str,
        mappable=contourf,
        cbar_title="Temperature",
        cbar_units='°C',
        cbar_ticks=np.arange(-50, 51, 5),
        footer_left=f"GFS 0.25° • INIT {run_date[0:10]} {run_date[11:-13]}z",
        save_path=composite_filename)

    
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
    fig, ax = map_builder(extent=[-119, -74, 23.5, 53.5], terrain=True, state_color='navy', border_color='navy', border_factor=1.5)

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
    

    cbar_x_start = 0.92 
    cbar_width = 0.02 
    cbar_height = 0.16
    cbar_spacing = 0.01
    cbar_y_start_top = 0.65

    ptype_order = ['Rain', 'FrzRain', 'Sleet', 'Snow'] 

    for i, name in enumerate(ptype_order):
        cbar_y_start = cbar_y_start_top - (i * (cbar_height + cbar_spacing))
        cax = fig.add_axes([cbar_x_start, cbar_y_start, cbar_width, cbar_height])
        cbar = fig.colorbar(cf_mappables[name],cax=cax,orientation='vertical',ticks=levels[::5],extendrect=True)
        cax.text(2.0, 0.5, f'{name}',ha='left',va='center',rotation=270,color='white',fontsize=12,fontweight='bold',transform=cax.transAxes)
        cbar.ax.tick_params(axis='y', labelcolor='white') 
        for t in cbar.ax.get_yticklabels():
            t.set_fontweight('bold')
            t.set_fontsize(9)
        cbar.ax.set_facecolor('black')
    

    # plot 500-1000 thickness 
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

    composite_filename = build_filename("staged_figures/conus_gfs_forecasts/", f"gfs_000b", utc_date, variant=fh)

    figure_builder(fig, ax,
        title=f"GFS Forecast • Surface",
        subtitle=f'MSLP (hPa), Composite Reflectivity Precip Type (dBZ), 1000-500 hPa Thickness (m), 10m Wind (kts)',
        valid=valid_date_str,
        mappable=None,
        cbar_title="Temperature",
        cbar_units='°C',
        cbar_ticks=np.arange(-50, 51, 5),
        footer_left=f"GFS 0.25° • INIT {run_date[0:10]} {run_date[11:-13]}z",
        save_path=composite_filename)
    
    print("    FINISHED SFC PTYPE MAP")
    #############################################################################################################################################################################
    #############################################################################################################################################################################
    #############################################################################################################################################################################

    
    
    print(f"    FINISHED F{fh} FIGURES")
    plt.close('all')
    gc.collect()

elapsed_time = comp_time.time() - st
print(f"############\nSCRIPT FINISHED: time: {comp_time.strftime('%H:%M:%S', comp_time.gmtime(elapsed_time))}\n############")
