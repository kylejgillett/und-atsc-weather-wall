##########################################################
#       SINGLE SITE NEXRAD ANALYSIS PLOT SCRIPT
#  (c) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2025
##########################################################


print("############\nSCRIPT RUNNING\n############")
import time as comp_time
st = comp_time.time()

import warnings
warnings.filterwarnings("ignore")

import os
from contextlib import contextmanager, redirect_stderr, redirect_stdout
@contextmanager
def suppress_stdout_stderr():
    """A context manager that redirects stdout and stderr to devnull"""
    with open(os.devnull, 'w') as fnull:
        with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
            yield (err, out)

with suppress_stdout_stderr():
    import pyart


import numpy as np
from datetime import datetime, timedelta
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from metpy.plots import USCOUNTIES
from cartopy.feature import ShapelyFeature
from cartopy.io.shapereader import Reader
from cartopy.io import img_tiles
import matplotlib.pyplot as plt
import sys
import fsspec
from matplotlib import colors
import matplotlib.cm as cm
import matplotlib


# get script dir
script_dir = os.path.dirname(os.path.abspath(__file__))

# get the parent dir 
project_root = os.path.abspath(os.path.join(script_dir, ".."))

if project_root not in sys.path:
    sys.path.append(project_root)

# import modules from sub dirs
from utils.colormaps import *
from get_data.get_metars import get_metar_data
from get_data.get_nexrad_from_aws import get_radar

hour = 1
minutes = [4, 10] #[24, 30, 36, 41, 47, 53, 58]

for minute in minutes:
    try:
        # get nexrad data
        radar_data, radar_scan_string = get_radar(date=datetime(2025, 7, 8, hour, minute))
        got_radar_data = True
    except: 
        got_radar_data = False
        pass



    def build_map(extent=[-122, -73, 21, 56], projection=ccrs.LambertConformal()):

        fig, axs = plt.subplots(1, 2, figsize=(24, 12), subplot_kw={'projection': projection})
        fig.set_facecolor("#4F5F68")

        for ax in axs:
            ax.set_extent(extent)
            ax.set_box_aspect(0.9)
            ax.add_feature(cfeature.STATES, edgecolor='white', alpha=0.2, linestyle='-', linewidth=2.5, zorder=10)
            ax.add_feature(cfeature.LAND, facecolor="#1a2637", alpha=0.3, zorder=0.1)
            ax.add_feature(USCOUNTIES.with_scale('5m'), alpha=0.7, edgecolor='white', linestyle=':', lw=1, zorder=9)
            satellite = img_tiles.GoogleTiles(style='satellite')
            ax.add_image(satellite, 10) 
        plt.tight_layout()
        
        return fig, axs


    #################################
    # BUILD MAP
    #################################
    center_lat = 47.7500 
    center_lon = -97.050
    box_size = 0.155
    fig, axs = build_map(extent=[center_lon-box_size, center_lon+box_size, center_lat-box_size, center_lat+box_size],
                        projection=ccrs.Mercator())


    ################################
    # PLOT RADAR DATA
    ################################  
    gatefilter = pyart.filters.GateFilter(radar_data)
    gatefilter.exclude_below('reflectivity', 10)
    gatefilter.exclude_equal('reflectivity', 2)
    gatefilter.exclude_equal('reflectivity', 1)
    display = pyart.graph.RadarMapDisplay(radar_data)
    rad_display = display.plot_ppi_map(field= 'reflectivity',
                    sweep=0,
                    ax=axs[0],
                    vmin=-32,
                    vmax=95,
                    title_flag = False,
                    colorbar_flag = False,
                    cmap=rs_expertreflect_cmap,
                    resolution='10m',
                    lat_lines=None, 
                    lon_lines=None,
                    add_grid_lines=False,
                    #gatefilter=gatefilter
                    zorder=5,
                    alpha=1)
        

    dealiased_vel = pyart.correct.dealias_region_based(radar_data,vel_field='velocity', nyquist_velocity=None,keep_original=False)
    radar_data.add_field('dealiased_velocity', dealiased_vel)
    storm_u = 11.0   # m/s (eastward)
    storm_v = -11.0    # m/s (northward)

    # Get dealiased velocity
    vel = radar_data.fields['dealiased_velocity']['data']

    # Radar geometry
    az = np.deg2rad(radar_data.azimuth['data'])      # shape (rays,)
    el = np.deg2rad(radar_data.elevation['data'])    # shape (rays,)

    # Expand to match (rays, gates)
    az_2d = np.repeat(az[:, None], vel.shape[1], axis=1)
    el_2d = np.repeat(el[:, None], vel.shape[1], axis=1)

    # Compute storm radial component TOWARD radar
    storm_radial = (
        storm_u * np.sin(az_2d) * np.cos(el_2d) +
        storm_v * np.cos(az_2d) * np.cos(el_2d)
    )

    # Storm-relative velocity (subtract storm motion)
    srv_data = vel - storm_radial

    # Build Py-ART field dictionary
    srv_field = {
        'data': srv_data,
        'units': 'm/s',
        'long_name': 'storm_relative_velocity',
        'standard_name': 'radial_velocity'
    }

    # Add it to radar object
    radar_data.add_field('storm_relative_velocity', srv_field)

    #################################
    # PLOT RADAR DATA
    #################################  
    gatefilter = pyart.filters.GateFilter(radar_data)
    display = pyart.graph.RadarMapDisplay(radar_data)
    rad_display = display.plot_ppi_map(field='storm_relative_velocity',
                    sweep=1,
                    ax=axs[1],
                    vmin=-30,
                    vmax=30,
                    title_flag = False,
                    colorbar_flag = False,
                    cmap="balance",
                    resolution='10m',
                    lat_lines=None, 
                    lon_lines=None,
                    add_grid_lines=False,
                    gatefilter=gatefilter,
                    zorder=5,
                    alpha=1)


    #################################
    # ADD MAP EXTRAS
    #################################
    # plot title, add one to the left with model name and data names, add another to the right with time info
    plt.figtext(0.0, 0.998, f'  KMVX | {radar_scan_string}z - ELEVATION: 0.5°', weight='bold', ha='left', fontsize=36, color='white')
    plt.figtext(0.0, 0.95, f'   BASE REFLECTIVITY', ha='left', fontsize=30, color='white')
    plt.figtext(0.73, 0.95, f' STORM RELATIVE VELOCITY', ha='left', fontsize=30, color='white')
    plt.figtext(0.915, 1.02, f' ', ha='left', fontsize=20)


    # Define the colormap and norm
    norm = colors.Normalize(vmin=-32, vmax=95)
    sm = cm.ScalarMappable(cmap=rs_expertreflect_cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=axs[0],  aspect=40, fraction=0.05, orientation='horizontal', pad=-0.01, extendrect=True)
    cbar.set_label("Equivalent Reflectivity Factor (dBZ)\n", fontsize=24, color='white') 
    cbar.ax.tick_params(color='white', labelcolor='white')
    cbar.ax.tick_params(labelsize=15)

    norm = colors.Normalize(vmin=-30, vmax=30)
    sm = cm.ScalarMappable(cmap="balance", norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=axs[1],  aspect=40, fraction=0.05, orientation='horizontal', pad=-0.01, extendrect=True)
    cbar.set_label("Storm Relative Base Velocity (m/s)\n", fontsize=24, color='white') 
    cbar.ax.tick_params(color='white', labelcolor='white')
    cbar.ax.tick_params(labelsize=15)

    plt.savefig(f"staged_figures/temp/nexrad_analysis__{hour}-{minute}.png", bbox_inches="tight")

    elapsed_time = comp_time.time() - st
    print(f"############\nSCRIPT FINISHED: time: {comp_time.strftime("%H:%M:%S", comp_time.gmtime(elapsed_time))}\n############")