##########################################################
#           ASOS/AWOS TIMESERIES PLOT SCRIPT
#  (c) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2025
##########################################################

print("############\nSCRIPT RUNNING\n############")
import time as comp_time
st = comp_time.time()

import warnings
warnings.filterwarnings("ignore")

import sys 
import os
import pandas as pd
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.lines as mlines
from matplotlib.patheffects import withStroke
import scipy.ndimage as ndimage


# get script dir
script_dir = os.path.dirname(os.path.abspath(__file__))

# get the parent dir 
project_root = os.path.abspath(os.path.join(script_dir, ".."))
if project_root not in sys.path:
    sys.path.append(project_root)

# import modules from sub dirs
from utils.utils import *
from get_data.get_asos_obs import get_asos_obs

# define stations to get data for and set some metadata 
stations = {
    "GFK": ["ASOS", "Grand Forks International Airport, ND", "47.94, -97.18", "843.0 ft"],
    "FAR": ["ASOS", "Hector International Airport, ND", "46.92, -96.81", "899ft ft"],
    "DVL": ["AWOS", "Devils Lake Municipal Airport, ND", "48.12, -98.92", "1457 ft"],
    "TVF": ["AWOS", "Thief River Falls Regional Airport, MN", "48.07, -96.18", "1115 ft"],
    "96D": ["AWOS", "Walhalla Airport, ND", "48.93, -97.90", "953 ft"]
            }





# get station data and build the plot
for station_id in stations:

    comp_time.sleep(10)

    hours = 96

    # get data
    df = get_asos_obs(station_id, hours=hours)

    # define a text outline stroke and a skip
    TEXT_OUTLINE = [withStroke(linewidth=3, foreground=(0, 0, 0, 0.3))]
    skip = int(len(df['drct'])/6)

    # build figure
    fig = plt.figure(figsize=(20, 12))
    fig.set_facecolor('white')
    fig.patch.set_edgecolor('#009946')
    fig.patch.set_linewidth(20) 
    

    # define all axes
    ax1 = plt.axes([0, 0.62, 0.91, 0.35])
    ax2 = plt.axes([0, 0.32, 0.91, 0.20])
    ax3 = plt.axes([0, 0.02, 0.91, 0.20])
    ax_text = plt.axes([0, -0.15, 0.91, 0.10]) 
    ax2b = ax2.twinx()
    ax3b = ax3.twinx()

    # apply base settings to all axes
    for ax in [ax1, ax2, ax3, ax2b, ax3b]:
        ax.set_facecolor('white')
        ax.tick_params(axis='x', labelsize=16, labelcolor='black', pad=10)
        ax.tick_params(axis='y', labelsize=20, labelcolor='black')
        for xlabel, ylabel in zip(ax.get_xticklabels(), ax.get_yticklabels()):
            xlabel.set_fontweight('bold')
            ylabel.set_fontweight('bold')
        ax.margins(0.02)
        ax.xaxis.set_major_formatter((mdates.DateFormatter(' %m-%d'+'\n '+'%H:%M'+'Z')))
        #ax.grid(which='major', axis='y', linestyle=(0, (5, 10)), color='black', linewidth=0.5)
        ax.grid(True, color='black', alpha=0.5)
        for spine in ax.spines.values():
            spine.set_linewidth(3)    
            spine.set_color("#000000") 




    # apply smoothing filter to make prettier
    # df['tmpf'] = ndimage.gaussian_filter1d(df['tmpf'],0.2)
    # df['feel'] = ndimage.gaussian_filter1d(df['feel'],0.2)
    # df['dwpf'] = ndimage.gaussian_filter1d(df['dwpf'],0.2)
    # df['vsby'] = ndimage.gaussian_filter1d(df['vsby'],0.2)
    # df['relh'] = ndimage.gaussian_filter1d(df['relh'],0.2)
    # df['dwpf'] = ndimage.gaussian_filter1d(df['sknt'],0.2)
    # df['vsby'] = ndimage.gaussian_filter1d(df['gust'],0.2)
    # df['relh'] = ndimage.gaussian_filter1d(df['mslp'],0.2)

    ###########################
    # TEMPERATURE AXIS
    ###########################
    # max and min
    ax1.set_ylim(df[['dwpf', 'feel']].min().min()-5, df['tmpf'].max()+12)

    # plot data
    ax1.plot(df['valid'], df['tmpf'], linestyle='-', zorder=1.5, color='orangered', linewidth=5, label='TMPF')
    ax1.plot(df['valid'], df['feel'], linestyle=':', zorder=1, color='blue', linewidth=3, label='FEEL')
    ax1.plot(df['valid'], df['dwpf'], linestyle='-', zorder=1.2, color='co', linewidth=5, label='DWPF')

    ax1.fill_between(df['valid'], df['tmpf'], -30,
                    color='orangered', alpha=0.3, interpolate=True)
    ax1.fill_between(df['valid'], df['dwpf'], -30,
                    color='lightblue', alpha=0.3, interpolate=True)

    # plot data annotations
    for x, y in zip(df['valid'][::-skip], df['tmpf'][::-skip]):
        ax1.text(x, y+3, f"{int(y)}°F", weight='bold', fontsize=15, color='orangered', 
                alpha=1, path_effects=TEXT_OUTLINE, ha='center', clip_on=True)
        ax1.add_line(mlines.Line2D([x, x], [y, y+3], alpha=0.3, color='black'))
    for x, y in zip(df['valid'][::-skip], df['dwpf'][::-skip]):
        ax1.text(x, y-3, f"{int(y)}°F", weight='bold', fontsize=15, color='lightblue', 
                alpha=1, path_effects=TEXT_OUTLINE, ha='center', clip_on=True)
        ax1.add_line(mlines.Line2D([x, x], [y, y-3], alpha=0.3, color='black'))

    # plot hline at 32F
    ax1.axhline(y=32, color='blue', alpha=0.6, zorder=0.2)

    # add legend
    leg = ax1.legend(ncol=3,fontsize=15, facecolor='gainsboro',edgecolor='black',labelcolor='black', loc='upper left')
    for text in leg.get_texts():
        text.set_fontweight('bold')
        text.set_color('black')



    ###########################
    # VIS AND RH AXIS
    ###########################
    # axis max/min
    ax2.set_ylim(0, 15)
    ax2b.set_ylim(0, 125)

    # plot data
    ax2a_line = ax2.plot(df['valid'], df['vsby'], linestyle='-', color='black', linewidth=5, alpha=0.7, label="VIS")
    ax2b_line = ax2b.plot(df['valid'], df['relh'], linestyle='--', color='mediumslateblue', linewidth=5, label="RELH")

    # set up legend (more complex for twin axes)
    ax2_lines = [ax2a_line, ax2b_line]
    ax2_labels = [l[0].get_label() for l in ax2_lines]
    leg = ax2.legend([l[0] for l in ax2_lines], ax2_labels, ncol=2, fontsize=15, facecolor='gainsboro',edgecolor='black',
                    labelcolor='black', loc='upper left')
    for text in leg.get_texts():
        text.set_fontweight('bold')
        text.set_color('black')

    # align axis grids
    n_ticks = 6
    left_ticks = np.linspace(*ax2.get_ylim(), n_ticks)
    right_ticks = np.linspace(*ax2b.get_ylim(), n_ticks)
    ax2.set_yticks(left_ticks)
    ax2b.set_yticks(right_ticks)
    ax2b.set_yticklabels([f"{t:.0f}" if i < len(right_ticks) - 1 else "" for i, t in enumerate(right_ticks)])
    ax2.grid(True, which='major', axis='y', linestyle='--', alpha=0.5)
    ax2b.grid(False)


    ###########################
    # WIND AXIS
    ###########################
    # axis max/min 
    ax3.set_ylim(0, 60)
    ax3b.set_ylim(980, 1055)

    # plot data
    ax3a_line1 = ax3.plot(df['valid'], df['sknt'], linestyle='-', color='lightblue', linewidth=5, label='SPD KT')
    ax3a_line2 = ax3.plot(df['valid'], df['gust'], linestyle='--', color='blue', linewidth=3, label="GST KT")
    ax3b_line = ax3b.plot(df['valid'], df['mslp'], linestyle=':', color='gray', linewidth=5, label="MSLP (mb)")

    # plot wind dir vectors 
    qvr = ax3.quiver(df['valid'][::-skip], df['y_arrow'][::-skip], df['u_norm'][::-skip], df['v_norm'][::-skip], pivot='middle',
            scale=50, width=0.002, color='black', alpha=0.5)

    # set up legend (more complex for twin axes)
    ax3_lines = [ax3a_line1, ax3a_line2, ax3b_line]
    ax3_labels = [l[0].get_label() for l in ax3_lines]
    leg = ax3.legend([l[0] for l in ax3_lines], ax3_labels, ncol=3, fontsize=15, facecolor='gainsboro',edgecolor='black',
                    labelcolor='black', loc='upper left')
    for text in leg.get_texts():
        text.set_fontweight('bold')
        text.set_color('black')

    # align axis grids
    n_ticks = 6
    left_ticks = np.linspace(*ax3.get_ylim(), n_ticks)
    right_ticks = np.linspace(*ax3b.get_ylim(), n_ticks)
    ax3.set_yticks(left_ticks)
    ax3b.set_yticks(right_ticks)
    ax3b.set_yticklabels([f"{t:.0f}" if i < len(right_ticks) - 1 else "" for i, t in enumerate(right_ticks)])
    ax3.grid(True, which='major', axis='y', linestyle='--', alpha=0.5)
    ax3b.grid(False)




    ###########################
    # LATEST INFO
    ###########################

    # turn off frame and ticks
    ax_text.axis("off")

    # get latest values
    latest_idx = -1 
    time_str = df['valid'].iloc[latest_idx].strftime('%Y-%m-%d %H:%M')
    tmp = df['tmpf'].iloc[latest_idx]
    dwp = df['dwpf'].iloc[latest_idx]
    rh = df['relh'].iloc[latest_idx]
    spd = df['sknt'].iloc[latest_idx]
    gst = df['gust'].iloc[latest_idx]
    vis = df['vsby'].iloc[latest_idx]
    slp = df['mslp'].iloc[latest_idx]
    p01 = df['p01i'].iloc[latest_idx]
    wxc = df['wxcodes'].iloc[latest_idx]
    card_dir = wind_to_dir(df['u'].iloc[-1], df['v'].iloc[-1])

    # if vals are nan, set to '--' to look nicer
    if math.isnan(gst):
        gst = '--'
    else:
        gst = f'{gst:.0f}'

    if isinstance(wxc, (int, float)):
        if math.isnan(wxc):
            wxc = '--'
        else:
            wxc = f'{wxc}'
    if math.isnan(slp):
        slp = '--'
    else:
        slp = f'{slp:.0f}'
    
    # define the text string to plot
    textstr = (
        f"Latest Observation ({time_str}z)\n"
        f"TMP: {tmp:.1f}°F,   DWP: {dwp:.1f}°F,   RH: {rh:.0f}%   |   WIND: {card_dir} {spd:.0f} kt (G{gst})   |   "
        f"PRESENT WX: {wxc},   VIS: {vis:.1f} mi,   MSLP: {slp} mb,    1HR ACCM: {p01:.1f} in"
    )

    # plot the text string and apply some formatting 
    ax_text.text(0.5, 0.5, textstr,
                ha='center', va='center',
                fontsize=16, fontweight='bold',
                color='black', family='monospace',
                bbox=dict(facecolor='gainsboro', edgecolor='black', boxstyle='round,pad=0.5'))


    ###########################
    # TITLES AND EXTRAS 
    ###########################
    plt.figtext(0.0, 1.01, f'  {hours}hr Surface Observations', weight='bold', ha='left', fontsize=32, color='black')
    plt.figtext(0.0, 0.98, f'   {station_id} {stations[station_id][0]} - {stations[station_id][1]}  |  {stations[station_id][2]} | Elev:  {stations[station_id][3]}', ha='left', fontsize=23, color='black')
    plt.figtext(-0.03, 1.038, f' ', ha='left', fontsize=20)
    img = Image.open('utils/images/und-logo.png')
    #                  side-side  up-down  size   size
    imgax = fig.add_axes([0.831, 0.995, 0.06, 0.06], anchor='SE', zorder=3)
    plt.figtext(0.81, 0.979, f'ATMOSPHERIC SCIENCES', ha='left', weight='bold', fontsize=10, color='black')
    imgax.imshow(img)
    imgax.axis('off')


    plt.savefig(f"staged_figures/asos_timeseries/{station_id}_timeseries.png", bbox_inches="tight")

    print(f"    FINISHED {station_id} timeseries")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################



elapsed_time = comp_time.time() - st
print(f"############\nSCRIPT FINISHED: time: {comp_time.strftime("%H:%M:%S", comp_time.gmtime(elapsed_time))}\n############")