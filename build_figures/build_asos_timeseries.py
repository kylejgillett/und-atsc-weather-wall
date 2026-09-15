##########################################################
#          ASOS/AWOS TIMESERIES PLOT SCRIPT
#  (c) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2025
##########################################################

print("############\nSCRIPT RUNNING\n############")

import time as comp_time
st = comp_time.time()

import warnings
warnings.filterwarnings("ignore")

import sys
import os
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.lines as mlines
from matplotlib.patheffects import withStroke

# get script dir
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(script_dir, ".."))

# get the parent dir
if project_root not in sys.path:
    sys.path.append(project_root)

# import modules from sub dirs
from utils.utils import *
from utils.figure import figure_builder
from get_data.get_asos_obs import get_asos_obs
from utils.style import *


TEXT_OUTLINE = [withStroke(linewidth=3, foreground=FIGURE_BG)]



def style_axis(ax, grid=True):
    ax.set_facecolor(PANEL_BG)
    ax.tick_params(axis="x", colors=TEXT_SECONDARY, labelsize=12, pad=8)
    ax.tick_params(axis="y", colors=TEXT_SECONDARY, labelsize=13)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight("bold")
    for spine in ax.spines.values():
        spine.set_color(SPINE_COLOR)
        spine.set_linewidth(1.3)
    ax.margins(x=0.01)

    if grid:
        ax.grid(True, color=GRID_COLOR, linewidth=0.8,alpha=0.35)


def style_legend(legend):
    legend.get_frame().set_facecolor(LEGEND_BG)
    legend.get_frame().set_edgecolor(LEGEND_EDGE)
    legend.get_frame().set_linewidth(1.0)
    legend.get_frame().set_alpha(0.95)
    for text in legend.get_texts():
        text.set_fontweight("bold")
        text.set_color(TEXT_PRIMARY)


def fmt_value(value, fmt=".0f"):
    if value is None or pd.isna(value):
        return "--"
    if isinstance(value, (int, float, np.number)):
        return format(value, fmt)
    return str(value)


def get_annotation_indices(df, interval="24h"):

    times = pd.to_datetime(df["valid"])
    targets = pd.date_range(
        times.iloc[0].ceil(interval),
        times.iloc[-1].floor(interval),
        freq=interval)
    indices = []
    for target in targets:
        idx = np.abs(times - target).argmin()
        indices.append(idx)
    indices.append(len(df) - 1)
    return sorted(set(indices))



# define stations to get data for and set some metadata 
stations = {
    "GFK": ["ASOS", "Grand Forks International Airport, ND", "47.94, -97.18", "843 ft"],
    "FAR": ["ASOS", "Hector International Airport, ND", "46.92, -96.81", "899 ft"],
    "DVL": ["AWOS", "Devils Lake Municipal Airport, ND", "48.12, -98.92", "1457 ft"],
    "TVF": ["AWOS", "Thief River Falls Regional Airport, MN", "48.07, -96.18", "1115 ft"],
    "96D": ["AWOS", "Walhalla Airport, ND", "48.93, -97.90", "953 ft"],
}

hours = 96


# get station data and build the plot
now_utc = datetime.now(timezone.utc)

for station_id, id_num in zip(stations, range(len(stations))):

    comp_time.sleep(10)

    # get data
    df = get_asos_obs(station_id, hours=hours)
    if df is None or len(df) == 0:
        print(f"    WARNING: NO DATA FOR {station_id}")
        continue
    df = df.copy()


    print(f"    {station_id} DATA LOADED.....{df['valid'].iloc[-1]}")
    skip = max(1, int(len(df) / 6))


    # build figure
    fig = plt.figure(figsize=(16, 10), dpi=250, facecolor=FIGURE_BG)

    # define all axes
    ax1 = fig.add_axes([0.055, 0.590, 0.890, 0.235])
    ax2 = fig.add_axes([0.055, 0.355, 0.890, 0.165], sharex=ax1)
    ax3 = fig.add_axes([0.055, 0.120, 0.890, 0.165], sharex=ax1)
    ax2b = ax2.twinx()
    ax3b = ax3.twinx()


    # apply base settings to all axes
    for ax in [ax1, ax2, ax3]:
        style_axis(ax)

    for ax in [ax2b, ax3b]:
        style_axis(ax, grid=False)
        ax.patch.set_visible(False)


    ax1.tick_params(axis="x", labelbottom=False)
    ax2.tick_params(axis="x", labelbottom=False)
    ax2b.tick_params(axis="x", bottom=False, labelbottom=False)
    ax3b.tick_params(axis="x", bottom=False, labelbottom=False)
    locator = mdates.AutoDateLocator(minticks=6, maxticks=9)
    ax3.xaxis.set_major_locator(locator)
    ax3.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d\n%H:%MZ"))


    ###########################
    # TEMPERATURE AXIS
    ###########################
    # max and min
    temp_min = df[["dwpf", "feel"]].min().min() - 5
    temp_max = df["tmpf"].max() + 10

    ax1.set_ylim(temp_min, temp_max)

    ax1.plot(df["valid"], df["tmpf"], color=TEMP_COLOR, linewidth=4.5, linestyle="-", label="TEMP", zorder=4)
    ax1.plot( df["valid"], df["feel"], color=FEEL_COLOR, linewidth=3, linestyle=":", label="FEELS LIKE", zorder=3)
    ax1.plot( df["valid"], df["dwpf"], color=DEWP_COLOR, linewidth=4.5, linestyle="-", label="DEWPOINT", zorder=4)
    ax1.fill_between(df["valid"], df["tmpf"], temp_min, color=TEMP_COLOR, alpha=0.14, interpolate=True)
    ax1.fill_between( df["valid"], df["feel"], temp_min, color=FEEL_COLOR, alpha=0.05, interpolate=True)
    ax1.axhline(32, color=FREEZE_COLOR, linewidth=1.5, linestyle="--", alpha=0.70, zorder=2)


    # plot data annotations
    annotation_idx = get_annotation_indices(df, interval="6h")
    for idx in annotation_idx:
        x = df["valid"].iloc[idx]
        tmp = df["tmpf"].iloc[idx]
        dwp = df["dwpf"].iloc[idx]

        if not pd.isna(tmp):
            ax1.scatter(x, tmp, s=5, color='k', zorder=6)
            ax1.annotate(  f"{tmp:.0f}°F", xy=(x, tmp), xytext=(0, 10), textcoords="offset points", ha="center",
                        va="bottom", fontsize=10, fontweight="bold", color=TEMP_COLOR, path_effects=TEXT_OUTLINE, clip_on=True)
        if not pd.isna(dwp):
            ax1.scatter(x, dwp, s=5, color='k', zorder=6)
            ax1.annotate(f"{dwp:.0f}°F", xy=(x, dwp),  xytext=(0, -10), textcoords="offset points", ha="center", 
                         va="top", fontsize=10, fontweight="bold", color=DEWP_COLOR,  path_effects=TEXT_OUTLINE, clip_on=True)


    leg = ax1.legend(ncol=3, fontsize=10.5, loc="lower left", bbox_to_anchor=(0.0, 1.015), borderaxespad=0, frameon=False, handlelength=2.5, columnspacing=1.8)
    for text in leg.get_texts():
        text.set_fontweight("bold")
        text.set_color(TEXT_PRIMARY)


    ###########################
    # VIS AND RH AXIS
    ###########################
    # axis max/min
    ax2.set_ylim(0, 15)
    ax2b.set_ylim(0, 125)

    vis_line = ax2.plot(df["valid"], df["vsby"], color=VIS_COLOR, linewidth=4.5, linestyle="-", alpha=0.95, label="VISIBILITY")
    rh_line = ax2b.plot(df["valid"], df["relh"], color=RH_COLOR, linewidth=4, linestyle="--", label="RH")

    ax2b.tick_params(axis="y", labelcolor=RH_COLOR)
    lines = [vis_line[0], rh_line[0]]
    labels = [line.get_label() for line in lines]

    leg = ax2.legend(lines, labels, ncol=2, fontsize=10.5, loc="lower left", bbox_to_anchor=(0.0, 1.025),
         borderaxespad=0, frameon=False, handlelength=2.5, columnspacing=1.8)
    for text in leg.get_texts():
        text.set_fontweight("bold")
        text.set_color(TEXT_PRIMARY)


    # align axis grids
    n_ticks = 6
    left_ticks = np.linspace(*ax2.get_ylim(), n_ticks)
    right_ticks = np.linspace(*ax2b.get_ylim(), n_ticks)
    ax2.set_yticks(left_ticks)
    ax2b.set_yticks(right_ticks)
    ax2b.set_yticklabels([f"{tick:.0f}" if i < len(right_ticks) - 1 else " " for i, tick in enumerate(right_ticks)])

    ax2.grid(True, which="major", axis="y", linestyle="--", color=GRID_COLOR, linewidth=0.8, alpha=0.45)
    ax2b.grid(False)


    ###########################
    # WIND & PRESSURE AXIS
    ###########################
    # axis max/min 
    ax3.set_ylim(0, 70)
    ax3b.set_ylim(980, 1055)

    # plot data
    wind_line = ax3.plot(df["valid"], df["sknt"], color=WIND_COLOR, linewidth=4.0, linestyle="-", label="WIND")
    gust_line = ax3.plot(df["valid"], df["gust"], color=GUST_COLOR, linewidth=2.8, linestyle="--", label="GUST")
    pressure_line = ax3b.plot(df["valid"], df["mslp"], color=MSLP_COLOR, linewidth=2.8, linestyle=":", alpha=0.85, label="MSLP")

    ax3b.tick_params(axis="y", labelcolor=MSLP_COLOR)

    # wind vectors
    wind_idx = get_annotation_indices(df, interval="6h")
    ax3.quiver(df["valid"].iloc[wind_idx], df["y_arrow"].iloc[wind_idx], df["u_norm"].iloc[wind_idx], df["v_norm"].iloc[wind_idx], 
               pivot="middle", scale=50, width=0.002, color=TEXT_MUTED, alpha=0.77, zorder=5)


    # set up legend (more complex for twin axes)
    lines = [wind_line[0], gust_line[0], pressure_line[0]]
    labels = [line.get_label() for line in lines]
    leg = ax3.legend(lines, labels, ncol=3, fontsize=10.5, loc="lower left", bbox_to_anchor=(0.0, 1.025), 
                     borderaxespad=0, frameon=False, handlelength=2.5, columnspacing=1.8)

    for text in leg.get_texts():
        text.set_fontweight("bold")
        text.set_color(TEXT_PRIMARY)

    # align axis grids
    left_ticks = np.linspace(*ax3.get_ylim(), n_ticks)
    right_ticks = np.linspace(*ax3b.get_ylim(), n_ticks)
    ax3.set_yticks(left_ticks)
    ax3b.set_yticks(right_ticks)
    ax3b.set_yticklabels([f"{tick:.0f}" if i < len(right_ticks) - 1 else "" for i, tick in enumerate(right_ticks)])
    ax3.grid(True,which="major", axis="y", linestyle="--", color=GRID_COLOR, linewidth=0.8, alpha=0.45)
    ax3b.grid(False)



    ###########################
    # LATEST OB TEXT
    ###########################
    latest_idx = -1
    latest_time = df["valid"].iloc[latest_idx]
    tmp = fmt_value(df["tmpf"].iloc[latest_idx], ".1f")
    dwp = fmt_value(df["dwpf"].iloc[latest_idx],".1f")
    rh = fmt_value(df["relh"].iloc[latest_idx])
    spd = fmt_value(df["sknt"].iloc[latest_idx])
    gst = fmt_value(df["gust"].iloc[latest_idx])
    vis = fmt_value(df["vsby"].iloc[latest_idx], ".1f")
    slp = fmt_value(df["mslp"].iloc[latest_idx])
    p01 = fmt_value(df["p01i"].iloc[latest_idx], ".1f")
    wxc = fmt_value(df["wxcodes"].iloc[latest_idx])
    card_dir = wind_to_dir(df["u"].iloc[latest_idx], df["v"].iloc[latest_idx])

    textstr = (F"Latest Observation:    "
        f"TMP {tmp}°F   •   "
        f"DWP {dwp}°F   •   "
        f"RH {rh}%   •   "
        f"WIND {card_dir} {spd} kt (G{gst})   •   "
        f"WX {wxc}   •   "
        f"VIS {vis} mi   •   "
        f"MSLP {slp} mb   •   "
        f"1HR PCPN {p01} in   ")


    ###########################
    # BUILD FILENAME
    ###########################
    filename = build_filename(
        "staged_figures/asos_timeseries/",
        "asos_timeseries",
        now_utc,
        variant=f"0{id_num}-{station_id}")



    ###########################
    # BUILD FIG
    ###########################
    figure_builder(
        fig,
        ax1,
        title=f"{hours} hr Surface Observations",
        subtitle=(
            f"{station_id} {stations[station_id][0]} • "
            f"{stations[station_id][1]} • "
            f"{stations[station_id][2]} • "
            f"Elev: {stations[station_id][3]}"),
        valid=(f"Valid • {latest_time.strftime('%Y-%m-%d %H:%MZ')}"),
        footer_left=textstr,
        footer_right=" ",
        save_path=filename,
    )


    plt.close(fig)
    print(f"    FINISHED {station_id} TIMESERIES")

elapsed_time = comp_time.time() - st

print("############\n"
    f"SCRIPT FINISHED: time: "
    f"{comp_time.strftime('%H:%M:%S', comp_time.gmtime(elapsed_time))}\n"
    "############")