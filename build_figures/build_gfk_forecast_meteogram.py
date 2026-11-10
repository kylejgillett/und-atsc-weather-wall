##########################################################
#      GFK 7-DAY NWS & GEFS FORECAST METEOGRAM
#  (c) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2026
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
from matplotlib.colors import to_rgba
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import FuncFormatter
from matplotlib.patheffects import withStroke

# get script dir
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(script_dir, ".."))

# get the parent dir
if project_root not in sys.path:
    sys.path.append(project_root)

# import modules from sub dirs
from utils.utils import *
from get_data.get_gfk_forecast import get_nws_forecast, get_gefs_plumes
from utils.style import *


TEXT_OUTLINE = [withStroke(linewidth=2.2, foreground=PANEL_BG)]
LOCAL_TZ = "America/Chicago"

# axes span (figure fraction); legends sit in the column to the right
AX_LEFT = 0.050
AX_WIDTH = 0.800
AX_RIGHT = AX_LEFT + AX_WIDTH

PRECIP_TYPE_COLORS = {
    "Rain": RAIN_COLOR,
    "Snow": SNOW_COLOR,
    "Wintry Mix": MIX_COLOR,
    "Freezing Rain": FZRA_COLOR,
    "Thunderstorms": TSTM_COLOR,
    "": POP_COLOR,
}

PRECIP_TYPE_LABELS = {
    "Rain": "NWS RAIN",
    "Snow": "NWS SNOW",
    "Wintry Mix": "NWS MIX/SLEET",
    "Freezing Rain": "NWS FZ RAIN",
    "Thunderstorms": "NWS T-STORMS",
}



def style_axis(ax, grid=True):
    ax.set_facecolor(PANEL_BG)
    ax.tick_params(axis="x", colors=TEXT_SECONDARY, labelsize=8, pad=4)
    ax.tick_params(axis="y", colors=TEXT_SECONDARY, labelsize=8.5, pad=3)
    for spine in ax.spines.values():
        spine.set_color(SPINE_COLOR)
        spine.set_linewidth(1.0)

    if grid:
        ax.grid(True, color=GRID_COLOR, linewidth=0.6, alpha=0.35)


def bold_ticks(ax):
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight("bold")


def side_legend(ax, handles, loc="upper left", anchor_y=1.02):
    leg = ax.legend(handles=handles, ncol=1, fontsize=7.2, loc=loc, bbox_to_anchor=(1.012, anchor_y),
                    borderaxespad=0, frameon=False, handlelength=2.0, handletextpad=0.6, labelspacing=0.45)
    for text in leg.get_texts():
        text.set_fontweight("bold")
        text.set_color(TEXT_PRIMARY)

    return leg


# label 00Z ticks with the date and 12Z ticks with the hour only
def utc_tick_label(x, pos):
    tick = mdates.num2date(x, tz=timezone.utc)
    if tick.hour == 0:
        return f"00Z\n{tick.strftime('%a %m/%d').upper()}"
    return f"{tick:%H}Z"


def format_utc_label(timestamp):
    return timestamp.strftime("%HZ %a %d %b %Y").upper()



# frame sized for the wide, short meteogram slot (same look as utils.figure.figure_builder,
# with positions and font sizes scaled for a 14 x 5.2 figure); header and footer
# are aligned with the left and right edges of the plot axes
def meteogram_frame(fig, *, title, subtitle, valid, footer_left, footer_right, save_path=None):

    # outer frame
    border = Rectangle((0.002, 0.002), 0.996, 0.996, transform=fig.transFigure,
        fill=False, edgecolor=UND_GREEN, linewidth=4, zorder=10000, clip_on=False)
    fig.add_artist(border)

    # header text
    fig.text(AX_LEFT, 0.968, title, ha="left", va="top", fontsize=16, fontweight="bold", color=TEXT_PRIMARY)
    fig.text(AX_LEFT, 0.902, subtitle, ha="left", va="top", fontsize=10, color=TEXT_SECONDARY)
    fig.text(AX_LEFT, 0.858, valid, ha="left", va="top", fontsize=8, fontweight="bold", color=TEXT_VALID)
    fig.add_artist(plt.Line2D([AX_LEFT, AX_LEFT + 0.073], [0.818, 0.818],
            transform=fig.transFigure, color=UND_GREEN, linewidth=2.2, solid_capstyle="round"))

    # und logo and text, right edge on the axes' right edge
    logo_width = 0.075
    logo_ax = fig.add_axes([AX_RIGHT - logo_width, 0.865, logo_width, 0.105], facecolor="none", zorder=100)
    logo_ax.imshow(plt.imread(os.path.join(project_root, "utils", "images", "und-logo.png")))
    logo_ax.set_anchor("E")
    logo_ax.axis("off")
    fig.text(AX_RIGHT - logo_width / 2, 0.855, "ATMOSPHERIC SCIENCES", ha="center", va="top", fontsize=7, fontweight="bold", color=TEXT_PRIMARY)

    # footer text and bar
    fig.add_artist(plt.Line2D([AX_LEFT, AX_RIGHT], [0.062, 0.062],
             transform=fig.transFigure, color=DIVIDER_COLOR, linewidth=0.7))
    fig.add_artist(plt.Line2D([AX_LEFT, AX_LEFT + 0.03], [0.062, 0.062],
            transform=fig.transFigure, color=UND_GREEN, linewidth=1.5))
    fig.text(AX_LEFT, 0.033, footer_left, ha="left", va="center", fontsize=8, color=TEXT_MUTED)
    fig.text(AX_RIGHT, 0.033, footer_right, ha="right", va="center", fontsize=8, color=TEXT_MUTED)

    # save settings
    if save_path is not None:
        fig.savefig(save_path, dpi=fig.dpi, facecolor=FIGURE_BG, pad_inches=0)



###########################
# LOAD FORECAST DATA
###########################
now_utc = datetime.now(timezone.utc)

try:
    nws, nws_extremes, nws_issued = get_nws_forecast()
    gefs, gefs_init = get_gefs_plumes()
except Exception as e:
    print(f"ERROR: FORECAST DATA UNAVAILABLE.....{e}")
    sys.exit(1)

# match the GEFS window to the NWS forecast window
start = nws["valid_time"].iloc[0]
end = nws["valid_time"].iloc[-1]
ens_temp = gefs["temperature_f"].loc[start:end]
# precipitation at each hour is the accumulation over the hour before it
ens_precip = gefs["precipitation_in"].loc[start + pd.Timedelta(hours=1):end]

print(f"    FORECAST WINDOW.....{start:%Y-%m-%d %HZ} TO {end:%Y-%m-%d %HZ}")



###########################
# ENSEMBLE STATISTICS
###########################
# member percentiles (actual ensemble distribution, not a normal approximation)
temp_stats = ens_temp.quantile([0.1, 0.5, 0.9], axis=1).T
temp_stats.columns = ["p10", "p50", "p90"]

ens_qpf = ens_precip.cumsum()
qpf_stats = ens_qpf.quantile([0.1, 0.5, 0.9], axis=1).T
qpf_stats.columns = ["p10", "p50", "p90"]

# fraction of members with measurable precip (>= 0.01 in) in the trailing 6 hr
ens_6hr = ens_precip.rolling(6, min_periods=1).sum()
ens_pop = (ens_6hr >= 0.01).mean(axis=1) * 100

# NWS cumulative accumulations
nws["qpf_total_in"] = nws["qpf_in"].cumsum()
nws_snow_total = nws["snow_in"].sum()
nws_ice_total = nws["ice_in"].sum()

print(f"    NWS 7-DAY QPF.....{nws['qpf_total_in'].iloc[-1]:.2f} IN, SNOW {nws_snow_total:.1f} IN")
print(f"    GEFS 7-DAY QPF.....MEDIAN {qpf_stats['p50'].iloc[-1]:.2f} IN, "
      f"10-90% {qpf_stats['p10'].iloc[-1]:.2f}-{qpf_stats['p90'].iloc[-1]:.2f} IN")



###########################
# BUILD FIGURE
###########################
fig = plt.figure(figsize=(14, 5.2), dpi=250, facecolor=FIGURE_BG)

# define all axes: temperature (with wind barbs), a thin sky / PoP / weather strip, and cumulative QPF
ax1 = fig.add_axes([AX_LEFT, 0.455, AX_WIDTH, 0.330])
ax2 = fig.add_axes([AX_LEFT, 0.310, AX_WIDTH, 0.132], sharex=ax1)
ax3 = fig.add_axes([AX_LEFT, 0.145, AX_WIDTH, 0.152], sharex=ax1)

# apply base settings to all axes
for ax in [ax1, ax2, ax3]:
    style_axis(ax)
ax2.grid(False)

ax1.tick_params(axis="x", labelbottom=False)
ax2.tick_params(axis="x", bottom=False, labelbottom=False)
ax3.xaxis.set_major_locator(mdates.HourLocator(byhour=[0, 12], tz=timezone.utc))
ax3.xaxis.set_major_formatter(FuncFormatter(utc_tick_label))
ax3.set_xlim(start, end)


# alternate shading for each local calendar day
local_days = pd.date_range(start.tz_convert(LOCAL_TZ).normalize(), end.tz_convert(LOCAL_TZ), freq="D")
for i, day in enumerate(local_days):
    day_start = max(day.tz_convert("UTC"), start)
    day_end = min((day + pd.Timedelta(days=1)).tz_convert("UTC"), end)
    for ax in [ax1, ax2, ax3]:
        if i % 2 == 0:
            ax.axvspan(day_start, day_end, color=TEXT_PRIMARY, alpha=0.035, linewidth=0, zorder=0)
        ax.axvline(day.tz_convert("UTC"), color=SPINE_COLOR, linewidth=0.8, alpha=0.6, zorder=1)

    # day name centered in the band (skip slivers too narrow to label)
    if (day_end - day_start) >= pd.Timedelta(hours=10):
        ax1.text(day_start + (day_end - day_start) / 2, 0.955, day.strftime("%a %m/%d").upper(), transform=ax1.get_xaxis_transform(),
                 ha="center", va="top", fontsize=8, fontweight="bold", color=TEXT_SECONDARY, zorder=9)



###########################
# TEMPERATURE AXIS
###########################
# GEFS plumes, 10-90% envelope, and median behind the NWS forecast
ax1.plot(ens_temp.index, ens_temp.values, color=ENS_COLOR, linewidth=0.5, alpha=0.20, zorder=2)
ax1.fill_between(temp_stats.index, temp_stats["p10"], temp_stats["p90"], color=ENS_COLOR, alpha=0.10, linewidth=0, zorder=2)
ax1.plot(temp_stats.index, temp_stats["p50"], color=ENS_MEDIAN_COLOR, linewidth=1.3, linestyle="--", alpha=0.75, zorder=3)

# NWS temperature, dewpoint, and wind chill (only where it is colder than the air temp)
wind_chill = nws["wind_chill_f"].where(nws["wind_chill_f"] < nws["temperature_f"] - 0.5)
ax1.plot(nws["valid_time"], nws["dewpoint_f"], color=DEWP_COLOR, linewidth=2.0, zorder=4)
ax1.plot(nws["valid_time"], wind_chill, color=FEEL_COLOR, linewidth=1.8, linestyle=":", zorder=4)
ax1.plot(nws["valid_time"], nws["temperature_f"], color=TEMP_COLOR, linewidth=4.2, zorder=5)
ax1.axhline(32, color=FREEZE_COLOR, linewidth=1.1, linestyle="--", alpha=0.70, zorder=1)

# axis max/min: the coldest line sits above the bottom BARB_BAND of the panel (wind barbs)
# and the warmest below TOP_BAND (day labels and high annotations)
BARB_BAND = 0.24
TOP_BAND = 0.73
temp_low = np.nanmin([nws["dewpoint_f"].min(), wind_chill.min(), temp_stats["p10"].min()])
temp_high = np.nanmax([nws["temperature_f"].max(), temp_stats["p90"].max()])
temp_range = (temp_high - temp_low) / (TOP_BAND - BARB_BAND)
temp_min = temp_low - BARB_BAND * temp_range
temp_max = temp_min + temp_range
ax1.set_ylim(temp_min, temp_max)

# only label temperatures above the barb band
ax1.yaxis.set_major_locator(plt.MaxNLocator(6, integer=True))
ax1.set_yticks([tick for tick in ax1.get_yticks() if temp_low - 3 <= tick <= temp_max])
ax1.set_ylim(temp_min, temp_max)
ax1.set_ylabel("Temperature °F", color=TEXT_PRIMARY, fontsize=8, fontweight="bold", labelpad=4)


# official NWS highs and lows, marked at the warmest/coldest hour of each period
for _, period in nws_extremes.iterrows():
    if period["start"] < start or period["end"] > end + pd.Timedelta(hours=1):
        continue

    in_period = nws.loc[(nws["valid_time"] >= period["start"]) & (nws["valid_time"] < period["end"])]
    if in_period.empty:
        continue

    if period["kind"] == "H":
        idx, color, marker, offset, va = in_period["temperature_f"].idxmax(), HIGH_COLOR, "^", (0, 8), "bottom"
    else:
        idx, color, marker, offset, va = in_period["temperature_f"].idxmin(), LOW_COLOR, "v", (0, -8), "top"

    x = nws.loc[idx, "valid_time"]
    y = nws.loc[idx, "temperature_f"]
    ax1.scatter(x, y, marker=marker, s=40, color=color, edgecolors=PANEL_BG, linewidths=0.9, zorder=7)
    ax1.annotate(f"{period['kind']} {period['temperature_f']:.0f}°", xy=(x, y), xytext=offset, textcoords="offset points",
                 ha="center", va=va, fontsize=9, fontweight="bold", color=TEXT_PRIMARY, path_effects=TEXT_OUTLINE, zorder=8)


# NWS wind barbs (kt) every 6 hr along the bottom of the panel; gusts >= 25 kt labeled
barb_hours = nws.loc[nws["valid_time"].dt.hour % 6 == 0].dropna(subset=["wind_speed_kt", "wind_dir_deg"])
barb_dir = np.deg2rad(barb_hours["wind_dir_deg"].to_numpy())
barb_u = -barb_hours["wind_speed_kt"].to_numpy() * np.sin(barb_dir)
barb_v = -barb_hours["wind_speed_kt"].to_numpy() * np.cos(barb_dir)
barb_x = mdates.date2num(barb_hours["valid_time"].dt.tz_localize(None))
ax1.barbs(barb_x, np.full(len(barb_x), BARB_BAND / 2), barb_u, barb_v, transform=ax1.get_xaxis_transform(),
          length=5.2, linewidth=0.8, pivot="middle", color=WIND_COLOR, zorder=6, clip_on=False)

for _, hour in barb_hours.iterrows():
    if hour["wind_gust_kt"] >= 25:
        ax1.text(hour["valid_time"] + pd.Timedelta(hours=1.5), BARB_BAND / 2, f"G{hour['wind_gust_kt']:.0f}",
                 transform=ax1.get_xaxis_transform(), ha="left", va="center", fontsize=6, fontweight="bold",
                 color=TEXT_PRIMARY, path_effects=TEXT_OUTLINE, zorder=7)

temp_handles = [
    Line2D([], [], color=TEMP_COLOR, linewidth=4.2, label="NWS TEMP"),
    Line2D([], [], color=DEWP_COLOR, linewidth=2.0, label="NWS DEWPOINT"),
    Line2D([], [], color=FEEL_COLOR, linewidth=1.8, linestyle=":", label="NWS WIND CHILL"),
    Line2D([], [], color=HIGH_COLOR, marker="^", linestyle="None", markersize=6, label="NWS HIGH"),
    Line2D([], [], color=LOW_COLOR, marker="v", linestyle="None", markersize=6, label="NWS LOW"),
    Line2D([], [], color=WIND_COLOR, linewidth=1.0, label="NWS WIND (KT)"),
    Line2D([], [], color=ENS_COLOR, linewidth=1.0, alpha=0.5, label="GEFS MEMBERS"),
    Line2D([], [], color=ENS_MEDIAN_COLOR, linewidth=1.3, linestyle="--", alpha=0.75, label="GEFS MEDIAN"),
    Patch(facecolor=ENS_COLOR, alpha=0.22, label="GEFS 10–90%"),
]
side_legend(ax1, temp_handles)



###########################
# PRECIP PROBABILITY STRIP
###########################
# four rows, top to bottom: NWS sky cover (brighter = cloudier), NWS non-precip weather
# (fog, haze, blowing snow...) as METAR-style codes, NWS PoP colored by forecast precip
# type, and GEFS PoP as the % of members with >= 0.01 in in the trailing 6 hr
# (PoP cell opacity scales with probability)
def pop_alpha(pop):
    return np.where(pop >= 5, 0.18 + 0.82 * np.clip(pop, 0, 100) / 100, 0.0)

hour_width = 1 / 24
nws_colors = [to_rgba(PRECIP_TYPE_COLORS.get(ptype, POP_COLOR), a)
              for ptype, a in zip(nws["precip_type"], pop_alpha(nws["pop_percent"].to_numpy()))]
ax2.bar(nws["valid_time"], 0.86, bottom=1.07, width=hour_width, color=nws_colors, align="edge", linewidth=0, zorder=3)

sky_alpha = 0.06 + 0.84 * np.clip(nws["sky_cover_percent"].fillna(0).to_numpy(), 0, 100) / 100
ax2.bar(nws["valid_time"], 0.86, bottom=3.07, width=hour_width, color=[to_rgba(SKY_COLOR, a) for a in sky_alpha],
        align="edge", linewidth=0, zorder=3)

# daily average sky cover labeled in each local day band
sky = nws.set_index("valid_time")["sky_cover_percent"]
for _, day_sky in sky.groupby(sky.index.tz_convert(LOCAL_TZ).date):
    if len(day_sky) >= 10 and day_sky.notna().any():
        ax2.text(day_sky.index[0] + (day_sky.index[-1] - day_sky.index[0]) / 2 + pd.Timedelta(hours=0.5), 3.5,
                 f"{day_sky.mean():.0f}%", ha="center", va="center", fontsize=6.5, fontweight="bold",
                 color=TEXT_PRIMARY, path_effects=TEXT_OUTLINE, zorder=6, clip_on=True)

ens_colors = [to_rgba(ENS_MEDIAN_COLOR, a) for a in pop_alpha(ens_pop.to_numpy())]
ax2.bar(ens_pop.index - pd.Timedelta(hours=1), 0.86, bottom=0.07, width=hour_width, color=ens_colors,
        align="edge", linewidth=0, zorder=3)

# non-precip weather: one cell per run of hours with the same codes, labeled when wide enough
wx = nws.set_index("valid_time")["obstruction"]
wx_runs = (wx != wx.shift()).cumsum()
for _, run in wx.groupby(wx_runs):
    if run.iloc[0] == "":
        continue
    run_start = run.index[0]
    run_hours = len(run)
    ax2.bar(run_start, 0.86, bottom=2.07, width=run_hours * hour_width, color=WX_COLOR, align="edge", linewidth=0, zorder=3)
    if run_hours >= 3:
        ax2.text(run_start + pd.Timedelta(hours=run_hours / 2), 2.5, run.iloc[0], ha="center", va="center", fontsize=5.8,
                 fontweight="bold", color=TEXT_PRIMARY, path_effects=TEXT_OUTLINE, zorder=6, clip_on=True)

if not (wx != "").any():
    ax2.text(0.5, 2.5, "NO FOG, HAZE, SMOKE, OR BLOWING SNOW FORECAST", transform=ax2.get_yaxis_transform(),
             ha="center", va="center", fontsize=5.8, fontweight="bold", color=TEXT_MUTED, zorder=6)

ax2.set_ylim(0, 4)
ax2.set_yticks([0.5, 1.5, 2.5, 3.5])
ax2.set_yticklabels(["GEFS PoP", "NWS PoP", "WX", "SKY"])
ax2.tick_params(axis="y", labelsize=7)
for row_edge in (1, 2, 3):
    ax2.axhline(row_edge, color=SPINE_COLOR, linewidth=0.6, zorder=2)

# label the peak probability in each row for each local day, keeping only the
# higher of two peaks less than 12 hr apart (one event spanning midnight)
pop_rows = ((nws.set_index("valid_time")["pop_percent"], 1.5, 0), (ens_pop, 0.5, -1))
for series, row_y, shift_hr in pop_rows:
    local_date = series.index.tz_convert(LOCAL_TZ).date
    peaks = [(day_pop.idxmax(), day_pop.max()) for _, day_pop in series.groupby(local_date) if day_pop.max() >= 20]

    kept = []
    for peak_time, peak_pop in peaks:
        if kept and peak_time - kept[-1][0] < pd.Timedelta(hours=12):
            if peak_pop > kept[-1][1]:
                kept[-1] = (peak_time, peak_pop)
            continue
        kept.append((peak_time, peak_pop))

    for peak_time, peak_pop in kept:
        x = peak_time + pd.Timedelta(hours=shift_hr + 0.5)
        ax2.text(x, row_y, f"{peak_pop:.0f}%", ha="center", va="center", fontsize=6.5, fontweight="bold",
                 color=TEXT_PRIMARY, path_effects=TEXT_OUTLINE, zorder=6, clip_on=True)

# legend lists only the precip types actually drawn (PoP >= 5%); faint gray low-PoP
# cells and the self-labeled WX row are left out so a busy week still fits the column
drawn_types = set(nws.loc[nws["pop_percent"] >= 5, "precip_type"]) - {""}
pop_handles = [Patch(facecolor=color, label=PRECIP_TYPE_LABELS[ptype])
               for ptype, color in PRECIP_TYPE_COLORS.items() if ptype in drawn_types]
pop_handles.append(Patch(facecolor=ENS_MEDIAN_COLOR, label="GEFS PoP (6H ≥0.01\")"))
side_legend(ax2, pop_handles, anchor_y=1.0)



###########################
# CUMULATIVE QPF AXIS
###########################
ax3.plot(ens_qpf.index, ens_qpf.values, color=ENS_COLOR, linewidth=0.6, alpha=0.35, zorder=2)
ax3.fill_between(qpf_stats.index, qpf_stats["p10"], qpf_stats["p90"], color=ENS_COLOR, alpha=0.14, linewidth=0, zorder=2)
ax3.plot(qpf_stats.index, qpf_stats["p50"], color=ENS_MEDIAN_COLOR, linewidth=1.6, linestyle="--", zorder=3)
ax3.plot(nws["valid_time"], nws["qpf_total_in"], color=QPF_COLOR, linewidth=2.6, drawstyle="steps-post", zorder=4)

# scale to the wettest member so no plume is clipped, with headroom for the totals row
qpf_max = max(0.25, float(ens_qpf.iloc[-1].max()), float(nws["qpf_total_in"].iloc[-1])) * 1.35
ax3.set_ylim(0, qpf_max)
ax3.yaxis.set_major_locator(plt.MaxNLocator(3))
ax3.set_ylabel("QPF (in)", color=TEXT_PRIMARY, fontsize=8, fontweight="bold", labelpad=4)

qpf_handles = [
    Line2D([], [], color=QPF_COLOR, linewidth=2.6, label="NWS TOTAL QPF"),
    Line2D([], [], color=ENS_COLOR, linewidth=1.0, alpha=0.7, label="GEFS MEMBERS"),
    Line2D([], [], color=ENS_MEDIAN_COLOR, linewidth=1.6, linestyle="--", label="GEFS MEDIAN"),
    Patch(facecolor=ENS_COLOR, alpha=0.3, label="GEFS 10–90%"),
]
side_legend(ax3, qpf_handles, loc="lower left", anchor_y=0.0)

# 7-day totals in one row across the top of the QPF axis
totals_str = (f"7-DAY TOTALS     NWS QPF {nws['qpf_total_in'].iloc[-1]:.2f}\"     NWS SNOW {nws_snow_total:.1f}\"     "
    + (f"NWS ICE {nws_ice_total:.2f}\"     " if nws_ice_total >= 0.01 else "")
    + f"GEFS MEDIAN {qpf_stats['p50'].iloc[-1]:.2f}\"     GEFS 10% {qpf_stats['p10'].iloc[-1]:.2f}\"     GEFS 90% {qpf_stats['p90'].iloc[-1]:.2f}\"")
ax3.text(0.006, 0.93, totals_str, transform=ax3.transAxes, ha="left", va="top", fontsize=7, fontweight="bold",
         color=TEXT_PRIMARY, path_effects=TEXT_OUTLINE, zorder=10)

for ax in [ax1, ax2, ax3]:
    bold_ticks(ax)



###########################
# FOOTER TEXT
###########################
gefs_str = f"{gefs_init:%HZ %d %b %Y}".upper() if gefs_init is not None else "LATEST"
nws_str = f"{nws_issued:%HZ %d %b %Y}".upper()
textstr = f"NWS FGF Forecast Issued {nws_str}   •   GEFS 0.25° Init {gefs_str} ({ens_temp.shape[1]} Members)"
source_str = "UND Atmospheric Sciences • Weather Wall"



###########################
# BUILD FILENAME
###########################
filename = build_filename("staged_figures/gfk_forecast/", "gfk_forecast_meteogram", now_utc)
os.makedirs(os.path.dirname(filename), exist_ok=True)



###########################
# BUILD FIG
###########################
meteogram_frame(
    fig,
    title="Grand Forks 7 Day Forecast",
    subtitle=f"KGFK Grand Forks Intl. • NWS FGF Forecast • GEFS {ens_temp.shape[1]}-Member Ensemble",
    valid=(f"Valid • {format_utc_label(start)} – {format_utc_label(end)}"),
    footer_left=textstr,
    footer_right=source_str,
    save_path=filename,
)

plt.close(fig)
print(f"    WROTE {filename}")


elapsed_time = comp_time.time() - st

print("############\n"
    f"SCRIPT FINISHED: time: "
    f"{comp_time.strftime('%H:%M:%S', comp_time.gmtime(elapsed_time))}\n"
    "############")
