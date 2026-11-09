##########################################################
#       GFK NWS / GEFS FORECAST METEOGRAM
#  UNIVERSITY OF NORTH DAKOTA - ATMOSPHERIC SCIENCES
##########################################################

import time as comp_time
st = comp_time.time()

import os
import sys
import warnings
from datetime import datetime, timezone

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.patheffects import withStroke
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=FutureWarning, module="cfgrib")

print("############\nSCRIPT RUNNING\n############")

# Find project root so this script can also be run directly.
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(script_dir, ".."))
if project_root not in sys.path:
    sys.path.append(project_root)

# Import project data and figure helpers.
from get_data.get_gfk_forecast import get_gefs_temperature_ensemble, get_nws_hourly_forecast, get_nws_qpf_forecast
from utils.style import (DEWP_COLOR, DIVIDER_COLOR, FIGURE_BG, GRID_COLOR, PANEL_BG, SPINE_COLOR, TEXT_PRIMARY,
                         TEXT_SECONDARY, TEXT_MUTED, TEXT_VALID, TEMP_COLOR, FREEZE_COLOR, UND_GREEN)
from utils.utils import build_filename


def add_meteogram_frame(fig, *, title, subtitle, valid, footer_left, footer_right):
    """Add a frame sized for the weather-wall's wide forecast-timeseries slot."""
    frame_left = 0.065
    frame_right = 0.88
    border = Rectangle( (0.002, 0.002), 0.996, 0.996, transform=fig.transFigure, fill=False, edgecolor=UND_GREEN, linewidth=4.0, zorder=10000, clip_on=False)
    fig.add_artist(border)

    fig.text(frame_left, 0.973, title, ha="left", va="top",  fontsize=17, fontweight="bold", color=TEXT_PRIMARY)
    fig.text(frame_left, 0.918, subtitle, ha="left", va="top", fontsize=10.5, color=TEXT_SECONDARY)
    fig.text(frame_left, 0.879, valid, ha="left", va="top", fontsize=8.5, fontweight="bold", color=TEXT_VALID)
    fig.add_artist(plt.Line2D([frame_left, frame_left + 0.063], [0.850, 0.850], transform=fig.transFigure, color=UND_GREEN, linewidth=2.2, solid_capstyle="round"))

    logo_path = os.path.join(project_root, "utils", "images", "und-logo.png")
    logo_width, logo_height = 0.095, 0.070
    logo_left = 0.75
    logo_center = logo_left + logo_width / 2
    logo_ax = fig.add_axes([logo_left, 0.905, logo_width, logo_height], zorder=100)
    logo_ax.imshow(plt.imread(logo_path))
    logo_ax.axis("off")
    fig.text(logo_center, 0.890, "ATMOSPHERIC SCIENCES", ha="center", va="top",
             fontsize=8.5, fontweight="bold", color=TEXT_PRIMARY)

    fig.add_artist(plt.Line2D([frame_left, frame_right], [0.078, 0.078], transform=fig.transFigure, color=DIVIDER_COLOR, linewidth=0.7))
    fig.add_artist(plt.Line2D([frame_left, frame_left + 0.03], [0.078, 0.078], transform=fig.transFigure, color=UND_GREEN, linewidth=1.5,))
    fig.text(frame_left, 0.043, footer_left, ha="left", va="center", fontsize=8.2, color=TEXT_MUTED)
    fig.text(frame_right, 0.043, footer_right, ha="right", va="center", fontsize=8.2, color=TEXT_MUTED)


now_utc = datetime.now(timezone.utc)

# ============================================================
# LOAD FORECAST DATA
# ============================================================

print("    LOADING NWS AND GEFS FORECASTS...")
nws = get_nws_hourly_forecast()
nws_qpf = get_nws_qpf_forecast()
gefs, gefs_cycle = get_gefs_temperature_ensemble()
print("    FORECAST DATA LOADED")


def format_utc_label(timestamp):
    """Format a valid time like ``18Z MON 09 NOV 2026``."""
    return timestamp.strftime("%HZ %a %d %b %Y").upper()

# Keep the ensemble within the NWS forecast period for a matched time series.
forecast_start = nws["valid_time"].min()
forecast_end = nws["valid_time"].max()
gefs = gefs.loc[(gefs["valid_time"] >= forecast_start) & (gefs["valid_time"] <= forecast_end)].copy()

# Accumulate the ensemble-mean QPF intervals. Estimate total spread by adding
# interval variances, assuming forecast-period precipitation errors are
# independent.
gefs = gefs.sort_values("valid_time").reset_index(drop=True)
gefs["cumulative_precipitation_in"] = gefs["precipitation_in"].fillna(0).cumsum()
gefs["cumulative_precipitation_spread_in"] = np.sqrt(gefs["precipitation_spread_in"].fillna(0).pow(2).cumsum())

# NWS gridded QPF values are interval totals. Summing them at their interval
# end times gives a cumulative forecast comparable to the GEFS cumulative QPF.
nws_qpf = nws_qpf.loc[(nws_qpf["valid_time"] >= forecast_start) & (nws_qpf["valid_time"] <= forecast_end)].sort_values("valid_time").reset_index(drop=True)
if not nws_qpf.empty:
    nws_qpf["cumulative_precipitation_in"] = nws_qpf["precipitation_in"].cumsum()
    print(f"    MAX NWS CUMULATIVE QPF: {nws_qpf['cumulative_precipitation_in'].iloc[-1]:.2f} IN")
else:
    print("    NWS QPF: NO VALID VALUES; PLOTTING GEFS QPF ONLY")







# ============================================================
# CALCULATE GEFS ENSEMBLE SUMMARIES
# ============================================================
temperature_stats = gefs[["valid_time", "temperature_f", "temperature_spread_f"]].copy()
temperature_stats["p10"] = (temperature_stats["temperature_f"] - 1.28155 * temperature_stats["temperature_spread_f"])
temperature_stats["p90"] = ( temperature_stats["temperature_f"] + 1.28155 * temperature_stats["temperature_spread_f"])

precipitation_stats = gefs[["valid_time", "cumulative_precipitation_in", "cumulative_precipitation_spread_in"]].copy()
precipitation_stats["p10"] = np.maximum(0, precipitation_stats["cumulative_precipitation_in"]- 1.28155 * precipitation_stats["cumulative_precipitation_spread_in"])
precipitation_stats["p90"] = (precipitation_stats["cumulative_precipitation_in"] + 1.28155 * precipitation_stats["cumulative_precipitation_spread_in"])

max_nws_pop = float(nws["pop_percent"].max())
max_gefs_qpf = float(precipitation_stats["p90"].max())
print(f"    MAX NWS PROBABILITY OF PRECIPITATION: {max_nws_pop:.0f}%")
if np.isfinite(max_gefs_qpf):
    print(f"    MAX GEFS 90TH-PERCENTILE CUMULATIVE QPF: {max_gefs_qpf:.2f} IN")
else:
    print("    GEFS QPF: NO VALID PRECIPITATION VALUES")






# ============================================================
# SET UP THE WIDE, SHORT FIGURE
# ============================================================
fig, (ax_temp, ax_pop) = plt.subplots(2,1,figsize=(14, 5.2), dpi=250, sharex=True, facecolor=FIGURE_BG, gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.12})
fig.subplots_adjust(left=0.065, right=0.87, top=0.835, bottom=0.17)
ax_temp.set_facecolor(PANEL_BG)
ax_pop.set_facecolor(PANEL_BG)

for ax in (ax_temp, ax_pop):
    ax.tick_params(axis="both", colors=TEXT_SECONDARY, labelsize=9, pad=4)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight("bold")
    for spine in ax.spines.values():
        spine.set_color(SPINE_COLOR)
        spine.set_linewidth(1.1)
    ax.grid(True, color=GRID_COLOR, linewidth=0.7, alpha=0.38)



# ============================================================
# TEMPERATURE, DEWPOINT, AND WIND CHILL
# ============================================================
ax_temp.fill_between(
    temperature_stats["valid_time"],
    temperature_stats["p10"],
    temperature_stats["p90"],
    color="#87AFC4",
    alpha=0.24,
    label="GEFS approx. 10–90%",
    zorder=3,
)
ax_temp.plot(
    temperature_stats["valid_time"], temperature_stats["temperature_f"],
    color="#8CC7E8", linewidth=2.5, label="GEFS mean", zorder=4,
)
nws_temp_line, = ax_temp.plot(
    nws["valid_time"], nws["temperature_f"],
    color=TEMP_COLOR, linewidth=3.4, label="NWS temperature", zorder=6,
)
nws_dew_line, = ax_temp.plot(
    nws["valid_time"], nws["dewpoint_f"],
    color=DEWP_COLOR, linewidth=2.2, linestyle="-.", label="NWS dewpoint", zorder=5,
)
if nws["wind_chill_f"].notna().any():
    ax_temp.plot(
        nws["valid_time"], nws["wind_chill_f"],
        color="#80C7E8", linewidth=2.1, linestyle=":", label="NWS wind chill", zorder=5,
    )
ax_temp.axhline(
    32, color=FREEZE_COLOR, linewidth=1.2, linestyle="--", alpha=0.32,
    label="32°F", zorder=1,
)
ax_temp.set_ylabel("Temperature (°F)", color=TEXT_PRIMARY, fontsize=10, fontweight="bold", labelpad=7)
ax_temp.tick_params(axis="x", labelbottom=False)
# Keep extra headroom above the warmest forecast for the high markers and legend.
ax_temp.margins(y=0.27)

# Annotate each local calendar day's NWS forecast high and low.
local_dates = nws["valid_time"].dt.tz_convert("America/Chicago").dt.date
outlined_text = [withStroke(linewidth=2.5, foreground=PANEL_BG)]
for _, day in nws.assign(local_date=local_dates).groupby("local_date", sort=True):
    for idx, label, color, offset, va, marker in (
        (day["temperature_f"].idxmax(), "H", "#FF554F", (0, 10), "bottom", "^"),
        (day["temperature_f"].idxmin(), "L", "#168BFF", (0, -10), "top", "v"),
    ):
        valid_time = nws.loc[idx, "valid_time"]
        temperature = nws.loc[idx, "temperature_f"]
        marker_edge = "#111111" if label == "H" else TEXT_PRIMARY
        ax_temp.scatter(
            valid_time, temperature, marker=marker, s=34,
            color=color, edgecolors=marker_edge, linewidths=1.1, zorder=8,
        )
        ax_temp.annotate(
            f"{label} {temperature:.0f}°",
            xy=(valid_time, temperature), xytext=offset,
            textcoords="offset points", ha="center", va=va,
            color=color, fontsize=8.0, fontweight="bold",
            path_effects=outlined_text, clip_on=True, zorder=9,
        )

temp_handles, temp_labels = ax_temp.get_legend_handles_labels()
temp_handles.extend([
    Line2D([], [], color="#FF554F", marker="^", linestyle="None", markersize=6, label="NWS high"),
    Line2D([], [], color="#168BFF", marker="v", linestyle="None", markersize=6, label="NWS low"),
])
legend = ax_temp.legend(
    temp_handles, temp_labels + ["NWS high", "NWS low"],
    loc="upper left", ncol=7, frameon=False, fontsize=7.0,
    handlelength=2.0, columnspacing=1.0, borderaxespad=0.5,
)
for label in legend.get_texts():
    label.set_color(TEXT_PRIMARY)
    label.set_fontweight("bold")

# ============================================================
# NWS PRECIPITATION PROBABILITY, TYPE, AND QPF
# ============================================================

pop_bars = ax_pop.bar(
    nws["valid_time"], nws["pop_percent"],
    width=0.034, color="#65A9D1", alpha=0.82,
    label="NWS probability", align="center", zorder=3,
)
ax_pop.set_ylim(0, 112)
ax_pop.set_ylabel("NWS PoP (%)", color=TEXT_PRIMARY, fontsize=9, fontweight="bold", labelpad=7)
ax_pop.set_yticks([0, 25, 50, 75, 100])

ax_qpf = ax_pop.twinx()
ax_qpf.patch.set_visible(False)
ax_qpf.tick_params(axis="y", colors=TEXT_SECONDARY, labelsize=8, pad=4)
for spine in ax_qpf.spines.values():
    spine.set_visible(False)
ax_qpf.fill_between(
    precipitation_stats["valid_time"],
    precipitation_stats["p10"],
    precipitation_stats["p90"],
    color="#E5B35A", alpha=0.16, zorder=1,
)
qpf_line, = ax_qpf.plot(
    precipitation_stats["valid_time"], precipitation_stats["cumulative_precipitation_in"],
    color="#F0C66E", linewidth=2.4, marker="o", markersize=2.7,
    label="GEFS mean cumulative QPF", zorder=4,
)
qpf_handles = [qpf_line]
qpf_labels = ["GEFS mean cumulative QPF"]
qpf_axis_max = float(precipitation_stats["p90"].max())
if not nws_qpf.empty:
    nws_qpf_line, = ax_qpf.plot(
        nws_qpf["valid_time"], nws_qpf["cumulative_precipitation_in"],
        color="#F28E8B", linewidth=2.5, drawstyle="steps-post",
        label="NWS cumulative QPF", zorder=5,
    )
    qpf_handles.append(nws_qpf_line)
    qpf_labels.append("NWS cumulative QPF")
    qpf_axis_max = max(
        qpf_axis_max,
        float(nws_qpf["cumulative_precipitation_in"].max()),
    )
ax_qpf.set_ylabel("Cumulative QPF (in)", color="#E5B35A", fontsize=9, fontweight="bold", labelpad=7)
max_qpf = qpf_axis_max
ax_qpf.set_ylim(0, max(0.1, max_qpf * 1.25))

# Mark the NWS text forecast's precipitation type when precipitation is
# mentioned. The forecast wording is the NWS categorical guidance.
precip_colors = {
    "Rain": "#5AA9E6",
    "Snow": "#F2F2EE",
    "Wintry Mix": "#C69AE8",
    "Thunderstorms": "#F2C14E",
}
type_handles = []
for precip_type, color in precip_colors.items():
    mask = (nws["precip_type"] == precip_type) & (nws["pop_percent"] > 0)
    if mask.any():
        ax_pop.scatter(
            nws.loc[mask, "valid_time"],
            np.full(mask.sum(), 106),
            s=14,
            marker="s",
            color=color,
            edgecolors=FIGURE_BG,
            linewidths=0.8,
            zorder=6,
        )
        type_handles.append(
            Line2D([], [], color=color, marker="s", linestyle="None", markersize=5, label=precip_type)
        )

ax_pop.xaxis.set_major_locator(mdates.HourLocator(byhour=[0, 12], tz=timezone.utc))
ax_pop.xaxis.set_major_formatter(mdates.DateFormatter("%a %m/%d\n%HZ", tz=timezone.utc))
ax_pop.set_xlim(forecast_start, forecast_end)

legend_handles = [pop_bars, *qpf_handles]
legend_labels = ["NWS probability", *qpf_labels]
legend_handles.extend(type_handles)
legend_labels.extend([handle.get_label() for handle in type_handles])
legend = ax_pop.legend( legend_handles, legend_labels, loc="upper left", ncol=min(6, len(legend_handles)), frameon=False, fontsize=7.6, handlelength=1.6, columnspacing=1.0)
for label in legend.get_texts():
    label.set_color(TEXT_PRIMARY)
    label.set_fontweight("bold")





# ============================================================
# SAVE THE METEOGRAM
# ============================================================
filename = build_filename("staged_figures/gfk_forecast/", "gfk_forecast_meteogram", now_utc)
os.makedirs(os.path.dirname(filename), exist_ok=True)
add_meteogram_frame(fig, title="Grand Forks 7 Day Forecast", 
                    subtitle="NWS FGF Forecast & GEFS Ensemble Mean/Spread for KGFK (Grand Forks Intl. Airport)",
                    valid=(f"Valid: {format_utc_label(forecast_start)} - {format_utc_label(forecast_end)}"),
                    footer_left=(f"GEFS 0.25° - INIT {pd.Timestamp(gefs_cycle).strftime('%d %b %Y %HZ').upper()} • KGFK NWS FGF Forecast - VALID {format_utc_label(forecast_start)}"),
                    footer_right="UND Atmospheric Sciences • Weather Wall")

fig.savefig(filename, dpi=fig.dpi, facecolor=FIGURE_BG, pad_inches=0)
plt.close(fig)
print(f"    WROTE {filename}")

elapsed_time = comp_time.time() - st
print("SCRIPT COMPLETE. Time elapsed:", comp_time.strftime("%H:%M:%S", comp_time.gmtime(elapsed_time)))
