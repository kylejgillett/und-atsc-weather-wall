import matplotlib.patches as mpatches
import cartopy.crs as ccrs

from get_data.get_nws_headlines import get_nws_headlines
from utils.nws_alert_colors import (
    NWS_ALERT_COLORS,
    NWS_ALERT_DEFAULT_COLOR,
)


# Radar-style colors for short-fused warnings
SBW_COLORS = {
    "Tornado Warning": "#FF0000",
    "Severe Thunderstorm Warning": "#FFFF00",
    "Flash Flood Warning": "#00FF00",
    "Special Marine Warning": "#00FFFF",
    "Snow Squall Warning": "#8A2BE2",
}


def _add_colors(gdf, colors):

    gdf = gdf.copy()
    gdf["_plot_color"] = (gdf["event"].map(colors).fillna(NWS_ALERT_DEFAULT_COLOR))
    return gdf


def add_nws_headlines(ax,bbox=None,show_wwa=True,show_sbw=True,wwa_alpha=0.35,sbw_alpha=0.65,edgecolor="black",linewidth=0.4,zorder=20,legend=False):

    if bbox is None:

        xmin, xmax, ymin, ymax = ax.get_extent(
            crs=ccrs.PlateCarree()
        )

        bbox = (
            xmin,
            ymin,
            xmax,
            ymax,
        )

    sbw, wwa = get_nws_headlines(bbox=bbox)


    wwa = _add_colors( wwa,NWS_ALERT_COLORS,)
    sbw = _add_colors(sbw,SBW_COLORS)

    if show_wwa and not wwa.empty:

        wwa.plot(
            ax=ax,
            color=wwa["_plot_color"],
            edgecolor=edgecolor,
            linewidth=linewidth,
            alpha=wwa_alpha,
            transform=ccrs.PlateCarree(),
            zorder=zorder)

    if show_sbw and not sbw.empty:

        sbw.plot(
            ax=ax,
            color=sbw["_plot_color"],
            edgecolor=edgecolor,
            linewidth=linewidth + 0.1,
            alpha=sbw_alpha,
            transform=ccrs.PlateCarree(),
            zorder=zorder + 0.1)

    if legend:

        handles = {}

        if show_wwa:

            for event in wwa["event"].dropna().unique():

                color = NWS_ALERT_COLORS.get(
                    event,
                    NWS_ALERT_DEFAULT_COLOR,
                )

                handles[event] = mpatches.Patch(
                    facecolor=color,
                    edgecolor=edgecolor,
                    label=event,
                    alpha=wwa_alpha,
                )

        if show_sbw:

            for event in sbw["event"].dropna().unique():

                color = SBW_COLORS.get(
                    event,
                    NWS_ALERT_DEFAULT_COLOR,
                )

                handles[event] = mpatches.Patch(
                    facecolor=color,
                    edgecolor=edgecolor,
                    label=event,
                    alpha=sbw_alpha,
                )

        if handles:

            ax.legend(
                handles=list(handles.values()),
                title="NWS Alerts",
                loc="lower left",
                fontsize=8,
                title_fontsize=9,
                framealpha=0.9,
            )

    return sbw, wwa