
##########################################################
#              BASIC FIGURE CREATION SCRIPT
#  (c) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2026
##########################################################

from pathlib import Path
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.colorbar import ColorbarBase

# get script dir
script_dir = os.path.dirname(os.path.abspath(__file__))
# get the parent dir 
project_root = os.path.abspath(os.path.join(script_dir, ".."))
if project_root not in sys.path:
    sys.path.append(project_root)
from utils.utils import *
from utils.style import *


def figure_builder(
    fig,
    ax,
    *,
    title,
    subtitle,
    valid=None,
    mappable=None,
    cbar_title=None,
    cbar_units=None,
    cbar_ticks=None,
    cbar_extend="both",
    category_colors=None,
    category_labels=None,
    category_title=None,
    footer_left=None,
    footer_right="UND Atmospheric Sciences • Weather Wall",
    logo_path="utils/images/und-logo.png",
    save_path=None):

    # outer frame
    border = Rectangle((0.002, 0.002), 0.996, 0.996, transform=fig.transFigure,
        fill=False, edgecolor=UND_GREEN, linewidth=6, zorder=10000, clip_on=False)
    fig.add_artist(border)


    # header text
    fig.text( 0.015, 0.973, title, ha="left", va="top", fontsize=24, fontweight="bold", color=TEXT_PRIMARY)
    fig.text(0.015, 0.937, subtitle, ha="left", va="top", fontsize=15, fontweight="medium", color=TEXT_SECONDARY)
    fig.text(0.015, 0.907, valid, ha="left", va="top", fontsize=11, fontweight="bold", color=TEXT_VALID)
    fig.add_artist(plt.Line2D([0.015, 0.115], [0.884, 0.884],
            transform=fig.transFigure, color=UND_GREEN, linewidth=3, solid_capstyle="round"))

    # und logo and text
    logo_path = Path(logo_path)
    LOGO_RECT = [0.875,0.905,0.095,0.070]
    logo_ax = fig.add_axes(LOGO_RECT,facecolor="none", zorder=100,)
    logo_ax.imshow(plt.imread(logo_path))
    logo_ax.axis("off")
    fig.text(0.86, 0.890, f'ATMOSPHERIC SCIENCES', ha='left', weight='bold', fontsize=10, color='white')



    # colorbar and legend settings
    cbar = None

    if mappable is not None:
        # if none assume colorbar
        COLORBAR_RECT = [0.915, 0.08, 0.016, 0.770]
        cax = fig.add_axes(COLORBAR_RECT, facecolor=FIGURE_BG)

        cbar = fig.colorbar(mappable,cax=cax, orientation="vertical", ticks=cbar_ticks, extend=cbar_extend, extendrect=True)
        cbar.ax.tick_params(axis="y", colors=TEXT_PRIMARY, labelsize=9, length=0, pad=6)

        for tick in cbar.ax.get_yticklabels():
            tick.set_fontweight("bold")

        cbar.outline.set_edgecolor(COLORBAR_EDGE)
        cbar.outline.set_linewidth(0.9)

        if cbar_title is not None:
            label = cbar_title

            if cbar_units is not None:
                label += f" ({cbar_units})"

            cax.text(2.8, 0.5, label, ha="left", va="center", rotation=270, color=TEXT_PRIMARY, fontsize=11, fontweight="bold", transform=cax.transAxes)



    elif category_colors is not None and category_labels is not None:
        # add legend if legend settings are passed
        LEGEND_RECT = [0.91, 0.15, 0.12, 0.62]
        lax = fig.add_axes(LEGEND_RECT, facecolor="none")
        lax.set_xlim(0, 1)
        lax.set_ylim(0, 1)
        lax.axis("off")

        bar_x = 0.06
        bar_y = 0.04
        bar_w = 0.34
        bar_h = 0.92

        n = len(category_labels)
        seg_h = bar_h / n

        for i, (color, label) in enumerate(zip(category_colors, category_labels)):

            y0 = bar_y + i * seg_h
            rect = Rectangle((bar_x, y0), bar_w, seg_h, facecolor=color, edgecolor=COLORBAR_EDGE, linewidth=0.8, transform=lax.transAxes,)
            lax.add_patch(rect)
            lax.text(bar_x + bar_w / 2, y0 + seg_h / 2, label, ha="center", va="center", fontsize=10.5, fontweight="bold", color="#111111", transform=lax.transAxes)

        outer = Rectangle((bar_x, bar_y), bar_w, bar_h, fill=False, edgecolor=COLORBAR_EDGE, linewidth=1.0, transform=lax.transAxes)
        lax.add_patch(outer)

        if category_title is not None:
            lax.text(0.96, 0.5, category_title, ha="center", va="center", rotation=270, color=TEXT_PRIMARY,
                fontsize=11, fontweight="bold", transform=lax.transAxes)


    # footer text and bar
    fig.add_artist(plt.Line2D([0.015, 0.905], [0.055, 0.055],
             transform=fig.transFigure, color=DIVIDER_COLOR, linewidth=0.8))
    fig.add_artist(plt.Line2D([0.015, 0.055],[0.055, 0.055],
            transform=fig.transFigure, color=UND_GREEN, linewidth=1.5,))
    fig.text(0.015, 0.04, footer_left, ha="left", va="center", fontsize=12, color=TEXT_MUTED)
    fig.text(0.905, 0.04, footer_right, ha="right", va="center", fontsize=12, color=TEXT_MUTED)


    # save settings
    if save_path is not None:
        fig.savefig(save_path, dpi=fig.dpi, facecolor=FIGURE_BG, pad_inches=0)

    return cbar