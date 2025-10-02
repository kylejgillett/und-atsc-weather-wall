##########################################################
#              WPC FRONTS RETREVIAL SCRIPT
#  (c) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2025
##########################################################

from siphon.catalog import TDSCatalog
from metpy.io import parse_wpc_surface_bulletin
from metpy.plots import (ColdFront, OccludedFront, StationaryFront,
                         StationPlot, WarmFront)
from urllib.request import urlopen
import cartopy.crs as ccrs


def plot_bulletin(ax):
    cat = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/noaaport/text/fronts/catalog.xml')
    bulletin_file = f"https://thredds.ucar.edu/thredds/fileServer/noaaport/text/fronts/{list(cat.datasets.values())[-1]}"
    bulletin_df = parse_wpc_surface_bulletin(urlopen(bulletin_file))

    """Plot a dataframe of surface features on a map."""
    size = 6
    fontsize = 9
    complete_style = {'HIGH': {'color': 'blue', 'fontsize': 12, 'weight': 'bold'},
                      'LOW': {'color': 'red', 'fontsize': 12, 'weight': 'bold'},
                      'WARM': {'linewidth': 1, 'path_effects': [WarmFront(size=size)]},
                      'COLD': {'linewidth': 1, 'path_effects': [ColdFront(size=size)]},
                      'OCFNT': {'linewidth': 1, 'path_effects': [OccludedFront(size=size)]},
                      'STNRY': {'linewidth': 1, 'path_effects': [StationaryFront(size=size)]},
                      'TROF': {'linewidth': 2, 'linestyle': 'dashed',
                               'edgecolor': 'darkorange'}}

    for field in ('HIGH', 'LOW'):
        rows = bulletin_df[bulletin_df.feature == field]
        x, y = zip(*((pt.x, pt.y) for pt in rows.geometry))
        sp = StationPlot(ax, x, y, transform=ccrs.PlateCarree(), clip_on=True)
        texts = sp.plot_text('C', [field[0]] * len(x), **complete_style[field])
        params = sp.plot_parameter('S', rows.strength, **complete_style[field])

    for field in ('WARM', 'COLD', 'STNRY', 'OCFNT', 'TROF'):
        rows = bulletin_df[bulletin_df.feature == field]
        geoms = ax.add_geometries(rows.geometry, crs=ccrs.PlateCarree(), **complete_style[field],
                                  facecolor='none')

    valid_time = bulletin_df['valid'][0].strftime('%H')

    return texts, params, geoms, valid_time