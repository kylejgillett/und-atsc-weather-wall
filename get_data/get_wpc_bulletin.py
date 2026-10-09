##########################################################
#              WPC FRONTS RETRIEVAL SCRIPT
#  (c) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2025
##########################################################

from siphon.catalog import TDSCatalog
from metpy.io import parse_wpc_surface_bulletin
from metpy.plots import (ColdFront, OccludedFront, StationaryFront, StationPlot, WarmFront)
from urllib.request import urlopen
import cartopy.crs as ccrs
import pandas as pd
import numpy as np


def plot_bulletin(ax):

    # Download and parse the latest bulletin
    try:
        cat = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/noaaport/text/fronts/catalog.xml')
        bulletin_file = (f'https://thredds.ucar.edu/thredds/fileServer/noaaport/text/fronts/{list(cat.datasets.values())[-1]}')

        with urlopen(bulletin_file, timeout=20) as response:
            bulletin_df = parse_wpc_surface_bulletin(response)
        if bulletin_df.empty:
            raise ValueError("WPC bulletin contains no features")
    except Exception as e:
        print(f"    WPC BULLETIN UNAVAILABLE.....{e}")
        return [], [], [], "N/A"

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
    

    # Track artists so they can be removed on failure
    texts = []
    params = []
    geoms = []
    created_artists = []

    try:
        # Plot HIGH and LOW centers
        for field in ('HIGH', 'LOW'):
            rows = bulletin_df[bulletin_df.feature == field].copy()
            if rows.empty:
                continue

            # Extract geometry coordinates
            rows['lon'] = rows.geometry.apply(lambda geom: geom.x if geom is not None
                and geom.geom_type == 'Point' and not geom.is_empty else np.nan)
            rows['lat'] = rows.geometry.apply(lambda geom: geom.y if geom is not None
                and geom.geom_type == 'Point' and not geom.is_empty else np.nan)
            
            # Convert pressure strengths to numeric
            rows['strength'] = pd.to_numeric(rows['strength'], errors='coerce')

            # Remove invalid geographic coordinates
            valid_coords = (np.isfinite(rows['lon']) & np.isfinite(rows['lat']) & rows['lon'].between(-180, 180) & rows['lat'].between(-90, 90))
            rows = rows.loc[valid_coords].copy()
            if rows.empty:
                continue
            # Validate the map projection
            xy = ax.projection.transform_points(ccrs.PlateCarree(), rows['lon'].to_numpy(dtype=float), rows['lat'].to_numpy(dtype=float))

            valid_xy = (np.isfinite(xy[:, 0]) & np.isfinite(xy[:, 1]))
            rows = rows.loc[valid_xy].reset_index(drop=True)
            if rows.empty:
                continue

            # Projected coordinates, already validated
            xy = xy[valid_xy]
            sp = StationPlot(ax, xy[:, 0], xy[:, 1], transform=ax.projection, clip_on=True)
            text_artist = sp.plot_text('C', [field[0]] * len(rows),  **complete_style[field])
            texts.append(text_artist)
            created_artists.append(text_artist)

            # Plot pressure only where pressure is finite
            pressure_valid = np.isfinite(rows['strength'].to_numpy())
            if pressure_valid.any():
                pressure_sp = StationPlot(ax, xy[pressure_valid, 0], xy[pressure_valid, 1], transform=ax.projection, clip_on=True)
                param_artist = pressure_sp.plot_parameter( 'S', rows.loc[pressure_valid, 'strength'].to_numpy(), **complete_style[field])
                params.append(param_artist)
                created_artists.append(param_artist)

        # Plot fronts and troughs
        for field in ('WARM', 'COLD', 'STNRY', 'OCFNT', 'TROF'):
            rows = bulletin_df[bulletin_df.feature == field]
            valid_geometries = [geom for geom in rows.geometry if geom is not None and not geom.is_empty]
            if not valid_geometries:
                continue

            geom_artist = ax.add_geometries( valid_geometries, crs=ccrs.PlateCarree(), **complete_style[field], facecolor='none')
            geoms.append(geom_artist)
            created_artists.append(geom_artist)

        # Extract bulletin valid time
        valid_dates = pd.to_datetime(bulletin_df['valid'], errors='coerce', utc=True).dropna()

        valid_time = (valid_dates.iloc[0].strftime('%H') if not valid_dates.empty else "N/A")

        print(f"    WPC BULLETIN LOADED.....{len(texts)} PRESSURE-CENTER LABEL GROUPS, {len(geoms)} FRONT GROUPS")

        return texts, params, geoms, valid_time

    except Exception as e:
        print(f"    WPC BULLETIN PLOTTING FAILED.....{e}")
        for artist in created_artists:
            if artist is not None:
                artist.remove()
        return [], [], [], "N/A"



# from siphon.catalog import TDSCatalog
# from metpy.io import parse_wpc_surface_bulletin
# from metpy.plots import (ColdFront, OccludedFront, StationaryFront,
#                          StationPlot, WarmFront)
# from urllib.request import urlopen
# import cartopy.crs as ccrs


# def plot_bulletin(ax):
#     cat = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/noaaport/text/fronts/catalog.xml')
#     bulletin_file = f"https://thredds.ucar.edu/thredds/fileServer/noaaport/text/fronts/{list(cat.datasets.values())[-1]}"
#     bulletin_df = parse_wpc_surface_bulletin(urlopen(bulletin_file))

#     """Plot a dataframe of surface features on a map."""
#     size = 6
#     fontsize = 9
#     complete_style = {'HIGH': {'color': 'blue', 'fontsize': 12, 'weight': 'bold'},
#                       'LOW': {'color': 'red', 'fontsize': 12, 'weight': 'bold'},
#                       'WARM': {'linewidth': 1, 'path_effects': [WarmFront(size=size)]},
#                       'COLD': {'linewidth': 1, 'path_effects': [ColdFront(size=size)]},
#                       'OCFNT': {'linewidth': 1, 'path_effects': [OccludedFront(size=size)]},
#                       'STNRY': {'linewidth': 1, 'path_effects': [StationaryFront(size=size)]},
#                       'TROF': {'linewidth': 2, 'linestyle': 'dashed',
#                                'edgecolor': 'darkorange'}}

#     for field in ('HIGH', 'LOW'):
#         rows = bulletin_df[bulletin_df.feature == field]
#         x, y = zip(*((pt.x, pt.y) for pt in rows.geometry))
#         sp = StationPlot(ax, x, y, transform=ccrs.PlateCarree(), clip_on=True)
#         texts = sp.plot_text('C', [field[0]] * len(x), **complete_style[field])
#         params = sp.plot_parameter('S', rows.strength, **complete_style[field])

#     for field in ('WARM', 'COLD', 'STNRY', 'OCFNT', 'TROF'):
#         rows = bulletin_df[bulletin_df.feature == field]
#         geoms = ax.add_geometries(rows.geometry, crs=ccrs.PlateCarree(), **complete_style[field],
#                                   facecolor='none')

#     valid_time = bulletin_df['valid'][0].strftime('%H')

#     return texts, params, geoms, valid_time