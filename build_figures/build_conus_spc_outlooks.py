print("############\nSCRIPT RUNNING\n############")
import time as comp_time
st = comp_time.time()

import warnings
warnings.filterwarnings("ignore")

import geopandas
from datetime import datetime, timedelta, timezone
import requests
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import calendar
from metpy.plots import USCOUNTIES

import sys
import os

# get script dir
script_dir = os.path.dirname(os.path.abspath(__file__))

# get the parent dir 
project_root = os.path.abspath(os.path.join(script_dir, ".."))

if project_root not in sys.path:
    sys.path.append(project_root)




# Current date (or target date)
target_dt = datetime.now(timezone.utc)

# SPC outlook issue times
valid_times = {
    1: ["0600", "1300", "1630", "2000"],  # Day 1 outlooks
    2: ["0600", "1300", "1630"],          # Day 2 outlooks
    3: ["0730"],                          # Day 3 outlooks
}

# Example for day 1
outlooks = []

for outlook_day in range(1,4):
    outlook_url = None
    for dt_try in [target_dt, target_dt - timedelta(days=1)]:  # try target date, then previous day
        year  = dt_try.strftime("%Y")
        month = dt_try.strftime("%m")
        day   = dt_try.strftime("%d")

        for t in valid_times[outlook_day][::-1]:
            url = f"https://www.spc.noaa.gov/products/outlook/archive/{year}/day{outlook_day}otlk_{year}{month}{day}_{t}_cat.lyr.geojson"
            try:
                # Quick check if file exists
                r = requests.head(url)
                if r.status_code == 200:
                    outlook_url = url
                    print(f"    SPC D{outlook_day} DATA LOADED.....{outlook_url}")
                    break
            except Exception as e:
                continue
        if outlook_url is not None:
            break

    if outlook_url is None:
        raise FileNotFoundError("No valid SPC outlook found for target or previous day.")

    # Read geojson
    outlook = geopandas.read_file(outlook_url)
    outlooks.append(outlook)


for outlook, outlook_day in zip(outlooks, range(1,4)):
    # find date & time objects of the outlooks for pretty titles 
    valid_year = outlook['VALID'][0][0:4]
    valid_month = outlook['VALID'][0][4:6]
    valid_day = outlook['VALID'][0][6:8]
    valid_day = outlook['VALID'][0][6:8]
    valid_hour = outlook['VALID'][0][8:10]
    valid_minute = outlook['VALID'][0][10:12]

    issue_year = outlook['ISSUE'][0][0:4]
    issue_month = outlook['ISSUE'][0][4:6]
    issue_day = outlook['ISSUE'][0][6:8]
    issue_day = outlook['ISSUE'][0][6:8]
    issue_hour = outlook['ISSUE'][0][8:10]
    issue_minute = outlook['ISSUE'][0][10:12]

    issue_date = datetime(int(issue_year), int(issue_month), int(issue_day), int(issue_hour))
    valid_date = datetime(int(valid_year), int(valid_month), int(valid_day), int(valid_hour))

    issue_dayint = issue_date.weekday()  # assigns an int to the local day of the week
    valid_dayint = valid_date.weekday()  # assigns an int to the local day of the week
    dayNames = {
        0: 'MON',
        1: 'TUE',
        2: 'WED',
        3: 'THU',
        4: 'FRI',
        5: 'SAT',
        6: 'SUN'}

    issue_dayName = dayNames[issue_dayint]
    valid_dayName = dayNames[valid_dayint]

    valid_monthName = calendar.month_abbr[int(valid_month)].upper()
    issue_monthName = calendar.month_abbr[int(issue_month)].upper()



    def build_map(extent=[-120, -73, 21, 53], projection=ccrs.LambertConformal(), style='light'):

        fig = plt.figure(figsize=(20, 12))
        fig.set_facecolor('#009946')
        ax = plt.axes(projection=projection)

        # apply the map extent (lat/lon bounding box)
        ax.set_extent(extent)
        # axis aspect ratio
        ax.set_box_aspect(0.7)
        # add map features
        
        if style == 'light':
            color = 'gray'
            alpha = 0.7
        else: 
            color = 'black'
            alpha = 0.8
        ax.add_feature(cfeature.STATES, edgecolor='navy', alpha=0.5, linestyle='-', linewidth=1, zorder=10)
        ax.add_feature(cfeature.LAND, facecolor=color, alpha=alpha, zorder=1)
        ax.add_feature(cfeature.OCEAN, facecolor=color, alpha=alpha+0.2, zorder=0)
        #ax.add_feature(cfeature.BORDERS, color='white', alpha=1, linestyle='-', linewidth=1, zorder=11)
        ax.add_feature(cfeature.COASTLINE, color='navy', alpha=0.5, linestyle='-', linewidth=1, zorder=11)

        # apply tight layout to the figure (keeps things tiddy)
        plt.tight_layout()

        # return the figure axis
        return fig, ax



    fig, ax = build_map(style='light')
    ax.add_feature(USCOUNTIES.with_scale('20m'), alpha=0.1, edgecolor='black', linestyle='-', lw=0.5, zorder=12.1)

    # this plots the outlook polygons, if they exist
    try:
        TSTM = ax.add_geometries(outlook.geometry[0], facecolor=outlook.fill[0], edgecolor='black', linewidth=1, alpha=0.6,
                                 zorder=6, crs=ccrs.PlateCarree(), label="THUNDERSTORM")
        MRGL = ax.add_geometries(outlook.geometry[1], facecolor=outlook.fill[1], edgecolor='black', linewidth=1, alpha=0.7,
                                 zorder=6, crs=ccrs.PlateCarree(), label="MARGINAL")
        SLGT = ax.add_geometries(outlook.geometry[2], facecolor=outlook.fill[2], edgecolor='black', linewidth=1, alpha=0.8,
                                 zorder=6, crs=ccrs.PlateCarree(), label="SLIGHT")
        ENH = ax.add_geometries(outlook.geometry[3], facecolor=outlook.fill[3], edgecolor='black', linewidth=1, alpha=0.9,
                                zorder=6, crs=ccrs.PlateCarree(), label="ENHANCED")
        MDT = ax.add_geometries(outlook.geometry[4], facecolor=outlook.fill[4], edgecolor='black', linewidth=1, alpha=0.9,
                                zorder=6, crs=ccrs.PlateCarree(), label="MODERATE")
        HIGH = ax.add_geometries(outlook.geometry[5], facecolor=outlook.fill[5], edgecolor='black', linewidth=1, alpha=0.9,
                                 zorder=6, crs=ccrs.PlateCarree(), label="HIGH")
    except:
        pass


    # this make dummy-polygons for the legend to look nice
    proxy_TSTM = mpatches.Rectangle((0, 0), 1, 0.1, linewidth=2, edgecolor='black', facecolor='#C1E9C1')
    proxy_MRGL = mpatches.Rectangle((0, 0), 1, 0.1, linewidth=2, edgecolor='black', facecolor='#66A366')
    proxy_SLGT = mpatches.Rectangle((0, 0), 1, 0.1, linewidth=2, edgecolor='black', facecolor='#FFE066')
    proxy_ENH = mpatches.Rectangle((0, 0), 1, 0.1, linewidth=2, edgecolor='black', facecolor='#FFA366')
    proxy_MDT = mpatches.Rectangle((0, 0), 1, 0.1, linewidth=2, edgecolor='black', facecolor='#E06666')
    proxy_HIGH = mpatches.Rectangle((0, 0), 1, 0.1, linewidth=2, edgecolor='black', facecolor='#EE99EE')




    #################################
    # ADD MAP EXTRAS
    #################################
    ax.legend([proxy_TSTM, proxy_MRGL, proxy_SLGT, proxy_ENH, proxy_MDT, proxy_HIGH],
              ['GENERAL THUNDER', '1: MARGINAL', '2: SLIGHT', '3: ENHANCED', '4: MODERATE', '5: HIGH'],
              loc='lower left', fontsize=12, facecolor='white', markerscale=8, framealpha=0.7, labelcolor='k', shadow=True,
              borderpad=0.7, title='SEVERE THUNDERSTORM\n      RISK CATEGORIES', title_fontsize=15).set_zorder(10)

    # plot title, add one to the left with model name and data names, add another to the right with time info
    plt.figtext(0.08, 1.03, f'   STORM PREDICTION CENTER DAY {outlook_day} CONVECTIVE OUTLOOK', weight='bold', ha='left', fontsize=20, color='white')
    plt.figtext(0.08, 1.00, f'    ISSUED: {issue_hour}{issue_minute}z {issue_dayName} {issue_monthName} {issue_day}, {issue_year}  |  VALID: {valid_dayName} {valid_monthName} {valid_day}, {valid_year}', ha='left', fontsize=18, color='white')
    plt.figtext(0.915, 1.04, f' ', ha='left', fontsize=20)

    from PIL import Image
    img = Image.open('utils/images/und-logo.png')
    #                  side-side  up-down  size   size
    imgax = fig.add_axes([0.85, 1, 0.06, 0.06], anchor='SE', zorder=3)
    imgax.imshow(img)
    imgax.axis('off')


    plt.savefig(f"staged_figures/middle_2_2_spc-d{outlook_day}outlook.png", bbox_inches="tight")


elapsed_time = comp_time.time() - st
print(f"############\nSCRIPT FINISHED: time: {comp_time.strftime("%H:%M:%S", comp_time.gmtime(elapsed_time))}\n############")