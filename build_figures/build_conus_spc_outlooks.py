##########################################################
#              SPC OUTLOOK PLOTS SCRIPT
#  (c) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2025
##########################################################


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

import io
import sys
import os

# get script dir
script_dir = os.path.dirname(os.path.abspath(__file__))

# get the parent dir 
project_root = os.path.abspath(os.path.join(script_dir, ".."))

if project_root not in sys.path:
    sys.path.append(project_root)

from utils.utils import *
from utils.map import map_builder
from utils.figure import figure_builder


# Current date (or target date)
target_dt = datetime.now(timezone.utc)

# SPC outlook issue times
valid_times = {
    1: ["0600", "1300", "1630", "2000"],  # Day 1 outlooks
    2: ["0600", "1730"],                  # Day 2 outlooks
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
                r = requests.get(url)
                if r.status_code == 200 and "json" in r.headers.get("Content-Type", "").lower():
                    #print('status is 200')
                    #print(url)
                    outlook_url = url
                    print(f"    SPC D{outlook_day} DATA LOADED.....{outlook_url}")
                    # Read directly from memory buffer
                    outlook = geopandas.read_file(io.BytesIO(r.content))
                    outlooks.append(outlook)
                    break
            except Exception as e:
                continue
        if outlook_url is not None:
            break

    if outlook_url is None:
        raise FileNotFoundError("No valid SPC outlook found for target or previous day.")


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


    fig, ax = map_builder(extent=[-119, -74, 23, 50], terrain=True, counties=True, county_alpha=0.9)

    # this plots the outlook polys
    try:
        TSTM = ax.add_geometries(outlook.geometry[0], facecolor=outlook.fill[0], edgecolor='black', linewidth=1, alpha=0.3,
                                 zorder=6, crs=ccrs.PlateCarree(), label="THUNDERSTORM")
        MRGL = ax.add_geometries(outlook.geometry[1], facecolor=outlook.fill[1], edgecolor='black', linewidth=1, alpha=0.6,
                                 zorder=6, crs=ccrs.PlateCarree(), label="MARGINAL")
        SLGT = ax.add_geometries(outlook.geometry[2], facecolor=outlook.fill[2], edgecolor='black', linewidth=1, alpha=0.7,
                                 zorder=6, crs=ccrs.PlateCarree(), label="SLIGHT")
        ENH = ax.add_geometries(outlook.geometry[3], facecolor=outlook.fill[3], edgecolor='black', linewidth=1, alpha=0.75,
                                zorder=6, crs=ccrs.PlateCarree(), label="ENHANCED")
        MDT = ax.add_geometries(outlook.geometry[4], facecolor=outlook.fill[4], edgecolor='black', linewidth=1, alpha=0.75,
                                zorder=6, crs=ccrs.PlateCarree(), label="MODERATE")
        HIGH = ax.add_geometries(outlook.geometry[5], facecolor=outlook.fill[5], edgecolor='black', linewidth=1, alpha=0.75,
                                 zorder=6, crs=ccrs.PlateCarree(), label="HIGH")
    except:
        pass


    spc_colors = [
        "#C1E9C1",  # TSTM
        "#66A366",  # MRGL
        "#FFE066",  # SLGT
        "#FFA366",  # ENH
        "#E06666",  # MDT
        "#EE99EE",  # HIGH
    ]

    spc_labels = [
        "TSTM",
        "MRGL",
        "SLGT",
        "ENH",
        "MDT",
        "HIGH",
    ]

    outlook_filename = build_filename("staged_figures/conus_spc_outlooks/", f"outlook", issue_date, variant=f"d{outlook_day}")

    figure_builder(fig, ax,
        title=f"Day {outlook_day} Convective Outlook",
        subtitle=f"NOAA Storm Prediction Center",
        valid=f"Issued: {issue_hour}{issue_minute}z {issue_dayName} {issue_monthName} {issue_day}, {issue_year}  •  Valid: {valid_dayName} {valid_monthName} {valid_day}, {valid_year}",
        category_colors=spc_colors,
        category_labels=spc_labels,
        category_title=None,
        footer_left=f"NOAA Storm Prediction Center Convective Outlook • https://www.spc.noaa.gov/products/outlook/day{outlook_day}otlk.html",
        save_path=outlook_filename)


elapsed_time = comp_time.time() - st
print(f"############\nSCRIPT FINISHED: time: {comp_time.strftime("%H:%M:%S", comp_time.gmtime(elapsed_time))}\n############")



    # # this make dummy-polygons for the legend to look nice
    # proxy_TSTM = mpatches.Rectangle((0, 0), 1, 0.1, linewidth=2, edgecolor='black', facecolor='#C1E9C1')
    # proxy_MRGL = mpatches.Rectangle((0, 0), 1, 0.1, linewidth=2, edgecolor='black', facecolor='#66A366')
    # proxy_SLGT = mpatches.Rectangle((0, 0), 1, 0.1, linewidth=2, edgecolor='black', facecolor='#FFE066')
    # proxy_ENH = mpatches.Rectangle((0, 0), 1, 0.1, linewidth=2, edgecolor='black', facecolor='#FFA366')
    # proxy_MDT = mpatches.Rectangle((0, 0), 1, 0.1, linewidth=2, edgecolor='black', facecolor='#E06666')
    # proxy_HIGH = mpatches.Rectangle((0, 0), 1, 0.1, linewidth=2, edgecolor='black', facecolor='#EE99EE')


    # #################################
    # # ADD MAP EXTRAS
    # #################################
    # ax.legend([proxy_TSTM, proxy_MRGL, proxy_SLGT, proxy_ENH, proxy_MDT, proxy_HIGH],
    #           ['GENERAL THUNDER', '1: MARGINAL', '2: SLIGHT', '3: ENHANCED', '4: MODERATE', '5: HIGH'],
    #           loc='lower left', fontsize=12, facecolor='white', markerscale=8, framealpha=0.7, labelcolor='k', shadow=True,
    #           borderpad=0.7, title='SEVERE THUNDERSTORM\n      RISK CATEGORIES', title_fontsize=15).set_zorder(10)

    # # plot title, add one to the left with model name and data names, add another to the right with time info
    # plt.figtext(0.08, 1.03, f'   STORM PREDICTION CENTER DAY {outlook_day} CONVECTIVE OUTLOOK', weight='bold', ha='left', fontsize=20, color='white')
    # plt.figtext(0.08, 1.00, f'    ISSUED: {issue_hour}{issue_minute}z {issue_dayName} {issue_monthName} {issue_day}, {issue_year}  |  VALID: {valid_dayName} {valid_monthName} {valid_day}, {valid_year}', ha='left', fontsize=18, color='white')
    # plt.figtext(0.915, 1.04, f' ', ha='left', fontsize=20)

    # from PIL import Image
    # img = Image.open('utils/images/und-logo.png')
    # #                  side-side  up-down  size   size
    # imgax = fig.add_axes([0.85, 1, 0.06, 0.06], anchor='SE', zorder=3)
    # imgax.imshow(img)
    # imgax.axis('off')
