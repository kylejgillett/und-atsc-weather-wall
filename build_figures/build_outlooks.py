##########################################################
#              NOAA / NWS OUTLOOK PLOTS
#  UND Atmospheric Sciences Weather Wall
#  KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2026
##########################################################

print("############\nSCRIPT RUNNING\n############")
import time as comp_time
st = comp_time.time()

import warnings
warnings.filterwarnings("ignore")
import os
import sys
from datetime import datetime, timezone
import cartopy.crs as ccrs
import matplotlib.pyplot as plt


# get script dir
script_dir = os.path.dirname(os.path.abspath(__file__))
# get the parent dir
project_root = os.path.abspath(os.path.join(script_dir, ".."))
if project_root not in sys.path:
    sys.path.append(project_root)
# import modules from sub dirs
from utils.utils import *
from utils.map import map_builder
from utils.figure import figure_builder
from get_data.get_outlooks import *


now_utc = datetime.now(timezone.utc)

CONUS_EXTENT = [-119, -74, 23.5, 53.5]
TROPICAL_EXTENT = [-140, -10, 5, 47]
#TROPICAL_PROJECTION = ccrs.Mercator(central_longitude=-75, min_latitude=0, max_latitude=55)
OUTPUT_DIR = "staged_figures/conus_outlooks/"
os.makedirs(OUTPUT_DIR, exist_ok=True)


# legend and color info
SPC_COLORS = {2: "#C1E9C1",  3: "#66A366",   4: "#FFE066",  5: "#FFA366",  6: "#E06666",  8: "#EE99EE"}
SPC_DN = [2, 3, 4, 5, 6, 8]
SPC_LABELS = ["TSTM", "MRGL", "SLGT", "ENH", "MDT", "HIGH"]

FIRE_COLORS = {("wind_rh", 5): "#E69800", ("wind_rh", 8): "#FF0000", ("wind_rh", 10): "#E600A9", 
               ("dry_thunder", 5): "#732600", ("dry_thunder", 8): "#B30000"}
FIRE_LEGEND = [ ("#E69800", "ELV"),  ("#FF0000", "CRT"), ("#E600A9", "EXT"), ("#732600", "IsoDT"), ("#B30000", "SctDT")]

ERO_COLORS = {1: "#38A800", 2: "#FFFE00", 3: "#F50000", 4: "#FF69C5"}
ERO_DN = [1, 2, 3, 4]
ERO_LABELS = ["MGRL\n5%+", "SLGT\n15%+", "MDT\n40%+", "HIGHT\n70%+"]

CPC_TEMP_LEGEND = [
    (("Above", 90), "#702100", "Abv\n90%"),
    (("Above", 80), "#912600", "Abv\n80%"),
    (("Above", 70), "#B32E05", "Abv\n70%"),
    (("Above", 60), "#C93B1A", "Abv\n60%"),
    (("Above", 50), "#DA5731", "Abv\n50%"),
    (("Above", 40), "#E38B4B", "Abv\n40%"),
    (("Above", 33), "#E7B168", "Abv\n33%"),
    (("Normal", 36), "#A0A0A0", "Nrml"),
    (("Below", 33), "#BFCBE4", "Blw\n33%"),
    (("Below", 40), "#A0C0DF", "Blw\n40%"),
    (("Below", 50), "#77B5E2", "Blw\n50%"),
    (("Below", 60), "#389FDC", "Blw\n60%"),
    (("Below", 70), "#005DA1", "Blw\n70%"),
    (("Below", 80), "#2E216F", "Blw\n80%"),
    (("Below", 90), "#221852", "Blw\n90%")]

CPC_PRCP_LEGEND = [
    (("Above", 90), "#285300", "Abv\n90%"),
    (("Above", 80), "#28600A", "Abv\n80%"),
    (("Above", 70), "#007814", "Abv\n70%"),
    (("Above", 60), "#009620", "Abv\n60%"),
    (("Above", 50), "#48B430", "Abv\n50%"),
    (("Above", 40), "#95CE7F", "Abv\n40%"),
    (("Above", 33), "#B3D9AB", "Abv\n33%"),
    (("Normal", 36), "#A0A0A0", "Nrml"),
    (("Below", 33), "#F0D493", "Blw\n33%"),
    (("Below", 40), "#D8A74F", "Blw\n40%"),
    (("Below", 50), "#BB6D33", "Blw\n50%"),
    (("Below", 60), "#9B5031", "Blw\n60%"),
    (("Below", 70), "#934639", "Blw\n70%"),
    (("Below", 80), "#804000", "Blw\n80%"),
    (("Below", 90), "#4F2F2F", "Blw\n90%")]

WSO_COLORS = {"10%": "#04BCCA", "30%": "#F7FF03", "50%": "#FF0000", "80%": "#990099"}
WSO_LABELS = ["10-\n30%", "30-\n50%", "50-\n80%", "80%+"]

NHC_COLORS = {"low": "#FFFF00", "medium": "#E69800", "high": "#E60000",}
NHC_STORM_COLORS = {"D": "#FFFF00",  "S": "#00BFFF",  "H": "#FF0000",  "M": "#FF00FF",}
NHC_STORM_LABELS = ["TD", "TS", "H", "MH"]

DROUGHT_COLORS = {0: "#FFFF00", 1: "#FCD37F", 2: "#FFAA00", 3: "#E60000", 4: "#730000"}
DROUGHT_LABELS = ["D0\nDRY", "D1\nMDT", "D2\nSVR", "D3\nEXT", "D4\nEXC"]





#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
#################################
# BUILD SPC DAY 1 CONVECTIVE OUTLOOK
#################################
outlook = get_spc_convective_outlook(day=1)

if outlook is not None:
    fig, ax = map_builder(extent=CONUS_EXTENT, terrain=True, counties=True, county_alpha=0.9)

    for dn in SPC_DN:
        data = outlook[outlook["dn"] == dn]
        if data.empty:
            continue

        ax.add_geometries(data.geometry, crs=ccrs.PlateCarree(), facecolor=SPC_COLORS[dn], edgecolor="black", 
                          linewidth=1.0, alpha=0.35 if dn == 2 else 0.75, zorder=6)

    issue_time = outlook.attrs.get("issue_time") or now_utc
    valid_start = outlook.attrs.get("valid_start")
    valid_end = outlook.attrs.get("valid_end")

    if valid_start is not None and valid_end is not None:
        valid_text = f'Issued: {issue_time.strftime("%H%MZ %a %b %d, %Y").upper()}  •  Valid: {valid_start.strftime("%H%MZ %a %b %d").upper()} - {valid_end.strftime("%H%MZ %a %b %d, %Y").upper()}'
    else:
        valid_text = f'Issued: {issue_time.strftime("%H%MZ %a %b %d, %Y").upper()}'

    outlook_filename = build_filename(OUTPUT_DIR, "outlook", now_utc, variant="00-spc-conv-d1")

    figure_builder(fig, ax,
        title="Day 1 Convective Outlook",
        subtitle="NOAA Storm Prediction Center",
        valid=valid_text,
        category_colors=[SPC_COLORS[dn] for dn in SPC_DN],
        category_labels=SPC_LABELS,
        category_title=None,
        footer_left="NOAA Storm Prediction Center Convective Outlook • https://www.spc.noaa.gov/products/outlook/",
        save_path=outlook_filename)

    plt.close(fig)
    print("    FINISHED SPC DAY 1 CONVECTIVE OUTLOOK")
else:
    print("    SKIPPING SPC DAY 1 CONVECTIVE OUTLOOK")
    
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################



#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
#################################
# BUILD SPC DAY 2 CONVECTIVE OUTLOOK
#################################
outlook = get_spc_convective_outlook(day=2)

if outlook is not None:
    fig, ax = map_builder(extent=CONUS_EXTENT, terrain=True, counties=True, county_alpha=0.9)

    for dn in SPC_DN:
        data = outlook[outlook["dn"] == dn]
        if data.empty:
            continue

        ax.add_geometries(data.geometry, crs=ccrs.PlateCarree(), facecolor=SPC_COLORS[dn], edgecolor="black", linewidth=1.0, alpha=0.35 if dn == 2 else 0.75, zorder=6)

    issue_time = outlook.attrs.get("issue_time") or now_utc
    valid_start = outlook.attrs.get("valid_start")
    valid_end = outlook.attrs.get("valid_end")

    if valid_start is not None and valid_end is not None:
        valid_text = f'Issued: {issue_time.strftime("%H%MZ %a %b %d, %Y").upper()}  •  Valid: {valid_start.strftime("%H%MZ %a %b %d").upper()} - {valid_end.strftime("%H%MZ %a %b %d, %Y").upper()}'
    else:
        valid_text = f'Issued: {issue_time.strftime("%H%MZ %a %b %d, %Y").upper()}'

    outlook_filename = build_filename(OUTPUT_DIR, "outlook", now_utc, variant="01-spc-conv-d2")

    figure_builder(fig, ax,
        title="Day 2 Convective Outlook",
        subtitle="NOAA Storm Prediction Center",
        valid=valid_text,
        category_colors=[SPC_COLORS[dn] for dn in SPC_DN],
        category_labels=SPC_LABELS,
        category_title=None,
        footer_left="NOAA Storm Prediction Center Convective Outlook • https://www.spc.noaa.gov/products/outlook/",
        save_path=outlook_filename)

    plt.close(fig)
    print("    FINISHED SPC DAY 2 CONVECTIVE OUTLOOK")
else:
    print("    SKIPPING SPC DAY 2 CONVECTIVE OUTLOOK")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################



#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
#################################
# BUILD SPC DAY 3 CONVECTIVE OUTLOOK
#################################
outlook = get_spc_convective_outlook(day=3)

if outlook is not None:
    fig, ax = map_builder(extent=CONUS_EXTENT, terrain=True, counties=True, county_alpha=0.9)

    for dn in SPC_DN:
        data = outlook[outlook["dn"] == dn]
        if data.empty:
            continue

        ax.add_geometries(data.geometry, crs=ccrs.PlateCarree(), facecolor=SPC_COLORS[dn], edgecolor="black", linewidth=1.0, alpha=0.35 if dn == 2 else 0.75, zorder=6)

    issue_time = outlook.attrs.get("issue_time") or now_utc
    valid_start = outlook.attrs.get("valid_start")
    valid_end = outlook.attrs.get("valid_end")

    if valid_start is not None and valid_end is not None:
        valid_text = f'Issued: {issue_time.strftime("%H%MZ %a %b %d, %Y").upper()}  •  Valid: {valid_start.strftime("%H%MZ %a %b %d").upper()} - {valid_end.strftime("%H%MZ %a %b %d, %Y").upper()}'
    else:
        valid_text = f'Issued: {issue_time.strftime("%H%MZ %a %b %d, %Y").upper()}'

    outlook_filename = build_filename(OUTPUT_DIR, "outlook", now_utc, variant="02-spc-conv-d3")

    figure_builder(fig, ax,
        title="Day 3 Convective Outlook",
        subtitle="NOAA Storm Prediction Center",
        valid=valid_text,
        category_colors=[SPC_COLORS[dn] for dn in SPC_DN],
        category_labels=SPC_LABELS,
        category_title=None,
        footer_left="NOAA Storm Prediction Center Convective Outlook • https://www.spc.noaa.gov/products/outlook/",
        save_path=outlook_filename)

    plt.close(fig)
    print("    FINISHED SPC DAY 3 CONVECTIVE OUTLOOK")
else:
    print("    SKIPPING SPC DAY 3 CONVECTIVE OUTLOOK")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################



#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
#################################
# BUILD SPC DAY 1 FIRE WEATHER OUTLOOK
#################################
outlook = get_spc_fire_outlook(day=1)

if outlook is not None:
    fig, ax = map_builder(extent=CONUS_EXTENT, terrain=True, counties=True, county_alpha=0.9)

    for risk_type, dn in [("wind_rh", 5), ("wind_rh", 8), ("wind_rh", 10), ("dry_thunder", 5), ("dry_thunder", 8)]:
        data = outlook[(outlook["risk_type"] == risk_type) & (outlook["dn"] == dn)]
        if data.empty:
            continue

        ax.add_geometries(data.geometry, crs=ccrs.PlateCarree(), facecolor=FIRE_COLORS[(risk_type, dn)],
                           edgecolor="black", linewidth=1.0, alpha=0.70, zorder=6 if risk_type == "wind_rh" else 7)

    issue_time = outlook.attrs.get("issue_time") or now_utc
    valid_start = outlook.attrs.get("valid_start")
    valid_end = outlook.attrs.get("valid_end")

    if valid_start is not None and valid_end is not None:
        valid_text = f'Issued: {issue_time.strftime("%H%MZ %a %b %d, %Y").upper()}  •  Valid: {valid_start.strftime("%H%MZ %a %b %d").upper()} - {valid_end.strftime("%H%MZ %a %b %d, %Y").upper()}'
    else:
        valid_text = f'Issued: {issue_time.strftime("%H%MZ %a %b %d, %Y").upper()}'

    outlook_filename = build_filename(OUTPUT_DIR, "outlook", now_utc, variant="03-spc-fire-d1")

    figure_builder(fig, ax,
        title="Day 1 Fire Weather Outlook",
        subtitle="NOAA Storm Prediction Center",
        valid=valid_text,
        category_colors=[item[0] for item in FIRE_LEGEND],
        category_labels=[item[1] for item in FIRE_LEGEND],
        category_title=None,
        footer_left="NOAA Storm Prediction Center Fire Weather Outlook • https://www.spc.noaa.gov/products/fire_wx/",
        save_path=outlook_filename)

    plt.close(fig)
    print("    FINISHED SPC DAY 1 FIRE WEATHER OUTLOOK")
else:
    print("    SKIPPING SPC DAY 1 FIRE WEATHER OUTLOOK")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################



#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
#################################
# BUILD SPC DAY 2 FIRE WEATHER OUTLOOK
#################################
outlook = get_spc_fire_outlook(day=2)

if outlook is not None:
    fig, ax = map_builder(extent=CONUS_EXTENT, terrain=True, counties=True, county_alpha=0.9)

    for risk_type, dn in [("wind_rh", 5), ("wind_rh", 8), ("wind_rh", 10), ("dry_thunder", 5), ("dry_thunder", 8)]:
        data = outlook[(outlook["risk_type"] == risk_type) & (outlook["dn"] == dn)]
        if data.empty:
            continue

        ax.add_geometries(data.geometry, crs=ccrs.PlateCarree(), facecolor=FIRE_COLORS[(risk_type, dn)],
                           edgecolor="black", linewidth=1.0, alpha=0.70, zorder=6 if risk_type == "wind_rh" else 7)

    issue_time = outlook.attrs.get("issue_time") or now_utc
    valid_start = outlook.attrs.get("valid_start")
    valid_end = outlook.attrs.get("valid_end")

    if valid_start is not None and valid_end is not None:
        valid_text = f'Issued: {issue_time.strftime("%H%MZ %a %b %d, %Y").upper()}  •  Valid: {valid_start.strftime("%H%MZ %a %b %d").upper()} - {valid_end.strftime("%H%MZ %a %b %d, %Y").upper()}'
    else:
        valid_text = f'Issued: {issue_time.strftime("%H%MZ %a %b %d, %Y").upper()}'

    outlook_filename = build_filename(OUTPUT_DIR, "outlook", now_utc, variant="04-spc-fire-d2")

    figure_builder(fig, ax,
        title="Day 2 Fire Weather Outlook",
        subtitle="NOAA Storm Prediction Center",
        valid=valid_text,
        category_colors=[item[0] for item in FIRE_LEGEND],
        category_labels=[item[1] for item in FIRE_LEGEND],
        category_title=None,
        footer_left="NOAA Storm Prediction Center Fire Weather Outlook • https://www.spc.noaa.gov/products/fire_wx/",
        save_path=outlook_filename)

    plt.close(fig)
    print("    FINISHED SPC DAY 2 FIRE WEATHER OUTLOOK")
else:
    print("    SKIPPING SPC DAY 2 FIRE WEATHER OUTLOOK")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################



#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
#################################
# BUILD WPC DAY 1 EXCESSIVE RAINFALL OUTLOOK
#################################
outlook = get_wpc_excessive_rainfall_outlook(day=1)

if outlook is not None:
    fig, ax = map_builder(extent=CONUS_EXTENT, terrain=True, counties=True, county_alpha=0.9)

    for dn in ERO_DN:
        data = outlook[outlook["dn"] == dn]
        if data.empty:
            continue

        ax.add_geometries(data.geometry, crs=ccrs.PlateCarree(), facecolor=ERO_COLORS[dn], edgecolor="black", linewidth=1.0, alpha=0.75, zorder=6)

    issue_time = outlook.attrs.get("issue_time") or now_utc
    valid_start = outlook.attrs.get("valid_start")
    valid_end = outlook.attrs.get("valid_end")

    if valid_start is not None and valid_end is not None:
        valid_text = f'Issued: {issue_time.strftime("%H%MZ %a %b %d, %Y").upper()}  •  Valid: {valid_start.strftime("%H%MZ %a %b %d").upper()} - {valid_end.strftime("%H%MZ %a %b %d, %Y").upper()}'
    else:
        valid_text = f'Issued: {issue_time.strftime("%H%MZ %a %b %d, %Y").upper()}'

    outlook_filename = build_filename(OUTPUT_DIR, "outlook", now_utc, variant="05-wpc-ero-d1")

    figure_builder(fig, ax,
        title="Day 1 Excessive Rainfall Outlook",
        subtitle="NOAA Weather Prediction Center",
        valid=valid_text,
        category_colors=[ERO_COLORS[dn] for dn in ERO_DN],
        category_labels=ERO_LABELS,
        category_title=None,
        footer_left="NOAA Weather Prediction Center Excessive Rainfall Outlook • https://www.wpc.ncep.noaa.gov/#page=ero",
        save_path=outlook_filename)

    plt.close(fig)
    print("    FINISHED WPC DAY 1 EXCESSIVE RAINFALL OUTLOOK")
else:
    print("    SKIPPING WPC DAY 1 EXCESSIVE RAINFALL OUTLOOK")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################



#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
#################################
# BUILD WPC DAY 2 EXCESSIVE RAINFALL OUTLOOK
#################################
outlook = get_wpc_excessive_rainfall_outlook(day=2)

if outlook is not None:
    fig, ax = map_builder(extent=CONUS_EXTENT, terrain=True, counties=True, county_alpha=0.9)

    for dn in ERO_DN:
        data = outlook[outlook["dn"] == dn]
        if data.empty:
            continue

        ax.add_geometries(data.geometry, crs=ccrs.PlateCarree(), facecolor=ERO_COLORS[dn], edgecolor="black", linewidth=1.0, alpha=0.75, zorder=6)

    issue_time = outlook.attrs.get("issue_time") or now_utc
    valid_start = outlook.attrs.get("valid_start")
    valid_end = outlook.attrs.get("valid_end")

    if valid_start is not None and valid_end is not None:
        valid_text = f'Issued: {issue_time.strftime("%H%MZ %a %b %d, %Y").upper()}  •  Valid: {valid_start.strftime("%H%MZ %a %b %d").upper()} - {valid_end.strftime("%H%MZ %a %b %d, %Y").upper()}'
    else:
        valid_text = f'Issued: {issue_time.strftime("%H%MZ %a %b %d, %Y").upper()}'

    outlook_filename = build_filename(OUTPUT_DIR, "outlook", now_utc, variant="06-wpc-ero-d2")

    figure_builder(fig, ax,
        title="Day 2 Excessive Rainfall Outlook",
        subtitle="NOAA Weather Prediction Center",
        valid=valid_text,
        category_colors=[ERO_COLORS[dn] for dn in ERO_DN],
        category_labels=ERO_LABELS,
        category_title=None,
        footer_left="NOAA Weather Prediction Center Excessive Rainfall Outlook • https://www.wpc.ncep.noaa.gov/#page=ero",
        save_path=outlook_filename)

    plt.close(fig)
    print("    FINISHED WPC DAY 2 EXCESSIVE RAINFALL OUTLOOK")
else:
    print("    SKIPPING WPC DAY 2 EXCESSIVE RAINFALL OUTLOOK")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################






#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
#################################
# BUILD WPC DAY 3 EXCESSIVE RAINFALL OUTLOOK
#################################
outlook = get_wpc_excessive_rainfall_outlook(day=3)

if outlook is not None:
    fig, ax = map_builder(extent=CONUS_EXTENT, terrain=True, counties=True, county_alpha=0.9)

    for dn in ERO_DN:
        data = outlook[outlook["dn"] == dn]
        if data.empty:
            continue

        ax.add_geometries(data.geometry, crs=ccrs.PlateCarree(), facecolor=ERO_COLORS[dn],
                           edgecolor="black", linewidth=1.0, alpha=0.75, zorder=6)

    issue_time = outlook.attrs.get("issue_time") or now_utc
    valid_start = outlook.attrs.get("valid_start")
    valid_end = outlook.attrs.get("valid_end")

    if valid_start is not None and valid_end is not None:
        valid_text = f'Issued: {issue_time.strftime("%H%MZ %a %b %d, %Y").upper()}  •  Valid: {valid_start.strftime("%H%MZ %a %b %d").upper()} - {valid_end.strftime("%H%MZ %a %b %d, %Y").upper()}'
    else:
        valid_text = f'Issued: {issue_time.strftime("%H%MZ %a %b %d, %Y").upper()}'

    outlook_filename = build_filename(OUTPUT_DIR, "outlook", now_utc, variant="07-wpc-ero-d3")

    figure_builder(fig, ax,
        title="Day 3 Excessive Rainfall Outlook",
        subtitle="NOAA Weather Prediction Center",
        valid=valid_text,
        category_colors=[ERO_COLORS[dn] for dn in ERO_DN],
        category_labels=ERO_LABELS,
        category_title=None,
        footer_left="NOAA Weather Prediction Center Excessive Rainfall Outlook • https://www.wpc.ncep.noaa.gov/#page=ero",
        save_path=outlook_filename)

    plt.close(fig)
    print("    FINISHED WPC DAY 3 EXCESSIVE RAINFALL OUTLOOK")
else:
    print("    SKIPPING WPC DAY 3 EXCESSIVE RAINFALL OUTLOOK")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################







#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
#################################
# BUILD CPC 6-10 DAY TEMPERATURE OUTLOOK
#################################
outlook = get_cpc_610_temperature_outlook()

if outlook is not None:
    fig, ax = map_builder(extent=CONUS_EXTENT, terrain=True, counties=True, county_alpha=0.9)

    for (cat, prob), color, label in CPC_TEMP_LEGEND:
        data = outlook[(outlook["cat"].astype(str).str.lower() == cat.lower()) & (outlook["prob"].astype(float).round().astype(int) == prob)]

        if data.empty:
            continue

        ax.add_geometries(data.geometry, crs=ccrs.PlateCarree(), facecolor=color,
                           edgecolor="#777777", linewidth=0.7, alpha=0.76, zorder=6)

    issue_time = outlook.attrs.get("issue_time") or now_utc
    valid_start = outlook.attrs.get("valid_start")
    valid_end = outlook.attrs.get("valid_end")

    if valid_start is not None and valid_end is not None:
        valid_text = f'Issued: {issue_time.strftime("%a %b %d, %Y").upper()}  •  Valid: {valid_start.strftime("%a %b %d").upper()} - {valid_end.strftime("%a %b %d, %Y").upper()}'
    else:
        valid_text = f'Issued: {issue_time.strftime("%a %b %d, %Y").upper()}'

    outlook_filename = build_filename(OUTPUT_DIR, "outlook", now_utc, variant="08-cpc-610-temp")

    figure_builder(fig, ax,
        title="6-10 Day Temperature Outlook",
        subtitle="NOAA Climate Prediction Center",
        valid=valid_text,
        category_colors=[item[1] for item in CPC_TEMP_LEGEND],
        category_labels=[item[2] for item in CPC_TEMP_LEGEND],
        category_title=None,
        footer_left="NOAA Climate Prediction Center Outlook • https://www.cpc.ncep.noaa.gov/",
        save_path=outlook_filename)

    plt.close(fig)
    print("    FINISHED CPC 6-10 DAY TEMPERATURE OUTLOOK")
else:
    print("    SKIPPING CPC 6-10 DAY TEMPERATURE OUTLOOK")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################







#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
#################################
# BUILD CPC 6-10 DAY PRECIPITATION OUTLOOK
#################################
outlook = get_cpc_610_precipitation_outlook()

if outlook is not None:
    fig, ax = map_builder(extent=CONUS_EXTENT, terrain=True, counties=True, county_alpha=0.9)

    for (cat, prob), color, label in CPC_PRCP_LEGEND:
        data = outlook[(outlook["cat"].astype(str).str.lower() == cat.lower()) & (outlook["prob"].astype(float).round().astype(int) == prob)]

        if data.empty:
            continue

        ax.add_geometries(data.geometry, crs=ccrs.PlateCarree(),
                           facecolor=color, edgecolor="#777777", linewidth=0.7, alpha=0.76, zorder=6)

    issue_time = outlook.attrs.get("issue_time") or now_utc
    valid_start = outlook.attrs.get("valid_start")
    valid_end = outlook.attrs.get("valid_end")

    if valid_start is not None and valid_end is not None:
        valid_text = f'Issued: {issue_time.strftime("%a %b %d, %Y").upper()}  •  Valid: {valid_start.strftime("%a %b %d").upper()} - {valid_end.strftime("%a %b %d, %Y").upper()}'
    else:
        valid_text = f'Issued: {issue_time.strftime("%a %b %d, %Y").upper()}'

    outlook_filename = build_filename(OUTPUT_DIR, "outlook", now_utc, variant="09-cpc-610-precip")

    figure_builder(fig, ax,
        title="6-10 Day Precipitation Outlook",
        subtitle="NOAA Climate Prediction Center",
        valid=valid_text,
        category_colors=[item[1] for item in CPC_PRCP_LEGEND],
        category_labels=[item[2] for item in CPC_PRCP_LEGEND],
        category_title=None,
        footer_left="NOAA Climate Prediction Center Outlook • https://www.cpc.ncep.noaa.gov/",
        save_path=outlook_filename)

    plt.close(fig)
    print("    FINISHED CPC 6-10 DAY PRECIPITATION OUTLOOK")
else:
    print("    SKIPPING CPC 6-10 DAY PRECIPITATION OUTLOOK")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################






TROPICAL_PROJECTION = ccrs.Mercator(
    central_longitude=-75,
    min_latitude=-15,
    max_latitude=65,
)
TROPICAL_EXTENT = [-150, -10, -10, 55]

#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
#################################
# BUILD NHC 7-DAY TROPICAL WEATHER OUTLOOK + ACTIVE STORMS
#################################
outlook = get_nhc_7day_tropical_outlook()

if outlook is not None:
    fig, ax = map_builder(extent=TROPICAL_EXTENT, projection=TROPICAL_PROJECTION, terrain=True, terrain_zoom=3, counties=False, map_scale="50m")

    regions = outlook.get("regions")
    motion = outlook.get("motion")
    points = outlook.get("points")

    storm_cone = outlook.get("storm_cone")
    storm_track = outlook.get("storm_track")
    storm_points = outlook.get("storm_points")

    # active storms 
    if storm_cone is not None:
        for _, row in storm_cone.iterrows():
            ax.add_geometries([row.geometry], crs=ccrs.PlateCarree(), facecolor="#E8E8E8",
                               edgecolor="#202020", linewidth=1.2, alpha=0.35, zorder=4)

    if storm_track is not None:
        for _, row in storm_track.iterrows():
            ax.add_geometries([row.geometry], crs=ccrs.PlateCarree(), facecolor="none",
                               edgecolor="black", linewidth=2.0, zorder=8)

    if storm_points is not None:
        for _, row in storm_points.iterrows():
            try:
                forecast_hour = float(row.get("fcstprd", 999.0))
            except Exception:
                forecast_hour = 999.0

            storm_type = str(row.get("dvlbl", "")).strip().upper()
            storm_color = NHC_STORM_COLORS.get(storm_type, "#FFFFFF")

            # fcst locations
            ax.scatter(row.geometry.x, row.geometry.y, s=24 if forecast_hour != 0 else 95, marker="o", c=storm_color, edgecolors="black",
                        linewidths=0.8 if forecast_hour != 0 else 1.4, transform=ccrs.PlateCarree(), zorder=10 if forecast_hour != 0 else 12, clip_on=True)

            # label only current storm
            if forecast_hour == 0:
                storm_name = str(row.get("stormname", "")).strip().upper()

                try:
                    max_wind = int(float(row.get("maxwind")))
                    storm_label = f"{storm_name}  {max_wind} kt"
                except Exception:
                    storm_label = storm_name

                ax.text(row.geometry.x + 0.8, row.geometry.y + 0.5, storm_label, fontsize=10, 
                        fontweight="bold", color="black", transform=ccrs.PlateCarree(), zorder=15, clip_on=True)

    # 7 day outlook
    if regions is not None:
        for _, row in regions.iterrows():
            try:
                probability = int(str(row["prob7day"]).replace("%", "").strip())
            except Exception:
                continue

            if probability <= 30:
                color = NHC_COLORS["low"]
            elif probability <= 60:
                color = NHC_COLORS["medium"]
            else:
                color = NHC_COLORS["high"]

            ax.add_geometries([row.geometry], crs=ccrs.PlateCarree(), facecolor=color, 
                              edgecolor=color, linewidth=1.5, alpha=0.28, hatch="///", zorder=5)

    if motion is not None:
        for _, row in motion.iterrows():
            try:
                probability = int(str(row["prob7day"]).replace("%", "").strip())
            except Exception:
                continue

            if probability <= 30:
                color = NHC_COLORS["low"]
            elif probability <= 60:
                color = NHC_COLORS["medium"]
            else:
                color = NHC_COLORS["high"]

            ax.add_geometries([row.geometry], crs=ccrs.PlateCarree(), 
                              facecolor="none", edgecolor=color, linewidth=2.5, zorder=7)

    if points is not None:
        for _, row in points.iterrows():
            try:
                probability = int(str(row["prob7day"]).replace("%", "").strip())
            except Exception:
                continue

            if probability <= 30:
                color = NHC_COLORS["low"]
            elif probability <= 60:
                color = NHC_COLORS["medium"]
            else:
                color = NHC_COLORS["high"]

            ax.scatter(row.geometry.x,row.geometry.y, s=85, marker="o", c=color, edgecolors="black", 
                       linewidths=1.0, transform=ccrs.PlateCarree(), zorder=9, clip_on=True)

            ax.text(row.geometry.x + 0.7, row.geometry.y + 0.4, f'{row["prob7day"]}', fontsize=10, 
                    fontweight="bold", color="black", transform=ccrs.PlateCarree(), zorder=15, clip_on=True)


    issue_time = outlook.get("issue_time") or now_utc
    valid_text = (f'Updated: {issue_time.strftime("%H%MZ %a %b %d, %Y").upper()}  •  7-Day Formation Probability + Active Tropical Cyclones')
    outlook_filename = build_filename(OUTPUT_DIR, "outlook", now_utc, variant="10-nhc-7d-tropical",)


    figure_builder(fig, ax,
        title="7-Day Tropical Weather Outlook",
        subtitle="NOAA National Hurricane Center • Atlantic and East Pacific Basins",
        valid=valid_text,
        category_colors=[
            NHC_COLORS["low"],
            NHC_COLORS["medium"],
            NHC_COLORS["high"],
            NHC_STORM_COLORS["D"],
            NHC_STORM_COLORS["S"],
            NHC_STORM_COLORS["H"],
            NHC_STORM_COLORS["M"]],
        category_labels=["LOW\n≤30%", "MED\n40-\n60%", "HIGH\n≥70%", "TD", "TS", "HUR", "MAJ\nHUR"],
        category_title=None,
        footer_left="NOAA National Hurricane Center • https://www.nhc.noaa.gov/gtwo.php",
        save_path=outlook_filename)

    plt.close(fig)
    print("    FINISHED NHC 7-DAY TROPICAL WEATHER OUTLOOK")
else:
    print("    SKIPPING NHC 7-DAY TROPICAL WEATHER OUTLOOK")


#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################











#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
#################################
# BUILD U.S. DROUGHT MONITOR
#################################
outlook = get_us_drought_monitor()

if outlook is not None:
    fig, ax = map_builder(extent=CONUS_EXTENT, terrain=True, counties=True, county_alpha=0.9)

    for drought_category in [0, 1, 2, 3, 4]:
        data = outlook[outlook["DM"] == drought_category]
        if data.empty:
            continue

        ax.add_geometries(data.geometry, crs=ccrs.PlateCarree(), facecolor=DROUGHT_COLORS[drought_category], edgecolor="black",
                           linewidth=0.45, alpha=0.55, zorder=5 + drought_category)

    issue_time = outlook.attrs.get("issue_time") or now_utc
    valid_start = outlook.attrs.get("valid_start")

    if valid_start is not None:
        valid_text = f'Valid: {valid_start.strftime("%a %b %d, %Y").upper()}'
    else:
        valid_text = f'Released: {issue_time.strftime("%a %b %d, %Y").upper()}'

    outlook_filename = build_filename(OUTPUT_DIR, "outlook", now_utc, variant="11-us-drought-monitor")

    figure_builder(fig, ax,
        title="U.S. Drought Monitor",
        subtitle="National Drought Mitigation Center • NOAA • USDA",
        valid=valid_text,
        category_colors=[DROUGHT_COLORS[i] for i in [0, 1, 2, 3, 4]],
        category_labels=DROUGHT_LABELS,
        category_title=None,
        footer_left="U.S. Drought Monitor • https://droughtmonitor.unl.edu/",
        save_path=outlook_filename)

    plt.close(fig)
    print("    FINISHED U.S. DROUGHT MONITOR")
else:
    print("    SKIPPING U.S. DROUGHT MONITOR")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################



#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
#################################
# BUILD WPC DAY 1 WINTER STORM OUTLOOK - SNOW
#################################
outlook = get_wpc_winter_storm_outlook(day=1, hazard="snow")

if outlook is not None:
    fig, ax = map_builder(extent=CONUS_EXTENT, terrain=True, counties=True, county_alpha=0.9)

    for threshold in ["10%", "30%", "50%", "80%"]:
        data = outlook[outlook["outlook"].astype(str).str.startswith(threshold)]
        if data.empty:
            continue

        ax.add_geometries(
            data.geometry,
            crs=ccrs.PlateCarree(),
            facecolor=WSO_COLORS[threshold],
            edgecolor="black",
            linewidth=0.8,
            alpha=0.76,
            zorder=6,
        )

    issue_time = outlook.attrs.get("issue_time") or now_utc
    valid_raw = str(outlook["valid_time"].iloc[0]) if "valid_time" in outlook.columns else ""
    valid_text = f'Issued: {issue_time.strftime("%H%MZ %a %b %d, %Y").upper()}'
    if valid_raw and valid_raw.lower() != "nan":
        valid_text += f'  •  Valid: {valid_raw}'

    outlook_filename = build_filename(OUTPUT_DIR, "outlook", now_utc, variant="12-wpc-wso-snow-d1")

    figure_builder(fig, ax,
        title="Day 1 Winter Storm Outlook • Snow",
        subtitle="NOAA Weather Prediction Center • Probability of Exceeding Warning Criteria",
        valid=valid_text,
        category_colors=[WSO_COLORS[item] for item in ["10%", "30%", "50%", "80%"]],
        category_labels=WSO_LABELS,
        category_title=None,
        footer_left="NOAA Weather Prediction Center Winter Storm Outlook • https://www.wpc.ncep.noaa.gov/wwd/wso/",
        save_path=outlook_filename)

    plt.close(fig)
    print("    FINISHED WPC DAY 1 WINTER STORM OUTLOOK - SNOW")
else:
    print("    SKIPPING WPC DAY 1 WINTER STORM OUTLOOK - SNOW")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################



#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
#################################
# BUILD WPC DAY 2 WINTER STORM OUTLOOK - SNOW
#################################
outlook = get_wpc_winter_storm_outlook(day=2, hazard="snow")

if outlook is not None:
    fig, ax = map_builder(extent=CONUS_EXTENT, terrain=True, counties=True, county_alpha=0.9)

    for threshold in ["10%", "30%", "50%", "80%"]:
        data = outlook[outlook["outlook"].astype(str).str.startswith(threshold)]
        if data.empty:
            continue

        ax.add_geometries(
            data.geometry,
            crs=ccrs.PlateCarree(),
            facecolor=WSO_COLORS[threshold],
            edgecolor="black",
            linewidth=0.8,
            alpha=0.76,
            zorder=6,
        )

    issue_time = outlook.attrs.get("issue_time") or now_utc
    valid_raw = str(outlook["valid_time"].iloc[0]) if "valid_time" in outlook.columns else ""
    valid_text = f'Issued: {issue_time.strftime("%H%MZ %a %b %d, %Y").upper()}'
    if valid_raw and valid_raw.lower() != "nan":
        valid_text += f'  •  Valid: {valid_raw}'

    outlook_filename = build_filename(OUTPUT_DIR, "outlook", now_utc, variant="13-wpc-wso-snow-d2")

    figure_builder(fig, ax,
        title="Day 2 Winter Storm Outlook • Snow",
        subtitle="NOAA Weather Prediction Center • Probability of Exceeding Warning Criteria",
        valid=valid_text,
        category_colors=[WSO_COLORS[item] for item in ["10%", "30%", "50%", "80%"]],
        category_labels=WSO_LABELS,
        category_title=None,
        footer_left="NOAA Weather Prediction Center Winter Storm Outlook • https://www.wpc.ncep.noaa.gov/wwd/wso/",
        save_path=outlook_filename)

    plt.close(fig)
    print("    FINISHED WPC DAY 2 WINTER STORM OUTLOOK - SNOW")
else:
    print("    SKIPPING WPC DAY 2 WINTER STORM OUTLOOK - SNOW")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################



#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
#################################
# BUILD WPC DAY 3 WINTER STORM OUTLOOK - SNOW
#################################
outlook = get_wpc_winter_storm_outlook(day=3, hazard="snow")

if outlook is not None:
    fig, ax = map_builder(extent=CONUS_EXTENT, terrain=True, counties=True, county_alpha=0.9)

    for threshold in ["10%", "30%", "50%", "80%"]:
        data = outlook[outlook["outlook"].astype(str).str.startswith(threshold)]
        if data.empty:
            continue

        ax.add_geometries(
            data.geometry,
            crs=ccrs.PlateCarree(),
            facecolor=WSO_COLORS[threshold],
            edgecolor="black",
            linewidth=0.8,
            alpha=0.76,
            zorder=6,
        )

    issue_time = outlook.attrs.get("issue_time") or now_utc
    valid_raw = str(outlook["valid_time"].iloc[0]) if "valid_time" in outlook.columns else ""
    valid_text = f'Issued: {issue_time.strftime("%H%MZ %a %b %d, %Y").upper()}'
    if valid_raw and valid_raw.lower() != "nan":
        valid_text += f'  •  Valid: {valid_raw}'

    outlook_filename = build_filename(OUTPUT_DIR, "outlook", now_utc, variant="14-wpc-wso-snow-d3")

    figure_builder(fig, ax,
        title="Day 3 Winter Storm Outlook • Snow",
        subtitle="NOAA Weather Prediction Center • Probability of Exceeding Warning Criteria",
        valid=valid_text,
        category_colors=[WSO_COLORS[item] for item in ["10%", "30%", "50%", "80%"]],
        category_labels=WSO_LABELS,
        category_title=None,
        footer_left="NOAA Weather Prediction Center Winter Storm Outlook • https://www.wpc.ncep.noaa.gov/wwd/wso/",
        save_path=outlook_filename)

    plt.close(fig)
    print("    FINISHED WPC DAY 3 WINTER STORM OUTLOOK - SNOW")
else:
    print("    SKIPPING WPC DAY 3 WINTER STORM OUTLOOK - SNOW")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################



#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
#################################
# BUILD WPC DAY 1 WINTER STORM OUTLOOK - FREEZING RAIN
#################################
outlook = get_wpc_winter_storm_outlook(day=1, hazard="freezing_rain")

if outlook is not None:
    fig, ax = map_builder(extent=CONUS_EXTENT, terrain=True, counties=True, county_alpha=0.9)

    for threshold in ["10%", "30%", "50%", "80%"]:
        data = outlook[outlook["outlook"].astype(str).str.startswith(threshold)]
        if data.empty:
            continue

        ax.add_geometries(
            data.geometry,
            crs=ccrs.PlateCarree(),
            facecolor=WSO_COLORS[threshold],
            edgecolor="black",
            linewidth=0.8,
            alpha=0.76,
            zorder=6,
        )

    issue_time = outlook.attrs.get("issue_time") or now_utc
    valid_raw = str(outlook["valid_time"].iloc[0]) if "valid_time" in outlook.columns else ""
    valid_text = f'Issued: {issue_time.strftime("%H%MZ %a %b %d, %Y").upper()}'
    if valid_raw and valid_raw.lower() != "nan":
        valid_text += f'  •  Valid: {valid_raw}'

    outlook_filename = build_filename(OUTPUT_DIR, "outlook", now_utc, variant="15-wpc-wso-ice-d1")

    figure_builder(fig, ax,
        title="Day 1 Winter Storm Outlook • Freezing Rain",
        subtitle="NOAA Weather Prediction Center • Probability of Exceeding Warning Criteria",
        valid=valid_text,
        category_colors=[WSO_COLORS[item] for item in ["10%", "30%", "50%", "80%"]],
        category_labels=WSO_LABELS,
        category_title=None,
        footer_left="NOAA Weather Prediction Center Winter Storm Outlook • https://www.wpc.ncep.noaa.gov/wwd/wso/",
        save_path=outlook_filename)

    plt.close(fig)
    print("    FINISHED WPC DAY 1 WINTER STORM OUTLOOK - FREEZING RAIN")
else:
    print("    SKIPPING WPC DAY 1 WINTER STORM OUTLOOK - FREEZING RAIN")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################



#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
#################################
# BUILD WPC DAY 2 WINTER STORM OUTLOOK - FREEZING RAIN
#################################
outlook = get_wpc_winter_storm_outlook(day=2, hazard="freezing_rain")

if outlook is not None:
    fig, ax = map_builder(extent=CONUS_EXTENT, terrain=True, counties=True, county_alpha=0.9)

    for threshold in ["10%", "30%", "50%", "80%"]:
        data = outlook[outlook["outlook"].astype(str).str.startswith(threshold)]
        if data.empty:
            continue

        ax.add_geometries(
            data.geometry,
            crs=ccrs.PlateCarree(),
            facecolor=WSO_COLORS[threshold],
            edgecolor="black",
            linewidth=0.8,
            alpha=0.76,
            zorder=6,
        )

    issue_time = outlook.attrs.get("issue_time") or now_utc
    valid_raw = str(outlook["valid_time"].iloc[0]) if "valid_time" in outlook.columns else ""
    valid_text = f'Issued: {issue_time.strftime("%H%MZ %a %b %d, %Y").upper()}'
    if valid_raw and valid_raw.lower() != "nan":
        valid_text += f'  •  Valid: {valid_raw}'

    outlook_filename = build_filename(OUTPUT_DIR, "outlook", now_utc, variant="16-wpc-wso-ice-d2")

    figure_builder(fig, ax,
        title="Day 2 Winter Storm Outlook • Freezing Rain",
        subtitle="NOAA Weather Prediction Center • Probability of Exceeding Warning Criteria",
        valid=valid_text,
        category_colors=[WSO_COLORS[item] for item in ["10%", "30%", "50%", "80%"]],
        category_labels=WSO_LABELS,
        category_title=None,
        footer_left="NOAA Weather Prediction Center Winter Storm Outlook • https://www.wpc.ncep.noaa.gov/wwd/wso/",
        save_path=outlook_filename)

    plt.close(fig)
    print("    FINISHED WPC DAY 2 WINTER STORM OUTLOOK - FREEZING RAIN")
else:
    print("    SKIPPING WPC DAY 2 WINTER STORM OUTLOOK - FREEZING RAIN")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################



#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################
#################################
# BUILD WPC DAY 3 WINTER STORM OUTLOOK - FREEZING RAIN
#################################
outlook = get_wpc_winter_storm_outlook(day=3, hazard="freezing_rain")

if outlook is not None:
    fig, ax = map_builder(extent=CONUS_EXTENT, terrain=True, counties=True, county_alpha=0.9)

    for threshold in ["10%", "30%", "50%", "80%"]:
        data = outlook[outlook["outlook"].astype(str).str.startswith(threshold)]
        if data.empty:
            continue

        ax.add_geometries(
            data.geometry,
            crs=ccrs.PlateCarree(),
            facecolor=WSO_COLORS[threshold],
            edgecolor="black",
            linewidth=0.8,
            alpha=0.76,
            zorder=6,
        )

    issue_time = outlook.attrs.get("issue_time") or now_utc
    valid_raw = str(outlook["valid_time"].iloc[0]) if "valid_time" in outlook.columns else ""
    valid_text = f'Issued: {issue_time.strftime("%H%MZ %a %b %d, %Y").upper()}'
    if valid_raw and valid_raw.lower() != "nan":
        valid_text += f'  •  Valid: {valid_raw}'

    outlook_filename = build_filename(OUTPUT_DIR, "outlook", now_utc, variant="17-wpc-wso-ice-d3")

    figure_builder(fig, ax,
        title="Day 3 Winter Storm Outlook • Freezing Rain",
        subtitle="NOAA Weather Prediction Center • Probability of Exceeding Warning Criteria",
        valid=valid_text,
        category_colors=[WSO_COLORS[item] for item in ["10%", "30%", "50%", "80%"]],
        category_labels=WSO_LABELS,
        category_title=None,
        footer_left="NOAA Weather Prediction Center Winter Storm Outlook • https://www.wpc.ncep.noaa.gov/wwd/wso/",
        save_path=outlook_filename)

    plt.close(fig)
    print("    FINISHED WPC DAY 3 WINTER STORM OUTLOOK - FREEZING RAIN")
else:
    print("    SKIPPING WPC DAY 3 WINTER STORM OUTLOOK - FREEZING RAIN")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################


elapsed_time = comp_time.time() - st
print(f"############\nSCRIPT FINISHED: time: {comp_time.strftime('%H:%M:%S', comp_time.gmtime(elapsed_time))}\n############")
