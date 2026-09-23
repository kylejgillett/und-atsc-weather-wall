##########################################################
#           REGIONAL OBS SOUNDING PLOTS SCRIPT
#  (c) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2025
##########################################################

print("############\nSCRIPT RUNNING\n############")
import time as comp_time
st = comp_time.time()

import warnings
warnings.filterwarnings("ignore")


import sounderpy as spy
from datetime import datetime, timedelta, timezone
import os
import sys

# get script dir
script_dir = os.path.dirname(os.path.abspath(__file__))

# get the parent dir 
project_root = os.path.abspath(os.path.join(script_dir, ".."))
if project_root not in sys.path:
    sys.path.append(project_root)
from utils.utils import *

    

utc_now = datetime.now(timezone.utc)
current_date = utc_now.date()

if utc_now.hour >= 18:
    set_date = utc_now.replace(
        hour=18, minute=0, second=0, microsecond=0)

elif utc_now.hour >= 12:
    set_date = utc_now.replace(
        hour=12, minute=0, second=0, microsecond=0)

elif utc_now.hour >= 0:
    set_date = utc_now.replace(
        hour=0, minute=0, second=0, microsecond=0)


year_str  = set_date.strftime("%Y")
month_str = set_date.strftime("%m")
day_str   = set_date.strftime("%d")
hour_str  = set_date.strftime("%H")


#################################
# OBS SOUNDINGS
#################################
ids = ["kbis", "kabr", "kggw", "kunr", "kinl", "kmpx"]

for id, i in zip(ids, range(3,len(ids)+3)):
    try:
        data = spy.get_obs_data(id, year_str, month_str, day_str, hour_str)

        snd_date = data['site_info']['valid-time']
        #snd_dt = datetime(int(snd_date[0]), int(snd_date[1]), int(snd_date[2]), int(snd_date[3]), tzinfo=timezone.utc)
        obs_filename = build_filename("staged_figures/soundings/", "sounding", set_date, variant=f"{i:02d}-obs-{id}")

        spy.build_sounding(data, special_parcels='simple', dark_mode=True, map_zoom=1, color_blind=True, radar='mosaic',
                    save=True, filename=obs_filename)
    except:
        print(f"    !!! NO DATA FOUND FOR {id.upper()}")


print(f"    FINISHED {id.upper()} OBS SOUNDING")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################


elapsed_time = comp_time.time() - st
print(f"############\nSCRIPT FINISHED: time: {comp_time.strftime('%H:%M:%S', comp_time.gmtime(elapsed_time))}\n############")
