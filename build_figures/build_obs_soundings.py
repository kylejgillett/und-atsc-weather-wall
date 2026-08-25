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

    

utc_now = datetime.now(timezone.utc)
current_date = utc_now.date()

# decide search time (00z or 12z)
if utc_now.hour >= 12:
    # after 12z, use 12z of the current day
    search_date = datetime(current_date.year, current_date.month, current_date.day, 12)
else:
    # between 00z and 12z, use 00z of the current day
    search_date = datetime(current_date.year, current_date.month, current_date.day, 0)


year_str   = search_date.strftime("%Y")
month_str  = search_date.strftime("%m")
day_str    = search_date.strftime("%d")
hour_str   = search_date.strftime("%H")


#################################
# OBS SOUNDINGS
#################################
ids = ["kggw", "kbis", "kabr", "kinl", "kmpx", "kunr"]

for id in ids:
    try:
        data = spy.get_obs_data(id, year_str, month_str, day_str, hour_str)

        spy.build_sounding(data, special_parcels='simple', map_zoom=1, color_blind=True, radar='mosaic',
                    save=True, filename=f"staged_figures/soundings/{id}_obs_sounding.png")
    except:
        print(f"    !!! NO DATA FOUND FOR {id.upper()}")


print(f"    FINISHED {id.upper()} OBS SOUNDING")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################


elapsed_time = comp_time.time() - st
print(f"############\nSCRIPT FINISHED: time: {comp_time.strftime('%H:%M:%S', comp_time.gmtime(elapsed_time))}\n############")
