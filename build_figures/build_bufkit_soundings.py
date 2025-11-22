##########################################################
#           BUFKIT F00 SOUNDING PLOTS SCRIPT
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



#################################
# KGFK BUFKIT SOUNDING
#################################
gfk_data = spy.get_bufkit_data('hrrr', "kgfk", 0)

spy.build_sounding(gfk_data, special_parcels='simple', map_zoom=1, color_blind=True, radar='mosaic',
                   save=True, filename="staged_figures/soundings/gfk_rap_sounding.png")
print("    FINISHED GFK BUFKIT SOUNDING")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################



#################################
# KGFK BUFKIT SOUNDING COMPOSITE
#################################

utc_now = datetime.now(timezone.utc)
offsets = [24, 18, 12, 6]

datas = []
for hrs in offsets:

    # format date/time
    t = utc_now - timedelta(hours=hrs)
    yy = t.strftime("%Y")
    mm = t.strftime("%m")
    dd = t.strftime("%d")
    hh = t.strftime("%H")

    # get data 
    datas.append(spy.get_bufkit_data('hrrr', "kgfk", 0, yy, mm, dd, hh))

datas.append(gfk_data)

spy.build_composite(datas, cmap='copper_r', lw_to_use=[4 for data in datas], alphas_to_use=[1, 0.9, 0.8, 0.7, 0.6][::-1],
                   save=True, filename="staged_figures/soundings/gfk_rap_composite_sounding.png")


print("    FINISHED GFK EVOLUTION COMPOSITE SOUNDING")
#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################


elapsed_time = comp_time.time() - st
print(f"############\nSCRIPT FINISHED: time: {comp_time.strftime('%H:%M:%S', comp_time.gmtime(elapsed_time))}\n############")
