##########################################################
#            UND WX MAPWALL UTILITES MODULE
#  (c) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2025
##########################################################

# imports
import math
import numpy as np
from datetime import timezone, datetime


# stanardized filename handler
def build_filename(storage, figtype, validtime, variant=None, order=None):

    time_str = validtime.astimezone(timezone.utc).strftime("%Y%m%d_%H%M%S")

    if variant is not None:
        filename = f"{storage}{figtype}_{time_str}_{variant}.png"
    else: 
        filename = f"{storage}{figtype}_{time_str}.png"

    return filename




# convert u and v components to cardinal direction strings 
def wind_to_dir(u, v):

    wind_dir_deg = (math.degrees(math.atan2(u, v)) + 360) % 360

    dirs = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
            "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"]

    # Each sector covers 360/16 = 22.5 degrees
    wind_dir_deg = (np.degrees(np.arctan2(-u, -v)) + 360) % 360
    idx = int((wind_dir_deg + 11.25) // 22.5) % 16

    return dirs[idx]




# convert datetime obj day to abbreviated day name str
def day_to_abbrev(date_obj):
    abbr = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]
    return abbr[date_obj.weekday()]