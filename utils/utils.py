##########################################################
#            UND WX MAPWALL UTILITES MODULE
#  (c) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2025
##########################################################

# imports
import math



# convert u and v components to cardinal direction strings 
def wind_to_dir(u, v):

    wind_dir_deg = (math.degrees(math.atan2(u, v)) + 360) % 360

    dirs = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
            "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"]

    # Each sector covers 360/16 = 22.5 degrees
    idx = int((wind_dir_deg + 11.25) // 22.5) % 16

    return dirs[idx]




# convert datetime obj day to abbreviated day name str
def day_to_abbrev(date_obj):
    abbr = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]
    return abbr[date_obj.weekday()]