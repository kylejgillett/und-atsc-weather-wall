# ================================================================
# UND WEATHER WALL STYLE COLORS
# KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2026
# ================================================================

UND_GREEN = "#009A44"

FIGURE_BG = "#2A2B28"

TEXT_PRIMARY = "#F5F4EF"
TEXT_SECONDARY = "#D8D5CB"
TEXT_VALID = "#E2E5DE"
TEXT_MUTED = "#979B91"

LAND_COLOR = "#242521"
WATER_COLOR = "#3D3E37"

STATE_COLOR = "#D2D0C7"
BORDER_COLOR = "#ECE9DF"
COUNTY_COLOR = "#9A9C90"
ROAD_COLOR = "#A7A79D"

DIVIDER_COLOR = "#3A3B36"
COLORBAR_EDGE = "#1A1B18"
LEGEND_EDGE = "#1A1B18"
LEGEND_BG = "#2C2C2C"

PANEL_BG = "#343530"
GRID_COLOR = "#666860"
SPINE_COLOR = "#7D8077"

TEMP_COLOR = "#F26B4A"
FEEL_COLOR = "#E5B35A"
DEWP_COLOR = "#67C587"

VIS_COLOR = "#E2DFD5"
RH_COLOR = "#7FBE72"

WIND_COLOR = "#D7D4C9"
WIND_FILL_COLOR = "#B98AF0"
WIND_FILL_DARK = "#4A3D5E"
WIND_FILL_LIGHT = "#F2EAFF"
GUST_COLOR = "#E99645"
MSLP_COLOR = "#C2A6D8"

FREEZE_COLOR = "#AFC7D5"

ENS_COLOR = "#A9B8C8"
ENS_MEDIAN_COLOR = "#D5E2EE"
QPF_COLOR = "#5CC98A"

HIGH_COLOR = "#FF6A5C"
LOW_COLOR = "#4FA3FF"

RAIN_COLOR = "#22B14C"
SNOW_COLOR = "#5CB8FF"
MIX_COLOR = "#FF8C1A"
FZRA_COLOR = "#FF6FA8"
TSTM_COLOR = "#F2C14E"
POP_COLOR = "#6F726B"
WX_COLOR = "#7A7466"
SKY_COLOR = "#C9CCC4"

# sounding background: between SounderPy's dark default (#161616) and FIGURE_BG, for more contrast
SOUNDING_BG = "#1C1D1B"

# sounding temperature: a true red (TEMP_COLOR is coral) that stays bright on SOUNDING_BG
SOUNDING_TEMP_COLOR = "#F0453A"

# SounderPy theme for the sounding figures, matching the Weather Wall maps
# (passed as build_sounding/build_composite(theme=...), also used by the winter sounding)
SOUNDING_THEME = {"text": TEXT_PRIMARY, "background": SOUNDING_BG, "border": UND_GREEN,
                  "marker": TEXT_PRIMARY, "barb": WIND_COLOR,
                  "adiabat": UND_GREEN, "skew_line": TEXT_MUTED, "isotherm": FREEZE_COLOR,
                  "shade": TEXT_MUTED, "pressure_line": GRID_COLOR,
                  "temperature": SOUNDING_TEMP_COLOR, "virtual": "#A8352D", "warm": SOUNDING_TEMP_COLOR,
                  "parcel": {"sb": HIGH_COLOR, "ml": GUST_COLOR, "mu": WIND_FILL_COLOR},
                  "special_parcel": TEXT_SECONDARY, "secondary_parcel": TEXT_MUTED,
                  "height_marker": TEXT_PRIMARY, "height_marker_edge": SOUNDING_BG,
                  "moisture": DEWP_COLOR, "moisture_strong": "#8FE388",
                  "composite_index": HIGH_COLOR,
                  "bar_edge": COLORBAR_EDGE, "watermark": TEXT_SECONDARY}

# SounderPy watermark text for the sounding figures (logo path is set in each builder)
SOUNDING_WATERMARK = {"title": "UND Atmospheric Science • Weather Map Wall",
                      "subtitle": "in the Leon F. Osborne Weather Center • Built by SounderPy (Gillett 2024)",
                      "alpha": 0.9}
