import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import ListedColormap
import numpy as np
import os





########### HIGH CONTRAST GREYS CMAP FOR SAT IR ###########
ir_colors = [(1, 1, 1), (0.6, 0.6, 0.6), (0.3, 0.3, 0.3), (0, 0, 0)]
ir_greys = mcolors.LinearSegmentedColormap.from_list("ir_greys", ir_colors, N=256)


########### DEWPOINT COLORMAP ###########
dwpt_colors = [
    (0.00, '#634235'), # -40: Dark Brown
    (0.20, '#96816D'), # -14: Muted Brown
    (0.30, '#EBE7D3'), # -1: Light Tan
    (0.40, '#D8E5C1'), # 12: Pale Green starts
    (0.55, '#6BB35F'), # 31: Vibrant Green
    (0.75, '#1E6E32'), # 57: Deep Green persists longer
    (0.82, '#3E5D61'), # 66: Teal (Transition begins later)
    (0.88, '#2E3A4A'), # 74: Deep Navy/Charcoal
    (0.94, '#5C4A8D'), # 82: Deep Purple
    (1.00, '#B28EA1')  # 90: Light Pink/Mauve
    ]

dwpt_cmap = LinearSegmentedColormap.from_list("dwpt_cmap", dwpt_colors, N=256)



########### PIVITOL WX STYLE TEMPERATURE CMAP ###########
temp_colors = [(0.000, '#B3ECE0'), (0.125, '#A08AC6'),(0.250, '#8C28AC'),(0.375, '#D4E2E8'),
               (0.4999,'#1450B4'), (0.5001,'#0F505F'), (0.625, '#F3F2A7'), (0.750, '#A56847'),
              (0.875, '#690B10'), (1.000, '#E8DFD6')]
temp_cmap = LinearSegmentedColormap.from_list("temp_cmap", temp_colors, N=256)



########### PIVITOL WX STYLE WIND SPEED CMAP ###########
wdsp_colors = [
    (0.000, "#bfbfbf05"), 
    (0.004, "#87CEFA1A"), 
    (0.166, "#695ACD3F"), 
    (0.250, '#E696DCFF'), 
    (0.333, '#C85ABEFF'),
    (0.416, '#A01496FF'), 
    (0.500, '#C80028FF'), 
    (0.583, '#DC283CFF'), 
    (0.666, '#F05050FF'),
    (0.750, '#FAF064FF'),
    (0.833, '#DCBE46FF'), 
    (0.916, '#BE8C28FF'),
    (1.000, '#A05A0AFF')
]
wdsp_cmap = LinearSegmentedColormap.from_list("wdsp_cmap", wdsp_colors, N=256)



########### PIVITOL WX STYLE PV CMAP ###########
pv_clevs  = np.append(np.arange(-1.4, 2, 0.2),np.arange(2,10.2,0.4))
pv_colors = ["blue","lightblue","lightblue","yellow","orange","red","darkred"]
pv_cmap   = LinearSegmentedColormap.from_list("pv_cmap",pv_colors)


########### RADARSCOPE STYLE REFLECTIVITY CMAP ###########
rs_colors = [(29/255, 37/255, 60/255),(89/255, 155/255, 171/255),(33/255, 186/255, 72/255),(5/255, 101/255, 1/255),
    [(251/255, 252/255, 0), (199/255, 176/255, 0)],[(253/255, 149/255, 2/255), (172/255, 92/255, 2/255)],
    [(253/255, 38/255, 0), (135/255, 43/255, 22/255)],[(193/255, 148/255, 179/255), (200/255, 23/255, 119/255)],
    [(165/255, 2/255, 215/255), (64/255, 0, 146/255)],[(135/255, 255/255, 253/255), (54/255, 120/255, 142/255)],
    (173/255, 99/255, 64/255),(105/255, 0, 4/255),(0, 0, 0)]
final_rs_colors = []
for c in rs_colors:
    if isinstance(c, list):
        final_rs_colors.extend(c)
    else:
        final_rs_colors.append(c)
rs_cmap = LinearSegmentedColormap.from_list("rs_reflect", final_rs_colors, N=256)


########### RADARSCOPE STYLE EXPERT MODE REFLECTIVITY CMAP ###########
script_dir = os.path.dirname(os.path.abspath(__file__))
txt_file_path = os.path.join(script_dir, "rs_basereflect_colormap.txt")
rs_data = np.loadtxt(txt_file_path, skiprows=3, usecols=(1, 2, 3, 4))
rgb = rs_data[:, 1:] / 255.0
levels = rs_data[:, 0]
rs_expertreflect_cmap = ListedColormap(rgb)

script_dir = os.path.dirname(os.path.abspath(__file__))
txt_file_path = os.path.join(script_dir,"rs_basevelocity_colormap.txt")
rs_data = np.loadtxt(txt_file_path, skiprows=3, usecols=(1, 2, 3, 4))
rgb = rs_data[:, 1:] / 255.0
levels = rs_data[:, 0]
rs_basevelocity_cmap = ListedColormap(rgb)



########### NWS STYLE REFLECTIVITY CMAP ###########
nws_cmap  = LinearSegmentedColormap.from_list('nws_reflect',['lightsteelblue','steelblue','lightgreen','forestgreen',
                                                               [255/255,255/255,77/255],[230/255,230/255,0],[255/255,195/255,77/255],
                                                               [230/255,153/255,0],[255/255,77/255,77/255],[230/255,0,0],[255/255,204/255,238/255],
                                                               [255/255,25/255,140/255],[212/255,0,255/255],[85/255,0,128/255]],N=256)



cape_colors = [
    # CAPE colormap (based on provided guide) with increasing opacity
    ("#6b8d97f1"),  # very light gray, mostly transparent
    ("#4c659bf9"),  # light blue
    ("#397efffc"),  # cyan
    ("#40ffc9ff"),  # green
    ("#40ff6dff"),  # yellow
    ("#f2ff00ff"),  # orange
    ("#ffac40ff"),  # red
    ("#c04000ff"),  # magenta
    ("#800037ff"),  # purple, fully opaque
]

cape_cmap = LinearSegmentedColormap.from_list("cape_cmap", cape_colors, N=256)
