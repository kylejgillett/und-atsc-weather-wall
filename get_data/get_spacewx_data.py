from urllib.request import urlopen
import urllib.request
import json
from datetime import datetime, timedelta
from numpy import loadtxt
import numpy as np



### GET ACE DATA ###
#################################################################################################
def get_ace_data(): 
    
    spd_data  = []
    dens_data = []
    temp_data = []
    bx_data   = []
    by_data   = []
    bz_data   = []
    bt_data   = []
    mag_dates = []
    pam_dates = []
    
    # MAG DATA ---------------------------------------------
    url ='https://services.swpc.noaa.gov/json/rtsw/rtsw_mag_1m.json'
    response = urllib.request.urlopen(url)
    mag_data = json.loads(response.read())
    mag_data = [item for item in mag_data if item.get('source') == 'ACE']

    for item in mag_data:
        mag_dates.append(datetime.strptime(str(item['time_tag']), '%Y-%m-%dT%H:%M:%S'))
        bz_data.append(item['bz_gsm'])
        bt_data.append(item['bt'])
        bx_data.append(item['bx_gsm'])
        by_data.append(item['by_gsm'])



    # SWEPAM DATA ---------------------------------------------
    url ='https://services.swpc.noaa.gov/json/rtsw/rtsw_wind_1m.json'
    response = urllib.request.urlopen(url)
    pam_data = json.loads(response.read())
    pam_data = [item for item in pam_data if item.get('source') == 'ACE']

    # parse out swepam data
    for item in pam_data:
        spd_data.append(item['proton_speed'])
        dens_data.append(item['proton_density']) 
        temp_data.append(item['proton_temperature'])
        pam_dates.append(datetime.strptime(str(item['time_tag']), '%Y-%m-%dT%H:%M:%S'))
        
        
    # FULL DATA DICT ----------------------------------------
    # add DSCOVR data to dictionary
    ace_data = {
    'pam_dates': np.array(pam_dates)[::-1],
    'mag_dates': np.array(mag_dates)[::-1],
    'density': np.array(dens_data, dtype=float)[::-1],
    'temperature': np.array(temp_data, dtype=float)[::-1],
    'speed': np.array(spd_data, dtype=float)[::-1],
    'bt': np.array(bt_data, dtype=float)[::-1],
    'bz': np.array(bz_data, dtype=float)[::-1],
    'bx': np.array(bx_data, dtype=float)[::-1],
    'by': np.array(by_data, dtype=float)[::-1],
    'meta': {'sat-name': 'ACE',
             'long-sat-name': 'NASA Advanced Composition Explorer',
             'data-name': 'MAG + SWEPAM',
             'long-data-name': '1-minute averaged Real-time Interplanetary Magnetic Field Values & Bulk Parameters of the Solar Wind Plasma',
             'missing': 'MAG: -999.9, DEN/SPD: -9999.9, TEMP: -1.00e+05',
             'source': 'Prepared by the U.S. Dept. of Commerce, NOAA, Space Weather Prediction Center',
        }
    }
    
    return ace_data



### GET KP-INDEX DATA ###
#################################################################################################
def get_kp_data():

    # get HTTP URL response for SWPC KP data 
    response = urllib.request.urlopen('https://services.swpc.noaa.gov/products/noaa-planetary-k-index.json')

    # load KP data
    data = json.loads(response.read())

    # perpare lists 
    kp_dates  = []
    kp_times = []
    kp_index = []

    # parse data lists 
    for kp_list in data[1:-1]:      
        kp_dates.append(datetime.strptime(f"{kp_list[0]}", '%Y-%m-%d %H:%M:%S.%f'))
        kp_index.append(kp_list[1])
        
    # build data dict
    kp_data = {
        'kp_dates': np.array(kp_dates),
        'kp_times': np.array([date.strftime('%H:%M') for date in np.array(kp_dates)]),
        'kp_index': np.array(kp_index),
        'current_kp': np.array(kp_index)[-1]
    }
    
    return kp_data
