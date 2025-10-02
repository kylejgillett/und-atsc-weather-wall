##########################################################
#             RADAR MOSAIC RETREVIAL SCRIPT
#  (c) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2025
##########################################################


from siphon.catalog import TDSCatalog
import numpy as np
from xarray.backends import NetCDF4DataStore
from xarray import open_dataset

# Download radar reflectivity data
north = 60
south = 10
west = -125
east = -50

def get_latest_mosaic(year, month, day):
    datestr = f'{year}{month}{day}'

    composite_url = 'https://thredds.ucar.edu/thredds/catalog/nexrad/composite/gini/dhr/1km/'+datestr+'/catalog.xml'
    best_radar = TDSCatalog(composite_url)
    radar_ds   = best_radar.datasets
    ncss1      = radar_ds[0].subset()
    query      = ncss1.query()

    query.lonlat_box(north=north+1,
                     south=south-1,
                     east=east+1,
                     west=west-1)

    query.add_lonlat(value=True)
    query.accept('netcdf4')
    query.variables('Reflectivity')
    radar_data = ncss1.get_data(query)
    radar_data = open_dataset(NetCDF4DataStore(radar_data))
    time = radar_data.time.values[0]
    radar_lat  = np.array(radar_data['lat'])
    radar_lon  = np.array(radar_data['lon'])
    radar_lon  = np.array(radar_data['lon'])

    dBz = np.array(radar_data['Reflectivity'])[0,:,:]
    dBz = np.ma.masked_array(dBz,dBz<10)

    print(f"    RADAR MOSAIC LOADED.....{time}z")

    return dBz, radar_lat, radar_lon, time

