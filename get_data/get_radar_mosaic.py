##########################################################
#             RADAR MOSAIC RETREVIAL SCRIPT
#  (c) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2025
##########################################################


from siphon.catalog import TDSCatalog
import numpy as np
from xarray.backends import NetCDF4DataStore
from xarray import open_dataset
import gzip
import os
import tempfile
from datetime import timezone
import pygrib
import requests


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




def get_latest_mrms(north=60, south=10, west=-125, east=-50):

    # MRMS QC Composite Reflectivity
    url = ("https://mrms.ncep.noaa.gov/2D/MergedReflectivityQCComposite/MRMS_MergedReflectivityQCComposite.latest.grib2.gz")

    # Download compressed GRIB2
    response = requests.get(url, timeout=90)
    response.raise_for_status()

    with tempfile.TemporaryDirectory() as tmpdir:

        grib_path = os.path.join(tmpdir, "mrms.grib2")

        # Decompress GRIB2
        with gzip.open if False else open(grib_path, "wb") as f:
            f.write(gzip.decompress(response.content))

        # Read MRMS field
        with pygrib.open(grib_path) as grbs:
            grb = grbs.message(1)
            dBz, radar_lat, radar_lon = grb.data()
            radar_lon = (radar_lon + 180) % 360 - 180

            mask = ((radar_lat >= south) & (radar_lat <= north) & (radar_lon >= west) & (radar_lon <= east))
            rows, cols = np.where(mask)
            if rows.size == 0:
                raise ValueError(
                    "MRMS geographic subset is empty.\n"
                    f"Requested: N={north}, S={south}, "
                    f"W={west}, E={east}\n"
                    f"Available latitude: "
                    f"{radar_lat.min():.2f} to {radar_lat.max():.2f}\n"
                    f"Available longitude: "
                    f"{radar_lon.min():.2f} to {radar_lon.max():.2f}")

            # Extract rectangular geographic subset
            r0, r1 = rows.min(), rows.max() + 1
            c0, c1 = cols.min(), cols.max() + 1

            dBz = dBz[r0:r1, c0:c1]
            radar_lat = radar_lat[r0:r1, c0:c1]
            radar_lon = radar_lon[r0:r1, c0:c1]

            # Extract radar valid time
            time = grb.validDate.replace(tzinfo=timezone.utc)

    # Mask missing data and weak echoes
    dBz = np.ma.masked_invalid(dBz)
    dBz = np.ma.masked_less(dBz, 10)

    # Diagnostics
    print(f"    MRMS COMPOSITE LOADED.....{time:%Y-%m-%d %H:%M}Z")
    print(f"    MRMS SUBSET: {dBz.shape} | LAT: {radar_lat.min():.2f} to {radar_lat.max():.2f} | LON: {radar_lon.min():.2f} to {radar_lon.max():.2f}")

    return dBz, radar_lat, radar_lon, time