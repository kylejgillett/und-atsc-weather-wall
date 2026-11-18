##########################################################
#             RADAR MOSAIC RETREVIAL SCRIPT
#  (c) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2025
##########################################################


from siphon.catalog import TDSCatalog
import numpy as np
from netCDF4 import Dataset
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

def fetch_latest_mosaic(year, month, day):
    """
    Download the latest radar mosaic from NCSS as raw NetCDF bytes.

    Network only, so it is safe to run in a background thread. Pass the
    result to get_latest_mosaic(..., raw_data=...) to open it on the
    main thread (HDF5 is not thread-safe).
    """
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

    return ncss1.get_data_raw(query)


def get_latest_mosaic(year, month, day, raw_data=None):

    # download now unless the raw bytes were already fetched
    if raw_data is None:
        raw_data = fetch_latest_mosaic(year, month, day)

    radar_data = Dataset('radar_mosaic.nc', memory=raw_data)
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
            dBz = grb.values

            # MRMS is a regular lat/lon grid, so build 1D lat/lon from the
            # GRIB header instead of decoding full 2D lat/lon arrays
            # (~24M points), which is much slower
            lat_step = grb['jDirectionIncrementInDegrees'] * (1 if grb['jScansPositively'] else -1)
            lats_1d = grb['latitudeOfFirstGridPointInDegrees'] + np.arange(grb.Nj) * lat_step
            lons_1d = grb['longitudeOfFirstGridPointInDegrees'] + np.arange(grb.Ni) * grb['iDirectionIncrementInDegrees']
            lons_1d = (lons_1d + 180) % 360 - 180

            rows = np.where((lats_1d >= south) & (lats_1d <= north))[0]
            cols = np.where((lons_1d >= west) & (lons_1d <= east))[0]
            if rows.size == 0 or cols.size == 0:
                raise ValueError(
                    "MRMS geographic subset is empty.\n"
                    f"Requested: N={north}, S={south}, "
                    f"W={west}, E={east}\n"
                    f"Available latitude: "
                    f"{lats_1d.min():.2f} to {lats_1d.max():.2f}\n"
                    f"Available longitude: "
                    f"{lons_1d.min():.2f} to {lons_1d.max():.2f}")

            # Extract rectangular geographic subset
            r0, r1 = rows.min(), rows.max() + 1
            c0, c1 = cols.min(), cols.max() + 1

            dBz = dBz[r0:r1, c0:c1]
            radar_lon, radar_lat = np.meshgrid(lons_1d[c0:c1], lats_1d[r0:r1])

            # Extract radar valid time
            time = grb.validDate.replace(tzinfo=timezone.utc)

    # Mask missing data and weak echoes
    dBz = np.ma.masked_invalid(dBz)
    dBz = np.ma.masked_less(dBz, 10)

    # Diagnostics
    print(f"    MRMS COMPOSITE LOADED.....{time:%Y-%m-%d %H:%M}Z")
    print(f"    MRMS SUBSET: {dBz.shape} | LAT: {radar_lat.min():.2f} to {radar_lat.max():.2f} | LON: {radar_lon.min():.2f} to {radar_lon.max():.2f}")

    return dBz, radar_lat, radar_lon, time