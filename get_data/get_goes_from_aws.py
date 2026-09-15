##########################################################
#              GOES-19 ABI RETRIEVAL SCRIPT
#  (c) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2026
##########################################################

from datetime import datetime, timedelta
from pathlib import Path
import cartopy.crs as ccrs
import s3fs
import numpy as np

# VISIBLE 
def download_goes19_visible(year, day_of_year, hour, output_dir="../temp_files", lookback_hours=2):

    fs = s3fs.S3FileSystem(anon=True)
    request_time = datetime.strptime(f"{int(year):04d}{int(day_of_year):03d}{int(hour):02d}", "%Y%j%H")

    latest_file = None

    for offset in range(lookback_hours + 1):
        search_time = request_time - timedelta(hours=offset)
        prefix = (f"noaa-goes19/ABI-L2-CMIPC/"
                f"{search_time:%Y}/"
                f"{search_time:%j}/"
                f"{search_time:%H}/")

        files = sorted(fs.glob(prefix + "*C02_G19*.nc"))

        if files:
            latest_file = files[-1]
            break

    if latest_file is None:
        print("    NO GOES-19 VISIBLE DATA FOUND")
        return None

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    local_file = output_dir / Path(latest_file).name

    if not local_file.exists():
        fs.get(latest_file, str(local_file))

    print(f"    GOES-19 VISIBLE DOWNLOADED.....{local_file.name}")
    return str(local_file)






# IR SAMMY FOR NIGHT TIME
def _get_scan_id(filepath):
    return Path(filepath).name.split("_s")[1].split("_e")[0]


def download_goes19_sandwich(year, day_of_year, hour, output_dir="../temp_files", lookback_hours=2):

    fs = s3fs.S3FileSystem(anon=True)
    request_time = datetime.strptime(f"{int(year):04d}{int(day_of_year):03d}{int(hour):02d}", "%Y%j%H")

    matched_files = None

    for offset in range(lookback_hours + 1):
        search_time = request_time - timedelta(hours=offset)
        prefix = (f"noaa-goes19/ABI-L2-CMIPC/"
                f"{search_time:%Y}/"
                f"{search_time:%j}/"
                f"{search_time:%H}/")

        files_c07 = sorted(fs.glob(prefix + "*C07_G19*.nc"))
        files_c13 = sorted(fs.glob(prefix + "*C13_G19*.nc"))

        if not files_c07 or not files_c13:
            continue

        # Map scan start time -> file
        scans_c07 = {_get_scan_id(f): f for f in files_c07}
        scans_c13 = {_get_scan_id(f): f for f in files_c13}

        # Find scans available for BOTH channels
        common_scans = sorted(set(scans_c07) & set(scans_c13))

        if not common_scans:
            continue

        # Latest matching scan
        latest_scan = common_scans[-1]
        matched_files = {"C07": scans_c07[latest_scan], "C13": scans_c13[latest_scan]}

        break

    if matched_files is None:
        print("    NO GOES-19 SANDWICH DATA FOUND")
        return None

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    local_files = {}

    for channel, remote_file in matched_files.items():
        local_file = output_dir / Path(remote_file).name
        if not local_file.exists():
            fs.get(remote_file, str(local_file))
        local_files[channel] = str(local_file)

    print(f"    GOES-19 SANDWICH DOWNLOADED.....C07 + C13")

    return local_files



# subset goes file to speed up plotting
def subset_goes_to_map(ds, west, east, south, north, pad_km=300):

    sat_crs = ds.FOV.crs
    map_crs = ccrs.PlateCarree()
    n = 50

    lons = np.concatenate([np.linspace(west, east, n),
                np.full(n, east),
                np.linspace(east, west, n),
                np.full(n, west)])

    lats = np.concatenate([np.full(n, south),
                np.linspace(south, north, n),
                np.full(n, north),
                np.linspace(north, south, n)])

    points = sat_crs.transform_points(map_crs, lons, lats)
    x_target = points[:, 0]
    y_target = points[:, 1]
    pad = pad_km * 1000.0
    xmin = np.nanmin(x_target) - pad
    xmax = np.nanmax(x_target) + pad
    ymin = np.nanmin(y_target) - pad
    ymax = np.nanmax(y_target) + pad
    x = np.asarray(ds.FOV.x.values)
    y = np.asarray(ds.FOV.y.values)

    x_idx = np.where((x >= xmin) &(x <= xmax))[0]

    y_idx = np.where((y >= ymin) &(y <= ymax))[0]

    x0 = x_idx.min()
    x1 = x_idx.max() + 1
    y0 = y_idx.min()
    y1 = y_idx.max() + 1

    visible = ds["CMI"].isel( x=slice(x0, x1), y=slice(y0, y1))

    x_sub = x[x0:x1]
    y_sub = y[y0:y1]

    dx = np.abs(np.nanmedian(np.diff(x_sub)))
    dy = np.abs(np.nanmedian(np.diff(y_sub)))

    extent = (np.nanmin(x_sub) - dx / 2,
              np.nanmax(x_sub) + dx / 2,
              np.nanmin(y_sub) - dy / 2,
              np.nanmax(y_sub) + dy / 2)

    print(f"    GOES SUBSET....."
        f"{visible.shape[1]}x{visible.shape[0]} "
        f"from {ds['CMI'].shape[1]}x{ds['CMI'].shape[0]}")

    return visible, sat_crs, extent










# ##########################################################
# #              GOES 19 ABI RETREVIAL SCRIPT
# #  (c) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2025
# ##########################################################
# import numpy as np
# import s3fs


# import warnings
# warnings.filterwarnings("ignore")



# def download_goes_file(year, day_of_year, hour, output_filename=None):
#     fs = s3fs.S3FileSystem(anon=True)

#     try:
#         files = np.array(fs.ls(f'noaa-goes19/ABI-L1b-RadF/{year}/{day_of_year}/{hour}/'))
#     except:
#         files = np.array(fs.ls(f'noaa-goes19/ABI-L1b-RadF/{year}/{day_of_year}/{str(int(hour)-1)}/'))
#         pass
#     if len(files) == 0:
#         print("     No files found for the specified date and time.")
#         return

#     latest_file = files[-1]

#     prefix = '../temp_files/'
#     if output_filename is None:
#         output_filename = f"{prefix}{latest_file.split('/')[-1]}"

#     fs.download(latest_file, output_filename)
#     print(f"    SATELLITE DATA LOADED.....{output_filename}")

#     return output_filename








# # def download_goes_file(year, day_of_year, hour):
# #     # Use the anonymous credentials to access public data
# #     fs = s3fs.S3FileSystem(anon=True)
# #
# #     # List contents of GOES-16 bucket.
# #     fs.ls('s3://noaa-goes16/')
# #
# #
# #     # List specific files of GOES-16/17 Full Disk/CONUS/Mesoscale sector data on a certain hour
# #     # Note: the `s3://` is not required
# #     # data structure is as such for goes 16
# #     # noaa-goes16/<Product>/<Year>/<Day of Year>/<Hour>
# #     # info at https://docs.opendata.aws/noaa-goes16/cics-readme.html
# #
# #     files = np.array(fs.ls(f'noaa-goes16/ABI-L1b-RadC/{year}/{day_of_year}/{hour}/'))
# #     print(files)
# #
# #     fs.download(files[-1], files[-1].split('/')[-1])
# #
# #     #fs.get(files[0:96], files[0:96].split('/')[-1])
# #
# #     # # Download the first file, and rename it the same name (without the directory structure). Downloads files into directory this script is run in.
# #     # for i in range(0, len(files[0:96])):
# #     #     fs.get(files[i], files[i].split('/')[-1])