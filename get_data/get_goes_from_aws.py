##########################################################
#              GOES 19 ABI RETREVIAL SCRIPT
#  (c) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2025
##########################################################
import numpy as np
import s3fs


def download_goes_file(year, day_of_year, hour, output_filename=None):
    fs = s3fs.S3FileSystem(anon=True)

    try:
        files = np.array(fs.ls(f'noaa-goes19/ABI-L1b-RadF/{year}/{day_of_year}/{hour}/'))
    except:
        files = np.array(fs.ls(f'noaa-goes19/ABI-L1b-RadF/{year}/{day_of_year}/{str(int(hour)-1)}/'))
        pass
    if len(files) == 0:
        print("     No files found for the specified date and time.")
        return

    latest_file = files[-1]

    prefix = '../temp_files/'
    if output_filename is None:
        output_filename = f"{prefix}{latest_file.split('/')[-1]}"

    fs.download(latest_file, output_filename)
    print(f"    SATELLITE DATA LOADED.....{output_filename}")

    return output_filename








# def download_goes_file(year, day_of_year, hour):
#     # Use the anonymous credentials to access public data
#     fs = s3fs.S3FileSystem(anon=True)
#
#     # List contents of GOES-16 bucket.
#     fs.ls('s3://noaa-goes16/')
#
#
#     # List specific files of GOES-16/17 Full Disk/CONUS/Mesoscale sector data on a certain hour
#     # Note: the `s3://` is not required
#     # data structure is as such for goes 16
#     # noaa-goes16/<Product>/<Year>/<Day of Year>/<Hour>
#     # info at https://docs.opendata.aws/noaa-goes16/cics-readme.html
#
#     files = np.array(fs.ls(f'noaa-goes16/ABI-L1b-RadC/{year}/{day_of_year}/{hour}/'))
#     print(files)
#
#     fs.download(files[-1], files[-1].split('/')[-1])
#
#     #fs.get(files[0:96], files[0:96].split('/')[-1])
#
#     # # Download the first file, and rename it the same name (without the directory structure). Downloads files into directory this script is run in.
#     # for i in range(0, len(files[0:96])):
#     #     fs.get(files[i], files[i].split('/')[-1])