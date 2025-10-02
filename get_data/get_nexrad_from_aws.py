import warnings
warnings.filterwarnings("ignore")

import os
from contextlib import contextmanager, redirect_stderr, redirect_stdout
@contextmanager
def suppress_stdout_stderr():
    """A context manager that redirects stdout and stderr to devnull"""
    with open(os.devnull, 'w') as fnull:
        with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
            yield (err, out)

with suppress_stdout_stderr():
    import pyart

import sys
import os
import fsspec
from datetime import datetime, timedelta, timezone


def get_radar(date=datetime.now(timezone.utc)):

    curr_dt = date
    time = curr_dt.strftime('%H%M')
    datestr = curr_dt.strftime('%Y%m%d')

    nexrad_site = "KMVX"

    fs = fsspec.filesystem('s3', anon=True)
    # try current time, then fallback to previous hour if needed
    time_candidates = [curr_dt + timedelta(minutes=delta) for delta in range(-60, 0)]

    files = []
    for t in time_candidates:
        prefix = t.strftime(f"s3://unidata-nexrad-level2/%Y/%m/%d/{nexrad_site}/{nexrad_site}%Y%m%d_%H%M")
        matches = fs.glob(prefix + "*")
        matches = [f for f in matches if not f.endswith('_MDM')]
        files.extend(matches)

    if files:
        files = sorted(files)
        file = files[-1]
    else:
        raise FileNotFoundError("    > RADAR ERROR: No single-site radar scan found")

    radar = pyart.io.read_nexrad_archive("s3://" + file)

    # extract timestamp from filename
    scan_time = str(file[51:53] + ':' + file[53:55])
    scan_date = str(file[42:46] + '-' + file[46:48] + '-' + file[48:50])
    scan_timestamp = str(scan_date + ' | ' + scan_time)
    print(f"    RADAR SCAN LOADED.....Site: {nexrad_site} Scan: {scan_timestamp}z")

    return radar, scan_timestamp