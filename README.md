# University of North Dakota - Department of Atmospheric Sciences - Leon F. Osborne Weather Center Mapwall Materials

Server side meteorological analysis figure generation and uploading scripts for the Leon F. Osborne Weather Center Mapwall. 

<br>

An electronic version of the mapwall may  be viewed at: --under construction--

<br>

Direct questions to Kyle Gillett (kyle.gillett@und.edu)

<br>

### Figure Update Schedule
| Time | Schedule | Scripts |
|---|---|---|
| `:00` | Every 15 minutes | NEXRAD Analysis, Regional RAP-surface analysis |
| `:04` | 03, 09, 15, 21 UTC | GFS Forecast |
| `:15` | Every 15 minutes | NEXRAD Analysis, Regional RAP-surface analysis |
| `:20` | 01, 13, 19 UTC | Regional Observed Soundings |
| `:30` | Every 15 minutes |NEXRAD Analysis, Regional RAP-surface analysis|
| `:34` | Every hour | CONUS RAP analysis, HRRR forecasts, NOAA Outlooks |
| `:45` | Every 15 minutes | NEXRAD Analysis, Regional RAP-surface analysis |
| `:50` | Every hour | Regional ASOS, BUFKIT analysis soundings |

<br>
<br>

## Production Deployment

The production Weather Wall VM tracks the `weatherwall-vm` branch. Development is generally performed on `master`. Changes intended for production are selectively merged or copied into `weatherwall-vm`.

On the production VM:

```bash
cd ~/und-atsc-weather-wall
git pull --ff-only origin weatherwall-vm
./systemd/install.sh
systemctl --user daemon-reload
```

<br>
<br>

## Figure Frame Counts

Number of figures ("frames") each builder produces per run, and per day at the schedule above. Counts marked *up to* depend on data being available (e.g. an outlook that isn't issued, or a sounding site with no launch).

| Builder | Frames per run | Runs per day | Frames per day |
|---|---|---|---|
| GFS forecast (`build_conus_gfs_forecast.py`) | 160 | 4 | 640 |
| HRRR forecast (`build_hrrr_forecast.py`) | 19 | 24 | 456 |
| NOAA outlooks (`build_outlooks.py`) | up to 19 | 24 | up to 456 |
| CONUS RAP analysis (`build_conus_rap_analysis.py`) | 11 | 24 | 264 |
| Regional ASOS timeseries (`build_asos_timeseries.py`) | 5 | 24 | 120 |
| Local NEXRAD analysis (`build_local_nexrad_analysis.py`) | 1 | 96 | 96 |
| Regional RAP surface analysis (`build_regional_sfc_analysis.py`) | 1 | 96 | 96 |
| BUFKIT soundings (`build_bufkit_soundings.py`) | 2 | 24 | 48 |
| Observed soundings (`build_obs_soundings.py`) | up to 6 | 3 | up to 18 |
| GFK 7-day meteogram (`build_gfk_forecast_meteogram.py`) | 1 | 4 | 4 |
| **Total** | | | **up to 2,198** |

**GFS forecast:** 8 map types x 20 forecast hours = 160 frames per run. Forecast hours are F006-F072 every 6 hours (12 frames per type) and F084-F168 every 12 hours (8 frames per type).

| Type | Map | Frames per run |
|---|---|---|
| `gfs_300a` | 300 hPa heights, wind | 20 |
| `gfs_300b` | 2 PVU pressure, potential temperature, wind | 20 |
| `gfs_500a` | 500 hPa heights, wind | 20 |
| `gfs_500b` | 500 hPa heights, relative vorticity, wind | 20 |
| `gfs_850a` | 850 hPa heights, temperature, wind | 20 |
| `gfs_850b` | 850 hPa heights, temperature advection, frontogenesis, wind | 20 |
| `gfs_000a` | MSLP, 2 m temperature, 10 m wind | 20 |
| `gfs_000b` | MSLP, precipitation type / reflectivity, 1000-500 hPa thickness, 10 m wind | 20 |

**HRRR forecast:** 1 surface map (reflectivity, cloud cover, MSLP, forecast station plot) for every hour F00-F18 = 19 frames per run.

**CONUS RAP analysis:** 11 frames per run.

| Type | Map |
|---|---|
| `300a` | 300 hPa heights, wind |
| `500a` | 500 hPa heights, wind |
| `500b` | 500 hPa heights, relative vorticity, wind |
| `700a` | 700 hPa heights, temperature, wind |
| `850a` | 850 hPa heights, wind |
| `850b` | 850 hPa heights, temperature, wind |
| `850c` | 850 hPa heights, temperature advection, frontogenesis, wind |
| `000a` | MSLP, temperature, wind |
| `000b` | MSLP, dewpoint, wind |
| `000c` | MSLP, MLCAPE, MLCIN |
| `000d` | MSLP, surface observations, WPC fronts, reflectivity mosaic, GOES visible |

**NOAA outlooks:** up to 19 frames per run, one per product.

| Source | Products | Frames |
|---|---|---|
| SPC convective outlook | Days 1-3 | 3 |
| SPC fire weather outlook | Days 1-2 | 2 |
| WPC excessive rainfall outlook | Days 1-3 | 3 |
| CPC 6-10 day outlook | Temperature, precipitation | 2 |
| NHC tropical weather outlook | 7 day | 1 |
| U.S. Drought Monitor | Current, flash drought | 2 |
| WPC winter storm outlook | Snow days 1-3, ice days 1-3 | 6 |

**Regional ASOS timeseries:** 1 frame per station, 5 stations (GFK, FAR, DVL, TVF, 96D).

**BUFKIT soundings:** KGFK HRRR analysis sounding and KGFK HRRR composite sounding, 2 frames per run.

**Observed soundings:** 1 frame per site, 6 sites (KBIS, KABR, KGGW, KUNR, KINL, KMPX).

**Not currently scheduled:** `build_conus_headlines.py` (1 frame) and `build_conus_spc_outlooks.py` (3 frames, Days 1-3) are not run by any upload script.
