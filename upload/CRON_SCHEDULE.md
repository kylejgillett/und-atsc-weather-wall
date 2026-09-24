# Weather Wall Cron Job Scheduling

All Weather Wall uploader jobs are executed with:

    upload/run_job.sh

Logs are found at:

    /home/kjgill/weather_wall/logs

<br>

<br>

<br>



## Schedule


| Time | Schedule | Scripts |
|---|---|---|
| `:00` | Every 15 minutes | `upload_local_nexrad_analysis.sh` + `upload_regional_rap_analysis.sh` |
| `:04` | 02, 08, 14, 20 UTC | `upload_gfs_forecasts.sh` |
| `:15` | Every 15 minutes | `upload_local_nexrad_analysis.sh` + `upload_regional_rap_analysis.sh` |
| `:20` | 01, 13, 19 UTC | `upload_soundings_obs.sh` |
| `:30` | Every 15 minutes | `upload_local_nexrad_analysis.sh` + `upload_regional_rap_analysis.sh` |
| `:34` | Every hour | `upload_conus_rap_analysis.sh` + `upload_hrrr_forecasts.sh` + `upload_outlooks.sh` |
| `:45` | Every 15 minutes | `upload_local_nexrad_analysis.sh` + `upload_regional_rap_analysis.sh` |
| `:50` | Every hour | `upload_asos.sh` + `upload_soundings_bufkit.sh` |

<br>

<br>

<br>

## Crontab

    SHELL=/bin/bash
    CRON_TZ=UTC

    REPO=/home/kjgill/weather_wall/und-atsc-weather-wall
    LOGDIR=/home/kjgill/weather_wall/logs


    # ------------------------------------------------------------
    # NEXRAD + REGIONAL SURFACE BUNDLE
    # Every 15 minutes
    # ~3:15 total runtime
    # ------------------------------------------------------------

    */15 * * * * $REPO/upload/run_job.sh upload_local_nexrad_analysis.sh upload_regional_rap_analysis.sh >> $LOGDIR/local_regional_analysis.log 2>&1


    # ------------------------------------------------------------
    # GFS FORECASTS
    # 02, 08, 14, 20 UTC
    # ~15 minute runtime
    # ------------------------------------------------------------

    4 3,9,15,21 * * * $REPO/upload/run_job.sh upload_gfs_forecasts.sh >> $LOGDIR/gfs_forecasts.log 2>&1


    # ------------------------------------------------------------
    # OBSERVED SOUNDINGS
    # 01, 13, 19 UTC
    # ~4-6 minute runtime
    # ------------------------------------------------------------

    20 1,13,19 * * * $REPO/upload/run_job.sh upload_soundings_obs.sh >> $LOGDIR/soundings_obs.log 2>&1


    # ------------------------------------------------------------
    # CONUS RAP + HRRR FORECAST + OUTLOOKS BUNDLE
    #
    # CONUS RAP  ~3:45
    # HRRR       ~2:30
    # Outlooks   ~2:40
    # Total      ~8:55
    # ------------------------------------------------------------

    34 * * * * $REPO/upload/run_job.sh upload_conus_rap_analysis.sh upload_hrrr_forecasts.sh upload_outlooks.sh >> $LOGDIR/hourly_forecasts.log 2>&1


    # ------------------------------------------------------------
    # ASOS + BUFKIT BUNDLE
    #
    # ASOS       ~0:35
    # BUFKIT     ~0:40
    # Total      ~1:15
    # ------------------------------------------------------------

    50 * * * * $REPO/upload/run_job.sh upload_asos.sh upload_soundings_bufkit.sh >> $LOGDIR/asos_bufkit.log 2>&1