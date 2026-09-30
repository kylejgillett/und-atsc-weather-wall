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
| `:04` | 03, 09, 15, 21 UTC | `upload_gfs_forecasts.sh` |
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

    # ============================================================
    # UND ATMOSPHERIC SCIENCES WEATHER WALL
    # OPERATIONAL SCHEDULE
    # ============================================================

    SHELL=/bin/bash
    CRON_TZ=UTC

    REPO=/home/kyle.gillett/und-atsc-weather-wall
    LOGDIR=/home/kyle.gillett/weather_wall_logs


    # Local NEXRAD -> Regional RAP
    */15 * * * * $REPO/upload/run_job.sh upload_local_nexrad_analysis.sh upload_regional_rap_analysis.sh >> $LOGDIR/local_regional_analysis.log 2>&1


    # GFS forecasts: 02:04, 08:04, 14:04, 20:04 UTC
    4 3,9,15,21 * * * $REPO/upload/run_job.sh upload_gfs_forecasts.sh >> $LOGDIR/gfs_forecasts.log 2>&1


    # Observed soundings: 01:20, 13:20, 19:20 UTC
    20 1,13,19 * * * $REPO/upload/run_job.sh upload_soundings_obs.sh >> $LOGDIR/soundings_obs.log 2>&1


    # CONUS RAP -> HRRR -> Outlooks
    34 * * * * $REPO/upload/run_job.sh upload_conus_rap_analysis.sh upload_hrrr_forecasts.sh upload_outlooks.sh >> $LOGDIR/hourly_forecasts.log 2>&1


    # ASOS -> BUFKIT
    50 * * * * $REPO/upload/run_job.sh upload_asos.sh upload_soundings_bufkit.sh >> $LOGDIR/asos_bufkit.log 2>&1