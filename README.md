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