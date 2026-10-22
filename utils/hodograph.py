
##########################################################
#              HODOGRAPH MAP PLOT UTILS
#  (c) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2026
##########################################################

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.collections import LineCollection




# hodograph colors 
HODO_COLORS = [
    "#7918B6",  # 0–1 km: purple
    "#DA161A",  # 1–3 km: red
    "#FF6730",  # 3–6 km: red-orange
    "#FFC400",  # 6–9 km: orange
]





def plot_hodograph(hax, u, v, height_agl, storm_motion=None, max_wind=50, ring_spacing=25, show_storm_motion=True):

    '''
    function that handles hodograph plotting and design
    '''

    # PERPARE AND CHECK DATA
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)
    z = np.asarray(height_agl, dtype=float)


    # HODOGRAPH WIND SEGMENTS
    points = np.column_stack([u, v])
    segments = np.stack([points[:-1], points[1:]], axis=1)
    zmid = 0.5 * (z[:-1] + z[1:]) / 1000.0

    colors = np.select([zmid < 1, zmid < 3, zmid < 6], HODO_COLORS[:3], default=HODO_COLORS[3])
    collection = LineCollection(segments, colors=colors, linewidths=2.25, zorder=5)
    hax.add_collection(collection)


    # STORM MOTION
    if show_storm_motion and storm_motion is not None:
        sm_u, sm_v = storm_motion
        if np.isfinite(sm_u) and np.isfinite(sm_v):
            hax.scatter(sm_u, sm_v, marker="x", s=13, color="white", linewidths=0.8, zorder=8)


    # HODOGRAPH SPEED RINGS, LINES, AND AXES PROPERTIES
    hax.set_xlim(-max_wind, max_wind)
    hax.set_ylim(-max_wind, max_wind)
    hax.set_aspect("equal", adjustable="box")
    ticks = np.arange(ring_spacing, max_wind + 1, ring_spacing)
    for radius in ticks:
        circle = plt.Circle((0, 0), radius,fill=False, color="#E6E0D8", linewidth=0.8,alpha=0.2, zorder=1)
        hax.add_patch(circle)

    hax.axhline(0, color="#E6E0D8", linewidth=0.65, alpha=0.2,zorder=2)
    hax.axvline(0, color="#E6E0D8", linewidth=0.65, alpha=0.2, zorder=2)
    
    for spine in hax.spines.values():
        spine.set_visible(False)

    hax.set_facecolor("none")
    hax.patch.set_alpha(0)
    hax.set_xticks([])
    hax.set_yticks([])

    for spine in hax.spines.values():
        spine.set_color("white")
        spine.set_linewidth(0.45)

    hax.set_clip_on(False)







def build_hybrid_wind_profile(u_hag, v_hag, z_hag, 
                              u_iso, v_iso, z_iso, t_iso, p_iso, 
                              p_sfc, t_sfc, max_height=9000):
    """
    Construct a wind profile from RAP's native iso vars to hag vars
    Returns: Profile height AGL (m), u wind (m/s), v wind (m/s).
    """

    # constants
    Rd = 287.05
    g = 9.80665



    # prep height-above-ground (hag) vars
    z_hag = np.asarray(z_hag, dtype=float).ravel()
    u_hag = np.asarray(u_hag, dtype=float).ravel()
    v_hag = np.asarray(v_hag, dtype=float).ravel()

    valid_hag = (np.isfinite(z_hag) & np.isfinite(u_hag) & np.isfinite(v_hag) & (z_hag >= 0) & (z_hag <= max_height))

    z_hag = z_hag[valid_hag]
    u_hag = u_hag[valid_hag]
    v_hag = v_hag[valid_hag]




    # estimate height of terrain with lowest p level and hypsometeric estimation
    p_iso = np.asarray(p_iso, dtype=float).ravel()
    z_iso = np.asarray(z_iso, dtype=float).ravel()
    u_iso = np.asarray(u_iso, dtype=float).ravel()
    v_iso = np.asarray(v_iso, dtype=float).ravel()
    t_iso = np.asarray(t_iso, dtype=float).ravel()

    valid_iso = (np.isfinite(p_iso) & np.isfinite(z_iso) & np.isfinite(u_iso) & np.isfinite(v_iso) & np.isfinite(t_iso) & (p_iso < p_sfc) & (p_iso > 0) & (t_iso > 150))

    if not np.any(valid_iso):
        return z_hag, u_hag, v_hag

    # The lowest valid pressure level is the largest pressure
    # below the local surface pressure.
    lowest_idx = np.where(valid_iso)[0][np.argmax(p_iso[valid_iso])]
    # Hypsometric estimate of the surface elevation.
    t_mean = 0.5 * (t_sfc + t_iso[lowest_idx])
    terrain_est = (z_iso[lowest_idx]- (Rd * t_mean / g) * np.log(p_sfc / p_iso[lowest_idx]))
    # Convert isobaric heights from MSL to estimated AGL.
    z_iso_agl = z_iso - terrain_est
    valid_iso &= (np.isfinite(z_iso_agl) & (z_iso_agl >= 0) & (z_iso_agl <= max_height))
    z_iso_agl = z_iso_agl[valid_iso]
    u_iso = u_iso[valid_iso]
    v_iso = v_iso[valid_iso]


    # combine and clean iso/hag profiles
    z = np.concatenate([z_hag, z_iso_agl])
    u = np.concatenate([u_hag, u_iso])
    v = np.concatenate([v_hag, v_iso])
    valid = (np.isfinite(z) & np.isfinite(u) & np.isfinite(v))
    z, u, v = z[valid], u[valid], v[valid]
    order = np.argsort(z)
    z, u, v = z[order], u[order], v[order]
    z_unique, unique_idx = np.unique(z, return_index=True)
    return z_unique, u[unique_idx], v[unique_idx]










def add_hodograph_network(ax, raw_data, lons, lats, cape, stride=35, inset_size=0.055, max_wind=60, 
                          max_height=9000, cape_threshold=50, show_storm_motion=True):

    '''
    function that actually adds all hodographs to map ax
    '''

    # Isobaric profiles: time, pressure, y, x.
    u_iso = np.asarray(raw_data["u-component_of_wind_isobaric"].values[0])
    v_iso = np.asarray(raw_data["v-component_of_wind_isobaric"].values[0])
    z_iso = np.asarray(raw_data["Geopotential_height_isobaric"].values[0])

    # Native height-above-ground winds
    u_hag_all = np.asarray(raw_data["u-component_of_wind_height_above_ground"].values[0])
    v_hag_all = np.asarray(raw_data["v-component_of_wind_height_above_ground"].values[0])
    z_hag_all = np.asarray(raw_data["height_above_ground"].values).ravel()
    
    # Isobaric temperature and pressure coordinates
    t_iso_all = np.asarray(raw_data["Temperature_isobaric"].values[0])
    p_iso_all = np.asarray(raw_data["isobaric"].values).ravel()

    # Surface pressure and the lowest available AGL temperature
    p_sfc_all = np.asarray(raw_data["Pressure_surface"].values[0]).squeeze()
    t_hag_all = np.asarray(raw_data["Temperature_height_above_ground"].values[0])

    # Bunkers storm motion, converted from m/s to knots.
    sm_u = np.asarray(raw_data["U-Component_Storm_Motion_height_above_ground_layer"].values[0, 0])
    sm_v = np.asarray(raw_data["V-Component_Storm_Motion_height_above_ground_layer"].values[0, 0])
    lon_values = np.asarray(lons)
    lat_values = np.asarray(lats)
    cape = np.asarray(cape)
    ny, nx = cape.shape

    # Convert map coordinates to axes-fraction coordinates.
    # This keeps inset placement aligned with the map projection.
    source_crs = ccrs.PlateCarree()

    fig = ax.figure
    added_axes = []

    for j in range(stride // 2, ny, stride):
        for i in range(stride // 2, nx, stride):
            if cape_threshold is not None:
                if cape[j, i] < cape_threshold:
                    continue
            lon = lon_values[j, i]
            lat = lat_values[j, i]

            # Transform geographic location into map coordinates.
            projected = ax.projection.transform_point(lon, lat, source_crs)
            display_xy = ax.transData.transform(projected)
            axes_xy = ax.transAxes.inverted().transform(display_xy)
            xfrac, yfrac = axes_xy

            # Skip locations outside the displayed map.
            if not (0.02 < xfrac < 0.98):
                continue
            if not (0.02 < yfrac < 0.98):
                continue

           # Construct the hybrid wind profile at this grid point.
            t_sfc = t_hag_all[0, j, i]

            z, u, v = build_hybrid_wind_profile(u_hag=u_hag_all[:, j, i], v_hag=v_hag_all[:, j, i], z_hag=z_hag_all,
                                                 u_iso=u_iso[:, j, i], v_iso=v_iso[:, j, i], z_iso=z_iso[:, j, i], t_iso=t_iso_all[:, j, i], p_iso=p_iso_all,
                                                   p_sfc=p_sfc_all[j, i], t_sfc=t_sfc, max_height=max_height)

            # Convert the combined profile to knots.
            u = u * 1.94384
            v = v * 1.94384

            # Interpolate only inside the observed profile coverage.
            z_regular = np.arange(0, max_height + 1, 250)
            inside = ((z_regular >= z.min()) & (z_regular <= z.max()))
            z_regular = z_regular[inside]
            u_regular = np.interp(z_regular, z, u)
            v_regular = np.interp(z_regular, z, v)

            # Convert the local Bunkers motion to knots.
            storm_motion = None
            if show_storm_motion:
                storm_motion = (float(sm_u[j, i] * 1.94384), float(sm_v[j, i] * 1.94384),)

            # Position the inset using axes-relative coordinates.
            left = xfrac - inset_size / 2
            bottom = yfrac - inset_size / 2

            hax = ax.inset_axes([left, bottom, inset_size, inset_size], transform=ax.transAxes, zorder=20)
            plot_hodograph( hax, u_regular,v_regular, z_regular, storm_motion=storm_motion, max_wind=max_wind, show_storm_motion=show_storm_motion)
            added_axes.append(hax)

    return added_axes
