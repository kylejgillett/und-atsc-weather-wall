import math

screen_w, screen_h = 2160, 3840
diaginal = 55


def compute_dpi(width_px, height_px, diagonal_in):
    return ( (width_px**2 + height_px**2) ** 0.5 ) / diagonal_in


def compute_figdims(screen_w, screen_h, cols, rows, dpi):
    fig_w = (screen_w / cols) / dpi
    fig_h = (screen_h / rows) / dpi
    return fig_w, fig_h


dpi = compute_dpi(2160, 3840, 55)