#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# plotting.py

import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from constant import GRID_X, GRID_Y

def plot_results(mesh_x, mesh_y, x_coords, y_coords, pressure, vel_x_center, vel_y_center, reynolds):
    plt.figure(figsize=(8, 6))
    contour_plot = plt.contourf(mesh_x, mesh_y, pressure, levels=25, cmap='jet', alpha=0.8)
    plt.colorbar(contour_plot, label='Pressure')
    plt.streamplot(x_coords, y_coords, vel_x_center, vel_y_center, density=2.0, color='black', linewidth=1)
    plt.title(f"Streamlines with Pressure Contours (Re = {reynolds})")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis('equal')
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(f"streamlines_Re{reynolds}.png")
    plt.close()

def plot_u_centerline(vel_x_center, y_coords, reynolds):
    center_x_idx = GRID_X // 2
    vx_centerline = vel_x_center[:, center_x_idx]
    
    # Ghia et al. (1982) data for different Reynolds numbers
    ghia_data = {
        100: {
            'y': [1.0000, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5000, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.0000],
            'vx': [1.0000, 0.8412, 0.7887, 0.7372, 0.6872, 0.2315, 0.0033, -0.1364, -0.2058, -0.2109, -0.1566, -0.1015, -0.0643, -0.0478, -0.0419, -0.0372, 0.0000]
        },
        400: {
            'y': [1.0000, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5000, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.0000],
            'vx': [1.00000, 0.75837, 0.68439, 0.61756, 0.55892, 0.29093, 0.16256, 0.02135, -0.11477, -0.17119, -0.32726, -0.24299, -0.14612, -0.10338, -0.09266, -0.08186, 0.00000]
        },
        1000: {
            'y': [1.0000, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5000, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.0000],
            'vx': [1.00000, 0.65928, 0.57492, 0.51117, 0.46604, 0.33304, 0.18719, 0.05702, -0.06080, -0.10648, -0.27805, -0.38289, -0.29730, -0.22220, -0.20196, -0.18109, 0.00000]
        }
    }
    
    ghia_y = ghia_data[reynolds]['y']
    ghia_vx = ghia_data[reynolds]['vx']
    
    interp_vx = interp1d(y_coords, vx_centerline, kind='cubic', bounds_error=False, fill_value="extrapolate")
    vx_interp = interp_vx(ghia_y)
    
    plt.figure(figsize=(6, 5))
    plt.plot(y_coords, vx_centerline, label='Current Simulation')
    plt.plot(ghia_y, ghia_vx, 'o', label='Ghia et al.', markersize=4)
    plt.title(f"x-velocity at x=0.5 (Re = {reynolds})")
    plt.xlabel("y")
    plt.ylabel("u")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"u_centerline_Re{reynolds}.png")
    plt.close()

def plot_v_centerline(vel_y_center, x_coords, reynolds):
    """Plot v-velocity at y=0.5 compared with Ghia et al. (1982)."""
    center_y_idx = GRID_Y // 2
    vy_centerline = vel_y_center[center_y_idx, :]
    
    # Ghia et al. (1982) data for different Reynolds numbers
    ghia_data = {
        100: {
            'x': [1.0000, 0.9688, 0.9609, 0.9531, 0.9453, 0.8516, 0.7344, 0.6172, 0.5000, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.0000],
            'vy': [0.0000, -0.0590, -0.0739, -0.0886, -0.1031, -0.1691, -0.2245, -0.2453, -0.2453, -0.2242, -0.1325, -0.0706, -0.0323, -0.0214, -0.0183, -0.0160, 0.0000]
        },
        400: {
            'x': [1.0000, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, 0.5000, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0.0000],
            'vy': [0.00000, -0.12146, -0.15663, -0.19254, -0.22847, -0.23827, -0.44993, -0.38598, 0.05186, 0.30174, 0.30203, 0.28124, 0.22965, 0.20920, 0.19713, 0.18360, 0.00000]
        },
        1000: {
            'x': [1.0000, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, 0.5000, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0.0000],
            'vy': [0.00000, -0.21388, -0.27669, -0.33714, -0.39188, -0.51550, -0.42665, -0.31966, 0.02526, 0.32325, 0.30305, 0.37095, 0.32627, 0.30353, 0.29012, 0.27485, 0.00000]
        }
    }
    
    ghia_x = ghia_data[reynolds]['x']
    ghia_vy = ghia_data[reynolds]['vy']
    
    interp_vy = interp1d(x_coords, vy_centerline, kind='cubic', bounds_error=False, fill_value="extrapolate")
    vy_interp = interp_vy(ghia_x)
    
    plt.figure(figsize=(6, 5))
    plt.plot(x_coords, vy_centerline, label='Current Simulation')
    plt.plot(ghia_x, ghia_vy, 'o', label='Ghia et al.', markersize=4)
    plt.title(f"y-velocity at y=0.5 (Re = {reynolds})")
    plt.xlabel("x")
    plt.ylabel("v")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"v_centerline_Re{reynolds}.png")
    plt.close()

def plot_u_contour(mesh_x, mesh_y, vel_x_center, reynolds):
    plt.figure(figsize=(8, 6))
    contour_plot = plt.contourf(mesh_x, mesh_y, vel_x_center, levels=20, cmap='jet')
    plt.colorbar(contour_plot, label='u-velocity')
    plt.title(f"u-velocity Contour (Re = {reynolds})")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis('equal')
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(f"u_contour_Re{reynolds}.png")
    plt.close()

def plot_v_contour(mesh_x, mesh_y, vel_y_center, reynolds):
    plt.figure(figsize=(8, 6))
    contour_plot = plt.contourf(mesh_x, mesh_y, vel_y_center, levels=20, cmap='jet')
    plt.colorbar(contour_plot, label='v-velocity')
    plt.title(f"v-velocity Contour (Re = {reynolds})")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis('equal')
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(f"v_contour_Re{reynolds}.png")
    plt.close()