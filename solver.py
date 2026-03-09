#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# solvers.py

import numpy as np
from constant import GRID_X, GRID_Y, RELAX_VEL_X, RELAX_VEL_Y, RELAX_PRESS, LID_VELOCITY, MAX_ITERATIONS, CONVERGENCE_TOL, PRESSURE_SWEEPS
from plotting import plot_results, plot_u_centerline, plot_u_contour, plot_v_contour, plot_v_centerline

def power_law_scheme(value, peclet):
    peclet = min(max(peclet, -10), 10)
    return value * max(0, (1 - 0.1 * abs(peclet))**5)

def solve_pressure_correction(press_corr, vel_x_star, vel_y_star, dx, dy):
    for _ in range(PRESSURE_SWEEPS):
        rhs = ((vel_x_star[:, 1:] - vel_x_star[:, :-1]) / dx + 
               (vel_y_star[1:, :] - vel_y_star[:-1, :]) / dy) * dx * dy
        for j in range(1, press_corr.shape[0] - 1):
            for i in range(1, press_corr.shape[1] - 1):
                press_corr[j, i] = 0.25 * (
                    press_corr[j, i+1] + press_corr[j, i-1] +
                    press_corr[j+1, i] + press_corr[j-1, i] - rhs[j, i])
    return press_corr

def update_velocities(vel_x, vel_y, pressure, dx, dy, coeffs):
    vel_x_star = np.copy(vel_x)
    vel_y_star = np.copy(vel_y)
    
    # update x-velocity
    for j in range(1, GRID_Y-1):
        for i in range(1, GRID_X-1):
            vx_east, vx_west, vx_center = vel_x[j, i+1], vel_x[j, i-1], vel_x[j, i]
            vy_east = 0.25 * (vel_y[j, i] + vel_y[j, i-1] + vel_y[j+1, i] + vel_y[j+1, i-1])
            vy_west = 0.25 * (vel_y[j, i-1] + vel_y[j, i-2] + vel_y[j+1, i-1] + vel_y[j+1, i-2])
            flux_east = vx_east * dy
            flux_west = vx_west * dy
            diff_east = diff_west = coeffs['viscosity'] * dy / dx
            coeff_east = max(flux_east, 0) + diff_east * power_law_scheme(1, flux_east/diff_east)
            coeff_west = max(-flux_west, 0) + diff_west * power_law_scheme(1, flux_west/diff_west)
            coeff_center = max(coeff_east + coeff_west + coeffs['x']['north'] + coeffs['x']['south'], 1e-8)
            diffusion = (coeff_west * vx_west + coeff_east * vx_east + 
                        coeffs['x']['north'] * vel_x[j+1, i] + coeffs['x']['south'] * vel_x[j-1, i])
            press_term = dy * (pressure[j, i-1] - pressure[j, i])
            vel_x_star[j, i] = (1 - RELAX_VEL_X) * vx_center + RELAX_VEL_X * (diffusion + press_term) / coeff_center
    
    # update y-velocity
    for j in range(2, GRID_Y-2):
        for i in range(2, GRID_X-2):
            vy_north, vy_south, vy_center = vel_y[j+1, i], vel_y[j-1, i], vel_y[j, i]
            vx_east = 0.25 * (vel_x[j, i] + vel_x[j, i+1] + vel_x[j-1, i] + vel_x[j-1, i+1])
            vx_west = 0.25 * (vel_x[j, i] + vel_x[j, i-1] + vel_x[j-1, i] + vel_x[j-1, i-1])
            flux_north = vy_north * dx
            flux_south = vy_south * dx
            diff_north = diff_south = coeffs['viscosity'] * dx / dy
            coeff_north = max(flux_north, 0) + diff_north * power_law_scheme(1, flux_north/diff_north)
            coeff_south = max(-flux_south, 0) + diff_south * power_law_scheme(1, flux_south/diff_south)
            coeff_center = max(coeff_north + coeff_south + coeffs['y']['east'] + coeffs['y']['west'], 1e-8)
            diffusion = (coeff_south * vy_south + coeff_north * vy_north + 
                        coeffs['y']['east'] * vel_y[j, i+1] + coeffs['y']['west'] * vel_y[j, i-1])
            press_term = dx * (pressure[j-1, i] - pressure[j, i])
            vel_y_star[j, i] = (1 - RELAX_VEL_Y) * vy_center + RELAX_VEL_Y * (diffusion + press_term) / coeff_center
    
    return vel_x_star, vel_y_star

def apply_boundary_conditions(vel_x, vel_y):
    vel_x[:, 0] = vel_x[:, -1] = 0.0
    vel_x[0, :] = 0.0
    vel_x[-1, :] = LID_VELOCITY
    vel_y[0, :] = vel_y[-1, :] = 0.0
    vel_y[:, 0] = vel_y[:, -1] = 0.0
    return vel_x, vel_y

def run_simple_loop(vel_x, vel_y, pressure, dx, dy, coeffs, x_coords, y_coords, mesh_x, mesh_y, reynolds):
    for iteration in range(MAX_ITERATIONS):
        vel_x_prev = np.copy(vel_x)
        vel_y_prev = np.copy(vel_y)
        pressure_prev = np.copy(pressure)
        
        vel_x_star, vel_y_star = update_velocities(vel_x, vel_y, pressure, dx, dy, coeffs)
        
        press_corr = np.zeros_like(pressure)
        press_corr = solve_pressure_correction(press_corr, vel_x_star, vel_y_star, dx, dy)
        
        vel_x[:, 1:-1] = vel_x_star[:, 1:-1] + dy * (press_corr[:, :-1] - press_corr[:, 1:]) / coeffs['x']['center']
        vel_y[1:-1, :] = vel_y_star[1:-1, :] + dx * (press_corr[:-1, :] - press_corr[1:, :]) / coeffs['y']['center']
        pressure += RELAX_PRESS * press_corr
        
        vel_x, vel_y = apply_boundary_conditions(vel_x, vel_y)
        
        if np.isnan(vel_x).any() or np.isnan(vel_y).any():
            print(f"NaN detected at iteration {iteration} for Re = {reynolds}")
            break
        
        residual_x = np.linalg.norm(vel_x - vel_x_prev)
        residual_y = np.linalg.norm(vel_y - vel_y_prev)
        residual_p = np.linalg.norm(pressure - pressure_prev)
        print(f"Re = {reynolds}, Iter {iteration}: x_res={residual_x:.2e}, y_res={residual_y:.2e}, p_res={residual_p:.2e}")
        
        if iteration > 0 and residual_x < CONVERGENCE_TOL and residual_y < CONVERGENCE_TOL:
            print(f"Converged in {iteration} iterations for Re = {reynolds}")
            break
    
    vel_x_center = 0.5 * (vel_x[:, :-1] + vel_x[:, 1:])
    vel_y_center = 0.5 * (vel_y[:-1, :] + vel_y[1:, :])
    
    plot_results(mesh_x, mesh_y, x_coords, y_coords, pressure, vel_x_center, vel_y_center, reynolds)
    plot_u_centerline(vel_x_center, y_coords, reynolds)
    plot_v_centerline(vel_y_center, x_coords, reynolds)
    plot_u_contour(mesh_x, mesh_y, vel_x_center, reynolds)
    plot_v_contour(mesh_x, mesh_y, vel_y_center, reynolds)