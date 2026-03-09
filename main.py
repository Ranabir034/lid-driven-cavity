#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 8 12:16:07 2025

@author: ranabir034
"""
# main.py

from constant import REYNOLDS_LIST
from mesh import initialize_mesh
from coefficients import compute_coefficients
from solver import run_simple_loop

def main():
    for reynolds in REYNOLDS_LIST:
        print(f"\nRunning simulation for Re = {reynolds}")
        dx, dy, x_coords, y_coords, mesh_x, mesh_y, vel_x, vel_y, pressure = initialize_mesh()
        coeffs = compute_coefficients(dx, dy, reynolds)
        run_simple_loop(vel_x, vel_y, pressure, dx, dy, coeffs, x_coords, y_coords, mesh_x, mesh_y, reynolds)

if __name__ == "__main__":
    main()