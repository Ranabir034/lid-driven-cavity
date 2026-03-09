#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# mesh.py

import numpy as np
from constant import GRID_X, GRID_Y, LENGTH_X, LENGTH_Y, LID_VELOCITY

def initialize_mesh():
    dx = LENGTH_X / GRID_X
    dy = LENGTH_Y / GRID_Y
    x_coords = np.linspace(dx/2, LENGTH_X - dx/2, GRID_X)
    y_coords = np.linspace(dy/2, LENGTH_Y - dy/2, GRID_Y)
    mesh_x, mesh_y = np.meshgrid(x_coords, y_coords)
    
    vel_x = np.zeros((GRID_Y, GRID_X + 1))
    vel_y = np.zeros((GRID_Y + 1, GRID_X))
    pressure = np.zeros((GRID_Y, GRID_X))
    vel_x[-1, :] = LID_VELOCITY
    
    return dx, dy, x_coords, y_coords, mesh_x, mesh_y, vel_x, vel_y, pressure