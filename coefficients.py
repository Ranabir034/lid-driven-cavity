#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# coefficients.py

# from constant import GRID_X, GRID_Y

def compute_coefficients(dx, dy, reynolds):
    viscosity = 1.0 / reynolds
    coeff_east_x = coeff_west_x = viscosity * dy / dx
    coeff_north_x = coeff_south_x = viscosity * dx / dy
    coeff_center_x = max(coeff_east_x + coeff_west_x + coeff_north_x + coeff_south_x, 1e-8)
    
    coeff_east_y = coeff_west_y = viscosity * dy / dx
    coeff_north_y = coeff_south_y = viscosity * dx / dy
    coeff_center_y = max(coeff_east_y + coeff_west_y + coeff_north_y + coeff_south_y, 1e-8)
    
    return {
        'x': {'east': coeff_east_x, 'west': coeff_west_x, 'north': coeff_north_x,
              'south': coeff_south_x, 'center': coeff_center_x},
        'y': {'east': coeff_east_y, 'west': coeff_west_y, 'north': coeff_north_y,
              'south': coeff_south_y, 'center': coeff_center_y},
        'viscosity': viscosity
    }