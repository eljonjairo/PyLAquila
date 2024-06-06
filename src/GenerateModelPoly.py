#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pickle
from pathlib import Path
from Classes.Model import Model
import numpy as np

if __name__ == '__main__':
    print()
    print(" ************************************************ ")
    print(" *        Starting Generate Model Program       * ")
    print(" ************************************************ ")
    print()

    # Fault name object
    name = 'LAquilaCirella03_eff_dhF1000m.pickle'

    # X and Y UTM limits (Km) for Central Italy Model 13 Stations
    x_min = 300.0
    x_max = 420.0
    y_min = 4610.0
    y_max = 4750.0
    z_min = -60.0

    # Horizontal and Topography steps (km)
    dh = 4.0
    dh_topo = 0.5

    # Distance to remove points around to the fault (km)
    dist_to_fault = 1.0

    # Z Coordinates of Velocity Model Layers (km) (distance between layers
    # no much bigger than dh)
    # Herrmann Velocity layers in z: 1.5, 4.5, 7.5, 14.5, 29.5, 35.5, 43.5
    # Last is the limit of the Domain
    z_layer = np.array([1.5, 4.5, 7.5, 11.0, 14.5, 19.5, 24.5, 29.5, 35.5,
                        39.5, 43.5, 48.5, 53.5, 60.0])

    src_dir = Path('')
    in_file = src_dir / 'Outputs/Fault/' / name

    with open(in_file, 'rb') as f:
        in_fault = pickle.load(f)

    LaquilaModel = Model(in_fault.name, x_min, x_max, y_min, y_max)
    print(LaquilaModel)

    LaquilaModel.generate_region_topo(33, 'N')
    LaquilaModel.generate_model_topo(dh_topo)
    LaquilaModel.generate_rectangular_model(dh, z_layer, in_fault, 
                                            dist_to_fault)
    # LaquilaModel.plot_topo()
    # LaquilaModel.plot_rectangular()
    LaquilaModel.triangulate_3d_model(in_fault)



