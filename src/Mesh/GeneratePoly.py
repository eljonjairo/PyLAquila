#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pickle

#
# Generate model's output files
#
#
# John Diaz June 2023
#

import Model
import numpy as np

if __name__ == '__main__':

    print()
    print(" ************************************************ ")
    print(" *        Starting Generate Model Program       * ")
    print(" ************************************************ ")
    print()

    # Fault Name
    name = 'LAquilaCirella03_eff_dhF1000m'

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

    fault_file = "../Fault/Outputs/" + name + ".pickle"
    # Load the Fault Object
    with open(fault_file, 'rb') as in_fault:
        Fault = pickle.load(in_fault)

    # Creates a new instance of the Model class
    LAquilaModel = Model.Model(Fault.name, x_min, x_max, y_min, y_max)
    print(LAquilaModel)
    LAquilaModel.generate_region_topo(zone_number=33, zone_letter='N')
    LAquilaModel.generate_model_topo(dh_topo)





